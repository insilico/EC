/* 
 * File:   ReliefFSeq.cpp
 * Author: Bill White
 * 
 * Created on: 7/21/12
 */

#include <iostream>
#include <iomanip>
#include <vector>

#include <gsl/gsl_cdf.h>

#include "ReliefF.h"
#include "ReliefFSeq.h"
#include "Dataset.h"
#include "DatasetInstance.h"
#include "DistanceMetrics.h"
#include "Insilico.h"
#include "Statistics.h"

using namespace std;

ReliefFSeq::ReliefFSeq(Dataset* ds) :
		ReliefF::ReliefF(ds, RNASEQ_ANALYSIS) {
	cout << Timestamp() << "ReliefFSeq initializing with data set only" << endl;
}

ReliefFSeq::ReliefFSeq(Dataset* ds, po::variables_map& vm) :
		ReliefF::ReliefF(ds, vm, RNASEQ_ANALYSIS) {
	cout << Timestamp() << "ReliefFSeq initializing with Boost parameters"
			<< endl;
	mode = "snr";
	if(vm.count("ec-seq-algorithm-mode")) {
		mode = vm["ec-seq-algorithm-mode"].as<string>();
		if((mode =="snr") || (mode == "tstat")) {
			cout << Timestamp() << "ReliefFSeq mode set to: " << mode << endl;
		}
		else {
			cerr << "ERROR: Unrecognized ReliefFSeq mode: " << mode << endl;
			exit(EXIT_FAILURE);
		}
	}
	s0 = 0.05;
	if(vm.count("ec-seq-algorithm-s0")) {
		s0 = vm["ec-seq-algorithm-s0"].as<double>();
		if((s0 >= 0) || (s0 <= 1.0)) {
			cout << Timestamp() << "ReliefFSeq s0 set to: " << s0 << endl;
		}
		else {
			cerr << "ERROR: ReliefFSeq s0 out of range (0, 1): " << s0 << endl;
			exit(EXIT_FAILURE);
		}

	}
}

ReliefFSeq::ReliefFSeq(Dataset* ds, ConfigMap& configMap) :
		ReliefF::ReliefF(ds, configMap, RNASEQ_ANALYSIS) {
	cout << Timestamp() << "ReliefFSeq initializing with config map" << endl;
}

ReliefFSeq::~ReliefFSeq() {
}

bool ReliefFSeq::ComputeAttributeScores() {
	// preconditions:
	// 1. case-control data
	// 2. all numeric variables

	W.resize(dataset->NumNumerics(), 0.0);

	// pre-compute all instance-to-instance distances and get nearest neighbors
	PreComputeDistances();

	// using pseudo-code notation from white board discussion - 7/21/12
	// changed to use Brett's email (7/21/12) equations - 7/23/12
	cout << Timestamp() << "Running ReliefFSeq algorithm" << endl;
//	vector<string> numNames;
//	numNames = dataset->GetNumericsNames();
	vector<unsigned int> numericIndices =
			dataset->MaskGetAttributeIndices(NUMERIC_TYPE);
	// DEBUG
//	string rawScoresFileName = dataset->GetNumericsFilename() + "_rawscores.tab";
//	ofstream outFile(rawScoresFileName.c_str());
//	outFile << "gene\tmuMiss\tmuHit\tsigmaMiss\tsigmaHit\tnum\tden\tdms0\tsnr" << endl;
#pragma omp parallel for
	for (unsigned int numIdx = 0; numIdx < numericIndices.size();
			++numIdx) {
		unsigned int alpha = numericIndices[numIdx];
		pair<double, double> muDeltaAlphas = MuDeltaAlphas(alpha);
		double muDeltaHitAlpha = muDeltaAlphas.first;
		double muDeltaMissAlpha = muDeltaAlphas.second;

		pair<double, double> sigmaDeltaAlphas = SigmaDeltaAlphas(alpha,
				muDeltaHitAlpha, muDeltaMissAlpha);
		double sigmaDeltaHitAlpha = sigmaDeltaAlphas.first;
		double sigmaDeltaMissAlpha = sigmaDeltaAlphas.second;

		double num = 0.0, den = 0.0;
		double alphaWeight = 0.0;
		if(mode == "snr") {
			// mode: snr (signal to noise ratio)
			num = fabs(muDeltaMissAlpha - muDeltaHitAlpha);
			den = sigmaDeltaMissAlpha + sigmaDeltaHitAlpha;
//			outFile << numNames[numIdx]
//					<< "\t" << muDeltaMissAlpha << "\t" << muDeltaHitAlpha
//					<< "\t" << sigmaDeltaMissAlpha << "\t" << sigmaDeltaHitAlpha
//					<< "\t" << num << "\t" << den << "\t" << (den + s0);
			alphaWeight = num / (den + s0);
		}
		else {
			// mode: tstat (t-statistic)
			// from Brett's email - 8/15/12
			// Also we could change the score to a real t-statistic:
			// (xbar1 – xbar2)/(Sp*sqrt(1/n1 + 1/n2)), where Sp = pooled standard
			// deviation=sqrt(((n1-1)*variance1 + (n2-1)*variance2)/(n1+n2-2)).
			double n1, n2;
			n1 = n2 = k;
			double variance1 = sigmaDeltaHitAlpha;
			double variance2 = sigmaDeltaMissAlpha;
			double pooledStdDev =
					sqrt(((n1 - 1) * variance1 + (n2 - 1) * variance2) / (n1 + n2 - 2));
			num = muDeltaMissAlpha - muDeltaHitAlpha;
			den = pooledStdDev * sqrt((1.0 / n1) + (1.0 / n2));
			double t = num / (den + s0);
			double df =  n1 + n2 - 2;
			double gslPval = 1.0;
			if(t < 0) {
				gslPval = gsl_cdf_tdist_P(-t, df);
			}
			else {
				gslPval = gsl_cdf_tdist_P(t, df);
			}
			alphaWeight = 1.0 - (2.0 * (1.0 - gslPval));
			// DEBUG
//			cout << setprecision(8);
//			cout << numIdx << "\t" << t << "\t" << gslPval << "\t" << alphaWeight << endl;
//			alphaWeight = fabs(t);
		}
		W[numIdx] = alphaWeight;
		// W[numIdx] = num; // should match ReliefF
		// DEBUG
//		outFile << "\t" << W[numIdx] << endl;
	}
	// DEBUG
//	outFile.close();

	return true;
}

AttributeScores ReliefFSeq::GetScores() {
	AttributeScores returnScores;
	vector<string> numNames;
	numNames = dataset->GetNumericsNames();
	for(unsigned int alpha = 0; alpha < dataset->NumNumerics(); ++alpha) {
		returnScores.push_back(make_pair(W[alpha], numNames[alpha]));
	}
	return returnScores;
}

pair<double, double> ReliefFSeq::MuDeltaAlphas(unsigned int alpha) {

	// for all instances
	double missSum = 0.0;
	double hitSum = 0.0;
	for(unsigned int i = 0; i < m; ++i) {
		DatasetInstance* S_i = dataset->GetInstance(i);

		// get hits and misses for this instance
		vector<unsigned int> hits(k);
		vector<unsigned int> misses(k);
		map<ClassLevel, vector<unsigned int> > allMisses;
		S_i->GetNNearestInstances(k, hits, allMisses);
		misses = allMisses.begin()->second;

		// sum over all nearest hits and misses neighbors
		for(unsigned int j = 0; j < k; ++j) {
			hitSum += diffManhattan(alpha, S_i, dataset->GetInstance(hits[j]));
			missSum += diffManhattan(alpha, S_i, dataset->GetInstance(misses[j]));
		}
	}

	// return the average of the hit and miss diffs
	pair<double, double> returnValues;
	double avgFactor = 1.0 / ((double) m * (double) k);
	returnValues.first = hitSum * avgFactor;
	returnValues.second = missSum * avgFactor;

	return returnValues;
}

pair<double, double> ReliefFSeq::SigmaDeltaAlphas(unsigned int alpha,
		double muDeltaHit, double muDeltaMiss) {

	double missSum = 0.0;
	double hitSum = 0.0;

	// for all instances
	for(unsigned int i = 0; i < m; ++i) {
		DatasetInstance* S_i = dataset->GetInstance(i);
		/// get hits and misses for this instance
		vector<unsigned int> hits(k);
		map<ClassLevel, vector<unsigned int> > allMisses;
		S_i->GetNNearestInstances(k, hits, allMisses);
		// for all nearest neighbor hits
		for(unsigned int j = 0; j < k; ++j) {
			hitSum +=
					pow(
							(diffManhattan(alpha, S_i, dataset->GetInstance(hits[j]))
									- muDeltaHit), 2);
		}
		// for all nearest neighbor misses
		// assume only one other miss class
		vector<unsigned int> misses(k);
		misses = allMisses.begin()->second;
		for(unsigned int j = 0; j < k; ++j) {
			missSum += pow(
					(diffManhattan(alpha, S_i, dataset->GetInstance(misses[j]))
							- muDeltaMiss), 2);
		}
	}

	// return the standard deviation of the hit and miss diffs
	pair<double, double> returnValues;
	double avgFactor = 1.0 / ((double) m * (double) k);
	returnValues.first = sqrt(hitSum * avgFactor);
	returnValues.second = sqrt(missSum * avgFactor);

	return returnValues;
}
