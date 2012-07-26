/* 
 * File:   ReliefFSeq.cpp
 * Author: Bill White
 * 
 * Created on: 7/21/12
 */

#include <iostream>
#include <vector>

#include "ReliefF.h"
#include "ReliefFSeq.h"
#include "Dataset.h"
#include "DatasetInstance.h"
#include "DistanceMetrics.h"
#include "Insilico.h"
#include "Statistics.h"

using namespace std;

ReliefFSeq::ReliefFSeq(Dataset* ds) :
		ReliefF::ReliefF(ds, SIGNAL_TO_NOISE_ANALYSIS) {
	cout << Timestamp() << "ReliefFSeq initializing" << endl;
}

ReliefFSeq::ReliefFSeq(Dataset* ds, po::variables_map& vm) :
		ReliefF::ReliefF(ds, vm, SIGNAL_TO_NOISE_ANALYSIS) {
	cout << Timestamp() << "ReliefFSeq initializing" << endl;
}

ReliefFSeq::ReliefFSeq(Dataset* ds, ConfigMap& configMap) :
		ReliefF::ReliefF(ds, configMap, SIGNAL_TO_NOISE_ANALYSIS) {
	cout << Timestamp() << "ReliefFSeq initialing" << endl;
}

ReliefFSeq::~ReliefFSeq() {
}

bool ReliefFSeq::ComputeAttributeScores() {
	// preconditions:
	// 1. case-control data
	// 2. all numeric variables

	W.resize(dataset->NumNumerics(), 0.0);

	// precompute all instance-to-instance distances and get nearest neighbors
	PreComputeDistances();

	// using pseudo-code notation from white board discussion - 7/21/12
	// changed to use Brett's email (7/21/12) equations - 7/23/12
	cout << Timestamp() << "Running Relief-F Seq algorithm" << endl;

#pragma omp parallel for
	for(unsigned int alpha = 0; alpha < dataset->NumNumerics(); ++alpha) {

		pair<double, double> muDeltaAlphas = MuDeltaAlphas(alpha);
		double muDeltaHitAlpha = muDeltaAlphas.first;
		double muDeltaMissAlpha = muDeltaAlphas.second;

		pair<double, double> sigmaDeltaAlphas = SigmaDeltaAlphas(alpha,
				muDeltaHitAlpha, muDeltaMissAlpha);
		double sigmaDeltaHitAlpha = sigmaDeltaAlphas.first;
		double sigmaDeltaMissAlpha = sigmaDeltaAlphas.second;

		double num = fabs(muDeltaMissAlpha - muDeltaHitAlpha);
		double den = sigmaDeltaMissAlpha + sigmaDeltaHitAlpha;

		W[alpha] = num / den;
	}

	return true;
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
