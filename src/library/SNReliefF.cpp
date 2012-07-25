/* 
 * File:   SNReliefF.cpp
 * Author: Bill White
 * 
 * Created on: 7/21/12
 */

#include <iostream>
#include <vector>

#include "ReliefF.h"
#include "SNReliefF.h"
#include "Dataset.h"
#include "DistanceMetrics.h"
#include "Insilico.h"
#include "Statistics.h"

using namespace std;

SNReliefF::SNReliefF(Dataset* ds) :
		ReliefF::ReliefF(ds, SIGNAL_TO_NOISE_ANALYSIS) {
	cout << Timestamp() << "SNReliefF initializing" << endl;
}

SNReliefF::SNReliefF(Dataset* ds, po::variables_map& vm) :
		ReliefF::ReliefF(ds, vm, SIGNAL_TO_NOISE_ANALYSIS) {
	cout << Timestamp() << "SNReliefF initializing" << endl;
}

SNReliefF::SNReliefF(Dataset* ds, ConfigMap& configMap) :
		ReliefF::ReliefF(ds, configMap, SIGNAL_TO_NOISE_ANALYSIS) {
	cout << Timestamp() << "SNReliefF initialing" << endl;
}

SNReliefF::~SNReliefF() {
}

bool SNReliefF::ComputeAttributeScores() {
	// preconditions:
	// 1. case-control data
	// 2. all numeric variables

	// precompute all instance-to-instance distances and get nearest neighbors
	PreComputeDistances();
	// precompute all nearest neighbor hit and miss variable averages and std's
	PreComputeNeighborGeneStats();
	// PrintNeighborStats();

	// variable weights
	std::vector<double> avgHitSum;
	std::vector<double> stdHitSum;
	std::vector<double> avgMissSum;
	std::vector<double> stdMissSum;

	avgHitSum.resize(dataset->NumVariables(), 0.0);
	stdHitSum.resize(dataset->NumVariables(), 0.0);
	avgMissSum.resize(dataset->NumVariables(), 0.0);
	stdMissSum.resize(dataset->NumVariables(), 0.0);

	// using pseudo-code notation from white board discussion - 7/21/12
	cout << Timestamp() << "Running SNRelief-F algorithm" << endl;

	for(unsigned int instanceIdx=0; instanceIdx < m; ++instanceIdx) {
		InstanceHitMissStats hitMissStats = neighborStats[instanceIdx];
		for(unsigned int varIdx=0; varIdx < dataset->NumVariables(); ++varIdx) {

			InstanceAttributeStats varIdxHitStats = hitMissStats.first;
			InstanceAttributeStats varIdxMissStats = hitMissStats.second;

			double avgHits = varIdxHitStats[varIdx].first;
			double stdHits = varIdxHitStats[varIdx].second;
			double avgMisses = varIdxMissStats[varIdx].first;
			double stdMisses = varIdxMissStats[varIdx].second;

			avgHitSum[varIdx] += avgHits;
			stdHitSum[varIdx] += stdHits;
			avgMissSum[varIdx] += avgMisses;
			stdMissSum[varIdx] += stdMisses;
		}
	}

	// divide all weights by the m-instance sums accumulated above
	W.resize(dataset->NumVariables(), 0.0);
	for(unsigned int i=0; i < dataset->NumVariables(); ++i) {
		W[i] = ((avgMissSum[i] -avgHitSum[i]) / (stdMissSum[i] + stdHitSum[i])) / m;
	}

	return true;
}

bool SNReliefF::PreComputeNeighborGeneStats() {

	cout << Timestamp() << "Precomputing nearest neighbor attribute stats" << endl;
	for(unsigned int i=0; i < dataset->NumInstances(); ++i) {
		DatasetInstance* M_i = dataset->GetInstance(i);
		// find k nearest hits and nearest misses
		vector<unsigned int> hits(k);
		map<ClassLevel, vector<unsigned int> > allMisses;
		bool canGetNeighbors = false;
		canGetNeighbors = M_i->GetNNearestInstances(k, hits, allMisses);
		if(allMisses.size() != 1) {
			cerr << "ERROR: SNReliefF requires case-control data" << endl;
			return false;
		}
		vector<unsigned int> misses(k);
		misses = allMisses.begin()->second;
		InstanceHitMissStats instanceHitMissStats;
		ComputeInstanceStats(M_i, hits, misses, instanceHitMissStats);
		neighborStats.push_back(instanceHitMissStats);
	}

	return true;
}

bool SNReliefF::ComputeInstanceStats(DatasetInstance* dsi,
		vector<unsigned int> hitIndicies, vector<unsigned int> missIndicies,
		InstanceHitMissStats& hitMissStats) {
	unsigned int numNumerics = dataset->NumNumerics();

	// get average/std for k values for each attribute in hits set
	InstanceAttributeStats hitStats;
	for(unsigned int numericIndex=0; numericIndex < numNumerics;
			++numericIndex) {
		vector<double> hitAttributeValues;
		vector<unsigned int>::const_iterator hitIt = hitIndicies.begin();
		for(; hitIt != hitIndicies.end(); ++hitIt) {
			NumericLevel thisValue =
					dataset->GetInstance(*hitIt)->GetNumeric(numericIndex);
			hitAttributeValues.push_back(thisValue);
		}
	  double average = 0.0;
	  pair<double, double> varStd = VarStd(hitAttributeValues, average);
	  pair<double, double> avgStd = make_pair(average, varStd.second);
	  hitStats.push_back(avgStd);
	}

	// misses
	InstanceAttributeStats missStats;
	for(unsigned int numericIndex=0; numericIndex < numNumerics;
			++numericIndex) {
		vector<double> missAttributeValues;
		vector<unsigned int>::const_iterator missIt = missIndicies.begin();
		for(; missIt != missIndicies.end(); ++missIt) {
			NumericLevel thisValue =
					dataset->GetInstance(*missIt)->GetNumeric(numericIndex);
			missAttributeValues.push_back(thisValue);
		}
	  double average = 0.0;
	  pair<double, double> varStd = VarStd(missAttributeValues, average);
	  pair<double, double> avgStd = make_pair(average, varStd.second);
	  // cout << average << ", " << varStd.second << endl;
	  missStats.push_back(avgStd);
	}
//	cout << "Miss stats:" << endl;
//	PrintInstanceAttributeStats(missStats);

	hitMissStats.first.resize(numNumerics);
	hitMissStats.first = hitStats;
	hitMissStats.second.resize(numNumerics);
	hitMissStats.second = missStats;

	return true;
}

void SNReliefF::PrintNeighborStats() {
	NeighborStatsCIt nit = neighborStats.begin();
	for(; nit != neighborStats.end(); ++nit) {
		cout << "********************************************************" << endl;
		InstanceHitMissStats thisInstanceHitMissStats = *nit;

		InstanceAttributeStats hitStats = thisInstanceHitMissStats.first;
		cout << "HITS" << endl;
		PrintInstanceAttributeStats(hitStats);
		cout << endl;

		InstanceAttributeStats missStats = thisInstanceHitMissStats.second;
		cout << "MISSES" << endl;
		PrintInstanceAttributeStats(missStats);
		cout << endl;
	}
}

void SNReliefF::PrintInstanceAttributeStats(InstanceAttributeStats stats) {

	InstanceAttributeStatsIt it = stats.begin();
	cout << "Avg" << "\t" << "Std" << endl;
	for(; it != stats.end(); ++it) {
		cout << it->first << "\t" << it->second << endl;
	}
	cout << endl;


}
