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

	// convenience variables
	unsigned int numVariables = dataset->NumNumerics();
	unsigned int numInstances = dataset->NumInstances();

	// precompute all instance-to-instance distances and get nearest neighbors
	PreComputeDistances();
	// precompute all nearest neighbor hit and miss variable averages and std's
	PreComputeNeighborGeneStats();

	// variable weights
	W.resize(numVariables, 0.0);

	// using pseudo-code notation from white board discussion - 7/21/12
	for(unsigned int varIdx=0; varIdx < numVariables; ++varIdx) {
		for(unsigned int instanceIdx=0; instanceIdx < numInstances; ++instanceIdx) {
			InstanceHitMissStats hitMissStats = neighborStats[instanceIdx];

			InstanceAttributeStats varIdxHitStats = hitMissStats.first;
			InstanceAttributeStats varIdxMissStats = hitMissStats.second;

			double avgHits = varIdxHitStats[varIdx].first;
			double stdHits = varIdxHitStats[varIdx].second;

			double avgMisses = varIdxMissStats[varIdx].first;
			double stdMisses = varIdxMissStats[varIdx].second;

			W[varIdx] += (avgMisses - avgHits) / (stdMisses + stdHits);
		}
	}

	// divide all weights by the m-instance sums accumulated above
	for(unsigned int i=0; i < W.size(); ++i) {
		W[i] /= m;
	}

	return true;
}

bool SNReliefF::PreComputeNeighborGeneStats() {

	for(unsigned int i=0; i < dataset->NumInstances(); ++i) {
		DatasetInstance* M_i = dataset->GetInstance(i);
		// find k nearest hits and nearest misses
		vector<unsigned int> hits;
		map<ClassLevel, vector<unsigned int> > allMisses;
		bool canGetNeighbors = false;
		canGetNeighbors = M_i->GetNNearestInstances(k, hits, allMisses);
		if(allMisses.size() != 1) {
			cerr << "ERROR: SNRelief requires case-control data" << endl;
			return false;
		}
		vector<unsigned int> misses = allMisses.begin()->second;
		InstanceHitMissStats instanceHitMissStats;
		ComputeInstanceStats(M_i, hits, misses, instanceHitMissStats);
		neighborStats.push_back(instanceHitMissStats);
	}

	return true;
}

bool SNReliefF::ComputeInstanceStats(DatasetInstance* dsi,
		vector<unsigned int> hitIndicies, vector<unsigned int> missIndicies,
		InstanceHitMissStats& hitMissStats) {
	unsigned int numAttributes = dataset->NumNumerics();

	InstanceAttributeStats hitStats;
	vector<unsigned int>::const_iterator hitIt = hitIndicies.begin();
	for(unsigned int attributeIndex=0; attributeIndex < numAttributes;
			++attributeIndex) {
		vector<double> hitAttributeValues;
		for(; hitIt != hitIndicies.end(); ++hitIt) {
			NumericLevel thisValue =
					dataset->GetInstance(*hitIt)->GetNumeric(attributeIndex);
			hitAttributeValues.push_back(thisValue);
		}
	  double average = 0.0;
	  pair<double, double> varStd = VarStd(hitAttributeValues, average);
	  pair<double, double> avgStd = make_pair(average, varStd.second);
	  hitStats.push_back(avgStd);
	}

	InstanceAttributeStats missStats;
	vector<unsigned int>::const_iterator missIt = missIndicies.begin();
	for(unsigned int attributeIndex=0; attributeIndex < numAttributes;
			++attributeIndex) {
		vector<double> missAttributeValues;
		for(; missIt != missIndicies.end(); ++missIt) {
			NumericLevel thisValue =
					dataset->GetInstance(*missIt)->GetNumeric(attributeIndex);
			missAttributeValues.push_back(thisValue);
		}
	  double average = 0.0;
	  pair<double, double> varStd = VarStd(missAttributeValues, average);
	  pair<double, double> avgStd = make_pair(average, varStd.second);
	  missStats.push_back(avgStd);
	}

	hitMissStats.first = hitStats;
	hitMissStats.second = missStats;

	return true;
}
