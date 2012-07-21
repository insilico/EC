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

using namespace std;

SNReliefF::SNReliefF(Dataset* ds) :
		ReliefF::ReliefF(ds, REGRESSION_ANALYSIS) {
	cout << Timestamp() << "SNReliefF initialization" << endl;
	if (!ds->HasContinuousPhenotypes()) {
		cerr << "ERROR: Attempting to construct SNReliefF object without a "
				"continuous phenotype data set" << endl;
		exit(EXIT_FAILURE);
	}
}

SNReliefF::SNReliefF(Dataset* ds, po::variables_map& vm) :
		ReliefF::ReliefF(ds, vm, REGRESSION_ANALYSIS) {
	cout << Timestamp() << "SNReliefF initialization" << endl;
	if (!ds->HasContinuousPhenotypes()) {
		cerr << "ERROR: Attempting to construct SNReliefF object without a "
				"continuous phenotype data set" << endl;
		exit(EXIT_FAILURE);
	}
}

SNReliefF::SNReliefF(Dataset* ds, ConfigMap& configMap) :
		ReliefF::ReliefF(ds, configMap, REGRESSION_ANALYSIS) {
	cout << Timestamp() << "SNReliefF initialization" << endl;
	if (!ds->HasContinuousPhenotypes()) {
		cerr << "ERROR: Attempting to construct SNReliefF object without a "
				"continuous phenotype data set" << endl;
		exit(EXIT_FAILURE);
	}
}

SNReliefF::~SNReliefF() {
}

bool SNReliefF::ComputeAttributeScores() {

	// precompute all instance-to-instance distances and get nearest neighbors
	PreComputeDistances();

	// results are stored in scores
	W.resize(dataset->NumVariables(), 0.0);

	// using pseudocode notation from whiteboard discussion - 7/21/12

	return true;
}
