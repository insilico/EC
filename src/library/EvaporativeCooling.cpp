/* 
 * File:   EvaporativeCooling.cpp
 * Author: billwhite
 * 
 * Created on July 14, 2011, 9:25 PM
 *
 * Implements the Evaporative Cooling algorithm in:
 * McKinney, et. al. "Capturing the Spectrum of Interaction Effects in Genetic
 * Association Studies by Simulated Evaporative Cooling Network Analysis."
 * PLoS Genetics, Vol 5, Issue 3, 2009.
 */

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <time.h>

#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <omp.h>
#include <gsl/gsl_rng.h>

#include "EvaporativeCooling.h"

// ReliefF project
#include "Dataset.h"
#include "Statistics.h"
#include "StringUtils.h"
#include "RandomJungle.h"
#include "ReliefF.h"
#include "RReliefF.h"
#include "Insilico.h"

using namespace std;
namespace po = boost::program_options;
using namespace boost;
using namespace insilico;

bool scoresSortAsc(const pair<double, string>& p1,
const pair<double, string>& p2) {
	return p1.first < p2.first;
}

bool scoresSortAscByName(const pair<double, string>& p1,
const pair<double, string>& p2) {
	return p1.second < p2.second;
}

bool scoresSortDesc(const pair<double, string>& p1,
const pair<double, string>& p2) {
	return p1.first > p2.first;
}

EvaporativeCooling::EvaporativeCooling(Dataset* ds, po::variables_map& vm,
		AnalysisType anaType) {
	cout << Timestamp() << "Evaporative Cooling initialization:" << endl;
	if (ds) {
		dataset = ds;
	} else {
		cerr << "ERROR: data set is not initialized" << endl;
		exit(EXIT_FAILURE);
	}
	paramsMap = vm;
	analysisType = anaType;

	reliefF = NULL;
	randomJungle = NULL;

	// set the number of target attributes
	numTargetAttributes = vm["ec-num-target"].as<unsigned int>();
	if (numTargetAttributes == 0) {
		numTargetAttributes = ds->NumVariables();
		numToRemovePerIteration = 0;
	}
	if (numTargetAttributes > dataset->NumVariables()) {
		cerr << "--ec-num-taget must be less than or equal to the "
				<< "number of attributes in the data set" << endl;
		exit(EXIT_FAILURE);
	}
	cout << Timestamp() << "EC is removing attributes until best "
			<< numTargetAttributes << " remain" << endl;

	if (paramsMap.count("ec-algorithm-steps")) {
		string ecAlgParam = to_upper(vm["ec-algorithm-steps"].as<string>());
		if (ecAlgParam == "ALL") {
			algorithmType = EC_ALL;
			cout << Timestamp()
					<< "Running EC in standard mode: Random Jungle + Relief-F" << endl;
		} else {
			if (ecAlgParam == "RJ") {
				algorithmType = EC_RJ;
				cout << Timestamp() << "Running EC in Random Jungle only mode" << endl;
			} else {
				if (ecAlgParam == "RF") {
					algorithmType = EC_RF;
					cout << Timestamp() << "Running EC in Relief-F only mode" << endl;
				} else {
					cerr << "ERROR: --ec-algorithm-steps must be one of: "
							<< "all, rj or rf" << endl;
					exit(EXIT_FAILURE);
				}
			}
		}
	}

	outFilesPrefix = paramsMap["out-files-prefix"].as<string>();

	// set the number of attributea to remove per iteration
	numToRemovePerIteration = 0;
	if (paramsMap.count("ec-iter-remove-n")) {
		numToRemovePerIteration = paramsMap["ec-iter-remove-n"].as<unsigned int>();
	}
	if (paramsMap.count("ec-iter-remove-percent")) {
		unsigned int iterPercentToRemove = paramsMap["ec-iter-remove-percent"].as<
				unsigned int>();
		numToRemovePerIteration = (unsigned int) (((double) iterPercentToRemove
				/ 100.0) * dataset->NumAttributes());
	}
	cout << Timestamp() << "EC will remove " << numToRemovePerIteration
			<< " attributes on first iteration" << endl;

	// multithreading setup
	unsigned int maxThreads = omp_get_num_procs();
	cout << Timestamp() << maxThreads << " OpenMP processors available to EC"
			<< endl;
	numRFThreads = maxThreads;
	cout << Timestamp() << "EC will use " << numRFThreads << " threads" << endl;

	// ------------------------------------------------------------- Random Jungle
	// initialize Random Jungle
	if ((algorithmType == EC_ALL) || (algorithmType == EC_RJ)) {
		randomJungle = new RandomJungle(dataset, paramsMap);
	}

	// ------------------------------------------------------------------ Relief-F
	// intialize Relief-F
	if ((algorithmType == EC_ALL) || (algorithmType == EC_RF)) {
		cout << Timestamp() << "Initializing Relief-F" << endl;
		if (dataset->HasContinuousPhenotypes()) {
			cout << Timestamp() << "RRelief-F" << endl;
			reliefF = new RReliefF(dataset, paramsMap);
		} else {
			cout << Timestamp() << "Relief-F" << endl;
			reliefF = new ReliefF(dataset, paramsMap, analysisType);
		}
	}

} // end of constructor

EvaporativeCooling::EvaporativeCooling(Dataset* ds, ConfigMap& configMap,
		AnalysisType anaType) {
	cout << Timestamp() << "Evaporative Cooling initialization:" << endl;
	if (ds) {
		dataset = ds;
	} else {
		cerr << "ERROR: data set is not initialized" << endl;
		exit(EXIT_FAILURE);
	}
	analysisType = anaType;

	reliefF = NULL;
	randomJungle = NULL;

	// set the number of target attributes
	string configValue;
	if(GetConfigValue(configMap, "ec-num-target", configValue)) {
		numTargetAttributes = lexical_cast<unsigned int>(configValue);
	}

	if (numTargetAttributes == 0) {
		numTargetAttributes = ds->NumVariables();
		numToRemovePerIteration = 0;
	}
	if (numTargetAttributes > dataset->NumVariables()) {
		cerr << "--ec-num-taget must be less than or equal to the "
				<< "number of attributes in the data set" << endl;
		exit(EXIT_FAILURE);
	}
	cout << Timestamp() << "EC is removing attributes until best "
			<< numTargetAttributes << " remain" << endl;

	string ecAlgParam = "ALL";
	algorithmType = EC_ALL;
	if (GetConfigValue(configMap, "ec-algorithm-steps", configValue)) {
		ecAlgParam = to_upper(configValue);
		if (ecAlgParam == "ALL") {
			algorithmType = EC_ALL;
			cout << Timestamp()
					<< "Running EC in standard mode: Random Jungle + Relief-F" << endl;
		} else {
			if (ecAlgParam == "RJ") {
				algorithmType = EC_RJ;
				cout << Timestamp() << "Running EC in Random Jungle only mode" << endl;
			} else {
				if (ecAlgParam == "RF") {
					algorithmType = EC_RF;
					cout << Timestamp() << "Running EC in Relief-F only mode" << endl;
				} else {
					cerr << "ERROR: --ec-algorithm-steps must be one of: "
							<< "all, rj or rf" << endl;
					exit(EXIT_FAILURE);
				}
			}
		}
	}

	if (GetConfigValue(configMap, "out-files-prefix", configValue)) {
		outFilesPrefix = configValue;
	}
	else {
		outFilesPrefix = "ec_run";
	}

	// set the number of attributea to remove per iteration
	numToRemovePerIteration = 0;
	if (GetConfigValue(configMap, "ec-iter-remove-n", configValue)) {
		numToRemovePerIteration = lexical_cast<unsigned int>(configValue);
	}
	if (GetConfigValue(configMap, "ec-iter-remove-percent", configValue)) {
		unsigned int iterPercentToRemove = lexical_cast<unsigned int>(configValue);
		numToRemovePerIteration = (unsigned int) (((double) iterPercentToRemove
				/ 100.0) * dataset->NumAttributes());
	}
	cout << Timestamp() << "EC will remove " << numToRemovePerIteration
			<< " attributes on first iteration" << endl;

	// multithreading setup
	unsigned int maxThreads = omp_get_num_procs();
	cout << Timestamp() << maxThreads << " OpenMP processors available to EC"
			<< endl;
	numRFThreads = maxThreads;
	cout << Timestamp() << "EC will use " << numRFThreads << " threads" << endl;

	// ------------------------------------------------------------- Random Jungle
	// initialize Random Jungle
	if ((algorithmType == EC_ALL) || (algorithmType == EC_RJ)) {
		randomJungle = new RandomJungle(dataset, configMap);
	}

	// ------------------------------------------------------------------ Relief-F
	// intialize Relief-F
	if ((algorithmType == EC_ALL) || (algorithmType == EC_RF)) {
		cout << Timestamp() << "Initializing Relief-F" << endl;
		if (dataset->HasContinuousPhenotypes()) {
			cout << Timestamp() << "RRelief-F" << endl;
			reliefF = new RReliefF(dataset, configMap);
		} else {
			cout << Timestamp() << "Relief-F" << endl;
			reliefF = new ReliefF(dataset, configMap, analysisType);
		}
	}

} // end of constructor

EvaporativeCooling::~EvaporativeCooling() {
	if (reliefF) {
		delete reliefF;
	}
	if (randomJungle) {
		delete randomJungle;
	}
}

bool EvaporativeCooling::ComputeECScores() {
	unsigned int numWorkingAttributes = dataset->NumVariables();
	if (numWorkingAttributes < numTargetAttributes) {
		cerr << "ERROR: The number of attributes in the data set "
				<< numWorkingAttributes
				<< " is less than the number of target attributes "
				<< numTargetAttributes << endl;
		return false;
	}

	// EC algorithm as in Figure 5, page 10 of the paper referenced
	// at top of this file. Modified per Brett's email to not do the
	// varying temperature and classifier accuracy optimization steps.
	unsigned int iteration = 1;
	clock_t t;
	float elapsedTime = 0.0;
	while (numWorkingAttributes >= numTargetAttributes) {
		cout << Timestamp()
				<< "----------------------------------------------------"
				<< "-------------------------" << endl;
		cout << Timestamp() << "EC algorithm...iteration: " << iteration
				<< ", working attributes: " << numWorkingAttributes
				<< ", target attributes: " << numTargetAttributes << endl;
		cout << fixed << setprecision(1);

		// -------------------------------------------------------------------------
		// run Random Jungle and get the normalized scores for use in EC
		if ((algorithmType == EC_ALL) || (algorithmType == EC_RJ)) {
			t = clock();
			cout << Timestamp() << "Running Random Jungle" << endl;
			if (randomJungle->ComputeAttributeScores()) {
				rjScores = randomJungle->GetScores();
			} else {
				cerr << "ERROR: In EC algorithm: Random Jungle failed" << endl;
				return false;
			}
			elapsedTime = (float) (clock() - t) / CLOCKS_PER_SEC;
			cout << Timestamp() << "Random Jungle finished in " << elapsedTime
					<< " secs" << endl;
			if ((algorithmType == EC_RJ)
					&& (numWorkingAttributes == numTargetAttributes)) {
				sort(rjScores.begin(), rjScores.end(), scoresSortDesc);
				ecScores.resize(numTargetAttributes);
				copy(rjScores.begin(), rjScores.begin() + numTargetAttributes,
						ecScores.begin());
				return true;
			}
		}

		// -------------------------------------------------------------------------
		// run Relief-F and get normalized score for use in EC
		if ((algorithmType == EC_ALL) || (algorithmType == EC_RF)) {
			t = clock();
			cout << Timestamp() << "Running ReliefF" << endl;
			if (!RunReliefF()) {
				cerr << "ERROR: In EC algorithm: ReliefF failed" << endl;
				return false;
			}
			cout << setprecision(1);
			elapsedTime = (float) (clock() - t) / CLOCKS_PER_SEC;
			cout << Timestamp() << "ReliefF finished in " << elapsedTime << " secs"
					<< endl;
			if ((algorithmType == EC_RF)
					&& (numWorkingAttributes == numTargetAttributes)) {
				sort(rfScores.begin(), rfScores.end(), scoresSortDesc);
				ecScores.resize(numTargetAttributes);
				copy(rfScores.begin(), rfScores.begin() + numTargetAttributes,
						ecScores.begin());
				return true;
			}
		}

		// -------------------------------------------------------------------------
		// compute free energy for all attributes
		t = clock();
		cout << Timestamp() << "Computing free energy" << endl;
		double temperature = 1.0;
		if (!ComputeFreeEnergy(temperature)) {
			cerr << "ERROR: In EC algorithm: ComputeFreeEnergy failed" << endl;
			return false;
		}
		elapsedTime = (float) (clock() - t) / CLOCKS_PER_SEC;
		cout << Timestamp() << "Free energy calculations complete in "
				<< elapsedTime << " secs" << endl;
		// PrintAllScoresTabular();
		// PrintKendallTaus();

		// -------------------------------------------------------------------------
		// remove the worst attributes and iterate
		t = clock();
		cout << Timestamp() << "Removing the worst attributes" << endl;
		unsigned int numToRemove = numToRemovePerIteration;
		if (paramsMap.count("ec-iter-remove-percent")) {
			unsigned int iterPercentToRemove = paramsMap["ec-iter-remove-percent"].as<
					unsigned int>();
			numToRemove = (int) (((double) iterPercentToRemove / 100.0)
					* dataset->NumVariables());
		}
		if ((numWorkingAttributes - numToRemove) < numTargetAttributes) {
			numToRemove = numWorkingAttributes - numTargetAttributes;
		}
		if (numToRemove < 1) {
//      cerr << "ERROR: Number of attributes to remove is less than one." << endl;
//      return false;
			break;
		}
		cout << Timestamp() << "Removing the worst " << numToRemove << " attributes"
				<< endl;
		if (!RemoveWorstAttributes(numToRemove)) {
			cerr << "ERROR: In EC algorithm: RemoveWorstAttribute failed" << endl;
			return false;
		}
		numWorkingAttributes -= numToRemove;
		elapsedTime = (float) (clock() - t) / CLOCKS_PER_SEC;
		cout << Timestamp() << "Attribute removal complete in " << elapsedTime
				<< " secs" << endl;

		++iteration;
	}

	cout << Timestamp() << "EC algorithm ran for " << iteration << " iterations"
			<< endl;

	// remaining free energy attributes are the ones we want to write as a
	// new dataset to be analyzed with (re)GAIN + SNPrank
	sort(freeEnergyScores.begin(), freeEnergyScores.end(), scoresSortDesc);
	ecScores.resize(numTargetAttributes);
	copy(freeEnergyScores.begin(), freeEnergyScores.begin() + numTargetAttributes,
			ecScores.begin());

	return true;
}

EcScores& EvaporativeCooling::GetRandomJungleScores() {
	return rjScores;
}

EcScores& EvaporativeCooling::GetReliefFScores() {
	return rfScores;
}

EcScores& EvaporativeCooling::GetECScores() {
	return ecScores;
}

EcAlgorithmType EvaporativeCooling::GetAlgorithmType() {
	return algorithmType;
}

void EvaporativeCooling::PrintAttributeScores(ofstream& outStream) {
	for (EcScoresCIt ecScoresIt = ecScores.begin(); ecScoresIt != ecScores.end();
			++ecScoresIt) {
		outStream << fixed << setprecision(8) << (*ecScoresIt).first << "\t"
				<< (*ecScoresIt).second << endl;
	}
}

void EvaporativeCooling::PrintRJAttributeScores(ofstream& outStream) {
	sort(rjScores.begin(), rjScores.end(), scoresSortDesc);
	for (EcScoresCIt rjScoresIt = rjScores.begin(); rjScoresIt != rjScores.end();
			++rjScoresIt) {
		outStream << fixed << setprecision(8) << (*rjScoresIt).first << "\t"
				<< (*rjScoresIt).second << endl;
	}
}

void EvaporativeCooling::PrintRFAttributeScores(ofstream& outStream) {
	sort(rfScores.begin(), rfScores.end(), scoresSortDesc);
	for (EcScoresCIt rfScoresIt = rfScores.begin(); rfScoresIt != rfScores.end();
			++rfScoresIt) {
		outStream << fixed << setprecision(8) << (*rfScoresIt).first << "\t"
				<< (*rfScoresIt).second << endl;
	}
}

void EvaporativeCooling::WriteAttributeScores(string baseFilename) {
	string resultsFilename = baseFilename;
	ofstream outFile;
	// added 9/26/11 for reflecting the fact that only parts of the
	// complete EC algorithm were performed
	switch (algorithmType) {
	case EC_ALL:
		resultsFilename = baseFilename + ".ec";
		outFile.open(resultsFilename.c_str());
		if (outFile.bad()) {
			cerr << "ERROR: Could not open scores file " << resultsFilename
					<< "for writing" << endl;
			exit(1);
		}
		PrintAttributeScores(outFile);
		outFile.close();

		resultsFilename = baseFilename + ".ec.rj";
		outFile.open(resultsFilename.c_str());
		if (outFile.bad()) {
			cerr << "ERROR: Could not open scores file " << resultsFilename
					<< "for writing" << endl;
			exit(1);
		}
		PrintRJAttributeScores(outFile);
		outFile.close();

		resultsFilename = baseFilename + ".ec.rf";
		outFile.open(resultsFilename.c_str());
		if (outFile.bad()) {
			cerr << "ERROR: Could not open scores file " << resultsFilename
					<< "for writing" << endl;
			exit(1);
		}
		PrintRFAttributeScores(outFile);
		outFile.close();
		break;
	case EC_RJ:
		resultsFilename += ".rj";
		outFile.open(resultsFilename.c_str());
		if (outFile.bad()) {
			cerr << "ERROR: Could not open scores file " << resultsFilename
					<< "for writing" << endl;
			exit(1);
		}
		PrintAttributeScores(outFile);
		outFile.close();
		break;
	case EC_RF:
		resultsFilename += ".rf";
		outFile.open(resultsFilename.c_str());
		if (outFile.bad()) {
			cerr << "ERROR: Could not open scores file " << resultsFilename
					<< "for writing" << endl;
			exit(1);
		}
		PrintAttributeScores(outFile);
		outFile.close();
		break;
	default:
		// we should not get here by the CLI front end but it is possible to call
		// this from other programs in the future or when used as a library
		// TODO: better message
		cerr << "ERROR: Attempting to write attribute scores before the analysis "
				<< "type was determined. " << endl;
		return;
	}
}

bool EvaporativeCooling::PrintAllScoresTabular() {
	// sanity checks
	if (rjScores.size() != rfScores.size()) {
		cerr
				<< "ERROR: Random Jungle and Relief-F scores lists are not the same size"
				<< endl;
		return false;
	}
	if (freeEnergyScores.size() != rfScores.size()) {
		cerr
				<< "ERROR: Random Jungle and Relief-F scores lists are not the same size"
				<< endl;
		return false;
	}

	sort(rjScores.begin(), rjScores.end(), scoresSortDesc);
	sort(rfScores.begin(), rfScores.end(), scoresSortDesc);
	sort(freeEnergyScores.begin(), freeEnergyScores.end(), scoresSortDesc);

	cout << "\t\t\tE (RF)\t\tS (RJ)\t\tF (free energy)\n";
	unsigned int numScores = freeEnergyScores.size();
	for (unsigned int i = 0; i < numScores; ++i) {
		pair<double, string> thisRJScores = rjScores[i];
		pair<double, string> thisRFScores = rfScores[i];
		pair<double, string> thisFEScores = freeEnergyScores[i];
		printf("\t\t\t%s\t%6.4f\t%s\t%6.4f\t%s\t%6.4f\n",
				thisRFScores.second.c_str(), thisRFScores.first,
				thisRJScores.second.c_str(), thisRJScores.first,
				thisFEScores.second.c_str(), thisFEScores.first);
	}

	return true;
}

bool EvaporativeCooling::PrintKendallTaus() {
	// sanity checks
	if (rjScores.size() != rfScores.size()) {
		cerr
				<< "ERROR: Random Jungle and Relief-F scores lists are not the same size"
				<< endl;
		return false;
	}
	if (freeEnergyScores.size() != rfScores.size()) {
		cerr
				<< "ERROR: Random Jungle and Relief-F scores lists are not the same size"
				<< endl;
		return false;
	}

	sort(rjScores.begin(), rjScores.end(), scoresSortDesc);
	sort(rfScores.begin(), rfScores.end(), scoresSortDesc);
	sort(freeEnergyScores.begin(), freeEnergyScores.end(), scoresSortDesc);

	vector<string> rjNames;
	vector<string> rfNames;
	vector<string> feNames;
	unsigned int numScores = freeEnergyScores.size();
	for (unsigned int i = 0; i < numScores; ++i) {
		pair<double, string> thisRJScores = rjScores[i];
		pair<double, string> thisRFScores = rfScores[i];
		pair<double, string> thisFEScores = freeEnergyScores[i];
		rjNames.push_back(thisRJScores.second);
		rfNames.push_back(thisRFScores.second);
		feNames.push_back(thisFEScores.second);
	}

	double tauRJRF = KendallTau(rjNames, rfNames);
	double tauRJFE = KendallTau(rjNames, feNames);
	double tauRFFE = KendallTau(rfNames, feNames);

	cout << "\t\t\tKendall tau's: " << "RJvRF: " << tauRJRF << ", RJvFE: "
			<< tauRJFE << ", RFvFE: " << tauRFFE << endl;

	return true;
}

bool EvaporativeCooling::RunReliefF() {

	if(reliefF->ComputeAttributeScores()) {
		rfScores = reliefF->GetScores();
		if(rfScores.size() == 0) {
			cerr << "ERROR: RunReliefF: No scores computed" << endl;
			return false;
		}
	}
	else {
		cerr << "ERROR: RunReliefF: ComputeAttributeScores failed" << endl;
		return false;
	}

	cout << Timestamp() << "Normalizing ReliefF scores to 0-1" << endl;
	pair<double, string> firstScore = rfScores[0];
	double minRFScore = firstScore.first;
	double maxRFScore = firstScore.first;
	EcScoresCIt rfScoresIt = rfScores.begin();
	for (; rfScoresIt != rfScores.end(); ++rfScoresIt) {
		pair<double, string> thisScore = *rfScoresIt;
		if (thisScore.first < minRFScore) {
			minRFScore = thisScore.first;
		}
		if (thisScore.first > maxRFScore) {
			maxRFScore = thisScore.first;
		}
	}

	// normalize attribute scores
	if (minRFScore == maxRFScore) {
		cout << Timestamp() << "WARNING: Relief-F min and max scores are the same. "
				<< "No normalization necessary" << endl;
		return true;
	}

	EcScores newRFScores;
	double rfRange = maxRFScore - minRFScore;
	for (EcScoresIt it = rfScores.begin(); it != rfScores.end(); ++it) {
		pair<double, string> thisScore = *it;
		double key = thisScore.first;
		string val = thisScore.second;
		newRFScores.push_back(make_pair((key - minRFScore) / rfRange, val));
	}

	rfScores.clear();
	rfScores = newRFScores;

	return true;
}

bool EvaporativeCooling::ComputeFreeEnergy(double temperature) {
	if (algorithmType == EC_ALL) {
		if (rjScores.size() != rfScores.size()) {
			cerr << "ERROR: EvaporativeCooling::ComputeFreeEnergy scores lists are "
					"unequal. RJ: " << rjScores.size() << " vs. RF: " << rfScores.size()
					<< endl;
			return false;
		}
	}

	freeEnergyScores.clear();
	EcScoresCIt rjIt = rjScores.begin();
	EcScoresCIt rfIt = rfScores.begin();
	switch (algorithmType) {
	case EC_ALL:
		sort(rjScores.begin(), rjScores.end(), scoresSortAscByName);
		sort(rfScores.begin(), rfScores.end(), scoresSortAscByName);
		for (; rjIt != rjScores.end(); ++rjIt, ++rfIt) {
			string val = rjIt->second;
			double key = rjIt->first;
			freeEnergyScores.push_back(
					make_pair((*rfIt).first + (temperature * key), val));
		}
		break;
	case EC_RJ:
		for (; rjIt != rjScores.end(); ++rjIt) {
			freeEnergyScores.push_back(make_pair(rjIt->first, rjIt->second));
		}
		break;
	case EC_RF:
		for (; rfIt != rfScores.end(); ++rfIt) {
			freeEnergyScores.push_back(make_pair(rfIt->first, rfIt->second));
		}
		break;
	default:
		cerr << "ERROR: EvaporativeCooling::ComputeFreeEnergy: "
				<< "could not determine EC algorithm type" << endl;
		return false;
	}

	return true;
}

bool EvaporativeCooling::RemoveWorstAttributes(unsigned int numToRemove) {
	unsigned int numToRemoveAdj = numToRemove;
	unsigned int numAttr = dataset->NumAttributes();
	if ((numAttr - numToRemove) < numTargetAttributes) {
		cout << Timestamp() << "WARNING: attempt to remove " << numToRemove
				<< " attributes which will remove more than target "
				<< "number of attributes " << numTargetAttributes << ". Adjusting"
				<< endl;
		numToRemoveAdj = numAttr - numTargetAttributes;
	}
	cout << Timestamp() << "Removing " << numToRemoveAdj << " attributes" << endl;
	sort(freeEnergyScores.begin(), freeEnergyScores.end(), scoresSortAsc);
	for (unsigned int i = 0; i < numToRemoveAdj; ++i) {

		// worst score and attribute name
		pair<double, string> worst = freeEnergyScores[i];
//    cout << "\t\t\t\tRemoving: "
//            << worst.second << " (" << worst.first << ")" << endl;

		// save worst
		evaporatedAttributes.push_back(worst);
		// remove the attribute from those under consideration
		if (!dataset->MaskRemoveVariable(worst.second)) {
			cerr << "ERROR: Could not remove worst attribute: " << worst.second
					<< endl;
			return false;
		}
	}

	return true;
}
