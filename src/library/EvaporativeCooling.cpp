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
  *
 * Made even more generic with main effects and interaction effects algorithms
 * in a class hierarchy from a AttributeRanker base. 8/12/12
 */

#include <cstdlib>
#include <iostream>
#include <iomanip>

#include <omp.h>

#include <gsl/gsl_rng.h>

#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/progress.hpp>

// class interface
#include "EvaporativeCooling.h"

// EC project
#include "Insilico.h"
#include "Dataset.h"
#include "Statistics.h"
#include "StringUtils.h"
#include "RandomJungle.h"
#include "Deseq.h"
#include "ReliefF.h"
#include "RReliefF.h"
#include "ReliefFSeq.h"

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

	optimizeTemperature = false;
	if(paramsMap.count("optimize-temp")) {
		optimizeTemperature = true;
	}
	optimalTemperature = 1.0;
	bestClassificationError = 1.0;

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

	/// set the EC steps to perform
	algorithmType = EC_ALG_ME_IT;
	if (paramsMap.count("ec-algorithm-steps")) {
		string ecAlgParam = to_upper(vm["ec-algorithm-steps"].as<string>());
		if (ecAlgParam == "ALL") {
			algorithmType = EC_ALG_ME_IT;
			cout << Timestamp()
					<< "Running EC in standard mode: Main effects + interaction effects"
					<< endl;
		} else {
			if (ecAlgParam == "ME") {
				algorithmType = EC_ALG_ME_ONLY;
				cout << Timestamp() << "Running EC in main effects only mode" << endl;
			} else {
				if (ecAlgParam == "IT") {
					algorithmType = EC_ALG_IT_ONLY;
					cout << Timestamp() << "Running EC in interactions effects only mode"
							<< endl;
				} else {
					cerr << "ERROR: --ec-algorithm-steps must be one of: "
							<< "all, me or it" << endl;
					exit(EXIT_FAILURE);
				}
			}
		}
	}

	/// set the main effects algorithm
	meAlgorithmType = EC_ME_ALG_RJ;
	if (paramsMap.count("ec-me-algorithm")) {
		string ecMeAlgParam = to_upper(vm["ec-me-algorithm"].as<string>());
		if (ecMeAlgParam == "RJ") {
			meAlgorithmType = EC_ME_ALG_RJ;
			cout << Timestamp()
					<< "EC main effects algorithm set to: Random Jungle"
					<< endl;
		} else {
			if (ecMeAlgParam == "DESEQ") {
				meAlgorithmType = EC_ME_ALG_DESEQ;
				cout << Timestamp()
						<< "Running EC in main effects algorithm set to: DESeq" << endl;
			} else {
				cerr << "ERROR: --ec-me-algorithm must be one of: "
						<< "rj or deseq" << endl;
				exit(EXIT_FAILURE);
			}
		}
	}
	maineffectAlgorithm = NULL;
	switch(meAlgorithmType) {
	case EC_ME_ALG_RJ:
		maineffectAlgorithm = new RandomJungle(ds, vm);
		break;
	case EC_ME_ALG_DESEQ:
		maineffectAlgorithm = new Deseq(ds);
		break;
	}

	/// set the interaction algorithm
	itAlgorithmType = EC_IT_ALG_RF;
	if (paramsMap.count("ec-it-algorithm")) {
		string ecItAlgParam = to_upper(vm["ec-it-algorithm"].as<string>());
		if (ecItAlgParam == "RF") {
			itAlgorithmType = EC_IT_ALG_RF;
			cout << Timestamp()
					<< "EC interaction effects algorithm set to: ReliefF"
					<< endl;
		} else {
			if (ecItAlgParam == "RFSEQ") {
				itAlgorithmType = EC_IT_ALG_RFSEQ;
				cout << Timestamp()
						<< "Running EC interaction effects algorithm set to: ReliefFSeq"
						<< endl;
			} else {
				cerr << "ERROR: --ec-it-algorithm must be one of: "
						<< "rf or rfseq" << endl;
				exit(EXIT_FAILURE);
			}
		}
	}
	interactionAlgorithm = NULL;
	switch(itAlgorithmType) {
	case EC_IT_ALG_RF:
		interactionAlgorithm = new ReliefF(ds, vm, anaType);
		break;
	case EC_IT_ALG_RFSEQ:
		interactionAlgorithm = new ReliefFSeq(ds, vm);
		break;
	}

	outFilesPrefix = paramsMap["out-files-prefix"].as<string>();

	// set the number of attributes to remove per iteration
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

	interactionAlgorithm = NULL;
	maineffectAlgorithm = NULL;

	string configValue;

	optimizeTemperature = false;
	if(GetConfigValue(configMap, "optimize-temp", configValue)) {
		optimizeTemperature = true;
	}
	optimalTemperature = 1.0;
	bestClassificationError = 1.0;

	// set the number of target attributes
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

	/// set the EC steps to perform
	algorithmType = EC_ALG_ME_IT;
	if (GetConfigValue(configMap, "ec-algorithm-steps", configValue)) {
		string ecAlgParam = to_upper(configValue);
		if (ecAlgParam == "ALL") {
			algorithmType = EC_ALG_ME_IT;
			cout << Timestamp()
					<< "Running EC in standard mode: Main effects + interaction effects"
					<< endl;
		} else {
			if (ecAlgParam == "ME") {
				algorithmType = EC_ALG_ME_ONLY;
				cout << Timestamp() << "Running EC in main effects only mode" << endl;
			} else {
				if (ecAlgParam == "IT") {
					algorithmType = EC_ALG_IT_ONLY;
					cout << Timestamp() << "Running EC in interactions effects only mode"
							<< endl;
				} else {
					cerr << "ERROR: --ec-algorithm-steps must be one of: "
							<< "all, me or it" << endl;
					exit(EXIT_FAILURE);
				}
			}
		}
	}

	/// set the main effects algorithm
	meAlgorithmType = EC_ME_ALG_RJ;
	if (GetConfigValue(configMap, "ec-me-algorithm", configValue)) {
		string ecMeAlgParam = to_upper(configValue);
		if (ecMeAlgParam == "RJ") {
			meAlgorithmType = EC_ME_ALG_RJ;
			cout << Timestamp()
					<< "EC main effects algorithm set to: Random Jungle"
					<< endl;
		} else {
			if (ecMeAlgParam == "DESEQ") {
				meAlgorithmType = EC_ME_ALG_DESEQ;
				cout << Timestamp()
						<< "Running EC in main effects algorithm set to: DESeq" << endl;
			} else {
				cerr << "ERROR: --ec-me-algorithm must be one of: "
						<< "rj or deseq" << endl;
				exit(EXIT_FAILURE);
			}
		}
	}
	maineffectAlgorithm = NULL;
	switch(meAlgorithmType) {
	case EC_ME_ALG_RJ:
		maineffectAlgorithm = new RandomJungle(ds, configMap);
		break;
	case EC_ME_ALG_DESEQ:
		maineffectAlgorithm = new Deseq(ds);
		break;
	}

	/// set the interaction algorithm
	itAlgorithmType = EC_IT_ALG_RF;
	if (GetConfigValue(configMap, "ec-it-algorithm", configValue)) {
		string ecItAlgParam = to_upper(configValue);
		if (ecItAlgParam == "RF") {
			itAlgorithmType = EC_IT_ALG_RF;
			cout << Timestamp()
					<< "EC interaction effects algorithm set to: ReliefF"
					<< endl;
		} else {
			if (ecItAlgParam == "RFSEQ") {
				itAlgorithmType = EC_IT_ALG_RFSEQ;
				cout << Timestamp()
						<< "Running EC in interaction effects algorithm set to: ReliefFSeq"
						<< endl;
			} else {
				cerr << "ERROR: --ec-it-algorithm must be one of: "
						<< "rf or rfseq" << endl;
				exit(EXIT_FAILURE);
			}
		}
	}
	interactionAlgorithm = NULL;
	switch(meAlgorithmType) {
	case EC_IT_ALG_RF:
		interactionAlgorithm = new ReliefF(ds, configMap, anaType);
		break;
	case EC_IT_ALG_RFSEQ:
		interactionAlgorithm = new ReliefFSeq(ds, configMap);
		break;
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

} // end of constructor

EvaporativeCooling::~EvaporativeCooling() {
	if (interactionAlgorithm) {
		delete interactionAlgorithm;
	}
	if (maineffectAlgorithm) {
		delete maineffectAlgorithm;
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
	// Added temperature optimization - 4/9/12
	unsigned int iteration = 1;
	boost::progress_timer t;
	float elapsedTime = 0.0;
	optimalTemperature = 1.0;
	while (numWorkingAttributes >= numTargetAttributes) {
		pair<unsigned int, unsigned int> titvCounts =
				dataset->GetAttributeTiTvCounts();
		double titvRatio = titvCounts.first;
		if(titvCounts.second) {
			titvRatio = (double) titvCounts.first / (double) titvCounts.second;
		}

		cout << Timestamp()
				<< "----------------------------------------------------"
				<< "-------------------------" << endl;
		cout << Timestamp() << "EC algorithm...iteration: " << iteration
				<< ", working attributes: " << numWorkingAttributes
				<< ", target attributes: " << numTargetAttributes
				<< ", temperature: " << optimalTemperature << endl;
		cout << Timestamp()
				<< "Ti/Tv: transitions: " << titvCounts.first
				<< " transversions: " << titvCounts.second
				<< " ratio: " << titvRatio
				<< endl;
		cout << fixed << setprecision(1);

		// -------------------------------------------------------------------------
		// run main effects algorithm and get the normalized scores for use in EC
		if ((algorithmType == EC_ALG_ME_IT) || (algorithmType == EC_ALG_ME_ONLY)) {
			cout << Timestamp() << "Running Random Jungle" << endl;
			maineffectScores = maineffectAlgorithm->ComputeScores();
			bestClassificationError = maineffectAlgorithm->GetClassificationError();
			cout << Timestamp() << "Main effects ranker finished in " << t.elapsed()
					<< " secs" << endl;
			// RJ standalone runs
			if ((algorithmType == EC_ALG_ME_ONLY)
					&& (numWorkingAttributes == numTargetAttributes)) {
				sort(maineffectScores.begin(), maineffectScores.end(), scoresSortDesc);
				ecScores.resize(numTargetAttributes);
				copy(maineffectScores.begin(), maineffectScores.begin() + numTargetAttributes,
						ecScores.begin());
				return true;
			}
		}

		// -------------------------------------------------------------------------
		// run interaction effects algorithm and get normalized score for use in EC
		if ((algorithmType == EC_ALG_ME_IT) || (algorithmType == EC_ALG_IT_ONLY)) {
			cout << Timestamp() << "Running ReliefF" << endl;
			if (!RunReliefF()) {
				cerr << "ERROR: In EC algorithm: ReliefF failed" << endl;
				return false;
			}
			cout << setprecision(1);
			cout << Timestamp() << "ReliefF finished in " << t.elapsed() << " secs"
					<< endl;
			// ReliefF standalone runs
			if ((algorithmType == EC_ALG_IT_ONLY)
					&& (numWorkingAttributes == numTargetAttributes)) {
				sort(interactionScores.begin(), interactionScores.end(), scoresSortDesc);
				ecScores.resize(numTargetAttributes);
				copy(interactionScores.begin(), interactionScores.begin() + numTargetAttributes,
						ecScores.begin());
				return true;
			}
		}

		// -------------------------------------------------------------------------
		// compute free energy for all attributes
		cout << Timestamp() << "Computing free energy" << endl;
		if (!ComputeFreeEnergy(optimalTemperature)) {
			cerr << "ERROR: In EC algorithm: ComputeFreeEnergy failed" << endl;
			return false;
		}
		cout << Timestamp() << "Free energy calculations complete in "
				<< t.elapsed() << " secs" << endl;
		// PrintAllScoresTabular();
		// PrintKendallTaus();

		// -------------------------------------------------------------------------
		// optimize the temperature by sampling a set of delta values around T
		if(optimizeTemperature) {
			cout << Timestamp() << "Optimizing coupling temperature T" << endl;
			vector<double> temperatureDeltas;
			temperatureDeltas.push_back(-0.2);
			temperatureDeltas.push_back(0.2);
			optimalTemperature = OptimizeTemperature(temperatureDeltas);
			cout << Timestamp() << "T optimization: " << optimalTemperature
					<< ", complete in " << t.elapsed() << " secs" << endl;
		}

		// -------------------------------------------------------------------------
		// remove the worst attributes and iterate
		cout << Timestamp() << "Removing the worst attributes" << endl;
		unsigned int numToRemove = numToRemovePerIteration;
		numToRemoveNextIteration = numToRemove - numToRemovePerIteration;
		if (paramsMap.count("ec-iter-remove-percent")) {
			unsigned int iterPercentToRemove = paramsMap["ec-iter-remove-percent"].as<
					unsigned int>();
			numToRemove = (int) (((double) iterPercentToRemove / 100.0)
					* dataset->NumVariables());
			numToRemoveNextIteration = (int) (((double) iterPercentToRemove / 100.0)
					* numToRemove);
			numToRemovePerIteration = numToRemove;
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
		cout << Timestamp() << "Attribute removal complete in " << t.elapsed()
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

AttributeScores& EvaporativeCooling::GetMaineffectScores() {
	return maineffectScores;
}

AttributeScores& EvaporativeCooling::GetInteractionScores() {
	return interactionScores;
}

AttributeScores& EvaporativeCooling::GetECScores() {
	return ecScores;
}

EcAlgorithmType EvaporativeCooling::GetAlgorithmType() {
	return algorithmType;
}

void EvaporativeCooling::PrintAttributeScores(ofstream& outStream) {
	for (AttributeScoresCIt ecScoresIt = ecScores.begin(); ecScoresIt != ecScores.end();
			++ecScoresIt) {
		outStream << fixed << setprecision(8) << (*ecScoresIt).first << "\t"
				<< (*ecScoresIt).second << endl;
	}
}

void EvaporativeCooling::PrintMaineffectAttributeScores(ofstream& outStream) {
	sort(maineffectScores.begin(), maineffectScores.end(), scoresSortDesc);
	for (AttributeScoresCIt rjScoresIt = maineffectScores.begin(); rjScoresIt != maineffectScores.end();
			++rjScoresIt) {
		outStream << fixed << setprecision(8) << (*rjScoresIt).first << "\t"
				<< (*rjScoresIt).second << endl;
	}
}

void EvaporativeCooling::PrintInteractionAttributeScores(ofstream& outStream) {
	sort(interactionScores.begin(), interactionScores.end(), scoresSortDesc);
	for (AttributeScoresCIt rfScoresIt = interactionScores.begin(); rfScoresIt != interactionScores.end();
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
	case EC_ALG_ME_IT:
		resultsFilename = baseFilename + ".ec";
		outFile.open(resultsFilename.c_str());
		if (outFile.bad()) {
			cerr << "ERROR: Could not open scores file " << resultsFilename
					<< "for writing" << endl;
			exit(1);
		}
		cout << Timestamp()
				<< "Writing EC scores to [" + resultsFilename + "]" << endl;
		PrintAttributeScores(outFile);
		outFile.close();

		resultsFilename = baseFilename + ".ec.me";
		outFile.open(resultsFilename.c_str());
		if (outFile.bad()) {
			cerr << "ERROR: Could not open scores file " << resultsFilename
					<< "for writing" << endl;
			exit(1);
		}
		cout << Timestamp()
				<< "Writing EC main effects scores to [" + resultsFilename + "]" << endl;
		PrintMaineffectAttributeScores(outFile);
		outFile.close();

		resultsFilename = baseFilename + ".ec.it";
		outFile.open(resultsFilename.c_str());
		if (outFile.bad()) {
			cerr << "ERROR: Could not open scores file " << resultsFilename
					<< "for writing" << endl;
			exit(1);
		}
		cout << Timestamp()
				<< "Writing EC interaction effects scores to [" + resultsFilename + "]"
				<< endl;
		PrintInteractionAttributeScores(outFile);
		outFile.close();
		break;
	case EC_ALG_ME_ONLY:
		resultsFilename += ".me";
		outFile.open(resultsFilename.c_str());
		if (outFile.bad()) {
			cerr << "ERROR: Could not open scores file " << resultsFilename
					<< "for writing" << endl;
			exit(1);
		}
		cout << Timestamp()
				<< "Writing EC main effects scores to [" + resultsFilename + "]" << endl;
		PrintAttributeScores(outFile);
		outFile.close();
		break;
	case EC_ALG_IT_ONLY:
		if(itAlgorithmType == EC_IT_ALG_RF) {
			resultsFilename += ".it";
		}
		if(itAlgorithmType == EC_IT_ALG_RFSEQ) {
			resultsFilename += ".itseq";
		}
		outFile.open(resultsFilename.c_str());
		if (outFile.bad()) {
			cerr << "ERROR: Could not open scores file " << resultsFilename
					<< "for writing" << endl;
			exit(1);
		}
		cout << Timestamp()
				<< "Writing EC interaction effects scores to [" + resultsFilename + "]"
				<< endl;
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
	if (maineffectScores.size() != interactionScores.size()) {
		cerr
				<< "ERROR: Random Jungle and Relief-F scores lists are not the same size"
				<< endl;
		return false;
	}
	if (freeEnergyScores.size() != interactionScores.size()) {
		cerr
				<< "ERROR: Random Jungle and Relief-F scores lists are not the same size"
				<< endl;
		return false;
	}

	sort(maineffectScores.begin(), maineffectScores.end(), scoresSortDesc);
	sort(interactionScores.begin(), interactionScores.end(), scoresSortDesc);
	sort(freeEnergyScores.begin(), freeEnergyScores.end(), scoresSortDesc);

	cout << "\t\t\tE (RF)\t\tS (RJ)\t\tF (free energy)\n";
	unsigned int numScores = freeEnergyScores.size();
	for (unsigned int i = 0; i < numScores; ++i) {
		pair<double, string> thisRJScores = maineffectScores[i];
		pair<double, string> thisRFScores = interactionScores[i];
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
	if (maineffectScores.size() != interactionScores.size()) {
		cerr
				<< "ERROR: Random Jungle and Relief-F scores lists are not the same size"
				<< endl;
		return false;
	}
	if (freeEnergyScores.size() != interactionScores.size()) {
		cerr
				<< "ERROR: Random Jungle and Relief-F scores lists are not the same size"
				<< endl;
		return false;
	}

	sort(maineffectScores.begin(), maineffectScores.end(), scoresSortDesc);
	sort(interactionScores.begin(), interactionScores.end(), scoresSortDesc);
	sort(freeEnergyScores.begin(), freeEnergyScores.end(), scoresSortDesc);

	vector<string> rjNames;
	vector<string> rfNames;
	vector<string> feNames;
	unsigned int numScores = freeEnergyScores.size();
	for (unsigned int i = 0; i < numScores; ++i) {
		pair<double, string> thisRJScores = maineffectScores[i];
		pair<double, string> thisRFScores = interactionScores[i];
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
	/// postcondition: interactionScores contains the newly-computed scores
	interactionScores = interactionAlgorithm->ComputeScores();
	if(interactionScores.size() == 0) {
		cerr << "ERROR: RunReliefF: No scores computed" << endl;
		return false;
	}

	// added for rnaSeq analysis - bcw - 8/21/12
	// normalizing "stretches" the distribution to many zeroes, thus
	// complicating ranking of many ties; don't do it!
	if(itAlgorithmType == EC_IT_ALG_RFSEQ) {
		cout << Timestamp() << "rnaSeq skipping normalization."	<< endl;
		return true;
	}

	cout << Timestamp() << "Normalizing ReliefF scores to 0-1" << endl;
	pair<double, string> firstScore = interactionScores[0];
	double minRFScore = firstScore.first;
	double maxRFScore = firstScore.first;
	AttributeScoresCIt rfScoresIt = interactionScores.begin();
	for (; rfScoresIt != interactionScores.end(); ++rfScoresIt) {
		pair<double, string> thisScore = *rfScoresIt;
		if (thisScore.first < minRFScore) {
			minRFScore = thisScore.first;
		}
		if (thisScore.first > maxRFScore) {
			maxRFScore = thisScore.first;
		}
	}

	// normalize attribute scores if necessary
	if (minRFScore == maxRFScore) {
		cout << Timestamp() << "WARNING: Relief-F min and max scores are the same. "
				<< "No normalization necessary" << endl;
		return true;
	}

	AttributeScores newRFScores;
	double rfRange = maxRFScore - minRFScore;
	for (AttributeScoresIt it = interactionScores.begin(); it != interactionScores.end(); ++it) {
		pair<double, string> thisScore = *it;
		double key = thisScore.first;
		string val = thisScore.second;
		newRFScores.push_back(make_pair((key - minRFScore) / rfRange, val));
	}

	interactionScores.clear();
	interactionScores = newRFScores;

	return true;
}

bool EvaporativeCooling::ComputeFreeEnergy(double temperature) {
	if (algorithmType == EC_ALG_ME_IT) {
		if (maineffectScores.size() != interactionScores.size()) {
			cerr << "ERROR: EvaporativeCooling::ComputeFreeEnergy scores lists are "
					"unequal. RJ: " << maineffectScores.size() << " vs. RF: " << interactionScores.size()
					<< endl;
			return false;
		}
	}

	freeEnergyScores.clear();
	AttributeScoresCIt rjIt = maineffectScores.begin();
	AttributeScoresCIt rfIt = interactionScores.begin();
	switch (algorithmType) {
	case EC_ALG_ME_IT:
		sort(maineffectScores.begin(), maineffectScores.end(), scoresSortAscByName);
		sort(interactionScores.begin(), interactionScores.end(), scoresSortAscByName);
		for (; rjIt != maineffectScores.end(); ++rjIt, ++rfIt) {
			string val = rjIt->second;
			double key = rjIt->first;
			freeEnergyScores.push_back(
					make_pair((*rfIt).first + (temperature * key), val));
		}
		break;
	case EC_ALG_ME_ONLY:
		for (; rjIt != maineffectScores.end(); ++rjIt) {
			freeEnergyScores.push_back(make_pair(rjIt->first, rjIt->second));
		}
		break;
	case EC_ALG_IT_ONLY:
		for (; rfIt != interactionScores.end(); ++rfIt) {
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

double EvaporativeCooling::OptimizeTemperature(vector<double> deltas) {
	cout << Timestamp() << "--- OPTIMIZER BEGIN: Classification error to beat: "
			<< bestClassificationError << endl;
	/// for each delta, run a classifier on the best attributes according
	/// to the free energy and update best temperature
	vector<double>::const_iterator deltaIt = deltas.begin();
	AttributeScores bestFreeEnergyScores = freeEnergyScores;
	for(; deltaIt != deltas.end(); ++deltaIt) {
		double thisTemp = optimalTemperature + *deltaIt;
		ComputeFreeEnergy(thisTemp);
		double thisClassificationError = ComputeClassificationErrorRJ();
		cout << Timestamp()
				<< "OPTIMIZER: Trying temperature: " << thisTemp
				<< " => Classification Error: " << thisClassificationError
				<< endl;
		/// if classification error is lower at this delta, update best temperature
		/// and best classification error
		if(thisClassificationError < bestClassificationError) {
			cout << Timestamp() << "OPTIMIZER: found better temperature: "
					<< thisTemp << endl;
			bestClassificationError = thisClassificationError;
			optimalTemperature = thisTemp;
			bestFreeEnergyScores = freeEnergyScores;
		}
	}

	cout << Timestamp() << "--- OPTIMIZER RESULTS: "
			<< "Temperature: " << optimalTemperature
			<< ", Classification error: " << bestClassificationError << endl;

	freeEnergyScores = bestFreeEnergyScores;

	return optimalTemperature;
}

double EvaporativeCooling::ComputeClassificationErrorRJ() {
	/// get the best attribute names based on free energy score
	unsigned int numBest = freeEnergyScores.size() - numToRemoveNextIteration;
	if(!numBest) {
		cerr << "ERROR: Best results calculation results in zero attributes" << endl;
		cerr << "Number of best to use in classifier: " << numBest << endl;
		cerr << "Free energy scores: " << freeEnergyScores.size() << endl;
		cerr << "Number to remove next iteration: " << numToRemoveNextIteration << endl;
		exit(EXIT_FAILURE);
	}
	cout << Timestamp() << "Getting best " << numBest
			<< " attributes for temporary CSV file" << endl;
	vector<string> bestAttributes;
	unsigned int numCopied = 0;
	AttributeScoresCIt scoreIt = freeEnergyScores.begin();
	for(; numCopied < numBest && scoreIt != freeEnergyScores.end();
			++numCopied, ++scoreIt) {
		bestAttributes.push_back(scoreIt->second);
	}
	if(!bestAttributes.size() || (bestAttributes.size() != numBest)) {
		cerr << "ERROR: could not get " << numBest << " attributes, got: "
				<< bestAttributes.size() << endl;
		exit(EXIT_FAILURE);
	}
	/// write new data set with worst attributes removed
	string newDatasetFilename = outFilesPrefix + "_CE.csv";
	bool newDatasetSuccess = dataset->WriteNewDataset(newDatasetFilename,
			bestAttributes, CSV_DELIMITED_DATASET);
	if(!newDatasetSuccess) {
		cerr << "ERROR: Could not write new data set: " << newDatasetFilename
				<< endl;
		exit(EXIT_FAILURE);
	}

	/// create a configuration map for RandomJungle constructor
	ConfigMap configMap;
	stringstream ss;
	if(paramsMap.count("rj-num-trees")) {
		unsigned int numTrees = paramsMap["rj-num-trees"].as<unsigned int>();
		ss << numTrees;
		configMap.insert(make_pair("rj-num-trees", ss.str()));
	}
	else {
		configMap.insert(make_pair("rj-num-trees", "1000"));
	}
	if(paramsMap.count("verbose")) {
		configMap.insert(make_pair("verbose", "true"));
	}
	configMap.insert(make_pair("out-files-prefix", outFilesPrefix));
	ss.str("");
	ss << omp_get_num_procs();
	configMap.insert(make_pair("num-threads", ss.str()));

	/// run Random Jungle classifier and read classification error
	double classifierError = 1.0;
	bool rjSuccess = RandomJungle::RunClassifier(newDatasetFilename,
			configMap, dataset->DetermineTreeType().first, classifierError);

	/// remove the temporary file
	cout << Timestamp() << "Removing temporary file for RJ: "
			<< newDatasetFilename << endl;
	unlink(newDatasetFilename.c_str());

	if(!rjSuccess) {
		cerr << "Error running Random Jungle classifier" << endl;
		exit(EXIT_FAILURE);
	}

	/// return the classification error on this data
	return classifierError;
}
