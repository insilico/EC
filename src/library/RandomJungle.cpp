/* 
 * RandomJungle.cpp
 * 
 * Created on October 16, 2011, 3:45 PM
 *
 * Adapter class to map EC call for Random Jungle importance scores
 * to Random Jungle library functions.
 */

// Random Jungle source distribution
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <iomanip>

#include <string.h>
#include <math.h>

#include <omp.h>

#include "rjungle/librjungle.h"
#include "rjungle/RJunglePar.h"
#include "rjungle/RJungleCtrl.h"
#include "rjungle/DataFrame.h"
#include "rjungle/FittingFct.h"
#include "rjungle/RJungleHelper.h"
#include "rjungle/Helper.h"

#include "boost/lexical_cast.hpp"

#include "RandomJungle.h"
#include "StringUtils.h"
#include "Insilico.h"
#include "AttributeRanker.h"

using namespace std;
using namespace insilico;
using namespace boost;

// static methods
bool RandomJungle::RunClassifier(string csvFile, ConfigMap& vm,
		RandomJungleTreeType treeType, double& classError) {
	string outPrefix(vm["out-files-prefix"]);
	string confusionFilename;
	if ((treeType == NUMERIC_NUMERIC_TREE)
			|| (treeType == NUMERIC_NOMINAL_TREE)) {
		confusionFilename = outPrefix + ".confusion";
	} else {
		confusionFilename = outPrefix + ".confusion2";
	}

	/// run rjungle through a system call to the shell
	stringstream rjCmd;
	rjCmd << "rjungle" << " -f " << csvFile << " -e ','" << " -D 'Class'"
			<< " -o " << outPrefix << " -U " << vm["num-threads"] << " -t "
			<< vm["rj-num-trees"] << " -B " << vm["rj-backsel"] << " -i "
			<< vm["rj-impmeasure"] << " -j " << vm["rj-nimpvar"] << " -y "
			<< treeType;
	if (vm.find("rj-mtry") != vm.end()) {
		rjCmd << " -m " << vm["rj-mtry"];
	}
	if (vm["verbose"] == "true") {
		rjCmd << " -v";
	}
	cout << Timestamp() << "Running RJ command: " << rjCmd.str() << endl;
	int systemCallReturnStatus = system(rjCmd.str().c_str());
	if (systemCallReturnStatus == -1) {
		cerr << "ERROR: Calling rjungle executable. -1 return code" << endl;
		return false;
	}

	/// loads rj classification error from confusion file
	cout << Timestamp() << "Loading RJ classification error " << "from ["
			<< confusionFilename << "]" << endl;
	if (!ReadClassificationError(confusionFilename, treeType, classError)) {
		cerr << "ERROR: Could not read Random Jungle classification error"
				<< endl;
		return false;
	}

	return true;
}

bool RandomJungle::ReadClassificationError(std::string confusionFilename,
		RandomJungleTreeType treeType, double& classifierError) {
	/// open the confusion file
	ifstream confusionStream(confusionFilename.c_str());
	if (!confusionStream.is_open()) {
		cerr << "ERROR: Could not open Random Jungle confusion file: "
				<< confusionFilename << endl;
		return false;
	}
	string line;
	/// strip the header line(s), read the error, cast to double
	getline(confusionStream, line);
	getline(confusionStream, line);
	vector<string> tokens;
	split(tokens, line, "\t");
	if (tokens.size() != 5) {
		cerr << "ERROR: RandomJungle::GetClassificationError: "
				<< "error parsing " << confusionFilename << "." << endl;
		cerr << "Read " << tokens.size()
				<< " columns from line 2, should be 5" << endl;
		return false;
	}
	classifierError = lexical_cast<double>(trim(tokens[4]));
	confusionStream.close();

	cout << Timestamp() << "Read classification error from "
			<< confusionFilename << ": " << setprecision(4) << classifierError
			<< endl;

	return true;
}

// standard methods
RandomJungle::RandomJungle(Dataset* ds, po::variables_map& vm):
	AttributeRanker::AttributeRanker(ds) {
	dataset = ds;

	cout << Timestamp()
			<< "Initializing Random Jungle with Boost program options" << endl;

	if (vm.count("rj-run-mode")) {
		runMode = (RandomJungleRunMode) vm["rj-run-mode"].as<unsigned int>();
	} else {
		runMode = LIBRARY_RUN_MODE;
	}

	/// set rjParams
	rjParams = initRJunglePar();
	rjParams.delimiter = '\t';

	if (vm.count("rj-num-trees")) {
		rjParams.ntree = vm["rj-num-trees"].as<uli_t>();
	} else {
		cerr << "RandomJungle constructor: Unexpected condition."
				<< "rj-num-trees should have a default" << endl;
		exit(1);
	}

	// added new RJ params 5/24/12 per Scott Dudek request/experience
	if (vm.count("rj-mtry")) {
		rjParams.mtry = vm["rj-mtry"].as<uli_t>();
		fixedMtry = true;
	} else {
		rjParams.mtry = (uli_t) sqrt((double) ds->NumVariables());
		fixedMtry = false;
	}
	if (vm.count("rj-impmeasure")) {
		rjParams.impMeasure = vm["rj-impmeasure"].as<unsigned int>();
	}
	if (vm.count("rj-nimpvar")) {
		rjParams.numOfImpVar = vm["rj-nimpvar"].as<uli_t>();
	}
	if (vm.count("rj-backsel")) {
		rjParams.backSel = vm["rj-backsel"].as<unsigned int>();
	}

	if (vm.count("rj-tree-type")) {
		rjParams.treeType = vm["rj-tree-type"].as<unsigned int>();
	} else {
		rjParams.treeType = NOMINAL_NUMERIC_TREE;
	}

	if (vm.count("rj-memory-mode")) {
		rjParams.memMode = vm["rj-memory-mode"].as<unsigned int>();
	} else {
		rjParams.memMode = 0;
	}

	// added user-specified random number generator seed - 7/8/12
	rjParams.rng = gsl_rng_alloc(gsl_rng_mt19937);
	rjParams.seed = vm["rj-rng-seed"].as<unsigned int>();
	gsl_rng_set(rjParams.rng, rjParams.seed);

	rjParams.nrow = dataset->NumInstances();
	rjParams.depVarName = (char *) "Class";
	rjParams.filename = (char*) "";

	rjParams.mpiId = 0;

	rjParams.verbose_flag = vm.count("verbose") ? true : false;

	string outFilesPrefix = vm["out-files-prefix"].as<string>();
	rjParams.outprefix = strdup(outFilesPrefix.c_str());

	int numProcs = omp_get_num_procs();
	cout << Timestamp() << "Using all " << numProcs
			<< " OpenMP processors available" << endl;
	rjParams.nthreads = numProcs;
}

RandomJungle::RandomJungle(Dataset* ds, ConfigMap& configMap):
			AttributeRanker::AttributeRanker(ds)  {

	cout << Timestamp() << "Initializing Random Jungle with ConfigMap" << endl;

	dataset = ds;

	rjParams = initRJunglePar();
	rjParams.delimiter = '\t';

	string configValue;

	// run mode: library = 1, or system call = 2
	if (GetConfigValue(configMap, "rj-run-mode", configValue)) {
		runMode = (RandomJungleRunMode) lexical_cast<unsigned int>(configValue);
	} else {
		runMode = LIBRARY_RUN_MODE;
	}

	if (GetConfigValue(configMap, "rj-num-trees", configValue)) {
		unsigned int numTrees = lexical_cast<unsigned int>(configValue);
		rjParams.ntree = numTrees;
	} else {
		cout << Timestamp() << "Setting RandomJungle number of trees to 500" << endl;
		rjParams.ntree = 500;
	}

	// added new RJ params 5/24/12 per Scott Dudek request/experience
	if (GetConfigValue(configMap, "rj-mtry", configValue)) {
		rjParams.mtry = lexical_cast<uli_t>(configValue);
		fixedMtry = true;
	} else {
		rjParams.mtry = (uli_t) sqrt((double) ds->NumVariables());
		fixedMtry = false;
	}
	if (GetConfigValue(configMap, "rj-nimpvar", configValue)) {
		rjParams.numOfImpVar = lexical_cast<uli_t>(configValue);
	}
	if (GetConfigValue(configMap, "rj-impmeasure", configValue)) {
		rjParams.impMeasure = lexical_cast<unsigned int>(configValue);
	}
	if (GetConfigValue(configMap, "rj-backsel", configValue)) {
		rjParams.backSel = lexical_cast<unsigned int>(configValue);
	}

	// added user-specified random number generator seed - 7/8/12
	rjParams.rng = gsl_rng_alloc(gsl_rng_mt19937);
	if (GetConfigValue(configMap, "rj-rng-seed", configValue)) {
		rjParams.seed = lexical_cast<unsigned int>(configValue);
	}
	gsl_rng_set(rjParams.rng, rjParams.seed);

	rjParams.nrow = dataset->NumInstances();
	strcpy(rjParams.depVarName, "Class");
	strcpy(rjParams.filename, "");

	rjParams.mpiId = 0;

	rjParams.verbose_flag =
			GetConfigValue(configMap, "verbose", configValue) ? true : false;

	string outFilesPrefix = "rj_run";
	if (GetConfigValue(configMap, "out-files-prefix", configValue)) {
		outFilesPrefix = configValue;
	}
	rjParams.outprefix = strdup(outFilesPrefix.c_str());

	int numProcs = omp_get_num_procs();
	int numThreads = omp_get_num_threads();
	cout << Timestamp() << numProcs << " OpenMP processors available" << endl;
	cout << Timestamp() << numThreads << " OpenMP threads running" << endl;
	rjParams.nthreads = numProcs;
}

RandomJungle::~RandomJungle() {
	cout << Timestamp() << "Removing temporary RandomJungle files" << endl;
	vector<string> tempFilenames;
	string outprefix(rjParams.outprefix);
	tempFilenames.push_back(outprefix + ".log");
	tempFilenames.push_back(outprefix + ".verbose");
	tempFilenames.push_back(outprefix + ".importance");
	tempFilenames.push_back(outprefix + ".confusion");
	tempFilenames.push_back(outprefix + ".confusion2");
	for(vector<string>::const_iterator it=tempFilenames.begin();
			it != tempFilenames.end(); ++it) {
		unlink((*it).c_str());
	}
	if (rjParams.rng) {
		gsl_rng_free(rjParams.rng);
	}
}

AttributeScores RandomJungle::ComputeScores() {

	cout << Timestamp() << "Computing Random Jungle variable importance scores"
			<< endl;

	if (runMode == SYSTEM_CALL_RUN_MODE) {
		ComputeAttributeScoresRjungle();
		return scores;
	}

	if (runMode == LIBRARY_FILE_RUN_MODE) {
		ComputeAttributeScoresFileIO();
		return scores;
	}

	cout << Timestamp() << "Running Random Jungle using C++ librjungle calls"
			<< endl;

	vector<uli_t>* colMaskVec = NULL;
	time_t start, end;
	clock_t startgrow, endgrow;

	/// set the default mtry
	if (!fixedMtry) {
		rjParams.mtry = (uli_t) sqrt((double) dataset->NumVariables());
	}

	// this needs to be done for every iteration of Random Jungle
	// so cannot be set once in the constructor
	rjParams.ncol = dataset->NumVariables() + 1;
	rjParams.depVar = rjParams.ncol - 1;

	string outPrefix(rjParams.outprefix);
	string importanceFilename = outPrefix + ".importance";
	string confusionFilename = outPrefix + ".confusion";

	unsigned int numInstances = dataset->NumInstances();
	vector<string> instanceIds = dataset->GetInstanceIds();
	//cout << "Variables names from the data set:" << endl;
	//copy(variableNames.begin(), variableNames.end(),
	//ostream_iterator<string>(cout, "\n"));
	vector<string> attributeNames = dataset->GetFileAttributeNames();
	vector<string> numericNames = dataset->GetFileNumericsNames();
	vector<string> variableNames(attributeNames.size() + numericNames.size() + 1);
	copy(attributeNames.begin(), attributeNames.end(), variableNames.begin());
	copy(numericNames.begin(), numericNames.end(),
			variableNames.begin() + attributeNames.size());
	variableNames[attributeNames.size() + numericNames.size()] = "Class";

	cout << Timestamp() << "Preparing Random Jungle type " << rjParams.treeType
			<< endl;

	// create controller
	RJungleCtrl<NumericLevel> rjCtrl;

	RJungleIO io;
	io.open(rjParams);

	time(&start);

	RJungleHelper<NumericLevel>::printHeader(rjParams, io, start);

	DataFrame<NumericLevel>* data = new DataFrame<NumericLevel>(rjParams);
	RJungleGen<NumericLevel> rjGen;
	data->setDim(rjParams.nrow, rjParams.ncol);
	data->initMatrix();
	data->setVarNames(variableNames);

	// load data frame
	// TODO: do not load data frame every time-- use column mask mechanism?
	cout << Timestamp() << "Loading RJ DataFrame with double values, "
			<< rjParams.nrow << " rows and " << rjParams.ncol << " columns."
			<< endl << Timestamp();
	cout.flush();
	for (unsigned int i = 0; i < numInstances; ++i) {
		unsigned int instanceIndex;
		dataset->GetInstanceIndexForID(instanceIds[i], instanceIndex);
		unsigned int j = 0;
		for (unsigned int aIdx = 0; aIdx < attributeNames.size(); aIdx++) {
			unsigned int attrIdx = dataset->GetAttributeIndexFromName(
					attributeNames[aIdx]);
			AttributeLevel A =
					dataset->GetInstance(instanceIndex)->attributes[attrIdx];
			data->set(i, j, static_cast<double>(A));
			++j;
		}
		for (unsigned int nIdx = 0; nIdx < numericNames.size(); nIdx++) {
			unsigned int numIdx = dataset->GetNumericIndexFromName(
					numericNames[nIdx]);
			NumericLevel N =
					dataset->GetInstance(instanceIndex)->numerics[numIdx];
			data->set(i, j, N);
			++j;
		}
		NumericLevel pheno = 0.0;
		if (dataset->HasContinuousPhenotypes()) {
			pheno = dataset->GetInstance(instanceIndex)->GetPredictedValueTau();
			data->set(i, j, pheno);
		} else {
			pheno =
				static_cast<NumericLevel>(dataset->GetInstance(instanceIndex)->GetClass());
//				cout << "Phenotype for instance " << i
//						<< ", ID " << instanceIds[i]
//						<< ", instanceIndex " << instanceIndex
//						<< ", class index " << j
//						<< ", class " << pheno << endl;
			data->set(i, j, pheno);
		}
		// happy lights
		if (i && ((i % 100) == 0)) {
			cout << i << "/" << numInstances << " ";
			cout.flush();
		}
		// happy lights
		if (i && ((i % 1000) == 0)) {
			cout << endl << Timestamp();
		}
	}
	cout << numInstances << "/" << numInstances << endl;

	data->setDepVarName(string(rjParams.depVarName));
	data->storeCategories();
	data->makeDepVecs();
	data->getMissings();

#ifdef __DEBUG__
	ofstream ofs("data_matrix_ec.txt");
	data->print(ofs);
	ofs.close();
#endif

	rjGen.init(rjParams, *data);

	RJungleHelper<NumericLevel>::printRJunglePar(rjParams, *io.outLog);

	startgrow = clock();
	cout << Timestamp() << "Running Random Jungle" << endl;
	rjCtrl.autoBuildInternal(rjParams, io, rjGen, *data, colMaskVec);
	classificationAccuracy = rjCtrl.getOobPredAcc();
	endgrow = clock();

	time(&end);

	// print info stuff
	RJungleHelper<NumericLevel>::printFooter(rjParams, io, start, end,
			startgrow, endgrow);
	delete data;

	// clean up
	if (colMaskVec != NULL) {
		delete colMaskVec;
	}

	// clean up Random Jungle run
	io.close();

	/// loads rjScores map
	cout << Timestamp() << "Loading RJ variable importance (VI) scores "
			<< "from [" << importanceFilename << "]" << endl;
	if (!ReadScores(importanceFilename)) {
		cerr << "ERROR: Could not read Random Jungle scores" << endl;
 		exit(-1);
	}

	/// rj classification error
	// from confusion file - 4/11/12
	// added the new getOob - 7/1/12
	cout << Timestamp() << "RJ classification accuracy: "
			<< classificationAccuracy << endl;

	return scores;
}

bool RandomJungle::ComputeAttributeScoresFileIO() {
	cout << Timestamp() << "Running rjungle through C++ library with file I/O"
			<< endl;
	vector<uli_t>* colMaskVec = NULL;

	string outPrefix(rjParams.outprefix);
	string importanceFilename = outPrefix + ".importance";
	string confusionFilename;
	if (dataset->HasContinuousPhenotypes()) {
		cerr << "File I/O run mode does not support continuous phenotypes yet."
				<< endl;
		return false;
		// confusionFilename = outPrefix + ".confusion";
	} else {
		confusionFilename = outPrefix + ".confusion2";
	}

	/// save the current data set to a temporary file for rjungle
	string tempFile = outPrefix + "_tmp.csv";
	cout << Timestamp() << "Writing temporary file for RJ: "
			<< tempFile	<< endl;
	dataset->WriteNewDataset(tempFile, TAB_DELIMITED_DATASET);

	// this needs to be done for every iteration of Random Jungle
	// so cannot be set once in the constructor
	rjParams.ncol = dataset->NumVariables() + 1;
	rjParams.depVar = rjParams.ncol - 1;
	rjParams.filename  = strdup(tempFile.c_str());
	/// set the default mtry
	if (!fixedMtry) {
		rjParams.mtry = (uli_t) sqrt((double) dataset->NumVariables());
	}

	// create controller
	switch (rjParams.memMode) {
	case 0:
	  RJungleCtrl<double>::autoBuild(rjParams);
	  break;
	case 1:
	  RJungleCtrl<float>::autoBuild(rjParams);
	  break;
	case 2:
	  RJungleCtrl<char>::autoBuild(rjParams);
	  break;
	default:
	  throw Exception(ERRORCODE_39);
	}

	unlink(tempFile.c_str());

	// clean up
	if (colMaskVec != NULL) {
		delete colMaskVec;
	}

	/// loads rjScores map
	cout << Timestamp() << "Loading RJ variable importance (VI) scores "
			<< "from [" << importanceFilename << "]" << endl;
	if (!ReadScores(importanceFilename)) {
		cerr << "ERROR: Could not read Random Jungle scores" << endl;
		return false;
	}

	/// rj classification error
	// from confusion file - 4/11/12
	// added the new getOob - 7/1/12
	classificationAccuracy = ReadClassificationError(confusionFilename);
	cout << Timestamp() << "RJ classification accuracy: "
			<< classificationAccuracy << endl;

	return true;
}

bool RandomJungle::ComputeAttributeScoresRjungle() {
	cout << Timestamp() << "Running rjungle through C system() call" << endl;

	string outPrefix(rjParams.outprefix);
	string importanceFilename = outPrefix + ".importance";
	string confusionFilename;
	if (dataset->HasContinuousPhenotypes()) {
		confusionFilename = outPrefix + ".confusion";
	} else {
		confusionFilename = outPrefix + ".confusion2";
	}

	/// save the current data set to a temporary file for rjungle
	string tempFile = outPrefix + "_tmp.csv";
	cout << Timestamp() << "Writing temporary file for RJ: " << tempFile
			<< endl;
	dataset->WriteNewDataset(tempFile, TAB_DELIMITED_DATASET);

	// base classifier: classification or regression trees?
	pair<RandomJungleTreeType, string> treeTypeResult;
	treeTypeResult = dataset->DetermineTreeType();
	cout << Timestamp() << treeTypeResult.second << ": " << treeTypeResult.first
			<< endl;
	rjParams.treeType = treeTypeResult.first;
	cout << Timestamp() << treeTypeResult.second << ": " << rjParams.treeType
			<< endl;

	/// run rjungle through a system call to the shell
	stringstream rjCmd;
	rjCmd << "rjungle" << " -f " << tempFile << " -e '\t'" << " -D 'Class'"
			<< " -o " << outPrefix << " -U " << rjParams.nthreads << " -t "
			<< rjParams.ntree << " -y " << rjParams.treeType << " -m "
			<< rjParams.mtry << " -B " << rjParams.backSel << " -i "
			<< rjParams.impMeasure << " -j " << rjParams.numOfImpVar;
	if (rjParams.verbose_flag) {
		rjCmd << " -v";
	}
	cout << Timestamp() << "Running RJ command: " << rjCmd.str() << endl;
	int systemCallReturnStatus = system(rjCmd.str().c_str());
	if (systemCallReturnStatus == -1) {
		cerr << "ERROR: Calling rjungle executable. -1 return code" << endl;
		return false;
	}

	/// loads rjScores map from importance file
	cout << Timestamp() << "Loading RJ variable importance (VI) scores "
			<< "from [" << importanceFilename << "]" << endl;
	if (!ReadScores(importanceFilename)) {
		cerr << "ERROR: Could not read Random Jungle scores" << endl;
		return false;
	}

	/// loads rj classification error from confusion file
	cout << Timestamp() << "Loading RJ classification error " << "from ["
			<< confusionFilename << "]" << endl;
	if (!ReadClassificationError(confusionFilename)) {
		cerr << "ERROR: Could not read Random Jungle classification error"
				<< endl;
		return false;
	}

	/// remove the temporary file
	cout << Timestamp() << "Removing temporary file for RJ: " << tempFile
			<< endl;
	unlink(tempFile.c_str());

	return true;
}

bool RandomJungle::ReadScores(string importanceFilename) {
	ifstream importanceStream(importanceFilename.c_str());
	if (!importanceStream.is_open()) {
		cerr << "ERROR: Could not open Random Jungle importance file: "
				<< importanceFilename << endl;
		return false;
	}
	string line;
	// strip the header line
	getline(importanceStream, line);
	// read and store variable name and gini index
	unsigned int lineNumber = 0;
	double minRJScore = 0; // initializations to keep compiler happy
	double maxRJScore = 0; // initializations are on line number 1 of file read
	scores.clear();
	while (getline(importanceStream, line)) {
		++lineNumber;
		vector<string> tokens;
		split(tokens, line, "\t");
		if (tokens.size() != 4) {
			cerr << "ERROR: EvaporativeCooling::ReadRandomJungleScores: "
					<< "error parsing line " << lineNumber << " of "
					<< importanceFilename << ". Read " << tokens.size()
					<< " columns. Should " << "be 4" << endl;
			return false;
		}
		string val = tokens[2];
		string keyVal = tokens[3];
		double key = boost::lexical_cast<double>(keyVal.c_str());
		// cout << setprecision(3);
		// cout << "Storing RJ: " << key << " => " << val << " form " << keyVal << endl;
		scores.push_back(make_pair(key, val));
		if (lineNumber == 1) {
			minRJScore = key;
			maxRJScore = key;
		} else {
			if (key < minRJScore) {
				minRJScore = key;
			} else {
				if (key > maxRJScore) {
					maxRJScore = key;
				}
			}
		}
	}
	importanceStream.close();
	cout << setprecision(3);
	cout << Timestamp() << "Read [" << scores.size() << "] scores from ["
			<< importanceFilename << "], min [" << minRJScore << "], max ["
			<< maxRJScore << "]" << endl;
	// normalize scores
	bool needsNormalization = true;
	if (minRJScore == maxRJScore) {
		cout << Timestamp() << "WARNING: Random Jungle min and max scores "
				<< "are the same" << endl;
		needsNormalization = false;
	}
	double rjRange = maxRJScore - minRJScore;
	vector<pair<double, string> > newRJScores;
	for (unsigned int i = 0; i < scores.size(); ++i) {
		pair<double, string> thisScore = scores[i];
		double key = thisScore.first;
		string val = thisScore.second;
		if (needsNormalization) {
			key = (key - minRJScore) / rjRange;
			newRJScores.push_back(make_pair(key, val));
		} else {
			newRJScores.push_back(make_pair(key, val));
		}
	}

	scores.clear();
	scores = newRJScores;

	return true;
}

bool RandomJungle::ReadClassificationError(std::string confusionFilename) {
	double classifierError = 1.0;
	if (!RandomJungle::ReadClassificationError(confusionFilename,
			(RandomJungleTreeType) rjParams.treeType, classifierError)) {
		return false;
	}
	classificationAccuracy = classifierError;

	return true;
}

bool RandomJungle::GetLibraryClassificationAccuracy() {
	double classifierError = 1.0;
	classificationAccuracy = classifierError;

	return true;
}
