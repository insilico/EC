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

#include <string.h>

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

using namespace std;
using namespace insilico;
using namespace boost;

// static methods
bool RandomJungle::RunClassifier(string csvFile, ConfigMap& vm,
		RandomJungleTreeType treeType, double& classError) {
	string outPrefix(vm["out-files-prefix"]);
	string confusionFilename;
	if((treeType == NUMERIC_NUMERIC_TREE) || (treeType == NUMERIC_NOMINAL_TREE)) {
		confusionFilename = outPrefix + ".confusion";
	}
	else {
		confusionFilename = outPrefix + ".confusion2";
	}

	/// run rjungle through a system call to the shell
	stringstream rjCmd;
	rjCmd << "rjungle"
			<< " -f " << csvFile
			<< " -e ','"
			<< " -D 'Class'"
		 	<< " -o " << outPrefix
		 	<< " -U " << vm["num-threads"]
		 	<< " -t " << vm["rj-num-trees"]
		 	<< " -y " << treeType;
	if(vm["verbose"] == "true") {
		rjCmd << " -v";
	}
	cout << Timestamp() << "Running RJ command: " << rjCmd.str() << endl;
	int systemCallReturnStatus = system(rjCmd.str().c_str());
	if(systemCallReturnStatus == -1) {
		cerr << "ERROR: Calling rjungle executable. -1 return code" << endl;
		return false;
	}

	/// loads rj classification error from confusion file
	cout << Timestamp() << "Loading RJ classification error "
			<< "from [" << confusionFilename << "]" << endl;
	if (!ReadClassificationError(confusionFilename, treeType, classError)) {
		cerr << "ERROR: Could not read Random Jungle classification error" << endl;
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
	if((treeType == NUMERIC_NUMERIC_TREE) ||
			(treeType == NUMERIC_NOMINAL_TREE) ||
			(treeType == NOMINAL_NUMERIC_TREE) ||
			(treeType == NOMINAL_NUMERIC_FLOATS)) {
		for(unsigned int i=0; i < 8; ++i) {
			getline(confusionStream, line);
		}
		vector<string> tokens;
		split(tokens, line, " ");
		if (tokens.size() != 3) {
			cout << Timestamp() << "WARNING: RandomJungle::GetClassificationError: "
					<< "error parsing " << confusionFilename << "." << endl;
			cout << Timestamp() << "WARNING: Read " << tokens.size()
					<< " columns from line 8, should be 3" << endl;
			cout << Timestamp() << "WARNING: Attempting to read alternate "
					<< "confusion file column" << endl;
		}
		if (tokens.size() == 3) {
			classifierError = lexical_cast<double>(trim(tokens[2]));
		}
		else {
			classifierError = lexical_cast<double>(trim(tokens[0]));
		}
	}
	else {
		getline(confusionStream, line);
		vector<string> tokens;
		split(tokens, line, ",");
		if (tokens.size() != 5) {
			cerr << "ERROR: RandomJungle::GetClassificationError: "
					<< "error parsing " << confusionFilename << "." << endl;
			cerr << "Read " << tokens.size()
					<< " columns from line 2, should be 5" << endl;
			return false;
		}
		classifierError = lexical_cast<double>(trim(tokens[4]));
	}
	confusionStream.close();

	cout << Timestamp() << "Read classification error from " << confusionFilename
			<< ": " << setprecision(4) << classifierError << endl;

	return true;
}

// standard methods
RandomJungle::RandomJungle(Dataset* ds, po::variables_map& vm) {
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

	if (vm.count("rj-num-trees")) {
		rjParams.ntree = vm["rj-num-trees"].as<uli_t>();
	} else {
		rjParams.ntree = 1000;
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

	rjParams.rng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rjParams.rng, rjParams.seed);

	rjParams.nrow = dataset->NumInstances();
	rjParams.depVarName = (char *) "Class";
	//  rjParams.verbose_flag = true;
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

RandomJungle::RandomJungle(Dataset* ds, ConfigMap& configMap) {

	cout << Timestamp() << "Initializing Random Jungle with ConfigMap" << endl;

	dataset = ds;

	rjParams = initRJunglePar();

	string configValue;

	// how should we set this default value?
	uli_t numTrees = 1000;
	if (GetConfigValue(configMap, "rj-num-trees", configValue)) {
		numTrees = lexical_cast<unsigned int>(configValue);
	}
	rjParams.ntree = numTrees;

	// fill in the parameters object for the RJ run
	rjParams.rng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rjParams.rng, rjParams.seed);

	rjParams.nrow = dataset->NumInstances();
	rjParams.depVarName = (char *) "Class";
	//  rjParams.verbose_flag = true;
	rjParams.filename = (char*) "";

	rjParams.mpiId = 0;

	rjParams.verbose_flag =
			GetConfigValue(configMap, "verbose", configValue)? true : false;

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
	if (rjParams.rng) {
		gsl_rng_free(rjParams.rng);
	}
}

bool RandomJungle::ComputeAttributeScores() {

	cout << Timestamp() << "Computing Random Jungle variable importance scores"
			<< endl;

	if(runMode == SYSTEM_CALL_RUN_MODE) {
		return ComputeAttributeScoresRjungle();
	}

	cout << Timestamp() << "Running Random Jungle using C++ librjungle calls"
			<< endl;

	vector<uli_t>* colMaskVec = NULL;
	time_t start, end;
	clock_t startgrow, endgrow;

	time(&start);

	rjParams.ncol = dataset->NumVariables() + 1;
	rjParams.depVar = rjParams.ncol - 1;
	rjParams.depVarCol = rjParams.ncol - 1;
	string outPrefix(rjParams.outprefix);
	string importanceFilename = outPrefix + ".importance";
	string confusionFilename = outPrefix + ".confusion";

	RJungleIO io;
	io.open(rjParams);
	unsigned int numInstances = dataset->NumInstances();
	vector<string> variableNames = dataset->GetVariableNames();
	variableNames.push_back(rjParams.depVarName);
	vector<string> instanceIds = dataset->GetInstanceIds();

//  cout << "Variables names from the data set:" << endl;
//  copy(variableNames.begin(), variableNames.end(), ostream_iterator<string>(cout, "\n"));

	vector<string> attributeNames = dataset->GetAttributeNames();
	vector<string> numericNames = dataset->GetNumericsNames();
	if ((rjParams.treeType == 1) || (rjParams.treeType == 3)
			|| (rjParams.treeType == 4) || (rjParams.treeType == 5)) {
		// regression
		cout << Timestamp() << "Preparing regression trees Random Jungle" << endl;
		DataFrame<NumericLevel>* data = new DataFrame<NumericLevel>(rjParams);
		data->setDim(rjParams.nrow, rjParams.ncol);
		data->setVarNames(variableNames);
		data->setDepVarName(rjParams.depVarName);
		data->setDepVar(rjParams.depVarCol);
		data->initMatrix();
		// load data frame
		// TODO: do not load data frame every time-- use column mask mechanism?
		cout << Timestamp() << "Loading RJ DataFrame with double values: " << endl
				<< Timestamp();
		cout.flush();
		for (unsigned int i = 0; i < numInstances; ++i) {
			unsigned int instanceIndex;
			dataset->GetInstanceIndexForID(instanceIds[i], instanceIndex);
			unsigned int j = 0;
			for (unsigned int a = 0; a < attributeNames.size(); ++a) {
//        cout << "Loading instance " << i << ", attribute: " << a
//                << " " << attributeNames[a] << endl;
				data->set(
						i,
						j,
						static_cast<NumericLevel>(dataset->GetAttribute(instanceIndex,
								attributeNames[a])));
				++j;
			}
			for (unsigned int n = 0; n < numericNames.size(); ++n) {
				data->set(i, j, dataset->GetNumeric(i, numericNames[n]));
				++j;
			}
			if (dataset->HasContinuousPhenotypes()) {
				data->set(i, j,
						dataset->GetInstance(instanceIndex)->GetPredictedValueTau());
			} else {
				data->set(
						i,
						j,
						static_cast<NumericLevel>(dataset->GetInstance(instanceIndex)->GetClass()));
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
		data->storeCategories();
		data->makeDepVecs();
		data->getMissings();
//   cout << "DEBUG data frame:" << endl;
//    data->print(cout);
//    data->printSummary();

		RJungleGen<NumericLevel> rjGen;
		rjGen.init(rjParams, *data);

		RJungleHelper<NumericLevel>::printRJunglePar(rjParams, *io.outLog);

		startgrow = clock();
		// create controller
		RJungleCtrl<NumericLevel> rjCtrl;
		cout << Timestamp() << "Running Random Jungle" << endl;
		rjCtrl.autoBuildInternal(rjParams, io, rjGen, *data, colMaskVec);
		endgrow = clock();

		time(&end);

		// print info stuff
		RJungleHelper<NumericLevel>::printFooter(rjParams, io, start, end,
				startgrow, endgrow);
		delete data;
	} else {
		cout << Timestamp() << "Preparing SNP classification trees Random Jungle"
				<< endl;
		DataFrame<char>* data = new DataFrame<char>(rjParams);
		data->setDim(rjParams.nrow, rjParams.ncol);
		data->setVarNames(variableNames);
		data->setDepVarName(rjParams.depVarName);
		data->setDepVar(rjParams.depVarCol);
		data->initMatrix();
		// TODO: do not load data frame every time-- use column mask mechanism?
		cout << Timestamp() << "Loading RJ DataFrame with character values: "
				<< endl << Timestamp();
		for (unsigned int i = 0; i < rjParams.nrow; ++i) {
			unsigned int instanceIndex;
			dataset->GetInstanceIndexForID(instanceIds[i], instanceIndex);
			// set attributes from data set
			for (unsigned int j = 0; j < rjParams.ncol - 1; ++j) {
				AttributeLevel intAttr = dataset->GetAttribute(instanceIndex,
						variableNames[j]);
				char attr = ' ';
				switch (intAttr) {
				case 0:
					attr = 0x0;
					break;
				case 1:
					attr = 0x1;
					break;
				case 2:
					attr = 0x2;
					break;
				case MISSING_ATTRIBUTE_VALUE:
					// Random Jungle says to avoid missing values. Duh!
					// What does this mean? Trying question mark.
					attr = '?';
					break;
				}
//        cout << "Setting (" << i << "," << j << ") => "
//                << intAttr << ", " << attr << endl;
				data->set(i, j, attr);
			}

			// set class from data set
			unsigned int intClass = dataset->GetInstance(instanceIndex)->GetClass();
			char classVal = ' ';
			if (intClass == 0) {
				classVal = 0x0;
			} else {
				classVal = 0x1;
			}
//      cout << "Setting class value (" << i << "," << (rjParams.ncol - 1)
//              << ") => " << intClass << ", " << classVal << endl;
			data->set(i, rjParams.ncol - 1, classVal);

			// happy lights
			if (i && ((i % 100) == 0)) {
				cout << i << "/" << numInstances << " ";
				cout.flush();
			}
			if (i && ((i % 1000) == 0)) {
				cout << i << endl << Timestamp();
			}
		}
		cout << numInstances << "/" << numInstances << endl;
		data->storeCategories();
		data->makeDepVecs();
		data->getMissings();
		// data->print(cout);

		RJungleGen<char> rjGen;
		rjGen.init(rjParams, *data);

		RJungleHelper<char>::printRJunglePar(rjParams, *io.outLog);

		startgrow = clock();
		// create controller
		RJungleCtrl<char> rjCtrl;
		cout << Timestamp() << "Running Random Jungle" << endl;
		rjCtrl.autoBuildInternal(rjParams, io, rjGen, *data, colMaskVec);
		endgrow = clock();

		time(&end);

		// print info stuff
		RJungleHelper<char>::printFooter(rjParams, io, start, end, startgrow,
				endgrow);
		delete data;
	}

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
		return false;
	}

	/// loads rj classification error from confusion file - 4/11/12
	cout << Timestamp() << "Loading RJ classification error "
			<< "from [" << confusionFilename << "]" << endl;
	if (!ReadClassificationError(confusionFilename)) {
		cerr << "ERROR: Could not read Random Jungle classification error" << endl;
		return false;
	}

	return true;
}

bool RandomJungle::ComputeAttributeScoresRjungle() {
	cout << Timestamp() << "Running rjungle through C system() call" << endl;

	string outPrefix(rjParams.outprefix);
	string importanceFilename = outPrefix + ".importance";
	string confusionFilename;
	if(dataset->HasContinuousPhenotypes()) {
		confusionFilename = outPrefix + ".confusion";
	}
	else {
		confusionFilename = outPrefix + ".confusion2";
	}

	/// save the current data set to a temporary file for rjungle
	string tempFile = outPrefix + "_tmp.csv";
	cout << Timestamp() << "Writing temporary file for RJ: " << tempFile << endl;
	dataset->WriteNewDataset(tempFile, CSV_DELIMITED_DATASET);

	// base classifier: classification or regression trees?
	pair<RandomJungleTreeType, string> treeTypeResult;
	treeTypeResult = dataset->DetermineTreeType();
	cout << Timestamp() << treeTypeResult.second << ": "
			<< treeTypeResult.first << endl;
	rjParams.treeType = treeTypeResult.first;
	cout << Timestamp() << treeTypeResult.second << ": "
			<< rjParams.treeType << endl;

	/// run rjungle through a system call to the shell
	stringstream rjCmd;
	rjCmd << "rjungle"
			<< " -f " << tempFile
			<< " -e ','"
			<< " -D 'Class'"
		 	<< " -o " << outPrefix
		 	<< " -U " << rjParams.nthreads
		 	<< " -t " << rjParams.ntree
		 	<< " -y " << rjParams.treeType;
	if(rjParams.verbose_flag) {
		rjCmd << " -v";
	}
	cout << Timestamp() << "Running RJ command: " << rjCmd.str() << endl;
	int systemCallReturnStatus = system(rjCmd.str().c_str());
	if(systemCallReturnStatus == -1) {
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
	cout << Timestamp() << "Loading RJ classification error "
			<< "from [" << confusionFilename << "]" << endl;
	if (!ReadClassificationError(confusionFilename)) {
		cerr << "ERROR: Could not read Random Jungle classification error" << endl;
		return false;
	}

	/// remove the temporary file
	cout << Timestamp() << "Removing temporary file for RJ: " << tempFile << endl;
	unlink(tempFile.c_str());

	return true;
}

vector<pair<double, string> > RandomJungle::GetScores() {
	return scores;
}

double RandomJungle::GetClassificationError() {
	return classificationError;
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
		split(tokens, line, " ");
		if (tokens.size() != 4) {
			cerr << "ERROR: EvaporativeCooling::ReadRandomJungleScores: "
					<< "error parsing line " << lineNumber << " of " << importanceFilename
					<< ". Read " << tokens.size() << " columns. Should " << "be 4"
					<< endl;
			return false;
		}
		string val = tokens[2];
		string keyVal = tokens[3];
		double key = strtod(keyVal.c_str(), NULL);
		// cout << "Storing RJ: " << key << " => " << val << endl;
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
	cout << Timestamp() << "Read [" << scores.size() << "] scores from ["
			<< importanceFilename << "]" << endl;
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
	if(!RandomJungle::ReadClassificationError(confusionFilename,
			(RandomJungleTreeType) rjParams.treeType, classifierError)) {
		return false;
	}
	classificationError = classifierError;

	return true;
}

