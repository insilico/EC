/*
 * File:   Insilico.cpp
 * Author: billwhite
 *
 * Created on October 13, 2011, 10:02 PM
 *
 * Common functions for Insilico Lab projects.
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <ctime>
#include <cmath>

#include "Dataset.h"
#include "DatasetInstance.h"
#include "ArffDataset.h"
#include "ArffDataset.h"
#include "PlinkDataset.h"
#include "PlinkRawDataset.h"
#include "PlinkBinaryDataset.h"
#include "StringUtils.h"
#include "Insilico.h"

using namespace std;
using namespace insilico;

RandomJungleTreeType DetermineRandomJungleTreeType(AttributeType attributeType,
		ClassType classType) {

	RandomJungleTreeType treeType = UNKNOWN_TREE_TYPE;
	bool treeTypeNotFound = true;
	if(treeTypeNotFound &&
			(classType == CASE_CONTROL_CLASS_TYPE) &&
			(attributeType == NUMERIC_TYPE)) {
		treeType = NOMINAL_NUMERIC_TREE;
		treeTypeNotFound = false;
	}
	if(treeTypeNotFound &&
			(classType == CASE_CONTROL_CLASS_TYPE) &&
			(attributeType == DISCRETE_TYPE)) {
		treeType = NOMINAL_NOMINAL_TREE;
		treeType = NOMINAL_NUMERIC_TREE;
		treeTypeNotFound = false;
	}
	if(treeTypeNotFound &&
			(classType == CONTINUOUS_CLASS_TYPE) &&
			(attributeType == NUMERIC_TYPE)) {
		treeType = NUMERIC_NUMERIC_TREE;
		treeTypeNotFound = false;
	}
	if(treeTypeNotFound &&
			(classType == CONTINUOUS_CLASS_TYPE) &&
			(attributeType == DISCRETE_TYPE)) {
		treeType = NUMERIC_NOMINAL_TREE;
		treeTypeNotFound = false;
	}

	if(treeTypeNotFound) {
		cout << "WARNING: Random Jungle tree type could not be determined" << endl;
	}

	return treeType;
}

string Timestamp() {

	time_t now = time(NULL);
	struct tm * ptm = localtime(&now);
	char buffer[32];
	// Format: 20091506 - 20:20:00 -
	strftime(buffer, 32, "%Y%m%d - %H:%M:%S - ", ptm);
	string ts(buffer);

	return ts;
}

Dataset* ChooseSnpsDatasetByType(string snpsFilename, string snpsFileType) {
	// determine the data set type from the passed SNP file type or from
	// the file extension if the file type is the empty string
	Dataset* ds = 0;

	if(datasetTypeToExt.size() == 0) {
		datasetTypeToExt["textwhitesp"] = "txt";
		datasetTypeToExt["plinkped"] = "ped";
		datasetTypeToExt["plinkbed"] = "bed";
		datasetTypeToExt["plinkraw"] = "raw";
		datasetTypeToExt["mayogeo"] = "csv";
		datasetTypeToExt["birdseed"] = "txt";
	}

	string fileExt = "";
	if(snpsFileType != "") {
		map<string, string>::const_iterator foundPos =
				datasetTypeToExt.find(snpsFileType);
		if(foundPos != datasetTypeToExt.end()) {
			fileExt = foundPos->second;
		}
		else {
			cerr << "ERROR: ChooseSnpsDatasetByType: cannot determine SNPs file "
					<< " type from type string: " << snpsFileType << endl;
			return ds;
		}
	}
	else {
		fileExt = GetFileExtension(snpsFilename);
	}
	// cout << "File extension: " << fileExt << endl;

	cout << Timestamp() << "Dataset detection for SNP file ["
			<< snpsFilename << "]" << endl;
	if (fileExt == "arff") {
		cout << Timestamp() << "ARFF" << endl;
		ds = new ArffDataset();
	} else {
		if ((fileExt == "dat") || (fileExt == "tab") || (fileExt == "txt")) {
			cout << Timestamp() << "Whitespace-delimited" << endl;
			ds = new Dataset();
		} else {
			if ((fileExt == "ped") || (fileExt == "map")) {
				cout << Timestamp() << "Plink map/ped" << endl;
				ds = new PlinkDataset();
			} else {
				if (fileExt == "raw") {
					cout << Timestamp() << "Plink raw" << endl;
					ds = new PlinkRawDataset();
				} else {
					if ((fileExt == "bed") || (fileExt == "bim") || (fileExt == "fam")) {
						cout << Timestamp() << "Plink binary" << endl;
						ds = new PlinkBinaryDataset();
					} else {
						cerr << endl;
						cerr << "ERROR: Cannot determine data set type by extension: "
								<< fileExt <<  " for " << snpsFilename << endl;
					}
				}
			}
		}
	}

	return ds;
}

bool LoadNumericIds(string filename, vector<string>& retIds) {
	ifstream dataStream(filename.c_str());
	if (!dataStream.is_open()) {
		cerr << "ERROR: Could not open ID file: " << filename << endl;
		return false;
	}

	// temporary string for reading file lines
	string line;

	// read the header
	getline(dataStream, line);
	vector<string> numNames;
	split(numNames, line);
	if (numNames.size() < 3) {
		cerr << "ERROR: ID file must have at least three columns: "
				<< "FID IID VAR1 . . . VARn" << endl;
		return false;
	}

	// read each line of the file and get the first tab-delimited field as the
	// individual's ID; insure each line has at least three tab-delimited fields
	retIds.clear();
	map<string, bool> idsSeen;
	unsigned int lineNumber = 0;
	retIds.clear();
	while (getline(dataStream, line)) {
		++lineNumber;
		vector<string> fieldsStringVector;
		split(fieldsStringVector, line);
		if (fieldsStringVector.size() < 3) {
			cerr << "ERROR: ID file must have at least three columns: "
					<< "FID IID VAR1 ... VARn" << endl;
			return false;
		}
    // use the PLINK IID as the instance ID - 9/18/12
		// use both FID and IID so all PLINK files will work - 4/10/13
		string ID = trim(fieldsStringVector[0]) + trim(fieldsStringVector[1]);
		if (idsSeen.find(ID) == idsSeen.end()) {
			idsSeen[ID] = true;
			retIds.push_back(ID);
		} else {
			cout << Timestamp() << "WARNING: Duplicate ID [" << ID
					<< "] detected and " << "skipped on line [" << lineNumber << "]"
					<< endl;
		}
	}
	dataStream.close();

	cout << Timestamp() << retIds.size() << " numeric IDs read" << endl;

	return true;
}

bool LoadPhenoIds(string filename, vector<string>& retIds) {
	ifstream dataStream(filename.c_str());
	if (!dataStream.is_open()) {
		cerr << "ERROR: Could not open ID file: " << filename << endl;
		return false;
	}

	// temporary string for reading file lines
	string line;

	// read each line of the file and get the first tab-delimited field as the
	// individual's ID; insure each line has three tab-delimited fields
	retIds.clear();
	map<string, bool> idsSeen;
	unsigned int lineNumber = 0;
	retIds.clear();
	int numPhenosAccepted = 0;
	while (getline(dataStream, line)) {
		++lineNumber;
		vector<string> fieldsStringVector;
		split(fieldsStringVector, line);
		if (fieldsStringVector.size() != 3) {
			cerr << "ERROR: ID file must have three columns: " << "FID IID PHENO"
					<< endl;
			return false;
		}
    // use the PLINK IID as the instance ID - 9/18/12
		// use both FID and IID so all PLINK files will work - 4/10/13
		// string ID = trim(fieldsStringVector[1]);
		string ID = trim(fieldsStringVector[0]) + trim(fieldsStringVector[1]);
		string pheno = trim(fieldsStringVector[2]);
		if (pheno == "-9") {
			continue;
		}
		if (idsSeen.find(ID) == idsSeen.end()) {
			idsSeen[ID] = true;
			retIds.push_back(ID);
			++numPhenosAccepted;
		} else {
			cout << Timestamp() << "WARNING: Duplicate ID [" << ID
					<< "] detected and " << "skipped on line [" << lineNumber << "]"
					<< endl;
		}
	}
	dataStream.close();

	cout << Timestamp() << numPhenosAccepted << " non-missing phenotype IDs read"
			<< endl;

	return true;
}

bool GetMatchingIds(string numericsFilename, string altPhenotypeFilename,
		vector<string> numericsIds, vector<string> phenoIds,
		vector<string>& matchingIds) {
	if (numericsFilename != "" && altPhenotypeFilename != "") {
		cout << Timestamp() << "IDs come from the numeric and the alternate "
				<< "phenotype files. Checking for intersection/matches" << endl;

		// find intersection of numeric and phenotype IDs
		unsigned int maxMatches = max(numericsIds.size(), phenoIds.size());
		unsigned int maxMismatches = numericsIds.size() + phenoIds.size();
		matchingIds.resize(maxMatches);
		vector<string>::iterator goodIdsIt;
		vector<string> skippedIds(maxMismatches);
		vector<string>::iterator badIdsIt;

		sort(numericsIds.begin(), numericsIds.end());
		sort(phenoIds.begin(), phenoIds.end());
		goodIdsIt = set_intersection(numericsIds.begin(), numericsIds.end(),
				phenoIds.begin(), phenoIds.end(), matchingIds.begin());
		badIdsIt = set_difference(numericsIds.begin(), numericsIds.end(),
				phenoIds.begin(), phenoIds.end(), skippedIds.begin());
		if (skippedIds.begin() != badIdsIt) {
			cerr << "\t\t\tWARNING: Covariates and phenotypes files do not contain "
					<< "the same IDs. These IDs differ: ";
			vector<string>::const_iterator skippedIt = skippedIds.begin();
			for (; skippedIt != badIdsIt; ++skippedIt) {
				cerr << *skippedIt << " ";
			}
			cerr << endl;
			return false;
		}
		matchingIds.resize(int(goodIdsIt - matchingIds.begin()));
	} else {
		if (numericsFilename != "") {
			cout << Timestamp() << "IDs come from the numerics file" << endl;
			matchingIds.resize(numericsIds.size());
			copy(numericsIds.begin(), numericsIds.end(), matchingIds.begin());
			// copy(numericsIds.begin(), numericsIds.end(), ostream_iterator<string > (cout, "\n"));
		} else {
			if (altPhenotypeFilename != "") {
				cout << Timestamp() << "IDs come from the alternate phenotype file"
						<< endl;
				// PrintVector(phenoIds, "phenoIds");
				matchingIds.resize(phenoIds.size());
				copy(phenoIds.begin(), phenoIds.end(), matchingIds.begin());
			} else {
				cout << Timestamp() << "IDs are not needed for this analysis" << endl;
			}
		}
	}

	return true;
}

ClassType DetectClassType(std::string filename, int classColumn,
		bool hasHeader) {
	ClassType detectedClass = NO_CLASS_TYPE;

	/// Open the file for reading
	ifstream phenoDataStream(filename.c_str());
	pair<map<string, unsigned int>::iterator, bool> retClassInsert;
	if (!phenoDataStream.is_open()) {
		cerr << "ERROR: DetectClassType: Could not open file: " << filename << endl;
		return detectedClass;
	}
	cout << Timestamp() << "Detecting class type from file: " << filename << endl;

	string line;

	/// Skip the header if it has one
	if(hasHeader) {
		getline(phenoDataStream, line);
	}

	/// Determine the phenotype type
	int lineNumber = 0;
	bool decimalFound = false;
	map<string, int> histogram;
	while (getline(phenoDataStream, line)) {
		++lineNumber;
		string trimmedLine = trim(line);
		if (trimmedLine[0] == '#') {
			continue;
		}
		if (trimmedLine.size() == 0) {
			continue;
		}
		vector<string> tokens;
		split(tokens, trimmedLine);
		int numColumns = tokens.size();
		if ((classColumn < 0) || (classColumn > numColumns)) {
			cerr << "ERROR: DetectClassType: reading file line " << lineNumber << ". "
					<< "Class column out of range: " << classColumn << endl;
			return detectedClass;
		}
		string thisClassString = tokens[classColumn - 1];
		if (thisClassString == "-9") {
			continue;
		}
		if (thisClassString.find('.') != string::npos) {
			decimalFound = true;
		}
		++histogram[thisClassString];
	}
	phenoDataStream.close();

	if (decimalFound) {
		detectedClass = CONTINUOUS_CLASS_TYPE;
	} else {
		if (histogram.size() == 2) {
			detectedClass = CASE_CONTROL_CLASS_TYPE;
		} else {
			if (histogram.size() > 2) {
				detectedClass = MULTI_CLASS_TYPE;
			}
		}
	}

	return detectedClass;
}

bool GetConfigValue(ConfigMap& configMap, std::string key, std::string& value) {
	if(configMap.find(key) != configMap.end()) {
		value = configMap[key];
		string trimmedKey = trim(value);
		if(trimmedKey == "") {
			return false;
		}
		return true;
	}
	return false;
}

string GetFileBasename(string fileName) {
  size_t pos = fileName.rfind('.');
  return fileName.substr(0, pos);
}

string GetFileExtension(string fileName) {
  size_t pos = fileName.rfind('.');
  return fileName.substr(pos + 1, fileName.size() - 1);
}

double ProtectedLog(double x) {
	if(x <= 0.0) {
		return 0.0;
	}
	else {
		return log(x);
	}
}
