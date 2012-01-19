/*
 * DgeData.cpp
 *
 *  Created on: Jan 18, 2012
 *      Author: billwhite
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <boost/lexical_cast.hpp>

#include "DgeData.h"
#include "Insilico.h"
#include "StringUtils.h"

using namespace std;
using namespace insilico;

DgeData::DgeData() {

}

DgeData::~DgeData() {
}

bool DgeData::LoadData(string countsFile, string phenoFile) {
	countsFilename = countsFile;
	ifstream countsStream(countsFilename.c_str());
	if (!countsStream.is_open()) {
		cerr << "ERROR: Could not open counts file: " << countsFilename << endl;
		return false;
	}
	cout << Timestamp() << "Reading CSV counts from [" << countsFilename << "]"
			<< endl;

	// temporary string for reading file lines
	string line;

	// read the header row - comma delimited sample names
	getline(countsStream, line);
	vector<string> tokens;
	split(tokens, line, ",");
	vector<string>::const_iterator it;
	unsigned int numSamples = 0;
	for (it = tokens.begin(); it != tokens.end(); ++it) {
		// TODO: remove quotes from sample names
		sampleNames.push_back(*it);
		++numSamples;
	}
	cout << Timestamp() << numSamples << " samples read" << endl;

	/// read gene counts, create dummy name for each gene
	unsigned int lineNumber = 0;
	while (getline(countsStream, line)) {
		++lineNumber;
		string trimmedLine = trim(line);
		// no blank lines in the data section
		if (!trimmedLine.size()) {
			cout << "WARNING: Blank line skipped at line number: "
					<< lineNumber << endl;
			continue;
		}
		string geneID = "GENE" + zeroPadNumber(lineNumber, 8);
		// split the line into counts vector
		vector<string> countsStringVector;
		split(countsStringVector, trimmedLine, ",");

		if ((countsStringVector.size() == 0)
				|| (countsStringVector.size() != numSamples)) {
			cout << Timestamp() << "ERROR: Skipping line " << lineNumber
					<< " samples read: " << countsStringVector.size() << " should be "
					<< numSamples << endl;
			return false;
		}

		/// load all counts for this gene as doubles
		geneNames.push_back(geneID);
		vector<string>::const_iterator it = countsStringVector.begin();
		vector<double> geneCounts;
		for (; it != countsStringVector.end(); ++it) {
			string thisStringCount = *it;
			double thisDoubleCount = boost::lexical_cast<double>(thisStringCount);
			geneCounts.push_back(thisDoubleCount);
		}

		/// save this gene's counts to the counts class member variable
		counts.push_back(geneCounts);
	}
	countsStream.close();

	if(!counts.size()) {
		cerr << "ERROR: No genes found" << endl;
		return false;
	}
	else {
		cout << Timestamp() << counts.size() << " genes read" << endl;
	}

	/// read phenotypes
	phenosFilename = phenoFile;
	ifstream phenosStream(phenosFilename.c_str());
	if (!phenosStream.is_open()) {
		cerr << "ERROR: Could not open phenotypes file: " << phenosFilename << endl;
		return false;
	}
	cout << Timestamp() << "Reading phenotypes from [" << phenosFilename << "]"
			<< endl;

	lineNumber = 0;
	while (getline(phenosStream, line)) {
		++lineNumber;
		string trimmedLine = trim(line);
		if (!trimmedLine.size()) {
			cout << "WARNING: Blank line skipped at line number: "
					<< lineNumber << endl;
			continue;
		}
		if((trimmedLine != "0") && (trimmedLine != "1")) {
			cerr << "ERROR: Phenotype is not 0 or 1 on line: " << lineNumber << endl;
			return false;
		}
		int thisPhenotype = boost::lexical_cast<int>(trimmedLine);
		phenotypes.push_back(thisPhenotype);
	}
	phenosStream.close();

	if(phenotypes.size() != numSamples) {
		cerr << "ERROR: Number of phenotypes read: " << phenotypes.size()
				<< " is not equal to the number of sample names read from the"
				<< " counts file: " << numSamples << endl;
		return false;
	}

	cout << Timestamp() << "Read " << numSamples << " samples with counts for "
			<< counts.size() << " genes" << endl;

	return true;
}
