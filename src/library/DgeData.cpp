/*
 * DgeData.cpp
 *
 *  Created on: Jan 18, 2012
 *      Author: billwhite
 */

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
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
//		string geneID = "GENE" + zeroPadNumber(lineNumber, 8);
		stringstream ssGeneID;
		ssGeneID << "GENE" << lineNumber;
		string geneID = ssGeneID.str();
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
		double minCount = 0;
		double maxCount = 0;
		int geneIndex = 0;
		for (; it != countsStringVector.end(); ++it, ++geneIndex) {
			string thisStringCount = *it;
			double thisDoubleCount = boost::lexical_cast<double>(thisStringCount);
			if(geneIndex == 0)  {
				minCount = maxCount = thisDoubleCount;
			}
			else {
				if(thisDoubleCount < minCount) {
					minCount = thisDoubleCount;
				}
				if(thisDoubleCount > maxCount) {
					maxCount = thisDoubleCount;
				}
			}
			geneCounts.push_back(thisDoubleCount);
		}
		minMaxGeneCounts.push_back(make_pair(minCount, maxCount));

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

	/// get min and max sample counts, and sample zeroes
	minMaxSampleCounts.resize(numSamples);
	sampleZeroes.resize(numSamples);
	for(int geneIndex=0; geneIndex < counts.size(); ++geneIndex) {
		vector<double> thisGeneCounts = counts[geneIndex];
		for(int sampleIndex=0; sampleIndex < thisGeneCounts.size(); ++sampleIndex) {
			double thisCount = thisGeneCounts[sampleIndex];
			if(geneIndex==0) {
				minMaxSampleCounts[sampleIndex].first = thisCount;
				minMaxSampleCounts[sampleIndex].second = thisCount;
			}
			else {
				if(thisCount < minMaxSampleCounts[sampleIndex].first) {
					minMaxSampleCounts[sampleIndex].first = thisCount;
				}
				if(thisCount > minMaxSampleCounts[sampleIndex].second) {
					minMaxSampleCounts[sampleIndex].second = thisCount;
				}
			}
			if(thisCount == 0) {
				sampleZeroes[sampleIndex].push_back(geneIndex);
			}
		}
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

vector<string> DgeData::GetSampleNames() {
	return sampleNames;
}

vector<string> DgeData::GetGeneNames() {
	return geneNames;
}

pair<double, double> DgeData::GetGeneMinMax(int geneIndex) {
	if((geneIndex >= 0) && (geneIndex < counts.size())) {
		return minMaxGeneCounts[geneIndex];
	}
	else {
		cerr << "ERROR: DgeData::GetGeneMinMax, index out of range: "
				<< geneIndex << endl;
		exit(EXIT_FAILURE);
	}
}

int DgeData::GetNumSamples() {
	return sampleNames.size();
}

int DgeData::GetNumGenes() {
	return geneNames.size();
}

vector<double> DgeData::GetSampleCounts(int sampleIndex) {
	if((sampleIndex < 0) || (sampleIndex >= sampleNames.size())) {
		cerr << "ERROR: DgeData::GetSampleCounts, index out of range: "
				<< sampleIndex << endl;
		exit(EXIT_FAILURE);
	}

	vector<double> returnVector;
	for(int i=0; i < geneNames.size(); ++i) {
		returnVector.push_back(counts[i][sampleIndex]);
	}

	return returnVector;
}

int DgeData::GetSamplePhenotype(int sampleIndex) {
	if((sampleIndex < 0) || (sampleIndex >= sampleNames.size())) {
		cerr << "ERROR: DgeData::GetSamplePhenotype, index out of range: "
				<< sampleIndex << endl;
		exit(EXIT_FAILURE);
	}

	return phenotypes[sampleIndex];
}
