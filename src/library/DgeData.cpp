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
#include <algorithm>

#include <boost/lexical_cast.hpp>

#include "DgeData.h"
#include "Insilico.h"
#include "StringUtils.h"

using namespace std;

DgeData::DgeData() {
	hasNormFactors = false;
}

DgeData::~DgeData() {
}

bool DgeData::LoadData(string countsFile, string normsFile) {

	// temporary string for reading file lines
	string line;

	if(normsFile != "") {
		normsFilename = normsFile;
		ifstream normsStream(normsFilename.c_str());
		if (!normsStream.is_open()) {
			cerr << "ERROR: Could not open normalization factors file: "
					<< normsFilename << endl;
			return false;
		}
		cout << Timestamp() << "Reading normalization factors from ["
				<< normsFilename << "]" << endl;
		int lineNumber = 0;
		while (getline(normsStream, line)) {
			++lineNumber;
			string trimmedLine = insilico::trim(line);
			if (!trimmedLine.size()) {
				cout << "WARNING: Blank line skipped at line number: "
						<< lineNumber << endl;
				continue;
			}
			double thisFactor = boost::lexical_cast<double>(trimmedLine);
			normFactors.push_back(thisFactor);
		}
		normsStream.close();
		hasNormFactors = true;
	}

	countsFilename = countsFile;
	ifstream countsStream(countsFilename.c_str());
	if (!countsStream.is_open()) {
		cerr << "ERROR: Could not open counts file: " << countsFilename << endl;
		return false;
	}
	cout << Timestamp() << "Reading CSV counts from [" << countsFilename << "]"
			<< endl;


	// read the header row - comma delimited sample phenotypes
	getline(countsStream, line);
	vector<string> tokens;
	insilico::split(tokens, line, ",");
	vector<string>::const_iterator it;
	unsigned int numPhenos = 0;
	unsigned int numSamples = 0;
	// skip the first column; the rest are phenotypes corresponding to subjects
	for (it = tokens.begin(); it != tokens.end(); ++it) {
		// Remove quotes from sample names
		string phenoString = *it;
		phenoString.erase(remove(phenoString.begin(), phenoString.end(), '"'), phenoString.end());
		phenotypes.push_back(boost::lexical_cast<int>(phenoString));
		++numPhenos;
		// use dummy subject/sample names
		ostringstream ss;
		ss << "Sample" << (numSamples + 1);
		sampleNames.push_back(ss.str());
		++numSamples;
	}
	cout << Timestamp() << phenotypes.size() << "/" << sampleNames.size()
			<< " samples/phenotypes read from file header" << endl;

	if(hasNormFactors && (normFactors.size() != phenotypes.size())) {
		cerr << "Number of phenotypes: " << phenotypes.size()
				<< " is not equal to the number of normalization factors: "
				<< normFactors.size() << endl;
		return false;
	}

	/// read gene counts
	unsigned int lineNumber = 0;
	while (getline(countsStream, line)) {
		++lineNumber;
		string trimmedLine = insilico::trim(line);
		// no blank lines in the data section
		if (!trimmedLine.size()) {
			cout << "WARNING: Blank line skipped at line number: "
					<< lineNumber << endl;
			continue;
		}
		// split the line into counts vector
		vector<string> countsStringVector;
		insilico::split(countsStringVector, trimmedLine, ",");
		unsigned int countsRead = countsStringVector.size() - 1;
		if (countsRead == 0) {
			cerr << "ERROR: Line: " << lineNumber << " could not be parsed" << endl;
			return false;
		}
		if(countsRead != numSamples) {
			cerr << "ERROR: Line: " << lineNumber
					<< " counts read: " << countsRead
					<< " should be " << numSamples << endl;
			return false;
		}

		// first column is gene ID
		string geneID = countsStringVector[0];
		geneNames.push_back(geneID);

		/// load all counts for this gene as doubles
		vector<string>::const_iterator it = countsStringVector.begin() + 1;
		vector<double> geneCounts;
		double minCount = 0;
		double maxCount = 0;
		int geneIndex = 0;
		for (; it != countsStringVector.end(); ++it, ++geneIndex) {
			string thisStringCount = *it;
			double thisDoubleCount = boost::lexical_cast<double>(thisStringCount);
			if(hasNormFactors) {
				thisDoubleCount *= normFactors[geneIndex];
			}
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
	for(size_t geneIndex=0; geneIndex < counts.size(); ++geneIndex) {
		vector<double> thisGeneCounts = counts[geneIndex];
		for(int sampleIndex=0; sampleIndex < (int) thisGeneCounts.size(); ++sampleIndex) {
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

	if(phenotypes.size() != numSamples) {
		cerr << "ERROR: Number of phenotypes read: " << phenotypes.size()
				<< " is not equal to the number of sample names read from the"
				<< " counts file: " << numSamples << endl;
		return false;
	}

	cout << Timestamp() << "Read " << numSamples << " samples with counts for "
			<< counts.size() << " genes";
	if(hasNormFactors) {
		cout << " (normalized)";
	}
	cout << endl;

	return true;
}

vector<string> DgeData::GetSampleNames() {
	return sampleNames;
}

vector<string> DgeData::GetGeneNames() {
	return geneNames;
}

pair<double, double> DgeData::GetGeneMinMax(int geneIndex) {
	if((geneIndex >= 0) && (geneIndex < (int) counts.size())) {
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
	if((sampleIndex < 0) || (sampleIndex >= (int) sampleNames.size())) {
		cerr << "ERROR: DgeData::GetSampleCounts, index out of range: "
				<< sampleIndex << endl;
		exit(EXIT_FAILURE);
	}

	vector<double> returnVector;
	for(int i=0; i < (int) geneNames.size(); ++i) {
		returnVector.push_back(counts[i][sampleIndex]);
	}

	return returnVector;
}

int DgeData::GetSamplePhenotype(int sampleIndex) {
	if((sampleIndex < 0) || (sampleIndex >= (int) sampleNames.size())) {
		cerr << "ERROR: DgeData::GetSamplePhenotype, index out of range: "
				<< sampleIndex << endl;
		exit(EXIT_FAILURE);
	}

	return phenotypes[sampleIndex];
}

vector<double> DgeData::GetNormalizationFactors() {
	return normFactors;
}

void DgeData::PrintSampleStats() {
	cout << Timestamp() << "DGE Sample Statistics" << endl;
	for(int sampleIndex = 0; sampleIndex < (int) sampleNames.size(); ++sampleIndex) {
		cout << Timestamp()
				 << sampleNames[sampleIndex]
		     << " min: " << minMaxSampleCounts[sampleIndex].first
		     << " max: " << minMaxSampleCounts[sampleIndex].second
		     << " zeroes: " << sampleZeroes[sampleIndex].size()
		     << endl;
	}
}

