/*
 * BirdseedData.cpp
 *
 *  Created on: Feb 12, 2012
 *      Author: billwhite
 */

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

#include <boost/lexical_cast.hpp>

#include "BirdseedData.h"
#include "Insilico.h"
#include "StringUtils.h"

using namespace std;
using namespace insilico;

BirdseedData::BirdseedData() {
	snpsFilename = "";
	subjectLabelsFilename = "";
	hasSubjectLabels = false;
	excludeSnpsFilename = "";
	includeSnpsFilename = "";
	hasIncludedSnps = false;
	hasExcludedSnps = false;
	phenosFilename = "";
	hasPhenotypes = false;
}

BirdseedData::~BirdseedData() {
}

bool BirdseedData::LoadData(string snpsFile, string phenoFile, string subjectsFile,
		string includeSnpsFile, string excludeSnpsFile) {

	// temporary string for reading file lines
	string line;

	/// read subjects file if specified
	if(subjectsFile != "") {
		subjectLabelsFilename = subjectsFile;
		ifstream subjectsStream(subjectLabelsFilename.c_str());
		if (!subjectsStream.is_open()) {
			cerr << "ERROR: Could not open subjects file: "
					<< subjectLabelsFilename << endl;
			return false;
		}
		cout << Timestamp() << "Reading subjects from ["
				<< subjectLabelsFilename << "]" << endl;
		int lineNumber = 0;
		while (getline(subjectsStream, line)) {
			++lineNumber;
			string trimmedLine = trim(line);
			if (!trimmedLine.size()) {
				cout << "WARNING: Blank line skipped at line number: "
						<< lineNumber << endl;
				continue;
			}
			subjectLabels.push_back(trimmedLine);
		}
		subjectsStream.close();
		hasSubjectLabels = true;
	}

	/// read SNP exclusion file
	if(excludeSnpsFile != "") {
		excludeSnpsFilename = excludeSnpsFile;
		ifstream excludeSnpsStream(excludeSnpsFilename.c_str());
		if (!excludeSnpsStream.is_open()) {
			cerr << "ERROR: Could not open subjects file: "
					<< excludeSnpsFilename << endl;
			return false;
		}
		cout << Timestamp() << "Reading excluded SNPs from ["
				<< excludeSnpsFilename << "]" << endl;
		int lineNumber = 0;
		while (getline(excludeSnpsStream, line)) {
			++lineNumber;
			string trimmedLine = trim(line);
			if (!trimmedLine.size()) {
				cout << "WARNING: Blank line skipped at line number: "
						<< lineNumber << endl;
				continue;
			}
			excludeSnps.push_back(trimmedLine);
		}
		excludeSnpsStream.close();
		hasExcludedSnps = true;
		cout << Timestamp() << excludeSnps.size()
				<< " SNPs in exclusion list" << endl;
	}

	/// read SNP inclusion file
	if(includeSnpsFile != "") {
		includeSnpsFilename = includeSnpsFile;
		ifstream includeSnpsStream(includeSnpsFilename.c_str());
		if (!includeSnpsStream.is_open()) {
			cerr << "ERROR: Could not open subjects file: "
					<< includeSnpsFilename << endl;
			return false;
		}
		cout << Timestamp() << "Reading excluded SNPs from ["
				<< includeSnpsFilename << "]" << endl;
		int lineNumber = 0;
		while (getline(includeSnpsStream, line)) {
			++lineNumber;
			string trimmedLine = trim(line);
			if (!trimmedLine.size()) {
				cout << "WARNING: Blank line skipped at line number: "
						<< lineNumber << endl;
				continue;
			}
			includeSnps.push_back(trimmedLine);
		}
		includeSnpsStream.close();
		hasIncludedSnps = true;
		cout << Timestamp() << includeSnps.size()
				<< " SNPs in inclusion list" << endl;
	}

	/// read SNPs data from the Birdseed file
	snpsFilename = snpsFile;
	ifstream genotypesStream(snpsFilename.c_str());
	if (!genotypesStream.is_open()) {
		cerr << "ERROR: Could not open SNPs file: " << snpsFilename << endl;
		return false;
	}
	cout << Timestamp() << "Reading SNPs from [" << snpsFilename << "]"
			<< endl;

	/// skip any header comment lines
	bool readingComments = true;
	while(readingComments) {
		if(getline(genotypesStream, line)) {
			if(line[0] != '#') {
				readingComments = false;
			}
		}
		else {
			cerr << "Unexpected end-of-file reading " << snpsFilename << endl;
			return false;
		}
	}

	/// assumption: past any comment rows and at the the header row
	/// 7 fields per subject
	const int COLUMNS_PER_SUBJECT = 7;
	// cout << "File header:" << endl << line << endl;
	vector<string> tokens;
	split(tokens, line, "\t");
	vector<string>::const_iterator it = tokens.begin();
	string probeId = *it;
	unsigned int numSubjects = 0;
	++it;
	for (; it != tokens.end()-1; it += COLUMNS_PER_SUBJECT) {
		string subjectName = *it;
		// cout << Timestamp() << "HEADER: Sample name: [" << subjectName << "]" << endl;
		// check that the subject name is in the subjKeys list
		// Remove quotes from subject names
		subjectName.erase(remove(subjectName.begin(), subjectName.end(), '"'), subjectName.end());
		// add to the subjectNames list
		subjectNames.push_back(subjectName);
		++numSubjects;
	}

	if(hasSubjectLabels) {
		if(subjectLabels.size() != subjectNames.size()) {
			cerr << "ERROR: Number of  subjects read from file " << subjectNames.size()
					<< " does not equal the number of labels specified in the subjects file "
					<< subjectLabelsFilename << ", which has " << subjectLabels.size()
					<< " labels" << endl;
			return false;
		}
	}
	if(subjectNames.size()) {
		cout << Timestamp() << subjectNames.size() << " subjects included" << endl;
	}
	else {
		cerr << "No subjects read from the file header" << endl;
		return false;
	}

	/// read SNP genotypes across all SNPs and all subjects
	std::vector<std::vector<std::string> > genotypes;
	cout << Timestamp()
			<< "Reading and encoding SNP data for all subjects" << endl;
	unsigned int lineNumber = 0;
	vector<string> tmpSnpNames;
	int numExcludedSnps = 0;
	int numIncludedSnps = 0;
	while (getline(genotypesStream, line)) {
		++lineNumber;
		if(lineNumber % 100000 == 0) {
			cout << Timestamp() << lineNumber << endl;
		}

		string trimmedLine = trim(line);
		// no blank lines in the data section
		if (!trimmedLine.size()) {
			cout << "WARNING: Blank line skipped at line number: "
					<< lineNumber << endl;
			continue;
		}

		/// split the line into genotypes
		vector<string> birdseedLineParts;
		split(birdseedLineParts, trimmedLine, "\t");

		/// first field is the Affymetrix SNP ID
		string snpID = birdseedLineParts[0];

		/// check for the inclusion/exclusion of this snpID
		if(hasExcludedSnps) {
			if(find(excludeSnps.begin(), excludeSnps.end(), snpID) != excludeSnps.end()) {
				/// skip this SNP
				++numExcludedSnps;
				continue;
			}
		}
		if(hasIncludedSnps) {
			if(find(includeSnps.begin(), includeSnps.end(), snpID) == includeSnps.end()) {
				/// skip this SNP
				continue;
			}
			else {
				++numIncludedSnps;
			}
		}

		tmpSnpNames.push_back(snpID);
//		cout << Timestamp() << "Splitting SNP [" << snpID
//				<< "] line into genotypes and allele counts" << endl;
		vector<string> genotypesStrings;
		vector<string>::const_iterator it =
				birdseedLineParts.begin() + 7;
		map<char, unsigned int> thisSnpAlleles;
		map<string, unsigned int> thisSnpGenotypes;
		char allele1, allele2;
		for (; it != birdseedLineParts.end()-1; it += 7) {
			string thisStringGenotype = *it;
			genotypesStrings.push_back(thisStringGenotype);
			if(thisStringGenotype == "---") {
				continue;
			}
			++thisSnpGenotypes[thisStringGenotype];
			allele1 = thisStringGenotype[0];
			allele2 = thisStringGenotype[1];
			++thisSnpAlleles[allele1];
			++thisSnpAlleles[allele2];
//			cout << "|" << thisStringGenotype << "("
//					<< allele1 << "," << allele2 << ")";
		}
//			cout << endl;
		/// save the allelic distribution for this SNP
		snpAlleleCounts.push_back(thisSnpAlleles);
		genotypeCounts.push_back(thisSnpGenotypes);

		/// save this gene's counts to the counts class member variable
		genotypes.push_back(genotypesStrings);
	}
	cout << Timestamp() << lineNumber << endl;
	genotypesStream.close();

	int numSnps = tmpSnpNames.size();
	cout << Timestamp() << "Read " << numSnps
			<< " SNPs from Birdseed file" << endl;

	/// for each SNP, map two-allele genotypes to integers using allele frequencies
	cout << Timestamp() << "Mapping genotype strings to integers" << endl;
	int monomorphs = 0;
	int newSnpIndex = 0;
	for(int snpIndex=0; snpIndex < numSnps; ++snpIndex) {
		map<char, unsigned int> thisAlleleCounts = snpAlleleCounts[snpIndex];
		map<string, unsigned int> thisGenotypeMap = genotypeCounts[snpIndex];
		if(thisGenotypeMap.size() == 1) {
			cout << Timestamp() << "WARNING: SNP " << snpNames[snpIndex]
						<< " is monomorphic - skipping" << endl;
			++monomorphs;
			continue;
		}
		snpNames.push_back(tmpSnpNames[snpIndex]);
		map<char, unsigned int>::const_iterator thisAlleleCountsIt = thisAlleleCounts.begin();
		char allele1 = thisAlleleCountsIt->first;
		int allele1Count = thisAlleleCountsIt->second;
		++thisAlleleCountsIt;
		char allele2 = thisAlleleCountsIt->first;
		int allele2Count = thisAlleleCountsIt->second;
		string majorAllele = " ";
		string minorAllele = " ";
		int majorAlleleCount = 0;
		if(allele1Count > allele2Count) {
			majorAllele[0] = allele1;
			minorAllele[0] = allele2;
			majorAlleleCount = allele1Count;
		}
		else {
			majorAllele[0] = allele2;
			minorAllele[0] = allele1;
			majorAlleleCount = allele2Count;
		}
		map<string, int> thisGenotypeStringMap;
		thisGenotypeStringMap[majorAllele + majorAllele] = 0;
		thisGenotypeStringMap[majorAllele + minorAllele] = 1;
		thisGenotypeStringMap[minorAllele + minorAllele] = 2;

		snpAlleleCounts.push_back(thisAlleleCounts);
		snpMajorMinorAlleles.push_back(make_pair(majorAllele[0], minorAllele[0]));

		/// map genotype string vector to genotype int vector
		vector<string> thisSnpGenotypes = genotypes[snpIndex];
		vector<int> thisSnpGenotypesInt;
		for(size_t sampleIndex=0; sampleIndex < thisSnpGenotypes.size(); ++sampleIndex) {
			string thisGenotypeString = thisSnpGenotypes[sampleIndex];
			if(thisGenotypeString == "---") {
				thisSnpGenotypesInt.push_back(MISSING_ATTRIBUTE_VALUE);
			}
			else {
				thisSnpGenotypesInt.push_back(thisGenotypeStringMap[thisGenotypeString]);
			}
		}
		snpGenotypes.push_back(thisSnpGenotypesInt);

		double majorAlleleFreq = majorAlleleCount / (thisSnpGenotypesInt.size() * 2);
		snpMajorAlleleFreq.push_back(majorAlleleFreq);

		++newSnpIndex;

		if(snpIndex && (snpIndex % 100000 == 0)) {
			cout << Timestamp() << snpIndex << endl;
		}
	}
	cout << Timestamp() << subjectNames.size() << " subjects read and encoded" << endl;
	cout << Timestamp() << snpGenotypes.size() << " SNPs read and encoded" << endl;
	if(monomorphs) {
		cout << Timestamp() << monomorphs << " monomorphic SNPs detected and skipped"
				<< endl;
	}

	/// read phenotypes
	if(phenoFile != "") {
		cout << Timestamp() << "Reading phenotypes from file: " << phenoFile << endl;
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

		if(phenotypes.size() != subjectNames.size()) {
			cerr << "ERROR: Number of phenotypes read: " << phenotypes.size()
					<< " is not equal to the number of sample names read from the"
					<< " counts file: " << subjectNames.size() << endl;
			return false;
		}
		hasPhenotypes = true;
	}
	else {
		for(int i=0; i < (int) subjectNames.size(); ++i) {
			phenotypes.push_back(MISSING_DISCRETE_CLASS_VALUE);
		}
		hasPhenotypes = false;
	}

	cout << Timestamp() << "Read " << subjectNames.size() << " samples with "
			<< snpGenotypes.size() << " SNPs each" << endl;

	return true;
}

vector<string> BirdseedData::GetSubjectNames() {
	return subjectNames;
}

vector<string> BirdseedData::GetSubjectLabels() {
	return subjectLabels;
}

bool BirdseedData::HasSubjectLabels() {
	return hasSubjectLabels;
}

vector<string> BirdseedData::GetSNPNames() {
	return snpNames;
}

int BirdseedData::GetNumSubjects() {
	return subjectNames.size();
}

int BirdseedData::GetNumSNPs() {
	return snpNames.size();
}

vector<int> BirdseedData::GetSubjectGenotypes(int subjectIndex) {
	if((subjectIndex < 0) || (subjectIndex >= (int) subjectNames.size())) {
		cerr << "ERROR: BirdseedData::GetSampleCounts, index out of range: "
				<< subjectIndex << endl;
		exit(EXIT_FAILURE);
	}

	vector<int> returnVector;
	for(int i=0; i < (int) snpNames.size(); ++i) {
		returnVector.push_back(snpGenotypes[i][subjectIndex]);
	}

	return returnVector;
}

int BirdseedData::GetSamplePhenotype(int subjectIndex) {
	if((subjectIndex < 0) || (subjectIndex >= (int) subjectNames.size())) {
		cerr << "ERROR: BirdseedData::GetSamplePhenotype, index out of range: "
				<< subjectIndex << endl;
		exit(EXIT_FAILURE);
	}

	return phenotypes[subjectIndex];
}

void BirdseedData::PrintInfo() {
	cout << Timestamp() << "Birdseed Statistics:" << endl;
	cout << Timestamp() << "Subjects: " << subjectNames.size( )<< endl;
	cout << Timestamp() << "SNPs:     " << snpNames.size( )<< endl;
	if(hasPhenotypes) {
		cout << Timestamp() << "Has Phenotypes" << endl;
	}
	else {
		cout << Timestamp() << "Does Not Have Phenotypes" << endl;
	}
	if(hasSubjectLabels) {
		cout << Timestamp() << "Using Subject Labels from Command-line File" << endl;
	}
	else {
		cout << Timestamp() << "Using Subject Labels from Original Data File" << endl;
	}
}

bool BirdseedData::HasPhenotypes() {
	return hasPhenotypes;
}

pair<char, char> BirdseedData::GetMajorMinorAlleles(int snpIndex) {
	pair<char, char> returnPair;
	if((snpIndex >= 0) && snpIndex < (int) snpMajorMinorAlleles.size()) {
		returnPair = snpMajorMinorAlleles[snpIndex];
	}
	else {
		cerr << "ERROR: SNP index out of range" << snpIndex << endl;
	}
	return returnPair;
}

double BirdseedData::GetMajorAlleleFrequency(int snpIndex) {
	double returnFreq = -1;
	if((snpIndex >= 0) && snpIndex < (int) snpMajorAlleleFreq.size()) {
		returnFreq = snpMajorAlleleFreq[snpIndex];
	}
	else {
		cerr << "ERROR: SNP index out of range" << snpIndex << endl;
	}
	return returnFreq;
}

map<char, unsigned int> BirdseedData::GetAlleleCounts(int snpIndex) {
	map<char, unsigned int> returnMap;
	if((snpIndex >= 0) && snpIndex < (int) snpAlleleCounts.size()) {
		returnMap = snpAlleleCounts[snpIndex];
	}
	else {
		cerr << "ERROR: SNP index out of range" << snpIndex << endl;
	}
	return returnMap;
}

map<string, unsigned int> BirdseedData::GetGenotypeCounts(int snpIndex) {
	map<string, unsigned int> returnMap;
	if((snpIndex >= 0) && snpIndex < (int) genotypeCounts.size()) {
		returnMap = genotypeCounts[snpIndex];
	}
	else {
		cerr << "ERROR: SNP index out of range" << snpIndex << endl;
	}
	return returnMap;
}
