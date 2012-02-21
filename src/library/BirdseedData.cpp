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
	vector<vector<string> > fileGenotypes;
	cout << Timestamp()
			<< "Reading and encoding SNP data for all subjects" << endl;
	unsigned int lineNumber = 0;
	vector<string> tmpSelectedSnpNames;
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
				/// found - skip this SNP
				++numExcludedSnps;
				continue;
			}
		}
		if(hasIncludedSnps) {
			if(find(includeSnps.begin(), includeSnps.end(), snpID) == includeSnps.end()) {
				/// not found - skip this SNP
				continue;
			}
			else {
				++numIncludedSnps;
			}
		}

		tmpSelectedSnpNames.push_back(snpID);
		//cout << Timestamp() << "Splitting SNP [" << snpID
		//		<< "] line into genotypes and allele counts" << endl;
		vector<string> thisSnpGenotypes;
		map<char, unsigned int> thisSnpAlleles;
		map<string, unsigned int> thisSnpGenotypeCounts;
		char allele1, allele2;
		int count = 0;
		for (unsigned int colIdx = COLUMNS_PER_SUBJECT;
				colIdx < birdseedLineParts.size();
				colIdx += COLUMNS_PER_SUBJECT) {
			string thisStringGenotype = birdseedLineParts[colIdx];
			thisSnpGenotypes.push_back(thisStringGenotype);
			++count;
			/// skip missing genotypes for allele updates
			if(thisStringGenotype != "---") {
				++thisSnpGenotypeCounts[thisStringGenotype];
				allele1 = thisStringGenotype[0];
				allele2 = thisStringGenotype[1];
				++thisSnpAlleles[allele1];
				++thisSnpAlleles[allele2];
//				if(snpID == "SNP_A-1978185") {
//					cout << colIdx << ": " << count << ": " << thisStringGenotype << "("
//							<< allele1 << "," << allele2 << ")" << endl;
//				}
			}
			else {
				/// missing data detected
				//cout << "MISSING!" << endl;
			}
		}
		// cout << endl;

		/// save the genotypic/allelic distribution for this SNP
		snpAlleleCounts.push_back(thisSnpAlleles);
		genotypeCounts.push_back(thisSnpGenotypeCounts);
		fileGenotypes.push_back(thisSnpGenotypes);
	}
	cout << Timestamp() << lineNumber << endl;
	genotypesStream.close();

	int numSnps = tmpSelectedSnpNames.size();
	if((int) snpAlleleCounts.size() != numSnps) {
		cerr << "Allele counts vector not equal to number of SNPs" << endl;
		return false;
	}
	if((int) genotypeCounts.size() != numSnps) {
		cerr << "Genotype counts vector not equal to number of SNPs" << endl;
		return false;
	}
	cout << Timestamp() << "Read " << numSnps << " SNPs from Birdseed file" << endl;

	// --------------------------------------------------------------------------

	/// for each SNP, map two-allele genotypes to integers using allele frequencies
	cout << Timestamp() << "Mapping genotype strings to integers" << endl;
	int monomorphs = 0;
	for(int snpIndex=0; snpIndex < numSnps; ++snpIndex) {
		// cout << endl << "Getting SNP info for index: " << snpIndex << endl;
		map<char, unsigned int> thisAlleleCounts = snpAlleleCounts[snpIndex];
		map<string, unsigned int> thisGenotypeMap = genotypeCounts[snpIndex];
		string majorAllele = " ";
		string minorAllele = " ";
		map<string, int> thisGenotypeStringMap;
		int majorAlleleCount = 0;
		// cout << "allele map size: " << thisAlleleCounts.size() << endl;
		// cout << "genotype map size: " << thisGenotypeMap.size() << endl;
		if(thisGenotypeMap.size() == 0) {
			// all missing SNPs
			cout << Timestamp() << "WARNING: SNP " << tmpSelectedSnpNames[snpIndex]
						<< " has all missing data - skipping" << endl;
			continue;
		}
		if(thisGenotypeMap.size() == 1) {
			cout << Timestamp() << "WARNING: SNP " << tmpSelectedSnpNames[snpIndex]
						<< " is monomorphic - keeping" << endl;
			map<char, unsigned int>::const_iterator monoIt = thisAlleleCounts.begin();
			char onlyAllele = monoIt->first;
			unsigned int onlyAlleleCount = monoIt->second;
			majorAllele[0] = onlyAllele;
			minorAllele[0] = onlyAllele;
			majorAlleleCount = onlyAlleleCount;
			++monomorphs;
			thisGenotypeStringMap[majorAllele + majorAllele] = 0;
			thisGenotypeStringMap[majorAllele + minorAllele] = 0;
			thisGenotypeStringMap[minorAllele + minorAllele] = 0;
		}
		else {
			map<char, unsigned int>::const_iterator thisAlleleCountsIt =
					thisAlleleCounts.begin();
			char allele1 = thisAlleleCountsIt->first;
			int allele1Count = thisAlleleCountsIt->second;
			++thisAlleleCountsIt;
			char allele2 = thisAlleleCountsIt->first;
			int allele2Count = thisAlleleCountsIt->second;
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
			thisGenotypeStringMap[majorAllele + majorAllele] = 0;
			thisGenotypeStringMap[majorAllele + minorAllele] = 1;
			thisGenotypeStringMap[minorAllele + majorAllele] = 1;
			thisGenotypeStringMap[minorAllele + minorAllele] = 2;
		}
		// cout << tmpSelectedSnpNames[snpIndex] << endl;
		snpNames.push_back(tmpSelectedSnpNames[snpIndex]);

		snpMajorMinorAlleles.push_back(make_pair(majorAllele[0], minorAllele[0]));
		// cout << majorAllele << "/" << minorAllele << endl;
		// cout << "mapping strings to ints through map" << endl;
		/// map genotype string vector to genotype int vector
		vector<string> thisSnpGenotypes = fileGenotypes[snpIndex];
		vector<int> thisSnpGenotypesInt;
		for(size_t sampleIndex=0; sampleIndex < thisSnpGenotypes.size(); ++sampleIndex) {
			string thisGenotypeString = thisSnpGenotypes[sampleIndex];
			if(thisGenotypeString == "---") {
				thisSnpGenotypesInt.push_back(MISSING_ATTRIBUTE_VALUE);
				missingValues[subjectNames[sampleIndex]].push_back(snpIndex);
			}
			else {
				// check to see if the genotype string is in the map
				if(thisGenotypeStringMap.find(thisGenotypeString) ==
						thisGenotypeStringMap.end()) {
					cerr << "ERROR: " << snpNames[snpIndex]
					     << ", genotype [" << thisGenotypeString
							 << "] not found in lookup map" << endl;
					exit(EXIT_FAILURE);
				}
				int thisGenotypeInt = thisGenotypeStringMap[thisGenotypeString];
				thisSnpGenotypesInt.push_back(thisGenotypeInt);
			}
		}
		snpGenotypes.push_back(thisSnpGenotypesInt);

		double majorAlleleFreq = ((double) majorAlleleCount) /
				(thisSnpGenotypesInt.size() * 2);
		// cout << "MAF: " << majorAlleleFreq << endl;
		snpMajorAlleleFreq.push_back(majorAlleleFreq);

		if(snpIndex && (snpIndex % 100000 == 0)) {
			cout << Timestamp() << snpIndex << endl;
		}
	}

	cout << Timestamp() << subjectNames.size() << " subjects read and encoded" << endl;
	cout << Timestamp() << snpGenotypes.size() << " SNPs read and encoded" << endl;
	cout << Timestamp() << missingValues.size() << " subjects had missing value(s)" << endl;
	if(monomorphs) {
		cout << Timestamp() << monomorphs << " monomorphic SNPs detected"
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

bool BirdseedData::GetMissingValues(std::string subjectName,
		std::vector<unsigned int>& missingValueIndices) {
	if(missingValues.find(subjectName) != missingValues.end()) {
		missingValueIndices = missingValues[subjectName];
	}
	return true;
}

void BirdseedData::PrintAlleleCounts() {
	for(unsigned int i=0; i < genotypeCounts.size(); ++i) {
		cout << "--------------------------------------------------" << endl;
		map<string, unsigned int>::const_iterator it = genotypeCounts[i].begin();
		for(; it != genotypeCounts[i].end(); ++it) {
			cout << it->first << " " << it->second << endl;
		}
	}
}
