/*
 * Edger.cpp - Bill White - 10/2/12
 * 
 * Edger algorithm implementation.
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "gsl/gsl_cdf.h"

#include "boost/lexical_cast.hpp"

#include "Edger.h"
#include "Dataset.h"
#include "Insilico.h"
#include "StringUtils.h"

using namespace std;
using namespace insilico;

Edger::Edger(Dataset* ds): AttributeRanker::AttributeRanker(ds) {
  if(ds) {
  	if(!ds->HasNumerics()) {
      cerr << "ERROR: Edger::constructor: data set must have numeric data"
      		<< endl;
      exit(-1);
  	}
    dataset = ds;
  } else {
    cerr << "ERROR: Edger::constructor: data set is NULL" << endl;
    exit(-1);
  }
}

Edger::~Edger() {
}

vector<pair<double, string> > Edger::ComputeScores() {
	cout << Timestamp() << "Running Edger through C system() call" << endl;

	/// save the current data set to a temporary file for Edger
	string tempFile = dataset->GetNumericsFilename() + "_tmp.csv";
	cout << Timestamp() << "Writing temporary file for Edger: " << tempFile
			<< endl;
	WriteDatasetInMayoFormat(tempFile);

	/// run edgeR through a system call to the shell
	stringstream edgerCmd;
	edgerCmd << "bash run_mayo_edge.sh " << tempFile << " 1";
	cout << Timestamp() << "Running edgeR command: " << edgerCmd.str() << endl;
	int systemCallReturnStatus = system(edgerCmd.str().c_str());
	if (systemCallReturnStatus == -1) {
		cerr << "ERROR: Calling R script for edgeR. -1 return code" << endl;
		exit(-1);
	}

	/// loads edgeR scores map from the output file
	// main1.csv.1.ranks.txt.sorted
	string resultsFilename = tempFile + ".1.ranks.txt.sorted";
	ReadEdgerScores(resultsFilename);

	/// remove the temporary files
	cout << Timestamp() << "Removing temporary file for edger: " << tempFile
			<< endl;
	vector<string> tempFiles;
	tempFiles.push_back(tempFile);
	//-rwxrwxr-x 1 bwhite bwhite   63850 2012-05-24 22:23 mainresults_edgeR/main1.csv.1.edger.viz.png
	//-rwxr-xr-x 1 bwhite bwhite     829 2012-03-22 20:57 mainresults_edgeR/main1.csv.1.normfactors
	//-rw-rw-r-- 1 bwhite bwhite 1064091 2012-10-01 13:38 mainresults_edgeR/main1.csv.1.ranks.txt
	//-rw-rw-r-- 1 bwhite bwhite 1064091 2012-10-01 13:38 mainresults_edgeR/main1.csv.1.ranks.txt.sorted
	//-rwxrwxr-x 1 bwhite bwhite   17439 2012-05-24 22:23 mainresults_edgeR/main1.csv.1.siggenes.dn.txt
	//-rwxrwxr-x 1 bwhite bwhite   15629 2012-05-24 22:23 mainresults_edgeR/main1.csv.1.siggenes.up.txt
	tempFiles.push_back(tempFile + ".1.edger." + "viz.png");
	tempFiles.push_back(tempFile + ".1.edger." + "normfactors");
	tempFiles.push_back(tempFile + ".1.edger." + "ranks.txt");
	tempFiles.push_back(tempFile + ".1.edger." + "ranks.txt.sorted");
	tempFiles.push_back(tempFile + ".1.edger." + "siggenes.dn.results.txt");
	tempFiles.push_back(tempFile + ".1.edger." + "siggenes.up.results.txt");
	for(unsigned int i=0; i < tempFiles.size(); ++i) {
		unlink(tempFiles[i].c_str());
	}

  return scores;
}

void Edger::PrintScores(ofstream& outFile, unsigned int topN) {
  if(topN == 0 || topN >= dataset->NumNumerics()) {
    topN = dataset->NumAttributes();
    cout << Timestamp() << "Edger values (1-pvalue) for all attributes:" << endl;
  } else {
    cout << Timestamp() << "Edger values (1-pvalue) for top [" << topN
            << "] attributes:" << endl;
  }
  vector<pair<double, string> >::const_iterator aIt = scores.begin();
  unsigned int n = 0;
  for(; aIt != scores.end() && n < topN; ++aIt, ++n) {
    outFile << (*aIt).first << "\t" << (*aIt).second << endl;
  }
}

void Edger::WriteScores(string outFilename, unsigned int topN) {
  if(topN == 0 || topN >= dataset->NumAttributes()) {
    topN = dataset->NumAttributes();
  }
  ostringstream resultsFilename;
  resultsFilename << outFilename << ".edger";
  ofstream outFile;
  outFile.open(resultsFilename.str().c_str());
  if(outFile.bad()) {
    cerr << "ERROR: Could not open scores file for writing" << endl;
    exit(-1);
  }
  PrintScores(outFile, topN);
  outFile.close();
}

void Edger::WriteDatasetInMayoFormat(string filename) {
  ofstream outFile;
  outFile.open(filename.c_str());
  if(outFile.bad()) {
    cerr << "ERROR: Edger::WriteDatasetInMayoFormat: " <<
    		"Could not open file for writing"	<< endl;
    exit(-1);
  }

  // write comma-delimited phenotypes, skipping first column
  vector<ClassLevel> classValues;
  dataset->GetClassValues(classValues);
  for(unsigned int i=0; i < classValues.size(); ++i) {
  	outFile << "," << classValues[i];
  }
  outFile << endl;

  // write one "gene" per line: counts for each individual for that gene
	vector<string> instanceIds = dataset->GetInstanceIds();
	vector<string> numericsNames =	dataset->GetNumericsNames();
	unsigned instanceIndex = 0;
	for (unsigned int nIdx = 0; nIdx < numericsNames.size(); nIdx++) {
		string thisNumericName = numericsNames[nIdx];
		outFile << thisNumericName;
		for (unsigned int iIdx = 0; iIdx < instanceIds.size(); iIdx++) {
			dataset->GetInstanceIndexForID(instanceIds[iIdx], instanceIndex);
			outFile << "," << dataset->GetNumeric(instanceIndex, thisNumericName);
		}
		outFile << endl;
	}

	outFile.close();
}

bool Edger::ReadEdgerScores(string resultsFilename) {
	ifstream resultsStream(resultsFilename.c_str());
	if (!resultsStream.is_open()) {
		cerr << "ERROR: Could not open edgeR results file: "
				<< resultsFilename << endl;
		return false;
	}
	string line;
	// strip the header line
	getline(resultsStream, line);
	// read and store variable name and gini index
	unsigned int lineNumber = 0;
	scores.clear();
	while (getline(resultsStream, line)) {
		++lineNumber;
		vector<string> tokens;
		split(tokens, line, "\t");
		// gene logConc logFC   p.value
		// 0    1       2       3
		if (tokens.size() != 4) {
			cerr << "ERROR: EvaporativeCooling::ReadEdgerScores: "
					<< "error parsing line " << lineNumber << " of "
					<< resultsFilename << ". Read " << tokens.size()
					<< " columns. Should " << "be 4" << endl;
			return false;
		}
		string gene = tokens[1];
		string genePval = tokens[3];
		double score = 1.0 - boost::lexical_cast<double>(genePval.c_str());
		// cout << "Storing Deseq: " << key << " => " << val << " form " << keyVal << endl;
		scores.push_back(make_pair(score, gene));
	}
	resultsStream.close();

	return true;
}
