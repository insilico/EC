/*
 * Deseq.cpp - Bill White - 8/8/12
 * 
 * Deseq algorithm implementation.
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "gsl/gsl_cdf.h"

#include "boost/lexical_cast.hpp"

#include "Deseq.h"
#include "Dataset.h"
#include "Insilico.h"
#include "StringUtils.h"

using namespace std;
using namespace insilico;

Deseq::Deseq(Dataset* ds): AttributeRanker::AttributeRanker(ds) {
  if(ds) {
  	if(!ds->HasNumerics()) {
      cerr << "ERROR: Deseq::constructor: data set must have numeric data"
      		<< endl;
      exit(-1);
  	}
    dataset = ds;
  } else {
    cerr << "ERROR: Deseq::constructor: data set is NULL" << endl;
    exit(-1);
  }
}

Deseq::~Deseq() {
}

vector<pair<double, string> > Deseq::ComputeScores() {
	cout << Timestamp() << "Running DESeq through C system() call" << endl;

	/// save the current data set to a temporary file for DESeq
	string tempFile = dataset->GetNumericsFilename() + "_tmp.csv";
	cout << Timestamp() << "Writing temporary file for DESeq: " << tempFile
			<< endl;
	WriteDatasetInMayoFormat(tempFile);

	/// run DESeq through a system call to the shell
	stringstream deseqCmd;
	deseqCmd << "bash run_mayo_deseq.sh " << tempFile << " 1";
	cout << Timestamp() << "Running DESeq command: " << deseqCmd.str() << endl;
	int systemCallReturnStatus = system(deseqCmd.str().c_str());
	if (systemCallReturnStatus == -1) {
		cerr << "ERROR: Calling R script for DESeq. -1 return code" << endl;
		exit(-1);
	}

	/// loads DESeq scores map from the output file
	string resultsFilename = tempFile + ".1.deseq.results.txt.sorted";
	ReadDeseqScores(resultsFilename);

	/// remove the temporary file
	cout << Timestamp() << "Removing temporary file for Deseq: " << tempFile
			<< endl;
	vector<string> tempFiles;
	tempFiles.push_back(tempFile);
	tempFiles.push_back(tempFile + ".1.deseq." + "viz.png");
	tempFiles.push_back(tempFile + ".1.deseq." + "hist.png");
	tempFiles.push_back(tempFile + ".1.deseq." + "sig.results.txt");
	tempFiles.push_back(tempFile + ".1.deseq." + "results.txt");
	tempFiles.push_back(tempFile + ".1.deseq." + "sig.dn.results.txt");
	tempFiles.push_back(tempFile + ".1.deseq." + "sig.up.results.txt");
	tempFiles.push_back(tempFile + ".1.deseq." + "ncctrl.png");
	tempFiles.push_back(tempFile + ".1.deseq." + "results.txt.sorted");
	for(unsigned int i=0; i < tempFiles.size(); ++i) {
		unlink(tempFiles[i].c_str());
	}

  return scores;
}

void Deseq::PrintScores(ofstream& outFile, unsigned int topN) {
  if(topN == 0 || topN >= dataset->NumNumerics()) {
    topN = dataset->NumAttributes();
    cout << Timestamp() << "Deseq values (1-pvalue) for all attributes:" << endl;
  } else {
    cout << Timestamp() << "Deseq values (1-pvalue) for top [" << topN
            << "] attributes:" << endl;
  }
  vector<pair<double, string> >::const_iterator aIt = scores.begin();
  unsigned int n = 0;
  for(; aIt != scores.end() && n < topN; ++aIt, ++n) {
    outFile << (*aIt).first << "\t" << (*aIt).second << endl;
  }
}

void Deseq::WriteScores(string outFilename, unsigned int topN) {
  if(topN == 0 || topN >= dataset->NumAttributes()) {
    topN = dataset->NumAttributes();
  }
  ostringstream resultsFilename;
  resultsFilename << outFilename << ".deseq";
  ofstream outFile;
  outFile.open(resultsFilename.str().c_str());
  if(outFile.bad()) {
    cerr << "ERROR: Could not open scores file for writing" << endl;
    exit(-1);
  }
  PrintScores(outFile, topN);
  outFile.close();
}

void Deseq::WriteDatasetInMayoFormat(string filename) {
  ofstream outFile;
  outFile.open(filename.c_str());
  if(outFile.bad()) {
    cerr << "ERROR: Deseq::WriteDatasetInMayoFormat: " <<
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

bool Deseq::ReadDeseqScores(string resultsFilename) {
	ifstream resultsStream(resultsFilename.c_str());
	if (!resultsStream.is_open()) {
		cerr << "ERROR: Could not open DESeq results file: "
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
		split(tokens, line, " ");
		// id baseMean baseMeanA baseMeanB foldChange log2FoldChange pval padj resVarA resVarB
		// 0  1        2         3         4          5              6    7    8
		if (tokens.size() != 9) {
			cerr << "ERROR: EvaporativeCooling::ReadDeseqScores: "
					<< "error parsing line " << lineNumber << " of "
					<< resultsFilename << ". Read " << tokens.size()
					<< " columns. Should " << "be 9" << endl;
			return false;
		}
		string gene = tokens[1];
		string geneCount = tokens[7];
		double key = 1.0 - boost::lexical_cast<double>(geneCount.c_str());
		// cout << "Storing Deseq: " << key << " => " << val << " form " << keyVal << endl;
		scores.push_back(make_pair(key, gene));
	}
	resultsStream.close();

	return true;
}
