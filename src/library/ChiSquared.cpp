/*
 * ChiSquared.cpp - Bill White - 6/15/05
 * 
 * ChiSquared algorithm implementation.
 * Reworked for McKinney Lab 2011.
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "gsl/gsl_cdf.h"

#include "AttributeRanker.h"
#include "ChiSquared.h"
#include "Dataset.h"

using namespace std;

ChiSquared::ChiSquared(Dataset* ds) : AttributeRanker::AttributeRanker(ds) {
  if(ds) {
    dataset = ds;
  } else {
    cerr << "ERROR: ChiSquared::constructor: dataset is NULL" << endl;
    exit(-1);
  }
  numClasses = dataset->NumClasses();
}

ChiSquared::~ChiSquared() {
}

const vector<pair<double, double> >& ChiSquared::ComputeScoresWithPValues() {
  // for each attribute
  for(unsigned int curAttribute = 0;
      curAttribute < dataset->NumAttributes();
      curAttribute++) {
    scoresPvalues.push_back(ComputeScore(curAttribute));
    // PrintTables();
  } // next attribute

  return scoresPvalues;
}

pair<double, double> ChiSquared::ComputeScore(unsigned int index) {

  PrepareForAttribute(index);

  // ------------ O B S E R V E D   V A L U E S  ----------------------------
  // go through each instance in the dataset and update
  // table counts based on class/attribute values
  unsigned int gtSum = 0;
  for(unsigned int curInst = 0;
      curInst < dataset->NumInstances();
      curInst++) {
    ClassLevel curClass = dataset->GetInstance(curInst)->GetClass();
    AttributeLevel curLevel =
            (dataset->GetInstance(curInst))->GetAttribute(index);
    if(curLevel != MISSING_ATTRIBUTE_VALUE) {
      observedFreqTable[curClass][curLevel]++;
      gtSum++;
    }
  }

  // -------------  R O W / C O L   C O U N T S  ----------------------------
  // calculate column totals of frequency table values
  double colSum = 0;
  vector<double> observedColSums(numLevels, 0.0);
  //   for each class
  for(unsigned int curLevel = 0; curLevel < numLevels; curLevel++) {
    for(unsigned int curClass = 0; curClass < numClasses; curClass++) {
      observedColSums[curLevel] +=
              observedFreqTable[curClass][curLevel];
    }
    colSum += observedColSums[curLevel];
  }

  // calculate row totals of contigency table values
  double rowSum = 0;
  vector<double> observedRowSums(numLevels, 0.0);
  for(unsigned int curClass = 0; curClass < numClasses; curClass++) {
    for(unsigned int curLevel = 0; curLevel < numLevels; curLevel++) {
      observedRowSums[curClass] +=
              observedFreqTable[curClass][curLevel];
    }
    rowSum += observedRowSums[curClass];
  }

  // --------------- D O U B L E  -  C H E C K   C O U N T S ----------------
  if((gtSum != rowSum) || (gtSum != colSum)) {
    cerr << endl
            << "ERROR: Attribute index: [" << index << "] => Grand Total: "
            << gtSum << " != (rowSum, colSum) "
            << " (" << rowSum << ", " << colSum << ")" << endl;
    PrintTables();
    exit(-1);
  }

  // ------------ E X P E C T E D   V A L U E S  ----------------------------
  //		cout << endl << "Expected values" << endl;
  // for each class
  for(unsigned int curClass = 0; curClass < numClasses; curClass++) {
    // for each level
    for(unsigned int curLevel = 0; curLevel < numLevels; curLevel++) {
      // calculate expected values
      expectedContingencyTable[curClass][curLevel] =
              (observedRowSums[curClass] * observedColSums[curLevel]) /
              (double) gtSum;
      //      cout << curClass << ", " << curLevel << " => "
      //           << observedRowSums[curClass] << " * "
      //           << observedColSums[curLevel] << " / "
      //           << gtSum << " => "
      //           << expectedContingencyTable[curClass][curLevel]
      //           << endl;
    }
  }

  // ------------ C H I - S Q U A R E D   V A L U E S  ----------------------
  //   for each class
  double curChiSquaredSum = 0.0;
  for(unsigned int curClass = 0; curClass < numClasses; curClass++) {
    //     for each level
    for(unsigned int curLevel = 0; curLevel < numLevels; curLevel++) {
      //        calculate chi-squared values
      double diff =
              observedFreqTable[curClass][curLevel] -
              expectedContingencyTable[curClass][curLevel];

      double denom = (double) expectedContingencyTable[curClass][curLevel];
      chiSquaredValues[curClass][curLevel] = (denom) ? (diff * diff) / denom : 0.0;

      curChiSquaredSum += chiSquaredValues[curClass][curLevel];
    }
  }

  // chi-squared p-value
  double degreesOfFreedom = (numClasses - 1) * (numLevels - 1);
  double pValue = gsl_cdf_chisq_Q(curChiSquaredSum, degreesOfFreedom);
  //		cout << "Attribute: " << curAttribute << " " << curChiSquaredSum << " "
  //            << p_value << endl;
  return make_pair(curChiSquaredSum, pValue);
}

void ChiSquared::PrintTables() {
  cout << endl << "------ CALCULATION TABLES ------" << endl;
  cout << endl << "Observed Frequencies Table:" << endl;
  for(unsigned int curClass = 0; curClass < numClasses; curClass++) {
    cout << "Class " << curClass << ": ";
    for(unsigned int curLevel = 0; curLevel < numLevels; curLevel++) {
      cout << observedFreqTable[curClass][curLevel] << " ";
    }
    cout << endl;
  }

  cout << endl << "Expected Contingency Table:" << endl;
  for(unsigned int curClass = 0; curClass < numClasses; curClass++) {
    cout << "Class " << curClass << ": ";
    for(unsigned int curLevel = 0; curLevel < numLevels; curLevel++) {
      cout << expectedContingencyTable[curClass][curLevel] << " ";
    }
    cout << endl;
  }

  cout << endl << "Chi_Squared Values Table:" << endl;
  for(unsigned int curClass = 0; curClass < numClasses; curClass++) {
    cout << "Class " << curClass << ": ";
    for(unsigned int curLevel = 0; curLevel < numLevels; curLevel++) {
      cout << chiSquaredValues[curClass][curLevel] << " ";
    }
    cout << endl;
  }
  cout << endl;
}

void ChiSquared::PrintScoresWithPValues(ofstream& outFile, unsigned int topN) {
  if(topN == 0 || topN >= dataset->NumAttributes()) {
    topN = dataset->NumAttributes();
    outFile << "Chi-squared values/p-values for all attributes:" << endl;
  } else {
    outFile << "Chi-squared values/p-values for top [" << topN
            << "] attributes:" << endl;
  }
  vector<string> attributeNames = dataset->GetAttributeNames();
  vector<pair<double, double> >::const_iterator aIt = scoresPvalues.begin();
  unsigned int n = 0;
  for(; aIt != scoresPvalues.end() && n < topN; ++aIt, ++n) {
    outFile << attributeNames[n] << " " << (*aIt).first << " "
            << (*aIt).second << endl;
  }
}

void ChiSquared::WriteScoresWithPValues(string outFilename, unsigned int topN) {
  if(topN == 0 || topN >= dataset->NumAttributes()) {
    topN = dataset->NumAttributes();
  }
  ostringstream resultsFilename;
  resultsFilename << outFilename << ".chisquared";
  ofstream outFile;
  outFile.open(resultsFilename.str().c_str());
  if(outFile.bad()) {
    cerr << "ERROR: Could not open scores file for writing" << endl;
    exit(-1);
  }
  PrintScoresWithPValues(outFile, topN);
  outFile.close();
}

void ChiSquared::PrepareForAttribute(unsigned int attributeIndex) {
  numLevels = dataset->NumLevels(attributeIndex);
  //  cout << "Preparing for attribute [" << attributeIndex << "] "
  //          << "with [" << numLevels << "] levels" << endl;

  observedFreqTable.resize(numClasses);
  expectedContingencyTable.resize(numClasses);
  chiSquaredValues.resize(numClasses);

  for(unsigned i = 0; i < numClasses; i++) {
    observedFreqTable[i].clear();
    observedFreqTable[i].resize(numLevels, 0.0);

    expectedContingencyTable[i].clear();
    expectedContingencyTable[i].resize(numLevels, 0.0);

    chiSquaredValues[i].clear();
    chiSquaredValues[i].resize(numLevels, 0.0);
  }
}

void ChiSquared::ClearTables() {
  //   for each class
  for(unsigned int curClass = 0; curClass < numClasses; curClass++) {
    //     for each level
    for(unsigned int curLevel = 0; curLevel < numLevels; curLevel++) {
      //       initialize count for observed values
      observedFreqTable[curClass][curLevel] = 0;
      expectedContingencyTable[curClass][curLevel] = 0;
      chiSquaredValues[curClass][curLevel] = 0.0;
    }
  }
}

AttributeScores ChiSquared::ComputeScores() {
	ComputeScoresWithPValues();
  vector<string> attributeNames = dataset->GetAttributeNames();
  vector<pair<double, double> >::const_iterator aIt = scoresPvalues.begin();
  scores.clear();
  for(unsigned int i =0; aIt != scoresPvalues.end(); ++aIt, ++i) {
  	scores.push_back(make_pair(aIt->first, attributeNames[i]));
  }
	return scores;
}
