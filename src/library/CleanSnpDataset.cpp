/*
 * CleanSnpDataset.C - Bill White - 9/22/11
 *
 * Minimalist data set assumes all-integer attributes and phenotypes.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "Dataset.h"
#include "CleanSnpDataset.h"
#include "StringUtils.h"
#include "Insilico.h"

using namespace std;
using namespace insilico;

/*****************************************************************************
 * Method: (constructor)
 *
 * IN:  none
 * OUT: none
 *
 * Constructs a dataset object
 ****************************************************************************/
CleanSnpDataset::CleanSnpDataset() : Dataset::Dataset() {
}

/*****************************************************************************
 * Method: LoadSnps
 *
 * IN:  SNPs filename
 * OUT: success
 *
 * Load SNPs assuming all tab-delimited integers 0/1/2 and binary 0/1 class.
 ****************************************************************************/
bool CleanSnpDataset::LoadSnps(string filename) {
  snpsFilename = filename;
  ifstream dataStream(snpsFilename.c_str());
  if(!dataStream.is_open()) {
    cerr << "ERROR: Could not open SNP dataset: " << snpsFilename << endl;
    exit(-1);
  }
  cout << Timestamp() << "\tReading whitespace-delimited SNP dataset lines from "
          << snpsFilename << ":" << endl;

  // temporary string for reading file lines
  string line;

  // read the header row - whitespace delimited attribute names
  // special attribute named "class" can be in any position
  // keep attribute names
  getline(dataStream, line);
  vector<string> tokens;
  split(tokens, line);
  vector<string>::const_iterator it;
  unsigned int numAttributes = 0;
  unsigned int classIndex = 0;
  for(it = tokens.begin(); it != tokens.end(); it++) {
    if(to_upper(*it) == "CLASS") {
      classColumn = classIndex;
    } else {
      attributeNames.push_back(*it);
      attributesMask[*it] = numAttributes;
      ++numAttributes;
    }
    ++classIndex;
  }
  // cout << "Class column is: " << classColumn << endl;

  levelCounts.resize(numAttributes);
  attributeLevelsSeen.resize(numAttributes);

  // read instance attributes from whitespace-delimited lines
  unsigned int instanceIndex = 0;
  unsigned int lineNumber = 0;
  cout << Timestamp();
  while(getline(dataStream, line)) {
    ++lineNumber;
    ostringstream ssLineNum;
    ssLineNum << lineNumber;
    string ID = ssLineNum.str();
    string trimmedLine = trim(line);
    vector<string> attributesStringVector;
    split(attributesStringVector, trimmedLine);
    unsigned int numAttributesRead = attributesStringVector.size() - 1;
    if(numAttributesRead == 0 || numAttributesRead != numAttributes) {
      cerr << "ERROR: Skipping line " << lineNumber
              << " instance has " << numAttributesRead
              << " should have " << numAttributes
              << endl;
      continue;
    }

    DatasetInstance* newInst = new DatasetInstance(this);
    vector<AttributeLevel> attributesIntVector;
    vector<string>::const_iterator it = attributesStringVector.begin();
    unsigned int attrIdx = 0;
    for(; it != attributesStringVector.end(); ++it) {
      string thisAttr = *it;
      unsigned int thisAttrLevel = (AttributeLevel) atoi(thisAttr.c_str());
      attributesIntVector.push_back(thisAttrLevel);
      if(attrIdx == classColumn) {
          continue;
      }
      attributeLevelsSeen[attrIdx].insert(thisAttr);
      ++levelCounts[attrIdx][thisAttrLevel];
      ++attrIdx;
    }
    newInst->LoadInstanceFromVector(attributesIntVector);
    
    // keep a map of classes and their instance indexes
    classIndexes[newInst->GetClass()].push_back(instanceIndex);
    instances.push_back(newInst);
    instanceIds.push_back(ID);
    instancesMask[ID] = instanceIndex;

    ++instanceIndex;

    // happy lights
    if(instanceIndex && ((instanceIndex % 100) == 0)) {
      cout << instanceIndex << " ";
      cout.flush();
    }
    if(instanceIndex && ((instanceIndex % 1000) == 0)) {
      cout << endl << Timestamp();
    }
  }
  cout << endl;
  dataStream.close();

  hasGenotypes = true;
  
  return true;
}
