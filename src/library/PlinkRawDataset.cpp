/*
 * PlinkRawDataset.cpp - Bill White - 2/24/11
 *
 * Collection class holding DatasetInstance from an Plink 
 * .raw format file
 */

#include <string>
#include <iostream>
#include <fstream>

#include <boost/lexical_cast.hpp>

#include "StringUtils.h"
#include "PlinkRawDataset.h"
#include "Insilico.h"

using namespace std;
using namespace insilico;
using namespace boost;

PlinkRawDataset::PlinkRawDataset() : Dataset::Dataset() {
}

bool PlinkRawDataset::LoadSnps(string filename) {
  snpsFilename = filename;
  ifstream dataStream(snpsFilename.c_str());
  if(!dataStream.is_open()) {
    cerr << "ERROR: Could not open dataset: " << snpsFilename << endl;
    exit(-1);
  }
  cout << Timestamp() << "Reading plink raw dataset lines from "
          << snpsFilename << endl;

  // temporary string for reading file lines
  string line;

  // read the header row - whitespace delimited attribute names
  // special attribute named "class" can be in any position
  // keep attribute names
  cout << Timestamp() << "Reading RAW file" << endl;

  getline(dataStream, line);
  string trimmedLine = trim(line);
  vector<string> tokens;
  split(tokens, trimmedLine);
  vector<string>::const_iterator it;
  unsigned int classIndex = 0;
  unsigned int numAttributes = 0;
  for(it = tokens.begin(); it != tokens.end(); it++) {
    string headerFieldName = to_upper(*it);
    if(headerFieldName == "PHENOTYPE") {
      cout << Timestamp() << "Class column detect at " << classIndex << endl;
      classColumn = classIndex;
    } else {
    	if(classIndex > 5) {
				attributeNames.push_back(*it);
				attributesMask[*it] = numAttributes;
				++numAttributes;
    	}
    }
    ++classIndex;
  }

  levelCounts.resize(numAttributes);
  levelCountsByClass.resize(numAttributes);
  attributeLevelsSeen.resize(numAttributes);

  attributeAlleleCounts.resize(numAttributes);
  attributeMinorAllele.resize(numAttributes);
  genotypeCounts.resize(numAttributes);
  attributeMutationTypes.resize(numAttributes);

  /// Detect the class type
 	bool classDetected = false;
 	switch (DetectClassType(filename, classColumn+1, true)) {
 	case CASE_CONTROL_CLASS_TYPE:
 		cout << Timestamp() << "Case-control phenotypes detected" << endl;
 		hasContinuousPhenotypes = false;
 		classDetected = true;
 		break;
 	case CONTINUOUS_CLASS_TYPE:
 		cout << Timestamp() << "Continuous phenotypes detected" << endl;
 		hasContinuousPhenotypes = true;
 		classDetected = true;
 		break;
 	case MULTI_CLASS_TYPE:
 		cout << "ERROR: more than two discrete phenotypes detected" << endl;
 		break;
 	case NO_CLASS_TYPE:
 		cout << "ERROR: phenotypes could not be detected" << endl;
 		break;
 	}
 	if (!classDetected) {
 		return false;
 	}

  // read instance attributes from whitespace-delimited lines
  unsigned int instanceIndex = 0;
  unsigned int lineNumber = 0;
  double minPheno = 0.0, maxPheno = 0.0;
  while(getline(dataStream, line)) {
    ++lineNumber;
    string trimmedLine = trim(line);
    vector<string> dataLineParts;
    split(dataLineParts, trimmedLine);
    // only load those instances with matching IDs
		// use both FID and IID so all PLINK files will work - 4/10/13
    string ID = dataLineParts[0] + dataLineParts[1];
    if(!IsLoadableInstanceID(ID)) {
      cout << Timestamp() << "WARNING: Dataset ID [" << ID << "] skipped. "
              << "Not found in numerics and/or phenotype file(s)"
              << endl;
      continue;
    }

    vector<string> attributesStringVector;
    attributesStringVector.resize(dataLineParts.size() - 5);
    copy(dataLineParts.begin() + 5, dataLineParts.end(),
         attributesStringVector.begin());

    vector<AttributeLevel> attributeVector;
    // assume genotype 0/1/2 and phenotype 0/1
    unsigned int attrIdx = 0;
    unsigned int vectorIdx = 0;
    vector<string>::const_iterator it = attributesStringVector.begin();
    ClassLevel discreteClassLevel = MISSING_DISCRETE_CLASS_VALUE;
    NumericLevel numericClassLevel = MISSING_NUMERIC_CLASS_VALUE;
    for(; it != attributesStringVector.end(); ++it, ++vectorIdx) {
      string thisAttr = *it;
      if(vectorIdx == (classColumn - 5)) {
        discreteClassLevel = MISSING_DISCRETE_CLASS_VALUE;
        numericClassLevel = MISSING_NUMERIC_CLASS_VALUE;
    		if (hasContinuousPhenotypes) {
    			if (thisAttr != "-9") {
    				numericClassLevel = lexical_cast<NumericLevel>(thisAttr);
    				if (lineNumber == 1) {
    					minPheno = maxPheno = numericClassLevel;
    				} else {
    					if (numericClassLevel < minPheno) {
    						minPheno = numericClassLevel;
    					}
    					if (numericClassLevel > maxPheno) {
    						maxPheno = numericClassLevel;
    					}
    				}
    			} else {
    				if (!hasAlternatePhenotypes) {
    					cout << Timestamp() << "Instance ID " << ID
    							<< " filtered out by missing value" << endl;
    					continue;
    				}
    			}
    		} else {
    			if (thisAttr != "-9") {
    				discreteClassLevel = lexical_cast<ClassLevel>(thisAttr) - 1;
    			} else {
    				if (!hasAlternatePhenotypes) {
    					cout << Timestamp() << "Instance ID " << ID
    							<< " filtered out by missing value" << endl;
    					continue;
    				}
    			}
    		}
     } else {
        AttributeLevel thisAttrLevel = MISSING_ATTRIBUTE_VALUE;
        if(thisAttr == "NA") {
          missingValues[ID].push_back(attrIdx);
        }
        else {
          thisAttrLevel = lexical_cast<AttributeLevel>(thisAttr);
          attributeLevelsSeen[attrIdx].insert(thisAttr);
        }
        attributeVector.push_back(thisAttrLevel);
        ++attrIdx;
      }
    }

    // create an instance from the vector of attribute and class values
		DatasetInstance * newInst = 0;
		if(attributeVector.size() != numAttributes) {
			cerr << "ERROR: Number of attributes parsed on line " << lineNumber
							<< ": " << attributesStringVector.size()
							<< " is not equal to the number of attributes "
							<< " read from the data file header: " << numAttributes
							<< endl;
			return false;
		}
		newInst = new DatasetInstance(this);
		if(newInst) {
			if(hasContinuousPhenotypes) {
				newInst->SetPredictedValueTau(numericClassLevel);
			} else {
				newInst->SetClass(discreteClassLevel);
				classIndexes[discreteClassLevel].push_back(instanceIndex);
			}
			newInst->LoadInstanceFromVector(attributeVector);
			instances.push_back(newInst);
			instanceIds.push_back(ID);
			instancesMask[ID] = instanceIndex;
		} else {
			cerr << "ERROR: loading PLINK RAW data set. "
							<< "Could not create dataset instance for line number "
							<< lineNumber << endl;
			return false;
		}
		++instanceIndex;

    // happy lights
    if(instanceIndex && ((instanceIndex % 100) == 0)) {
      cout << Timestamp() << instanceIndex << endl;
    }
  }
  cout << Timestamp() << instanceIndex << " lines read" << endl;

  dataStream.close();

  cout << Timestamp() << "There are " << NumInstances()
          << " instances in the data set" << endl;
  cout << Timestamp() << "There are " << instancesMask.size()
          << " instances in the instance mask" << endl;
  if(instancesMask.size() == 0) {
    cerr << "ERROR: no instances in the instance mask" << endl;
    return false;
  }

  if(hasContinuousPhenotypes) {
    continuousPhenotypeMinMax = make_pair(minPheno, maxPheno);
    cout << Timestamp() << "Continuous phenotypes." << endl;
  } else {
    cout << Timestamp() << "There are " << classIndexes.size()
            << " classes in the data set" << endl;
  }

  hasGenotypes = true;
  UpdateAllLevelCounts();

  CreateDummyAlleles();
  hasAllelicInfo = true;

  return true;
}
