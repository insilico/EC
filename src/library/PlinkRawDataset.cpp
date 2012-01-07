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
  vector<string> tokens;
  split(tokens, line);
  vector<string>::const_iterator it;
  unsigned int classIndex = 0;
  unsigned int numAttributes = 0;
  for(it = tokens.begin() + 5; it != tokens.end(); it++) {
    string headerFieldName = to_upper(*it);
    if(headerFieldName == "PHENOTYPE") {
      cout << Timestamp() << "Class column detect at " << classIndex << endl;
      classColumn = classIndex;
    } else {
      attributeNames.push_back(*it);
      attributesMask[*it] = numAttributes;
      ++numAttributes;
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

  vector<string> missingValuesToCheck;
  missingValuesToCheck.push_back("?");
  missingValuesToCheck.push_back("-9");
  missingValuesToCheck.push_back("NA");

  // read instance attributes from whitespace-delimited lines
  unsigned int instanceIndex = 0;
  unsigned int lineNumber = 0;
  cout << Timestamp();
  ValueType classType = NO_VALUE;
  double minPheno = 0.0, maxPheno = 0.0;
  while(getline(dataStream, line)) {
    ++lineNumber;
    string trimmedLine = trim(line);
    vector<string> dataLineParts;
    split(dataLineParts, trimmedLine);
    // only load those instances with matching IDs
    string ID = dataLineParts[0];
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
    ClassLevel discreteClassLevel = MISSING_DISCRETE_CLASS_VALUE;
    NumericLevel numericClassLevel = MISSING_NUMERIC_CLASS_VALUE;
    bool makeLineIntoInstance = false;
    vector<string>::const_iterator it = attributesStringVector.begin();
    for(; it != attributesStringVector.end(); ++it, ++vectorIdx) {
      string thisAttr = *it;
      if(vectorIdx == classColumn) {
        if(lineNumber == 1) {
          classType = GetClassValueType(thisAttr, missingValuesToCheck);
          switch(classType) {
            case DISCRETE_VALUE:
              hasContinuousPhenotypes = false;
              // cout << Timestamp() << "Detected DISCRETE phenotype" << endl;
              break;
            case NUMERIC_VALUE:
              hasContinuousPhenotypes = true;
              // cout << Timestamp() << "Detected NUMERIC phenotype" << endl;
              break;
            case MISSING_VALUE:
              cout << Timestamp()
                      << "WARNING: missing phenotype - skipping line: "
                      << lineNumber << endl;
              continue;
            default:
              cerr << "Could not determine class type on line 1" << endl;
              return false;
          }
        }
        discreteClassLevel = MISSING_DISCRETE_CLASS_VALUE;
        numericClassLevel = MISSING_NUMERIC_CLASS_VALUE;
        switch(classType) {
          case DISCRETE_VALUE:
            discreteClassLevel = MISSING_DISCRETE_CLASS_VALUE;
            if(!GetDiscreteClassLevel(thisAttr, missingValuesToCheck,
                                      discreteClassLevel)) {
              cerr << "ERROR: Could not get class level on line: "
                      << lineNumber << endl;
              return false;
            } else {
              if(discreteClassLevel == MISSING_DISCRETE_CLASS_VALUE) {
                cout << Timestamp()
                        << "WARNING: missing discrete phenotype skipped on line: "
                        << lineNumber << endl;
                continue;
              }
              makeLineIntoInstance = true;
            }
            break;
          case NUMERIC_VALUE:
            numericClassLevel = MISSING_NUMERIC_CLASS_VALUE;
            if(!GetNumericClassLevel(thisAttr, missingValuesToCheck,
                                     numericClassLevel)) {
              cerr << "ERROR: Could not get class level on line: "
                      << lineNumber << endl;
              return false;
            } else {
              if(numericClassLevel == MISSING_NUMERIC_CLASS_VALUE) {
                cout << Timestamp()
                        << "WARNING: missing numeric phenotype skipped on line: "
                        << lineNumber << endl;
                continue;
              }
              makeLineIntoInstance = true;
              if(lineNumber == 1) {
                minPheno = maxPheno = numericClassLevel;
              } else {
                if(numericClassLevel < minPheno) {
                  minPheno = numericClassLevel;
                }
                if(numericClassLevel > maxPheno) {
                  maxPheno = numericClassLevel;
                }
              }

            }
            break;
          case NO_VALUE:
            cerr << "ERROR: class type could not be determined on line: "
                    << lineNumber << endl;
            return false;
            break;
          case MISSING_VALUE:
            cout << "WARNING: missing phenotype - skipping line: "
                    << lineNumber << endl;
            continue;
        }
      } else {
        AttributeLevel thisAttrLevel = MISSING_ATTRIBUTE_VALUE;
        if(!GetAttributeLevel(thisAttr, missingValuesToCheck,
                              thisAttrLevel)) {
          cout << "ERROR: reading SNP on line: "
                  << lineNumber << endl;
          return false;
        }
        if(thisAttrLevel == MISSING_ATTRIBUTE_VALUE) {
          missingValues[ID].push_back(attrIdx);
        }
        attributeLevelsSeen[attrIdx].insert(thisAttr);
        attributeVector.push_back(thisAttrLevel);
        ++attrIdx;
      }
    }

    // create an instance from the vector of attribute and class values
    if(makeLineIntoInstance) {
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
        // instanceIdsToLoad.push_back(ID);
        instancesMask[ID] = instanceIndex;
      } else {
        cerr << "ERROR: loading PLINK RAW data set. "
                << "Could not create dataset instance for line number "
                << lineNumber << endl;
        return false;
      }
      ++instanceIndex;
    }

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

  UpdateAllLevelCounts();

  hasGenotypes = true;

  return true;
}

ValueType PlinkRawDataset::GetClassValueType(string value,
                                             vector<string> missingValues) {
  ValueType returnValueType = NO_VALUE;
  if(find(missingValues.begin(), missingValues.end(), value) !=
     missingValues.end()) {
    return MISSING_VALUE;
  } else {
    if((value == "1") || (value == "2")) {
      return DISCRETE_VALUE;
    } else {
      return NUMERIC_VALUE;
    }
  }
  return returnValueType;
}

bool PlinkRawDataset::GetDiscreteClassLevel(string inLevel,
                                            vector<string> missingValues,
                                            ClassLevel& outLevel) {
  if(find(missingValues.begin(), missingValues.end(), inLevel) !=
     missingValues.end()) {
    outLevel = MISSING_DISCRETE_CLASS_VALUE;
  } else {
    if((inLevel == "1") || (inLevel == "2")) {
      outLevel = lexical_cast<ClassLevel > (inLevel) - 1;
    } else {
      return false;
    }
  }

  return true;
}

bool PlinkRawDataset::GetNumericClassLevel(string inLevel,
                                           vector<string> missingValues,
                                           NumericLevel& outLevel) {
  if(find(missingValues.begin(), missingValues.end(), inLevel) !=
     missingValues.end()) {
    outLevel = MISSING_NUMERIC_CLASS_VALUE;
  } else {
    outLevel = lexical_cast<NumericLevel > (inLevel);
  }

  return true;
}
