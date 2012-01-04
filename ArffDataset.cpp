/**
 * \class ArffDataset
 *
 * \brief * Collection class holding DatasetInstance from an ARFF format file.
 * http://www.cs.waikato.ac.nz/ml/weka/arff.html
 *
 * \sa Dataset
 *
 * \author Bill White
 * \version 1.0
 *
 * Contact: bill.c.white@gmail.com
 * Created on: 2/24/11
 */

#include <string>
#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <sstream>

#include <boost/lexical_cast.hpp>

#include "Dataset.h"
#include "DatasetInstance.h"
#include "StringUtils.h"
#include "ArffDataset.h"
#include "Insilico.h"

using namespace std;
using namespace insilico;
using namespace boost;

ArffDataset::ArffDataset() {
  ::Dataset();
  missingAttributeValuesToCheck.push_back("?");
  missingClassValuesToCheck.push_back("?");
}

bool ArffDataset::LoadSnps(string filename) {
  snpsFilename = filename;
  ifstream dataStream(snpsFilename.c_str());
  if(!dataStream.is_open()) {
    cerr << "ERROR: Could not open dataset: " << snpsFilename << endl;
    return false;
  }
  cout << Timestamp() << "ArffDataset: Reading lines from "
          << snpsFilename << endl;
  string line;
  int firstSpace = -1, secondSpace = -1;
  string attributeName = "";
  int attributeIndex = 0;
  string attributeType = "";
  string classTypeString = "";
  ValueType classType = NO_VALUE;
  unsigned int lineNumber = 0;
  double minPheno = 0.0, maxPheno = 0.0;
  while(getline(dataStream, line)) {
    ++lineNumber;
    string trimmedLine = trim(line);
    // skip blank lines
    if(trimmedLine.size() < 1) {
      continue;
    }
    // check the first character of the line
    switch(trimmedLine.at(0)) {
      case '%':
        // skip comment lines
        continue;
      case '@':
        // relation, attribute or data
        firstSpace = trimmedLine.find(" ");
        string keyword = to_upper(trimmedLine.substr(1, firstSpace - 1));
        // cout << "keyword => " << keyword << endl;
        if(keyword == "RELATION") {
          relationName = keyword;
        }
        if(keyword == "ATTRIBUTE") {
          secondSpace = trimmedLine.find(" ", firstSpace + 1);
          attributeName = trimmedLine.substr(firstSpace + 1,
                                             secondSpace - firstSpace - 1);
          // cout << "\tname: [" << attributeName << "]" << endl;
          if(to_upper(attributeName) == "CLASS") {
            classColumn = attributeIndex;
            classTypeString = to_upper(trimmedLine.substr(secondSpace + 1));
            if(attributeType == "NUMERIC") {
              hasContinuousPhenotypes = true;
              classType = NUMERIC_VALUE;
            } else {
              hasContinuousPhenotypes = false;
              classType = DISCRETE_VALUE;
            }
          } else {
            attributeNames.push_back(attributeName);
            attributeType = to_upper(trimmedLine.substr(secondSpace + 1));
            if(attributeType == "NUMERIC") {
              cerr << "ERROR: NUMERIC attributes are not yet supported" << endl;
              return false;
              attributeTypes.push_back(ARFF_NUMERIC_TYPE);
            }
            if(attributeType == "STRING") {
              cerr << "ERROR: STRING attributes are not yet supported" << endl;
              return false;
              attributeTypes.push_back(ARFF_STRING_TYPE);
            }
            if(attributeType == "DATE") {
              cerr << "ERROR: DATE attributes are not yet supported" << endl;
              return false;
              attributeTypes.push_back(ARFF_DATE_TYPE);
            }
            // must be nominal type - add nominal values to map
            attributeTypes.push_back(ARFF_NOMINAL_TYPE);
            vector<string> tokens;
            split(tokens, trimmedLine, "{");
            vector<string>::const_iterator it = tokens.end() - 1;
            string nominalsListWithCurly = *it;
            string nominalsList = trim(nominalsListWithCurly.
                                       substr(0, nominalsListWithCurly.size() - 1));
            vector<string> nominals;
            split(nominals, nominalsList, ",");
            // GENETICS CHECK HERE for plink recodeA encoding
            if(nominals.size() != 3) {
              cerr << "ERROR: This dataset is currently unsupported. SNP data "
                      << "must be encoded with {0, 1, 2} for {homozygous1, "
                      << "heterzygote, homozygous2} respectively. The following "
                      << "attributes were read successfully" << endl;
              PrintNominalsMapping();
              return false;
            }
            // GENETICS CHECK HERE
            if((nominals[0] == "0") &&
               (nominals[1] == "1") &&
               (nominals[2] == "2")) {
              nominalValues[attributeName] = nominals;
              attributesMask[attributeName] = attributeIndex;
              ++attributeIndex;
            } else {
              cerr << "ERROR: This dataset is currently unsupported. SNP data "
                      << "must be encoded with {0, 1, 2} for {homozygous1, "
                      << "heterzygote, homozygous2} respectively. The following "
                      << "attributes were read successfully" << endl;
              PrintNominalsMapping();
              return false;
            }
          }
        }
        // the rest of the file is instances
        if(keyword == "DATA") {
          if(attributeNames.size() != attributeTypes.size()) {
            cerr << "ERROR: The number of attribute names: "
                    << attributeNames.size()
                    << " is not equal to the number ot attribute types: "
                    << attributeTypes.size() << endl;
            exit(0);
          }
          unsigned int numAttributes = attributeNames.size();

          levelCounts.resize(numAttributes);
          levelCountsByClass.resize(numAttributes);
          attributeLevelsSeen.resize(numAttributes);

          attributeAlleleCounts.resize(numAttributes);
          attributeMinorAllele.resize(numAttributes);
          genotypeCounts.resize(numAttributes);
          attributeMutationTypes.resize(numAttributes);

          lineNumber = 0;
          bool makeLineIntoInstance = true;
          unsigned int instanceIndex = 0;
          cout << Timestamp();

          while(getline(dataStream, line)) {
            ++lineNumber;
            string trimmedLine = trim(line);
            // skip blank lines in the data section (usually end of file)
            if(!trimmedLine.size()) {
              continue;
            }
            // only load matching IDs, line numbers for non-plink files
            ostringstream ssLineNum;
            ssLineNum << zeroPadNumber(lineNumber, 8);
            string ID = ssLineNum.str();
            // filter out IDs
            if(!IsLoadableInstanceID(ID)) {
              cout << Timestamp() << "WARNING: "
                      << "Dataset instance ID [" << ID << "] skipped. "
                      << "Not found in list of loadable IDs. Numerics and/or "
                      << "phenotype file(s) matching filtered out this ID"
                      << endl;
              continue;
            }

            vector<string> attributesStringVector;
            split(attributesStringVector, trimmedLine, ",");
            vector<AttributeLevel> attributeVector;
            unsigned int attrIdx = 0;
            unsigned int vectorIdx = 0;
            makeLineIntoInstance = true;
            vector<string>::const_iterator it = attributesStringVector.begin();
            ClassLevel discreteClassLevel = MISSING_DISCRETE_CLASS_VALUE;
            NumericLevel numericClassLevel = MISSING_NUMERIC_CLASS_VALUE;
            for(; it != attributesStringVector.end(); ++it, ++vectorIdx) {
              string thisAttr = *it;
              if(vectorIdx == classColumn) {
                string thisClassString = *it;
                if(lineNumber == 1) {
                  classType = GetClassValueType(thisClassString,
                                                missingClassValuesToCheck);
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
                /// assign class level
                discreteClassLevel = MISSING_DISCRETE_CLASS_VALUE;
                numericClassLevel = MISSING_NUMERIC_CLASS_VALUE;
                switch(classType) {
                  case DISCRETE_VALUE:
                    discreteClassLevel = MISSING_DISCRETE_CLASS_VALUE;
                    if(!GetDiscreteClassLevel(thisClassString,
                                              missingClassValuesToCheck,
                                              discreteClassLevel)) {
                      cerr << "ERROR: Could not get discrete class level on line: "
                              << lineNumber << endl;
                      return false;
                    } else {
                      if(discreteClassLevel == MISSING_DISCRETE_CLASS_VALUE) {
                        cout << Timestamp()
                                << "WARNING: missing phenotype skipped on line: "
                                << lineNumber << endl;
                        makeLineIntoInstance = false;
                        continue;
                      }
                    }
                    break;
                  case NUMERIC_VALUE:
                    numericClassLevel = MISSING_NUMERIC_CLASS_VALUE;
                    if(!GetNumericClassLevel(thisClassString,
                                             missingClassValuesToCheck,
                                             numericClassLevel)) {
                      cerr << "ERROR: Could not get numeric class level on line: "
                              << lineNumber << endl;
                      return false;
                    } else {
                      if(numericClassLevel == MISSING_NUMERIC_CLASS_VALUE) {
                        cout << Timestamp()
                                << "WARNING: missing phenotype skipped on line: "
                                << lineNumber << endl;
                        makeLineIntoInstance = false;
                        continue;
                      }
                      numericClassLevel = numericClassLevel;
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
                if(!GetAttributeLevel(thisAttr, missingAttributeValuesToCheck,
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
                cerr << "ERROR: loading ARFF @data section. "
                        << "Could not create dataset instance for line number "
                        << lineNumber << endl;
                return false;
              }
              ++instanceIndex;
            }

            // happy lights
            if((lineNumber - 1) && ((lineNumber % 100) == 0)) {
              cout << lineNumber << " ";
              cout.flush();
            }
            if((lineNumber - 1) && ((lineNumber % 1000) == 0)) {
              cout << endl << Timestamp();
            }
          }
        }
        break;
    } // end switch

  }// end while
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

  return true;
}

ArffAttributeType ArffDataset::GetTypeOf(unsigned int columnIndex) {
  if(columnIndex < attributeTypes.size()) {
    return attributeTypes[columnIndex];
  } else
    return(ARFF_ERROR_TYPE);
}

void ArffDataset::PrintNominalsMapping() {
  cout << Timestamp() << "Nominals and their accepted values:" << endl;
  map<string, vector<string> >::const_iterator mit = nominalValues.begin();
  for(; mit != nominalValues.end(); ++mit) {
    cout << (*mit).first << ":";
    vector<string>::const_iterator it = (*mit).second.begin();
    for(; it != (*mit).second.end(); ++it) {
      cout << " " << *it;
    }
    cout << endl;
  }
}

bool ArffDataset::GetAttributeLevel(string inLevel,
                                    vector<string> missingValues,
                                    AttributeLevel& outLevel) {
  if(find(missingValues.begin(), missingValues.end(), inLevel) !=
     missingValues.end()) {
    outLevel = MISSING_DISCRETE_CLASS_VALUE;
  } else {
    if((inLevel == "0") || (inLevel == "1") || (inLevel == "2")) {
      outLevel = lexical_cast<AttributeLevel > (inLevel);
    } else {
      return false;
    }
  }

  return true;
}

bool ArffDataset::GetDiscreteClassLevel(string inLevel,
                                        vector<string> missingValues,
                                        ClassLevel& outLevel) {
  if(find(missingValues.begin(), missingValues.end(), inLevel) !=
     missingValues.end()) {
    outLevel = MISSING_DISCRETE_CLASS_VALUE;
  } else {
    if((inLevel == "0") || (inLevel == "1")) {
      outLevel = lexical_cast<ClassLevel > (inLevel);
    } else {
      return false;
    }
  }

  return true;
}

bool ArffDataset::GetNumericClassLevel(string inLevel,
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
