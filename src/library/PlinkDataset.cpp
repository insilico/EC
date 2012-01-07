/*
 * PlinkDataset.cpp - Bill White - 2/24/11
 *
 * Collection class holding DatasetInstance from Plink format files
 * .map and .ped pairs
 */

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <time.h>

#include <boost/lexical_cast.hpp>

#include "StringUtils.h"
#include "FilesystemUtils.h"
#include "PlinkDataset.h"
#include "Insilico.h"

using namespace std;
using namespace insilico;
using boost::lexical_cast;

PlinkDataset::PlinkDataset() : Dataset::Dataset() {
  missingClassValuesToCheck.push_back("0");
  missingClassValuesToCheck.push_back("-9");
}

bool PlinkDataset::LoadSnps(string filename) {
  snpsFilename = filename;

  cout << Timestamp() << "PlinkDataset loading" << endl;
  filenameBase = GetFileBasename(snpsFilename);
  cout << Timestamp() << "Plink filename prefix for map and ped files: "
          << filenameBase << endl;

  // temporary string for reading file lines
  string line;

  /// read attribute information from the map file
  string mapFilename = filenameBase + ".map";
  ifstream mapDataStream(mapFilename.c_str());
  if(!mapDataStream.is_open()) {
    cerr << "ERROR: Could not open plink map file: " << mapFilename << endl;
    return false;
  }
  cout << Timestamp() << "Reading plink map/attribute metadata from "
          << mapFilename << endl;
  unsigned int mapLineNumber = 0;
  MapFileType mapFileType = ERROR_FILE;
  map<MapFileType, string> mapFileString;
  mapFileString[ERROR_FILE] = "ERROR: Unindentified file type";
  mapFileString[MAP3_FILE] = "MAP3 File Tye";
  mapFileString[MAP4_FILE] = "MAP4 File Type";
  unsigned int attrIdx = 0;
  while(getline(mapDataStream, line)) {
    ++mapLineNumber;
    string trimmedLine = trim(line);
    if(trimmedLine[0] == '#') {
      continue;
    }
    if(trimmedLine.size() == 0) {
      continue;
    }
    vector<string> tokens;
    split(tokens, trimmedLine);
    if(tokens.size() == 3) {
      mapFileType = MAP3_FILE;
      classColumn = 2;
    } else {
      if(tokens.size() == 4) {
        mapFileType = MAP4_FILE;
        classColumn = 5;
      } else {
        cerr << "ERROR: reading plink map file line "
                << mapLineNumber << ". "
                << "Each row should have three (map3) or four "
                << "(standard map) whitespace-separated columns"
                << endl;
        return false;
      }
    }
    attributeNames.push_back(tokens[1]);
    attributesMask[tokens[1]] = attrIdx;
    ++attrIdx;
  }
  mapDataStream.close();
  cout << Timestamp() << mapFileString[mapFileType] << endl;

  unsigned int numAttributes = attributeNames.size();
  cout << Timestamp() << "Setting up attribute metadata structures for "
          << "[" << numAttributes << "] attributes" << endl;
  vector<map<string, AttributeLevel> > attributeStringToInt;
  attributeStringToInt.resize(numAttributes);
  levelCounts.resize(numAttributes);
  levelCountsByClass.resize(numAttributes);
  attributeLevelsSeen.resize(numAttributes);

  attributeAlleles.resize(numAttributes);
  attributeAlleleCounts.resize(numAttributes);
  attributeMinorAllele.resize(numAttributes);
  genotypeCounts.resize(numAttributes);
  attributeMutationTypes.resize(numAttributes);
  
  vector<vector<string> > genotypeMatrix;
  /// read attribute values from the ped file
  string pedFilename = filenameBase + ".ped";
  ifstream pedDataStream(pedFilename.c_str());
  if(!pedDataStream.is_open()) {
    cerr << "ERROR: Could not open plink ped file: " << pedFilename << endl;
    return false;
  }
  cout << Timestamp() << "Reading plink attribute values from "
          << pedFilename << endl;
  unsigned int pedLineNumber = 0;
  unsigned int instanceIndex = 0;
  cout << Timestamp();
  ValueType classType = NO_VALUE;
  vector<pair<string, string> > alleles;
  double minPheno = 0.0, maxPheno = 0.0;
  while(getline(pedDataStream, line)) {
    ++pedLineNumber;
    string trimmedLine = trim(line);
    if(trimmedLine[0] == '#') {
      continue;
    }
    if(trimmedLine.size() == 0) {
      continue;
    }
    vector<string> pedColumnsParsed;
    split(pedColumnsParsed, trimmedLine);
    /// determine the MAP file type
    if(mapFileType == MAP4_FILE) {
      if(pedColumnsParsed.size() != numAttributes * 2 + 6) {
        cerr << "ERROR: readling line " << pedLineNumber << " from the ped file"
                << endl
                << pedColumnsParsed.size() << " columns read, "
                << (numAttributes + 6) << " expected" << endl;
        exit(1);
      }
    } else {
      if(pedColumnsParsed.size() != numAttributes * 2 + 2) {
        cerr << "ERROR: readling line " << pedLineNumber << " from the ped file"
                << endl
                << pedColumnsParsed.size() << " columns read, "
                << (numAttributes + 1) << " expected" << endl;
        return false;
      }
    }

    /// get ID for matching between PLINK data, numeric and pheno files
    string ID = pedColumnsParsed[0];
    if(!IsLoadableInstanceID(ID)) {
      cout << Timestamp() << "WARNING: Dataset instance ID [" << ID << "] skipped. "
              << "Not found in numerics and/or phenotype file(s)"
              << endl;
      continue;
    }

    string thisClassString = pedColumnsParsed[classColumn];
    /// determine class data type
    if(pedLineNumber == 1) {
      classType = GetClassValueType(thisClassString, missingClassValuesToCheck);
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
                  << pedLineNumber << endl;
          continue;
        default:
          cerr << "Could not determine class type on line 1" << endl;
          return false;
      }
    }
    /// assign class level
    ClassLevel discreteClassLevel = MISSING_DISCRETE_CLASS_VALUE;
    NumericLevel numericClassLevel = MISSING_NUMERIC_CLASS_VALUE;
    switch(classType) {
      case DISCRETE_VALUE:
        discreteClassLevel = MISSING_DISCRETE_CLASS_VALUE;
        if(!GetDiscreteClassLevel(thisClassString, missingClassValuesToCheck,
                                  discreteClassLevel)) {
          cerr << "ERROR: Could not get class level on line: "
                  << pedLineNumber << endl;
          return false;
        } else {
          if(discreteClassLevel == MISSING_DISCRETE_CLASS_VALUE) {
            cout << Timestamp()
                    << "WARNING: missing phenotype skipped on line: "
                    << pedLineNumber << endl;
            continue;
          }
        }
        break;
      case NUMERIC_VALUE:
         numericClassLevel = MISSING_NUMERIC_CLASS_VALUE;
        if(!GetNumericClassLevel(thisClassString, missingClassValuesToCheck,
                                 numericClassLevel)) {
          cerr << "ERROR: Could not get class level on line: "
                  << pedLineNumber << endl;
          return false;
        } else {
          if(numericClassLevel == MISSING_NUMERIC_CLASS_VALUE) {
            cout << Timestamp()
                    << "WARNING: missing phenotype skipped on line: "
                    << pedLineNumber << endl;
            continue;
          }
          numericClassLevel = numericClassLevel;
          if(pedLineNumber == 1) {
            minPheno = maxPheno = numericClassLevel;
          }
          else {
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
                << pedLineNumber << endl;
        return false;
        break;
      case MISSING_VALUE:
        cout << "WARNING: missing phenotype - skipping line: "
                << pedLineNumber << endl;
        continue;
    }

    DatasetInstance* newInst = new DatasetInstance(this);
    if(hasContinuousPhenotypes) {
      newInst->SetPredictedValueTau(numericClassLevel);
    } else {
      newInst->SetClass(discreteClassLevel);
      classIndexes[discreteClassLevel].push_back(instanceIndex);
    }
    instances.push_back(newInst);
    instanceIds.push_back(ID);
    // instanceIdsToLoad.push_back(ID);
    instancesMask[ID] = instanceIndex;

    // the remaining columns in the line are gentoypes for the instance/subject
    // as allele duets separated by spaces
    vector<string> attributeDuets;
    attributeDuets.resize(numAttributes * 2);
    if(mapFileType == MAP4_FILE) {
      copy(pedColumnsParsed.begin() + 6, pedColumnsParsed.end(),
           attributeDuets.begin());
    } else {
      copy(pedColumnsParsed.begin() + 2, pedColumnsParsed.end(),
           attributeDuets.begin());
    }

    // collapse the allelic encodings into genotypes
    string littleBuff = "  ";
    unsigned int attrIdx = 0;
    vector<string> attributesStringVector;
    vector<string>::const_iterator it = attributeDuets.begin();
    // cout << "Parsing duets..." << endl;
    for(unsigned int i = 0; it != attributeDuets.end(); ++i, ++it) {
      // cout << i << " => " << *it << endl;
      if(i % 2 == 0) {
        littleBuff[0] = (*it)[0];
      } else {
        littleBuff[1] = (*it)[0];
        string uppercaseGenotype = to_upper(littleBuff);
        if(uppercaseGenotype != "00") {
          ++attributeAlleleCounts[attrIdx][uppercaseGenotype[0]];
          ++attributeAlleleCounts[attrIdx][uppercaseGenotype[1]];
        }
//        cout << "\t" << attrIdx << ": " << littleBuff
//                << " (" << uppercaseGenotype << ")" << endl;
        ++genotypeCounts[attrIdx][uppercaseGenotype];
        attributesStringVector.push_back(uppercaseGenotype);
        ++attrIdx;
      }
    }
    genotypeMatrix.push_back(attributesStringVector);

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

  pedDataStream.close();

  // make genotype to integer map by determining the minor allele 
  // for each attribute
  attributeStringToInt.resize(numAttributes);
  for(attrIdx = 0; attrIdx < numAttributes; ++attrIdx) {
    map<char, unsigned int> thisAttrMap = attributeAlleleCounts[attrIdx];
    if(thisAttrMap.size() != 2) {
      cerr << "ERROR: Only biallelic genotypes are supported" << endl;
      return false;
    }
    map<char, unsigned int>::const_iterator mapIt = thisAttrMap.begin();
    char allele1 = mapIt->first;
    unsigned int allele1Count = mapIt->second;
    ++mapIt;
    char allele2 = mapIt->first;
    unsigned int allele2Count = mapIt->second;
    string majorAllele = " ";
    string minorAllele = " ";
    double attributeMaf = 0.0;
    if(allele1Count <= allele2Count) {
      minorAllele[0] = allele1;
      attributeMaf = ((double) allele1Count) / (NumInstances() * 2.0);
      majorAllele[0] = allele2;
    }
    else {
      minorAllele[0] = allele2;
      attributeMaf = ((double) allele2Count) / (NumInstances() * 2.0);
      majorAllele[0] = allele1;
    }

    /// set the mutation type
    attributeMutationTypes[attrIdx] =
            attributeMutationMap[make_pair(minorAllele[0],
                                           majorAllele[0])];

//    cout << "Attribute: " << attrIdx
//            << ", A1: " << minorAllele
//            << ", A2: " << majorAllele
//            << ", MAF: " << attributeMaf
//            << endl;
    attributeMinorAllele[attrIdx] = make_pair(minorAllele[0], attributeMaf);
    
    string genotype0 = majorAllele + majorAllele;
    string genotype1 = minorAllele + majorAllele;
    // all heterozygotes regardless of allele order are coded 1
    string genotype11 = majorAllele + minorAllele;
    string genotype2 = minorAllele + minorAllele;

//    cout << "0: " << genotype0
//            << " 1: " << genotype1
//            << " 2: " << genotype2 << endl;

    attributeStringToInt[attrIdx].insert(make_pair(genotype0, 0));
    attributeStringToInt[attrIdx].insert(make_pair(genotype1, 1));
    attributeStringToInt[attrIdx].insert(make_pair(genotype11, 1));
    attributeStringToInt[attrIdx].insert(make_pair(genotype2, 2));
    levelCounts[attrIdx][0] = 0;
    levelCounts[attrIdx][1] = 0;
    levelCounts[attrIdx][2] = 0;
  }

  // map all genotypes to integers to populate the data set
  // PrintStringToIntMap();
  for(instanceIndex = 0; instanceIndex < instances.size(); ++instanceIndex) {
    for(unsigned int attributeIndex = 0;
        attributeIndex < genotypeMatrix[instanceIndex].size();
        ++attributeIndex) {
      string thisAttr = genotypeMatrix[instanceIndex][attributeIndex];
      AttributeLevel thisAttrLevel;
      if(thisAttr == "00") {
        thisAttrLevel = MISSING_ATTRIBUTE_VALUE;
      } else {
        thisAttrLevel = attributeStringToInt[attributeIndex][thisAttr];
//        cout << "(" << instanceIndex << "," << attributeIndex << ") -> "
//                << thisAttr << " (" << thisAttrLevel << ")" << endl;
        attributeLevelsSeen[attributeIndex].insert(thisAttr);
      }
      instances[instanceIndex]->attributes.push_back(thisAttrLevel);
    }
  }

  cout << Timestamp() << "There are " << NumInstances()
          << " instances in the data set" << endl;
  cout << Timestamp() << "There are " << instancesMask.size()
          << " instances in the instance mask" << endl;
  if(hasContinuousPhenotypes) {
    continuousPhenotypeMinMax = make_pair(minPheno, maxPheno);
    cout << Timestamp() << "Continuous phenotypes." << endl;
  } else {
    cout << Timestamp() << "There are " << classIndexes.size()
            << " classes in the data set" << endl;
  }

  UpdateAllLevelCounts();

  hasGenotypes = true;

  cout << Timestamp() << "Dataset read and transformed into integer encoding"
          << endl;

  return true;
}

ValueType PlinkDataset::GetClassValueType(string value,
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

bool PlinkDataset::GetDiscreteClassLevel(string inLevel,
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

bool PlinkDataset::GetNumericClassLevel(string inLevel,
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

pair<char, double> PlinkDataset::GetAttributeMAF(unsigned int attributeIndex) {
  pair<char, double> returnPair = make_pair(' ', 0.0);
  if(attributeIndex < NumAttributes()) {
    returnPair = attributeMinorAllele[attributeIndex];
  }
  return returnPair;
}

AttributeMutationType
PlinkDataset::GetAttributeMutationType(unsigned int attributeIndex) {
  AttributeMutationType returnType = UNKNOWN_MUTATION;
  if(attributeIndex < NumAttributes()) {
    returnType = attributeMutationTypes[attributeIndex];
  }
  return returnType;
}

