/*
 * Dataset.cpp - Bill White - 6/14/05
 *
 * Collection class holding DatasetInstances
 * Added interaction infomation week of 4/18-26/06
 * Reworked entirely for McKinney Lab work - 2/28/11
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <iterator>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <sstream>

#include <limits.h>
#include <sys/types.h>
#include <unistd.h>
#include <assert.h>
#include <time.h>

#include <boost/lexical_cast.hpp>

#include "gsl/gsl_cdf.h"
#include "GSLRandomFlat.h"

#include "ChiSquared.h"
#include "Dataset.h"
#include "DatasetInstance.h"
#include "StringUtils.h"
#include "Statistics.h"
#include "Debugging.h"
#include "Insilico.h"
#include "DgeData.h"

using namespace std;
using namespace insilico;
using namespace boost;

Dataset::Dataset() {
  /// Set private data to defaults.
  snpsFilename = "";
  numericsFilename = "";
  alternatePhenotypesFilename = "";

  hasGenotypes = false;
  hasNumerics = false;
  hasAlternatePhenotypes = false;
  hasContinuousPhenotypes = false;

  classColumn = 0;

  maskIsPushed = false;

  /// Load attribute mutation map for transitions/transversions.
  attributeMutationMap[make_pair('A', 'G')] = TRANSITION_MUTATION;
  attributeMutationMap[make_pair('G', 'A')] = TRANSITION_MUTATION;
  attributeMutationMap[make_pair('C', 'T')] = TRANSITION_MUTATION;
  attributeMutationMap[make_pair('T', 'C')] = TRANSITION_MUTATION;
  attributeMutationMap[make_pair('A', 'C')] = TRANSVERSION_MUTATION;
  attributeMutationMap[make_pair('C', 'A')] = TRANSVERSION_MUTATION;
  attributeMutationMap[make_pair('G', 'T')] = TRANSVERSION_MUTATION;
  attributeMutationMap[make_pair('T', 'G')] = TRANSVERSION_MUTATION;
  attributeMutationMap[make_pair('A', 'T')] = TRANSVERSION_MUTATION;
  attributeMutationMap[make_pair('T', 'A')] = TRANSVERSION_MUTATION;
  attributeMutationMap[make_pair('C', 'G')] = TRANSVERSION_MUTATION;
  attributeMutationMap[make_pair('G', 'C')] = TRANSVERSION_MUTATION;

  rng = NULL;
}

Dataset::~Dataset() {
  vector<DatasetInstance*>::const_iterator it;
  for(it = instances.begin(); it != instances.end(); it++) {
    if(*it) {
      delete *it;
    }
  }
  if(rng) {
    delete rng;
  }
}

bool Dataset::LoadDataset(string snpsFilename,
                          string numericsFilename,
                          string altPhenoFilename,
                          vector<string> ids) {

  // instance IDs to use when loading from data, numeric and
  // alternate phenotype files
  if(ids.size()) {
    instanceIdsToLoad.resize(ids.size());
    copy(ids.begin(), ids.end(), instanceIdsToLoad.begin());
  }

  // load SNPs
  if(snpsFilename != "") {
    if(!LoadSnps(snpsFilename)) {
      cerr << "ERROR in LoadDataset. Could not load SNPs file" << endl;
      return false;
    }
    if(instancesMask.size() == 0) {
      cerr << "ERROR: No instances for analysis" << endl;
      return false;
    }
    if(attributesMask.size() == 0) {
      cerr << "ERROR: No attributes for analysis" << endl;
      return false;
    }
    hasGenotypes = true;

    // remove the instances that don't match the covariate and/or phenotype files
    unsigned int nextInstanceIndex = 0;
    if(instanceIdsToLoad.size()) {
      cout << Timestamp() << "Finding matching IDs" << endl;
      vector<DatasetInstance*> tempInstances;
      vector<string> tempInstanceIds;
      map<string, unsigned int> tempInstancesMask;
      for(unsigned int i = 0; i < instanceIdsToLoad.size(); ++i) {
        unsigned int instanceIndex = 0;
        string ID = instanceIdsToLoad[i];
        if(!GetInstanceIndexForID(ID, instanceIndex)) {
          cerr << "ERROR: Could not find ID in data set: " << ID << endl;
          return false;
        }
        tempInstances.push_back(instances[instanceIndex]);
        tempInstanceIds.push_back(ID);
        tempInstancesMask[ID] = nextInstanceIndex;
        ++nextInstanceIndex;
      }
      instances = tempInstances;
      instanceIds = tempInstanceIds;
      instancesMask = tempInstancesMask;
    }
    cout << Timestamp() << NumInstances()
            << " instances remain after covariate/phenotype matching" << endl;
  }

  // load numerics
  if(numericsFilename != "") {
    if(!LoadNumerics(numericsFilename)) {
      cerr << "ERROR in LoadDataset. Could not load numerics file" << endl;
      return false;
    }
  }

  // load alternate phenotypes
  if(altPhenoFilename != "") {
    if(!LoadAlternatePhenotypes(altPhenoFilename)) {
      cerr << "ERROR in LoadDataset. Could not load alternate phenotypes file" << endl;
      return false;
    }
  }

  // create and seed a random number generator for random sampling
  rng = new GSLRandomFlat(getpid() * time((time_t*) 0), 0.0, NumInstances());

  return true;
}

bool Dataset::LoadDataset(DgeData* dgeData) {
  	// TODO: check for genotypes already loaded; assume DGE only for now

	numericsFilename = "DGE CLASS";

	cout << Timestamp() << "Reading numerics from " << numericsFilename << endl;

	// populate numericNames, numericsMinMax and numericsMask
	vector<string> geneNames = dgeData->GetGeneNames();
	for(int i=0; i < geneNames.size(); ++i) {
		numericsNames.push_back(geneNames[i]);
		numericsMinMax.push_back(dgeData->GetGeneMinMax(i));
		numericsMask[geneNames[i]] = i;
	}

	// load the data set instances: set the instance numerics,
	// instance IDs, instance mask and phenotype
	vector<string> sampleNames = dgeData->GetSampleNames();
	for(int instanceIndex = 0; instanceIndex < sampleNames.size(); ++instanceIndex) {
		vector<double> sampleValues = dgeData->GetSampleCounts(instanceIndex);
		DatasetInstance* dsi = new DatasetInstance(this);
		for(int numericIndex=0; numericIndex < geneNames.size(); ++numericIndex) {
			dsi->AddNumeric(sampleValues[numericIndex]);
		}
		instances.push_back(dsi);
		string ID = sampleNames[instanceIndex];
		instanceIds.push_back(ID);
		numericsIds.push_back(ID);
		instancesMask[ID] = instanceIndex;
		ClassLevel thisClass = dgeData->GetSamplePhenotype(instanceIndex);
		dsi->SetClass(thisClass);
    classIndexes[thisClass].push_back(instanceIndex);
	}

	cout << Timestamp() << "Read " << NumNumerics() << " numeric attributes"
			<< " from DGE counts data" << endl;

	hasNumerics = true;

  // create and seed a random number generator for random sampling
  rng = new GSLRandomFlat(getpid() * time((time_t*) 0), 0.0, NumInstances());

  return true;
}

bool Dataset::GetAttributeRowCol(unsigned int row, unsigned int col,
                                 AttributeLevel& attrVal) {
  unsigned long numInstances = instances.size();
  unsigned long numAttributes = instances[0]->NumAttributes();

  if((row < numInstances) && (col < numAttributes)) {
    attrVal = instances[row]->GetAttribute(col);
    return true;
  }

  return false;
}

bool Dataset::GetNumericRowCol(unsigned int row, unsigned int col,
                               NumericLevel& numVal) {
  unsigned long numInstances = instances.size();
  unsigned long numNumerics = instances[0]->NumNumerics();

  if((row < numInstances) && (col < numNumerics)) {
    numVal = instances[row]->GetNumeric(col);
    return true;
  }

  return false;
}

bool Dataset::WriteNewDataset(string newDatasetFilename,
                              OutputDatasetType outputDatasetType) {

  switch(outputDatasetType) {
    case TAB_DELIMITED_DATASET:
    case CSV_DELIMITED_DATASET:
    case ARFF_DATASET:
      break;
    default:
      cerr << "ERROR: Output data set file type not recognized." << endl;
      return false;
  }

  ofstream newDatasetStream(newDatasetFilename.c_str());
  if(!newDatasetStream.is_open()) {
    cerr << "ERROR: Could not open new dataset file: "
            << newDatasetFilename << endl;
    return false;
  }

  /// write the attribute names header
  if(outputDatasetType == ARFF_DATASET) {
    newDatasetStream << "@RELATION dataset" << endl << endl;
  }
  map<string, unsigned int>::const_iterator ait = attributesMask.begin();
  for(; ait != attributesMask.end(); ++ait) {
    switch(outputDatasetType) {
      case TAB_DELIMITED_DATASET:
        newDatasetStream << ait->first << "\t";
        break;
      case CSV_DELIMITED_DATASET:
        newDatasetStream << ait->first << ",";
        break;
      case ARFF_DATASET:
        newDatasetStream << "@ATTRIBUTE " << ait->first << " {0,1,2}" << endl;
        break;
      case NO_OUTPUT_DATASET:
      default:
        cerr << "ERROR: Unrecognized output data set type: " <<
                outputDatasetType << endl;
        return false;
    }
  }
  map<string, unsigned int>::const_iterator nit = numericsMask.begin();
  for(; nit != numericsMask.end(); ++nit) {
    switch(outputDatasetType) {
      case TAB_DELIMITED_DATASET:
        newDatasetStream << nit->first << "\t";
        break;
      case CSV_DELIMITED_DATASET:
        newDatasetStream << nit->first << ",";
        break;
      case ARFF_DATASET:
        newDatasetStream << "@ATTRIBUTE " << nit->first << " numeric" << endl;
        break;
      case NO_OUTPUT_DATASET:
      default:
        cerr << "ERROR: Unrecognized output data set type: " <<
                outputDatasetType << endl;
        return false;
    }
  }
  switch(outputDatasetType) {
    case TAB_DELIMITED_DATASET:
    case CSV_DELIMITED_DATASET:
      newDatasetStream << "Class" << endl;
      break;
    case ARFF_DATASET:
      if(hasContinuousPhenotypes) {
        newDatasetStream << "@ATTRIBUTE Class numeric" << endl;
      } else {
        newDatasetStream << "@ATTRIBUTE Class {0,1}" << endl;
      }
      break;
    case NO_OUTPUT_DATASET:
    default:
      cerr << "ERROR: Unrecognized output data set type: " <<
              outputDatasetType << endl;
      return false;
  }

  /// write the data, respecting the masked attributes, numerics
  /// and masked instances - 10/28/11
  /// write the attribute names header
  if(outputDatasetType == ARFF_DATASET) {
    newDatasetStream << endl << "@DATA" << endl;
  }
  vector<string> instanceIds = GetInstanceIds();
  for(unsigned int iIdx = 0; iIdx < NumInstances(); iIdx++) {
    unsigned instanceIndex = 0;
    GetInstanceIndexForID(instanceIds[iIdx], instanceIndex);
    // write discrete attribute values
    vector<unsigned int> attrIndices = MaskGetAttributeIndices(DISCRETE_TYPE);
    for(unsigned int aIdx = 0; aIdx < attrIndices.size(); aIdx++) {
      AttributeLevel A = instances[instanceIndex]->GetAttribute(attrIndices[aIdx]);
      switch(outputDatasetType) {
        case TAB_DELIMITED_DATASET:
          if(A == MISSING_ATTRIBUTE_VALUE) {
            newDatasetStream << "?" << "\t";
          } else {
            newDatasetStream << A << "\t";
          }
          break;
        case CSV_DELIMITED_DATASET:
        case ARFF_DATASET:
          if(A == MISSING_ATTRIBUTE_VALUE) {
            newDatasetStream << "?" << ",";
          } else {
            newDatasetStream << A << ",";
          }
          break;
        case NO_OUTPUT_DATASET:
        default:
          cerr << "ERROR: Unrecognized output data set type: " <<
                  outputDatasetType << endl;
          return false;
      }
    }
    /// write continuous attribute values
    vector<unsigned int> numIndices = MaskGetAttributeIndices(NUMERIC_TYPE);
    for(unsigned int nIdx = 0; nIdx < numIndices.size(); nIdx++) {
      NumericLevel N = instances[instanceIndex]->GetNumeric(numIndices[nIdx]);
      switch(outputDatasetType) {
        case TAB_DELIMITED_DATASET:
          if(N == MISSING_NUMERIC_VALUE) {
            newDatasetStream << "?" << "\t";
          } else {
            newDatasetStream << N << "\t";
          }
          break;
        case CSV_DELIMITED_DATASET:
        case ARFF_DATASET:
          if(N == MISSING_NUMERIC_VALUE) {
            newDatasetStream << "?" << ",";
          } else {
            newDatasetStream << N << ",";
          }
          break;
        case NO_OUTPUT_DATASET:
        default:
          cerr << "ERROR: Unrecognized output data set type: " <<
                  outputDatasetType << endl;
          return false;
      }
    }
    // class/phenotype
    ClassLevel discreteClassLevel = MISSING_DISCRETE_CLASS_VALUE;
    NumericLevel continuousClassLevel = MISSING_NUMERIC_CLASS_VALUE;
    if(hasContinuousPhenotypes) {
      continuousClassLevel = instances[instanceIndex]->GetPredictedValueTau();
      newDatasetStream << continuousClassLevel << endl;
    } else {
      discreteClassLevel = instances[instanceIndex]->GetClass();
      newDatasetStream << discreteClassLevel << endl;
    }
  }

  newDatasetStream.close();

  return true;
}

bool Dataset::ExtractAttributes(string scoresFilename,
                                unsigned int topN,
                                string newDatasetFilename) {
  // read attribute scores from file
  ifstream scoresStream(scoresFilename.c_str());
  if(!scoresStream.is_open()) {
    cerr << "ERROR: Could not open scores file: " << scoresFilename << endl;
    return false;
  }
  // build a "rank map": pair<ranker score, attribute index>
  // for all scores in the score file
  string line;
  vector<pair<double, string> > rankMap;
  vector<string> tokens;
  while(getline(scoresStream, line)) {
//    cout << line << endl;
    tokens.clear();
    split(tokens, line);
//    cout << "t0: " << tokens[0] << ", t1: " << tokens[1] << endl;
    rankMap.push_back(make_pair(lexical_cast<double>(tokens[0]), tokens[1]));
  }
  scoresStream.close();

  // sort the rank map in ascending score order
  sort(rankMap.begin(), rankMap.end());

  // print the top N sorted index
//  vector<pair<double, string> >::const_iterator rankMapIt;
//  for(rankMapIt = rankMap.end() - 1;
//      rankMapIt != (rankMap.end() - topN - 1);
//      rankMapIt--) {
//    cout << rankMapIt->first << " => " << rankMapIt->second << endl;
//  }

  // write top N ranked attributes to new file
  // honor MDR format by using header row infomation - bcw - 9/13/05
  ofstream newDatasetStream(newDatasetFilename.c_str());
  if(!newDatasetStream.is_open()) {
    cerr << "ERROR: Could not open new dataset file: "
            << newDatasetFilename << endl;
    return false;
  }

  // write MDR header - bcw - 9/13/05
  vector<pair<double, string> >::const_iterator rankHdrIt;
  for(rankHdrIt = rankMap.end() - 1;
      rankHdrIt != rankMap.end() - topN - 1;
      rankHdrIt--) {
    newDatasetStream << rankHdrIt->second << "\t";
  }
  newDatasetStream << "Class" << endl;

  // write dataset instances
  // MDR format -> tab-delimited attributes + class - bcw - 9/13/05
  // adjusted for integrated data and continuous class
  vector<DatasetInstance*>::const_iterator dsiIt;
  vector<pair<double, string> >::const_iterator rankIt;
  for(dsiIt = instances.begin(); dsiIt != instances.end(); dsiIt++) {
    for(rankIt = rankMap.end() - 1;
        rankIt != rankMap.end() - topN - 1;
        rankIt--) {
      if(MaskSearchVariableType(rankIt->second, DISCRETE_TYPE)) {
        newDatasetStream << (*dsiIt)->GetAttribute(GetAttributeIndexFromName(rankIt->second)) << "\t";
      } else {
        if(MaskSearchVariableType(rankIt->second, NUMERIC_TYPE)) {
          newDatasetStream << (*dsiIt)->GetNumeric(GetNumericIndexFromName(rankIt->second)) << "\t";
        } else {
          cerr << "ERROR: Dataset::ExtractAttributes: requested variable name ["
                  << rankIt->second << "] is not in the data set" << endl;
          return false;
        }
      }
    }
    if(hasContinuousPhenotypes) {
      newDatasetStream << (*dsiIt)->GetPredictedValueTau() << endl;
    } else {
      newDatasetStream << (*dsiIt)->GetClass() << endl;
    }
  }

  newDatasetStream.close();

  return true;
}

bool Dataset::SwapAttributes(unsigned int a1, unsigned int a2) {
  vector<DatasetInstance*>::const_iterator it;
  for(it = instances.begin(); it != instances.end(); it++) {
    (*it)->SwapAttributes(a1, a2);
  }

  return true;
}

unsigned int Dataset::NumVariables() {
  return(NumAttributes() + NumNumerics());
}

vector<string> Dataset::GetVariableNames() {
  vector<string> variableNames(attributesMask.size() + numericsMask.size());
  vector<string> attrNames = GetAttributeNames();
  vector<string> numNames = GetNumericsNames();

  copy(attrNames.begin(), attrNames.end(), variableNames.begin());
  copy(numNames.begin(), numNames.end(),
       variableNames.begin() + attrNames.size());

  return variableNames;
}

unsigned int Dataset::NumInstances() {
  return instancesMask.size();
}

DatasetInstance* Dataset::GetInstance(unsigned int index) {
  if(index < instances.size()) {
    return instances[index];
  }
  cout << "ERROR: Instance index out of range: " << index << endl;
  exit(-1);
}

DatasetInstance* Dataset::GetRandomInstance() {
  unsigned int index = (unsigned int) rng->nextRandVal();
  return instances[index];
}

vector<string> Dataset::GetInstanceIds() {
  vector<string> idsToReturn;
  map<string, unsigned int>::const_iterator it = instancesMask.begin();
  for(; it != instancesMask.end(); ++it) {
    idsToReturn.push_back(it->first);
  }
  return idsToReturn;
}

bool Dataset::GetInstanceIndexForID(string ID, unsigned int& instanceIndex) {

  if(instancesMask.find(ID) != instancesMask.end()) {
    instanceIndex = instancesMask[ID];
    return true;
  }

  cerr << "ERROR: Could not find instance ID: " << ID << endl;

  return false;
}

unsigned int Dataset::NumAttributes() {
  return attributesMask.size();
}

vector<string> Dataset::GetAttributeNames() {
  vector<string> names;
  map<string, unsigned int>::const_iterator it = attributesMask.begin();
  for(; it != attributesMask.end(); ++it) {
    names.push_back(it->first);
  }
  return names;
}

bool Dataset::GetAttributeValues(unsigned int attributeIndex,
                                 vector<AttributeLevel>& attributeValues) {
  if((attributeIndex < 0) || (attributeIndex > attributeNames.size())) {
    cerr << "ERROR: Dataset::GetAttributeValues: attribute index out of range: "
            << attributeIndex << endl;
    return false;
  }
  if(hasGenotypes) {
    attributeValues.clear();
    map<string, unsigned int>::const_iterator it;
    for(it = instancesMask.begin(); it != instancesMask.end(); it++) {
      AttributeLevel thisAttribute =
              instances[it->second]->GetAttribute(attributeIndex);
      attributeValues.push_back(thisAttribute);
    }
  } else {
    cerr << "ERROR: attempting to access SNP data when none "
            << "have been loaded" << endl;
    return false;
  }

  return true;
}

bool Dataset::GetAttributeValues(string attributeName,
                                 vector<AttributeLevel>& attributeValues) {
  if(MaskSearchVariableType(attributeName, DISCRETE_TYPE)) {
    GetAttributeValues(GetAttributeIndexFromName(attributeName), attributeValues);
  } else {
    cerr << "ERROR: Dataset::GetAttributeValues cannot get attribute values for: "
            << attributeName << ". Either doesn't exist or is excluded" << endl;
    return false;
  }
  return true;
}

std::string Dataset::GetSnpsFilename() {
  return snpsFilename;
}

bool Dataset::HasGenotypes() {
  return hasGenotypes;
}

AttributeLevel Dataset::GetAttribute(unsigned instanceIndex, string name) {
  if((instanceIndex) < 0 || (instanceIndex >= instances.size())) {
    cerr << "ERROR: Dataset::GetAttribute: instance index " << instanceIndex
            << " out of range" << endl;
    exit(1);
  }
  map<string, unsigned int>::iterator pos = attributesMask.find(name);
  if(pos != attributesMask.end()) {
    return instances[instanceIndex]->GetAttribute(pos->second);
  } else {
    cerr << "ERROR: Dataset::GetAttribute: " << name << " at instance index: "
            << instanceIndex << " not found" << endl;
    exit(1);
  }
}

pair<char, double> Dataset::GetAttributeMAF(unsigned int attributeIndex) {
  /// An Intriduction to Genetic Analysis by Griffiths, Miller, Suzuki,
  /// Lewontin and Gelbart, 2000, page 715.
  pair<char, double> returnPair = make_pair(' ', 0.0);
  double p = 0.0;
  double q = 0.0;
  if(attributeIndex < NumAttributes()) {
    map<string, unsigned int> genotypes = genotypeCounts[attributeIndex];
    if(genotypes.size() == 3) {
      double f_AA = ((double) genotypes["0"]) / NumInstances();
      double f_Aa = ((double) genotypes["1"]) / NumInstances();
      double f_aa = ((double) genotypes["2"]) / NumInstances();
      p = f_AA + 0.5 * f_Aa;
      q = f_aa + 0.5 * f_Aa;
      if(p < q) {
        returnPair = make_pair('p', p);
      } else {
        returnPair = make_pair('q', q);
      }
    }
  }
  return returnPair;
}

AttributeMutationType
Dataset::GetAttributeMutationType(unsigned int attributeIndex) {
  return UNKNOWN_MUTATION;
}

unsigned int Dataset::NumLevels(unsigned int index) {
  if((index >= 0) && (index < attributeLevelsSeen.size())) {
    return attributeLevelsSeen[index].size();
  }
  cerr << "ERROR: Attempt to access number of levels of an attribute "
          << "index that is out of range: " << index << " out of "
          << attributeLevelsSeen.size() << std::endl;
  return 0;
}

unsigned int Dataset::GetAttributeIndexFromName(string attributeName) {
  for(unsigned int i = 0; i < attributeNames.size(); i++) {
    if(attributeNames[i] == attributeName) {
      return i;
    }
  }
  return INVALID_INDEX;
}

unsigned int Dataset::NumNumerics() {
  return numericsMask.size();
}

vector<string> Dataset::GetNumericsNames() {
  vector<string> names;
  map<string, unsigned int>::const_iterator it = numericsMask.begin();
  for(; it != numericsMask.end(); ++it) {
    names.push_back(it->first);
  }
  return names;
}

pair<NumericLevel, NumericLevel>
Dataset::GetMinMaxForNumeric(unsigned int numericIdx) {
  return numericsMinMax[numericIdx];
}

double Dataset::GetMeanForNumeric(unsigned int numericIdx) {
  double sum = 0.0;
  vector<unsigned int> instanceIndicies = MaskGetInstanceIndices();
  for(unsigned int i = 0; i < instanceIndicies.size(); ++i) {
    sum += instances[instanceIndicies[i]]->numerics[numericIdx];
  }

  return(sum / ((double) instanceIndicies.size()));
}

bool Dataset::HasNumerics() {
  return hasNumerics;
}

NumericLevel Dataset::GetNumeric(unsigned instanceIndex, string name) {
  if((instanceIndex < 0) || (instanceIndex >= instances.size())) {
    cerr << "ERROR: Dataset::GetNumeric: instance index " << instanceIndex
            << " out of range" << endl;
    exit(1);
  }
  map<string, unsigned int>::iterator pos = numericsMask.find(name);
  if(pos != numericsMask.end()) {
    return instances[instanceIndex]->GetNumeric(pos->second);
  } else {
    cerr << "ERROR: Dataset::GetNumeric" << endl;
    exit(1);
  }
}

bool Dataset::GetNumericValues(string numericName,
                               vector<NumericLevel>& numericValues) {
  if(MaskSearchVariableType(numericName, NUMERIC_TYPE)) {
    GetNumericValues(GetNumericIndexFromName(numericName), numericValues);
  } else {
    cerr << "ERROR: Dataset::GetNumericValues cannot get numeric values for: "
            << numericName << ". Either doesn't exist or is excluded" << endl;
    return false;
  }
  return true;
}

std::string Dataset::GetNumericsFilename() {
  return numericsFilename;
}

unsigned int Dataset::GetNumericIndexFromName(string numericName) {
  for(unsigned int i = 0; i < numericsNames.size(); i++) {
    if(numericsNames[i] == numericName) {
      return i;
    }
  }
  return INVALID_INDEX;
}

unsigned int Dataset::NumClasses() {
  return classIndexes.size();
}

unsigned int Dataset::GetClassColumn() {
  return classColumn;
}

bool Dataset::GetClassValues(vector<ClassLevel>& classValues) {
  classValues.clear();
  map<string, unsigned int>::const_iterator it;
  for(it = instancesMask.begin(); it != instancesMask.end(); it++) {
    classValues.push_back(instances[it->second]->GetClass());
  }
  return true;
}

const std::map<ClassLevel, std::vector<unsigned int> >&
Dataset::GetClassIndexes() {
  return classIndexes;
}

bool Dataset::HasAlternatePhenotypes() {
  return hasAlternatePhenotypes;
}

string Dataset::GetAlternatePhenotypesFilename() {
  return alternatePhenotypesFilename;
}

bool Dataset::HasContinuousPhenotypes() {
  return hasContinuousPhenotypes;
}

pair<double, double> Dataset::GetMinMaxForContinuousPhenotype() {
  return continuousPhenotypeMinMax;
}

void Dataset::Print() {
  PrintStats();
  map<string, unsigned int>::const_iterator it;
  for(it = instancesMask.begin(); it != instancesMask.end(); it++) {
    instances[it->second]->Print();
  }
}

void Dataset::PrintStats() {
  unsigned int numInstances = NumInstances();
  unsigned int numClasses = NumClasses();
  unsigned int numAttributes = NumAttributes();
  unsigned int numNumerics = NumNumerics();
  unsigned int numElements =
          (numInstances * (numAttributes + numNumerics)) + numInstances;

  // TODO: add class stats
  cout << Timestamp() << "Dataset has:" << endl
          << Timestamp() << "instances:      " << numInstances << endl;
  if(hasGenotypes) {
    cout << Timestamp() << "SNPs:           " << numAttributes << endl;
  }
  if(hasNumerics) {
    cout << Timestamp() << "numerics:       " << numNumerics << endl;
  }
  if(hasContinuousPhenotypes) {
    cout << Timestamp() << "continuous phenotype, min: " <<
            continuousPhenotypeMinMax.first
            << ", max: " << continuousPhenotypeMinMax.second << endl;
  } else {
    cout << Timestamp() << "classes:        " << numClasses << endl;
    PrintClassIndexInfo();
  }
  cout << Timestamp() << "total elements: " << numElements << endl;

  //PrintIntToStringMap();
  //PrintStringToIntMap();
  //  PrintLevelCounts();
  // PrintAttributeLevels();
  PrintMissingValuesStats();
}

void Dataset::PrintNumericsStats() {
  unsigned int numInstances = NumInstances();
  unsigned int numClasses = NumClasses();
  unsigned int numAttributes = NumAttributes();
  unsigned int numNumerics = NumNumerics();
  unsigned int numElements =
          (numInstances * (numAttributes + numNumerics)) + numInstances;

  // TODO: add class stats
  cout << Timestamp() << "Dataset has:" << endl
          << Timestamp() << "instances:      " << numInstances << endl;
  if(hasGenotypes) {
    cout << Timestamp() << "SNPs:           " << numAttributes << endl;
  }
  if(hasNumerics) {
    cout << Timestamp() << "numerics:       " << numNumerics << endl;
  }
//  vector< pair<double, double> >::const_iterator minMaxIt = numericsMinMax.begin();
//  for(unsigned int i = 0; minMaxIt != numericsMinMax.end(); ++minMaxIt, ++i) {
//    cout << Timestamp()
//            << (*minMaxIt).first << " <= "
//            << numericsNames[i] << " <= "
//            << (*minMaxIt).second << endl;
//  }
  if(hasContinuousPhenotypes) {
    cout << Timestamp() << "continuous phenotype, min: "
            << continuousPhenotypeMinMax.first
            << ", max: " << continuousPhenotypeMinMax.second << endl;
  } else {
    cout << Timestamp() << "classes:        " << numClasses << endl;
    PrintClassIndexInfo();
  }
  cout << Timestamp() << "total elements: " << numElements << endl;

  PrintMissingValuesStats();
}

void Dataset::PrintStatsSimple() {
  unsigned int numInstances = NumInstances();
  unsigned int numClasses = NumClasses();
  unsigned int numAttributes = NumAttributes();
  unsigned int numNumerics = NumNumerics();
  unsigned int numElements =
          (numInstances * (numAttributes + numNumerics)) + numInstances;

  cout << Timestamp() << "Dataset has:" << endl
          << Timestamp() << "instances:      " << numInstances << endl;
  if(hasGenotypes) {
    cout << Timestamp() << "SNPs:           " << numAttributes << endl;
  }
  if(hasNumerics) {
    cout << Timestamp() << "numerics:       " << numNumerics << endl;
  }

  if(hasContinuousPhenotypes) {
    cout << Timestamp() << "continuous phenotype: " << endl;
    cout << Timestamp() << "continuous phenotype, min: " << continuousPhenotypeMinMax.first
            << ", max: " << continuousPhenotypeMinMax.second << endl;
  } else {
    cout << Timestamp() << "classes:        " << numClasses << endl;
    PrintClassIndexInfo();
  }

  cout << Timestamp() << "total elements: " << numElements << endl;
}

void Dataset::PrintClassIndexInfo() {
  cout << Timestamp() << "Data Set Class Index" << endl;
  cout << Timestamp() << "Index has [" << classIndexes.size() << "] entries:" << endl;
  map<ClassLevel, vector<unsigned int> >::const_iterator mit = classIndexes.begin();
  for(; mit != classIndexes.end(); ++mit) {
    cout << Timestamp() << (*mit).first << ": " << (*mit).second.size() << endl;
  }
}

void Dataset::PrintMissingValuesStats() {
  if(missingValues.size()) {
    cout << Timestamp() << "Missing Attributes Values Detected" << endl;
    cout << Timestamp() << "Instance # Missing" << endl;
    map<string, vector<unsigned int> >::const_iterator mit;
    for(mit = missingValues.begin(); mit != missingValues.end(); ++mit) {
      cout << Timestamp()
              << mit->first
              << setw(10) << mit->second.size() << endl;
    }
  } else {
    cout << Timestamp() << "0 missing attribute values detected" << endl;
  }
  if(missingNumericValues.size()) {
    cout << Timestamp() << "Missing Numeric Values Detected" << endl;
    cout << Timestamp() << "Instance # Missing" << endl;
    map<string, vector<unsigned int> >::const_iterator mit;
    for(mit = missingNumericValues.begin();
        mit != missingNumericValues.end(); ++mit) {
      cout << Timestamp()
              << mit->first
              << setw(10) << mit->second.size() << endl;
    }
  } else {
    cout << Timestamp() << "0 missing numeric values detected" << endl;
  }
}

void Dataset::PrintLevelCounts() {
  cout << Timestamp() << "Dataset attribute level counts:" << endl;
  vector<map<AttributeLevel, unsigned int> >::const_iterator levelCountsIt =
          levelCounts.begin();
  for(unsigned int attrIdx = 0;
      levelCountsIt != levelCounts.end();
      ++levelCountsIt, ++attrIdx) {
    cout << Timestamp() << "Attribute [" << attrIdx << "] => [ ";
    map<AttributeLevel, unsigned int>::const_iterator itsIt =
            (*levelCountsIt).begin();
    for(; itsIt != (*levelCountsIt).end(); ++itsIt) {
      cout << (*itsIt).first << "->" << (*itsIt).second << " ";
    }
    cout << " ]" << endl;
  }

}

void Dataset::WriteLevelCounts(std::string levelsFilename) {
  ofstream outFile;
  outFile.open(levelsFilename.c_str());
  if(outFile.bad()) {
    cerr << "ERROR: Could not open levels counts file for writing" << endl;
    exit(1);
  }

  vector<map<AttributeLevel, unsigned int> >::const_iterator levelCountsIt =
          levelCounts.begin();
  for(unsigned int attrIdx = 0;
      levelCountsIt != levelCounts.end();
      ++levelCountsIt, ++attrIdx) {
    map<AttributeLevel, unsigned int>::const_iterator itsIt =
            (*levelCountsIt).begin();
    outFile << attributeNames[attrIdx];
    for(; itsIt != (*levelCountsIt).end(); ++itsIt) {
      outFile << "\t" << (*itsIt).second;
    }
    outFile << endl;
  }

  outFile.close();
}

void Dataset::PrintRecodeMap(vector<map<unsigned int, unsigned int> > recodeMap) {
  cout << Timestamp() << "Recoding maps for heterozygosity in level 1" << endl;
  vector<map<unsigned int, unsigned int> >::const_iterator levelCountsIt = recodeMap.begin();
  for(unsigned int attrIdx = 0; levelCountsIt != recodeMap.end(); ++levelCountsIt, ++attrIdx) {
    cout << Timestamp() << "Attribute [" << attrIdx << "] => [ | ";
    map<unsigned int, unsigned int>::const_iterator itsIt = (*levelCountsIt).begin();
    for(; itsIt != (*levelCountsIt).end(); ++itsIt) {
      cout << (*itsIt).first << " => " << (*itsIt).second << " | ";
    }
    cout << "]" << endl;
  }

}

void Dataset::PrintAttributeLevelsSeen() {
  cout << Timestamp() << "Dataset attribute levels seen:" << endl;
  vector<set<string> >::const_iterator attributeIt = attributeLevelsSeen.begin();
  for(unsigned int i = 0; attributeIt != attributeLevelsSeen.end();
      ++attributeIt, ++i) {
    cout << Timestamp() << "Attribute [" << i << "] => [ ";
    set<string>::const_iterator sIt = (*attributeIt).begin();
    for(; sIt != (*attributeIt).end(); ++sIt) {
      cout << *sIt << " ";
    }
    cout << "]" << endl;
  }
}

bool Dataset::MaskRemoveVariable(string variableName) {
	if(MaskSearchVariableType(variableName, DISCRETE_TYPE)) {
		MaskRemoveVariableType(variableName, DISCRETE_TYPE);
		return true;
	}
	if(MaskSearchVariableType(variableName, NUMERIC_TYPE)) {
		MaskRemoveVariableType(variableName, NUMERIC_TYPE);
		return true;
	}
	return false;
}

bool Dataset::MaskRemoveVariableType(string variableName,
                                     AttributeType varType) {
  map<string, unsigned int>::iterator pos;
  if(varType == DISCRETE_TYPE) {
    pos = attributesMask.find(variableName);
    if(pos != attributesMask.end()) {
      attributesMask.erase(pos);
    } else {
      cerr << "ERROR: Dataset::MaskRemoveVariable failed for SNP attribute name: "
              << variableName << ". name not found" << endl;
      return false;
    }
  } else {
    pos = numericsMask.find(variableName);
    if(pos != numericsMask.end()) {
      numericsMask.erase(pos);
    } else {
      cerr << "ERROR: Dataset::MaskRemoveVariable failed for numerics"
              << " attribute name: " << variableName << ". name not found"
              << endl;
      return false;
    }
  }
  return true;
}

bool Dataset::MaskSearchVariableType(string variableName,
                                 AttributeType varType) {
  map<string, unsigned int>::iterator pos;
  if(varType == DISCRETE_TYPE) {
    pos = attributesMask.find(variableName);
    if(pos != attributesMask.end()) {
      return true;
    } else {
      return false;
    }
  } else {
    pos = numericsMask.find(variableName);
    if(pos != numericsMask.end()) {
      return true;
    } else {
      return false;
    }
  }
}

bool Dataset::MaskIncludeAllAttributes(AttributeType attrType) {
  if(attrType == DISCRETE_TYPE) {
    attributesMask.clear();
    unsigned int attributeIndex = 0;
    vector<string>::const_iterator it = attributeNames.begin();
    for(; it != attributeNames.end(); ++it) {
      attributesMask[*it] = attributeIndex;
      ++attributeIndex;
    }
    return true;
  } else {
    numericsMask.clear();
    unsigned int attributeIndex = 0;
    vector<string>::const_iterator it = numericsNames.begin();
    for(; it != numericsNames.end(); ++it) {
      numericsMask[*it] = attributeIndex;
      ++attributeIndex;
    }
    return true;
  }
}

vector<unsigned int> Dataset::MaskGetAttributeIndices(AttributeType attrType) {
  vector<unsigned int> indices;
  if(attrType == DISCRETE_TYPE) {
    map<string, unsigned int>::const_iterator it = attributesMask.begin();
    for(; it != attributesMask.end(); ++it) {
      indices.push_back(it->second);
    }
  } else {
    map<string, unsigned int>::const_iterator it = numericsMask.begin();
    for(; it != numericsMask.end(); ++it) {
      indices.push_back(it->second);
    }
  }
  return indices;
}

const map<string, unsigned int>&
Dataset::MaskGetAttributeMask(AttributeType attrType) {
  if(attrType == DISCRETE_TYPE) {
    return attributesMask;
  } else {
    return numericsMask;
  }
}

vector<string> Dataset::MaskGetAllVariableNames() {
  vector<string> names;
  map<string, unsigned int>::const_iterator ait = attributesMask.begin();
  for(; ait != attributesMask.end(); ++ait) {
    names.push_back(ait->first);
  }
  map<string, unsigned int>::const_iterator nit = numericsMask.begin();
  for(; nit != numericsMask.end(); ++nit) {
    names.push_back(nit->first);
  }
  return names;
}

bool Dataset::MaskRemoveInstance(std::string instanceId) {
  map<string, unsigned int>::iterator
  pos = instancesMask.find(instanceId);
  if(pos != instancesMask.end()) {
    instancesMask.erase(pos);
  } else {
    cerr << "ERROR: Dataset::MaskRemoveInstance failed for instance ID: "
            << instanceId << endl;
    return false;
  }
  return true;
}

bool Dataset::MaskSearchInstance(string instanceId) {
  map<string, unsigned int>::iterator
  pos = instancesMask.find(instanceId);
  if(pos != instancesMask.end()) {
    return true;
  } else {
    return false;
  }
}

bool Dataset::MaskIncludeAllInstances() {
  instancesMask.clear();
  vector<string>::const_iterator it = instanceIds.begin();
  unsigned int instanceIndex = 0;
  for(; it != instanceIds.end(); ++it) {
    instancesMask[*it] = instanceIndex;
    ++instanceIndex;
  }

  return true;
}

vector<unsigned int> Dataset::MaskGetInstanceIndices() {
  vector<unsigned int> indices;
  map<string, unsigned int>::const_iterator it = instancesMask.begin();
  for(; it != instancesMask.end(); ++it) {
    indices.push_back(it->second);
  }
  return indices;
}

vector<string> Dataset::MaskGetInstanceIds() {
  vector<string> ids;
  map<string, unsigned int>::const_iterator it = instancesMask.begin();
  for(; it != instancesMask.end(); ++it) {
    ids.push_back(it->first);
  }
  return ids;
}

const map<string, unsigned int>& Dataset::MaskGetInstanceMask() {
  return instancesMask;
}

bool Dataset::MaskPushAll() {
  if(!maskIsPushed) {
    attributesMaskPushed = attributesMask;
    numericsMaskPushed = numericsMask;
    instancesMaskPushed = instancesMask;
    maskIsPushed = true;
    return true;
  } else {
    cerr << "ERROR: only one level of pushing is allowed for "
            << "attribute masks" << endl;
  }
  return false;
}

bool Dataset::MaskPopAll() {
  if(maskIsPushed) {
    attributesMask = attributesMaskPushed;
    numericsMask = numericsMaskPushed;
    instancesMask = instancesMaskPushed;
    maskIsPushed = false;
    return true;
  } else {
    cerr << "ERROR: attempt to pop an unpushed attribute mask" << endl;
    return false;
  }
}

bool Dataset::MaskWriteNewDataset(string newDatasetFilename) {

  ofstream outFile;
  outFile.open(newDatasetFilename.c_str());
  if(outFile.bad()) {
    cerr << "ERROR: Could not open new data set file [" << newDatasetFilename
            << "] for writing" << endl;
    return false;
  }

  if(HasNumerics()) {
    vector<string> numericNames = GetNumericsNames();
    outFile << join(numericNames.begin(), numericNames.end(), "\t");
    outFile << "Class" << endl;
    for(unsigned int i = 0; i < NumInstances(); ++i) {
      ostringstream newInstanceString;
      for(unsigned int j = 0; j < numericNames.size() - 1; ++j) {
        newInstanceString << GetNumeric(i, numericNames[j]) << "\t";
      }
      if(hasContinuousPhenotypes) {
        outFile << newInstanceString.str() << instances[i]->GetPredictedValueTau() << endl;
      } else {
        outFile << newInstanceString.str() << instances[i]->GetClass() << endl;
      }
    }
  } else {
    vector<string> attributeNames = GetAttributeNames();
    outFile << join(attributeNames.begin(), attributeNames.end(), "\t");
    outFile << "Class" << endl;
    for(unsigned int i = 0; i < NumInstances(); ++i) {
      ostringstream newInstanceString;
      for(unsigned int j = 0; j < attributeNames.size() - 1; ++j) {
        newInstanceString << GetAttribute(i, attributeNames[j]) << "\t";
      }
      if(hasContinuousPhenotypes) {
        outFile << newInstanceString.str() << instances[i]->GetPredictedValueTau() << endl;
      } else {
        outFile << newInstanceString.str() << instances[i]->GetClass() << endl;
      }
    }
  }
  outFile.close();

  return true;
}

void Dataset::PrintMaskStats() {
  cout << endl << Timestamp() << "Dataset Mask Statistics" << endl
          << Timestamp() << "Attributes mask size: " << attributesMask.size() << endl
          << Timestamp() << "Numerics mask size: " << numericsMask.size() << endl
          << Timestamp() << "Instances mask size: " << instancesMask.size() << endl
          << Timestamp() << "Mask is pushed? "
          << (maskIsPushed ? "true" : "false") << endl;
}

void Dataset::RunSnpDiagnosticTests(string logFilename,
                                    double globalGenotypeThreshold,
                                    unsigned int cellThreshold) {
  cout << endl << "BEGIN disgnostic tests BEGIN" << endl << endl;
  cout << "Diagnostic Information for: [" << snpsFilename << "]" << endl;
  PrintStatsSimple();
  cout << endl;

  if(!hasGenotypes) {
    cout << "WARNING: This dataset does not have any SNPs. "
            << "Cannnot perform diagnostic tests" << endl;
    return;
  }

  map<string, vector<string> > screwySnps;
  vector<string> badSnps;

  // ---------------------------------------------------------------------------
  // check for missing data in instances (subjects) - bcw - 7/12/11
  cout << "-missing values check... ";
  unsigned int totalMissing = 0;
  map<string, vector<unsigned int> >::const_iterator mit =
          missingValues.begin();
  for(; mit != missingValues.end(); ++mit) {
    totalMissing += (*mit).second.size();
  }
  cout << missingValues.size() << " individuals had a total of "
          << totalMissing << " missing values" << endl;

  // ---------------------------------------------------------------------------
  cout << "-global frequency count threshold check... ";
  unsigned int attributeIndex = 0;
  unsigned int instanceThreshold = (unsigned int) (NumInstances() * globalGenotypeThreshold);
  vector<map<AttributeLevel, unsigned int> >::const_iterator lcIt;
  // PrintLevelCounts();
  unsigned int freqCountBad = 0;
  for(lcIt = levelCounts.begin(); lcIt != levelCounts.end(); ++lcIt, ++attributeIndex) {
    vector<unsigned int> genotypeCounts;
    map<AttributeLevel, unsigned int>::const_iterator countsIt = (*lcIt).begin();
    for(unsigned int levelIdx = 0; countsIt != (*lcIt).end(); ++countsIt, ++levelIdx) {
      string attributeName = attributeNames[attributeIndex];
      if((*countsIt).second < instanceThreshold) {
        //        cout << "Attribute [" << attributeName
        //                << "] does not pass global threshold check of [" << instanceThreshold << "]. "
        //                << "Level: [" << (*countsIt).first
        //                << "], count : [" << (*countsIt).second
        //                << "]" << endl;
        ostringstream gtmsg;
        gtmsg << "does not pass global threshold check of ["
                << instanceThreshold << "]. "
                << "Level: [" << (*countsIt).first
                << "], count : [" << (*countsIt).second
                << "]";

        screwySnps[attributeName].push_back(gtmsg.str());
        ++freqCountBad;
      }
      if((*countsIt).second == 0) {
        //        cout << "Attribute [" << attributeName
        //                << "] Level: [" << levelIdx
        //                << "] is not represented in the dataset" << endl << endl;
        ostringstream mgmsg;
        mgmsg << "Level: [" << levelIdx
                << "] is not represented in the dataset";
        screwySnps[attributeName].push_back(mgmsg.str());
        ++freqCountBad;
      }
      genotypeCounts.push_back((*countsIt).second);
    }

    //    if(!CheckHardyWeinbergEquilibrium(genotypeCounts)) {
    //      cerr << "Attribute [" << attributeNames[attributeIndex]
    //              << "] fails the HWE test" << endl << endl;
    //      screwySnps[attributeNames[attributeIndex]].push_back("HWE");
    //    }
  }
  cout << freqCountBad << " bad frequency counts" << endl;

  // ---------------------------------------------------------------------------
  cout << "-chi-squared cell count threshold check... ";
  ChiSquared chiSq(this);
  unsigned int chiSqBad = 0;
  for(unsigned int attributeIndex = 0;
      attributeIndex < NumAttributes();
      attributeIndex++) {
    chiSq.ComputeScore(attributeIndex);
    vector<vector<double> > frequencyCounts = chiSq.GetFrequencyCounts();
    for(unsigned int curClass = 0; curClass < NumClasses(); curClass++) {
      for(unsigned int curLevel = 0; curLevel < NumLevels(attributeIndex); curLevel++) {
        if(frequencyCounts[curClass][curLevel] < cellThreshold) {
          //          cout << "Attribute [" << attributeNames[attributeIndex]
          //                  << "] does not pass x^2 cell threshold check of [" << cellThreshold << "]. "
          //                  << "Level: [" << curLevel << "], Class: [" << curClass << "] => "
          //                  << "count: [" << frequencyCounts[curClass][curLevel]
          //                  << "]" << endl;
          ostringstream chimsg;
          chimsg << "does not pass x^2 cell threshold check of [" << cellThreshold << "]. "
                  << "Level: [" << curLevel << "], Class: [" << curClass << "] => "
                  << "count: [" << frequencyCounts[curClass][curLevel]
                  << "]";
          screwySnps[attributeNames[attributeIndex]].push_back(chimsg.str());
          ++chiSqBad;
        }
      }
    }
  }
  cout << chiSqBad << " bad x^2 threshold checks" << endl << endl;

  ofstream outFile;
  outFile.open(logFilename.c_str());
  if(outFile.bad()) {
    cerr << "Could not open levels counts file for writing" << endl;
    exit(1);
  }
  map<string, vector<string> >::const_iterator it = screwySnps.begin();
  for(; it != screwySnps.end(); ++it) {
    outFile << (*it).first << "\n"
            << "\t" << join((*it).second.begin(), (*it).second.end(), "\n\t") << endl;
  }
  outFile.close();

  cout << "END diagnostic tests END" << endl << endl;
}

bool
Dataset::CheckHardyWeinbergEquilibrium(vector<unsigned int> genotypeCounts) {
  /// observered counts
  unsigned int AA = genotypeCounts[0];
  unsigned int Aa = genotypeCounts[1];
  unsigned int aa = genotypeCounts[2];
  unsigned int sum = AA + Aa + aa;
  if(sum < 1) {
    cerr << "ERROR: Hardy-Weinberg test failed to find any observed values"
            << endl;
    return false;
  }
  if(sum != NumInstances()) {
    cerr << "ERROR: genotype counts sum [" << sum
            << "] does not add up to number of instances ["
            << NumInstances() << "]" << endl;
    PrintLevelCounts();
    return false;
  }

  /// HWE probabilities
  double fAA = (double) AA / sum;
  double fAa = (double) Aa / sum;
  double p = fAA + 0.5 * fAa;
  double q = 1.0 - p;
  //  cout << p << " " << q << " " << endl;

  /// expected values
  double eAA = p * p * sum;
  double eAa = 2.0 * p * q * sum;
  double eaa = q * q * sum;
  if(eAA == 0 || eAa == 0 || eaa == 0) {
    cout << "ERROR: HWE chi-square could not be computed due to zeroes "
            << "in expected counts" << endl;
    return false;
  }
  unsigned int expectedSum = (unsigned int) (eAA + eAa + eaa + 0.5);
  if(expectedSum != NumInstances()) {
    cerr << "ERROR: HWE expected genotype counts sum [" << expectedSum
            << "] does not add up to number of instances ["
            << NumInstances() << "]" << endl;
    return false;
  }

  /// perform Pearson's chi-squared test
  double chiSquaredSum =
          ((AA - eAA) / eAA) +
          ((Aa - eAa) / eAa) +
          ((aa - eaa) / eaa);

  //  cout << eAA << " " << eAa << " " << eaa << endl;

  /// one degree of freedom (# genotypes - # alleles), 5% significance level
  double pValue = 1.0 - gsl_cdf_chisq_Q(chiSquaredSum, 1);
  if(pValue > 0.1) {
    cerr << "ERROR: x^2 p-value too high: [" << pValue << "]" << endl;
    cerr << AA << " / " << Aa << " / " << aa << endl;
    return false;
  }

  return true;
}

double Dataset::SNPHWE(int obs_hets, int obs_hom1, int obs_hom2) {
  if(obs_hom1 < 0 || obs_hom2 < 0 || obs_hets < 0) {
    printf("FATAL ERROR - SNP-HWE: Current genotype configuration "
           "(%d  %d %d ) includes a negative count",
           obs_hets, obs_hom1, obs_hom2);
    exit(EXIT_FAILURE);
  }

  int obs_homc = obs_hom1 < obs_hom2 ? obs_hom2 : obs_hom1;
  int obs_homr = obs_hom1 < obs_hom2 ? obs_hom1 : obs_hom2;

  int rare_copies = 2 * obs_homr + obs_hets;
  int genotypes = obs_hets + obs_homc + obs_homr;

  double * het_probs = (double *) malloc((size_t) (rare_copies + 1) * sizeof(double));
  if(het_probs == NULL) {
    printf("FATAL ERROR - SNP-HWE: Unable to allocate array for "
           "heterozygote probabilities");
    exit(EXIT_FAILURE);
  }

  int i;
  for(i = 0; i <= rare_copies; i++)
    het_probs[i] = 0.0;

  /* start at midpoint */
  int mid = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes);

  /* check to ensure that midpoint and rare alleles have same parity */
  if((rare_copies & 1) ^ (mid & 1))
    mid++;

  int curr_hets = mid;
  int curr_homr = (rare_copies - mid) / 2;
  int curr_homc = genotypes - curr_hets - curr_homr;

  het_probs[mid] = 1.0;
  double sum = het_probs[mid];
  for(curr_hets = mid; curr_hets > 1; curr_hets -= 2) {
    het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets *
            (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0));
    sum += het_probs[curr_hets - 2];

    /* 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote */
    curr_homr++;
    curr_homc++;
  }

  curr_hets = mid;
  curr_homr = (rare_copies - mid) / 2;
  curr_homc = genotypes - curr_hets - curr_homr;
  for(curr_hets = mid; curr_hets <= rare_copies - 2; curr_hets += 2) {
    het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc
            / ((curr_hets + 2.0) * (curr_hets + 1.0));
    sum += het_probs[curr_hets + 2];

    /* add 2 heterozygotes for next iteration -> subtract one rare,
     * one common homozygote */
    curr_homr--;
    curr_homc--;
  }

  for(i = 0; i <= rare_copies; i++)
    het_probs[i] /= sum;

  /* alternate p-value calculation for p_hi/p_lo
  double p_hi = het_probs[obs_hets];
  for (i = obs_hets + 1; i <= rare_copies; i++)
    p_hi += het_probs[i];

  double p_lo = het_probs[obs_hets];
  for (i = obs_hets - 1; i >= 0; i--)
     p_lo += het_probs[i];


  double p_hi_lo = p_hi < p_lo ? 2.0 * p_hi : 2.0 * p_lo;
   */

  double p_hwe = 0.0;
  /*  p-value calculation for p_hwe  */
  for(i = 0; i <= rare_copies; i++) {
    if(het_probs[i] > het_probs[obs_hets])
      continue;
    p_hwe += het_probs[i];
  }

  p_hwe = p_hwe > 1.0 ? 1.0 : p_hwe;

  free(het_probs);

  return p_hwe;
}

double Dataset::GetClassProbability(ClassLevel thisClass) {
  if(NumInstances()) {
    return((double) classIndexes[thisClass].size() / (double) NumInstances());
  }
  return 0.0;
}

double Dataset::GetProbabilityValueGivenClass(unsigned int attributeIndex,
                                              AttributeLevel A,
                                              ClassLevel classValue) {
  map<pair<AttributeLevel, ClassLevel>, unsigned int> thisAttrLevelCounts =
          levelCountsByClass[attributeIndex];
  unsigned int instancesInThisClass = classIndexes[classValue].size();
  pair<AttributeLevel, ClassLevel> key = make_pair(A, classValue);
  if(thisAttrLevelCounts.find(key) != thisAttrLevelCounts.end()) {
    //    cout << attributeIndex << " [" << A << "|" << classValue << "]: "
    //            << thisAttrLevelCounts[key] << "/" << instancesInThisClass << endl;
    return((double) thisAttrLevelCounts[key] / (double) instancesInThisClass);
  } else {
    return 0.0;
  }
}

void Dataset::AttributeInteractionInformation() {
  map<pair<unsigned int, unsigned int>, map<string, double> > results;
  CalculateInteractionInformation(results);

  cout << endl
          << "---------------------------------" << endl
          << "Attribute Interaction Information" << endl
          << "---------------------------------" << endl
          << endl;

  // display results header line
  cout << "A\tB\t";
  map<pair<unsigned int, unsigned int>, map<string, double> >::const_iterator
  rIt = results.begin();
  map<string, double> hdrResults = (*rIt).second;
  map<string, double>::const_iterator rdIt;
  for(rdIt = hdrResults.begin(); rdIt != hdrResults.end(); rdIt++) {
    cout << setw(10) << (*rdIt).first;
  }
  cout << endl;

  /// get the column sum
  double metricSum = 0.0;
  for(rIt = results.begin(); rIt != results.end(); rIt++) {
    map<string, double> thisResults = (*rIt).second;
    for(rdIt = thisResults.begin(); rdIt != thisResults.end(); rdIt++) {
      string metric = (*rdIt).first;
      if(metric == "I(A;B|C)") {
        metricSum += (*rdIt).second;
      }
    }
  }

  /// display results detail; selected column as percentages
  for(rIt = results.begin(); rIt != results.end(); rIt++) {
    pair<unsigned int, unsigned int> thisCombo = (*rIt).first;
    map<string, double> thisResults = (*rIt).second;
    cout << thisCombo.first << "\t" << thisCombo.second << "\t";
    for(rdIt = thisResults.begin(); rdIt != thisResults.end(); rdIt++) {
      string metric = (*rdIt).first;
      if(metric == "I(A;B|C)") {
        cout << fixed << setw(10) << setprecision(3) << (*rdIt).second;
        cout << fixed << setw(10) << setprecision(3)
                << "( " << (((*rdIt).second / metricSum) * 100) << "% )";
      } else {
        cout << fixed << setw(10) << setprecision(3) << (*rdIt).second;
      }
    }
    cout << endl;
  }
  cout << endl;
}

void Dataset::CalculateInteractionInformation(
                                              map<pair<unsigned int, unsigned int>, map<string, double> >& results) {

  vector<AttributeLevel> a;
  vector<AttributeLevel> b;
  vector<ClassLevel> c;
  vector<AttributeLevel> ab;

  if(!HasGenotypes() || HasNumerics()) {
    cerr << "ERROR: Dataset::CalculateInteractionInformation "
            << "only works on SNP data" << endl;
    return;
  }

  GetClassValues(c);
  //copy(c.begin(), c.end(), ostream_iterator<InstanceClass>(cout, "\n"));

  // for all possible (unique) interactions, ie nCk
  // #pragma omp parallel for
  for(int i = 0; i < (int) NumAttributes(); i++) {
    for(int j = i + 1; j < (int) NumAttributes(); j++) {

      // cout << "Calculating attributes: " << i << "," << j << endl;
      // load attribute values (columns) into vectors for math routines
      GetAttributeValues(attributeNames[i], a);
      GetAttributeValues(attributeNames[j], b);
      //      cout << "A:" << endl;
      //      copy(a.begin(), a.end(), ostream_iterator<Attribute>(cout, "\n"));
      //      cout << "B:" << endl;
      //      copy(b.begin(), b.end(), ostream_iterator<Attribute>(cout, "\n"));
      //      cout << endl;

      // construct a new attribute with a and b
      // based on passed constructWith parameter
      ConstructAttributeCart(a, b, ab);

      // compute and save result for this i,j combination
      map<string, double> r;

      r["H(A)"] = Entropy(a);
      r["H(B)"] = Entropy(b);
      r["H(C)"] = Entropy(c);
      r["H(AB)"] = Entropy(ab);

      r["H(A|C)"] = ConditionalEntropy(a, c);
      r["H(B|C)"] = ConditionalEntropy(b, c);
      r["H(AB|C)"] = ConditionalEntropy(ab, c);
      r["H(C|A)"] = ConditionalEntropy(c, a);
      r["H(C|B)"] = ConditionalEntropy(c, b);

      r["I(A;C)"] = r["H(C)"] - r["H(C|A)"];
      r["I(B;C)"] = r["H(C)"] - r["H(C|B)"];
      r["I(A;B)"] = r["H(A)"] + r["H(B)"] - r["H(AB)"];
      r["I(A;B|C)"] = r["H(A|C)"] + r["H(B|C)"] - r["H(AB|C)"];
      r["I(A;B;C)"] = r["I(A;B|C)"] - r["I(A;B)"];
      r["I(AB;C)"] = r["I(A;B;C)"] + r["I(A;C)"] + r["I(B;C)"];

      //			cout << "Adding pair: " << i << ", " << j << endl;
      //      cout << attributeNames[i] << "\t" << attributeNames[j]
      //              << "\t" << r["I(A;B|C)"] << "\t" << r["I(A;C)"] << endl;

      results[make_pair(i, j)] = r;
    }
  }

}

bool Dataset::CalculateGainMatrix(double** gainMatrix) {
  if(!HasGenotypes() || HasNumerics()) {
    cerr << "Dataset::CalculateInteractionInformation only works on SNP data"
            << endl;
    return false;
  }

  vector<ClassLevel> c;
  GetClassValues(c);

  /// for all possible (unique) interactions, ie nCk
  cout << Timestamp() << "Constructing attribute interaction matrix in parallel:"
          << endl << Timestamp();
  int numAttributes = NumAttributes();
#pragma omp parallel for schedule(dynamic, 1)
  for(int i = 0; i < numAttributes; i++) {
    for(int j = i; j < (int) numAttributes; j++) {

      vector<AttributeLevel> a;
      vector<AttributeLevel> b;
      vector<AttributeLevel> ab;

      // cout << "Calculating attributes: " << i << "," << j << endl;
      // TODO: optimization here since new a is old b after first iteration
      GetAttributeValues(attributeNames[i], a);
      GetAttributeValues(attributeNames[j], b);

      ConstructAttributeCart(a, b, ab);

      /// compute gainMatrix[i][j]
      map<string, double> r;

      r["H(A)"] = Entropy(a);
      r["H(B)"] = Entropy(b);
      r["H(C)"] = Entropy(c);
      r["H(AB)"] = Entropy(ab);

      r["H(A|C)"] = ConditionalEntropy(a, c);
      r["H(B|C)"] = ConditionalEntropy(b, c);
      r["H(AB|C)"] = ConditionalEntropy(ab, c);
      r["H(C|A)"] = ConditionalEntropy(c, a);
      r["H(C|B)"] = ConditionalEntropy(c, b);

      r["I(A;C)"] = r["H(C)"] - r["H(C|A)"];
      r["I(B;C)"] = r["H(C)"] - r["H(C|B)"];
      r["I(A;B)"] = r["H(A)"] + r["H(B)"] - r["H(AB)"];
      r["I(A;B|C)"] = r["H(A|C)"] + r["H(B|C)"] - r["H(AB|C)"];
      r["I(A;B;C)"] = r["I(A;B|C)"] - r["I(A;B)"];
      r["I(AB;C)"] = r["I(A;B;C)"] + r["I(A;C)"] + r["I(B;C)"];

      //      map<string, double>::const_iterator it = r.begin();
      //      for(; it != r.end(); ++it) {
      //        cout << (*it).first << " = " << (*it).second << endl;
      //      }
      //      exit(1);

      //      def interaction_information(self, attrA, attrB):
      //        # I(A;B;C)=I(A;B|C)-I(A;B)
      //        # I(A;B;C)=H(AB)+H(BC)+H(AC)-H(A)-H(B)-H(C)-H(ABC)
      //        # where C is the class
      //
      //        H_ABC   = self.entropy(attrA, attrB, self.class_idx)
      //        H_AB    = self.entropy(attrA, attrB)
      //        H_AC    = self.entropy(attrA, self.class_idx)
      //        H_BC    = self.entropy(attrB, self.class_idx)
      //        H_A     = self.entropy(attrA)
      //        H_B     = self.entropy(attrB)
      //        H_C     = self.entropy(self.class_idx)
      //        return H_AB + H_BC + H_AC - H_A - H_B - H_C - H_ABC
      //
      //      def autointeraction(self, attr):
      //              return self.entropy(attr, self.class_idx) - self.entropy(attr) - self.entropy(self.class_idx)

      // I_3 below from:
      // McKinney BA, Crowe JE Jr, Guo J, Tian D (2009) Capturing the Spectrum
      // of Interaction Effects in Genetic Association Studies by Simulated
      // Evaporative Cooling Network Analysis. PLoS Genet 5(3): e1000432.
      // doi:10.1371/journal.pgen.1000432
//      double I_3_pygain = r["H(AB)"] + r["H(B|C)"] + r["H(A|C)"] -
//              r["H(A)"] - r["H(B)"] - r["H(C)"] - r["H(AB|C)"];
      double I_3_paper = r["I(AB;C)"] - r["I(A;C)"] - r["I(B;C)"];
      // are the two I_3 values not close?
      //      if(fabs(I_3_pygain - I_3_paper) > 1e-6) {
      //        cerr << "WARNING: I_3 values are not the same. "
      //                << "I_3 from pygain: " << I_3_pygain
      //                << " I_3 from paper:  " << I_3_paper
      //                << " Using formula from pygain" << endl;
      //      }
      double I_3 = I_3_paper;
      double I_2 = r["H(A|C)"] - r["H(A)"] - r["H(C)"];

      if(i == j) {
        gainMatrix[i][j] = I_2;
      } else {
        gainMatrix[i][j] = gainMatrix[j][i] = I_3;
      }

      // DEBUG
      //      cout << attributeNames[i] << "\t" << attributeNames[j]
      //              << "\t" << r["I(A;B|C)"] << "\t" << r["I(A;C)"] << endl;
    }
    // happy lights
    if(i && (i % 100 == 0)) {
       cout << i << "/" << numAttributes << " ";
       cout.flush();
     }
     if(i && ((i % 1000) == 0)) {
       cout << endl << Timestamp();
     }
   }
   cout << numAttributes << "/" << numAttributes << " done" << endl;

  return true;
}

bool Dataset::LoadSnps(std::string filename) {
  snpsFilename = filename;
  ifstream dataStream(snpsFilename.c_str());
  if(!dataStream.is_open()) {
    cerr << "ERROR: Could not open SNP data set: " << snpsFilename << endl;
    exit(-1);
  }
  cout << Timestamp() << "Reading whitespace-delimited SNP data set lines from "
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
  for(it = tokens.begin(); it != tokens.end(); ++it) {
    if(to_upper(*it) == "CLASS") {
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
  genotypeCounts.resize(numAttributes);

  vector<string> missingValuesToCheck;
  missingValuesToCheck.push_back("?");
  missingValuesToCheck.push_back("-9");

  // read instance attributes from whitespace-delimited lines
  unsigned int instanceIndex = 0;
  unsigned int lineNumber = 0;
  cout << Timestamp();
  ValueType classType = NO_VALUE;
  double minPheno = 0.0, maxPheno = 0.0;
  while(getline(dataStream, line)) {
    ++lineNumber;
    // only load matching IDs from numerics and phenotype files
    // the delimited text file uses line numbers for instance IDs
    string trimmedLine = trim(line);
    // skip blank lines in the data section (usually end of file)
    if(!trimmedLine.size()) {
      continue;
    }

    string ID = zeroPadNumber(lineNumber, 8);
    if(!IsLoadableInstanceID(ID)) {
      cout << Timestamp() << "WARNING: Dataset ID [" << ID << "] skipped. "
              << " Line number: " << lineNumber
              << ". Not found in numerics and/or phenotype file(s)"
              << endl;
      continue;
    }

    vector<string> attributesStringVector;
    split(attributesStringVector, trimmedLine);
    unsigned int numAttributesRead = attributesStringVector.size() - 1;
    if((numAttributesRead == 0) || (numAttributesRead != numAttributes)) {
      cout << Timestamp() << "WARNING: Skipping line " << lineNumber
              << " instance has " << numAttributesRead
              << " should have " << numAttributes
              << endl;
      continue;
    }

    vector<AttributeLevel> attributeVector;
    unsigned int attrIdx = 0;
    unsigned int thisCol = 0;
    vector<string>::const_iterator it = attributesStringVector.begin();
    ClassLevel discreteClassLevel = MISSING_DISCRETE_CLASS_VALUE;
    NumericLevel numericClassLevel = MISSING_NUMERIC_CLASS_VALUE;
    for(; it != attributesStringVector.end(); ++it, ++thisCol) {
      string thisAttr = *it;
      if(thisCol == classColumn) {
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
              cerr << "ERROR: Could not get class discrete level on line: "
                      << lineNumber << endl;
              return false;
            } else {
              if(discreteClassLevel == MISSING_DISCRETE_CLASS_VALUE) {
                cout << Timestamp()
                        << "WARNING: missing discrete phenotype skipped on line: "
                        << lineNumber << endl;
                continue;
              }
            }
            break;
          case NUMERIC_VALUE:
            numericClassLevel = MISSING_NUMERIC_CLASS_VALUE;
            if(!GetNumericClassLevel(thisAttr, missingValuesToCheck,
                                     numericClassLevel)) {
              cerr << "ERROR: Could not get numeric class level on line: "
                      << lineNumber << endl;
              return false;
            } else {
              if(numericClassLevel == MISSING_NUMERIC_CLASS_VALUE) {
                cout << Timestamp()
                        << "WARNING: missing numeric phenotype skipped on line: "
                        << lineNumber << endl;
                continue;
              }
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
          missingValues[ID].push_back(attrIdx);
        } else {
          attributeLevelsSeen[attrIdx].insert(thisAttr);
          ++genotypeCounts[attrIdx][thisAttr];
        }
        attributeVector.push_back(thisAttrLevel);
        ++attrIdx;
      }
    }

    // create an instance from the vector of attribute and class values
    if(attributeVector.size() != numAttributes) {
      cerr << "ERROR: Number of attributes parsed on line " << lineNumber
              << ": " << attributesStringVector.size()
              << " is not equal to the number of attributes"
              << " read from the data file header: " << numAttributes
              << endl;
      return false;
    }
    if(classType == NO_VALUE) {
      cerr << "Could not determine class type on line 1" << endl;
      return false;
    }
    DatasetInstance* newInst = new DatasetInstance(this);
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
      cerr << "ERROR: loading tab-delimited data set. "
              << "Could not create dataset instance for line number "
              << lineNumber << endl;
      return false;
    }
    ++instanceIndex;

    // happy lights
    if((lineNumber - 1) && ((instanceIndex % 100) == 0)) {
      cout << instanceIndex << " ";
      cout.flush();
    }
    if((lineNumber - 1) && ((instanceIndex % 1000) == 0)) {
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

bool Dataset::GetAttributeLevel(string inLevel,
                                vector<string> missingValues,
                                AttributeLevel& outLevel) {
  if(find(missingValues.begin(), missingValues.end(), inLevel) !=
     missingValues.end()) {
    outLevel = MISSING_DISCRETE_CLASS_VALUE;
  } else {
    if((inLevel == "0") || (inLevel == "1") || (inLevel == "2")) {
      outLevel = lexical_cast<AttributeLevel > (inLevel);
    } else {
      outLevel = INVALID_ATTRIBUTE_VALUE;
      return false;
    }
  }

  return true;
}

bool Dataset::GetDiscreteClassLevel(string inLevel,
                                    vector<string> missingValues,
                                    ClassLevel& outLevel) {
  if(find(missingValues.begin(), missingValues.end(), inLevel) !=
     missingValues.end()) {
    outLevel = MISSING_DISCRETE_CLASS_VALUE;
  } else {
    if((inLevel == "0") || (inLevel == "1")) {
      outLevel = lexical_cast<ClassLevel > (inLevel);
    } else {
      outLevel = INVALID_DISCRETE_CLASS_VALUE;
      return false;
    }
  }

  return true;
}

bool Dataset::GetNumericClassLevel(string inLevel,
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

ValueType Dataset::GetAttributeValueType(string value,
                                         vector<string> missingValues) {
  ValueType returnValueType = NO_VALUE;
  if(find(missingValues.begin(), missingValues.end(), value) !=
     missingValues.end()) {
    return MISSING_VALUE;
  } else {
    if((value == "0") || (value == "1") || (value == "2")) {
      return DISCRETE_VALUE;
    } else {
      return NUMERIC_VALUE;
    }
  }
  return returnValueType;
}

ValueType Dataset::GetClassValueType(string value,
                                     vector<string> missingValues) {
  ValueType returnValueType = NO_VALUE;
  if(find(missingValues.begin(), missingValues.end(), value) !=
     missingValues.end()) {
    return MISSING_VALUE;
  } else {
    if((value == "0") || (value == "1")) {
      return DISCRETE_VALUE;
    } else {
      return NUMERIC_VALUE;
    }
  }
  return returnValueType;
}

void Dataset::UpdateAllLevelCounts() {
  cout << Timestamp() << "Updating all level counts:" << endl << Timestamp();
  levelCounts.clear();
  levelCounts.resize(NumAttributes());
  if(!hasContinuousPhenotypes) {
    levelCountsByClass.clear();
    levelCountsByClass.resize(NumAttributes());
  }
  unsigned int instanceCount = 0;
  map<string, unsigned int>::const_iterator it = instancesMask.begin();
  for(; it != instancesMask.end(); ++it) {
    UpdateLevelCounts(instances[it->second]);
    if(instanceCount && ((instanceCount % 100) == 0)) {
      cout << instanceCount << " ";
      cout.flush();
    }
    if(instanceCount && ((instanceCount % 1000) == 0)) {
      cout << endl << Timestamp();
    }
    ++instanceCount;
  }
  cout << instanceCount << endl;
}

void Dataset::UpdateLevelCounts(DatasetInstance* dsi) {
  ClassLevel thisClassLevel = dsi->GetClass();
  map<string, unsigned int>::const_iterator it = attributesMask.begin();
  for(; it != attributesMask.end(); ++it) {
    AttributeLevel thisAttributeLevel = dsi->GetAttribute(it->second);
    if(thisAttributeLevel != MISSING_ATTRIBUTE_VALUE) {
      ++levelCounts[it->second][thisAttributeLevel];
      if(!hasContinuousPhenotypes) {
        ++levelCountsByClass[it->second][make_pair(thisAttributeLevel,
                                                   thisClassLevel)];
      }
    }
  }
}

bool Dataset::LoadNumerics(string filename) {
  numericsFilename = filename;
  ifstream dataStream(numericsFilename.c_str());
  if(!dataStream.is_open()) {
    cerr << "ERROR: Could not open numerics file: "
            << numericsFilename << endl;
    return false;
  }
  cout << Timestamp() << "Reading numerics from " << numericsFilename << endl;
  //  PrintStats();

  // temporary string for reading file lines
  unsigned int lineNumber = 0;
  string line;
  // read the header for covariate names
  getline(dataStream, line);
  ++lineNumber;
  vector<string> numNames;
  split(numNames, line);
  if(numNames.size() < 3) {
    cerr << "ERROR: Covariate file must have at least three columns: "
            << "FID IID COV1 ... COVN" << endl;
    return false;
  }
  numericsNames.resize(numNames.size() - 2);
  copy(numNames.begin() + 2, numNames.end(), numericsNames.begin());
  vector<string>::const_iterator it = numericsNames.begin();
  unsigned int numIdx = 0;
  for(; it != numericsNames.end(); ++it) {
    numericsMask[*it] = numIdx;
    ++numIdx;
  }

  // if no snp data then need to create instances in this loop - 6/19/11
  // read each new set of numerics
  map<string, bool> idsSeen;
  unsigned int newInstanceIdx = 0;
  DatasetInstance* tempInstance = 0;
  unsigned int instanceIndex = 0;
  while(getline(dataStream, line)) {
    ++lineNumber;

    // cout << lineNumber << ": " << line << endl;
    vector<string> numericsStringVector;
    split(numericsStringVector, line);
    // cout << line << endl;
    if(numericsStringVector.size() < 3) {
      cerr << "ERROR: Covariate file must have at least three columns: "
              << "FID IID COV1 ... COVN" << endl;
      return false;
    }
    if(numericsNames.size() != (numericsStringVector.size() - 2)) {
      cerr << "ERROR: Number of numeric values read from the covariate file header: ["
              << numericsNames.size() << "] is not equal to the number of numeric"
              << " values read [" << (numericsStringVector.size() - 2)
              << "] on line: " << lineNumber << " of " << numericsFilename << endl;
      return false;
    }

    // first column is the ID for matching data set rows - bcw - 8/3/11
    string ID = numericsStringVector[0];
    if(!IsLoadableInstanceID(ID)) {
      cout << Timestamp() << "WARNING: Skipping Numeric ID [" << ID << "]. "
              << "It does not match the data set and/or phenotype file"
              << endl;
      continue;
    }

    if(idsSeen.find(ID) == idsSeen.end()) {
      idsSeen[ID] = true;
    } else {
      cout << Timestamp() << "WARNING: Duplicate ID [" << ID << "] detected and "
              << "skipped on line [" << lineNumber << "]" << endl;
      continue;
    }

    // cout << "Numerics ID string from file: " << thisID << endl;
    numericsIds.push_back(ID);
    if(!hasGenotypes) {
      instancesMask[ID] = newInstanceIdx++;
      tempInstance = new DatasetInstance(this);
    }
    // skip the first two columns: familiy and individual IDs
    vector<string>::const_iterator it = numericsStringVector.begin() + 2;
    // add new numeric columns to this instance
    unsigned int numericsIndex = 0;
    for(; it != numericsStringVector.end(); it++) {
      //      cout << "instance " << instanceIndex << ", read from file: " << *it
      //              << ", numeric value: " << strtod((*it).c_str(), NULL) << endl;
      NumericLevel thisValue = 0.0;
      if((*it == "-9") || (*it == "?")) {
        thisValue = MISSING_NUMERIC_VALUE;
        missingNumericValues[ID].push_back(numericsIndex);
      } else {
        thisValue = lexical_cast<NumericLevel > (*it);
      }
      if(!hasGenotypes) {
        tempInstance->AddNumeric(thisValue);
      } else {
        unsigned int lookupIdIndex = 0;
        if(GetInstanceIndexForID(ID, lookupIdIndex)) {
          // cout << "ID: " << thisID << ", Lookup index: " << lookupIdIndex << ", numInstances = " << instances.size() << endl;
          instances[lookupIdIndex]->AddNumeric(thisValue);
        }
      }
      ++numericsIndex;
      // cout << "-----------------------------------------------------" << endl;
    }
    if(!hasGenotypes) {
      instances.push_back(tempInstance);
    }
    ++instanceIndex;
  }

  hasNumerics = true;

  // find the min and max values for each numeric attribute
  // used in diff/distance calculation metrics
  vector<NumericLevel> numericColumn;
  for(unsigned int i = 0; i < NumNumerics(); ++i) {
    GetNumericValues(i, numericColumn);
    double minElement = *numericColumn.begin();
    double maxElement = *numericColumn.begin();
    for(vector<NumericLevel>::const_iterator it = numericColumn.begin();
        it != numericColumn.end(); ++it) {
      if((*it != MISSING_NUMERIC_VALUE) && (*it < minElement)) {
        minElement = *it;
      }
      if((*it != MISSING_NUMERIC_VALUE) && (*it > maxElement)) {
        maxElement = *it;
      }
    }
    numericsMinMax.push_back(make_pair<double, double>(minElement, maxElement));
  }

  cout << Timestamp() << "Read " << NumNumerics() << " numeric attributes" << endl;

  return true;
}

bool Dataset::GetNumericValues(unsigned int numericIndex,
                               vector<double>& numericValues) {
  if(hasNumerics) {
    if((numericIndex < 0) || (numericIndex > numericsNames.size())) {
      cerr << "ERROR: Dataset::GetNumericValues: numeric index out of range: "
              << numericIndex << endl;
      return false;
    }
    numericValues.clear();
    map<string, unsigned int>::const_iterator it;
    for(it = instancesMask.begin(); it != instancesMask.end(); it++) {
      double thisNumeric = instances[it->second]->GetNumeric(numericIndex);
      numericValues.push_back(thisNumeric);
    }
  } else {
    cout << Timestamp() << "WARNING: attempting to access numeric data when none "
            << "have been loaded" << endl;
    return false;
  }

  return true;
}

bool Dataset::LoadAlternatePhenotypes(string phenotypesFilename) {
  if((!hasGenotypes) && (!hasNumerics)) {
    cerr << "ERROR: Dataset::LoadAlternatePhenotypes: SNP and/or numeric data "
            << "must be loaded before alternate phenotypes" << endl;
    return false;
  }
  alternatePhenotypesFilename = phenotypesFilename;
  ifstream dataStream(alternatePhenotypesFilename.c_str());
  if(!dataStream.is_open()) {
    cerr << endl << "ERROR: Could not open alternate phenotype file: "
            << alternatePhenotypesFilename << endl;
    return false;
  }
  cout << Timestamp() << "Reading alternate phenotypes from "
          << alternatePhenotypesFilename << "... " << endl;

  // temporary string for reading file lines
  string line;

  // read each new phenotype value
  // remove header? decided we will not use a header
  // getline(dataStream, line);
  unsigned int lineNumber = 0;
  // PrintClassIndexInfo();
  // clear any existing class info for reading alternate phenotype file
  classIndexes.clear();
  unsigned int instancesRead = 0;
  map<string, bool> idsSeen;
  double minPheno = 0.0, maxPheno = 0.0;
  vector<string> idsToDelete;
  while(getline(dataStream, line)) {
    ++lineNumber;
    string trimmedLine = trim(line);
    if(trimmedLine == "") {
      cout << Timestamp() << "WARNING: Line [" << lineNumber
              << "] in phenotype file is blank. Skipping" << endl;
      continue;
    }
    vector<string> lineParts;
    split(lineParts, trimmedLine);
    if(lineParts.size() != 3) {
      cerr << "ERROR: Phenotype file must have at three columns: "
              << "FID IID PHENOTYPE" << endl;
      return false;
    }

    string ID = lineParts[0];

    // skip IDs that don't match common IDs between numerics, phenotypes and data
    if(!IsLoadableInstanceID(ID)) {
      cout << Timestamp() << "WARNING: Skipping alternate phenotype ID [" << ID << "]. "
              << "It does not match the data set and/or numerics file(s)"
              << endl;
      continue;
    }

    // skip duplicate IDs
    if(idsSeen.find(ID) == idsSeen.end()) {
      idsSeen[ID] = true;
    } else {
      cout << Timestamp() << "WARNING: Duplicate ID [" << ID << "] detected and "
              << "skipped on line [" << lineNumber << "]" << endl;
      continue;
    }

    // lookup instance index for this ID; cannot assume in same order
    unsigned int instanceIndex = 0;
    if(!GetInstanceIndexForID(ID, instanceIndex)) {
      cerr << "ERROR: on lookup in GetInstanceIndexForID: "
              << ID << endl;
      return false;
    }

    ClassLevel classValue = 0;
    NumericLevel predictedValue = 0.0;
    string classString = lineParts[2];
    if((lineParts[2] == "-9") || (lineParts[2] == "?")) {
      if(hasContinuousPhenotypes) {
        predictedValue = (NumericLevel) MISSING_DISCRETE_CLASS_VALUE;
      } else {
        classValue = (ClassLevel) MISSING_DISCRETE_CLASS_VALUE;
      }
      cout << Timestamp() << "WARNING: missing phenotype value read from " <<
              "alternate phenotype file line: " << lineNumber << endl;
      idsToDelete.push_back(ID);
      continue;
    } else {
      if(hasContinuousPhenotypes) {
        predictedValue = lexical_cast<NumericLevel > (classString);
        instances[instanceIndex]->SetPredictedValueTau(predictedValue);
        if(lineNumber == 1) {
          minPheno = maxPheno = predictedValue;
        } else {
          if(predictedValue < minPheno) {
            minPheno = predictedValue;
          }
          if(predictedValue > maxPheno) {
            maxPheno = predictedValue;
          }
        }
      } else {
        classValue = lexical_cast<AttributeLevel > (classString);
        instances[instanceIndex]->SetClass(classValue);
        classIndexes[classValue].push_back(instanceIndex);
      }
    }

    ++instancesRead;
  }

  // if we found missing genotypes in the alternate phenotype file load,
  // we need to remove the already loaded instances with the missing phenotype
  // instance IDs - also update meta data - 10/27/11
  if(idsToDelete.size()) {
    for(vector<string>::const_iterator delIt = idsToDelete.begin();
        delIt != idsToDelete.end(); ++delIt) {

      string delId = *delIt;
      unsigned int delIdIndex = instancesMask[delId];
      ClassLevel delClass = instances[delIdIndex]->GetClass();

      // remove instanceIndex from classIndexes
      vector<unsigned int>::iterator ciIt = find(classIndexes[delClass].begin(),
                                                 classIndexes[delClass].end(),
                                                 delIdIndex);
      if(ciIt != classIndexes[delClass].end()) {
        classIndexes[delClass].erase(ciIt);
      }

      // remove ID from instanceIds
      vector<string>::iterator iiIt = find(instanceIds.begin(),
                                           instanceIds.end(),
                                           delId);
      if(iiIt != instanceIds.end()) {
        instanceIds.erase(iiIt);
      }

      // remove from instancesMask
      map<string, unsigned int>::iterator imIt = instancesMask.find(delId);
      if(imIt != instancesMask.end()) {
        instancesMask.erase(imIt);
      }
    }
    classIndexes.erase(MISSING_DISCRETE_CLASS_VALUE);
  }

  // update all counts
  UpdateAllLevelCounts();

  if(hasContinuousPhenotypes) {
    continuousPhenotypeMinMax = make_pair(minPheno, maxPheno);
  }

  hasAlternatePhenotypes = true;

  cout << Timestamp() << "Read " << instancesRead << " phenotypes. "
          << NumInstances() << " instances" << endl;

  return true;
}

bool Dataset::IsLoadableInstanceID(std::string ID) {
  if(!instanceIdsToLoad.size()) {
    return true;
  }
  if(std::find(instanceIdsToLoad.begin(), instanceIdsToLoad.end(), ID)
     == instanceIdsToLoad.end()) {
    return false;
  } else {
    return true;
  }
}
