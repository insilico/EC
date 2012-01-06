/*
 * DatasetInstance.C - Bill White - 6/14/05
 * 
 * Class to hold dataset instances (rows)
 * Reworked entirely for McKinney Lab work - 2/28/11
 */

#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "Dataset.h"
#include "DatasetInstance.h"
#include "StringUtils.h"
#include "best_n.h"
#include "Debugging.h"

using namespace std;
using namespace insilico;

/// functor for T comparison
typedef DistancePair T;

class deref_less_bcw : public std::binary_function<T, T, bool>
{
public:

  bool operator()(const T a, const T b) const {
    return(a.first < b.first);
  }
};

DatasetInstance::DatasetInstance(Dataset* ds) {
  dataset = ds;
  classLabel = MISSING_DISCRETE_CLASS_VALUE;
  predictedValueTau = MISSING_DISCRETE_CLASS_VALUE;
}

DatasetInstance::~DatasetInstance() {
}

Dataset* DatasetInstance::GetDatasetPtr() {
  return dataset;
}

bool
DatasetInstance::LoadInstanceFromVector(vector<AttributeLevel> newAttributes) {
  if(!newAttributes.size()) {
    return false;
  }
  attributes.clear();
  vector<AttributeLevel>::const_iterator it;
  for(it = newAttributes.begin(); it != newAttributes.end(); it++) {
    attributes.push_back(*it);
  }
  return true;
}

unsigned int DatasetInstance::NumAttributes() {
  return(attributes.size());
}

AttributeLevel DatasetInstance::GetAttribute(unsigned int index) {
  if(attributes.size()) {
    if(index < attributes.size()) {
      return attributes[index];
    } else {
      cerr << "ERROR: Attribute index is out of range: " << index << endl;
      exit(1);
    }
  } else {
    cerr << "ERROR: Attempting to access attribute value when none "
            << "have been loaded" << endl;
    exit(1);
  }

}

unsigned int DatasetInstance::NumNumerics() {
  return(numerics.size());
}

double DatasetInstance::GetNumeric(unsigned int index) {
  if(numerics.size()) {
    if(index < numerics.size()) {
      return numerics[index];
    } else {
      cerr << "ERROR: Numeric index out of range: " << index << endl;
      exit(1);
    }
  } else {
    cerr << "ERROR: Attempting to access numeric value when none "
            << "have been loaded" << endl;
    exit(1);
  }
}

bool DatasetInstance::AddNumeric(NumericLevel newNum) {
  //  cout << "DatasetInstance::AddNumeric(double newNum): " << newNum << endl;
  numerics.push_back(newNum);
  return true;
}

ClassLevel DatasetInstance::GetClass() {
  return classLabel;
}

void DatasetInstance::SetClass(ClassLevel classValue) {
  classLabel = classValue;
}

double DatasetInstance::GetPredictedValueTau() {
  return predictedValueTau;
}

void DatasetInstance::SetPredictedValueTau(double newValue) {
  predictedValueTau = newValue;
}

double DatasetInstance::GetInfluenceFactorD(unsigned int neighborIndex) {
  return neighborInfluenceFactorDs[neighborIndex];
}

void DatasetInstance::ClearInfluenceFactors() {
  neighborInfluenceFactorDs.clear();
}

bool DatasetInstance::AddInfluenceFactorD(double factor) {
  neighborInfluenceFactorDs.push_back(factor);
  return true;
}

void DatasetInstance::Print() {
  vector<AttributeLevel>::const_iterator it = attributes.begin();
  for(; it != attributes.end(); ++it) {
    cout << *it << " ";
  }
  if(numerics.size()) {
    cout << " | ";
    vector<double>::const_iterator dit = numerics.begin();
    for(; dit != numerics.end(); ++dit) {
      cout << *dit << " ";
    }
  }
  // for RReliefF - bcw - 9/30/11
  if(dataset->HasContinuousPhenotypes()) {
    cout << "=> [" << predictedValueTau << "]" << endl;
  } else {
    cout << "=> [" << classLabel << "]" << endl;
  }
}

bool DatasetInstance::SwapAttributes(unsigned int a1, unsigned int a2) {
  if(a1 >= attributes.size()) {
    return false;
  }
  if(a2 >= attributes.size()) {
    return false;
  }
  // hahaha
  if(a1 == a2) {
    return true;
  }

  AttributeLevel temp = attributes[a1];
  attributes[a1] = attributes[a2];
  attributes[a2] = temp;

  return true;
}

void DatasetInstance::SetDistanceSums(unsigned int kNearestNeighbors,
                                      DistancePairs& sameClassSums,
                                      map<ClassLevel, DistancePairs>& diffClassSums) {
  // added 9/22/11 for iterative Relief-F and EC
  bestNeighborIdsSameClass.clear();
  bestNeighborIdsDiffClass.clear();

  // use Nate's best_n.h algorithm
  // cout << "Same class sums:" << endl;
  // PrintDistancePairs(sameClassSum);
  DistancePairs bestInstancesHits;
  best_n(sameClassSums.begin(), sameClassSums.end(),
         back_insert_iterator<DistancePairs > (bestInstancesHits),
         kNearestNeighbors, deref_less_bcw());
  // cout << "Hits:" << endl;
  DistancePairsIt hit;
  for(hit = bestInstancesHits.begin(); hit != bestInstancesHits.end(); ++hit) {
    DistancePair thisHit = *hit;
    // cout << thisHit.first << " => " << thisHit.second << endl;
    bestNeighborIdsSameClass.push_back(thisHit.second);
  }

  map<ClassLevel, DistancePairs>::const_iterator it = diffClassSums.begin();
  for(; it != diffClassSums.end(); ++it) {
    ClassLevel thisClass = it->first;
    DistancePairs thisDiffSums = it->second;
    DistancePairs bestInstancesMisses;
    best_n(thisDiffSums.begin(), thisDiffSums.end(),
           back_insert_iterator<DistancePairs > (bestInstancesMisses),
           kNearestNeighbors, deref_less_bcw());
    DistancePairsIt mit;
    // cout << "Class " << thisClass << ", Different class sums:" << endl;
    // PrintDistancePairs(bestInstanceMisses);
    for(mit = bestInstancesMisses.begin(); mit != bestInstancesMisses.end(); ++mit) {
      DistancePair thisMiss = *mit;
      // cout << thisMiss.first << " => " << thisMiss.second << endl;
      bestNeighborIdsDiffClass[thisClass].push_back(thisMiss.second);
    }
  }
  // cout << "----------------------------------------------------------" << endl;
}

void DatasetInstance::SetDistanceSums(unsigned int kNearestNeighbors,
                                      DistancePairs instanceSums) {
  bestNeighborIds.clear();

  //  cout << "Instance sums:" << endl;
  //  PrintDistancePairs(instanceSums);

  // use Nate's best_n.h algorithm to select the nearest neighbors
  DistancePairs bestInstances;
  best_n(instanceSums.begin(), instanceSums.end(),
         back_insert_iterator<DistancePairs > (bestInstances),
         kNearestNeighbors, deref_less_bcw());
  // sort(bestInstances.begin(), bestInstances.end());
  //  cout << "Best instances:" << endl;
  //  PrintDistancePairs(bestInstances);

  DistancePairsIt it = bestInstances.begin();
  for(; it != bestInstances.end(); ++it) {
    //    cout << it->first << " => " << it->second << endl;
    bestNeighborIds.push_back(it->second);
  }
  //  PrintVector(bestNeighborIds, "Best neighbor IDs");
  //  cout << "----------------------------------------------------------" << endl;
}

void DatasetInstance::PrintDistancePairs(const DistancePairs& distPairs) {
  for(DistancePairsIt dpit = distPairs.begin(); dpit != distPairs.end(); ++dpit) {
    cout << (*dpit).first << "\t" << (*dpit).second << endl;
  }
}

bool DatasetInstance::GetNNearestInstances(unsigned int n,
                                           vector<unsigned int>& sameClassInstances,
                                           vector<unsigned int>& diffClassInstances) {
  if((bestNeighborIdsSameClass.size() < n) || (bestNeighborIdsDiffClass.size() < n)) {
    cerr << endl << "ERROR: GetNNearestInstances: N: [" << n
            << "] is larger than the number of neighbors: "
            << "Same: " << bestNeighborIdsSameClass.size()
            << ", Different: " << bestNeighborIdsDiffClass.size() << endl;
    return false;
  }

  sameClassInstances.clear();
  diffClassInstances.clear();
  for(unsigned int i = 0; i < n; ++i) {
    unsigned int sameIdx;
    dataset->GetInstanceIndexForID(bestNeighborIdsSameClass[i], sameIdx);
    sameClassInstances.push_back(sameIdx);
    unsigned int diffIdx;
    dataset->GetInstanceIndexForID(bestNeighborIdsSameClass[i], diffIdx);
    sameClassInstances.push_back(diffIdx);
  }

  return true;
}

bool
DatasetInstance::GetNNearestInstances
(
 unsigned int n,
 vector<unsigned int>& sameClassInstances,
 map<ClassLevel, vector<unsigned int> >& diffClassInstances
 ) {

  if(bestNeighborIdsSameClass.size() < n) {
    cerr << endl << "ERROR: GetNNearestInstances: N: [" << n
            << "] is larger than the number of neighbors "
            << " in same class: " << bestNeighborIdsSameClass.size() << endl;
    return false;
  }

  sameClassInstances.clear();
  for(unsigned int i = 0; i < n; ++i) {
    unsigned int sameIdx;
    dataset->GetInstanceIndexForID(bestNeighborIdsSameClass[i], sameIdx);
    //    cout << bestNeighborIdsSameClass[i] << ", " << sameIdx << endl;
    sameClassInstances.push_back(sameIdx);
  }
  //  cout << "------" << endl;
  // diffClassInstances.clear();
  map<ClassLevel, std::vector<std::string> >::const_iterator it;
  for(it = bestNeighborIdsDiffClass.begin();
      it != bestNeighborIdsDiffClass.end(); ++it) {
    ClassLevel thisClass = it->first;
    vector<string> ids = it->second;
    if(ids.size() < n) {
      cerr << endl << "ERROR: GetNNearestInstances: N: [" << n
              << "] is larger than the number of neighbors for class "
              << thisClass << ": " << bestNeighborIdsDiffClass.size() << endl;
      return false;
    }
    for(unsigned int i = 0; i < n; ++i) {
      unsigned int diffIdx;
      dataset->GetInstanceIndexForID(ids[i], diffIdx);
      //      cout << ids[i] << ", " << diffIdx << endl;
      diffClassInstances[thisClass].push_back(diffIdx);
    }
  }
  return true;
}

bool DatasetInstance::GetNNearestInstances(unsigned int n,
                                           vector<unsigned int>& closestInstances) {

  if(bestNeighborIds.size() < n) {
    cerr << "ERROR: GetNNearestInstances: k: [" << n
            << "] is larger than the number of neighbors" << endl;
    return false;
  }

  //  cout << "Same sums (" << sameSums.size() << ")" << endl;
  //  copy(neighborSums.begin(), neighborSums.end(), ostream_iterator<double>(cout, "\n"));

  closestInstances.clear();
  for(unsigned int i = 0; i < n; ++i) {
    unsigned int instanceIndex;
    dataset->GetInstanceIndexForID(bestNeighborIds[i], instanceIndex);
    closestInstances.push_back(instanceIndex);
  }

  return true;
}

