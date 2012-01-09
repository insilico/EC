/**
 * \file example2.cpp
 * \author billwhite
 *
 * Created on December 21, 2011, 11:21 AM
 *
 * Compile with: make example2
 */

#include <iostream>
#include <string>
#include <vector>

/// Our Dataset class is in the parent directory.
#include "../Dataset.h"
#include "../DatasetInstance.h"

using namespace std;

/** 
 * Print the genotype frequencies on the console.
 */
int main(int argc, char** argv) {

  /// Create a pointer to default data set object.
  Dataset* example1Dataset = new Dataset();

  /// Load the data set with the example1.txt file - other parameters will
  /// be explained in a later tutorial example.
  string example1Filename("example1.txt");
  vector<string> ids;
  if(!example1Dataset->LoadDataset(example1Filename, false, "", "", ids)) {
    cerr << "ERROR: Could not load data set." << endl;
    exit(1);
  }

  /// Loop through the attribute values to get genotype counts.
  cout << endl << "Dataset attribute level counts and frequencies:" << endl;
  unsigned int numInstances = example1Dataset->NumInstances();
  for(unsigned int attributeIndex=0;
      attributeIndex < example1Dataset->NumAttributes();
      ++attributeIndex) {

    /// Get the attribute values for all instances at the attribute index.
    vector<AttributeLevel> attributeValues;
    if(!example1Dataset->GetAttributeValues(attributeIndex, attributeValues)) {
      cerr << "Could not get attribute values for attribute index: "
              << attributeIndex << endl;
      exit(1);
    }

    /// Loop through the attribute values, updating the attribute counts.
    map<AttributeLevel, unsigned int> attributeCounts;
    vector<AttributeLevel>::const_iterator countsIt = attributeValues.begin();
    for(; countsIt != attributeValues.end(); ++countsIt) {
      ++attributeCounts[*countsIt];
    }

    /// Print the attribute counts and frequencies.
    map<AttributeLevel, unsigned int>::const_iterator mapIt;
    for(mapIt = attributeCounts.begin();
        mapIt != attributeCounts.end();
        ++mapIt) {
      AttributeLevel thisGenotype = mapIt->first;
      unsigned int thisCount = mapIt->second;
      cout << "Attribute index: " << attributeIndex
              << ", genotype: " << thisGenotype
              << ", count: " << thisCount
              << ", frequency: " << (thisCount / (double) numInstances)
              << endl;
    }
    cout << endl;

  }

  /// Sanity check to see if level counts match.
  example1Dataset->PrintLevelCounts();

  return 0;
}

