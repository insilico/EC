/**
 * \file example1.cpp
 * \author billwhite
 *
 * Created on December 20, 2011, 11:37 PM
 *
 * Compile with: make example1
 */

#include<cstdlib>
#include <iostream>
#include <string>
#include <vector>

/// Our Dataset class is in the parent directory.
#include "ec/Dataset.h"
#include "ec/DatasetInstance.h"

using namespace std;

/** 
 * Print a data set to the console.
 */
int main(int argc, char** argv) {

  /// Create a pointer to default data set object.
  Dataset* example1Dataset = new Dataset();

  /// Load the data set with the example1.txt file - other parameters will
  /// be explained in a later tutorial example.
  string example1Filename("example1.txt");
  vector<string> ids;
  if(!example1Dataset->LoadDataset(example1Filename, "", "", ids)) {
    cerr << "ERROR: Could not load data set." << endl;
    exit(EXIT_FAILURE);
  }

  /// Loop through the instances printing the attributes and class for each.
  for(unsigned int instanceIndex = 0;
      instanceIndex < example1Dataset->NumInstances();
      ++instanceIndex) {
    DatasetInstance* thisInstance = example1Dataset->GetInstance(instanceIndex);
    for(unsigned int attributeIndex=0;
        attributeIndex < thisInstance->NumAttributes();
        ++attributeIndex) {
      cout << thisInstance->GetAttribute(attributeIndex) << " ";
    }
    cout << thisInstance->GetClass() << endl;
  }

  return 0;
}

