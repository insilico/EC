/**
 * \file example3.cpp
 * \author billwhite
 *
 * Created on December 22, 2011, 7:35 PM
 *
 * Compile with: make example3
 */

#include <iostream>
#include <string>
#include <vector>

/// Our Dataset class is in the parent directory.
#include "../Dataset.h"
#include "../DatasetInstance.h"
#include "../PlinkDataset.h"
#include "../PlinkBinaryDataset.h"
#include "../Insilico.h"

using namespace std;

/** 
 * Print the attribute allele frequencies on the console.
 */
int main(int argc, char** argv) {

  /// Use PLINK and tab-delimited formats to compare results.
  vector<string> datasetNames;
  datasetNames.push_back("example3.ped");
  datasetNames.push_back("example3.bed");
  datasetNames.push_back("example3.txt");

  /// Loop through the data set types.
  vector<string>::const_iterator fileIt = datasetNames.begin();
  for(; fileIt != datasetNames.end(); ++fileIt) {

    /// Create a pointer to a data set object by calling the library function
    /// to determine the data set class to use by its filename extension.
    Dataset* example3Dataset = ChooseSnpsDatasetByExtension(*fileIt);
    vector<string> ids;
    if(!example3Dataset->LoadDataset(*fileIt, false, "", "", ids)) {
      cerr << "ERROR: Could not load data set." << endl;
      exit(1);
    }

    /// Loop through the attribute values to get the minor allele
    /// and its frequency.
    cout << endl << "Dataset: " << *fileIt
            << ", minor allele frequencies:" << endl;
    for(unsigned int attributeIndex=0;
        attributeIndex < example3Dataset->NumAttributes();
        ++attributeIndex) {

      /// Print the attribute minor allele and its frequency.
      pair<char, double> attributeMAF =
        example3Dataset->GetAttributeMAF(attributeIndex);
      if(attributeMAF.first == ' ') {
        cout << "WARNING: Attribute Index: " << attributeIndex
                << " minor allele could not be determined." << endl;
      }
      else {
        cout << "Attribute Index: " << attributeIndex
                << ", Allele: " << attributeMAF.first
                << ", Frequency: " << attributeMAF.second << endl;
      }
    }

    /// Clean up dynamically allocated memory for the data set.
    delete example3Dataset;
  }

  return 0;
}

