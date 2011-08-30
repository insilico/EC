/* 
 * File:   EvaporativeCooling.cpp
 * Author: billwhite
 * 
 * Created on July 14, 2011, 9:25 PM
 */

#include "EvaporativeCooling.h"
#include "../cpprelieff/Dataset.h"
#include <boost/program_options.hpp>
#include <omp.h>

using namespace std;
namespace po = boost::program_options;

EvaporativeCooling::EvaporativeCooling(Dataset* ds, po::variables_map& vm) {
  cout << "\tEvaporative Cooling initialization:" << endl;
  if(ds) {
    dataset = ds;
  } else {
    cerr << "ERROR: dataset is not initialized." << endl;
    exit(-1);
  }
  if(vm.count("iter-remove-n")) {
    removePerIteration = vm["iter-remove-n"].as<unsigned int>();
    if((removePerIteration < 1) ||
       (removePerIteration >= dataset->NumAttributes())) {
      cerr << "ERROR: Number to remove per iteratopn ["
              << removePerIteration << "] not in valid range." << endl;
      exit(-1);
    }
   cout << "\t\tIteratively removing " << removePerIteration << endl;
  }
  if(vm.count("iter-remove-percent")) {
    double percentage = vm["iter-remove-percent"].as<unsigned int>() / 100.0;
    removePerIteration = (unsigned int)
            ((double) dataset->NumAttributes() * percentage + 0.5);
    if((removePerIteration < 1) ||
       (removePerIteration >= dataset->NumAttributes())) {
      cerr << "ERROR: Number to remove per iteratopn ["
              << removePerIteration << "] not in valid range." << endl;
      exit(-1);
    }
    cout << "\t\tIteratively removing " << (percentage * 100)
            << "% = " << removePerIteration << endl;
  }
  cout << "\t\t" << omp_get_num_procs() << " OpenMP processors available" << endl;
}

EvaporativeCooling::~EvaporativeCooling() {
}

bool EvaporativeCooling::GetECScores(std::map<std::string, double>& scores) {
  // run random jungle (76)
  // get attribute importance scores to map (257)
  // store normalize scores to map igMap (282)

  // declare and initialize iterative relieff runtime variables (150)
  // while number of attributes grater than the target number (214)
    // run relieff on current set of attributes (220)
    // get relieff ranked attributes (238)
    // find max and min (256)
    // store scaled attribute values in map rfMap (299)

      // do local search for best attribute based on classification accuracy
      // search around a best temp estimate
          // Given at trial temperature, trialTemp
          // 1. compute the information free energy (389)
          // 2. remove the attribute with the lowest free energy
          // 3. compute classification accuracy of remaining attributes
          // 4. pick trial temperature (and attribute to remove) that
          // gives best accuracy

  // set return values in the passed map reference

  return true;
}