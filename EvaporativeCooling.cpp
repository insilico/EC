/* 
 * File:   EvaporativeCooling.cpp
 * Author: billwhite
 * 
 * Created on July 14, 2011, 9:25 PM
 */

#include <iostream>
#include <iomanip>

#include <boost/program_options.hpp>
#include "EvaporativeCooling.h"
#include "../cpprelieff/Dataset.h"
#include <omp.h>

#include "librjungle.h"
#include "RJunglePar.h"
#include "RJungleCtrl.h"
#include "DataFrame.h"
#include "FittingFct.h"
#include "RJungleHelper.h"

#undef bool

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

  // ------------------------------------------------------------- Random Jungle

  // load Random Jungle parameters object from
  // boost program options variable map
  rjParams = initRJunglePar();
  rjParams.mpiId = 0;
  // set the number of threads that will be used
  if((rjParams.nthreads < omp_get_max_threads()) && (rjParams.nthreads != 0)) {
    omp_set_num_threads(rjParams.nthreads);
  }

  // ------------------------------------------------------------------ Relief-F

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

bool EvaporativeCooling::ComputeECScores() {
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
  std::vector<uli_t> *colMaskVec = NULL;
  time_t start, end;
  clock_t startgrow, endgrow;
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus);
  uli_t ntree = 1000;
  rjParams.rng = rng;
  rjParams.ntree = ntree;
  
  // prepare input/output
  if(strcmp(rjParams.outprefix, "") == 0) {
    rjParams.outprefix = (char *) "rjungle";
  }
  RJungleIO io;
  io.open(rjParams);

  if(dataset->HasNumerics()) {
    // numeric data
    // create controller
    RJungleCtrl<double> rjCtrl;

    time(&start);

    // load data frame
    DataFrame<double>* data = new DataFrame<double>(rjParams);
    data->setDim(dataset->NumInstances(), dataset->NumAttributes() + 1);
    data->setVarNames(dataset->GetAttributeNames());
    data->setDepVar(dataset->GetClassColumn());
    data->initMatrix();
    for(unsigned int i=0; i < dataset->NumInstances(); ++i) {
      for(unsigned int j=0; j < dataset->NumNumerics(); ++j) {
        data->set(i, j, dataset->GetInstance(i)->GetNumeric(j));
      }
      data->set(i, dataset->NumNumerics(), dataset->GetInstance(i)->GetClass());
    }

    RJungleGen<double> rjGen;
    rjGen.init(rjParams, *data);

    startgrow = clock();
    TIMEPROF_START("RJungleCtrl~~RJungleCtrl::autoBuildInternal");
    rjCtrl.autoBuildInternal(rjParams, io, rjGen, *data, colMaskVec);
    TIMEPROF_STOP("RJungleCtrl~~RJungleCtrl::autoBuildInternal");
    endgrow = clock();

    time(&end);

    RJungleHelper<double>::printFooter(rjParams, io, start, end, startgrow, endgrow);

  }
  else {
    // SNP data
    if(dataset->HasGenotypes()) {
      
    }
    else {
      cerr << "ERROR: Dataset is no loaded or of unknown data type." << endl;
      return false;
    }
  }

  io.close();

  gsl_rng_free(rjParams.rng);

  return true;
}

/*****************************************************************************
 * Method: PrintAttributeScores
 *
 * IN:  output file stresm onto which scores will be sent
 * OUT: none
 *
 * Send score and attribute name to output filestream passed.
 ****************************************************************************/
void EvaporativeCooling::PrintAttributeScores(ofstream& outFile) {
  vector<string> attributeNames = dataset->GetAttributeNames();
  map<string, double>::const_iterator ecScoresIt;
  for(ecScoresIt = ecScores.begin(); ecScoresIt != ecScores.end(); ++ecScoresIt) {
      outFile << fixed << setprecision(8) << (*ecScoresIt).second << "\t"
              << (*ecScoresIt).first << endl;
  }
}

/*****************************************************************************
 * Method: WriteAttributeScores
 *
 * IN:  output file name into which scores will be printed
 * OUT: none
 *
 * Send score and attribute name to output filename passed.
 ****************************************************************************/
void EvaporativeCooling::WriteAttributeScores(string baseFilename) {
  string resultsFilename = baseFilename + ".ec";
  ofstream outFile;
  outFile.open(resultsFilename.c_str());
  if(outFile.bad()) {
    cerr << "Could not open scores file " << resultsFilename
            << "for writing." << endl;
    exit(1);
  }
  PrintAttributeScores(outFile);
  outFile.close();
}
