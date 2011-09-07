/* 
 * File:   EvaporativeCooling.cpp
 * Author: billwhite
 * 
 * Created on July 14, 2011, 9:25 PM
 */

#include <cstdlib>
#include <iostream>
#include <iomanip>

#include <boost/program_options.hpp>
#include <omp.h>
#include <gsl/gsl_rng.h>

#include "librjungle.h"
#include "RJunglePar.h"
#include "RJungleCtrl.h"
#include "DataFrame.h"
#include "FittingFct.h"
#include "RJungleHelper.h"

#include "EvaporativeCooling.h"
#include "../cpprelieff/Dataset.h"
#include "../cpprelieff/StringUtils.h"

using namespace std;
namespace po = boost::program_options;
using namespace insilico;

EvaporativeCooling::EvaporativeCooling(Dataset* ds, po::variables_map& vm) {
  cout << "\t\tEvaporative Cooling initialization:" << endl;
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
  cout << "\t\t\t" << omp_get_num_procs() << " OpenMP processors available" << endl;

  // Relief-F parameters
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

}

EvaporativeCooling::~EvaporativeCooling() {
}

bool EvaporativeCooling::ComputeECScores() {
  // run random jungle (76)
  // get attribute importance scores to map (257)
  // store normalize scores to map igMap (282)
  // set return values in the passed map reference

  cout << "\t\tRunning EC algorithm..." << endl;

  cout << "\t\t\tRunning Random Jungle..." << endl;
  std::vector<uli_t> *colMaskVec = NULL;
  time_t start, end;
  clock_t startgrow, endgrow;
  uli_t ntree = 1000;

  // fill in the parameters object for the RJ run
  rjParams.rng = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rjParams.rng, rjParams.seed);
  rjParams.ntree = ntree;
  rjParams.nrow = dataset->NumInstances();
  rjParams.depVarName = (char *) "Class";
  rjParams.verbose_flag = true;
  rjParams.filename = (char*) "";
  if(dataset->HasNumerics()) {
    rjParams.outprefix = (char *) dataset->GetNumericsFilename().c_str();
    rjParams.ncol = dataset->NumNumerics() + 1;
  }
  else {
    rjParams.outprefix = (char *) dataset->GetSnpsFilename().c_str();
    rjParams.ncol = dataset->NumAttributes() + 1;
  }
  rjParams.depVar = rjParams.ncol - 1;
  rjParams.depVarCol = rjParams.ncol - 1;
  string outPrefix(rjParams.outprefix);
  string importanceFilename = outPrefix + ".importance";
  RJungleIO io;
  io.open(rjParams);

  if(dataset->HasNumerics()) {
    // numeric data
    cout << "\t\t\t\tPreparing numeric version of Random Jungle." << endl;
    time(&start);

    // load data frame
    cout << "\t\t\t\tLoading RJ DataFrame with double values." << endl;
    DataFrame<double>* data = new DataFrame<double>(rjParams);
    data->setDim(rjParams.nrow, rjParams.ncol);
    vector<string> numericNames = dataset->GetNumericsNames();
    numericNames.push_back(rjParams.depVarName);
    data->setVarNames(numericNames);
    data->setDepVarName(rjParams.depVarName);
    data->setDepVar(rjParams.depVarCol);
    data->initMatrix();
    for(unsigned int i=0; i < dataset->NumInstances(); ++i) {
      for(unsigned int j=0; j < dataset->NumNumerics(); ++j) {
        data->set(i, j, dataset->GetInstance(i)->GetNumeric(j));
      }
      data->set(i, rjParams.depVarCol, (double) dataset->GetInstance(i)->GetClass());
    }
    data->storeCategories();
    data->makeDepVecs();
    data->getMissings();

    RJungleGen<double> rjGen;
    rjGen.init(rjParams, *data);

    startgrow = clock();
    TIMEPROF_START("RJungleCtrl~~RJungleCtrl::autoBuildInternal");
    // create controller
    RJungleCtrl<double> rjCtrl;
    cout << "\t\t\t\tRunning Random Jungle" << endl;
    rjCtrl.autoBuildInternal(rjParams, io, rjGen, *data, colMaskVec);
    TIMEPROF_STOP("RJungleCtrl~~RJungleCtrl::autoBuildInternal");
    endgrow = clock();

    // print info stuff
    RJungleHelper<double>::printRJunglePar(rjParams, *io.outLog);

    // clean up
    if(data != NULL) {
      delete data;
    }
    if (colMaskVec != NULL) {
      delete colMaskVec;
    }

    time(&end);
  }
  else {
    // SNP data
    if(dataset->HasGenotypes()) {
      cout << "\t\t\tPreparing discrete version of Random Jungle" << endl;
      time(&start);

      // load data frame
      cout << "\t\t\t\tLoading RJ DataFrame with SNP (character) values." << endl;
      DataFrame<int>* data = new DataFrame<int>(rjParams);
      data->setDim(rjParams.nrow, rjParams.ncol);
      vector<string> attributeNames = dataset->GetAttributeNames();
      attributeNames.push_back(rjParams.depVarName);
      data->setVarNames(attributeNames);
      data->setDepVarName(rjParams.depVarName);
      data->setDepVar(rjParams.depVarCol);
      data->initMatrix();
      for(unsigned int i=0; i < dataset->NumInstances(); ++i) {
        for(unsigned int j=0; j < dataset->NumAttributes(); ++j) {
          data->set(i, j, (int) dataset->GetInstance(i)->GetAttribute(j));
        }
        data->set(i, rjParams.depVarCol, (int) dataset->GetInstance(i)->GetClass());
      }
      data->storeCategories();
      data->makeDepVecs();
      data->getMissings();

      RJungleGen<int> rjGen;
      rjGen.init(rjParams, *data);

      startgrow = clock();
      TIMEPROF_START("RJungleCtrl~~RJungleCtrl::autoBuildInternal");
      // create controller
      RJungleCtrl<int> rjCtrl;
      cout << "\t\t\t\tRunning Random Jungle" << endl;
      rjCtrl.autoBuildInternal(rjParams, io, rjGen, *data, colMaskVec);
      TIMEPROF_STOP("RJungleCtrl~~RJungleCtrl::autoBuildInternal");
      endgrow = clock();

      // clean up
      if(data != NULL) {
        delete data;
      }
      if(colMaskVec != NULL) {
        delete colMaskVec;
      }

      time(&end);
    }
    else {
      cerr << "ERROR: Dataset is no loaded or of unknown data type." << endl;
      return false;
    }
  }

  // print info stuff
  RJungleHelper<double>::printRJunglePar(rjParams, *io.outLog);

  
  RJungleHelper<double>::printFooter(rjParams, io, start, end,
                                     startgrow, endgrow);

#ifdef HAVE_TIMEPROF
  timeProf.writeToFile(*io.outTimeProf);
#endif

  // clean up Random Jungle run
  io.close();
  gsl_rng_free(rjParams.rng);

  // set RJ scores for use in EC
  cout << "\t\t\t\tLoading RJ variable importance (VI) scores" << endl;
  ifstream importanceStream(importanceFilename.c_str());
    if(!importanceStream.is_open()) {
    cerr << "ERROR: Could not open Random Jungle importnace file: "
            << importanceFilename << endl;
    return false;
  }
  string line;
  // strip the header line
  getline(importanceStream, line);
  // read and store variable name and gini index
  unsigned int lineNumber = 0;
  while(getline(importanceStream, line)) {
    ++lineNumber;
    vector<string> tokens;
    split(tokens, line);
    rjScores[tokens[2]] = strtod(tokens[3].c_str(), NULL);
  }
  importanceStream.close();
  
  cout << "\t\t\tRandom Jungle finished." << endl;
  
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
