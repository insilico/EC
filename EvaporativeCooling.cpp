/* 
 * File:   EvaporativeCooling.cpp
 * Author: billwhite
 * 
 * Created on July 14, 2011, 9:25 PM
 *
 * Implements the Evaporative Cooling algorithm in:
 * McKinney, et. al. "Capturing the Spectrum of Interaction Effects in Genetic
 * Association Studies by Simulated Evaporative Cooling Network Analysis."
 * PLoS Genetics, Vol 5, Issue 3, 2009.
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
#include "Helper.h"

#include "EvaporativeCooling.h"
#include "../cpprelieff/Dataset.h"
#include "../cpprelieff/StringUtils.h"
#include "../cpprelieff/ReliefF.h"

using namespace std;
namespace po = boost::program_options;
using namespace insilico;

bool freeEnergyScoreSort(const pair<string, double>& p1,
                         const pair<string, double>& p2) {
  return p1.second < p2.second;
}

/*****************************************************************************
 * Method: constructor
 *
 * IN:  pointer to a Dataset or subclass, reference to program options map
 * OUT: none
 *
 * Initialize an EvaporativeCooling object.
 ****************************************************************************/
EvaporativeCooling::EvaporativeCooling(Dataset* ds, po::variables_map& vm) {
  cout << "\t\tEvaporative Cooling initialization:" << endl;
  if(ds) {
    dataset = ds;
  } else {
    cerr << "ERROR: dataset is not initialized." << endl;
    exit(-1);
  }

  paramsMap = vm;

  numTargetAttributes = vm["ec-num-target"].as<unsigned int>();
  if(numTargetAttributes < 1) {
    cerr << "Use --ec-num-target to the number of cooled attributes desired." << endl;
    exit(-1);
  }
  if(numTargetAttributes > dataset->NumAttributes()) {
    cerr << "--ec-num-taget must be less than or equal to the "
            << "number of attributes in the data set." << endl;
    exit(-1);
  }
 
  rjParams.verbose_flag = vm["verbose"].as<bool>();

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


}

/*****************************************************************************
 * Method: destructor
 *
 * IN:  none
 * OUT: none
 *
 * Destroy object.
 ****************************************************************************/
EvaporativeCooling::~EvaporativeCooling() {
}

/*****************************************************************************
 * Method: ComputeECScores
 *
 * IN:  none
 * OUT: success
 *
 * Runs the Evaporative Cooling algorithm.
 ****************************************************************************/
bool EvaporativeCooling::ComputeECScores() {
  // initialize EC algorithm working variables
  unsigned int iteration = 0;

  // all attributes are considered in Step 0
  unsigned int numWorkingAttributes = dataset->NumAttributes();
  vector<string> attributeNames = dataset->GetAttributeNames();
  vector<string>::const_iterator it = attributeNames.begin();
  for(; it != attributeNames.end(); ++it) {
    attributesToConsider[*it] = true;
  }
  
  // Relief-F algorithm
  reliefF = new ReliefF(dataset, paramsMap);

  // EC algorithm as in Figure 5, page 10 of the paper referenced
  // at top of this file. Modified per Brett's email to not do the
  // varying temperature and classifier accuracy optimization steps.
  while(numWorkingAttributes > numTargetAttributes) {
    cout << "\t\tRunning EC algorithm...iteration: " << iteration << endl;

    // -------------------------------------------------------------------------
    // run Random Jungle and get the normalized scores for use in EC
    cout << "\t\t\tRunning Random Jungle..." << endl;
    if(!RunRandomJungle()) {
      cerr << "Random Jungle failed." << endl;
      return false;
    }
    cout << "\t\t\tRandom Jungle finished." << endl;

    // -------------------------------------------------------------------------
    // run Relief-F and get normalized score for use in EC
    cout << "\t\t\tRunning ReliefF..." << endl;
    if(!RunReliefF()) {
      cerr << "ReliefF failed." << endl;
      return false;
    }
    cout << "\t\t\tReliefF finished." << endl;

    // -------------------------------------------------------------------------
    // compute free energy for all attributes
    cout << "\t\t\tComputing free energy..." << endl;
    double temperature = 1.0;
    if(!ComputeFreeEnergy(temperature)) {
      cerr << "ComputeFreeEnergy failed." << endl;
      return false;
    }

    // -------------------------------------------------------------------------
    // remove the worst attribute and iterate
    cout << "\t\t\tRemoving the worst attribute..." << endl;
    if(!RemoveWorstAttributes()) {
      cerr << "RemoveWorstAttribute failed." << endl;
      return false;
    }

    // go to next iteration
    --numWorkingAttributes;
    ++iteration;
  }

  // remaining free energy attributes are the ones we want to write as a
  // new dataset to be analyzed with (re)GAIN + SNPrank

  // clean up
  delete reliefF;
  reliefF = NULL;
  
  return true;
}

/*****************************************************************************
 * Method: RemoveWorstAttribute
 *
 * IN:  none
 * OUT: success
 *
 * Remove the worst of the attributes based on free energy score.
 ****************************************************************************/
bool EvaporativeCooling::RemoveWorstAttributes() {
  // find the minumum score
  EcScoresMapCIt it = freeEnergyScores.begin();
  string minAttr = it->first;
  double minScore = it->second;
  for(; it != freeEnergyScores.end(); ++it) {
    if(it->second < minScore) {
      minAttr = it->first;
      minScore = it->second;
    }
  }

  // save worst
  ecScores.insert(make_pair(minAttr, minScore));

  // remove the attribute from those under consideration
  attributesToConsider.erase(minAttr);
  if(!reliefF->RemoveAttributeFromDistanceCalc(minAttr)) {
    cerr << "EC could not remove attribute: " << minAttr 
            << " from Relief-F." << endl;
    return false;
  }
//  if(!dataset->RemoveAttribute(minAttr)) {
//    cerr << "EC could not remove attribute: " << minAttr
//            << " from Dataset." << endl;
//    return false;
//  }

  int iterNumToRemove = 0;
  int iterPercentToRemove = 0;
  if(paramsMap.count("iter-remove-n")) {
    iterNumToRemove = paramsMap["iter-remove-n"].as<int>();
  }
  if(paramsMap.count("iter-remove-percent")) {
    iterPercentToRemove = paramsMap["iter-remove-percent"].as<int>();
  }
  vector<double> scores;
  if((iterNumToRemove == 0) && (iterPercentToRemove == 0)) {
    reliefF->PreComputeDistances();
  } else {
    reliefF->PreComputeDistancesIterative();
  }

  return true;
}

/*****************************************************************************
 * Method: ComputeFreeEnergy
 *
 * IN:  none
 * OUT: success
 *
 * Compute the free energey for each attribute based on Random Jungle
 * and ReliefF attribute scores. F = E - TS
 * where E=Relief-F, S=Random Jungle, T=temperature
 ****************************************************************************/
bool EvaporativeCooling::ComputeFreeEnergy(double temperature) {
  if(rjScores.size() != rfScores.size()) {
    cerr << "EvaporativeCooling::ComputeFreeEnergy scores lists are " 
            "unequal. RJ: " << rjScores.size() << " vs. RF: " <<
            rfScores.size() << endl;
    return false;
  }
  freeEnergyScores.clear();
  EcScoresMapCIt rjIt = rjScores.begin();
  for(; rjIt != rjScores.end(); ++rjIt) {
    string key = rjIt->first;
    double val = rjIt->second;
    freeEnergyScores.insert(make_pair(key, rfScores[key] - (temperature * val)));
  }
  return true;
}

/*****************************************************************************
 * Method: RunReliefF
 *
 * IN:  none
 * OUT: success
 *
 * Run the ReliefF algorithm to get interaction ranked variables.
 ****************************************************************************/
bool EvaporativeCooling::RunReliefF() {
  int iterNumToRemove = 0;
  int iterPercentToRemove = 0;
  if(paramsMap.count("iter-remove-n")) {
    iterNumToRemove = paramsMap["iter-remove-n"].as<int>();
  }
  if(paramsMap.count("iter-remove-percent")) {
    iterPercentToRemove = paramsMap["iter-remove-percent"].as<int>();
  }
  if((iterNumToRemove == 0) && (iterPercentToRemove == 0)) {
    cout << "\t\t\t\tRunning standard ReliefF..." << endl;
    rfScores = reliefF->ComputeAttributeScores();
  } else {
    // determine the number to remove per iteration
    cout << "\t\t\t\tRunning Iterative ReliefF..." << endl;
    rfScores = reliefF->ComputeAttributeScoresIteratively();
  }
  double minRFScore = rfScores.begin()->second;
  double maxRFScore = rfScores.begin()->second;
  EcScoresMapCIt rfScoresIt = rfScores.begin();
  for(; rfScoresIt != rfScores.end(); ++rfScoresIt) {
    if(rfScoresIt->second < minRFScore) {
      minRFScore = rfScoresIt->second;
    }
    if(rfScoresIt->second > maxRFScore) {
      maxRFScore = rfScoresIt->second;
    }
  }

  // get and store normalized attribute values
  bool needsNormalization = true;
  if(minRFScore == maxRFScore) {
    cerr << "WARNING: Relief-F min and max scores are the same." << endl;
    needsNormalization = false;;
  }
  double rfRange = maxRFScore - minRFScore;
  for(EcScoresMapIt it = rfScores.begin(); it != rfScores.end(); ++it) {
    string key = (*it).first;
    double val = (*it).second;
    if(needsNormalization) {
      rfScores[key] = (val - minRFScore) / rfRange;
    }
    else {
      rfScores[key] = val;
    }
  }

  return true;
}

/*****************************************************************************
 * Method: RunRandomJungle
 *
 * IN:  none
 * OUT: success
 *
 * Run the Random Jungle algorithm to get main effects ranked variables.
 ****************************************************************************/
bool EvaporativeCooling::RunRandomJungle() {
  vector<uli_t>* colMaskVec = NULL;
  time_t start, end;
  clock_t startgrow, endgrow;
  uli_t ntree = 1000;

  // fill in the parameters object for the RJ run
  rjParams.rng = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rjParams.rng, rjParams.seed);


  if(paramsMap.count("rj-num-trees")) {
    rjParams.ntree = paramsMap["rj-num-trees"].as<uli_t>();
  } else {
    rjParams.ntree = ntree;
  }

  rjParams.nrow = dataset->NumInstances();
  rjParams.depVarName = (char *) "Class";
  //  rjParams.verbose_flag = true;
  rjParams.filename = (char*) "";
  if(dataset->HasNumerics()) {
    rjParams.outprefix = (char *) dataset->GetNumericsFilename().c_str();
    rjParams.ncol = dataset->NumNumerics() + 1;
  } else {
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
    for(unsigned int i = 0; i < dataset->NumInstances(); ++i) {
      for(unsigned int j = 0; j < dataset->NumNumerics(); ++j) {
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
    if(colMaskVec != NULL) {
      delete colMaskVec;
    }

    time(&end);
  } else {
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
      for(unsigned int i = 0; i < dataset->NumInstances(); ++i) {
        for(unsigned int j = 0; j < dataset->NumAttributes(); ++j) {
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
    } else {
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

  cout << "\t\t\t\tLoading RJ variable importance (VI) scores" << endl;
  if(!ReadRandomJungleScores(importanceFilename)) {
    cerr << "Could not read Random Jungle scores." << endl;
    return false;
  }

  return true;
}

/*****************************************************************************
 * Method: ReadRandomJungleScores
 *
 * IN:  filename of the ranked attributes from Random Jungle
 * OUT: success
 *
 * Read and normalizae the ranked attributes from the Random Jungle run.
 * Side effect: sets rjScores map.
 ****************************************************************************/
bool EvaporativeCooling::ReadRandomJungleScores(string importanceFilename) {
  ifstream importanceStream(importanceFilename.c_str());
  if(!importanceStream.is_open()) {
    cerr << "ERROR: Could not open Random Jungle importance file: "
            << importanceFilename << endl;
    return false;
  }
  string line;
  // strip the header line
  getline(importanceStream, line);
  // read and store variable name and gini index
  unsigned int lineNumber = 0;
  double minRJScore = 0; // initializations to keep compiler happy
  double maxRJScore = 0; // initializations are on line number 1 of file read
  while(getline(importanceStream, line)) {
    ++lineNumber;
    vector<string> tokens;
    split(tokens, line, " ");
    if(tokens.size() != 4) {
      cerr << "EvaporativeCooling::ReadRandomJungleScores: error parsing line "
              << lineNumber << ". Read " << tokens.size() << " columns. Should "
              << "be 4." << endl;
      return false;
    }
    string key = tokens[2];
    double val = strtod(tokens[3].c_str(), NULL);
    rjScores[key] = val;
    if(lineNumber == 1) {
      minRJScore = val;
      maxRJScore = val;
    } else {
      if(val < minRJScore) {
        minRJScore = val;
      } else {
        if(val > maxRJScore) {
          maxRJScore = val;
        }
      }
    }
  }
  importanceStream.close();
  // normalize map scores
  bool needsNormalization = true;
  if(minRJScore == maxRJScore) {
    cerr << "WARNING: Random Jungle min and max scores are the same." << endl;
    needsNormalization = false;
  }
  double rjRange = maxRJScore - minRJScore;
  for(EcScoresMapIt it = rjScores.begin(); it != rjScores.end(); ++it) {
    string key = (*it).first;
    double val = (*it).second;
    if(needsNormalization) {
      rjScores[key] = (val - minRJScore) / rjRange;
    }
    else {
      rjScores[key] = val;
    }
  }
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

  for(EcScoresMapCIt ecScoresIt = ecScores.begin();
      ecScoresIt != ecScores.end(); ++ecScoresIt) {
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
