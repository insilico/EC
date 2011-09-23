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
#include "../cpprelieff/Statistics.h"
#include "../cpprelieff/StringUtils.h"
#include "../cpprelieff/ReliefF.h"

using namespace std;
namespace po = boost::program_options;
using namespace insilico;

bool scoresSortAsc(const pair<double, string>& p1,
                   const pair<double, string>& p2) {
  return p1.first < p2.first;
}

bool scoresSortAscByName(const pair<double, string>& p1,
                         const pair<double, string>& p2) {
  return p1.second < p2.second;
}

bool scoresSortDesc(const pair<double, string>& p1,
                    const pair<double, string>& p2) {
  return p1.first > p2.first;
}

/*****************************************************************************
 * Method: constructor
 *
 * IN:  pointer to a Dataset or subclass, reference to program options map
 * OUT: none
 *
 * Initialize an EvaporativeCooling object.
 ****************************************************************************/
EvaporativeCooling::EvaporativeCooling(Dataset* ds, po::variables_map& vm,
                                       AnalysisType anaType) {
  cout << "\t\tEvaporative Cooling initialization:" << endl;
  if(ds) {
    dataset = ds;
  } else {
    cerr << "ERROR: dataset is not initialized." << endl;
    exit(-1);
  }
  paramsMap = vm;
  analysisType = anaType;

  // set the number of attributea to remove per iteration
  if(paramsMap.count("ec-iter-remove-n")) {
    numToRemovePerIteration = paramsMap["ec-iter-remove-n"].as<unsigned int>();
  }
  if(paramsMap.count("ec-iter-remove-percent")) {
    int iterPercentToRemove = paramsMap["ec-iter-remove-percent"].as<unsigned int>();
    numToRemovePerIteration = (int) (((double) iterPercentToRemove / 100.0) *
                                     dataset->NumAttributes());
  }
  cout << "\t\t\tEC will remove " << numToRemovePerIteration
          << " attributes per iteration." << endl;

  // set the number of attributea to remove per iteration - ReliefF
  rfNumToRemovePerIteration = 0;
  if(paramsMap.count("iter-remove-n")) {
    rfNumToRemovePerIteration = paramsMap["iter-remove-n"].as<unsigned int>();
  }
  if(paramsMap.count("iter-remove-percent")) {
    unsigned int iterPercentToRemove = paramsMap["iter-remove-percent"].as<unsigned int>();
    rfNumToRemovePerIteration = (int) (((double) iterPercentToRemove / 100.0) *
                                       dataset->NumAttributes());
  }
  cout << "\t\t\tRelief-F will remove " << rfNumToRemovePerIteration
          << " attributes per iteration." << endl;

  // set the number of target attributes
  numTargetAttributes = vm["ec-num-target"].as<unsigned int>();
  if(numTargetAttributes < 1) {
    cerr << "Use --ec-num-target to the number of best attributes desired."
            << endl;
    exit(-1);
  }
  if(numTargetAttributes > dataset->NumAttributes()) {
    cerr << "--ec-num-taget must be less than or equal to the "
            << "number of attributes in the data set." << endl;
    exit(-1);
  }
  cout << "\t\t\tEC is removing attributes until best " << numTargetAttributes
          << " remain." << endl;

  // multithreading setup
  unsigned int maxThreads = omp_get_num_procs();
  cout << "\t\t\t" << maxThreads << " OpenMP processors available." << endl;

  numRJThreads = vm["rj-num-threads"].as<unsigned int>();
  if((numRJThreads < 1) || (numRJThreads > maxThreads)) {
    numRJThreads = maxThreads;
  }
  cout << "\t\t\tRandom Jungle will use " << numRJThreads << " threads." << endl;

  numRFThreads = vm["rf-num-threads"].as<unsigned int>();
  if((numRFThreads < 1) || (numRFThreads > maxThreads)) {
    numRFThreads = maxThreads;
  }
  cout << "\t\t\tRelief-F will use " << numRFThreads << " threads." << endl;

  // ------------------------------------------------------------- Random Jungle
  // initialize Random Jungle
  uli_t numTrees = vm["rj-num-trees"].as<uli_t > ();
  cout << "\t\t\tInitializing Random Jungle with " << numTrees << " trees." << endl;
  if(!InitializeRandomJungle(numTrees)) {
    cerr << "ERROR: could not initialize Random Jungle." << endl;
    exit(1);
  }

  // ------------------------------------------------------------------ Relief-F
  // intialize Relief-F
  cout << "\t\t\tInitializing Relief-F." << endl;
  reliefF = new ReliefF(dataset, paramsMap, analysisType);

  // end of constructor
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
  FinalizeRandomJungle();
  if(reliefF) {
    delete reliefF;
  }
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
  unsigned int numWorkingAttributes = dataset->NumAttributes();
  if(numWorkingAttributes <= numTargetAttributes) {
    cerr << "ERROR: The number of attributes in the data set "
            << numWorkingAttributes
            << " is less than the number of target attributes "
            << numTargetAttributes << endl;
    return false;
  }

  // EC algorithm as in Figure 5, page 10 of the paper referenced
  // at top of this file. Modified per Brett's email to not do the
  // varying temperature and classifier accuracy optimization steps.
  unsigned int iteration = 1;
  while(numWorkingAttributes > numTargetAttributes) {
    cout << "\t\t----------------------------------------------------"
            << "-------------------------" << endl;
    cout << "\t\tEC algorithm...iteration: " << iteration
            << ", working attributes: " << numWorkingAttributes
            << ", target attributes: " << numTargetAttributes
            << endl;

    // -------------------------------------------------------------------------
    // run Random Jungle and get the normalized scores for use in EC
    cout << "\t\t\tRunning Random Jungle..." << endl;
    if(!RunRandomJungle()) {
      cerr << "ERROR: In EC algorithm: Random Jungle failed." << endl;
      return false;
    }
    cout << "\t\t\tRandom Jungle finished." << endl;

    // -------------------------------------------------------------------------
    // run Relief-F and get normalized score for use in EC
    cout << "\t\t\tRunning ReliefF..." << endl;
    if(!RunReliefF()) {
      cerr << "ERROR: In EC algorithm: ReliefF failed." << endl;
      return false;
    }
    cout << "\t\t\tReliefF finished." << endl;

    // -------------------------------------------------------------------------
    // compute free energy for all attributes
    cout << "\t\t\tComputing free energy..." << endl;
    double temperature = 1.0;
    if(!ComputeFreeEnergy(temperature)) {
      cerr << "ERROR: In EC algorithm: ComputeFreeEnergy failed." << endl;
      return false;
    }
    cout << "Free energy calculations complete." << endl;
    // PrintAllScoresTabular();
    // PrintKendallTaus();

    // -------------------------------------------------------------------------
    // remove the worst attributes and iterate
    cout << "\t\t\tRemoving the worst attributes..." << endl;
    unsigned int numToRemove = numToRemovePerIteration;
    if(paramsMap.count("ec-iter-remove-percent")) {
      unsigned int iterPercentToRemove =
        paramsMap["ec-iter-remove-percent"].as<unsigned int>();
      numToRemove = (int) (((double) iterPercentToRemove / 100.0) *
                                       dataset->NumAttributes());
    }
    if((numWorkingAttributes - numToRemove) < numTargetAttributes) {
      numToRemove = numWorkingAttributes - numTargetAttributes;
    }
    if(numToRemove < 1) {
      break;
    }
    cout << "\t\t\t\tRemoving the worst " << numToRemove << " attributes..." << endl;
    if(!RemoveWorstAttributes(numToRemove)) {
      cerr << "ERROR: In EC algorithm: RemoveWorstAttribute failed." << endl;
      return false;
    }
    numWorkingAttributes -= numToRemove;
    cout << "Attribute removal complete." << endl;

    ++iteration;
  }

  cout << "\t\tEC algorithm ran for " << iteration << " iterations." << endl;

  // remaining free energy attributes are the ones we want to write as a
  // new dataset to be analyzed with (re)GAIN + SNPrank
  sort(freeEnergyScores.begin(), freeEnergyScores.end(), scoresSortDesc);
  ecScores.resize(numTargetAttributes);
  copy(freeEnergyScores.begin(),
       freeEnergyScores.begin() + numTargetAttributes,
       ecScores.begin());

  return true;
}

/*****************************************************************************
 * Methods: scores getters
 *
 * Get scores from phases of the EC algorithm.
 ****************************************************************************/
EcScores& EvaporativeCooling::GetRandomJungleScores() {
  return rjScores;
}

EcScores& EvaporativeCooling::GetReliefFScores() {
  return rfScores;
}

EcScores& EvaporativeCooling::GetECScores() {
  return ecScores;
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

  for(EcScoresCIt ecScoresIt = ecScores.begin();
      ecScoresIt != ecScores.end(); ++ecScoresIt) {
    outFile << fixed << setprecision(8) << (*ecScoresIt).first << "\t"
            << (*ecScoresIt).second << endl;
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

/*****************************************************************************
 * Method: PrintAllScoresTabular
 *
 * IN:  none
 * OUT: success
 *
 * Print all scores to stdout.
 ****************************************************************************/
bool EvaporativeCooling::PrintAllScoresTabular() {
  // sanity checks
  if(rjScores.size() != rfScores.size()) {
    cerr << "Random Jungle and Relief-F scores lists are not the same size."
            << endl;
    return false;
  }
  if(freeEnergyScores.size() != rfScores.size()) {
    cerr << "Random Jungle and Relief-F scores lists are not the same size."
            << endl;
    return false;
  }

  sort(rjScores.begin(), rjScores.end(), scoresSortDesc);
  sort(rfScores.begin(), rfScores.end(), scoresSortDesc);
  sort(freeEnergyScores.begin(), freeEnergyScores.end(), scoresSortDesc);

  cout << "\t\t\tE (RF)\t\tS (RJ)\t\tF (free energy)\n";
  unsigned int numScores = freeEnergyScores.size();
  for(unsigned int i = 0; i < numScores; ++i) {
    pair<double, string> thisRJScores = rjScores[i];
    pair<double, string> thisRFScores = rfScores[i];
    pair<double, string> thisFEScores = freeEnergyScores[i];
    printf("\t\t\t%s\t%6.4f\t%s\t%6.4f\t%s\t%6.4f\n",
           thisRFScores.second.c_str(), thisRFScores.first,
           thisRJScores.second.c_str(), thisRJScores.first,
           thisFEScores.second.c_str(), thisFEScores.first);
  }

  return true;
}

/*****************************************************************************
 * Method: PrintKendallTaus
 *
 * IN:  none
 * OUT: success
 *
 * Print all ranked scores list combination Kendal Tau correlation coefficient.
 ****************************************************************************/
bool EvaporativeCooling::PrintKendallTaus() {
  // sanity checks
  if(rjScores.size() != rfScores.size()) {
    cerr << "Random Jungle and Relief-F scores lists are not the same size."
            << endl;
    return false;
  }
  if(freeEnergyScores.size() != rfScores.size()) {
    cerr << "Random Jungle and Relief-F scores lists are not the same size."
            << endl;
    return false;
  }

  sort(rjScores.begin(), rjScores.end(), scoresSortDesc);
  sort(rfScores.begin(), rfScores.end(), scoresSortDesc);
  sort(freeEnergyScores.begin(), freeEnergyScores.end(), scoresSortDesc);

  vector<string> rjNames;
  vector<string> rfNames;
  vector<string> feNames;
  unsigned int numScores = freeEnergyScores.size();
  for(unsigned int i = 0; i < numScores; ++i) {
    pair<double, string> thisRJScores = rjScores[i];
    pair<double, string> thisRFScores = rfScores[i];
    pair<double, string> thisFEScores = freeEnergyScores[i];
    rjNames.push_back(thisRJScores.second);
    rfNames.push_back(thisRFScores.second);
    feNames.push_back(thisFEScores.second);
  }

  double tauRJRF = KendallTau(rjNames, rfNames);
  double tauRJFE = KendallTau(rjNames, feNames);
  double tauRFFE = KendallTau(rfNames, feNames);

  cout << "\t\t\tKendall tau's: "
          << "RJvRF: " << tauRJRF
          << ", RJvFE: " << tauRJFE
          << ", RFvFE: " << tauRFFE
          << endl;

  return true;
}

/*****************************************************************************
 * Method: InitializeRandomJungle
 *
 * IN:  none
 * OUT: success
 *
 * Run the Random Jungle algorithm to get main effects ranked variables.
 ****************************************************************************/
bool EvaporativeCooling::InitializeRandomJungle(uli_t ntree) {
  rjParams = initRJunglePar();
  rjParams.mpiId = 0;
  rjParams.nthreads = numRJThreads;
  rjParams.verbose_flag = paramsMap["verbose"].as<bool>();

  // fill in the parameters object for the RJ run
  rjParams.rng = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(rjParams.rng, rjParams.seed);

  if(paramsMap.count("rj-num-trees")) {
    rjParams.ntree = paramsMap["rj-num-trees"].as<uli_t > ();
  } else {
    rjParams.ntree = ntree;
  }

  rjParams.nrow = dataset->NumInstances();
  rjParams.depVarName = (char *) "Class";
  //  rjParams.verbose_flag = true;
  rjParams.filename = (char*) "";

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

  // only support numeric OR SNP data
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

  unsigned int numInstances = dataset->NumInstances();

  if(dataset->HasNumerics()) {
    // numeric data
    cout << "\t\t\t\tPreparing numeric version of Random Jungle." << endl;
    time(&start);

    // load data frame
    // TODO: do not load data frame every time-- use column mask mechanism?
    cout << "\t\t\t\tLoading RJ DataFrame with double values: ";
    DataFrame<double>* data = new DataFrame<double>(rjParams);
    data->setDim(rjParams.nrow, rjParams.ncol);
    vector<string> numericNames = dataset->GetNumericsNames();
    numericNames.push_back(rjParams.depVarName);
    data->setVarNames(numericNames);
    data->setDepVarName(rjParams.depVarName);
    data->setDepVar(rjParams.depVarCol);
    data->initMatrix();
    for(unsigned int i = 0; i < numInstances; ++i) {
      for(unsigned int j = 0; j < numericNames.size() - 1; ++j) {
        data->set(i, j, dataset->GetNumeric(i, numericNames[j]));
      }
      data->set(i, rjParams.depVarCol,
                (double) dataset->GetInstance(i)->GetClass());
      // happy lights
      if(i && ((i % 100) == 0)) {
        cout << i << "/" << numInstances << " ";
        cout.flush();
      }
    }
    cout << numInstances << "/" << numInstances << endl;
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
      cout << "\t\t\t\tLoading RJ DataFrame with SNP (character) values: ";
      DataFrame<int>* data = new DataFrame<int>(rjParams);
      data->setDim(rjParams.nrow, rjParams.ncol);
      vector<string> attributeNames = dataset->GetAttributeNames();
      attributeNames.push_back(rjParams.depVarName);
      data->setVarNames(attributeNames);
      data->setDepVarName(rjParams.depVarName);
      data->setDepVar(rjParams.depVarCol);
      data->initMatrix();
      for(unsigned int i = 0; i < dataset->NumInstances(); ++i) {
        for(unsigned int j = 0; j < attributeNames.size() - 1; ++j) {
          data->set(i, j, (int) dataset->GetAttribute(i, attributeNames[j]));
        }
        data->set(i, rjParams.depVarCol, (int) dataset->GetInstance(i)->GetClass());
        if(i && ((i % 100) == 0)) {
          cout << i << "/" << numInstances << " ";
          cout.flush();
        }
      }
      cout << numInstances << "/" << numInstances << endl;
      ;
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
      cerr << "ERROR: Dataset is not loaded or of unknown data type." << endl;
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
  rjScores.clear();
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
    string val = tokens[2];
    string keyVal = tokens[3];
    double key = strtod(keyVal.c_str(), NULL);
    // cout << "Storing RJ: " << key << " => " << val << endl;
    rjScores.push_back(make_pair(key, val));
    if(lineNumber == 1) {
      minRJScore = key;
      maxRJScore = key;
    } else {
      if(key < minRJScore) {
        minRJScore = key;
      } else {
        if(key > maxRJScore) {
          maxRJScore = key;
        }
      }
    }
  }
  importanceStream.close();
  cout << "\t\t\tRead " << rjScores.size() << " scores from "
          << importanceFilename << endl;
  // normalize map scores
  bool needsNormalization = true;
  if(minRJScore == maxRJScore) {
    cerr << "WARNING: Random Jungle min and max scores are the same." << endl;
    needsNormalization = false;
  }
  double rjRange = maxRJScore - minRJScore;
  EcScores newRJScores;
  for(unsigned int i = 0; i < rjScores.size(); ++i) {
    pair<double, string> thisScore = rjScores[i];
    double key = thisScore.first;
    string val = thisScore.second;
    if(needsNormalization) {
      key = (key - minRJScore) / rjRange;
      newRJScores.push_back(make_pair(key, val));
    } else {
      newRJScores.push_back(make_pair(key, val));
    }
  }

  rjScores.clear();
  rjScores = newRJScores;

  return true;
}

/*****************************************************************************
 * Method: FinalizeRandomJungle
 *
 * IN:  none
 * OUT: success
 *
 * Clean up Random Jungle setup.
 ****************************************************************************/
bool EvaporativeCooling::FinalizeRandomJungle() {
  if(rjParams.rng) {
    gsl_rng_free(rjParams.rng);
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
  //  int saveThreads = omp_get_num_threads();
  //  omp_set_num_threads(numRFThreads);

  if(rfNumToRemovePerIteration) {
    cout << "\t\t\t\tRunning Iterative ReliefF..." << endl;
    reliefF->ComputeAttributeScoresIteratively();
  } else {
    if(analysisType == SNP_ONLY_ANALYSIS) {
      cout << "\t\t\t\tRunning standard ReliefF..." << endl;
      reliefF->ComputeAttributeScores();
    } else {
      cout << "\t\t\t\tRunning CLEAN SNPS ReliefF..." << endl;
      reliefF->ComputeAttributeScoresCleanSnps();
    }
  }
  rfScores = reliefF->GetScores();

  cout << "\t\t\t\tNormalizing ReliefF scores to 0-1..." << endl;
  pair<double, string> firstScore = rfScores[0];
  double minRFScore = firstScore.first;
  double maxRFScore = firstScore.first;
  EcScoresCIt rfScoresIt = rfScores.begin();
  for(; rfScoresIt != rfScores.end(); ++rfScoresIt) {
    pair<double, string> thisScore = *rfScoresIt;
    if(thisScore.first < minRFScore) {
      minRFScore = thisScore.first;
    }
    if(thisScore.first > maxRFScore) {
      maxRFScore = thisScore.first;
    }
  }

  // normalize attribute scores
  if(minRFScore == maxRFScore) {
    cerr << "\t\t\t\t\tWARNING: Relief-F min and max scores are the same."
            << "No normalization necessary." << endl;
    return true;
  }

  EcScores newRFScores;
  double rfRange = maxRFScore - minRFScore;
  for(EcScoresIt it = rfScores.begin(); it != rfScores.end(); ++it) {
    pair<double, string> thisScore = *it;
    double key = thisScore.first;
    string val = thisScore.second;
    newRFScores.push_back(make_pair((key - minRFScore) / rfRange, val));
  }

  rfScores.clear();
  rfScores = newRFScores;

  // omp_set_num_threads(saveThreads);

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
  sort(rjScores.begin(), rjScores.end(), scoresSortAscByName);
  sort(rfScores.begin(), rfScores.end(), scoresSortAscByName);
  freeEnergyScores.clear();
  EcScoresCIt rjIt = rjScores.begin();
  EcScoresCIt rfIt = rfScores.begin();
  for(; rjIt != rjScores.end(); ++rjIt, ++rfIt) {
    string val = rjIt->second;
    double key = rjIt->first;
    freeEnergyScores.push_back(make_pair((*rfIt).first + (temperature * key), val));
  }
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
bool EvaporativeCooling::RemoveWorstAttributes(unsigned int numToRemove) {
  unsigned int numToRemoveAdj = numToRemove;
  unsigned int numAttr = dataset->NumAttributes();
  if((numAttr - numToRemove) < numTargetAttributes) {
    cerr << "WARNING: attempt to remove " << numToRemove
            << " attributes which will remove more than target "
            << "number of attributes " << numTargetAttributes
            << ". Adjusting." << endl;
    numToRemoveAdj = numAttr - numTargetAttributes;
  }
  cout << "\t\t\tRemoving " << numToRemoveAdj << " attributes..." << endl;
  sort(freeEnergyScores.begin(), freeEnergyScores.end(), scoresSortAsc);
  for(unsigned int i = 0; i < numToRemoveAdj; ++i) {

    // worst score and attribute name
    pair<double, string> worst = freeEnergyScores[i];
    //    cout << "\t\t\t\tRemoving: " << worst.second
    //            << " (" << worst.first << ")" << endl;

    // save worst
    evaporatedAttributes.push_back(worst);
    // remove the attribute from those under consideration
    dataset->MaskRemoveAttribute(worst.second);
  }

  return true;
}
