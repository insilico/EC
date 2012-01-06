/*
 * ReliefF.C - Bill White - 7/16/05
 * 
 * ReliefF algorithm implementation.
 *
 * Totally redone for the McKinney insilico lab in 2011.
 * Using OpenMP for multi-core parallelization - April 2011.
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <cmath>
#include <sstream>

#include <omp.h>
#include <boost/program_options.hpp>

#include "ReliefF.h"
#include "Dataset.h"
#include "DatasetInstance.h"
#include "StringUtils.h"
#include "DistanceMetrics.h"
#include "best_n.h"
#include "Insilico.h"

namespace po = boost::program_options;
using namespace std;
using namespace insilico;

/// scores map: score->attribute index
typedef vector<pair<double, unsigned int> > ScoresMap;
/// scores map iterator
typedef vector<pair<double, unsigned int> >::iterator ScoresMapIt;
/// attribute index map: attribute index->score
typedef vector<pair<unsigned int, double> > AttributeIndex;
/// attribute index map iterator
typedef vector<pair<unsigned int, double> >::const_iterator AttributeIndexIt;

/// attribute score sorting functor
bool scoreSort(const pair<double, string>& p1,
               const pair<double, string>& p2) {
  return p1.first < p2.first;
}

/// attribute index sorting functor
bool attributeSort(const pair<unsigned int, double>& p1,
                   const pair<unsigned int, double>& p2) {
  return p1.first < p2.first;
}

/// functor for T comparison
typedef pair<unsigned int, DatasetInstance*> T;

class deref_less : public std::binary_function<T, T, bool>
{
public:

  bool operator()(const T a, const T b) const {
    return(a.first < b.first);
  }
};

ReliefF::ReliefF(Dataset* ds, AnalysisType anaType) {
  cout << Timestamp() << "ReliefF initialization" << endl;
  if(ds) {
    dataset = ds;
  } else {
    cerr << "ERROR: dataset is not initialized" << endl;
    exit(-1);
  }
  analysisType = anaType;
  m = dataset->NumInstances();
  k = 10;
  weightByDistanceMethod = "equal";
  snpMetric = "gm";
  numMetric = "manhattan";
  removePerIteration = 0;

  cout << Timestamp() << "Number of samples: m = " << m << endl;
  randomlySelect = true;
  if(m == 0 || m == ds->NumInstances()) {
    // sample deterministically unless a sample size has been set
    cout << Timestamp() << "Sampling all instances deterministically" << endl;
    randomlySelect = false;
    m = ds->NumInstances();
  } else {
    cout << Timestamp() << "Sampling instances randomly" << endl;
    randomlySelect = true;
  }

  cout << Timestamp() << "Number of nearest neighbors: k = "
          << k << endl;

  one_over_m_times_k =
          1.0 / ((double) m * (double) k);
  if(to_upper(snpMetric) == "GM") {
    snpDiff = diffGMM;
  } else {
    if(to_upper(snpMetric) == "AM") {
      snpDiff = diffAMM;
    } else {
      cerr << "ERROR: [" << snpMetric << "] is not a valid SNP metric type" << endl;
      exit(1);
    }
  }

  if(to_upper(numMetric) == "MANHATTAN") {
    numDiff = diffManhattan;
  } else {
    cerr << "ERROR: [" << numMetric << "] is not a valid numeric metric type" << endl;
    exit(1);
  }

  cout << Timestamp() << "SNP distance metric: " << snpMetric << endl;
  cout << Timestamp() << "Continuous distance metric: " << numMetric << endl;

  int numProcs = omp_get_num_procs();
  int numThreads = omp_get_num_threads();
  cout << Timestamp() << numProcs << " OpenMP processors available" << endl;
  cout << Timestamp() << numThreads << " OpenMP threads running" << endl;

  vector<string> atrNames = dataset->GetAttributeNames();
  vector<string> numNames = dataset->GetNumericsNames();
  scoreNames.resize(atrNames.size() + numNames.size());
  copy(atrNames.begin(), atrNames.end(), scoreNames.begin());
  copy(numNames.begin(), numNames.end(), scoreNames.begin() + atrNames.size());
}

ReliefF::ReliefF(Dataset* ds, po::variables_map& vm, AnalysisType anaType) {
  cout << Timestamp() << "ReliefF initialization:" << endl;
  if(ds) {
    dataset = ds;
  } else {
    cerr << "ERROR: dataset is not initialized" << endl;
    exit(-1);
  }
  analysisType = anaType;

  // remove any attributes from distance calculations that
  // are in the attribute exclusion file
  if(vm.count("snp-exclusion-file")) {
    string snpExclusionFile = vm["snp-exclusion-file"].as<string > ();
    cout << Timestamp() << "Processing exclusion file [" << snpExclusionFile << "]" << endl;
    if(!ProcessExclusionFile(snpExclusionFile)) {
      cerr << "ERROR: Could not process exclusion file" << endl;
      exit(1);
    }
    cout << Timestamp() << dataset->NumVariables()
            << " attributes after processing exclusions file" << endl;
  }
  if(vm.count("number-random-samples")) {
    m = vm["number-random-samples"].as<unsigned int>();
    if(!m) {
      m = dataset->NumInstances();
    }
  } else {
    m = dataset->NumInstances();
  }
  if(vm.count("k-nearest-neighbors")) {
    k = vm["k-nearest-neighbors"].as<unsigned int>();
  } else {
    k = 10;
  }
  if(vm.count("snp-metric")) {
    snpMetric = vm["snp-metric"].as<string > ();
  } else {
    snpMetric = "gm";
  }
  if(vm.count("numeric-metric")) {
    numMetric = vm["numeric-metric"].as<string > ();
  } else {
    numMetric = "manhattan";
  }
  removePerIteration = 0;
  if(vm.count("iter-remove-n")) {
    removePerIteration = vm["iter-remove-n"].as<unsigned int>();
    if((removePerIteration < 1) ||
       (removePerIteration >= dataset->NumAttributes())) {
      cerr << "ERROR: Number to remove per iteratopn ["
              << removePerIteration << "] not in valid range" << endl;
      exit(-1);
    }
    cout << Timestamp() << "Iteratively removing " << removePerIteration << endl;
  }
  if(vm.count("iter-remove-percent")) {
    doRemovePercent = true;
    removePercentage = vm["iter-remove-percent"].as<unsigned int>() / 100.0;
    removePerIteration = (unsigned int)
            ((double) dataset->NumAttributes() * removePercentage + 0.5);
    if((removePerIteration < 1) ||
       (removePerIteration >= dataset->NumAttributes())) {
      cerr << "ERROR: Number to remove per iteratopn ["
              << removePerIteration << "] not in valid range" << endl;
      exit(-1);
    }
    cout << Timestamp() << "Iteratively removing " << (removePercentage * 100)
            << "% = " << removePerIteration << endl;
  }

  cout << Timestamp() << "Number of samples: m = " << m << endl;
  randomlySelect = true;
  if(m == 0 || m == ds->NumInstances()) {
    // sample deterministically unless a sample size has been set
    cout << Timestamp() << "Sampling all instances deterministically" << endl;
    randomlySelect = false;
    m = ds->NumInstances();
  } else {
    cout << Timestamp() << "Sampling instances randomly" << endl;
    randomlySelect = true;
  }

  cout << Timestamp() << "Number of nearest neighbors: k = " << k << endl;

  // k nearest neighbors and m randomly selected instances
  // spread differences and thus weight updates
  // over (m x k) iterations
  one_over_m_times_k =
          1.0 / ((double) m * (double) k);
  //                                  m           *                  k

  if(to_upper(snpMetric) == "GM") {
    snpDiff = diffGMM;
  } else {
    if(to_upper(snpMetric) == "AM") {
      snpDiff = diffAMM;
    } else {
      cerr << "ERROR: [" << snpMetric << "] is not a valid SNP metric type" << endl;
      exit(1);
    }
  }

  if(to_upper(numMetric) == "MANHATTAN") {
    numDiff = diffManhattan;
  } else {
    cerr << "ERROR: [" << numMetric << "] is not a valid numeric metric type" << endl;
    exit(1);
  }

  cout << Timestamp() << "SNP distance metric: " << snpMetric << endl;
  cout << Timestamp() << "Continuous distance metric: " << numMetric << endl;

  weightByDistanceMethod = vm["weight-by-distance-method"].as<string > ();
  if((weightByDistanceMethod != "exponential") &&
     (weightByDistanceMethod != "equal")) {
    cerr << "ERROR: Invalid --weight-by-distance-method: "
            << weightByDistanceMethod << endl;
    exit(1);
  }
  weightByDistanceSigma = vm["weight-by-distance-sigma"].as<double>();
  cout << Timestamp() << "Weight by distance method: " << weightByDistanceMethod;
  if(weightByDistanceMethod == "exponential") {
    cout << Timestamp() << ", using sigma = " << weightByDistanceSigma << endl;
  } else {
    cout << endl;
  }

  int numProcs = omp_get_num_procs();
  int numThreads = omp_get_num_threads();
  cout << Timestamp() << numProcs << " OpenMP processors available" << endl;
  cout << Timestamp() << numThreads << " OpenMP threads in work team" << endl;

  vector<string> atrNames = dataset->GetAttributeNames();
  vector<string> numNames = dataset->GetNumericsNames();
  scoreNames.resize(atrNames.size() + numNames.size());
  copy(atrNames.begin(), atrNames.end(), scoreNames.begin());
  copy(numNames.begin(), numNames.end(), scoreNames.begin() + atrNames.size());
}

ReliefF::~ReliefF() {
}

bool ReliefF::ComputeAttributeScores() {

  // changed from matric to map for ID matching - November 2011
  PreComputeDistancesByMap();

  // algorithm line 1
  W.resize(dataset->NumVariables(), 0.0);

  // pointer to the instance being sampled
  DatasetInstance* R_i;
  int i = 0;
  cout << Timestamp() << "Running Relief-F algorithm" << endl;
  cout << Timestamp() << "Averaging factor 1/(m*k): " << one_over_m_times_k << endl;

  vector<string> instanceIds = dataset->GetInstanceIds();
  cout << Timestamp();
  // algorithm line 2
  for(i = 0; i < (int) m; i++) {
    // algorithm line 3
    if(randomlySelect) {
      // randomly sample an instance (without replacement?)
      R_i = dataset->GetRandomInstance();
    } else {
      // deterministic/indexed instance sampling, ie, every instance against
      // every other instance
      unsigned int instanceIndex;
      dataset->GetInstanceIndexForID(instanceIds[i], instanceIndex);
      R_i = dataset->GetInstance(instanceIndex);
    }
    if(!R_i) {
      cerr << "ERROR: Random or indexed instance count not be found for index: ["
              << i << "]" << endl;
      return false;
    }
    ClassLevel class_R_i = R_i->GetClass();

    // algorithm lines 4, 5 and 6
    // find k nearest hits and nearest misses
    vector<unsigned int> hits;
    map<ClassLevel, vector<unsigned int> > misses;
    bool canGetNeighbors = false;
    canGetNeighbors = R_i->GetNNearestInstances(k, hits, misses);
    //    cout << "Instance class: " << R_i->GetClass() << ", hits: ";
    //    for(unsigned int ii = 0; ii < hits.size(); ++ii) {
    //      cout << hits[ii] << " ";
    //    }
    //    cout << endl << "Misses:" << endl;
    //    map<ClassLevel, vector<unsigned int> >::const_iterator iit;
    //    for(iit = misses.begin(); iit != misses.end(); ++iit) {
    //      cout << "Class: " << iit->first << ", misses: ";
    //      vector<unsigned int> ids = iit->second;
    //      for(unsigned int jj = 0; jj < ids.size(); ++jj) {
    //        cout << ids[jj] << " ";
    //      }
    //      cout << endl;
    //    }

    if(!canGetNeighbors) {
      cerr << "ERROR: relieff cannot get " << k
              << " nearest neighbors" << endl;
      return false;
    }

    // check algorithm preconditions
    if(hits.size() < 1) {
      cerr << "ERROR: No nearest hits found" << endl;
      return false;
    }
    if(hits.size() < k) {
      cerr << "ERROR: Could not find enough neighbors that are hits" << endl;
      exit(1);
    }
    map<ClassLevel, vector<unsigned int> >::const_iterator it;
    for(it = misses.begin(); it != misses.end(); ++it) {
      vector<unsigned int> missIds = it->second;
      if(missIds.size() < 1) {
        cerr << "ERROR: No nearest misses found" << endl;
        return false;
      }
      if(missIds.size() < k) {
        cerr << "ERROR: Could not find enough neighbors that are misses" << endl;
        return false;
      }
      if(missIds.size() != hits.size()) {
        cerr << "ERROR: Could not find equal number of neighbors for hits and misses:"
                << hits.size() << " vs. " << misses.size() << endl;
        return false;
      }
    }

    // UPDATE WEIGHTS FOR ATTRIBUTE 'A' BASED ON THIS AND NEIGHBORING INSTANCES
    // update weights/relevance scores for each attribute averaged
    // across k nearest neighbors and m (possibly randomly) selected instances
    unsigned int A = 0;
    unsigned int scoresIdx = 0;
    if(dataset->HasGenotypes()) {
      vector<unsigned int> attributeIndicies = 
        dataset->MaskGetAttributeIndices(DISCRETE_TYPE);
      // algorithm line 7
      for(unsigned int attrIdx = 0; attrIdx < attributeIndicies.size(); ++attrIdx) {
        A = attributeIndicies[attrIdx];
        double hitSum = 0.0, missSum = 0.0;
        // algorithm line 8
        for(unsigned int j = 0; j < k; j++) {
          DatasetInstance* H_j = dataset->GetInstance(hits[j]);
          hitSum += (snpDiff(A, R_i, H_j) * one_over_m_times_k);
        }
        // algorithm line 9
        map<ClassLevel, vector<unsigned int> >::const_iterator mit;
        for(mit = misses.begin(); mit != misses.end(); ++mit) {
          ClassLevel C = mit->first;
          vector<unsigned int> missIds = mit->second;
          double P_C = dataset->GetClassProbability(C);
          double P_C_R = dataset->GetClassProbability(class_R_i);
          double adjustmentFactor = P_C / (1.0 - P_C_R);
          double tempSum = 0.0;
          for(unsigned int j = 0; j < k; j++) {
            DatasetInstance* M_j = dataset->GetInstance(missIds[j]);
            tempSum += (snpDiff(A, R_i, M_j) * one_over_m_times_k);
          } // nearest neighbors
          missSum += (adjustmentFactor * tempSum);
        }
        W[scoresIdx] = W[scoresIdx] - hitSum + missSum;
        ++scoresIdx;

      } // all attributes
    } // has genotypes

    // loop here for numeric attributes if they exist - 6/19/11
    if(dataset->HasNumerics()) {
      vector<unsigned int> numericIndices = 
        dataset->MaskGetAttributeIndices(NUMERIC_TYPE);
      for(unsigned int numIdx = 0; numIdx < numericIndices.size(); ++numIdx) {
        A = numericIndices[numIdx];
        double hitSum = 0.0, missSum = 0.0;
        for(unsigned int j = 0; j < k; j++) {
          DatasetInstance* H_j = dataset->GetInstance(hits[j]);
          hitSum += (numDiff(A, R_i, H_j) * one_over_m_times_k);
        }

        map<ClassLevel, vector<unsigned int> >::const_iterator mit;
        for(mit = misses.begin(); mit != misses.end(); ++mit) {
          ClassLevel C = mit->first;
          vector<unsigned int> missIds = mit->second;
          double P_C = dataset->GetClassProbability(C);
          double P_C_R = dataset->GetClassProbability(class_R_i);
          double adjustmentFactor = P_C / (1.0 - P_C_R);
          double tempSum = 0.0;
          for(unsigned int j = 0; j < k; j++) {
            DatasetInstance* M_j = dataset->GetInstance(missIds[j]);
            tempSum += (numDiff(A, R_i, M_j) * one_over_m_times_k);
          } // nearest neighbors
          missSum += (adjustmentFactor * tempSum);
        }
        W[scoresIdx] = W[scoresIdx] - hitSum + missSum;
        ++scoresIdx;

      }
    } // has numerics

    // happy lights
    if(i && ((i % 100) == 0)) {
      cout << i << "/" << m << " ";
      cout.flush();
    }
    if(i && ((i % 1000) == 0)) {
      cout << endl << Timestamp();
    }

  } // number to randomly select

  cout << i << "/" << m << endl;

  return true;
}

bool ReliefF::ComputeAttributeScoresCleanSnps() {

  PreComputeDistances();

  cout << Timestamp() << "Running Relief-F algorithm" << endl;
  cout.flush();

  W.resize(dataset->NumVariables(), 0.0);

  // pointer to the instance being sampled
  DatasetInstance* R;
  int i = 0;

  // K NEAREST NEIGHBORS
  // find k nearest hits and nearest misses
  vector<unsigned int> hits;
  vector<unsigned int> misses;

  //#pragma omp parallel default(shared)
  // {
  //#pragma omp parallel for
  vector<string> instanceIds = dataset->GetInstanceIds();
  cout << Timestamp();
  for(i = 0; i < (int) m; i++) {
    // INSTANCE SAMPLING
    // sample an instance from all instances
    if(randomlySelect) {
      // randomly sample an instance (without replacement?)
      R = dataset->GetRandomInstance();
    } else {
      // deterministic/indexed instance sampling, ie, every instance against
      // every other instance
      unsigned int instanceIndex;
      dataset->GetInstanceIndexForID(instanceIds[i], instanceIndex);
      R = dataset->GetInstance(instanceIndex);
    }
    if(!R) {
      cerr << "ERROR: Random or indexed instance count not be found for index: ["
              << i << "]" << endl;
      exit(1);
    }

    //      cout << "Getting nearest hits and misses for instance index: "
    //              << i << ", which points to: " << R << endl;
    bool canGetNeighbors = R->GetNNearestInstances(k, hits, misses);
    if(!canGetNeighbors) {
      cerr << "\nERROR: relieff cannot get " << k
              << " nearest neighbors" << endl;
      exit(1);
    }

    // check algorithm preconditions
    if(hits.size() < 1) {
      cerr << "\nERROR: No nearest hits found" << endl;
      exit(1);
    }
    if(hits.size() < k) {
      cerr << "\nERROR: Could not find enough neighbors that are hits" << endl;
      exit(1);
    }
    if(misses.size() < 1) {
      cerr << "\nERROR: No nearest misses found" << endl;
      exit(1);
    }
    if(misses.size() < k) {
      cerr << "\nERROR: Could not find enough neighbors that are misses" << endl;
      exit(1);
    }
    if(misses.size() != hits.size()) {
      cerr << "\nERROR: Could not find equal number of neighbors for hits and misses:"
              << hits.size() << " vs. " << misses.size() << endl;
      exit(1);
    }

    // UPDATE WEIGHTS FOR ATTRIBUTE 'A' BASED ON THIS AND NEIGHBORING INSTANCES
    // update weights/relevance scores for each attribute averaged
    // across k nearest neighbors and m (possibly randomly) selected instances
    // dataset->PrintMaskStats();
    unsigned int A = 0;
    vector<unsigned int> attributeIndicies = 
      dataset->MaskGetAttributeIndices(DISCRETE_TYPE);
    for(unsigned int attrIdx = 0; attrIdx < attributeIndicies.size(); ++attrIdx) {
      A = attributeIndicies[attrIdx];
      //        cout << attrIdx << " (" << A << ") " << endl;
      //        cout << "Hits:" << endl;
      //        copy(hits.begin(), hits.end(), ostream_iterator<unsigned int>(cout, "\n"));
      //        cout << "Misses STL iterators:" << endl;
      //        copy(misses.begin(), misses.end(), ostream_iterator<unsigned int>(cout, "\n"));
      //        cout << "Misses for loop with [] operator:" << endl;
      //        for(unsigned int hi=0; hi < hits.size(); ++hi) {
      //          cout << hi << " => " << hits[hi] << endl;
      //        }
      double hitSum = 0.0, missSum = 0.0;
      for(unsigned int j = 0; j < k; j++) {
        //          cout << j << ", " << hits[j] << endl;
        DatasetInstance* H = dataset->GetInstance(hits[j]);
        hitSum += snpDiff(A, R, H);
        DatasetInstance* M = dataset->GetInstance(misses[j]);
        missSum += snpDiff(A, R, M);
        // "normalize" and update weights/relevance scores
        W[attrIdx] += ((missSum - hitSum) / one_over_m_times_k);
      }
    }
    //      cout << endl;

    // happy lights
    if(i && ((i % 100) == 0)) {
      cout << i << "/" << m << " ";
      cout.flush();
    }
    if(i && ((i % 1000) == 0)) {
      cout << endl << Timestamp();
    }

  }
  //  } // end omp parallel
  cout << m << "/" << m << endl;

  return true;
}

bool ReliefF::ComputeAttributeScoresIteratively() {

  // save the current dataset mask
  dataset->MaskPushAll();

  // IterativeReliefF or TuRF (Tuned Relief-F)
  unsigned int iterations = 1;
  while(dataset->NumVariables() > 0) {

    cout << Timestamp() << "------------------------------------------------------------"
            << "-----------------------------------------" << endl;
    cout << Timestamp() << "[" << iterations << "] Working attributes: "
            << dataset->NumVariables() << endl;

    ComputeAttributeScores();
    vector<pair<double, string> > attributeScores = GetScores();

    // save worst attributes and remove from consideration on next iteration
    sort(attributeScores.begin(), attributeScores.end(), scoreSort);
    unsigned int removeThisIteration = 0;
    if(dataset->NumVariables() < removePerIteration) {
      removeThisIteration = dataset->NumVariables();
    } else {
      if(doRemovePercent) {
        removeThisIteration = (unsigned int)
                ((double) dataset->NumAttributes() * removePercentage + 0.5);
      } else {
        removeThisIteration = removePerIteration;
      }
    }
    for(unsigned int i = 0; i < removeThisIteration; ++i) {
      string attributeToDelete = attributeScores[i].second;
      //      cout << "\t\t\tremoving attribute: " << attributeToDelete << endl;
      if(dataset->MaskSearchAttribute(attributeToDelete, DISCRETE_TYPE)) {
        dataset->MaskRemoveAttribute(attributeToDelete, DISCRETE_TYPE);
        //        scores[dataset->GetAttributeIndexFromName(attributeToDelete)] =
        //                attributeScores[i].first;
      } else {
        if(dataset->MaskSearchAttribute(attributeToDelete, NUMERIC_TYPE)) {
          dataset->MaskRemoveAttribute(attributeToDelete, NUMERIC_TYPE);
          //          scores[dataset->NumAttributes() +
          //                 dataset->GetNumericIndexFromName(attributeToDelete)] =
          //                attributeScores[i].first;
        } else {
          cerr << "ERROR: ReliefF::ComputeAttributeScoresIteratively: "
                  << "could not find attribute name in data set: "
                  << attributeToDelete << endl;
          return false;
        }
      }
      finalScores[attributeToDelete] = attributeScores[i].first;
    }

    ++iterations;

    ResetForNextIteration();
  } // iterate

  // populate finalScores with remaining scores
  vector<double>::const_iterator scoresIt;
  vector<string> attrNames = dataset->GetAttributeNames();
  //  cout << "Scores: " << scores.size()
  //          << ", Attribute names: " << attrNames.size()
  //          << ", Final Scores: " << finalScores.size() << endl;
  for(unsigned int i = 0; i < attrNames.size(); ++i) {
    cout << attrNames[i] << " => " << scoresIt[i] << endl;
    finalScores[attrNames[i]] = W[i];
  }

  W.resize(scoreNames.size());
  for(unsigned int i = 0; i < scoreNames.size(); ++i) {
    if(finalScores.find(scoreNames[i]) == finalScores.end()) {
      cerr << "ERROR: Logic error. See Bill" << endl;
      exit(1);
    }
    W[i] = finalScores[scoreNames[i]];
  }

  // restore the dataset attribute mask
  dataset->MaskPopAll();

  return true;
}

bool ReliefF::ResetForNextIteration() {

  PreComputeDistancesByMap();

  return true;
}

void ReliefF::PrintAttributeScores(ofstream & outFile) {
  vector<double>::const_iterator scoresIt = W.begin();
  unsigned int nameIdx = 0;
  for(; scoresIt != W.end(); ++scoresIt) {
    outFile << fixed << setprecision(6) << *scoresIt << "\t"
            << scoreNames[nameIdx] << endl;
    //    outFile << fixed << *scoresIt << "\t" << scoreNames[nameIdx] << endl;
    ++nameIdx;
  }
}

void ReliefF::WriteAttributeScores(string baseFilename) {

  string resultsFilename = baseFilename;
  if(dataset->HasContinuousPhenotypes()) {
    resultsFilename += ".rrelieff";
  } else {
    resultsFilename += ".relieff";
  }

  ofstream outFile;
  outFile.open(resultsFilename.c_str());
  if(outFile.bad()) {
    cerr << "ERROR: Could not open scores file " << resultsFilename
            << "for writing" << endl;
    exit(1);
  }
  PrintAttributeScores(outFile);
  outFile.close();
}

bool ReliefF::ProcessExclusionFile(string exclusionFilename) {
  ifstream dataStream(exclusionFilename.c_str());
  if(!dataStream.is_open()) {
    cerr << "ERROR: Could not open exclusion file: " << exclusionFilename << endl;
    return false;
  }

  // temporary string for reading file lines
  string line;
  unsigned int lineNumber = 0;
  while(getline(dataStream, line)) {
    ++lineNumber;
    string attributeName = trim(line);
    if(!dataset->MaskRemoveAttribute(attributeName, DISCRETE_TYPE)) {
      cerr << "ERROR: attribute to exclude [" << attributeName
              << "] on line [" << lineNumber
              << "] in the exclusion file. It is not in the data set" << endl;
      return false;
    }
  }
  dataStream.close();

  return true;
}

double
ReliefF::ComputeInstanceToInstanceDistance(DatasetInstance* dsi1,
                                           DatasetInstance * dsi2) {
  double distance = 0;

  if(dataset->HasGenotypes()) {
    vector<unsigned int> attributeIndices = 
      dataset->MaskGetAttributeIndices(DISCRETE_TYPE);
    for(unsigned int i = 0; i < attributeIndices.size(); ++i) {
      distance += snpDiff(attributeIndices[i], dsi1, dsi2);
    }
    // cout << "SNP distance = " << distance << endl;
  }

  // added 6/16/11
  // compute numeric distances
  if(dataset->HasNumerics()) {
    vector<unsigned int> numericIndices = 
      dataset->MaskGetAttributeIndices(NUMERIC_TYPE);
    for(unsigned int i = 0; i < numericIndices.size(); ++i) {
      double numDistance = numDiff(numericIndices[i], dsi1, dsi2);
      // cout << "Numeric distance " << i << " => " << numDistance << endl;
      distance += numDistance;
    }
  }

  return distance;
}

bool ReliefF::PreComputeDistances() {
  cout << Timestamp() << "Precomputing instance distances" << endl;
  map<string, unsigned int> instanceMask = dataset->MaskGetInstanceMask();
  vector<string> instanceIds = dataset->MaskGetInstanceIds();
  int numInstances = instanceIds.size();

  // create a distance matrix
  cout << Timestamp() << "Allocating distance matrix";
  double** distanceMatrix;
  distanceMatrix = new double*[numInstances];
  for(int i = 0; i < numInstances; ++i) {
    distanceMatrix[i] = new double[numInstances];
    for(int j = 0; j < numInstances; ++j) {
      distanceMatrix[i][j] = 0.0;
    }
  }
  cout << " done" << endl;

  // populate the matrix - upper triangular
  // NOTE: make complete symmetric matrix for neighbor-to-neighbor sums
  cout << Timestamp() << "1) Computing instance-to-instance distances... ";
  //  omp_set_nested(1);
#pragma omp parallel for schedule(dynamic, 1)
  for(int i = 0; i < numInstances; ++i) {
    // cout << "Computing instance to instance distances. Row: " << i << endl;
    // #pragma omp parallel for
    for(int j = i + 1; j < numInstances; ++j) {
      unsigned int dsi1Index;
      dataset->GetInstanceIndexForID(instanceIds[i], dsi1Index);
      unsigned int dsi2Index;
      dataset->GetInstanceIndexForID(instanceIds[j], dsi2Index);
      distanceMatrix[i][j] =
              ComputeInstanceToInstanceDistance(dataset->GetInstance(dsi1Index),
                                                dataset->GetInstance(dsi2Index));
      // cout << i << ", " << j << " => " << distanceMatrix[i][j] << endl;
      distanceMatrix[j][i] = distanceMatrix[i][j];
    }
    if(i % 100 == 0) {
      cout << i << "/" << numInstances << " ";
      cout.flush();
    }
    if(i && ((i % 1000) == 0)) {
      cout << endl << Timestamp();
    }

  }
  cout << numInstances << "/" << numInstances << " done" << endl;

  //  DEBUG
  //  ofstream outFile;
  //  outFile.open("distanceMatrix.csv");
  //  for(unsigned int i=0; i < dataset->NumInstances(); ++i) {
  //    for(unsigned int j=0; j < dataset->NumInstances(); ++j) {
  //      if(j)
  //        outFile << "," << distanceMatrix[i][j];
  //      else
  //        outFile << distanceMatrix[i][j];
  //    }
  //    outFile << endl;
  //  }
  //  outFile.close();
  //  DEBUG

  // for each instance: if discrete class, store the distance sums for same
  // and different classes, else store distances to all other instances
  // (regression ReliefF)
  if(dataset->HasContinuousPhenotypes()) {
    cout << Timestamp() << "2) Calculating continuous phenotype nearest neighbors... ";
  } else {
    cout << Timestamp() << "2) Calculating same and different class nearest neighbors... ";
  }

  pair<string, string> key;
  DistancePair nnInfo;
  for(int i = 0; i < (int) numInstances; ++i) {
    unsigned int thisInstanceIndex = instanceMask[instanceIds[i]];
    DatasetInstance* thisInstance = dataset->GetInstance(thisInstanceIndex);

    if(dataset->HasContinuousPhenotypes()) {
      DistancePairs instanceDistances;
      for(int j = 0; j < (int) numInstances; ++j) {
        if(i == j) continue;
        double instanceToInstanceDistance = distanceMatrix[i][j];
        DistancePair nearestNeighborInfo;
        nearestNeighborInfo = make_pair(instanceToInstanceDistance, instanceIds[j]);
        instanceDistances.push_back(nearestNeighborInfo);
      }
      dataset->GetInstance(i)->SetDistanceSums(k, instanceDistances);
    } else {
      ClassLevel thisClass = thisInstance->GetClass();
      DistancePairs sameSums;
      // changed to an array for multiclass - 12/1/11
      map<ClassLevel, DistancePairs> diffSums;
      for(int j = 0; j < numInstances; ++j) {
        if(i == j) continue;
        if(j < i) {
          key = make_pair(instanceIds[j], instanceIds[i]);
        } else {
          key = make_pair(instanceIds[i], instanceIds[j]);
        }
        double instanceToInstanceDistance = distanceMatrix[i][j];
        unsigned int otherInstanceIndex = instanceMask[instanceIds[j]];
        DatasetInstance* otherInstance =
                dataset->GetInstance(otherInstanceIndex);
        nnInfo = make_pair(instanceToInstanceDistance, instanceIds[j]);
        if(otherInstance->GetClass() == thisClass) {
          sameSums.push_back(nnInfo);
        } else {
          ClassLevel otherClass = otherInstance->GetClass();
          diffSums[otherClass].push_back(nnInfo);
        }
      }
      thisInstance->SetDistanceSums(k, sameSums, diffSums);
    }

    if(i % 100 == 0) {
      cout << i << "/" << numInstances << " ";
      cout.flush();
    }
    if(i && ((i % 1000) == 0)) {
      cout << endl << Timestamp();
    }
  }
  cout << Timestamp() << numInstances << "/" << numInstances << " done" << endl;

  cout << Timestamp() << "3) Calculating weight by distance factors for "
          << "nearest neighbors... " << endl;
  ComputeWeightByDistanceFactors();

  // release the dynamically-allocated distance matrix
  cout << Timestamp() << "Freeing distance matrix memory";
  for(int i = 0; i < numInstances; ++i) {
    delete [] distanceMatrix[i];
  }
  delete [] distanceMatrix;
  cout << Timestamp() << "done" << endl;

  return true;
}

bool ReliefF::PreComputeDistancesByMap() {

  cout << Timestamp() << "Precomputing instance distances by map" << endl;
  map<string, unsigned int> instanceMask = dataset->MaskGetInstanceMask();
  vector<string> instanceIds = dataset->MaskGetInstanceIds();
  int numInstances = instanceIds.size();

  cout << Timestamp()
          << "1) Computing instance-to-instance distances in parallel... ";
  map<pair<string, string>, double > distanceMatrix;
  //        ID1     ID2     dist
  int i = 0;
#pragma omp parallel for schedule(dynamic, 1)
  for(i = 0; i < numInstances; ++i) {
    for(int j = i + 1; j < numInstances; ++j) {
      string instanceId1 = instanceIds[i];
      string instanceId2 = instanceIds[j];
      unsigned int dsi1Index = instanceMask[instanceId1];
      unsigned int dsi2Index = instanceMask[instanceId2];
      pair<string, string> key = make_pair(instanceId1, instanceId2);
      pair<string, string> symkey = make_pair(instanceId2, instanceId1);
      double distance =
              ComputeInstanceToInstanceDistance(dataset->GetInstance(dsi1Index),
                                                dataset->GetInstance(dsi2Index));
#pragma omp critical
      {
        distanceMatrix[key] = distance;
        //        cout << i << "," << j << " -> (" << key.first << "," << key.second
        //             << ") => " << distanceMatrix[key] << endl;
      }
    }
    if(i && (i % 100 == 0)) {
      cout << i << "/" << numInstances << " ";
      cout.flush();
    }
    if(i && ((i % 1000) == 0)) {
      cout << endl << Timestamp();
    }
  }
  cout << numInstances << "/" << numInstances << " done" << endl;

  if(dataset->HasContinuousPhenotypes()) {
    cout << Timestamp() << "2) Calculating continuous phenotype nearest neighbors... ";
  } else {
    // multiclass - 12/1/11
    if(dataset->NumClasses() > 2) {
      cout << Timestamp() << "2) Calculating same and different classes nearest neighbors... ";
    } else {
      cout << Timestamp() << "2) Calculating same and different class nearest neighbors... ";
    }
  }

  pair<string, string> key;
  DistancePair nnInfo;
  double instanceToInstanceDistance;
  for(i = 0; i < numInstances; ++i) {
    unsigned int thisInstanceIndex = instanceMask[instanceIds[i]];
    DatasetInstance* thisInstance = dataset->GetInstance(thisInstanceIndex);
    if(dataset->HasContinuousPhenotypes()) {
      DistancePairs instanceDistances;
      for(int j = 0; j < numInstances; ++j) {
        if(i == j) continue;
        if(j < i) {
          key = make_pair(instanceIds[j], instanceIds[i]);
        } else {
          key = make_pair(instanceIds[i], instanceIds[j]);
        }
        instanceToInstanceDistance = distanceMatrix[key];
        nnInfo = make_pair(instanceToInstanceDistance, instanceIds[j]);
        instanceDistances.push_back(nnInfo);
      }
      thisInstance->SetDistanceSums(k, instanceDistances);
    } else {
      ClassLevel thisClass = thisInstance->GetClass();
      DistancePairs sameSums;
      // changed to an array for multiclass - 12/1/11
      map<ClassLevel, DistancePairs> diffSums;
      for(int j = 0; j < numInstances; ++j) {
        if(i == j) continue;
        if(j < i) {
          key = make_pair(instanceIds[j], instanceIds[i]);
        } else {
          key = make_pair(instanceIds[i], instanceIds[j]);
        }
        instanceToInstanceDistance = distanceMatrix[key];
        unsigned int otherInstanceIndex = instanceMask[instanceIds[j]];
        DatasetInstance* otherInstance =
                dataset->GetInstance(otherInstanceIndex);
        nnInfo = make_pair(instanceToInstanceDistance, instanceIds[j]);
        if(otherInstance->GetClass() == thisClass) {
          sameSums.push_back(nnInfo);
        } else {
          ClassLevel otherClass = otherInstance->GetClass();
          diffSums[otherClass].push_back(nnInfo);
        }
      }
      thisInstance->SetDistanceSums(k, sameSums, diffSums);
    }

    if(i % 100 == 0) {
      cout << i << "/" << numInstances << " ";
      cout.flush();
    }
    if(i && ((i % 1000) == 0)) {
      cout << endl << Timestamp();
    }
  }
  cout << numInstances << "/" << numInstances << " done" << endl;

  cout << Timestamp() << "3) Calculating weight by distance factors for "
          << "nearest neighbors... " << endl;
  ComputeWeightByDistanceFactors();

  return true;
}

vector<pair<double, string> > ReliefF::GetScores() {

  vector<pair<double, string> > returnScores;
  vector<double>::const_iterator scoresIt = W.begin();
  unsigned int nameIdx = 0;
  vector<string> maskNames = dataset->MaskGetAllVariableNames();
  for(; scoresIt != W.end(); ++scoresIt) {
    returnScores.push_back(make_pair(*scoresIt, maskNames[nameIdx]));
    ++nameIdx;
  }

  return returnScores;
}

bool ReliefF::ComputeWeightByDistanceFactors() {

  vector<string> instanceIds = dataset->GetInstanceIds();
  for(unsigned int i = 0; i < dataset->NumInstances(); ++i) {

    // this instance
    unsigned int instanceIndex;
    dataset->GetInstanceIndexForID(instanceIds[i], instanceIndex);
    DatasetInstance* dsi = dataset->GetInstance(instanceIndex);

    vector<double> d1_ij;
    double d1_ij_sum = 0.0;
    for(unsigned int rank_j = 1; rank_j <= k; ++rank_j) {
      double d1_ij_value = 0.0;
      if(weightByDistanceMethod == "exponential") {
        double exponentArg = (double) rank_j / weightByDistanceSigma;
        d1_ij_value = exp(-(exponentArg * exponentArg));
      } else {
        if(weightByDistanceMethod == "one_over_k") {
          d1_ij_value = 1.0 / (double) rank_j;
        } else {
          // equal
          d1_ij_value = 1.0 / (double) k;
        }
      }
      d1_ij.push_back(d1_ij_value);
      d1_ij_sum += d1_ij_value;
    }

    // "normalize" the factors - divide through by the total/sum
    dsi->ClearInfluenceFactors();
    for(unsigned int neighborIdx = 0; neighborIdx < k; ++neighborIdx) {
      double influenceFactorD = d1_ij[neighborIdx] / d1_ij_sum;
      //      cout << "d_ij: " << d1_ij[neighborIdx]
      //              << ", cummulative sum: " << d1_ij_sum
      //              << ", normalized value: " << influenceFactorD << endl;
      dsi->AddInfluenceFactorD(influenceFactorD);
    }
    //    cout << "---------------------------------------------------------" << endl;
  } // end all instances

  return true;
}

void librelieff_is_present(void) {
  ;
}