/*
 * EvaporativeCoolingCLI.C - Bill White - 8/29/11
 *
 * Evaporative Cooling (EC) Optimization of Information Free Energy -
 * Command Line Interface (CLI)
 *
 * Implements the Evaporative Cooling algorithm as described in:
 * McKinney, et. al. "Capturing the Spectrum of Interaction Effects in Genetic
 * Association Studies by Simulated Evaporative Cooling Network Analysis."
 * PLoS Genetics, Vol 5, Issue 3, 2009.
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <map>
#include <sstream>
#include <vector>
#include <time.h>

#include <boost/program_options.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/parsers.hpp>

#include "Insilico.h"
#include "EvaporativeCooling.h"
#include "Dataset.h"

using namespace std;
using namespace insilico;
namespace po = boost::program_options;

int main(int argc, char** argv) {

  // ---------------------------------------------------------------------------
  cout << argv[0] << " starting..." << endl;
  clock_t t;
  t = clock();

  // ---------------------------------------------------------------------------
  cout << "\tProcessing command line arguments..." << endl;

  // command line processing variables: defaults and storage for boost
  bool verbose = false;
  // data set files
  string snpsFilename = "";
  string snpExclusionFile = "";
  bool doRecodeA = false;
  string cleanSnpsFilename = "";
  string numericsFilename = "";
  bool continuousPhenotype = false;
  string altPhenotypeFilename = "";
  string outputDatasetFilename = "";
  string outputFilesPrefix = "ec_run";
 // Random Jungle
  uli_t rjNumTrees = 100;
  unsigned int rjNumThreads = 0;
  // ReliefF
  unsigned int k = 10;
  unsigned int m = 0;
  string snpMetric = "gm";
  string numMetric = "manhattan";
  string weightByDistanceMethod = "equal";
  double weightByDistanceSigma = 2.0;
  unsigned int iterNumToRemove = 0;
  unsigned int iterPercentToRemove = 0;
  unsigned int rfNumThreads = 0;
  // diagnostic
  string diagnosticLogFilename = "";
  string diagnosticLevelsCountsFilename = "";
  // EC parameters
  string ecAlgorithmSteps = "all";
  unsigned int ecNumTarget = 0;
  unsigned int ecIterNumToRemove = 1;
  unsigned int ecIterPercentToRemove = 0;
  
  // declare the supported options
  po::options_description desc("Allowed options");
  desc.add_options()
          ("help", "produce help message")
          (
           "snp-data",
           po::value<string>(&snpsFilename),
           "read SNP attributes from genotype filename: txt, ARFF, plink (map/ped, binary, raw)"
           )
          (
           "snp-exclusion-file",
           po::value<string>(&snpExclusionFile),
           "file of SNP names to be excluded"
           )
          (
           "snp-data-clean",
           po::value<string>(&cleanSnpsFilename),
           "read SNP attributes from genotype filename - assumes no missing data, recodeA encoding"
           )
          (
           "recode-a",
           po::value<bool>(&doRecodeA)->default_value(doRecodeA),
           "do a plink recodeA encoding to insure genotype data values (0=no|1=yes)"
           )
          (
           "numeric-data",
           po::value<string>(&numericsFilename),
           "read SNP attributes from genotype filename: txt, ARFF, plink (map/ped, binary, raw)"
           )
          (
           "continuous-phenotype",
           po::value<bool>(&continuousPhenotype)->default_value(continuousPhenotype),
           "phenotype is continuous? (0=no|1=yes)"
           )
          (
           "alternate-pheno-file",
           po::value<string>(&altPhenotypeFilename),
           "specifies an alternative phenotype/class label file; one value per line"
           )
          (
           "out-dataset-filename",
           po::value<string>(&outputDatasetFilename),
           "write a new tab-delimited data set with EC filtered attributes"
           )
          (
           "out-files-prefix",
           po::value<string>(&outputFilesPrefix)->default_value(outputFilesPrefix),
           "use prefix for all output files"
           )
          (
           "verbose",
           po::value<bool>(&verbose)->default_value(verbose),
           "verbose output to stdout (0=no|1=yes)"
           )
          (
           "rj-num-trees",
           po::value<uli_t>(&rjNumTrees)->default_value(rjNumTrees),
           "Random Jungle number of trees to grow"
           )
          (
           "rj-num-threads",
           po::value<unsigned int>(&rjNumThreads)->default_value(rjNumThreads),
           "number of threads to use in Random Jungle, 0=all"
           )
          (
           "k-nearest-neighbors",
           po::value<unsigned int>(&k)->default_value(k),
           "set k nearest neighbors"
           )
          (
           "number-random-samples",
           po::value<unsigned int>(&m)->default_value(m),
           "number of random samples (0=all|1 <= n <= number of samples)"
           )
          (
           "snp-metric",
           po::value<string>(&snpMetric)->default_value(snpMetric),
           "metric for determining the difference between SNPs (gm|am)"
           )
            (
           "numeric-metric",
           po::value<string>(&numMetric)->default_value(numMetric),
           "metric for determining the difference between numeric attributes (manhattan=|euclidean)"
           )
          (
           "weight-by-distance-method",
           po::value<string>(&weightByDistanceMethod)->default_value(weightByDistanceMethod),
           "weight-by-distance method (equal|one_over_k|exponential)"
           )
          (
           "weight-by-distance-sigma",
           po::value<double>(&weightByDistanceSigma)->default_value(weightByDistanceSigma),
           "weight by distance sigma"
           )
          (
           "iter-remove-n",
           po::value<unsigned int>(&iterNumToRemove),
           "iterative ReliefF number of attributes to remove per iteration"
           )
          (
           "iter-remove-percent",
           po::value<unsigned int>(&iterPercentToRemove),
           "iterative ReliefF precentage of attributes to remove per iteration"
           )
          (
           "rf-num-threads",
           po::value<unsigned int>(&rfNumThreads)->default_value(rfNumThreads),
           "number of threads to use in Relief-F, 0=all"
           )
          (
           "diagnostic-tests",
           po::value<string> (&diagnosticLogFilename),
           "performs diagnostic tests and sends output to filename without running Relief-F"
           )
          (
           "diagnostic-levels-file",
           po::value<string > (&diagnosticLevelsCountsFilename),
           "write diagnostic attribute level counts to filename"
           )
          (
           "ec-algorithm-steps",
           po::value<string>(&ecAlgorithmSteps)->default_value(ecAlgorithmSteps),
           "EC steps to run (all|rj|rf)"
           )
          (
           "ec-num-target",
           po::value<unsigned int>(&ecNumTarget)->default_value(ecNumTarget),
           "EC N_target - target number of attributes to keep"
           )
          (
           "ec-iter-remove-n",
           po::value<unsigned int>(&ecIterNumToRemove)->default_value(ecIterNumToRemove),
           "Evaporative Cooling number of attributes to remove per iteration"
           )
          (
           "ec-iter-remove-percent",
           po::value<unsigned int>(&ecIterPercentToRemove),
           "Evaporative Cooling precentage of attributes to remove per iteration"
           )
          ;

  // parse the command line into a map
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  // po::store(po::command_line_parser(argc, argv).options(desc).positional(pd).run(), vm);
  po::notify(vm);

  if(vm.count("help")) {
    cerr << desc << "\n";
    exit(1);
  }

  // -------------------------------------------------------------------------
  // determine the analysis type
  cout << "\tDetermining analysis type..." << endl;
  AnalysisType analysisType = NO_ANALYSIS;
  if(vm.count("diagnostic-tests") || vm.count("diagnostic-levels-file")) {
    cout << "\t\tDiagnostic test requested." << endl;
    analysisType = DIAGNOSTIC_ANALYSIS;
  } else {
    if(vm.count("snp-data") && vm.count("numeric-data")) {
      cout << "\t\tIntegrated analysis requested." << endl;
      analysisType = INTEGRATED_ANALYSIS;
    }
    if(vm.count("snp-data-clean") && vm.count("numeric-data")) {
      cout << "\t\tIntegrated analysis requested." << endl;
      analysisType = INTEGRATED_ANALYSIS;
    }
    if((vm.count("snp-data") || vm.count("snp-data-clean")) &&
       !vm.count("numeric-data")) {
      if(vm.count("snp-data")) {
        cout << "\t\tSNP-only analysis requested." << endl;
        analysisType = SNP_ONLY_ANALYSIS;
      }
      else {
        cout << "\t\tClean SNP analysis requested." << endl;
        analysisType = SNP_CLEAN_ANALYSIS;     
      }
    } else {
      if(!vm.count("snp-data") && vm.count("numeric-data")) {
        cout << "\t\tNumeric-only analysis requested." << endl;
        analysisType = NUMERIC_ONLY_ANALYSIS;
        // must have an alternate phenotype file for numeric only
        if(!vm.count("alternate-pheno-file")) {
          cerr << "An alternate phenotype file must be specified with the "
                  << "--alternate-pheno-file option for numeric only data."
                  << endl;
          exit(1);
        }
      } else {
        if(vm.count("snp-data") && vm.count("numeric-data")) {
          cout << "\t\tIntegrated analysis requested." << endl;
          analysisType = INTEGRATED_ANALYSIS;
        }
        else {
          cerr << "ERROR: Could not determine the analysis to do based on "
                  << "command line options: " << endl << desc << endl;
          exit(1);
        }
      }
    }
  }

  if(continuousPhenotype && altPhenotypeFilename == "") {
    cerr << "ERROR: Continuous phenotype option --continuous-phenotype "
            << "requires --alternate-pheno-file option." << endl;
    exit(1);
  }

  // -------------------------------------------------------------------------
  // added bcw 7/15/11 - number of numerics or phenotypes might not be the
  // same as the data set - only load those in the numerics/phenotype file
  // if covariate or alternate phenotype file is present, read it first to
  // get the keys for reading the data set instances
  cout << "\tChecking for covariates and/or alternate phenotype files." << endl;
  vector<string> numericsIds;
  vector<string> phenoIds;
  if(analysisType == SNP_ONLY_ANALYSIS ||
     analysisType == NUMERIC_ONLY_ANALYSIS ||
     analysisType == INTEGRATED_ANALYSIS) {
    if(numericsFilename != "") {
      cout << "\t\tLoading individual IDs from covar file: "
              << numericsFilename << endl;
      if(!LoadIndividualIds(numericsFilename, numericsIds, true)) {
        exit(1);
      }
      // copy (numericsIds.begin(), numericsIds.end(), ostream_iterator<string> (cout, "\n"));
    }
    if(altPhenotypeFilename != "") {
      cout << "\t\tLoading individual IDs from alternate phenotype file: "
              << altPhenotypeFilename << endl;
      if(!LoadIndividualIds(altPhenotypeFilename, phenoIds, false)) {
        exit(1);
      }
      // copy(phenoIds.begin(), phenoIds.end(), ostream_iterator<string> (cout, "\n"));
    }
  }
  else {
    cout << "\t\tCovariate and alternate phenotype files not used for "
            << "this analysis type." << endl;
  }

  // -------------------------------------------------------------------------
  // find IDs for loading from the dataset
  cout << "\tDetermining the IDs to be read from the dataset." << endl;
  vector<string> indIds;
  if(numericsFilename != "" && altPhenotypeFilename != "") {
    cout << "\t\tIDs come from the numeric and the alternate "
            << "phenotype files. Checking for intersection/matches." << endl;

    // find intersection of numeric and phenotype IDs
    unsigned int maxMatches = max(numericsIds.size(), phenoIds.size());
    unsigned int maxMismatches = numericsIds.size() + phenoIds.size();
    indIds.resize(maxMatches);
    vector<string>::iterator goodIdsIt;
    vector<string> skippedIds(maxMismatches);
    vector<string>::iterator badIdsIt;

    sort(numericsIds.begin(), numericsIds.end());
    sort(phenoIds.begin(), phenoIds.end());
    goodIdsIt = set_intersection(numericsIds.begin(), numericsIds.end(),
                                 phenoIds.begin(), phenoIds.end(),
                                 indIds.begin());
    badIdsIt = set_difference(numericsIds.begin(), numericsIds.end(),
                              phenoIds.begin(), phenoIds.end(),
                              skippedIds.begin());
    if(skippedIds.begin() != badIdsIt) {
      cerr << "\t\t\tWARNING: Covariates and phenotypes files do not contain "
              << "the same IDs. These IDs differ: ";
      vector<string>::const_iterator skippedIt = skippedIds.begin();
      for(; skippedIt != badIdsIt; ++skippedIt) {
        cerr << *skippedIt << " ";
      }
      cerr << endl;
      // chaged to warning - 8/15/11
      // exit(1);
    }
    indIds.resize(int(goodIdsIt - indIds.begin()));
  } else {
    if(numericsFilename != "") {
      cout << "\t\tIDs come from the numerics file." << endl;
      indIds.resize(numericsIds.size());
      copy(numericsIds.begin(), numericsIds.end(), indIds.begin());
      // copy(numericsIds.begin(), numericsIds.end(), ostream_iterator<string > (cout, "\n"));
    } else {
      if(altPhenotypeFilename != "") {
        cout << "\t\tIDs come from the alternate phenotype file." << endl;
        // PrintVector(phenoIds, "phenoIds");
        indIds.resize(phenoIds.size());
        copy(phenoIds.begin(), phenoIds.end(), indIds.begin());
      } else {
        cout << "\t\tIDs are not needed for this analysis." << endl;
      }
    }
  }
  cout << "\t\t" << indIds.size()
          << " individual IDs read from numeric and/or phenotype file(s)."
          << endl;

  // -------------------------------------------------------------------------
  // prepare data for running EC
  cout << "\tPreparing data set for EC analysis..." << endl;
  Dataset* ds = 0;
  bool datasetLoaded = false;
  switch(analysisType) {
    case SNP_ONLY_ANALYSIS:
      cout << "\tReading SNPs data set" << endl;
      ds = ChooseSnpsDatasetByExtension(snpsFilename);
      datasetLoaded = ds->LoadDataset(snpsFilename, doRecodeA, "",
                                      altPhenotypeFilename, indIds,
                                      continuousPhenotype);
      break;
    case SNP_CLEAN_ANALYSIS:
      cout << "\tReading CLEAN SNPs data set" << endl;
      ds = ChooseSnpsDatasetByExtension(cleanSnpsFilename, true);
      datasetLoaded = ds->LoadDataset(cleanSnpsFilename, false, "",
                                      "", indIds, continuousPhenotype);
      break;
    case NUMERIC_ONLY_ANALYSIS:
      cout << "\tReading numerics only data set" << endl;
      ds = new Dataset();
      datasetLoaded = ds->LoadDataset("", doRecodeA, numericsFilename,
                                      altPhenotypeFilename, indIds,
                                      continuousPhenotype);
      break;
    case INTEGRATED_ANALYSIS:
      cout << "\tReading datasets for integrated analysis" << endl;
      ds = ChooseSnpsDatasetByExtension(snpsFilename);
      datasetLoaded = ds->LoadDataset(snpsFilename, doRecodeA, numericsFilename,
                                      altPhenotypeFilename, indIds,
                                      continuousPhenotype);
      break;
    case DIAGNOSTIC_ANALYSIS:
      cout << "\tPerforming SNP diagnostics on the data set" << endl;
      if(snpsFilename == "" || cleanSnpsFilename == "") {
        cerr << "Cannot run diagnostics without a SNP file specified with "
                << "--snp-data or --snp-data-clean flag." << endl;
        exit(1);
      }
      if(cleanSnpsFilename == "") {
        ds = ChooseSnpsDatasetByExtension(snpsFilename);
      }
      else {
        ds = ChooseSnpsDatasetByExtension(snpsFilename, true);
      }
      ds->LoadDataset(snpsFilename, doRecodeA, numericsFilename,
                      altPhenotypeFilename, indIds, continuousPhenotype);
      ds->RunSnpDiagnosticTests(diagnosticLevelsCountsFilename);
      if(diagnosticLevelsCountsFilename != "") {
        ds->WriteLevelCounts(diagnosticLevelsCountsFilename + ".counts");
      }
      // brutal exit!
      cout << argv[0] << " done." << endl;
      exit(0);
      break;
    case NO_ANALYSIS:
      cerr << "Analysis type could not be determined." << endl;
      exit(1);
    default:
      cerr << "Undefined analysis type: " << analysisType << endl;
      exit(1);
  }

  if(!datasetLoaded) {
    cerr << "ERROR: Failure to load dataset for analysis." << endl << endl;
    exit(1);
  }

  // happy lights
  if((analysisType == SNP_ONLY_ANALYSIS) ||
     (analysisType == SNP_CLEAN_ANALYSIS)) {
    ds->PrintStats();  
  }
  else {
    if(analysisType == NUMERIC_ONLY_ANALYSIS ||
       analysisType == INTEGRATED_ANALYSIS) {
      ds->PrintNumericsStats();
    }
  }

  // ---------------------------------------------------------------------------
  // FINALLY! run EC algorithm
  cout << "\tRunning EC..." << endl;
  EvaporativeCooling ec(ds, vm, analysisType);
  if(!ec.ComputeECScores()) {
    cerr << "ERROR: Failed to calculate EC scores." << endl;
    exit(1);
  }
  cout << "\tEC done." << endl;

  // ---------------------------------------------------------------------------
  // write the scores to the same name as the dataset with
  // <metric>.relieff suffix
  cout << "\tWriting EC scores to [" + outputFilesPrefix << ".ec]" << endl;
  ec.WriteAttributeScores(outputFilesPrefix);

  // write the EC filtered attributes as a new data set
  if(outputDatasetFilename != "") {
    cout << "\tWriting EC filtered data set to [" << outputDatasetFilename
            << "]" << endl;
    ds->WriteNewDataset(outputDatasetFilename);
  }

  // ---------------------------------------------------------------------------
  cout << "\tClean up and shutdown." << endl;
  // delete ds;

  float elapsedTime = (float)(clock() - t) / CLOCKS_PER_SEC;
  cout << "\tEC elapsed time " << elapsedTime << " secs." << endl;

  cout << argv[0] << " done." << endl;

  return 0;
}

