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
#include "FilesystemUtils.h"

using namespace std;
namespace po = boost::program_options;

int main(int argc, char** argv) {

  // ---------------------------------------------------------------------------
  cout << Timestamp() << argv[0] << " starting" << endl;
  clock_t t;
  t = clock();

  // ---------------------------------------------------------------------------
  cout << Timestamp() << "Processing command line arguments" << endl;

  // command line processing variables: defaults and storage for boost
  bool verbose = false;
  // data set files
  string snpsFilename = "";
  string snpExclusionFile = "";
  string cleanSnpsFilename = "";
  string numericsFilename = "";
  string altPhenotypeFilename = "";
  string outputDatasetFilename = "";
  string outputFilesPrefix = "ec_run";
 // Random Jungle
  uli_t rjNumTrees = 1000;
  // ReliefF
  unsigned int k = 10;
  unsigned int m = 0;
  string snpMetric = "gm";
  string numMetric = "manhattan";
  string weightByDistanceMethod = "equal";
  double weightByDistanceSigma = 2.0;
  // diagnostic
  string diagnosticLogFilename = "";
  string diagnosticLevelsCountsFilename = "";
  // EC parameters
  string ecAlgorithmSteps = "all";
  unsigned int ecNumTarget = 1;
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
           "numeric-data",
           po::value<string>(&numericsFilename),
           "read SNP attributes from genotype filename: txt, ARFF, plink (map/ped, binary, raw)"
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

  if(vm.count("help") || (argc == 0)) {
    cerr << desc << endl;
    exit(COMMAND_LINE_ERROR);
  }

  /// determine the output data set type
  OutputDatasetType outputDatasetType = NO_OUTPUT_DATASET;
  if(outputDatasetFilename != "") {
    // determine the dataset type
    string outFileExtension = GetFileExtension(outputDatasetFilename);
    if(outFileExtension == "txt") {
      outputDatasetType = TAB_DELIMITED_DATASET;
    } else {
      if(outFileExtension == "csv") {
        outputDatasetType = CSV_DELIMITED_DATASET;
      } else {
        if(outFileExtension == "arff") {
          outputDatasetType = ARFF_DATASET;
        } else {
          cerr << "Unrecognized output data set filename extension: ["
                  << outFileExtension << "]" << endl;
          exit(COMMAND_LINE_ERROR);
        }
      }
    }
  }
  cout << Timestamp() << "Writing ReliefF filtered data set to ["
          << outputDatasetFilename
          << "]" << endl;

  /// determine the analysis type
  cout << Timestamp() << "Determining analysis type" << endl;
  AnalysisType analysisType = NO_ANALYSIS;
  if(vm.count("diagnostic-tests") || vm.count("diagnostic-levels-file")) {
    cout << Timestamp() << "Diagnostic test requested" << endl;
    analysisType = DIAGNOSTIC_ANALYSIS;
  } else {
    if(vm.count("snp-data") && vm.count("numeric-data")) {
      cout << Timestamp() << "Integrated analysis requested" << endl;
      analysisType = INTEGRATED_ANALYSIS;
    }
    if(vm.count("snp-data-clean") && vm.count("numeric-data")) {
      cout << Timestamp() << "Integrated analysis requested" << endl;
      analysisType = INTEGRATED_ANALYSIS;
    }
    if((vm.count("snp-data") || vm.count("snp-data-clean")) &&
       !vm.count("numeric-data")) {
      if(vm.count("snp-data")) {
        cout << Timestamp() << "SNP-only analysis requested" << endl;
        analysisType = SNP_ONLY_ANALYSIS;
      }
      else {
        cout << Timestamp() << "Clean SNP analysis requested" << endl;
        analysisType = SNP_CLEAN_ANALYSIS;     
      }
    } else {
      if(!vm.count("snp-data") && vm.count("numeric-data")) {
        cout << Timestamp() << "Numeric-only analysis requested" << endl;
        analysisType = NUMERIC_ONLY_ANALYSIS;
        // must have an alternate phenotype file for numeric only
        if(!vm.count("alternate-pheno-file")) {
          cerr << "An alternate phenotype file must be specified with the "
                  << "--alternate-pheno-file option for numeric only data"
                  << endl;
          exit(COMMAND_LINE_ERROR);
        }
      } else {
        if(vm.count("snp-data") && vm.count("numeric-data")) {
          cout << "\t\tIntegrated analysis requested" << endl;
          analysisType = INTEGRATED_ANALYSIS;
        }
        else {
          cerr << "ERROR: Could not determine the analysis to perform based on "
                  << "command line options: " << endl << desc << endl;
          exit(COMMAND_LINE_ERROR);
        }
      }
    }
  }

  // -------------------------------------------------------------------------
  // added bcw 7/15/11 - number of numerics or phenotypes might not be the
  // same as the data set - only load those in the numerics/phenotype file
  // if covariate or alternate phenotype file is present, read it first to
  // get the keys for reading the data set instances
  cout << Timestamp()
          << "Checking for covariates and/or alternate phenotype files" << endl;
  vector<string> numericsIds;
  vector<string> phenoIds;
  if(analysisType == SNP_ONLY_ANALYSIS ||
     analysisType == NUMERIC_ONLY_ANALYSIS ||
     analysisType == INTEGRATED_ANALYSIS) {
    if(numericsFilename != "") {
      cout << Timestamp() << "Loading individual IDs from covar file: "
              << numericsFilename << endl;
      if(!LoadIndividualIds(numericsFilename, numericsIds, true)) {
        exit(COMMAND_LINE_ERROR);
      }
      // copy (numericsIds.begin(), numericsIds.end(), ostream_iterator<string> (cout, "\n"));
    }
    if(altPhenotypeFilename != "") {
      cout << Timestamp() << "Loading individual IDs from alternate phenotype file: "
              << altPhenotypeFilename << endl;
      if(!LoadIndividualIds(altPhenotypeFilename, phenoIds, false)) {
        exit(COMMAND_LINE_ERROR);
      }
      // copy(phenoIds.begin(), phenoIds.end(), ostream_iterator<string> (cout, "\n"));
    }
  }
  else {
    cout << Timestamp() << "Covariate and alternate phenotype files not used for "
            << "this analysis type" << endl;
  }

  // -------------------------------------------------------------------------
  // find IDs for loading from the dataset
  cout << Timestamp() << "Determining the IDs to be read from the dataset" << endl;
  vector<string> indIds;
  if(!GetMatchingIds(numericsFilename, altPhenotypeFilename,
                     numericsIds, phenoIds, indIds)) {
    cerr << "ERROR: could not get matching IDs from numeric " <<
            " and/or phenotype files." << endl;
    exit(COMMAND_LINE_ERROR);
  }
  cout << Timestamp() << indIds.size()
          << " individual IDs read from numeric and/or phenotype file(s)"
          << endl;

  // -------------------------------------------------------------------------
  // prepare data for running EC
  cout << Timestamp() << "Preparing data set for EC analysis" << endl;
  Dataset* ds = 0;
  bool datasetLoaded = false;
  switch(analysisType) {
    case SNP_ONLY_ANALYSIS:
      cout << Timestamp() << "Reading SNPs data set" << endl;
      ds = ChooseSnpsDatasetByExtension(snpsFilename);
      datasetLoaded = ds->LoadDataset(snpsFilename, "",
                                      altPhenotypeFilename, indIds);
      break;
    case SNP_CLEAN_ANALYSIS:
      cout << Timestamp() << "Reading CLEAN SNPs data set" << endl;
      ds = ChooseSnpsDatasetByExtension(cleanSnpsFilename, true);
      datasetLoaded = ds->LoadDataset(cleanSnpsFilename, "",
                                      "", indIds);
      break;
    case NUMERIC_ONLY_ANALYSIS:
      cout << Timestamp() << "Reading numerics only data set" << endl;
      ds = new Dataset();
      datasetLoaded = ds->LoadDataset("", numericsFilename,
                                      altPhenotypeFilename, indIds);
      break;
    case INTEGRATED_ANALYSIS:
      cout << Timestamp() << "Reading datasets for integrated analysis" << endl;
      ds = ChooseSnpsDatasetByExtension(snpsFilename);
      datasetLoaded = ds->LoadDataset(snpsFilename, numericsFilename,
                                      altPhenotypeFilename, indIds);
      break;
    case DIAGNOSTIC_ANALYSIS:
      cout << Timestamp() << "Performing SNP diagnostics on the data set" << endl;
      if(snpsFilename == "" || cleanSnpsFilename == "") {
        cerr << "Cannot run diagnostics without a SNP file specified with "
                << "--snp-data or --snp-data-clean flag" << endl;
        exit(COMMAND_LINE_ERROR);
      }
      if(cleanSnpsFilename == "") {
        ds = ChooseSnpsDatasetByExtension(snpsFilename);
      }
      else {
        ds = ChooseSnpsDatasetByExtension(snpsFilename, true);
      }
      ds->LoadDataset(snpsFilename, numericsFilename,
                      altPhenotypeFilename, indIds);
      ds->RunSnpDiagnosticTests(diagnosticLevelsCountsFilename);
      if(diagnosticLevelsCountsFilename != "") {
        ds->WriteLevelCounts(diagnosticLevelsCountsFilename + ".counts");
      }
      // brutal exit!
      cout << argv[0] << " done" << endl;
      exit(0);
      break;
    case NO_ANALYSIS:
      cerr << "Analysis type could not be determined" << endl;
      exit(COMMAND_LINE_ERROR);
    default:
      cerr << "Undefined analysis type: " << analysisType << endl;
      exit(COMMAND_LINE_ERROR);
  }

  if(!datasetLoaded) {
    cerr << "ERROR: Failure to load dataset for analysis" << endl << endl;
    exit(COMMAND_LINE_ERROR);
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
  cout << Timestamp() << "Running EC" << endl;
  EvaporativeCooling ec(ds, vm, analysisType);
  if(!ec.ComputeECScores()) {
    cerr << "ERROR: Failed to calculate EC scores" << endl;
    exit(COMMAND_LINE_ERROR);
  }
  cout << Timestamp() << "EC done" << endl;

  // ---------------------------------------------------------------------------
  // write the scores to the same name as the dataset with
  // <metric>.relieff suffix
  cout << Timestamp() << "Writing EC scores to [" + outputFilesPrefix << ".ec]" << endl;
  ec.WriteAttributeScores(outputFilesPrefix);

  /// write the ReliefF filtered attributes as a new data set
  switch(outputDatasetType) {
    case TAB_DELIMITED_DATASET:
      ds->WriteNewDataset(outputDatasetFilename, TAB_DELIMITED_DATASET);
      break;
    case CSV_DELIMITED_DATASET:
      ds->WriteNewDataset(outputDatasetFilename, CSV_DELIMITED_DATASET);
      break;
    case ARFF_DATASET:
      ds->WriteNewDataset(outputDatasetFilename, ARFF_DATASET);
      break;
    case NO_OUTPUT_DATASET:
    default:
      break;
  }

  // ---------------------------------------------------------------------------
  cout << Timestamp() << "Clean up and shutdown" << endl;
  // delete ds;

  float elapsedTime = (float)(clock() - t) / CLOCKS_PER_SEC;
  cout << Timestamp() << "EC elapsed time " << elapsedTime << " secs" << endl;

  cout << Timestamp() << argv[0] << " done" << endl;

  return 0;
}