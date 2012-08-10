/*
 * EvaporativeCoolingCLINewAPI.cpp - Bill White - 1/30/11
 *
 * Evaporative Cooling (EC) Optimization of Information Free Energy -
 * Command Line Interface (CLI) Using the New Map Parameters API
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

#include <boost/program_options.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/parsers.hpp>

#include "EvaporativeCooling.h"
#include "Insilico.h"
#include "Dataset.h"
#include "DgeData.h"

using namespace std;
namespace po = boost::program_options;

int main(int argc, char** argv) {

  // command line processing variables: defaults and storage for boost
  // data set files
  string configFilename = "";
  string snpsFilename = "";
  string snpExclusionFile = "";
  string numericsFilename = "";
  string dgeCountsFilename = "";
  string dgePhenosFilename = "";
  string dgeNormsFilename = "";
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
  unsigned int ecNumTarget = 0;
  unsigned int ecIterNumToRemove = 0;
  unsigned int ecIterPercentToRemove = 0;

  // declare the supported options
  po::options_description desc("Allowed options");
  desc.add_options()
          ("help", "produce help message")
          ("verbose", "verbose output")
          (
           "config-file,c",
           po::value<string > (&configFilename),
           "read configuration options from file - command line overrides these"
           )
          (
           "snp-data,s",
           po::value<string > (&snpsFilename),
           "read SNP attributes from genotype filename: txt, ARFF, plink (map/ped, binary, raw)"
           )
          (
           "numeric-data,n",
           po::value<string > (&numericsFilename),
           "read continuous attributes from PLINK-style covar file"
           )
          (
           "alternate-pheno-file,a",
           po::value<string > (&altPhenotypeFilename),
           "specifies an alternative phenotype/class label file; one value per line"
           )
           (
            "ec-algorithm-steps,g",
            po::value<string > (&ecAlgorithmSteps)->default_value(ecAlgorithmSteps),
            "EC steps to run (all|rj|rf)"
            )
           (
            "ec-num-target,t",
            po::value<unsigned int>(&ecNumTarget)->default_value(ecNumTarget),
            "EC N_target - target number of attributes to keep"
            )
           (
            "ec-iter-remove-n,r",
            po::value<unsigned int>(&ecIterNumToRemove)->default_value(ecIterNumToRemove),
            "Evaporative Cooling number of attributes to remove per iteration"
            )
           (
            "ec-iter-remove-percent,p",
            po::value<unsigned int>(&ecIterPercentToRemove),
            "Evaporative Cooling precentage of attributes to remove per iteration"
            )
          (
           "out-dataset-filename,O",
           po::value<string > (&outputDatasetFilename),
           "write a new tab-delimited data set with EC filtered attributes"
           )
          (
           "out-files-prefix,o",
           po::value<string > (&outputFilesPrefix)->default_value(outputFilesPrefix),
           "use prefix for all output files"
           )
           (
            "snp-metric,S",
            po::value<string > (&snpMetric)->default_value(snpMetric),
            "metric for determining the difference between SNPs (gm|am)"
            )
           (
            "numeric-metric,N",
            po::value<string > (&numMetric)->default_value(numMetric),
            "metric for determining the difference between numeric attributes (manhattan=|euclidean)"
            )
          (
           "rj-num-trees,j",
           po::value<uli_t > (&rjNumTrees)->default_value(rjNumTrees),
           "Random Jungle number of trees to grow"
           )
           (
            "snp-exclusion-file,x",
            po::value<string > (&snpExclusionFile),
            "file of SNP names to be excluded"
            )
          (
           "k-nearest-neighbors,k",
           po::value<unsigned int>(&k)->default_value(k),
           "set k nearest neighbors"
           )
          (
           "number-random-samples,m",
           po::value<unsigned int>(&m)->default_value(m),
           "number of random samples (0=all|1 <= n <= number of samples)"
           )
          (
           "weight-by-distance-method,b",
           po::value<string > (&weightByDistanceMethod)->default_value(weightByDistanceMethod),
           "weight-by-distance method (equal|one_over_k|exponential)"
           )
          (
           "weight-by-distance-sigma",
           po::value<double>(&weightByDistanceSigma)->default_value(weightByDistanceSigma),
           "weight by distance sigma"
           )
          (
           "diagnostic-tests,d",
           po::value<string > (&diagnosticLogFilename),
           "performs diagnostic tests and sends output to filename without running EC"
           )
          (
           "diagnostic-levels-file,D",
           po::value<string > (&diagnosticLevelsCountsFilename),
           "write diagnostic attribute level counts to filename"
           )
           (
            "dge-counts-data",
            po::value<string > (&dgeCountsFilename),
            "read digital gene expression counts from text file"
            )
            (
             "dge-phenos-data",
             po::value<string > (&dgePhenosFilename),
             "read digital gene expression phenotypes from text file"
             )
             (
              "dge-norm-factors",
              po::value<string > (&dgeNormsFilename),
              "read digital gene expression normalization factors from text file"
              )
          ;

  // parse the command line into a map
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  // po::store(po::command_line_parser(argc, argv).options(desc).positional(pd).run(), vm);
  po::notify(vm);

  if(vm.count("help") || (argc == 1)) {
    cerr << desc << endl;
    exit(COMMAND_LINE_ERROR);
  }

  // ---------------------------------------------------------------------------
  cout << Timestamp() << argv[0] << " starting" << endl;
  clock_t t;
  t = clock();

  // ---------------------------------------------------------------------------
  cout << Timestamp() << "Processing command line arguments" << endl;

  // read config file if specified
  if(vm.count("config-file")) {
  	ifstream configStream(configFilename.c_str());
  	if (!configStream.is_open()) {
  		cerr << "ERROR: Could not open configuration file: "
  				<< configFilename << endl;
  		exit(EXIT_FAILURE);
  	}
  	cout << Timestamp() << "Reading configuration options from: "
  			<< configFilename << endl;
    po::store(po::parse_config_file(configStream, desc), vm);
    po::notify(vm);
    configStream.close();
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
    cout << Timestamp() << "Writing ReliefF filtered data set to ["
            << outputDatasetFilename
            << "]" << endl;
  }

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
    if(vm.count("snp-data") && !vm.count("numeric-data")) {
			cout << Timestamp() << "SNP-only analysis requested" << endl;
			analysisType = SNP_ONLY_ANALYSIS;
    } else {
      if(!vm.count("snp-data") && (vm.count("numeric-data"))) {
        cout << Timestamp() << "Numeric-only analysis requested" << endl;
        analysisType = NUMERIC_ONLY_ANALYSIS;
        // must have an alternate phenotype file for numeric only
        if(!vm.count("alternate-pheno-file")) {
          cerr << "An alternate phenotype file must be specified with the "
                  << "--alternate-pheno-file option for numeric-only or "
                  << "digital gene expression data"
                  << endl;
          exit(COMMAND_LINE_ERROR);
        }
      } else {
        if(vm.count("snp-data") && vm.count("numeric-data")) {
          cout << Timestamp() << "Integrated analysis requested" << endl;
          analysisType = INTEGRATED_ANALYSIS;
        } else {
          if(vm.count("dge-counts-data") && vm.count("dge-phenos-data")) {
            cout << Timestamp() << "DGE analysis requested" << endl;
            analysisType = DGE_ANALYSIS;
          }
          else {
						cerr << "ERROR: Could not determine the analysis to perform based on "
										<< "command line options: " << endl << desc << endl;
						exit(COMMAND_LINE_ERROR);
          }
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
      if(!LoadNumericIds(numericsFilename, numericsIds)) {
        exit(COMMAND_LINE_ERROR);
      }
      // copy (numericsIds.begin(), numericsIds.end(), ostream_iterator<string> (cout, "\n"));
    }
    if(altPhenotypeFilename != "") {
      cout << Timestamp() << "Loading individual IDs from alternate phenotype file: "
              << altPhenotypeFilename << endl;
      if(!LoadPhenoIds(altPhenotypeFilename, phenoIds)) {
        exit(COMMAND_LINE_ERROR);
      }
      // copy(phenoIds.begin(), phenoIds.end(), ostream_iterator<string> (cout, "\n"));
    }
  } else {
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
	DgeData* dge;
  bool datasetLoaded = false;
  switch(analysisType) {
    case SNP_ONLY_ANALYSIS:
      cout << Timestamp() << "Reading SNPs data set" << endl;
      ds = ChooseSnpsDatasetByType(snpsFilename);
      datasetLoaded = ds->LoadDataset(snpsFilename, "",
                                      altPhenotypeFilename, indIds);
      break;
    case NUMERIC_ONLY_ANALYSIS:
      cout << Timestamp() << "Reading numerics only data set" << endl;
      ds = new Dataset();
      datasetLoaded = ds->LoadDataset("", numericsFilename,
                                      altPhenotypeFilename, indIds);
      break;
    case INTEGRATED_ANALYSIS:
      cout << Timestamp() << "Reading datasets for integrated analysis" << endl;
      ds = ChooseSnpsDatasetByType(snpsFilename);
      datasetLoaded = ds->LoadDataset(snpsFilename, numericsFilename,
                                      altPhenotypeFilename, indIds);
      break;
    case DGE_ANALYSIS:
      cout << Timestamp() << "Reading numerics data set from digital gene "
      << "expression (DGE) data" << endl;
    	dge = new DgeData();
    	if(dge->LoadData(dgeCountsFilename, dgeNormsFilename)) {
    		ds = new Dataset();
    		datasetLoaded = ds->LoadDataset(dge);
    	}
    	break;
    case DIAGNOSTIC_ANALYSIS:
      cout << Timestamp() << "Performing SNP diagnostics on the data set" << endl;
      if(snpsFilename == "") {
        cerr << "Cannot run diagnostics without a SNP file specified with "
                << "--snp-data or --snp-data-clean flag" << endl;
        exit(COMMAND_LINE_ERROR);
      }
      ds = ChooseSnpsDatasetByType(snpsFilename);
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
  if(analysisType == SNP_ONLY_ANALYSIS) {
    ds->PrintStats();
  } else {
    if(analysisType == NUMERIC_ONLY_ANALYSIS ||
       analysisType == INTEGRATED_ANALYSIS ||
       analysisType == DGE_ANALYSIS) {
      ds->PrintNumericsStats();
    }
  }

  // ---------------------------------------------------------------------------
  // FINALLY! run EC algorithm
  cout << Timestamp() << "Running EC" << endl;
  /// Test the new API with a ConfigMap rather than Boost Program Options
  /// ConfigMap defined in Insilico.h: std::map<std::string, std::string>
  ConfigMap configMap;
  string tempKey;
  stringstream ss;
  if(vm.count("verbose")) {
  	configMap.insert(make_pair("verbose", "true"));
  }
  else {
  	configMap.insert(make_pair("verbose", ""));
  }
  configMap.insert(make_pair("ec-algorithm-steps", ecAlgorithmSteps));
  ss << ecNumTarget;
  configMap.insert(make_pair("ec-num-target", ss.str()));
  ss.str("");
  ss << ecIterNumToRemove;
  configMap.insert(make_pair("ec-iter-remove-n",ss.str()));
  ss.str("");
  ss << ecIterPercentToRemove;
  configMap.insert(make_pair("ec-iter-remove-percent",ss.str()));
  ss.str("");
  configMap.insert(make_pair("out-files-prefix", outputFilesPrefix));
  configMap.insert(make_pair("snp-metric", snpMetric));
  configMap.insert(make_pair("numeric-metric", numMetric));
  ss << rjNumTrees;
  configMap.insert(make_pair("rj-num-trees",ss.str()));
  ss.str("");
  configMap.insert(make_pair("snp-exclusion-file", snpExclusionFile));
  ss << k;
  configMap.insert(make_pair("k-nearest-neighbors",ss.str()));
  ss.str("");
  ss << m;
  configMap.insert(make_pair("number-random-samples", ss.str()));
  ss.str("");
  configMap.insert(make_pair("weight-by-distance-method",weightByDistanceMethod));
  ss << weightByDistanceSigma;
  configMap.insert(make_pair("weight-by-distance-sigma", ss.str()));
  EvaporativeCooling ec(ds, configMap, analysisType);

  if(!ec.ComputeECScores()) {
    cerr << "ERROR: Failed to calculate EC scores" << endl;
    exit(COMMAND_LINE_ERROR);
  }

  cout << Timestamp() << "EC done" << endl;

  // ---------------------------------------------------------------------------
  // write the scores to the same name as the dataset with
  // <metric>.relieff suffix
  string resultsFilename = outputFilesPrefix;
  switch(ec.GetAlgorithmType()) {
    case EC_ALG_ALL:
      resultsFilename += ".ec";
      break;
    case EC_ALG_RJ:
      resultsFilename += ".ec.rj";
      break;
    case EC_ALG_RF:
      resultsFilename += ".ec.rf";
      break;
    default:
      // we should not get here by the CLI front end but it is possible to call
      // this from other programs in the future or when used as a library
      cerr << "ERROR: Attempting to write attribute scores before the algorithm "
              << "type was determined. " << endl;
      return false;
  }

  cout << Timestamp() << "Writing EC scores to [" + resultsFilename + "]" << endl;
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

  /// Remove temporary Random Jungle files if they exist
  if((ec.GetAlgorithmType() == EC_ALG_ALL) || (ec.GetAlgorithmType() == EC_ALG_RJ)) {
  	cout << Timestamp() << "Removing temporary RandomJungle files" << endl;
		vector<string> tempFilenames;
		tempFilenames.push_back(outputFilesPrefix + ".log");
		tempFilenames.push_back(outputFilesPrefix + ".verbose");
		tempFilenames.push_back(outputFilesPrefix + ".importance");
		tempFilenames.push_back(outputFilesPrefix + ".confusion");
		tempFilenames.push_back(outputFilesPrefix + ".confusion2");
		for(vector<string>::const_iterator it=tempFilenames.begin();
				it != tempFilenames.end(); ++it) {
			unlink((*it).c_str());
		}
  }

  float elapsedTime = (float) (clock() - t) / CLOCKS_PER_SEC;
  cout << Timestamp() << "EC elapsed time " << elapsedTime << " secs" << endl;

  cout << Timestamp() << argv[0] << " done" << endl;

  return 0;
}
