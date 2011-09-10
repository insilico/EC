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

#include <boost/program_options.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/parsers.hpp>

#include "EvaporativeCooling.h"
#include "ArffDataset.h"
#include "StringUtils.h"
#include "FilesystemUtils.h"
#include "ArffDataset.h"
#include "PlinkDataset.h"
#include "PlinkRawDataset.h"
#include "PlinkBinaryDataset.h"
#include "Debugging.h"

using namespace std;
using namespace insilico;
namespace po = boost::program_options;

// forward declarations of main program helper functions
Dataset* ChooseSnpsDatasetByExtension(string snpsFilename);
bool LoadIndividualIds(string filename, vector<string>& retIds, bool hasHeader);

int main(int argc, char** argv) {

  // ---------------------------------------------------------------------------
  cout << argv[0] << " starting..." << endl;

  // ---------------------------------------------------------------------------
  cout << "\tProcessing command line arguments..." << endl;

  // command line processing variables: defaults and storage for boost
  string snpsFilename = "";
  string numericsFilename = "";
  unsigned int k = 10;
  unsigned int m = 0;
  bool wbd = false;
  double sigma = 0.0;
  string snpMetric = "gm";
  string numMetric = "manhattan";
  string diagnosticLogFilename = "";
  string diagnosticLevelsCountsFilename = "";
  unsigned int iterNumToRemove = 0;
  unsigned int iterPercentToRemove = 0;
  unsigned int ecIterNumToRemove = 1;
  unsigned int ecIterPercentToRemove = 0;
  string snpExclusionFile = "";
  bool doRecodeA = false;
  string altPhenotypeFilename = "";
  bool verbose = false;
  unsigned int ecNumTarget = 0;
  unsigned int rjNumTrees = 1000;

  // declare the supported options
  po::options_description desc("Allowed options");
  desc.add_options()
          ("help", "produce help message")
          (
           "ec-num-target",
           po::value<unsigned int>(&ecNumTarget)->default_value(ecNumTarget),
           "EC N_target - target number of attributes to keep"
           )
          (
           "rj-num-trees",
           po::value<unsigned int>(&rjNumTrees)->default_value(rjNumTrees),
           "Random Jungle number of trees to grow"
           )
          (
           "alternate-pheno-file,a",
           po::value<string> (&altPhenotypeFilename),
           "specifies an alternative phenotype/class label file; one value per line"
           )
          (
           "snp-data,d",
           po::value<string> (&snpsFilename),
           "read SNP attributes from genotype filename: txt, ARFF, plink (map/ped, binary, raw)"
           )
          (
           "snp-metric",
           po::value<string> (&snpMetric)->default_value(snpMetric),
           "metric for determining the difference between SNPs (gm=default|am)"
           )
          (
           "numeric-data,n",
           po::value<string> (&numericsFilename),
           "read SNP attributes from genotype filename: txt, ARFF, plink (map/ped, binary, raw)"
           )
            (
           "numeric-metric",
           po::value<string> (&numMetric)->default_value(numMetric),
           "metric for determining the difference between numeric attributes (manhattan=default|euclidean)"
           )
          (
           "diagnostic-tests,g",
           po::value<string> (&diagnosticLogFilename),
           "performs diagnostic tests and sends output to filename without running Relief-F"
           )
          (
           "iter-remove-n,i",
           po::value<unsigned int>(&iterNumToRemove),
           "iterative ReliefF number of attributes to remove per iteration"
           )
          (
           "ec-iter-remove-n",
           po::value<unsigned int>(&ecIterNumToRemove)->default_value(ecIterNumToRemove),
           "Evaporative Cooling number of attributes to remove per iteration"
           )
          (
           "k-nearest-neighbors,k",
           po::value<unsigned int>(&k)->default_value(k),
           "set k nearest neighbors"
           )
          (
           "diagnostic-levels-file,l",
           po::value<string > (&diagnosticLevelsCountsFilename),
           "write diagnostic attribute level counts to filename"
           )
          (
           "number-random-samples,m",
           po::value<unsigned int>(&m)->default_value(m),
           "number of random samples (default=0=all|1 <= n <= number of samples)"
           )
          (
           "iter-remove-percent,p",
           po::value<unsigned int>(&iterPercentToRemove),
           "iterative ReliefF precentage of attributes to remove per iteration"
           )
          (
           "ec-iter-remove-percent,p",
           po::value<unsigned int>(&ecIterPercentToRemove),
           "Evaporative Cooling precentage of attributes to remove per iteration"
           )
          (
           "recode-a,r",
           po::value<bool>(&doRecodeA)->default_value(doRecodeA),
           "do a plink recodeA encoding to insure genotype data values (0=no=default|1=yes)"
           )
          (
           "weight-by-distance-sigma,s",
           po::value<double>(&sigma),
           "weight by distance sigma (default=20.0)"
           )
          (
           "weight-by-distance,w",
           po::value<bool>(&wbd)->default_value(wbd),
           "weight differences by distance (0=no|1=yes)"
           )
          (
           "verbose,v",
           po::value<bool>(&verbose)->default_value(verbose),
           "verbose output to stdout (0=no|1=yes)"
           )
          (
           "snp-exclusion-file,x",
           po::value<string>(&snpExclusionFile),
           "file of SNP names to be excluded"
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
      cerr << "Integrated data not supported in this version of EC." << endl;
      exit(1);
    }
    if(vm.count("snp-data") && !vm.count("numeric-data")) {
      cout << "\t\tSNP-only analysis requested." << endl;
      analysisType = SNP_ONLY_ANALYSIS;
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
        cerr << "ERROR: Could not determine the analysis to do based on "
                << "command line options: " << endl << desc << endl;
        exit(1);
      }
    }
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
     analysisType == NUMERIC_ONLY_ANALYSIS) {
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
                                      altPhenotypeFilename, indIds);
      break;
    case NUMERIC_ONLY_ANALYSIS:
      cout << "\tReading numerics only data set" << endl;
      ds = new Dataset();
      datasetLoaded = ds->LoadDataset("", doRecodeA, numericsFilename,
                                      altPhenotypeFilename, indIds);
      break;
    case DIAGNOSTIC_ANALYSIS:
      cout << "\tPerforming SNP diagnostics on the data set" << endl;
      if(snpsFilename == "") {
        cerr << "Cannot run diagnostics without a SNP file specified with "
                << "--snp-data flag." << endl;
        exit(1);
      }
      ds = ChooseSnpsDatasetByExtension(snpsFilename);
      ds->LoadDataset(snpsFilename, doRecodeA, numericsFilename,
                      altPhenotypeFilename, indIds);
      ds->RunSnpDiagnosticTests(diagnosticLevelsCountsFilename);
      if(diagnosticLevelsCountsFilename != "") {
        ds->WriteLevelCounts(diagnosticLevelsCountsFilename + ".counts");
      }
      // brutal exit!
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
  if(analysisType == SNP_ONLY_ANALYSIS) {
    ds->PrintStats();  
  }
  else {
    if(analysisType == NUMERIC_ONLY_ANALYSIS) {
      ds->PrintNumericsStats();
    }
  }

  // ---------------------------------------------------------------------------
  // FINALLY! run EC algorithm
  cout << "\tRunning EC..." << endl;
  EvaporativeCooling ec(ds, vm);
  map<string, double> ecScores;
  if(!ec.ComputeECScores()) {
    cerr << "ERROR: Failed to calculate EC scores." << endl;
    exit(1);
  }

  // ---------------------------------------------------------------------------
  // write the scores to the same name as the dataset with
  // <metric>.relieff suffix
  string fileToWriteOutput;
  if(snpsFilename != "") {
    fileToWriteOutput = snpsFilename + "." + snpMetric + "." + numMetric;
  } else {
    fileToWriteOutput = numericsFilename + "." + snpMetric + "." + numMetric;
  }
  cout << "\tWriting EC scores to [" + fileToWriteOutput << ".ec]" << endl;
  ec.WriteAttributeScores(fileToWriteOutput);

  // ---------------------------------------------------------------------------
  cout << "\tClean up and shutdown." << endl;
  // delete ds;

  cout << argv[0] << " done." << endl;

  return 0;
}

/*****************************************************************************
 * Function: ChooseSnpsDatasetByExtension
 *
 * IN:    SNP data set filename
 * OUT:   Pointer to a newly-instantiated dataset
 *        or NULL if could not find a matching data set type
 *
 * Determines the data set type to instantiate based on the
 * data set filenames's extension.
 ****************************************************************************/
Dataset* ChooseSnpsDatasetByExtension(string snpsFilename) {
  string fileExt = "";
  fileExt = GetFileExtension(snpsFilename);
  // cout << "File extension: " << fileExt << endl;
  Dataset* ds = 0;
  if(fileExt == "arff") {
    cout << "\t\tARFF";
    ds = new ArffDataset();
  } else {
    if(fileExt == "tab" || fileExt == "txt") {
      cout << "\t\tWhitespace-delimited";
      ds = new Dataset();
    } else {
      if(fileExt == "ped" || fileExt == "map") {
        cout << "\t\tPlink map/ped";
        ds = new PlinkDataset();
      } else {
        if(fileExt == "raw") {
          cout << "\t\tPlink raw";
          ds = new PlinkRawDataset();
        } else {
          if(fileExt == "bed" || fileExt == "bim" || fileExt == "fam") {
            cout << "\t\tPlink binary";
            ds = new PlinkBinaryDataset();
          } else {
            cerr << endl;
            cerr << "ERROR: Cannot determine data set type by extension: "
                    << fileExt << endl;
            return 0;
          }
        }
      }
    }
  }
  cout << "." << endl;

  return ds;
}

/*****************************************************************************
 * Function: LoadIndividualIds
 *
 * IN:    filename that contains ID in second column (covar or pheno file)
 * INOUT: vector of individual (instance) IDs (strings)
 * OUT:   success
 *
 * Loads the individual (instance) IDs from the numerics or alternate 
 * phenotype file. Returns the IDs through reference parameter retIds.
 ****************************************************************************/
bool LoadIndividualIds(string filename, vector<string>& retIds,
                       bool hasHeader) {
  ifstream dataStream(filename.c_str());
  if(!dataStream.is_open()) {
    cerr << "ERROR: Could not open ID file: "
            << filename << endl;
    return false;
  }

  // temporary string for reading file lines
  string line;

  // strip the header if there is one and validate three tab-delimited fields
  if(hasHeader) {
    // read the header
    getline(dataStream, line);
    vector<string> numNames;
    split(numNames, line);
    if(numNames.size() < 3) {
      cerr << "ERROR: ID file must have at least three columns: "
              << "FID IID VAR1 . . . VARn" << endl;
      return false;
    }
  }

  // read each line of the file and get the first tab-delimited field as the
  // individual's ID; insure each line has three tab-delimited fields
  retIds.clear();
  map<string, bool> idsSeen;
  unsigned int lineNumber = 0;
  retIds.clear();
  while(getline(dataStream, line)) {
    ++lineNumber;
    vector<string> fieldsStringVector;
    split(fieldsStringVector, line);
    if(fieldsStringVector.size() < 3) {
      cerr << "ERROR: ID file must have at least three columns: "
              << "FID IID VAR1 ... VARn" << endl;
      return false;
    }
    string ID = trim(fieldsStringVector[0]);
    if(idsSeen.find(ID) == idsSeen.end()) {
      idsSeen[ID] = true;
      retIds.push_back(ID);
    } else {
      cerr << "\t\t\tWARNING: Duplicate ID [" << ID << "] detected and "
              << "skipped on line [" << lineNumber << "]." << endl;
    }
  }
  dataStream.close();

  return true;
}
