/**
 * \file Insilico.h
 *
 * \brief Common functions for Insilico Lab projects.
 *
 * \author: Bill White
 * \version 1.0
 *
 * Contact: bill.c.white@gmail.com
 * Created on 10/13/11
 */

#ifndef INSILICO_H
#define	INSILICO_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <iterator>

class Dataset;

/// T Y P E D E F S

/// type of discrete attribute values
typedef int AttributeLevel;
/// type of continuous attributes
typedef double NumericLevel;
/// type of instance class labels
typedef int ClassLevel;

/// distance pair type: distance, instance ID
typedef std::pair<double, std::string> DistancePair;
/// vector of distance pairs represents distances to nearest neighbors
typedef std::vector<DistancePair> DistancePairs;
/// distance pairs iterator
typedef DistancePairs::const_iterator DistancePairsIt;

/// Configuration map as an alternative to Boost::program_options
typedef std::map<std::string, std::string> ConfigMap;

/// C O N S T A N T S

/// Error codes.
const static int COMMAND_LINE_ERROR = EXIT_FAILURE;
const static int DATASET_LOAD_ERROR = EXIT_FAILURE;

/// return value for invalid distance
const static int INVALID_DISTANCE = INT_MAX;
/// return value for invalid index into attributes
const static int INVALID_INDEX = INT_MAX;
/// return value for invalid index into attributes
const static unsigned int INVALID_INT_VALUE = UINT_MAX;

/// invalid attribute value
const static AttributeLevel INVALID_ATTRIBUTE_VALUE = INT_MIN;
/// invalid attribute value
const static NumericLevel INVALID_NUMERIC_VALUE = INT_MIN;
/// stored value for missing discrete class
const static ClassLevel INVALID_DISCRETE_CLASS_VALUE = INT_MIN;
/// stored value for missing numeric class
const static NumericLevel INVALID_NUMERIC_CLASS_VALUE = INT_MIN;

/// stored value for missing discrete attribute
const static AttributeLevel MISSING_ATTRIBUTE_VALUE = -9;
/// stored value for missing numeric attribute
const static NumericLevel MISSING_NUMERIC_VALUE = -9;
/// stored value for missing discrete class
const static ClassLevel MISSING_DISCRETE_CLASS_VALUE = -9;
/// stored value for missing numeric class
const static NumericLevel MISSING_NUMERIC_CLASS_VALUE = -9;

/// E N U M S

/**
 * \enum OutputDatasetType.
 * Type of data set to write filtered output.
 */
enum OutputDatasetType {
	TAB_DELIMITED_DATASET, /**< tab-delimited .txt file */
	CSV_DELIMITED_DATASET, /**< comma separated values .csv file */
	ARFF_DATASET, /**< WEKA ARFF format .arff file */
	PLINK_PED_DATASET, /**< PLINK ped/map format */
	PLINK_BED_DATASET, /**< PLINK bed/bim/fam format */
	NO_OUTPUT_DATASET /**< no output data set specified */
};

/**
 * \enum AnalysisType.
 * Type of analysis to perform.
 */
enum AnalysisType
{
  SNP_ONLY_ANALYSIS, /**< discrete analysis */
  NUMERIC_ONLY_ANALYSIS, /**< continuous attributes */
  INTEGRATED_ANALYSIS, /**< discrete and continuous analysis  */
  DIAGNOSTIC_ANALYSIS, /**< diagnostic mode - no ReliefF analysis */
  REGRESSION_ANALYSIS, /**< regression ReliefF analysis */
  DGE_ANALYSIS, /**< digital gene expression (DGE) analysis */
  BIRDSEED_ANALYSIS, /**< Birdseed called SNPs analysis */
  DISTANCE_MATRIX_ANALYSIS, /**< distance matrix calculation */
  DATASET_CONVERSION, /**< convert data set format types */
  NO_ANALYSIS /**< no analysis specified */
};

/**
 * \enum ValueType.
 * Return types for determing a value's type.
 */
enum ValueType
{
  NUMERIC_VALUE, /**< continuous numeric value */
  DISCRETE_VALUE, /**< discrete genotype value */
  MISSING_VALUE, /**< missing value */
  NO_VALUE /**< default no value type */
};

/**
 * \enum AttributeType.
 * Type of attributes that are stored in data set instances.
 */
enum AttributeType
{
  NUMERIC_TYPE, /**< continuous numeric type */
  DISCRETE_TYPE, /**< discrete genotype type */
  NO_TYPE /**< default no type */
};

/**
 * \enum ClassType.
 * Type of classes that are stored in data set instances.
 */
enum ClassType
{
  CONTINUOUS_CLASS_TYPE, /**< continuous numeric type */
  CASE_CONTROL_CLASS_TYPE, /**< discrete case-control type */
  MULTI_CLASS_TYPE, /**< multiclass type */
  NO_CLASS_TYPE /**< default no type */
};

/**
 * \enum AttributeMutationType.
 * Type of attribute mutation.
 */
enum AttributeMutationType
{
  TRANSITION_MUTATION, /**< transition within family */
  TRANSVERSION_MUTATION, /**< transversion between families */
  UNKNOWN_MUTATION /**< unknown - no allele information */
};

/**
 * \enum RandomJungleTreeType.
 * Type random jungle trees.
 */
enum RandomJungleTreeType
{
  UNKNOWN_TREE_TYPE=0, /**< place holder = 0 */
  NOMINAL_NUMERIC_TREE, /**< classification trees, numeric attributes (integers) */
  NOMINAL_NOMINAL_TREE, /**< classification trees, discrete attributes (0/1/2) */
  NUMERIC_NUMERIC_TREE, /**< regression trees, numeric attributes (doubles) */
  NUMERIC_NOMINAL_TREE, /**< regression trees, discrete attributes (0/1/2) */
  NOMINAL_NUMERIC_FLOATS /**< classification trees, numeric attributes (doubles) */
};

/**
 * \enum RandomJungleRunMode.
 * Run mode for random jungle.
 */
enum RandomJungleRunMode
{
	UNKNOWN_RUN_MODE, /**< unknown run mode */
	LIBRARY_RUN_MODE, /**< call Random Jungle through C++ library calls */
	SYSTEM_CALL_RUN_MODE /**< call Random Jungle through C system() call */
};

static std::map<std::string, std::string> datasetTypeToExt;

/***************************************************************************//**
 * Return random jungle tree type from the class and attribute types
 * \param [in] attributeType attribute data type
 * \param [in] classType class data type
 * \return Random Jungle tree type
 ******************************************************************************/
RandomJungleTreeType DetermineRandomJungleTreeType(AttributeType attributeType,
		ClassType classType);
/***************************************************************************//**
 * Return a timestamp string for logging purposes.
 * \return fixed-length, formatted timestamp as a string
 ******************************************************************************/
std::string Timestamp();
/***************************************************************************//**
 * Determines the data set type to instantiate based on the
 * passed type string or data set filenames's extension.
 * \param [in] snpsFilename SNP data set filename
 * \param [in] snpsFileType SNP data set file type (overrides detect by extension)
 * \return pointer to new data set or NULL if could not match a data set type
 ******************************************************************************/
Dataset* ChooseSnpsDatasetByType(std::string snpsFilename,
		std::string snpsFileType="");
/***************************************************************************//**
 * Loads the individual (instance) IDs from the numerics file.
 * Returns the IDs through reference parameter retIds.
 * \param [in] filename filename that contains numerics IDs
 * \param [out] vector of individual (instance) IDs (strings)
 * \return success
 ******************************************************************************/
bool LoadNumericIds(std::string filename,
                    std::vector<std::string>& retIds);
/***************************************************************************//**
 * Loads the individual (instance) IDs from the numerics file.
 * Returns the IDs through reference parameter retIds.
 * \param [in] filename filename that contains numerics IDs
 * \param [out] vector of individual (instance) IDs (strings)
 * \return success
 ******************************************************************************/
bool LoadPhenoIds(std::string filename,
                  std::vector<std::string>& retIds);
/***************************************************************************//**
 * Return matching IDs from numeric and/or phenotype file IDs
 * \param [in] numericsFilename name of the PLINK covar format file
 * \param [in] altPhenotypeFilename name of the alternate pheno file PLINK
 * \param [in] numericsIds covar format file ids
 * \param [in] phenoIds alternate phenotype file ids
 * \param [out] matchingIds ids that match between numerics and phenotypes
 * \return success
 ******************************************************************************/
bool GetMatchingIds(std::string numericsFilename,
                    std::string altPhenotypeFilename,
                    std::vector<std::string> numericsIds,
                    std::vector<std::string> phenoIds,
                    std::vector<std::string>& matchingIds);
/***************************************************************************//**
 * Detect the class type by reading the specified column from a whitespace-
 * delimited text file.
 * \param [in] filename whitespace-delimited text file name
 * \param [in] classColumn the column containing the class values
 * \param [in] heasHeader does the file have a header line?
 * \return ClassType defined in Dataset.h
 ******************************************************************************/
ClassType DetectClassType(std::string filename, int classColumn, bool hasHeader);
/***************************************************************************//**
 * Get the parameter value from the configuration map key.
 * \param [in] configMap reference to a configuration map
 * \param [in] key parameter name
 * \param [out] parameter value
 * \return true if key found, false if not found
 ******************************************************************************/
bool GetConfigValue(ConfigMap& configMap, std::string key, std::string& value);
/***************************************************************************//**
 * Get the full filename without the extension.
 * \param [in] fullFilename complete filename
 * \return path/filename without extension
 ******************************************************************************/
std::string GetFileBasename(std::string fullFilename);
/***************************************************************************//**
 * Get the filename extension.
 * \param [in] fullFilename complete filename
 * \return filename extension
 ******************************************************************************/
std::string GetFileExtension(std::string fullFilename);
/***************************************************************************//**
 * Print a vector of T values with optional title.
 * \param [in] vec vector of T type values
 * \param [in] title optional title to print before the vector
 ******************************************************************************/
template <class T> void PrintVector(std::vector<T> vec, std::string title="");
template <class T> void PrintVector(std::vector<T> vec, std::string title) {
  if(title != "") {
    std::cout << title << ": ";
  }
  std::cout << "[ ";
  std::copy(vec.begin(), vec.end(), std::ostream_iterator<T>(std::cout, " "));
  std::cout << "]" << std::endl;
}

/// protected log function returns 0 for 0
double ProtectedLog(double x);

#endif	/* INSILICO_H */
