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
#include <string>
#include <vector>

/// Forward reference to Dataset class.
#include "Dataset.h"

typedef std::map<std::string, std::string> ConfigMap;

/**
 * \enum AnalysisType.
 * Type of analysis to perform.
 */
typedef enum
{
  SNP_ONLY_ANALYSIS, /**< discrete analysis */
  SNP_CLEAN_ANALYSIS, /**< discrete analysis - no filtering */
  NUMERIC_ONLY_ANALYSIS, /**< continuous attributes */
  INTEGRATED_ANALYSIS, /**< discrete and continuous analysis  */
  DIAGNOSTIC_ANALYSIS, /**< diagnostic mode - no ReliefF analysis */
  REGRESSION_ANALYSIS, /**< regression ReliefF analysis */
  DGE_ANALYSIS, /**< digital gene expression (DGE) analysis */
  BIRDSEED_ANALYSIS, /**< Birdseed called SNPs analysis */
  DISTANCE_MATRIX_ANALYSIS, /**< distance matrix calculation */
  NO_ANALYSIS /**< no analysis specified */
} AnalysisType;

/// Error codes.
const static int COMMAND_LINE_ERROR = EXIT_FAILURE;

/***************************************************************************//**
 * Return a timestamp string for logging purposes.
 * \return fixed-length, formatted timestamp as a string
 ******************************************************************************/
std::string Timestamp();
/***************************************************************************//**
 * Determines the data set type to instantiate based on the
 * data set filenames's extension.
 * \param [in] snpsFilename SNP data set filename
 * \return pointer to new dataset or NULL if could not match filename extension
 ******************************************************************************/
Dataset* ChooseSnpsDatasetByExtension(std::string snpsFilename);
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
///
bool GetConfigValue(ConfigMap& configMap, std::string key, std::string& value);

#endif	/* INSILICO_H */
