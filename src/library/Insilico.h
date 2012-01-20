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
class Dataset;

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
 * \param isCleanSnps is this a CleanSnpsDataset
 * \return pointer to new dataset or NULL if could not match filename extension
 ******************************************************************************/
Dataset* ChooseSnpsDatasetByExtension(std::string snpsFilename,
                                      bool isCleanSnps = false);
/***************************************************************************//**
 * Loads the individual (instance) IDs from the numerics or alternate
 * phenotype file. Returns the IDs through reference parameter retIds.
 * \param [in] filename filename that contains covar or pheno file IDs
 * \param [out] vector of individual (instance) IDs (strings)
 * \param [in] does the file have a header row?
 * \return success
 ******************************************************************************/
bool LoadIndividualIds(std::string filename,
                       std::vector<std::string>& retIds,
                       bool hasHeader);
/***************************************************************************//**
 * Return matching IDs from numeric and/or phenotype file IDs
 * \param [in] numericsFilename
 * \param [in] altPhenotypeFilename
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

#endif	/* INSILICO_H */
