/**
 * \file DistanceMetrics.h
 *
 * \brief Distance metrics for ReliefF.
 *
 * \author: Bill White
 * \version 1.0
 *
 * Contact: bill.c.white@gmail.com
 * Created on 3/29/11
 */

#ifndef DISTANCEMETRICS_H
#define	DISTANCEMETRICS_H

/// Forward reference to a DatasetInstance class.
class DatasetInstance;

/***************************************************************************//**
 * Check for a missing discrete value and return value.
 * \param [in] attributeIndex index into the vector of attributes
 * \param [in] dsi1 data set instance 1
 * \param [in] dsi2 data set instance 2
 * \return pair: has missing value?, value for missing (true) or 0.0 (false)
 ******************************************************************************/
std::pair<bool, double> CheckMissing(unsigned int attributeIndex,
                                     DatasetInstance* dsi1,
                                     DatasetInstance* dsi2);
/***************************************************************************//**
 * Check for a missing continuous value and return value.
 * \param [in] attributeIndex index into the vector of attributes
 * \param [in] dsi1 data set instance 1
 * \param [in] dsi2 data set instance 2
 * \return pair: has missing value?, value for missing (true) or 0.0 (false)
 ******************************************************************************/
std::pair<bool, double> CheckMissingNumeric(unsigned int numericIndex,
                                            DatasetInstance* dsi1,
                                            DatasetInstance* dsi2);
/***************************************************************************//**
 * Normalizes a given value of a numeric attribute.
 * Borrowed from Weka 8/18/11
 * \param [in] x value
 * \param [in] minX minimum value for x
 * \param [in] maxX maximum value for x
 * \return normalized value
 ******************************************************************************/
double norm(double x, double minX, double maxX);
/***************************************************************************//**
 * Allele mismatch metric.
 * \param [in] attributeIndex index into the vector of attributes
 * \param [in] dsi1 data set instance 1
 * \param [in] dsi2 data set instance 2
 * \return diff(erence) between attribute values: 0.0, 0.5, 1.0
 ******************************************************************************/
double diffAMM(unsigned int attributeIndex,
               DatasetInstance* dsi1,
               DatasetInstance* dsi2);
/***************************************************************************//**
 * Genotype mismatch metric.
 * \param [in] attributeIndex index into the vector of attributes
 * \param [in] dsi1 data set instance 1
 * \param [in] dsi2 data set instance 2
 * \return diff(erence) between attribute values: 0.0 (same) or 1.0 (not same)
 ****************************************************************************/
double diffGMM(unsigned int attributeIndex,
               DatasetInstance* dsi1,
               DatasetInstance* dsi2);
/***************************************************************************//**
 * Nucleotide count array (NCA) metric.
 * \param [in] attributeIndex index into the vector of attributes
 * \param [in] dsi1 data set instance 1
 * \param [in] dsi2 data set instance 2
 * \return diff(erence) considering nucleotide counts
 ****************************************************************************/
double diffNCA(unsigned int attributeIndex,
               DatasetInstance* dsi1,
               DatasetInstance* dsi2);
/***************************************************************************//**
 * "Manhattan" distance between continuous attributes.
 * \param [in] attributeIndex index into the vector of attributes
 * \param [in] dsi1 data set instance 1
 * \param [in] dsi2 data set instance 2
 * \return absolute value of difference divided by attribute's range
 ******************************************************************************/
double diffManhattan(unsigned int attributeIndex,
                     DatasetInstance* dsi1,
                     DatasetInstance* dsi2);
/***************************************************************************//**
 * Same as "Manhattan" distance but uses method calls versus public variables.
 * \param [in] attributeIndex index into the vector of attributes
 * \param [in] dsi1 data set instance 1
 * \param [in] dsi2 data set instance 2
 * \return absolute value of difference divided by attribute's range
 ******************************************************************************/
double diffPredictedValueTau(DatasetInstance* dsi1,
                             DatasetInstance* dsi2);

#endif	/* DISTANCEMETRICS_H */

