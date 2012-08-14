/**
 * \class ChiSquared
 *
 * \brief Chi-squared attribute ranking algorithm.
 *
 * ChiSquared algorithm interface. For performing chi-squared tests of
 * association between an attribute and its class across all instances
 * in a data set.
 **
 * \author Bill White
 * \version 1.0
 *
 * Contact: bill.c.white@gmail.com
 * Created on: 6/15/05
 *
 * Modified to implement new AttributeRanker interface.
 */

#ifndef CHI_SQUARED_H
#define CHI_SQUARED_H

#include <vector>
#include <fstream>

#include "AttributeRanker.h"
#include "Dataset.h"

class ChiSquared: public AttributeRanker
{
public:
  /*************************************************************************//**
   * Construct an chi-squared algorithm object.
   * \param [in] ds pointer to a Dataset object
   ****************************************************************************/
  ChiSquared(Dataset* ds);
  ~ChiSquared();
  /*************************************************************************//**
   * For each attribute, calculate chi-squared and associated p-value. Return
   * in a vector of pairs indexed by attribute index.
   * \return vector of pairs of chi-squared scores and associated p-values
   ****************************************************************************/
  const std::vector<std::pair<double, double> >& ComputeScoresWithPValues();
  /*************************************************************************//**
   * For the attribute at the specified index, calculate the chi-squared and
   * associated p-value. Return as a pair.
   * \param [in] index index into the attributes of the data set
   * \return pairs of chi-squared score and associated p-value for the attribute
   ****************************************************************************/
  std::pair<double, double> ComputeScore(unsigned int index);
  /// Print calculation tables
  void PrintTables();
  /*************************************************************************//**
   * Print the scores to a stream.
   * \param [in] outStream reference to an output stream
   * \param [in] topN top number of attributes to print
   ****************************************************************************/
  void PrintScoresWithPValues(std::ofstream& outStream, unsigned int topN = 0);
  /*************************************************************************//**
   * Print the scores to a stream.
   * \param [in] outFilename filename to write scores to
   * \param [in] topN top number of attributes to print
   ****************************************************************************/
  void WriteScoresWithPValues(std::string outFilename, unsigned int topN = 0);
  /// Get the observed frequencies table as a vector of vector of doubles
  std::vector<std::vector<double> > GetFrequencyCounts() {
    return observedFreqTable;
  }
  // added to support AttributeRanker interface - 8/13/12
  AttributeScores ComputeScores();
private:
  /*************************************************************************//**
   * Private method to setup the chi-squared contingency tables for a particular
   * attribute.
   * \param [in] attributeIndex attribute index
   ****************************************************************************/
  void PrepareForAttribute(unsigned int attributeIndex);
   /// Clear calculation tables
  void ClearTables();

  /// number of levels in the attributes
  unsigned int numLevels;
  /// number of classes in the instances
  unsigned int numClasses;
  /// observed frequencies
  std::vector<std::vector<double> > observedFreqTable;
  //// expected frequencies
  std::vector<std::vector<double> > expectedContingencyTable;
  /// chi squared computed values
  std::vector<std::vector<double> > chiSquaredValues;
  /// chi-squared value, p-value for each attribute
  std::vector<std::pair<double, double> > scoresPvalues;
};

#endif
