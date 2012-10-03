/**
 * \class Edger
 *
 * \brief Edger attribute ranking algorithm.
 *
 * Edger algorithm interface.
 **
 * \author Bill White
 * \version 1.0
 *
 * Contact: bill.c.white@gmail.com
 * Created on: 10/2/12
 */

#ifndef Edger_H
#define Edger_H

#include <vector>
#include <fstream>

#include "AttributeRanker.h"
#include "Dataset.h"

class Edger : public AttributeRanker
{
public:
  /*************************************************************************//**
   * Construct an Edger algorithm object.
   * \param [in] ds pointer to a Dataset object
   ****************************************************************************/
  Edger(Dataset* ds);
  ~Edger();
  /*************************************************************************//**
   * For each attribue, calculate Edger and associated p-value. Return
   * in a vector of pairs indexed by attribute index.
   * \return vector of pairs of Edger scores and associated p-values
   ****************************************************************************/
  std::vector<std::pair<double, std::string> > ComputeScores();
  /*************************************************************************//**
   * Print the scores to a stream.
   * \param [in] outStream reference to an output stream
   * \param [in] topN top number of attributes to print
   ****************************************************************************/
  void PrintScores(std::ofstream& outStream, unsigned int topN = 0);
  /*************************************************************************//**
   * Print the scores to a file.
   * \param [in] outFilename filename to write scores to
   * \param [in] topN top number of attributes to print
   ****************************************************************************/
  void WriteScores(std::string outFilename, unsigned int topN = 0);
private:
  Dataset* dataset;
  /// vector of pairs: scores, attribute names
  /// Edger 1.0 - pvalue for each attribute
  std::vector<std::pair<double, std::string> > scores;
  /// Write the data set in Mayo GEO format for reading into my Edger script.
  void WriteDatasetInMayoFormat(std::string filename);
  /// Read Edger scores into scores map from Edger results file.
  bool ReadEdgerScores(std::string resultsFilename);
};

#endif
