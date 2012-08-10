/**
 * \class Deseq
 *
 * \brief DESeq attribute ranking algorithm.
 *
 * Deseq algorithm interface.
 **
 * \author Bill White
 * \version 1.0
 *
 * Contact: bill.c.white@gmail.com
 * Created on: 8/8/12
 */

#ifndef DESEQ_H
#define DESEQ_H

#include <vector>
#include <fstream>

#include "Dataset.h"

class Deseq
{
public:
  /*************************************************************************//**
   * Construct an deseq algorithm object.
   * \param [in] ds pointer to a Dataset object
   ****************************************************************************/
  Deseq(Dataset* ds);
  ~Deseq();
  /*************************************************************************//**
   * For each attribue, calculate deseq and associated p-value. Return
   * in a vector of pairs indexed by attribute index.
   * \return vector of pairs of deseq scores and associated p-values
   ****************************************************************************/
  const std::vector<std::pair<double, std::string> >& ComputeScores();
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
  /// deseq 1.0 - pvalue for each attribute
  std::vector<std::pair<double, std::string> > scores;
  /// Write the data set in Mayo GEO format for reading into my DESeq script.
  void WriteDatasetInMayoFormat(std::string filename);
  /// Read DESeq scores into scores map from DESeq results file.
  bool ReadDeseqScores(std::string resultsFilename);
};

#endif
