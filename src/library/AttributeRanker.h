/**
 * \class AttributeRanker
 *
 * \brief Abstract base class for attribute rankers.
 *
 * \author Bill White
 * \version 1.0
 *
 * Contact: bill.c.white@gmail.com
 * Created on: 8/13/12
 */

#ifndef ATTRIBUTERANKER_H
#define ATTRIBUTERANKER_H

#include <string>
#include <fstream>

#include "Dataset.h"
#include "Insilico.h"

class AttributeRanker
{
public:
  /// Construct a default data set.
  AttributeRanker(Dataset* ds);
  /// Destruct all dynamically allocated memory.
  virtual ~AttributeRanker();

  /// Compute the attribute scores for the current set of attributes.
  virtual AttributeScores ComputeScores() = 0;
  /*************************************************************************//**
   * Get the (importance) scores as a vector of pairs: score, attribute name
   * \return vector of pairs
   ****************************************************************************/
  virtual AttributeScores GetScores();
	/*************************************************************************//**
	 * Write the scores and attribute names to file.
	 * \param [in] baseFilename filename to write score-attribute name pairs
	 ****************************************************************************/
	virtual void WriteScores(std::string baseFilename);
	/*************************************************************************//**
	 * Write the EC scores and attribute names to stream.
	 * \param [in] outStream stream to write score-attribute name pairs
	 ****************************************************************************/
	virtual void PrintScores(std::ofstream& outStream);
  /// Error from using ranked attributes in a classifier.
	virtual double GetClassificationError();
protected:
  /// The Dataset on which the ranking algorithm is working.
  Dataset* dataset;
  /// attribute scores and names
  AttributeScores scores;
  /// Error from using ranked attributes in a classifier.
  double classificationAccuracy;
};

#endif
