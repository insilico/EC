/**
 * \class RandomJungle
 *
 * \brief RandomJungle attribute ranking algorithm.
 *
 * Adapter class to map EC call for Random Jungle importance scores
 * to Random Jungle library functions.
 *
 * \author Bill White
 * \version 1.0
 *
 * Contact: bill.c.white@gmail.com
 * Created on: 10/16/11
 */

#ifndef RANDOMJUNGLE_H
#define	RANDOMJUNGLE_H

#include <vector>
#include <fstream>

#include "Dataset.h"
#include <boost/program_options.hpp>
#include "rjungle/RJunglePar.h"

namespace po = boost::program_options;

class RandomJungle
{
public:
  /*************************************************************************//**
   * Construct an RandomJungle algorithm object.
   * \param [in] ds pointer to a Dataset object
   * \param [in] vm reference to a Boost map of command line options
   ****************************************************************************/
  RandomJungle(Dataset* ds, po::variables_map& vm);
  virtual ~RandomJungle();
  bool ComputeAttributeScores();
  /*************************************************************************//**
   * Get the (importance) scores as a vector of pairs: score, attribute name
   * \return vector of pairs
   ****************************************************************************/
  std::vector<std::pair<double, std::string> > GetScores();
private:
  /// Read the importance scores as attribute rankings from file.
  bool ReadScores(std::string importanceFilename);
  /// RandomJungle parameters object
  RJunglePar rjParams;
  /// pointer to a Dataset object
  Dataset* dataset;
  /// vector of pairs: scores, attribute names
  std::vector<std::pair<double, std::string> > scores;
};

#endif	/* RANDOMJUNGLE_H */

