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

#include "Insilico.h"
#include "Dataset.h"
#include <boost/program_options.hpp>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "rjungle/RJunglePar.h"

namespace po = boost::program_options;

class RandomJungle
{
public:
	/// Run Random jungle as a classifier without instantiating a Random Jungle
  static bool RunClassifier(std::string csvFile, ConfigMap& vm,
  		RandomJungleTreeType treeType, double& classError);
  /// Read the classification error from file into variable
  static bool ReadClassificationError(std::string confusionFilename,
  		RandomJungleTreeType treeType, double& classifierError);

  /*************************************************************************//**
   * Construct an RandomJungle algorithm object.
   * \param [in] ds pointer to a Dataset object
   * \param [in] vm reference to a Boost map of command line options
   ****************************************************************************/
  RandomJungle(Dataset* ds, po::variables_map& vm);
  /*************************************************************************//**
   * Construct an RandomJungle algorithm object.
   * \param [in] ds pointer to a Dataset object
   * \param [in] configMap reference ConfigMap (map<string, string>)
   ****************************************************************************/
  RandomJungle(Dataset* ds, ConfigMap& vm);
  virtual ~RandomJungle();
  /// Score attributes by getting Random Jungle importance scores
  bool ComputeAttributeScores();
  /// Score attributes by getting Random Jungle importance scores
  bool ComputeAttributeScoresRjungle();
  /*************************************************************************//**
   * Get the (importance) scores as a vector of pairs: score, attribute name
   * \return vector of pairs
   ****************************************************************************/
  std::vector<std::pair<double, std::string> > GetScores();
  /// Get the classification error of the last classifier run
  double GetClassificationError();
private:
  /// Read the importance scores as attribute rankings from file into member
  /// vector scores: pair<double, string>
  bool ReadScores(std::string importanceFilename);
  /// Read classification error from file into member variable classificationError
  bool ReadClassificationError(std::string confusionFilename);

  /// RandomJungle parameters object
  RJunglePar rjParams;
  /// RandomJungle calling style
  RandomJungleRunMode runMode;
  /// pointer to a Dataset object
  Dataset* dataset;
  /// vector of pairs: scores, attribute names
  std::vector<std::pair<double, std::string> > scores;
  /// last classification error
  double classificationError;
};

#endif	/* RANDOMJUNGLE_H */

