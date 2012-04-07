/**
 * \class ReliefF
 *
 * \brief ReliefF attribute ranking algorithm.
 *
 * Totally redone for the McKinney insilico lab in 2011.
 * Large refactoring to move all attribute elimination handling to the
 * Dataset and its subclasses. 9/11/11
 *
 * \sa RReliefF
 *
 * \author Bill White
 * \version 1.0
 *
 * Contact: bill.c.white@gmail.com
 * Created on: 7/16/05
 */

#ifndef RELIEFF_H
#define RELIEFF_H

#include <vector>
#include <fstream>

#include <boost/program_options.hpp>

#include "Dataset.h"
#include "Insilico.h"

namespace po = boost::program_options;

class ReliefF
{
public:
  /*************************************************************************//**
   * Construct an ReliefF algorithm object.
   * \param [in] ds pointer to a Dataset object
   * \param [in] anaType analysis type
   ****************************************************************************/
  ReliefF(Dataset* ds, AnalysisType anaType);
  /*************************************************************************//**
   * Construct an ReliefF algorithm object.
   * \param [in] ds pointer to a Dataset object
   * \param [in] vm reference to a Boost map of command line options
   * \param [in] anaType analysis type
   ****************************************************************************/
  ReliefF(Dataset* ds, po::variables_map& vm, AnalysisType anaType);
  /*************************************************************************//**
   * Construct an ReliefF algorithm object.
   * \param [in] ds pointer to a Dataset object
   * \param [in] configMap reference to a ConfigMap (map<string, string>)
   * \param [in] anaType analysis type
   ****************************************************************************/
  ReliefF(Dataset* ds, ConfigMap& vm, AnalysisType anaType);
  virtual ~ReliefF();
  /**
   * Compute the ReliefF scores for the current set of attributes.
   * Implements ReliefF algorithm:
   * Marko Robnik-Sikonja, Igor Kononenko: Theoretical and Empirical Analysis of
   * ReliefF and RReliefF. Machine Learning Journal, 53:23-69, 2003
   * http://lkm.fri.uni-lj.si/rmarko/papers/robnik03-mlj.pdf
   */
  virtual bool ComputeAttributeScores();
  /// Compute the ReliefF scores by iteratively removing worst attributes.
  bool ComputeAttributeScoresIteratively();
  /// Resets some data structures for the next iteration of ReliefF
  bool ResetForNextIteration();
  /*************************************************************************//**
   * Write the scores and attribute names to stream.
   * \param [in] outStream stream to write score-attribute name pairs
   ****************************************************************************/
  void PrintAttributeScores(std::ofstream& outStream);
  /*************************************************************************//**
   * Write the scores and attribute names to file.
   * \param [in] baseFIlename filename to write score-attribute name pairs
   ****************************************************************************/
  void WriteAttributeScores(std::string baseFilename);
  /// Precompute all pairwise instance-to-instance distances.
  bool PreComputeDistances();
  /// Precompute all pairwise distances homoring excluded instances.
  bool PreComputeDistancesByMap();
  /// Get the last computed ReliefF scores.
  std::vector<std::pair<double, std::string> > GetScores();
protected:
  /// Compute the weight by distance factors for nearest neighbors.
  bool ComputeWeightByDistanceFactors();
  /// type of analysis to perform
  AnalysisType analysisType;
  /*************************************************************************//**
   * Compute the discrete difference in an attribute between two instances.
   * \param [in] attributeIndex index into vector of all attributes
   * \param [in] dsi1 pointer to DatasetInstance 1
   * \param [in] dsi2 pointer to DatasetInstance 2
   * \return diff(erence)
   ****************************************************************************/
  double (*snpDiff)(unsigned int attributeIndex,
                    DatasetInstance* dsi1,
                    DatasetInstance* dsi2);
  /*************************************************************************//**
   * Compute the continuous difference in an attribute between two instances.
   * \param [in] attributeIndex index into vector of all attributes
   * \param [in] dsi1 pointer to DatasetInstance 1
   * \param [in] dsi2 pointer to DatasetInstance 2
   * \return diff(erence)
   ****************************************************************************/
  double (*numDiff)(unsigned int attributeIndex,
                    DatasetInstance* dsi1,
                    DatasetInstance* dsi2);
  /// the name of discrete diff(erence) function
  std::string snpMetric;
  /// the name of continuous diff(erence) function
  std::string numMetric;
  /// the dataset on which the algorithm is working
  Dataset* dataset;
  /// nomalizing factor for ReliefF m * k loop
  double one_over_m_times_k;
  /// number of instances to sample
  unsigned int m;
  /// are instances being randomly selected?
  bool randomlySelect;
  /// k nearest neighbors
  unsigned int k;
  /// number of attributes to remove each iteration if running iteratively
  unsigned int removePerIteration;
  /// are we removing a percentage per iteration?
  bool doRemovePercent;
  /// percentage of attributes to remove per iteration if running iteratively
  double removePercentage;
  /// name of the weight-by-distance method
  std::string weightByDistanceMethod;
  /// sigma value used in exponential decay weight-by-distance
  double weightByDistanceSigma;

  /// attribute scores/weights
  std::vector<double> W;
  /// attribute names associated with scores
  std::vector<std::string> scoreNames;
  /// final scores after all iterations
  std::map<std::string, double> finalScores;
};

#endif
