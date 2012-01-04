/**
 * \class EvaporativeCooling
 *
 * \brief Evaporative Cooling attribute ranking algorithm.
 *
 * Implements the Evaporative Cooling algorithm in:
 * McKinney, et. al. "Capturing the Spectrum of Interaction Effects in Genetic
 * Association Studies by Simulated Evaporative Cooling Network Analysis."
 * PLoS Genetics, Vol 5, Issue 3, 2009.
 *
 * \sa ReliefF
 * \sa RReliefF
 * \sa RandomJungle
 *
 * \author Bill White
 * \version 1.0
 *
 * Contact: bill.c.white@gmail.com
 * Created on: 7/14/11
 */

#ifndef EVAPORATIVECOOLING_H
#define	EVAPORATIVECOOLING_H

#include <vector>
#include <boost/program_options.hpp>

#include "Dataset.h"
#include "RandomJungle.h"
#include "ReliefF.h"

namespace po = boost::program_options;

/// evaporative cooling scores - sorted by score key
typedef std::vector<std::pair<double, std::string> > EcScores;
/// evaporative cooling scores iterator - sorted by score key
typedef std::vector<std::pair<double, std::string> >::iterator EcScoresIt;
/// evaporative cooling scores constant iterator - sorted by score key
typedef std::vector<std::pair<double, std::string> >::const_iterator EcScoresCIt;

/**
 * \enum EcAlgorithmType.
 * Type of algorithm steps to perform.
 */
typedef enum
{
  EC_ALL, /**< Run RandomJungle and ReliefF */
  EC_RJ,  /**< Run only RandomJungle */
  EC_RF   /**< Run only ReliefF */
} EcAlgorithmType;

class EvaporativeCooling
{
public:
  /*************************************************************************//**
   * Construct an EC algorithm object.
   * \param [in] ds pointer to a Dataset object
   * \param [in] vm reference to a Boost map of command line options
   * \param [in] anaType analysis type
   ****************************************************************************/
  EvaporativeCooling(Dataset* ds, po::variables_map& vm,
                     AnalysisType anaType=SNP_ONLY_ANALYSIS);
  virtual ~EvaporativeCooling();
  /// Compute the EC scores based on the current set of attributes.
  bool ComputeECScores();
  /// Get the last computed RandomJungle scores.
  EcScores& GetRandomJungleScores();
  /// Get the last computed ReliefF scores.
  EcScores& GetReliefFScores();
  /// Get the last computed EC scores.
  EcScores& GetECScores();
  /*************************************************************************//**
   * Write the scores and attribute names to file.
   * \param [in] baseFilename filename to write score-attribute name pairs
   ****************************************************************************/
  void WriteAttributeScores(std::string baseFilename);
  /*************************************************************************//**
   * Write the scores and attribute names to stream.
   * \param [in] outStream stream to write score-attribute name pairs
   ****************************************************************************/
  void PrintAttributeScores(std::ofstream& outStream);
  /// Print the current attributes scores to stdout in tab-delimited format.
  bool PrintAllScoresTabular();
  /// Print the kendall taus between the ReliefF and RandomJungle scores.
  bool PrintKendallTaus();
private:
  /// Run the ReliefF algorithm.
  bool RunReliefF();
  /*************************************************************************//**
   * Compute the attributes' free energy using the couple temperature.
   * \param [in] tempreatire coupling temperature T
   * \return distance
   ****************************************************************************/
  bool ComputeFreeEnergy(double temperature);
  /*************************************************************************//**
   * Remove the worst attribute based on free energy scores.
   * \param [in] numToRemove number of attributes to remove/evaporate
   * \return distance
   ****************************************************************************/
  bool RemoveWorstAttributes(unsigned int numToRemove=1);
  
  /// pointer to a Dataset object
  Dataset* dataset;
  /// command line parameters map
  po::variables_map paramsMap;
  /// prefix for all output files
  std::string outFilesPrefix;

  /// type of analysis to perform
  /// \sa ReliefF
  AnalysisType analysisType;
  /// algorithm steps to perform
  EcAlgorithmType algorithmType;

  /// pointer to a ReliefF or RReliefF algorithm object
  ReliefF* reliefF;
  /// pointer to a RandomJungle algorithm onject
  RandomJungle* randomJungle;

  /// current random jungle scores
  EcScores rjScores;
  /// current relieff scores
  EcScores rfScores;
  /// current free energy scores
  EcScores freeEnergyScores;

  // number of threads to use for random jungle
  unsigned int numRFThreads;
  /// number of attributes to remove per iteration
  unsigned int numToRemovePerIteration;
  /// number of target attributes
  unsigned int numTargetAttributes;
  /// attributes that have been evaporated so far
  EcScores evaporatedAttributes;
  /// current set of ec scores
  EcScores ecScores;
};

#endif	/* EVAPORATIVECOOLING_H */
