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
#include "Deseq.h"
#include "Insilico.h"

namespace po = boost::program_options;

/// evaporative cooling scores - sorted by score key
typedef std::vector<std::pair<double, std::string> > EcScores;
/// evaporative cooling scores iterator - sorted by score key
typedef std::vector<std::pair<double, std::string> >::iterator EcScoresIt;
/// evaporative cooling scores constant iterator - sorted by score key
typedef std::vector<std::pair<double, std::string> >::const_iterator EcScoresCIt;

class EvaporativeCooling {
public:
	/*************************************************************************//**
	 * Construct an EC algorithm object.
	 * \param [in] ds pointer to a Dataset object
	 * \param [in] vm reference to a Boost map of command line options
	 * \param [in] anaType analysis type
	 ****************************************************************************/
	EvaporativeCooling(Dataset* ds, po::variables_map& vm,
			AnalysisType anaType =SNP_ONLY_ANALYSIS);
	/*************************************************************************//**
	 * Construct an EC algorithm object.
	 * \param [in] ds pointer to a Dataset object
	 * \param [in] configMap reference to a ConfigMap (map<string, string>)
	 * \param [in] anaType analysis type
	 ****************************************************************************/
	EvaporativeCooling(Dataset* ds, ConfigMap& configMap,
			AnalysisType anaType = SNP_ONLY_ANALYSIS);
	virtual ~EvaporativeCooling();
	/// Compute the EC scores based on the current set of attributes.
	bool ComputeECScores();
	/// Get the last computed RandomJungle scores.
	EcScores& GetRandomJungleScores();
	/// Get the last computed ReliefF scores.
	EcScores& GetReliefFScores();
	/// Get the last computed EC scores.
	EcScores& GetECScores();
	/// Return the algorithm type: EC_ALL, EC_RJ or EC_RF.
	EcAlgorithmType GetAlgorithmType();
	/*************************************************************************//**
	 * Write the scores and attribute names to file.
	 * \param [in] baseFilename filename to write score-attribute name pairs
	 ****************************************************************************/
	void WriteAttributeScores(std::string baseFilename);
	/*************************************************************************//**
	 * Write the EC scores and attribute names to stream.
	 * \param [in] outStream stream to write score-attribute name pairs
	 ****************************************************************************/
	void PrintAttributeScores(std::ofstream& outStream);
	/*************************************************************************//**
	 * Write the RJ scores and attribute names to stream.
	 * \param [in] outStream stream to write score-attribute name pairs
	 ****************************************************************************/
	void PrintRJAttributeScores(std::ofstream& outStream);
	/*************************************************************************//**
	 * Write the Deseq scores and attribute names to stream.
	 * \param [in] outStream stream to write score-attribute name pairs
	 ****************************************************************************/
	void PrintDeseqAttributeScores(std::ofstream& outStream);
	/*************************************************************************//**
	 * Write the RF scores and attribute names to stream.
	 * \param [in] outStream stream to write score-attribute name pairs
	 ****************************************************************************/
	void PrintRFAttributeScores(std::ofstream& outStream);
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
	bool RemoveWorstAttributes(unsigned int numToRemove = 1);
	/// optimize the temperature coupling constant
	double OptimizeTemperature(std::vector<double> deltas);
	/// use Random Jungle to compute the classification error of the current
	/// set of attributes with numToRemovePerIteration attributes removed
	double ComputeClassificationErrorRJ();

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
	/// pointer to a RandomJungle algorithm object
	RandomJungle* randomJungle;
	/// pointer to a DESeq algorithm object
	Deseq* deseq;

	bool optimizeTemperature;
	double optimalTemperature;
	double bestClassificationError;

	/// current random jungle scores
	EcScores rjScores;
	/// current deseq scores
	EcScores deseqScores;
	/// current relieff scores
	EcScores rfScores;
	/// current free energy scores
	EcScores freeEnergyScores;

	// number of threads to use for random jungle
	unsigned int numRFThreads;
	/// number of attributes to remove per iteration
	unsigned int numToRemovePerIteration;
	/// number of attributes to remove next iteration
	unsigned int numToRemoveNextIteration;

	/// number of target attributes
	unsigned int numTargetAttributes;
	/// attributes that have been evaporated so far
	EcScores evaporatedAttributes;
	/// current set of ec scores
	EcScores ecScores;
};

/// HACK FOR AUTOTOOLS LIBRARY DETECTION
extern "C" {
void libec_is_present(void);
}

#endif	/* EVAPORATIVECOOLING_H */
