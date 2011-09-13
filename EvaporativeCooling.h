/* 
 * File:   EvaporativeCooling.h
 * Author: billwhite
 *
 * Created on July 14, 2011, 9:25 PM
 *
 * Implements the Evaporative Cooling algorithm in:
 * McKinney, et. al. "Capturing the Spectrum of Interaction Effects in Genetic
 * Association Studies by Simulated Evaporative Cooling Network Analysis."
 * PLoS Genetics, Vol 5, Issue 3, 2009.
 */

#ifndef EVAPORATIVECOOLING_H
#define	EVAPORATIVECOOLING_H

#include <vector>
#include <boost/program_options.hpp>
#include "rjungle/RJunglePar.h"
#include "../cpprelieff/Dataset.h"
#include "../cpprelieff/ReliefF.h"

namespace po = boost::program_options;

typedef std::vector<std::pair<double, std::string> > EcScores;
typedef std::vector<std::pair<double, std::string> >::iterator EcScoresIt;
typedef std::vector<std::pair<double, std::string> >::const_iterator EcScoresCIt;

class EvaporativeCooling
{
public:
  EvaporativeCooling(Dataset* ds, po::variables_map& vm);
  virtual ~EvaporativeCooling();
  bool ComputeECScores();
  EcScores& GetRandomJungleScores();
  EcScores& GetReliefFScores();
  EcScores& GetECScores();
  void WriteAttributeScores(std::string baseFilename);
  void PrintAttributeScores(std::ofstream& outFile);
  bool PrintAllScoresTabular();
private:
  bool InitializeRandomJungle(uli_t ntree=100);
  bool RunRandomJungle();
  bool ReadRandomJungleScores(std::string filename);
  bool FinalizeRandomJungle();
  bool RunReliefF();
  bool ComputeFreeEnergy(double temperature);
  bool RemoveWorstAttributes(unsigned int numToRemove=1);
  
  Dataset* dataset;
  po::variables_map paramsMap;
  
  RJunglePar rjParams;

  ReliefF* reliefF;
  unsigned int rfNumToRemovePerIteration;

  EcScores rjScores;
  EcScores rfScores;
  EcScores freeEnergyScores;

  unsigned int numToRemovePerIteration;
  unsigned int numTargetAttributes;
  EcScores evaporatedAttributes;
  EcScores ecScores;
};

#endif	/* EVAPORATIVECOOLING_H */
