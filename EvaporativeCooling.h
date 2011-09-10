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

#include <map>
#include <boost/program_options.hpp>
#include "rjungle/RJunglePar.h"
#include "../cpprelieff/Dataset.h"
#include "../cpprelieff/ReliefF.h"

namespace po = boost::program_options;

typedef std::map<std::string, double> EcScoresMap;
typedef std::map<std::string, double>::iterator EcScoresMapIt;
typedef std::map<std::string, double>::const_iterator EcScoresMapCIt;

class EvaporativeCooling
{
public:
  EvaporativeCooling(Dataset* ds, po::variables_map& vm);
  bool ComputeECScores();
  EcScoresMap& GetRandomJungleScores() {
    return rjScores;
  }
  EcScoresMap& GetReliefFScores() {
    return rfScores;
  }
  EcScoresMap& GetECScores() {
    return ecScores;
  }
  void WriteAttributeScores(std::string baseFilename);
  void PrintAttributeScores(std::ofstream& outFile);
  virtual ~EvaporativeCooling();
private:
  bool RunRandomJungle();
  bool ReadRandomJungleScores(std::string filename);
  bool RunReliefF();
  bool ComputeFreeEnergy(double temperature);
  bool RemoveWorstAttributes();
  
  Dataset* dataset;
  po::variables_map paramsMap;
  
  RJunglePar rjParams;

  ReliefF* reliefF;

  EcScoresMap rjScores;
  EcScoresMap rfScores;
  EcScoresMap freeEnergyScores;

  unsigned int numTargetAttributes;
  EcScoresMap evaporatedAttributes;
  EcScoresMap ecScores;
  
  // attributes to be considered
  std::map<std::string, bool> attributesToConsider;
};

#endif	/* EVAPORATIVECOOLING_H */

