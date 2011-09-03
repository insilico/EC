/* 
 * File:   EvaporativeCooling.h
 * Author: billwhite
 *
 * Created on July 14, 2011, 9:25 PM
 */

#ifndef EVAPORATIVECOOLING_H
#define	EVAPORATIVECOOLING_H

#include "../cpprelieff/Dataset.h"
#include <boost/program_options.hpp>
#include "rjungle/RJunglePar.h"

namespace po = boost::program_options;

class EvaporativeCooling
{
public:
  EvaporativeCooling(Dataset* ds, po::variables_map& vm);
  bool ComputeECScores();
  void WriteAttributeScores(std::string baseFilename);
  void PrintAttributeScores(std::ofstream& outFile);
  virtual ~EvaporativeCooling();
private:
  bool RunRandomJungle();
  Dataset* dataset;

  std::map<std::string, double> rjScores;
  std::map<std::string, double> rfScores;
  std::map<std::string, double> ecScores;
  
  RJunglePar rjParams;

  unsigned int removePerIteration;
};

#endif	/* EVAPORATIVECOOLING_H */

