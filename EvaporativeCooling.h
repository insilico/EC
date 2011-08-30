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

namespace po = boost::program_options;

class EvaporativeCooling
{
public:
  EvaporativeCooling(Dataset* ds, po::variables_map& vm);
  bool GetECScores(std::map<std::string, double>& scores);
  virtual ~EvaporativeCooling();
private:
  Dataset* dataset;
  unsigned int removePerIteration;
};

#endif	/* EVAPORATIVECOOLING_H */

