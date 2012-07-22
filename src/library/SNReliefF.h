/**
 * \class SNReliefF
 *
 * \brief Signal-to-Noise ReliefF attribute ranking algorithm.
 *
 * Designed to handle digital gene expression (DGE) data sets, particularly
 * rna-Seq high throughput count data, by accounting for variable-specific
 * variance in counts. Data is known to follow a Poisson or negative binomial
 * distribution. This algorithm is a more computationally practical approach
 * than others that use more sophisticated methods. Our approach keeps the
 * ReliefF algorithm general while addressing the variance "dispersion" problem
 * as a special case.
 *
 * \sa ReliefF
 *
 * \author Bill White
 * \version 1.0
 *
 * Contact: bill.c.white@gmail.com
 * Created on: 7/21/12
 */

#ifndef SNRELIEFF_H
#define	SNRELIEFF_H

#include <vector>
#include <map>
#include <string>
#include <fstream>

#include "ReliefF.h"
#include "Dataset.h"
#include "DatasetInstance.h"
#include "Insilico.h"
#include <boost/program_options.hpp>

/// a pair of vectors for hit and miss statistics for each instance
typedef std::vector<std::pair<double, double> > InstanceAttributeStats;
typedef std::pair<InstanceAttributeStats, InstanceAttributeStats>
	InstanceHitMissStats;

/// a map from instance ID to neighbor statistics
typedef std::vector<InstanceHitMissStats> NeighborStats;
typedef std::vector<InstanceHitMissStats>::iterator NeighborStatsIt;
typedef std::vector<InstanceHitMissStats>::const_iterator NeighborStatsCIt;

namespace po = boost::program_options;

class SNReliefF : public ReliefF
{
public:
  /*************************************************************************//**
   * Construct an ReliefF algorithm object.
   * \param [in] ds pointer to a Dataset object
   ****************************************************************************/
  SNReliefF(Dataset* ds);
  /*************************************************************************//**
   * Construct an ReliefF algorithm object.
   * \param [in] ds pointer to a Dataset object
   * \param [in] vm reference to a Boost map of command line options
   ****************************************************************************/
  SNReliefF(Dataset* ds, po::variables_map& vm);
  /*************************************************************************//**
   * Construct an ReliefF algorithm object.
   * \param [in] ds pointer to a Dataset object
   * \param [in] configMap reference to a ConfigMap (map<string, string>)
   ****************************************************************************/
  SNReliefF(Dataset* ds, ConfigMap& configMap);
  bool ComputeAttributeScores();
  /// Precompute nearest neighbor gene statistics for all instances.
  bool PreComputeNeighborGeneStats();
  virtual ~SNReliefF();
private:
  bool ComputeInstanceStats(DatasetInstance* dsi,
  		std::vector<unsigned int> hitIndicies,
  		std::vector<unsigned int> missIndicies,
  		InstanceHitMissStats& hitMissStats);
  /// nearest neighbor attribute averages and standard deviations
  NeighborStats neighborStats;
};

#endif	/* SNReliefF_H */
