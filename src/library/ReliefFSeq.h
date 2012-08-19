/**
 * \class ReliefFSeq
 *
 * \brief ReliefFSeq attribute ranking algorithm.
 *
 * Designed to handle digital gene expression (DGE) data sets, particularly
 * RNA-Seq high-throughput count data, by accounting for variable-specific
 * variance in counts. Data is known to follow a Poisson or negative binomial
 * distribution. This algorithm is a more computationally practical approach
 * than others that use more sophisticated statistical methods and models.
 * Our approach keeps the ReliefF algorithm general while addressing the
 * variance "dispersion" problem as a special case.
 *
 * \sa ReliefF
 *
 * \author Bill White
 * \version 1.0
 *
 * Contact: bill.c.white@gmail.com
 * Created on: 7/23/12
 */

#ifndef RELIEFFSEQ_H
#define	RELIEFFSEQ_H

#include <vector>
#include <map>
#include <string>
#include <fstream>

#include "ReliefF.h"
#include "Dataset.h"
#include "DatasetInstance.h"
#include "Insilico.h"
#include <boost/program_options.hpp>

namespace po = boost::program_options;

class ReliefFSeq : public ReliefF
{
public:
  /*************************************************************************//**
   * Construct an ReliefFSeq algorithm object.
   * \param [in] ds pointer to a Dataset object
   ****************************************************************************/
  ReliefFSeq(Dataset* ds);
  /*************************************************************************//**
   * Construct an ReliefFSeq algorithm object.
   * \param [in] ds pointer to a Dataset object
   * \param [in] vm reference to a Boost map of command line options
   ****************************************************************************/
  ReliefFSeq(Dataset* ds, po::variables_map& vm);
  /*************************************************************************//**
   * Construct an ReliefFSeq algorithm object.
   * \param [in] ds pointer to a Dataset object
   * \param [in] configMap reference to a ConfigMap (map<string, string>)
   ****************************************************************************/
  ReliefFSeq(Dataset* ds, ConfigMap& configMap);
  bool ComputeAttributeScores();
  // average hit and miss diffs for gene alpha
  std::pair<double, double> MuDeltaAlphas(unsigned int alpha);
  /// standard deviations of hit and miss diffs for gene alpha
  std::pair<double, double> SigmaDeltaAlphas(unsigned int alpha,
  		double muDeltaHit, double muDeltaMiss);
  virtual ~ReliefFSeq();
private:
  std::string mode;
  double s0;
};

#endif	/* SNReliefF_H */
