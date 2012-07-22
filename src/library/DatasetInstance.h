/**
 * \class DatasetInstance
 *
 * \brief Class to hold dataset instances (rows of attributes).
 *
 * Reworked entirely for McKinney Lab work - 2/28/11
 *
 * \author Bill White
 * \version 1.0
 *
 * Contact: bill.c.white@gmail.com
 * Created on: 6/14/05
 */

#ifndef DATASET_INSTANCE_H
#define DATASET_INSTANCE_H

#include <string>
#include <vector>
#include <map>
#include <algorithm>

#include "Insilico.h"

/// forward reference to avoid circular include problems
class Dataset;

class DatasetInstance
{
public:
  /*************************************************************************//**
   * Construct an data set instance object.
   * \param [in] ds pointer to a Dataset object
   ****************************************************************************/
  DatasetInstance(Dataset* ds);
  ~DatasetInstance();
  /// return the Dataset pointer associated with this instance
  Dataset* GetDatasetPtr();
  /*************************************************************************//**
   * Load this instance with the attributes and class value from the
   * newAttributes vector.
   * \param [in] newAttributes vector of new attribute values
   * \return success
   ****************************************************************************/
  bool LoadInstanceFromVector(std::vector<AttributeLevel> newAttributes);
  /// return the number of discrete attributes
  unsigned int NumAttributes();
  /*************************************************************************//**
   * Get and return an attribute value at index.
   * \param [in] index attribute index
   * \return attribute value at index
   ****************************************************************************/
  AttributeLevel GetAttribute(unsigned int index);
  /// return the number of continuous attributes
  unsigned int NumNumerics();
  /*************************************************************************//**
   * Get and return numeric value at index.
   * \param [in] index numeric index
   * \return numeric value at index
   ****************************************************************************/
  NumericLevel GetNumeric(unsigned int index);
  /*************************************************************************//**
   * Add a numeric value to the instance's numerics vector
   * \param [in] newNum new numeric value
   * \return success
   ****************************************************************************/
  bool AddNumeric(NumericLevel newNum);
  /// Get the discrete class value.
  ClassLevel GetClass();
  /// Set the discrete class value.
  void SetClass(ClassLevel classValue);
  /// Get the continuous class value.
  double GetPredictedValueTau();
  /// Set the continuous class value.
  void SetPredictedValueTau(double newValue);
  /// Get the nearest neighbor value at neighborIndex.
  double GetInfluenceFactorD(unsigned int neighborIndex);
  /// Clear all nearest neighbor values.
  void ClearInfluenceFactors();
  /// Add the next nearest neighbor influence factor.
  bool AddInfluenceFactorD(double factor);
  /// Print the attributes, numerics and class name of this instance to stdout
  void Print();
  /*************************************************************************//**
   * Swap attribute/column values in this instance.
   * \param [in] a1 attribue index 1
   * \param [in] a2 attribue index 2
   * \return bool success
   ****************************************************************************/
  bool SwapAttributes(unsigned int a1, unsigned int a2);
  /*************************************************************************//**
   * Set the best kNearestNeighbors from the same and different classes
   * SIDE_EFFECT: Sorts and loads class the vairables: sameSums snd diffSums
   *              from the neighbors
   * \param [in] kNearestNeighbors k nearest nerighbors,
   * \param [in] sameCLassSums vectors of pairs <instance, sum> of same class
   * \param [in] diffClassSums vectors of pairs <instance, sum> of other classes
   * \return nothing
   ****************************************************************************/
  void SetDistanceSums(unsigned int kNearestNeighbors,
                       DistancePairs& sameClassSums,
                       std::map<ClassLevel,
                       DistancePairs>& diffClassSums);
  /*************************************************************************//**
   * Set the best kNearestNeighbors from all other instances/neighbors.
   * SIDE_EFFECT: Sorts and loads neighborSums from the instanceSums
   * \param [in] kNearestNeighbors k nearest neighbors
   * \param [in] instanceSums vectors of k pairs <instance, sum> for neighbors
   * \return nothing
   ****************************************************************************/
  void SetDistanceSums(unsigned int kNearestNeighbors,
                       DistancePairs instancesSums);
  /*************************************************************************//**
   * Prints passed distance pairs.
   * \param [in] distPairs distance pairs
   ****************************************************************************/
  void PrintDistancePairs(const DistancePairs& distPairs);
  /*************************************************************************//**
   * Returns N closest instances using the sameSums and diffSums class variables
   * \param [in] n n nearest nerighbors
   * \param [in] sameCLassInstances vector of same class instances indices
   * \param [in] diffClassInstances vector of different class instance indices
   * \return success
   ****************************************************************************/
  bool GetNNearestInstances(unsigned int n,
                            std::vector<unsigned int>& sameClassInstances,
                            std::vector<unsigned int>& diffClassInstances);
  /*************************************************************************//**
   * Returns N closest instances using the sameSums and diffSums class variables
   * \param [in] n n nearest nerighbors
   * \param [in] sameCLassInstances vector of same class instances indices
   * \param [in] diffClassInstances vector of different classes instance indices
   * \return success
   ****************************************************************************/
  bool GetNNearestInstances(unsigned int n,
                            std::vector<unsigned int>& sameClassInstances,
                            std::map<ClassLevel, std::vector<unsigned int> >& diffClassInstances);
  /*************************************************************************//**
   * Returns N closest instances to this instance.
   * \param [in] n n nearest neighbors
   * \param [in] closestInstances reference to a vector of instance indices
   * \return success
   ****************************************************************************/
  bool GetNNearestInstances(unsigned int n,
                            std::vector<unsigned int>& closestInstances);
  /// discrete attributes
  std::vector<AttributeLevel> attributes;
  /// continuous attributes
  std::vector<NumericLevel> numerics;
private:
  /// pointer to a Dataset object
  Dataset* dataset;
  /// the class value for this instance
  ClassLevel classLabel;
  /// vector of instance IDs for the best neighbors in this instance's class
  std::vector<std::string> bestNeighborIdsSameClass;
  /// vector of instance IDs for the best neighbors of different class(es)
  std::map<ClassLevel, std::vector<std::string> > bestNeighborIdsDiffClass;
  /// best neighbor IDs for continuous class
  std::vector<std::string> bestNeighborIds;
  /// nearest neighbor weighting factors
  std::vector<double> neighborInfluenceFactorDs;
  /// continuous value for this class
  double predictedValueTau;
};

#endif
