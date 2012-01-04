/**
 * \file Statstics.h
 *
 * \brief Statistical utilities.
 *
 * \author Bill White
 * \version 1.0
 *
 * Contact: bill.c.white@gmail.com
 * Created on: 2/13/06
 */

#ifndef STATISTICS_H
#define STATISTICS_H

#include <vector>
#include <map>
#include <numeric>
#include <iterator>
#include <cmath>
#include <algorithm>

#include "Dataset.h"

/// vector of doubles type
typedef std::vector<double> VectorDouble;
/// vector of doubles iterator
typedef std::vector<double>::const_iterator VectorDoubleIt;
/// histogram type is a map: value->count
typedef std::map<unsigned int, unsigned int> Histogram;
/// historgram iterator
typedef std::map<unsigned int, unsigned int>::const_iterator HistogramIt;

/***************************************************************************//**
 * ZTransform input values.
 * \param [in] inputValues const vector of double input values
 * \param [out] outputValues transformed input values to z-scores with mean=0, stddev=1
 * \return success
 ******************************************************************************/
bool ZTransform(const VectorDouble& inputValues, VectorDouble& outputValues);
/***************************************************************************//**
 * Calculates the entropy of a sequence of unsigned integers.
 * \param [in] attributeValues vector of sequence values - unsigned ints - positive categorical
 * \return entropy as a double-precision float
 ******************************************************************************/
double Entropy(const std::vector<AttributeLevel>& attributeValues);
/***************************************************************************//**
 * Calculates the conditional entropy of a sequence of unsigned integers
 * based (conditioned) on another sequence of unsigned integers (the givens).
 * P(sequenceValues | givenValues)
 * \param [in] attributeValues vector of values
 * \param [in] givenValues vector of givens
 * \return conditiional entropy as a double-precision float
 ******************************************************************************/
double ConditionalEntropy(const std::vector<AttributeLevel>& attributeValues,
                          const std::vector<AttributeLevel>& givenValues);
/***************************************************************************//**
 * Create a new attribute that is the cartesian product of a and b.
 * NOTE: works for genotypes; need to verify for missign data levels, etc.
 * \param [in] a attributes vector a
 * \param [in] b attributes vector b
 * \param [out] vector ab, the cartesian product of a and b
 * \return success
 ******************************************************************************/
bool ConstructAttributeCart(const std::vector<AttributeLevel>& a, 
                            const std::vector<AttributeLevel>& b,
                            std::vector<AttributeLevel>& ab);
/***************************************************************************//**
 * Compute KendallTau for two ranked vectors of strings.
 * Why Kenall Tau - G. E. NOETHER
 * http://www.rsscse-edu.org.uk/tsj/bts/noether/text.html
 * \param [in] X ranked attribute vector X
 * \param [in] Y ranked attribute vector Y
 * \return Kendall Tau value (-1, 1)
 ******************************************************************************/
double KendallTau(std::vector<std::string> X, std::vector<std::string> Y);
/***************************************************************************//**
 * Compute KendallTau for two ranked vectors of doubles.
 * Why Kenall Tau - G. E. NOETHER
 * http://www.rsscse-edu.org.uk/tsj/bts/noether/text.html
 * \param [in] X ranked attribute vector X
 * \param [in] Y ranked attribute vector Y
 * \return Kendall Tau value (-1, 1)
 ******************************************************************************/
double KendallTau(std::vector<double> X, std::vector<double> Y);
/***************************************************************************//**
 * Compute KendallTau for two ranked vectors of integers.
 * Why Kenall Tau - G. E. NOETHER
 * http://www.rsscse-edu.org.uk/tsj/bts/noether/text.html
 * \param [in] X ranked attribute vector X
 * \param [in] Y ranked attribute vector Y
 * \return Kendall Tau value (-1, 1)
 ******************************************************************************/
double KendallTau(std::vector<int> X, std::vector<int> Y);
/***************************************************************************//**
 * Calculate variance and standard deviation of a vector of values.
 * \param [in]  ranked attribute lists X and Y
 * \return Kendall Tau value (-1, 1)
 ******************************************************************************/
template<class T>
std::pair<double, double> VarStd(std::vector<T>& values) {
  double average = (double) std::accumulate(values.begin(), values.end(), 0) /
          values.size();

  // calculate and return variance
  double SSE = 0.0;
  typename std::vector<T>::const_iterator it;
  for(it = values.begin(); it != values.end(); it++) {
    SSE += ((*it - average) * (*it - average));
  }
  double var = SSE / (values.size() - 1);
  double std = sqrt(var);

  return std::make_pair(var, std);
}

#endif
