/*
  Statistics.C - Bill White - 2/13/06

  Statistical routines: implementations.
*/

#include <iostream>
#include <vector>
#include <cmath>

#include "Dataset.h"
#include "DatasetInstance.h"
#include "Statistics.h"

using namespace std;

#define DEBUG_Z 0
#define DEBUG_E 0

bool ZTransform(const VectorDouble& inputValues, 
								VectorDouble& outputValues)
{
	// sum values in input vector
	double inputSum = 0.0;
	double inputMin = 0, inputMax = 0;
	VectorDoubleIt sumIt;
	for(sumIt = inputValues.begin(); sumIt != inputValues.end(); sumIt++) {
		if(sumIt == inputValues.begin()) {
			inputMin = *sumIt;
			inputMax = *sumIt;
		}
		else {
			if(*sumIt < inputMin) { inputMin = *sumIt; }
			if(*sumIt > inputMax) { inputMax = *sumIt; }
		}
		inputSum += (*sumIt);
	}

	// calculate the mean
	double inputMean = inputSum / inputValues.size();

	// get the sum-squared deviations from the mean
	double sumSqDev = 0.0;
	VectorDoubleIt sumsqIt;
	for(sumsqIt = inputValues.begin(); sumsqIt != inputValues.end(); sumsqIt++) {
		double inputX = *sumsqIt;
		sumSqDev += ((inputX - inputMean) * (inputX - inputMean));
	}

	// compute variance and standard deviation 
	double inputVar = sumSqDev / (inputValues.size()-1);
	double inputStd = sqrt(inputVar);

#if DEBUG_Z
	cout << "z-transform stats: " << endl;
	cout << "\tN:        " << inputValues.size() << endl;
	cout << "\tMinimum:  " << inputMin << endl;
	cout << "\tMaximum:  " << inputMax << endl;
	cout << "\tSum:      " << inputSum << endl;
	cout << "\tMean:     " << inputMean << endl;
	cout << "\tSumSqDev: " << sumSqDev << endl;
	cout << "\tVar:      " << inputVar << endl;
	cout << "\tStdDev:   " << inputStd << endl;
#endif
	
	// compute and return the z-transformed input values, ie z-scores
	outputValues.clear();
	VectorDoubleIt ztransIt;
	for(ztransIt = inputValues.begin(); ztransIt != inputValues.end(); ztransIt++) {
		double curZScore = ((*ztransIt) - inputMean) / inputStd;
#if DEBUG_Z
		cout << curZScore << endl;
#endif
		outputValues.push_back(curZScore);
	}

	return true;
}

double Entropy(const vector<AttributeLevel>& sequenceValues)
{
	// get the counts for each level in the sequence 
	// (histogram of sequence values)
	Histogram countsByLevel;
	for(vector<AttributeLevel>::const_iterator aItKey=sequenceValues.begin();
			aItKey != sequenceValues.end(); 
			aItKey++)	{
		++countsByLevel[*aItKey];
	}
	
	// calculate the entropy using prior probabilities (p)
	double entropy = 0.0;
	unsigned int sequenceSize = sequenceValues.size();
	for(HistogramIt cblItKey=countsByLevel.begin(); 
			cblItKey != countsByLevel.end(); 
			cblItKey++)	{
		unsigned int thisCount = (*cblItKey).second;
		double p = ((double) thisCount) / ((double) sequenceSize);
		entropy -= p * log(p) / log(2.0);
	}

	return entropy;
}

double ConditionalEntropy(const vector<AttributeLevel>& sequenceValues,
													const vector<AttributeLevel>& givenValues)
{
	// check sequences for proper format
	unsigned int sequenceSize = sequenceValues.size();
	unsigned int givensSize   = givenValues.size();
	if(sequenceSize != givensSize) {
		cerr << "Statistics::ConditionalEntropy Failed: sequence sizes are not equal: " 
				 << sequenceSize << " vs. " << givensSize << endl;
		exit(-1); // this is rather drastic!
	}

	// get a map of given levels to sequence level counts, ie P(sequence | givens)
	map<unsigned int, Histogram*> counts;
	Histogram* subCounts;
	for(unsigned int sequenceIdx=0; sequenceIdx < sequenceSize; sequenceIdx++)	{
		unsigned int curSequenceValue = sequenceValues[sequenceIdx];
		unsigned int curGivenValue = givenValues[sequenceIdx];
		if(counts.find(curGivenValue) != counts.end()) {
			subCounts = counts[curGivenValue];
			if(subCounts->find(curSequenceValue) != subCounts->end()) {
				(*subCounts)[curSequenceValue]++;
			}
			else {
				(*subCounts)[curSequenceValue] = 1;
			}
		}
		else {
			subCounts = new Histogram();
			(*subCounts)[curSequenceValue] = 1;
			counts[curGivenValue] = subCounts;
		}
	}

#if DEBUG_E
	cout << "Counts matrix:" << endl << endl;
	map<unsigned int, Histogram*>::const_iterator printCountsIt;
	for(printCountsIt=counts.begin(); printCountsIt != counts.end(); printCountsIt++) {
		cout << (*printCountsIt).first << ":" << endl;
		Histogram* levels =(*printCountsIt).second;
		for(HistogramIt hIt = levels->begin(); hIt != levels->end(); hIt++) {
			cout << "\t" << (*hIt).first << "\t" << (*hIt).second << endl;
		}
	}
	cout << endl;
#endif

	// calculate the entropy using probabilities from frequency counts
	double entropy = 0.0;
	double subEntropy = 0.0;

	map<unsigned int, Histogram*>::const_iterator countsIt;
	for(countsIt=counts.begin(); countsIt != counts.end(); countsIt++) {
		Histogram* subCounts = (*countsIt).second;
		HistogramIt subCountsIt;
		unsigned int subTotal = 0;
		for(subCountsIt = subCounts->begin(); subCountsIt != subCounts->end(); subCountsIt++) {
			subTotal += (*subCountsIt).second;
		}
		for(subCountsIt = subCounts->begin(); subCountsIt != subCounts->end(); subCountsIt++) {
			double p = (double) (*subCountsIt).second / (double) subTotal;;
			subEntropy -= p * log(p) / log(2.0);
		}
		entropy += subTotal * subEntropy / givensSize;
	}

	return entropy;
}

bool ConstructAttributeCart(const vector<AttributeLevel>& a,
                            const vector<AttributeLevel>& b,
														vector<AttributeLevel>& ab)
{
	ab.clear();
	vector<AttributeLevel>::const_iterator aIt;
	vector<AttributeLevel>::const_iterator bIt;
	for(aIt = a.begin(),bIt  = b.begin(); aIt != a.end(); aIt++, bIt++) {
		ab.push_back((*aIt) * 3 + (*bIt));
	}
	return true;
}

double KendallTau(vector<string> X, vector<string> Y)
{
  unsigned int numX = X.size();
  unsigned int numY = Y.size();
  if(numX != numY) {
    cerr << "ERROR: KendallTau: lists must be the same size" << endl;
    return -2.0;
  }
  double n = (double) numX;
  double C = 0;
  double D = 0;
  double S = 0;
  for(unsigned int i=0; i < n; ++i) {
    for(unsigned int j=i+1; j < n; ++j) {
      string Xi = X[i];
      string Xj = X[j];
      string Yi = Y[i];
      string Yj = Y[j];
      if((Xi == Yi) && (Xj == Yj)) {
        ++C;
      }
      else {
        ++D;
      }
    }
  }
  S = C - D;
  double tau = 2.0 * S / (n * (n - 1));

  return tau;
}

double KendallTau(vector<double> X, vector<double> Y)
{
  unsigned int numX = X.size();
  unsigned int numY = Y.size();
  if(numX != numY) {
    cerr << "ERROR: KendallTau: lists must be the same size" << endl;
    return -2.0;
  }
  double n = (double) numX;
  double C = 0;
  double D = 0;
  double S = 0;
  for(unsigned int i=0; i < n; ++i) {
    for(unsigned int j=i+1; j < n; ++j) {
      double Xi = X[i];
      double Xj = X[j];
      double Yi = Y[i];
      double Yj = Y[j];
      double xDiff = Xj - Xi;
      double yDiff = Yj - Yi;
      if( xDiff * yDiff > 0) {
        ++C;
      }
      else {
        ++D;
      }
    }
  }
  S = C - D;
  double tau = 2.0 * S / (n * (n - 1));

  return tau;
}

double KendallTau(vector<int> X, vector<int> Y)
{
  unsigned int numX = X.size();
  unsigned int numY = Y.size();
  if(numX != numY) {
    cerr << "ERROR: KendallTau: lists must be the same size" << endl;
    return -2.0;
  }
  double n = (double) numX;
  double C = 0;
  double D = 0;
  double S = 0;
  for(unsigned int i=0; i < n; ++i) {
    for(unsigned int j=i+1; j < n; ++j) {
      double Xi = (double) X[i];
      double Xj = (double) X[j];
      double Yi = (double) Y[i];
      double Yj = (double) Y[j];
      double xDiff = Xj - Xi;
      double yDiff = Yj - Yi;
      if( xDiff * yDiff > 0) {
        ++C;
      }
      else {
        ++D;
      }
    }
  }
  S = C - D;
  double tau = 2.0 * S / (n * (n - 1));

  return tau;
}
