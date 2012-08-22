/*
 * File:   DistanceMetrics.cpp
 * Author: billwhite
 *
 * Created on March 29, 2011, 5:23 PM
 */

#include <cmath>
#include <iostream>
#include <map>
#include <utility>
#include <cmath>

#include "Dataset.h"
#include "DistanceMetrics.h"
#include "DatasetInstance.h"
#include "Statistics.h"

using namespace std;

pair<bool, double> CheckMissing(unsigned int attributeIndex,
                                DatasetInstance* dsi1,
                                DatasetInstance* dsi2) {
  int numMissing = 0;
  pair<double, double> hasMissing = make_pair(false, false);
  if(dsi1->attributes[attributeIndex] == MISSING_ATTRIBUTE_VALUE) {
    hasMissing.first = true;
    ++numMissing;
  }
  if(dsi2->attributes[attributeIndex] == MISSING_ATTRIBUTE_VALUE) {
    hasMissing.second = true;
    ++numMissing;
  }

  // RELIEF-D
  if(!numMissing) {
    return make_pair(false, 0.0);
  } else {
    return make_pair(true, 2.0 / 3.0);
  }

  pair<bool, double> retValue;
  double diff = 0.0;
  //  unsigned int numLevels = 0;
  //  unsigned int V = 0;
  //  Dataset* ds = dsi1->GetDatasetPtr();
  switch(numMissing) {
    case 0:
      retValue = make_pair(false, 0.0);
      break;
    case 1:
      if(hasMissing.first == true) {
        diff = 1.0 / 3.0;
        //        diff = 1.0 / (double) ds->NumLevels(attributeIndex);
        //        diff = ds->GetProbabilityValueGivenClass(attributeIndex,
        //                                                 dsi2->attributes[attributeIndex],
        //                                                 dsi1->GetClass());
      } else {
        diff = 1.0 / 3.0;
        //        diff = 1.0 / (double) ds->NumLevels(attributeIndex);
        //        diff = ds->GetProbabilityValueGivenClass(attributeIndex,
        //                                                 dsi1->attributes[attributeIndex],
        //                                                 dsi2->GetClass());
      }
      // cout << "One missing value, diff = " << diff << endl;
      retValue = make_pair(true, 1.0 - diff);
      break;
    case 2:
      diff = 1.0 / 3.0;
      //      diff = 1.0 / (double) ds->NumLevels(attributeIndex);
      //      numLevels = ds->NumLevels(attributeIndex);
      //      diff = 0.0;
      //      double PVCI1 = 0.0, PVCI2 = 0.0;
      //      for(V = 0; V < numLevels; ++V) {
      //        PVCI1 = ds->GetProbabilityValueGivenClass(attributeIndex, V, dsi1->GetClass());
      //        PVCI2 = ds->GetProbabilityValueGivenClass(attributeIndex, V, dsi2->GetClass());
      //        cout << PVCI1 << " * " << PVCI2 << " = " << (PVCI1 * PVCI2) << endl;
      //        diff += (PVCI1 * PVCI2);
      //      }
      // cout << "Two missing values, diff = " << diff << endl;
      retValue = make_pair(true, 1.0 - diff);
      break;
  }

  return retValue;
}

pair<bool, double> CheckMissingNumeric(unsigned int numericIndex,
                                       DatasetInstance* dsi1,
                                       DatasetInstance* dsi2) {
  int numMissing = 0;
  pair<double, double> hasMissing = make_pair(false, false);
  if(dsi1->numerics[numericIndex] == MISSING_NUMERIC_VALUE) {
    hasMissing.first = true;
    ++numMissing;
  }
  if(dsi2->numerics[numericIndex] == MISSING_NUMERIC_VALUE) {
    hasMissing.second = true;
    ++numMissing;
  }

  if(numMissing) {
    // ripped from Weka ReliefFAttributeEval.java:difference(), lines 828-836
    double diff = 0.0;
    if(numMissing == 2) {
      return make_pair(true, 1.0);
    } else {
      if(hasMissing.first) {
        pair<double, double> thisMinMax =
                dsi2->GetDatasetPtr()->GetMinMaxForNumeric(numericIndex);
        diff = norm(dsi2->numerics[numericIndex], thisMinMax.first, thisMinMax.second);
      } else {
        pair<double, double> thisMinMax =
                dsi1->GetDatasetPtr()->GetMinMaxForNumeric(numericIndex);
        diff = norm(dsi1->numerics[numericIndex], thisMinMax.first, thisMinMax.second);
      }
      if(diff < 0.5) {
        diff = 1.0 - diff;
      }
      return make_pair(true, diff);
    }
  } else {
    return make_pair(false, 0.0);
  }
}

double norm(double x, double minX, double maxX) {
  if(minX == maxX) {
    return 0;
  } else {
    return(x - minX) / (maxX - minX);
  }
}

double diffAMM(unsigned int attributeIndex,
               DatasetInstance* dsi1,
               DatasetInstance* dsi2) {
  double distance = 0.0;
  pair<bool, double> checkMissing = CheckMissing(attributeIndex, dsi1, dsi2);
  if(checkMissing.first) {
    distance = checkMissing.second;
  } else {
    distance = (double)
            abs((int) dsi1->attributes[attributeIndex] -
                (int) dsi2->attributes[attributeIndex]) * 0.5;
  }
  return distance;
}

double diffGMM(unsigned int attributeIndex,
               DatasetInstance* dsi1,
               DatasetInstance* dsi2) {
  double distance = 0.0;
  pair<bool, double> checkMissing = CheckMissing(attributeIndex, dsi1, dsi2);
  if(checkMissing.first) {
    distance = checkMissing.second;
  } else {
    distance = (dsi1->attributes[attributeIndex] !=
                dsi2->attributes[attributeIndex]) ? 1.0 : 0.0;
  }
  return distance;
}

double diffNCA(unsigned int attributeIndex,
               DatasetInstance* dsi1,
               DatasetInstance* dsi2) {
  double distance = 0.0;
  // TODO: need special missing value checks for NCA metrics
  pair<bool, double> checkMissing = CheckMissing(attributeIndex, dsi1, dsi2);
  if(checkMissing.first) {
    distance = checkMissing.second;
  } else {
  	pair<char, char> alleles =
  			dsi1->GetDatasetPtr()->GetAttributeAlleles(attributeIndex);
  	string a1 = " ";
  	a1[0] = alleles.first;
  	string a2 = " ";
  	a2[0] = alleles.second;
  	map<AttributeLevel, string> genotypeMap;
  	genotypeMap[0] = a1 + a1;
  	genotypeMap[1] = a1 + a2;
  	genotypeMap[2] = a2 + a2;
  	AttributeLevel attrLevel1 = dsi1->attributes[attributeIndex];
  	AttributeLevel attrLevel2 = dsi2->attributes[attributeIndex];
  	string genotype1 = genotypeMap[attrLevel1];
  	string genotype2 = genotypeMap[attrLevel2];
  	map<char, unsigned int> nca1;
  	nca1['A'] = 0; nca1['T'] = 0; nca1['C'] = 0; nca1['G'] = 0;
  	++nca1[genotype1[0]];
  	++nca1[genotype1[1]];
  	map<char, unsigned int> nca2;
  	nca2['A'] = 0; nca2['T'] = 0; nca2['C'] = 0; nca2['G'] = 0;
  	++nca2[genotype2[0]];
  	++nca2[genotype2[1]];
  	map<char, unsigned int>::const_iterator nca1It = nca1.begin();
  	map<char, unsigned int>::const_iterator nca2It = nca2.begin();
  	for(; nca1It != nca1.end(); ++nca1It, ++nca2It) {
  		double nucleotideCount1 = (double) nca1It->second;
  		double nucleotideCount2 = (double) nca2It->second;
  		distance += abs(nucleotideCount1 - nucleotideCount2);
  	}
  }
  return distance;
}

double diffNCA6(unsigned int attributeIndex,
                DatasetInstance* dsi1,
                DatasetInstance* dsi2) {
  double distance = 0.0;
  // TODO: need special missing value checks for NCA metrics
  pair<bool, double> checkMissing = CheckMissing(attributeIndex, dsi1, dsi2);
  if(checkMissing.first) {
    distance = checkMissing.second;
  } else {
  	pair<char, char> alleles =
  			dsi1->GetDatasetPtr()->GetAttributeAlleles(attributeIndex);
  	string a1 = " ";
  	a1[0] = alleles.first;
  	string a2 = " ";
  	a2[0] = alleles.second;
  	map<AttributeLevel, string> genotypeMap;
  	genotypeMap[0] = a1 + a1;
  	genotypeMap[1] = a1 + a2;
  	genotypeMap[2] = a2 + a2;
  	AttributeLevel attrLevel1 = dsi1->attributes[attributeIndex];
  	AttributeLevel attrLevel2 = dsi2->attributes[attributeIndex];
  	string genotype1 = genotypeMap[attrLevel1];
  	string genotype2 = genotypeMap[attrLevel2];
  	map<char, unsigned int> nca1;
  	nca1['A'] = 0; nca1['T'] = 0; nca1['C'] = 0; nca1['G'] = 0;
  	nca1['X'] = 0; nca1['Y'] = 0;
  	++nca1[genotype1[0]];
  	++nca1[genotype1[1]];
  	nca1['X'] = nca1['A'] + nca1['G'];
  	nca1['Y'] = nca1['C'] + nca1['T'];
  	map<char, unsigned int> nca2;
  	nca2['A'] = 0; nca2['T'] = 0; nca2['C'] = 0; nca2['G'] = 0;
  	nca2['X'] = 0; nca2['Y'] = 0;
  	++nca2[genotype2[0]];
  	++nca2[genotype2[1]];
  	nca2['X'] = nca2['A'] + nca2['G'];
  	nca2['Y'] = nca2['C'] + nca2['T'];
  	map<char, unsigned int>::const_iterator nca1It = nca1.begin();
  	map<char, unsigned int>::const_iterator nca2It = nca2.begin();
  	for(; nca1It != nca1.end(); ++nca1It, ++nca2It) {
  		double nucleotideCount1 = (double) nca1It->second;
  		double nucleotideCount2 = (double) nca2It->second;
  		distance += abs(nucleotideCount1 - nucleotideCount2);
  	}
  }
  return distance;
}

double diffKM(unsigned int attributeIndex,
               DatasetInstance* dsi1,
               DatasetInstance* dsi2) {
  double distance = 0.0;

	AttributeLevel dsi1Al =  dsi1->GetAttribute(attributeIndex);
	AttributeLevel dsi2Al =  dsi2->GetAttribute(attributeIndex);
	if(dsi1Al != dsi2Al) {
		if(dsi1->GetDatasetPtr()->GetAttributeMutationType(attributeIndex) ==
				TRANSITION_MUTATION) {
			distance = 1.0;
		}
		if(dsi1->GetDatasetPtr()->GetAttributeMutationType(attributeIndex) ==
				TRANSVERSION_MUTATION) {
			distance = 2.0;
		}
	}

  return distance;
}

double diffManhattan(unsigned int attributeIndex,
                     DatasetInstance* dsi1,
                     DatasetInstance* dsi2) {

  //  double wekaNormDiff = fabs(norm(dsi1->numerics[attributeIndex],
  //                                  minMax.first, minMax.second) -
  //                             norm(dsi2->numerics[attributeIndex],
  //                                  minMax.first, minMax.second));
  // the above code is equivalent; why???
  double distance = 0.0;
  pair<bool, double> checkMissing =
          CheckMissingNumeric(attributeIndex, dsi1, dsi2);
  if(checkMissing.first) {
    distance = checkMissing.second;
  } else {
    pair<double, double> minMax =
            dsi1->GetDatasetPtr()->GetMinMaxForNumeric(attributeIndex);
    distance =
            fabs(dsi1->numerics[attributeIndex] -
                 dsi2->numerics[attributeIndex]) /
            (minMax.second - minMax.first);
  }
  return distance;
}

double diffEuclidean(unsigned int attributeIndex,
                     DatasetInstance* dsi1,
                     DatasetInstance* dsi2) {
  double distance = 0.0;
  pair<bool, double> checkMissing =
          CheckMissingNumeric(attributeIndex, dsi1, dsi2);
  if(checkMissing.first) {
    distance = checkMissing.second;
  } else {
    distance =
    		hypot(dsi1->numerics[attributeIndex], dsi2->numerics[attributeIndex]);
  }
  return distance;
}

double diffPredictedValueTau(DatasetInstance* dsi1, DatasetInstance* dsi2) {

  pair<double, double> minMax =
          dsi1->GetDatasetPtr()->GetMinMaxForContinuousPhenotype();

  //  double diff = fabs(dsi1->GetPredictedValueTau() - dsi2->GetPredictedValueTau());
  double diff =
          fabs(dsi1->GetPredictedValueTau() - dsi2->GetPredictedValueTau()) /
          (minMax.second - minMax.first);

  return diff;
}
