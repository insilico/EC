/*
 * CleanSnpDataset.h - Bill White - 9/22/11
 * 
* Minimalist data set assumes all integer data.
 */

#ifndef CLEAN_SNP_DATASET_H
#define CLEAN_SNP_DATASET_H

#include"Dataset.h"

class CleanSnpDataset : public Dataset
{
public:
  CleanSnpDataset();
  ~CleanSnpDataset() { ; }

  bool LoadSnps(std::string filename, bool doRecodeA);
};

#endif
