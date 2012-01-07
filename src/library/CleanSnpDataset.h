/**
 * \class CleanSnpDataset
 *
 * \brief Experiemnetal data set reader for large/GWAS data.
 *
 * \sa Dataset
 *
 * \author Bill White
 * \version 1.0
 *
 * Contact: bill.c.white@gmail.com
 * Created on: 9/22/11
 */

#ifndef CLEAN_SNP_DATASET_H
#define CLEAN_SNP_DATASET_H

#include"Dataset.h"

class CleanSnpDataset : public Dataset
{
public:
  CleanSnpDataset();
  ~CleanSnpDataset() { ; }

  bool LoadSnps(std::string filename);
};

#endif
