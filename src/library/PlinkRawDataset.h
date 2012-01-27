/**
 * \class PlinkRawDataset
 *
 * \brief Plink recodeA/RAW file format reader.
 *
 * \sa Dataset
 *
 * \author Bill White
 * \version 1.0
 *
 * Contact: bill.c.white@gmail.com
 * Created on: 2/24/11
 */

#ifndef PLINKRAWDATASET_H
#define	PLINKRAWDATASET_H

#include <string>
#include <vector>

#include "Dataset.h"

class PlinkRawDataset : public Dataset
{
public:
  PlinkRawDataset();
  ~PlinkRawDataset() { ; }
private:
  bool LoadSnps(std::string filename);
};

#endif	/* PLINKDATASET_H */

