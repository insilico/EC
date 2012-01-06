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
  ~PlinkRawDataset();
private:
  bool LoadSnps(std::string filename);
  ValueType GetClassValueType(std::string value,
                              std::vector<std::string> missingValues);
  bool GetDiscreteClassLevel(std::string inLevel,
                             std::vector<std::string> missingValues,
                             ClassLevel& outLevel);
  bool GetNumericClassLevel(std::string inLevel,
                            std::vector<std::string> missingValues,
                            NumericLevel& outLevel);
};

#endif	/* PLINKDATASET_H */

