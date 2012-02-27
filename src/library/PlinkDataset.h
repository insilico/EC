/**
 * \class PlinkDataset
 *
 * \brief Plink MAP/PED file format reader.
 *
 * \sa Dataset
 *
 * \author Bill White
 * \version 1.0
 *
 * Contact: bill.c.white@gmail.com
 * Created on: 2/1/11
 */

#ifndef PLINKDATASET_H
#define	PLINKDATASET_H

#include <set>
#include <vector>
#include <string>

#include "Dataset.h"
#include "Insilico.h"

 /**
 * \enum MapFileType.
 * PLINK map file types.
 */
typedef enum
{
  MAP3_FILE, /**< map 3 simplified format */
  MAP4_FILE, /**< map 4 standard format */
  ERROR_FILE /**< default */
} MapFileType;

class PlinkDataset : public Dataset
{
public:
  /// Construct a PLINK data set reader. Calls Dataset base class constructor.
  PlinkDataset();
  ~PlinkDataset() { ; }
private:
  bool LoadSnps(std::string filename);
  std::pair<char, double> GetAttributeMAF(unsigned int attributeIndex);
  AttributeMutationType GetAttributeMutationType(unsigned int attributeIndex);
  
  /// base filename for auxiliary files
  std::string filenameBase;
  /// missing class values
  std::vector<std::string> missingClassValuesToCheck;
};

#endif	/* PLINKDATASET_H */

