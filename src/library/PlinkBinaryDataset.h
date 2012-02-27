/**
 * \class PlinkBinaryDataset
 *
 *
 * \brief Plink binary PED/BED file format reader.
 *
 * \sa Dataset
 *
 * \author Bill White
 * \version 1.0
 *
 * Contact: bill.c.white@gmail.com
 * Created on: 3/10/11
 */

#ifndef PLINKBINARYDATASET_H
#define	PLINKBINARYDATASET_H

#include "Dataset.h"
#include "Insilico.h"

class PlinkBinaryDataset : public Dataset
{
public:
  PlinkBinaryDataset();
  ~PlinkBinaryDataset() { ; }
private:
  /*************************************************************************//**
   * Load attribute information.
   * \param [in] PLINK bim filename
   * \return success
   ****************************************************************************/
  bool ReadBimFile(std::string bimFilename);
  /*************************************************************************//**
   * Load individual information.
   * \param [in] PLIN fam filename
   * \return success
   ****************************************************************************/
  bool ReadFamFile(std::string famFilename);
  bool LoadSnps(std::string filename);
  std::pair<char, double> GetAttributeMAF(unsigned int attributeIndex);
  AttributeMutationType GetAttributeMutationType(unsigned int attributeIndex);

  unsigned int numInstancesRead;
  unsigned int numAttributesRead;
  unsigned int numClassesRead;

  std::vector<int> instanceIndicesToKeep;
  std::vector<int> missingPhenoLines;

  std::string filenameBase;

  /// for checking attribute values
  std::vector<std::string> validAttributeValues;
  /// missing class values
  std::vector<std::string> missingClassValuesToCheck;
  /// missing attribute values
  std::vector<std::string> missingAttributeValuesToCheck;
};

#endif	/* PLINKBINARYDATASET_H */
