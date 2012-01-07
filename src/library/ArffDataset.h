/**
 * \class ArffDataset
 *
 * \brief ARFF file format reader.
 * http://www.cs.waikato.ac.nz/ml/weka/arff.html
 *
 * \sa Dataset
 *
 * \author Bill White
 * \version 1.0
 *
 * Contact: bill.c.white@gmail.com
 * Created on: 2/24/11
 */

#ifndef ARFFDATASET_H
#define	ARFFDATASET_H

#include <vector>
#include <map>
#include <string>

#include "Dataset.h"

 /**
 * \enum ArffAttributeType.
 * ARFF attribute types.
 */
typedef enum
{
  ARFF_NUMERIC_TYPE, /**< continuous levels */
  ARFF_NOMINAL_TYPE, /**< discrete levels */
  ARFF_STRING_TYPE,  /**< string levels */
  ARFF_DATE_TYPE,    /**< date levels */
  ARFF_ERROR_TYPE    /**< unknown type */
} ArffAttributeType;

class ArffDataset : public Dataset
{
public:
  ArffDataset();
  /*************************************************************************\\**
   * Get the ARFF data type of the attribute in columnIndex
   * NOTE: only NOMINAL_TYPE is currently supported.
   *
   * \param [in] columnIndex SNP column index
   * \return data type of the indexed variable
   ****************************************************************************/
  ArffAttributeType GetTypeOf(unsigned int columnIndex);
  /*************************************************************************\\**
   * Print the nominals levels to stdout.
   ****************************************************************************/
  void PrintNominalsMapping();
  ~ArffDataset();
private:
  bool LoadSnps(std::string filename);
  bool GetAttributeLevel(std::string inLevel,
                         std::vector<std::string> missingValues,
                         AttributeLevel& outLevel);
  bool GetDiscreteClassLevel(std::string inLevel,
                             std::vector<std::string> missingValues,
                             ClassLevel& outLevel);
  bool GetNumericClassLevel(std::string inLevel,
                            std::vector<std::string> missingValues,
                            NumericLevel& outLevel);

  /// ARFF relation name
  std::string relationName;
  /// vector of attribute types
  std::vector<ArffAttributeType> attributeTypes;
  /// map of attribute names to valid nominal values
  std::map<std::string, std::vector<std::string> > nominalValues;

  /// missing class values
  std::vector<std::string> missingClassValuesToCheck;
  /// missing attribute values
  std::vector<std::string> missingAttributeValuesToCheck;
};

#endif	/* ARFFDATASET_H */
