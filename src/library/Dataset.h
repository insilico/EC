/**
 * \class Dataset
 *
 * \brief Base class for collections of instances containing attributea and class.
 *
 * Added interaction infomation week of 4/18-26/06
 * Totally redone for McKinney Lab. February 2011.
 *
 * \author Bill White
 * \version 1.0
 *
 * Contact: bill.c.white@gmail.com
 * Created on: 6/14/05
 */

#ifndef DATASET_H
#define DATASET_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <climits>

#include "DatasetInstance.h"
#include "Insilico.h"

// GSL random number generator base class
#include "GSLRandomFlat.h"

class DgeData;
class BirdseedData;

class Dataset
{
public:
  /// Construct a default data set.
  Dataset();
  /// Destruct all dynamically allocated memory.
  virtual ~Dataset();
  /*************************************************************************//**
   * Load the data set from "raw data".
   * \param [in] dataMatrix reference to a matrix of SNP values 0-1-2-missing
   * \param [in] classLabels reference to a vector of case-control labels
   * \param [in] attrNames reference to a vector of attribute names
   * \return success
   ****************************************************************************/
  bool LoadDataset(std::vector<std::vector<int> >& dataMatrix,
  		std::vector<int>& classLabels, std::vector<std::string>& attrNames);
  /*************************************************************************//**
   * Load the data set from files passed as parameters.
   * \param [in] snpFilename discrete values (SNPs) filename
   * \param [in] doRecodeA perform recodeA encoding after reading
   * \param [in] numericsFilename continuous values (numerics) filename or empty string
   * \param [in] altPhenoFilename alternate class (phenotype) filename or empty string
   * \param [in] ids vector of possibly empty IDs to match in auxiliary files
   * \return success
   ****************************************************************************/
  bool LoadDataset(std::string snpsFilename,
                   std::string numericsFilename,
                   std::string altPhenoFilename,
                   std::vector<std::string> ids);
  /*************************************************************************//**
   * Load the data set from DGE data.
   * \param [in] dgeData pointer to a digital gene expression (DGE) data object
   * \return success
   ****************************************************************************/
  bool LoadDataset(DgeData* dgeData);
  /*************************************************************************//**
   * Load the data set from Birdseed data.
   * \param [in] birdseedData pointer to a Birdseed-called SNP data object
   * \return success
   ****************************************************************************/
  bool LoadDataset(BirdseedData* birdseedData);
  /*************************************************************************//**
   * Get the attribute value at row, column.
   * Same as instance index, attribute index.
   * \param [in] row instance row
   * \param [in] col attribute column
   * \param [out] attrVal attribute value
   * \return success
   ****************************************************************************/
  bool GetAttributeRowCol(unsigned int row, unsigned int col,
                          AttributeLevel& attrVal);
  /*************************************************************************//**
   * Get the numeric value at row, column.
   * Same as instance index, numeric index.
   * \param [in] row instance row
   * \param [in] col numeric column
   * \param [out] numVal numeric value
   * \return success
   ****************************************************************************/
  bool GetNumericRowCol(unsigned int row, unsigned int col,
                        NumericLevel& numVal);
  /*************************************************************************//**
   * Write the dataset to a new filename, respecting masked attributes
   * and numerics and class/phenotype data type.
   * \param [in] newDatasetFilename new dataset filename
   * \return success
   ****************************************************************************/
  bool WriteNewDataset(std::string newDatasetFilename,
                       OutputDatasetType outputDatasetType);
  /*************************************************************************//**
   * Extracts top N attributes based on a file of attribute scores
   * and writes a new dataset.
   * Revised 10/3/11 for numerics and continuous class/phenotypes.
   * \param  [in] scoresFilename filename of attribute scores and names
   * \param  [in] topN top N attributes
   * \param  [in] newDatasetFilename filename of new dataset
   * \return success
   *
   ****************************************************************************/
  bool ExtractAttributes(std::string scoresFilename,
                         unsigned int topN,
                         std::string newDatasetFilename);
  /*************************************************************************//**
   * Swap two attributes/columns in the dataset.
   * \param [in] a1 attribue index 1
   * \param [in] a2 attribue index 2
   * \return success
   ****************************************************************************/
  bool SwapAttributes(unsigned int a1, unsigned int a2);
  /*************************************************************************//**
   * Return the number of discrete plus continuous variables in the data set.
   * The number does not include masked variables removed.
   * \return number of discrete plus continuous variables
   ****************************************************************************/
  unsigned int NumVariables();
  /*************************************************************************//**
   * Returns the names of discrete and continuous variables in the data set.
   * \return vector of names as strings
   ****************************************************************************/
  std::vector<std::string> GetVariableNames();
  /// Returns the number of instances in the data set.
  virtual unsigned int NumInstances();
  /*************************************************************************//**
   * Returns a pointer to a dataset instance selected by index.
   * \param [in] index index of instance
   * \return pointer to an instance
   ****************************************************************************/
  DatasetInstance* GetInstance(unsigned int index);
  /*************************************************************************//**
   * Returns a pointer to a randomly chosen data set instance.
   * The random number generator is set to give values in range
   * of instance indexes.
   * \return pointer to a data set instance
   ****************************************************************************/
  DatasetInstance* GetRandomInstance();
  /*************************************************************************//**
   * Get all instance IDs.
   * \return vector of instance IDs
   ****************************************************************************/
  std::vector<std::string> GetInstanceIds();
  /*************************************************************************//**
   * Get the instance index from the instance ID.
   * \param [in] ID string ID
   * \param [out] instanceIndex instance index
   * \return success
   ****************************************************************************/
  bool GetInstanceIndexForID(std::string ID, unsigned int& instanceIndex);
  /// Return the number of unmasked discrete attributes in the data set.
  virtual unsigned int NumAttributes();
  /*************************************************************************//**
   * Return the discrete (SNP) attribute names.
   * \return vector of attribute names
   ****************************************************************************/
  std::vector<std::string> GetAttributeNames();
  /*************************************************************************//**
   * Loads the referenced vector with an attribute's values (column).
   * from the dataset
   * \param [in] attributeIndex attribute index
   * \param [out] attributeValues reference to a a vector allocated by the caller
   * \return success
   ****************************************************************************/
  bool GetAttributeValues(unsigned int attributeIndex,
                          std::vector<AttributeLevel>& attributeValues);
  /*************************************************************************//**
   * Loads the referenced vector with an attribute's values (column)
   * from the dataset.
   * \param [in] attributeName attribute name
   * \param [out] attributeValues reference to a a vector allocated by the caller
   * \return success
   ****************************************************************************/
  bool GetAttributeValues(std::string attributeName,
                          std::vector<AttributeLevel>& attributeValues);
  /// Get the filename SNPs were read from.
  std::string GetSnpsFilename();
  /*************************************************************************//**
   * Looks up original attribute index from attribute name.
   * \param [in] attributeName attribute name
   * \return attribute index or INVALID_INDEX
   ****************************************************************************/
  unsigned int GetAttributeIndexFromName(std::string attributeName);
  /// Does the data set have genotype variables?
  bool HasGenotypes();
  /// Does the data set have allelic information for genotypes?
  bool HasAllelicInfo();
  /*************************************************************************//**
   * Get attribute value for attribute name at instance index.
   * \param [in] instanceIndex instance index
   * \param [in] name attribute name
   * \return attributevalue
   ****************************************************************************/
  AttributeLevel GetAttribute(unsigned instanceIndex, std::string name);
  /*************************************************************************//**
   * Get attribute major and minor alleles
   * \param [in] attribute index
   * \return pair (major allele, minor allele frequency)
   ****************************************************************************/
  std::pair<char, char> GetAttributeAlleles(unsigned int attributeIndex);
  /*************************************************************************//**
   * Get attribute minor allele and frequency.
   * \param [in] attribute index
   * \return pair (minor allele, minor allele frequency)
   ****************************************************************************/
  virtual std::pair<char, double> GetAttributeMAF(unsigned int attributeIndex);
  /*************************************************************************//**
   * Get attribute mutation type.
   * \param [in] attribute index
   * \return mutation type (transition, transversion, unknown)
   ****************************************************************************/
  virtual AttributeMutationType
  GetAttributeMutationType(unsigned int attributeIndex);
  /*************************************************************************//**
   * Get integer value for string genotype.
   * \param [in] genotype genotype string
   * \param [out] newAttr new attribute value
   * \return success
   ****************************************************************************/
  bool GetIntForGenotype(std::string genotype, AttributeLevel& newAttr);
  /*************************************************************************//**
   * Returns the number of levels in a given attribute index.
   * \param [in] index attribute index
   * \return number of levels
   ****************************************************************************/
  unsigned int NumLevels(unsigned int index);
  /// Return the number of unmasked discrete attributes in the data set.
  unsigned int NumNumerics();
  /*************************************************************************//**
   * Return the numeric attribute names.
   * \return vector of attribute names
   ****************************************************************************/
  std::vector<std::string> GetNumericsNames();
  /*************************************************************************//**
   * Get the minimum and maximum values for a numeric at index.
   * \param [in] numericIdx numeric index
   * \return minimum/maximum pair
   ****************************************************************************/
  std::pair<double, double> GetMinMaxForNumeric(unsigned int numericIdx);
  /*************************************************************************//**
   * Get the mean/average of numeric at index.
   * \param [in] numericIdx numeric index
   * \return average value of numeric attribute at index
   ****************************************************************************/
  double GetMeanForNumeric(unsigned int numericIdx);
  /// Does the data set have numeric variables? setter/getter
  bool HasNumerics();
  void HasNumerics(bool setHasNumerics);
  /*************************************************************************//**
   * Get numeric value for numeric name at instance index.
   * \param [in] instanceIndex instance index
   * \param [in] name numeric name
   * \return numeric value at index
   ****************************************************************************/
  NumericLevel GetNumeric(unsigned int instanceIndex, std::string name);
  /*************************************************************************//**
   * Loads the referenced vector with a numeric's values (column)
   * from the dataset.
   * \param [in] numericName numeric name
   * \param [out] numericValues reference to a a vector allocated by the caller
   * \return success
   ****************************************************************************/
  bool GetNumericValues(std::string numericName,
                        std::vector<NumericLevel>& numericValues);
  /// Get the filename numerics were read from.
  std::string GetNumericsFilename();
  /*************************************************************************//**
   * Looks up original numeric index from numeric name.
   * \param [in] numericName numeric name
   * \return attribute index or INVALID_INDEX
   ****************************************************************************/
  unsigned int GetNumericIndexFromName(std::string numericName);
  /// Get the number of classes in the data set.
  unsigned int NumClasses();
  /// Get the class column as read from the file.
  unsigned int GetClassColumn();
  /*************************************************************************//**
   * Loads the referenced vector with the dataset's class labels.
   * \param [out] classValues reference to a a vector allocated by the caller
   * \return success
   ****************************************************************************/
  bool GetClassValues(std::vector<ClassLevel>& classValues);
  /*************************************************************************//**
   * Get a map from class levels to a vector of instance indices.
   * \return map of class => instance indices
   ****************************************************************************/
  const std::map<ClassLevel, std::vector<unsigned int> >& GetClassIndexes();
  /// Does the data set have alternate phenotypes loaded?
  bool HasAlternatePhenotypes();
  void HasAlternatePhenotypes(bool setHasAlternatePhenotypes);
  /// Get the alternate phenotype filename.
  std::string GetAlternatePhenotypesFilename();
  /// Does the data set have continuous phenotypes?
  bool HasContinuousPhenotypes();
  /// Does the data set have any valid phenotypes?
  bool HasPhenotypes();
  /*************************************************************************//**
   * Get the minumum and maximum values for the continuous phenotype.
   * \return minimum/maximum pair
   ****************************************************************************/
  std::pair<double, double> GetMinMaxForContinuousPhenotype();
  /// Print the entire data set in compact format.
  void Print();
  /// Print basic statstics abou the data set - discrete/SNPs only.
  void PrintStats();
  /// Print statistics about the data set including numerics.
  void PrintNumericsStats();
  /// Print very simple statistics abou the data set with no formatting.
  void PrintStatsSimple(std::ostream& outStream=std::cout);
  /// Print class index information.
  void PrintClassIndexInfo(std::ostream& outStream=std::cout);
  /// Print missing value statistics.
  void PrintMissingValuesStats();
  /// Prit attribute level counts.
  void PrintLevelCounts();
  /// 
  /*************************************************************************//**
   * Write attribute level counts to a text file.
   * \param [in] levelsFilename filename to write levels to
   ****************************************************************************/
  void WriteLevelCounts(std::string levelsFilename);
  /// Print unique attribute levels seen.
  void PrintAttributeLevelsSeen();
  /*************************************************************************//**
   * Removes the variable name from consideration in any data set operations.
   * \param [in] variableName variable name
   * \return success
   ****************************************************************************/
  bool MaskRemoveVariable(std::string variableName);
  /*************************************************************************//**
   * Removes the attribute name from consideration in any data set operations.
   * \param [in] attributeName attribute name
   * \param [in] attrType attribute type
   * \return success
   ****************************************************************************/
  bool MaskRemoveVariableType(std::string variableName,
  		AttributeType varType);
  /*************************************************************************//**
   * Determines if the named variable is in the current masked data set.
   * \param [in] attributeName attribute name
   * \param [in] attributeType attribute type
   * \return true if discrete attribute name is being considered in operations.
   ****************************************************************************/
  bool MaskSearchVariableType(std::string variableName,
  		AttributeType attrType);
  /*************************************************************************//**
   * Mark all attributes for inclusion in data set operations.
   * \param [in] attrType attribute type
   * \return success
   ****************************************************************************/
  bool MaskIncludeAllAttributes(AttributeType attrType);
  /*************************************************************************//**
   * Return a vector of all the attribute indices under consideration.
   * \param attrType attribute type
   * \return vector of indices into currently considered discrete attributes
   ****************************************************************************/
  std::vector<unsigned int> MaskGetAttributeIndices(AttributeType attrType);
  /*************************************************************************//**
   * Return a map of attribute name to attribute index of attributes to include.
   * \param [in] attrType attribute type
   * \return attributes mask: name->index
   ****************************************************************************/
  const std::map<std::string, unsigned int>&
  MaskGetAttributeMask(AttributeType attrType);
  /*************************************************************************//**
   * Return a vector of all the variable names under consideration.
   * \return vector of discrete and numeric variable
   ****************************************************************************/
  std::vector<std::string> MaskGetAllVariableNames();
  /*************************************************************************//**
   * Removes the instance from consideration in any data set operations.
   * \param [in] instanceId instance ID
   * \return success
   ****************************************************************************/
  bool MaskRemoveInstance(std::string instanceId);
  /*************************************************************************//**
   * Determines if the names Instance is in the current masked dataaset.
   * \param [in] instanceID instance ID
   * \return true if instance ID is in the dataset, considering instance mask
   ****************************************************************************/
  bool MaskSearchInstance(std::string instanceId);
  /*************************************************************************//**
   * Mark all instances for inclusion in algorithms.
   * \return success
   ****************************************************************************/
  bool MaskIncludeAllInstances();
  /*************************************************************************//**
   * Return a vector of all the instance indices under consideration.
   * \retrun vector of indices into current instances
   ****************************************************************************/
  std::vector<unsigned int> MaskGetInstanceIndices();
  /*************************************************************************//**
   * Return a vector of all the instance ids under consideration.
   * \return vector of ids of currently included instances
   ****************************************************************************/
  std::vector<std::string> MaskGetInstanceIds();
  /*************************************************************************//**
   * Return a map of instance name to instance index of instances to include.
   * \return instances mask: instance ID=>vector of instance indices
   ****************************************************************************/
  const std::map<std::string, unsigned int>& MaskGetInstanceMask();
  /*************************************************************************//**
   * Save the current masks for later restore.
   * \return success
   ****************************************************************************/
  bool MaskPushAll();
  /*************************************************************************//**
   * Restore the masks previously pushed.
   * \return success
   ****************************************************************************/
  bool MaskPopAll();
  /*************************************************************************//**
   * Saved the unmasked attributes as a tab-delimited text file.
   * \param [in] newDatasetFilename new data set filename
   * \return success
   ****************************************************************************/
  bool MaskWriteNewDataset(std::string newDatasetFilename);
  /// Print mask statistics.
  void PrintMaskStats();
  /*************************************************************************//**
   * Perform and report SNP diagnostic test information.
   * \param [in] logFilename log filename
   * \param [in] globalGenotypeThreshold genotype count threshold
   * \param [in] cellThreshold x^2 cell count threshold
   ****************************************************************************/
  void RunSnpDiagnosticTests(std::string logFilename,
                             double globalGenotypeThreshold = 0.01,
                             unsigned int cellThreshold = 5);
  /*************************************************************************//**
   * Calculate whether passed genotype counts are in HWE.
   * \param genotypeCounts vector of genotype counts: AA, Aa, aa
   * \return counts are in HWE?
   ****************************************************************************/
  bool CheckHardyWeinbergEquilibrium(std::vector<unsigned int>& chkGenotypeCounts);
  /*************************************************************************//**
   * This code implements an exact SNP test of Hardy-Weinberg Equilibrium.
   * As described in Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note
   * on Exact Tests of Hardy-Weinberg Equilibrium. American Journal of Human
   * Genetics: 76. Written by Jan Wigginton.
   * \param [in] obs_hets observed heterozygotes
   * \param [in] obs_hom1 observed homozygotes type 1
   * \param [in] obs_hom2 homozygotes type 2
   * \return HWE value
   ****************************************************************************/
  double SNPHWE(int obs_hets, int obs_hom1, int obs_hom2);
  /*************************************************************************//**
   * Get the probability of a class value in the data set.
   * \param  thisClass class value
   * \return probability
   ****************************************************************************/
  double GetClassProbability(ClassLevel thisClass);
  /*************************************************************************//**
   * Get the probability of an attribute value at an attribute index.
   * \param [in] attributeIndex attribute index
   * \param [in] A attribute value
   * \param [in] classValue class value
   * \return probability of the value in attribute given class
   ****************************************************************************/
  double GetProbabilityValueGivenClass(unsigned int attributeIndex,
                                       AttributeLevel A, ClassLevel classValue);
  /*************************************************************************//**
   * Calculate and display interaction information for all attribute combinations.
   ****************************************************************************/
  void AttributeInteractionInformation();
  /*************************************************************************//**
   * Calculate all the information needed to construct the interaction diagram.
   * \param [out] results map of attribute combinations to results
   ****************************************************************************/
  void
  CalculateInteractionInformation(std::map<std::pair<int, int>,
  		std::map<std::string, double> > & results);
  /*************************************************************************//**
   * Calculate the GAIN matrix to run snprank on this data set.
   * Uses OpenMP to calculate matrix entries in parallel threads.
   * \param [out] gainMatrix pointer to an allocated n x n matrix,
   *                         n = number of attributes
   * \return success
   ****************************************************************************/
  bool CalculateGainMatrix(double** gainMatrix,
  		std::string matrixFilename="");
  /*************************************************************************//**
   * Calculate the instance-to-instance distance matrix for this data set.
   * Uses OpenMP to calculate matrix entries in parallel threads.
   * \param [out] distanceMatrix pointer to an allocated m x m matrix,
   *                         m = number of instances
   * \param [in] distanceMatrixFilename filename to write matrix
   * \return success
   ****************************************************************************/
  bool CalculateDistanceMatrix(double** distanceMatrix,
  		std::string matrixFilename="");
  /*************************************************************************//**
   * Calculate the instance-to-instance distance matrix for this data set.
   * Uses OpenMP to calculate matrix entries in parallel threads.
   * \param [out] distanceMatrix vector of vectors of double: m x m matrix,
   *                         m = number of instances
   * \return success
   ****************************************************************************/
  bool CalculateDistanceMatrix(std::vector<std::vector<double> >& distanceMatrix);
  /*************************************************************************//**
   * Compute the distance between two DatasetInstances.
   * \param [in] dsi1 pointer to DatasetInstance 1
   * \param [in] dsi2 pointer to DatasetInstance 2
   * \return distance
   ****************************************************************************/
  double ComputeInstanceToInstanceDistance(DatasetInstance* dsi1,
                                           DatasetInstance* dsi2);
  /*************************************************************************//**
   * Set the he distance metrics used to compute instance-to-instance distances.
   * \param [in] snpMetric name of SNP metric
   * \param [in] numMetric name of the numeric metric
   * \return distance
   ****************************************************************************/
  bool SetDistanceMetrics(std::string newSnpMetric,
  		std::string newNumMetric="manhattan");
protected:
  /*************************************************************************//**
   * Load SNPs from file using the data set filename.
   * \param [in] filename SNPs filename
   * \param [in] deRecodeA perform a recodeA operation after reading raw data?
   * \return success
   ****************************************************************************/
  virtual bool LoadSnps(std::string filename);
  /// Update level counts for all instances by calling UpdateLevelCounts(inst)
  void UpdateAllLevelCounts();
  /// Create dummy alleles from genotypes for data sets that have no allele info
  void CreateDummyAlleles();
  /*************************************************************************//**
   * Update all attribute level counts from one data set instance.
   * Updates levelCountsByClass.
   * \param [in] dsi pointer to a data set instance
   ****************************************************************************/
  void UpdateLevelCounts(DatasetInstance* dsi);
  /*************************************************************************//**
   * Load numerics (continuous attributes) from a file set in the constructor.
   * \param [in] filename numerics data filename in PLINK covar format
   * \return success
   ****************************************************************************/
  bool LoadNumerics(std::string filename);
  /*************************************************************************//**
   * Loads the referenced vector with an numeric's values (column).
   * from the dataset
   * \param [in] numericIndex numeric index
   * \param [out] numericValues reference to a a vector allocated by the caller
   * \return success
   ****************************************************************************/
  bool GetNumericValues(unsigned int numericIndex,
                        std::vector<NumericLevel>& numericValues);
  /*************************************************************************//**
   * Load alternate phenotype/class values from a plink covariate .cov file.
   * Format described here:
   * http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#covar
   * MAJOR CHANGES: for continuous phenotypes/class - 9/29/11
   * \param [in] filename alternate phenotype data filename in PLINK covar format
   * \return success
   ****************************************************************************/
  bool LoadAlternatePhenotypes(std::string filename);
  /*************************************************************************//**
   * Is the passed instance ID loadable (not filtered).
   * \param [in] ID instance ID
   * \return [out] success
   ****************************************************************************/
  bool IsLoadableInstanceID(std::string ID);

  /*************************************************************************//**
   * Compute the discrete difference in an attribute between two instances.
   * \param [in] attributeIndex index into vector of all attributes
   * \param [in] dsi1 pointer to DatasetInstance 1
   * \param [in] dsi2 pointer to DatasetInstance 2
   * \return diff(erence)
   ****************************************************************************/
  double (*snpDiff)(unsigned int attributeIndex,
                    DatasetInstance* dsi1,
                    DatasetInstance* dsi2);
  /*************************************************************************//**
   * Compute the continuous difference in an attribute between two instances.
   * \param [in] attributeIndex index into vector of all attributes
   * \param [in] dsi1 pointer to DatasetInstance 1
   * \param [in] dsi2 pointer to DatasetInstance 2
   * \return diff(erence)
   ****************************************************************************/
  double (*numDiff)(unsigned int attributeIndex,
                    DatasetInstance* dsi1,
                    DatasetInstance* dsi2);

  /// the name of discrete diff(erence) function
  std::string snpMetric;
  /// the name of continuous diff(erence) function
  std::string numMetric;

  /// file from which the discrete attributes (SNPSs) were read
  std::string snpsFilename;
  /// does the data set contain any genotypes?
  bool hasGenotypes;
  /// discrete attribute names read from file
  std::vector<std::string> attributeNames;
  /// attribute values/levels counts
  std::vector<std::map<AttributeLevel, unsigned int> > levelCounts;
  /// attribute values/levels counts by discrete class
  std::vector<std::map<std::pair<AttributeLevel, ClassLevel>, unsigned int> > levelCountsByClass;
  /// unique attribute values/levels read from file
  std::vector<std::set<std::string> > attributeLevelsSeen;
  /// allele1, allele2
  std::vector<std::pair<char, char> > attributeAlleles;
  /// allele->count
  std::vector<std::map<char, unsigned int> > attributeAlleleCounts;
  /// minor allele, minor allele frequency
  std::vector<std::pair<char, double> > attributeMinorAllele;
  /// Does this data set have alelelic information?
  bool hasAllelicInfo;
  /// genotype->count
  std::vector<std::map<std::string, unsigned int> > genotypeCounts;
  /// Keep mutation type for all attributes.
  std::vector<AttributeMutationType> attributeMutationTypes;
  /// Lookup table for mutation type.
  std::map<std::pair<char, char>, AttributeMutationType> attributeMutationMap;

  /// file from which the continuous attributes were read
  std::string numericsFilename;
  /// does the data set contain any continuous attributes?
  bool hasNumerics;
  /// IDs associated with the numerics read from file
  std::vector<std::string> numericsIds;
  /// the minimum and maximum value for each continuous attribute
  std::vector< std::pair<NumericLevel, NumericLevel> > numericsMinMax;
  /// continuous attribute names read from file
  std::vector<std::string> numericsNames;

  /// Does the data set contain phenotypes?
  bool hasPhenotypes;
  /// file from which the alternate phenotypes (class labels) were read
  std::string alternatePhenotypesFilename;
  /// does the data set contain alternate phenotypes?
  bool hasAlternatePhenotypes;
  /// IDs associated with the phenotypes/classes read from file
  std::vector<std::string> phenotypesIds;
  /// does the data set contain continuous phenotypes?
  bool hasContinuousPhenotypes;
  /// the minimum and maximum value for each continuous phenotype
  std::pair<NumericLevel, NumericLevel> continuousPhenotypeMinMax;

  /// vector of pointers to all instances in the data set
  std::vector<DatasetInstance*> instances;
  /// IDs associated with the instances read from file
  std::vector<std::string> instanceIds;
  /// IDs of instances to load from numeric and/or phenotype files
  std::vector<std::string> instanceIdsToLoad;
  /// missing discrete values and their instance indices
  std::map<std::string, std::vector<unsigned int> > missingValues;
  /// missing continuous values and their instance indices
  std::map<std::string, std::vector<unsigned int> > missingNumericValues;

  /// class column from the original data set
  unsigned int classColumn;
  /// class values mapped to instance indices
  std::map<ClassLevel, std::vector<unsigned int> > classIndexes;

  /***
   * Masks specify the columns from the data set being considered
   * when algorithms call methods on this object:
   * key = attribute name, value = original index into all.
   */
  std::map<std::string, unsigned int> attributesMask;
  std::map<std::string, unsigned int> numericsMask;
  std::map<std::string, unsigned int> instancesMask;
  /// masks can be temporarily pushed and popped
  std::map<std::string, unsigned int> attributesMaskPushed;
  std::map<std::string, unsigned int> numericsMaskPushed;
  std::map<std::string, unsigned int> instancesMaskPushed;
  bool maskIsPushed;

  /// random number generator classes use GNU Scienitifc Library (GSL)
  GSLRandomFlat* rng;
};

#endif // DATASET_H
