/**
 * \class BirdseedData
 *
 * \brief Read Broad's Birdsuite Birdseed-called SNP data
 *
 * \author Bill White
 * \version 1.0
 *
 * Contact: bill.c.white@gmail.com
 * Created on: 2/12/12
 */

#ifndef BIRDSEEDDATA_H_
#define BIRDSEEDDATA_H_

class BirdseedData {
public:
	BirdseedData();
	virtual ~BirdseedData();
	/// Create a new set of Birdseed data with a SNPs file and
	/// optional phenotype file and optional subject names file
	bool LoadData(std::string snpsFile, std::string phenoFile="",
			std::string subjsFile="", std::string includeSnpsFile="",
			std::string excludeSnpsFile="");
	/// Get the subject names/IDs
	std::vector<std::string> GetSubjectNames();
	/// Get the subject labels
	std::vector<std::string> GetSubjectLabels();
	/// Do the subjects have labels?
	bool HasSubjectLabels();
	/// Get the SNP names/IDs
	std::vector<std::string> GetSNPNames();
	/// Get the number of subjects
	int GetNumSubjects();
	/// Get the number of SNPs
	int GetNumSNPs();
	/// Get SNPs for sample at index
	std::vector<int> GetSubjectGenotypes(int subjectIndex);
	/// Get the phenotype at sample index
	int GetSamplePhenotype(int subjectIndex);
	/// Print basic statistics to the console
	void PrintInfo();
	/// Does this data have phenotypes?
	bool HasPhenotypes();
private:
	/// Filename containing birdseed-called SNPs
	std::string snpsFilename;

	/// Filename containing subject names
	std::string subjectLabelsFilename;
	std::vector<std::string> subjectLabels;
	bool hasSubjectLabels;
	std::vector<std::string> subjectNames;

	std::string excludeSnpsFilename;
	std::vector<std::string> excludeSnps;
	bool hasExcludedSnps;
	std::string includeSnpsFilename;
	std::vector<std::string> includeSnps;
	bool hasIncludedSnps;
	/// SNP names
	std::vector<std::string> snpNames;
	/// SNP genotypes
	std::vector<std::vector<int> > snpGenotypes;
	/// SNP genotypes alleles
	std::vector<std::string> snpMajorAllele;
	/// SNP genotypes major allele frequency
	std::vector<double> snpMajorAlleleFreq;

	/// Sample phenotypes
	/// Filename containing subject phenotypes
	std::string phenosFilename;
	/// vector of phenotypes (case-control)
	std::vector<int> phenotypes;
	/// has phenotypes?
	bool hasPhenotypes;
};

#endif /* DGEDATA_H_ */
