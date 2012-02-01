/**
 * \class DgeData
 *
 * \brief Digital gene expression data
 *
 * \author Bill White
 * \version 1.0
 *
 * Contact: bill.c.white@gmail.com
 * Created on: 1/18/12
 */

#ifndef DGEDATA_H_
#define DGEDATA_H_

class DgeData {
public:
	DgeData();
	virtual ~DgeData();
	/// Create a new set of DGE data with a counts file and a phenotype file
	bool LoadData(std::string countsFile, std::string phenoFile, std::string normsFile="");
	/// Get the sample names/IDs
	std::vector<std::string> GetSampleNames();
	/// Get the gene names/IDs
	std::vector<std::string> GetGeneNames();
	/// Get the min and max values for gene at index
	std::pair<double, double> GetGeneMinMax(int geneIndex);
	/// Get the number of samples
	int GetNumSamples();
	/// Get the number of genes
	int GetNumGenes();
	/// Get sample counts for sample at index
	std::vector<double> GetSampleCounts(int sampleIndex);
	/// Get the phenotype at sample index
	int GetSamplePhenotype(int sampleIndex);
	/// Get the normalization factors
	std::vector<double> GetNormalizationFactors();
	/// Print the Sample statistics to the console
	void PrintSampleStats();
private:
	/// Filename containing DGE counts
	std::string countsFilename;
	/// Filename containing DGE phenotypes
	std::string phenosFilename;
	/// Filename containing DGE normalization factors
	std::string normsFilename;
	/// Are we using normalization?
	bool hasNormFactors;
	/// Vector of (optional) normalization factors for each sample
	std::vector<double> normFactors;
	/// Gene names
	std::vector<std::string> geneNames;
	/// Digital gene expression counts
	std::vector<std::vector<double> > counts;
	/// Sample names
	std::vector<std::string> sampleNames;
	/// Sample phenotypes
	std::vector<int> phenotypes;
	/// Min and max count for genes
	std::vector<std::pair<double, double> > minMaxGeneCounts;
	/// Min and max values for samples
	std::vector<std::pair<double, double> > minMaxSampleCounts;
	/// Zero count sample indices
	std::vector<std::vector<int> > sampleZeroes;
};

#endif /* DGEDATA_H_ */
