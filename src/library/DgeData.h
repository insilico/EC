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
	bool LoadData(std::string countsFile, std::string phenoFile);
private:
	/// Filename containing DGE counts
	std::string countsFilename;
	/// Filename containing DGE phenotypes
	std::string phenosFilename;
	/// Gene names
	std::vector<std::string> geneNames;
	/// Digital gene expression counts
	std::vector<std::vector<double> > counts;
	/// Sample names
	std::vector<std::string> sampleNames;
	/// Sample phenotypes
	std::vector<int> phenotypes;
};

#endif /* DGEDATA_H_ */
