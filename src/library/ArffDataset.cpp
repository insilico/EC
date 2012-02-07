/*
 * ArffDataset.cpp - Bill White
 * see ArffDataset.h and http://www.cs.waikato.ac.nz/ml/weka/arff.html
 */

#include <string>
#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <sstream>

#include <boost/lexical_cast.hpp>

#include "Dataset.h"
#include "DatasetInstance.h"
#include "StringUtils.h"
#include "ArffDataset.h"
#include "Insilico.h"

using namespace std;
using namespace insilico;
using namespace boost;

ArffDataset::ArffDataset() {
	::Dataset();
	missingAttributeValuesToCheck.push_back("?");
	missingClassValuesToCheck.push_back("?");
}

bool ArffDataset::LoadSnps(string filename) {
	snpsFilename = filename;
	ifstream dataStream(snpsFilename.c_str());
	if (!dataStream.is_open()) {
		cerr << "ERROR: Could not open dataset: " << snpsFilename << endl;
		return false;
	}
	cout << Timestamp() << "ArffDataset: Reading lines from " << snpsFilename
			<< endl;
	string line;
	int firstSpace = -1, secondSpace = -1;
	string attributeName = "";
	int attributeIndex = 0;
	int numericsIndex = 0;
	string attributeType = "";
	string classTypeString = "";
	ClassType classType = NO_CLASS_TYPE;
	unsigned int lineNumber = 0;
	double minPheno = 0.0, maxPheno = 0.0;
	while (getline(dataStream, line)) {
		++lineNumber;
		string trimmedLine = trim(line);
		// skip blank lines
		if (trimmedLine.size() < 1) {
			continue;
		}
		// check the first character of the line
		switch (trimmedLine.at(0)) {
		case '%':
			// skip comment lines
			continue;
		case '@':
			// relation, attribute or data
			firstSpace = trimmedLine.find(" ");
			string keyword = to_upper(trimmedLine.substr(1, firstSpace - 1));
			// cout << "keyword => " << keyword << endl;
			if (keyword == "RELATION") {
				relationName = keyword;
			}
			if (keyword == "ATTRIBUTE") {
				secondSpace = trimmedLine.find(" ", firstSpace + 1);
				attributeName = trimmedLine.substr(firstSpace + 1,
						secondSpace - firstSpace - 1);
				// cout << "DEBUG attribute name: [" << attributeName << "]" << endl;
				if (to_upper(attributeName) == "CLASS") {
					classColumn = attributeIndex;
					cout << Timestamp() << "Class column detect: " << classColumn << endl;
					classTypeString = to_upper(trimmedLine.substr(secondSpace + 1));
					if (classTypeString == "NUMERIC") {
						hasContinuousPhenotypes = true;
						classType = CONTINUOUS_CLASS_TYPE;
						cout << Timestamp() << "Detected continuous phenotype" << endl;
					} else {
						hasContinuousPhenotypes = false;
						classType = CASE_CONTROL_CLASS_TYPE;
						cout << Timestamp() << "Detected case-control phenotype" << endl;
					}
				} else {
					attributeNames.push_back(attributeName);
					attributeType = to_upper(trimmedLine.substr(secondSpace + 1));
					if (attributeType == "STRING") {
						cerr << "ERROR: STRING attributes are not yet supported" << endl;
						return false;
						attributeTypes.push_back(ARFF_STRING_TYPE);
					}
					if (attributeType == "DATE") {
						cerr << "ERROR: DATE attributes are not yet supported" << endl;
						return false;
						attributeTypes.push_back(ARFF_DATE_TYPE);
					}
					if (attributeType == "NUMERIC") {
						attributeTypes.push_back(ARFF_NUMERIC_TYPE);
						numericsMask[attributeName] = numericsIndex;
						numericsNames.push_back(attributeName);
						++numericsIndex;
					} else {
						// must be nominal type - add nominal values to map
						attributeTypes.push_back(ARFF_NOMINAL_TYPE);
						vector<string> tokens;
						split(tokens, trimmedLine, "{");
						vector<string>::const_iterator it = tokens.end() - 1;
						string nominalsListWithCurly = *it;
						string nominalsList = trim(
								nominalsListWithCurly.substr(0,
										nominalsListWithCurly.size() - 1));
						vector<string> nominals;
						split(nominals, nominalsList, ",");
						// GENETICS CHECK HERE for plink recodeA encoding
						if (nominals.size() != 3) {
							cerr << "ERROR: This dataset is currently unsupported. SNP data "
									<< "must be encoded with {0, 1, 2} for {homozygous1, "
									<< "heterzygote, homozygous2} respectively. The following "
									<< "attributes were read successfully" << endl;
							PrintNominalsMapping();
							return false;
						}
						// GENETICS CHECK HERE
						if ((nominals[0] == "0") && (nominals[1] == "1")
								&& (nominals[2] == "2")) {
							nominalValues[attributeName] = nominals;
							attributesMask[attributeName] = attributeIndex;
						} else {
							cerr << "ERROR: This dataset is currently unsupported. SNP data "
									<< "must be encoded with {0, 1, 2} for {homozygous1, "
									<< "heterzygote, homozygous2} respectively. The following "
									<< "attributes were read successfully" << endl;
							PrintNominalsMapping();
							return false;
						} // end genetics check
					} // end nominal
				} // end class or attribute
				++attributeIndex;
			} // keyword = attribute

			// the rest of the file is instances
			if (keyword == "DATA") {
				int numAttributes = attributesMask.size();
				if (numAttributes) {
					hasGenotypes = true;
					levelCounts.resize(numAttributes);
					levelCountsByClass.resize(numAttributes);
					attributeLevelsSeen.resize(numAttributes);

					attributeAlleleCounts.resize(numAttributes);
					attributeMinorAllele.resize(numAttributes);
					genotypeCounts.resize(numAttributes);
					attributeMutationTypes.resize(numAttributes);
				} else {
					hasGenotypes = false;
				}
				int numNumerics = numericsMask.size();
				if (numNumerics) {
					hasNumerics = true;
				} else {
					hasNumerics = false;
				}

				lineNumber = 0;
				bool makeLineIntoInstance = true;
				unsigned int instanceIndex = 0;
				int numericsAdded = 0;
				cout << Timestamp();

				while (getline(dataStream, line)) {
					++lineNumber;
					string trimmedLine = trim(line);
					// skip blank lines in the data section (usually end of file)
					if (!trimmedLine.size()) {
						continue;
					}
					// only load matching IDs, line numbers for non-plink files
					ostringstream ssLineNum;
					ssLineNum << zeroPadNumber(lineNumber, 8);
					string ID = ssLineNum.str();
					// filter out IDs
					if (!IsLoadableInstanceID(ID)) {
						cout << Timestamp() << "WARNING: " << "Dataset instance ID [" << ID
								<< "] skipped. "
								<< "Not found in list of loadable IDs. Numerics and/or "
								<< "phenotype file(s) matching filtered out this ID" << endl;
						continue;
					}

					vector<string> attributesStringVector;
					split(attributesStringVector, trimmedLine, ",");
					vector<AttributeLevel> attributeVector;
					vector<NumericLevel> numericsVector;
					unsigned int attrIdx = 0;
					unsigned int vectorIdx = 0;
					makeLineIntoInstance = true;
					vector<string>::const_iterator it = attributesStringVector.begin();
					ClassLevel discreteClassLevel = MISSING_DISCRETE_CLASS_VALUE;
					NumericLevel numericClassLevel = MISSING_NUMERIC_CLASS_VALUE;
					for (; it != attributesStringVector.end(); ++it, ++vectorIdx) {
						string thisAttr = *it;
						if (vectorIdx == classColumn) {
			        discreteClassLevel = MISSING_DISCRETE_CLASS_VALUE;
			        numericClassLevel = MISSING_NUMERIC_CLASS_VALUE;
			    		if (hasContinuousPhenotypes) {
			    			if (thisAttr != "-9") {
			    				numericClassLevel = lexical_cast<NumericLevel>(thisAttr);
			    				if (lineNumber == 1) {
			    					minPheno = maxPheno = numericClassLevel;
			    				} else {
			    					if (numericClassLevel < minPheno) {
			    						minPheno = numericClassLevel;
			    					}
			    					if (numericClassLevel > maxPheno) {
			    						maxPheno = numericClassLevel;
			    					}
			    				}
			    			} else {
			    				if (!hasAlternatePhenotypes) {
			    					cout << Timestamp() << "Instance ID " << ID
			    							<< " filtered out by missing value" << endl;
			    					continue;
			    				}
			    			}
			    		} else {
			    			if (thisAttr != "-9") {
			    				discreteClassLevel = lexical_cast<ClassLevel>(thisAttr);
			    			} else {
			    				if (!hasAlternatePhenotypes) {
			    					cout << Timestamp() << "Instance ID " << ID
			    							<< " filtered out by missing value" << endl;
			    					continue;
			    				}
			    			}
			    		}
						} else {
							if (attributeTypes[attrIdx] == ARFF_NUMERIC_TYPE) {
								if (thisAttr == "?") {
									numericsVector.push_back(MISSING_NUMERIC_VALUE);
									missingNumericValues[ID].push_back(attrIdx);
								} else {
									double thisNumericValue = lexical_cast<NumericLevel>(
											thisAttr);
									numericsVector.push_back(thisNumericValue);
								}
								++numericsAdded;
							} else {
								if (attributeTypes[attrIdx] == ARFF_NOMINAL_TYPE) {
									AttributeLevel thisAttrLevel = MISSING_ATTRIBUTE_VALUE;
									if (thisAttr == "?") {
										missingValues[ID].push_back(attrIdx);
									} else {
										thisAttrLevel = lexical_cast<AttributeLevel>(thisAttr);
										attributeLevelsSeen[attrIdx].insert(thisAttr);
									}
									attributeVector.push_back(thisAttrLevel);
									++attrIdx;
								} else {
									cout << Timestamp() << "Unrecognized attribute type!" << endl;
									return false;
								}
							}
						}
					}

					// create an instance from the vector of attribute and class values
					if (makeLineIntoInstance) {
						DatasetInstance * newInst = 0;
						if ((int) attributeVector.size() != numAttributes) {
							cerr << "ERROR: Number of attributes parsed on line "
									<< lineNumber << ": " << attributesStringVector.size()
									<< " is not equal to the number of attributes "
									<< " read from the data file header: " << numAttributes
									<< endl;
							return false;
						}
						if ((int) numericsVector.size() != numNumerics) {
							cerr << "ERROR: Number of numerics parsed on line " << lineNumber
									<< ": " << numericsVector.size()
									<< " is not equal to the number of attributes "
									<< " read from the data file header: " << numNumerics << endl;
							return false;
						}
						newInst = new DatasetInstance(this);
						if (newInst) {
							if (hasContinuousPhenotypes) {
								newInst->SetPredictedValueTau(numericClassLevel);
							} else {
								newInst->SetClass(discreteClassLevel);
								classIndexes[discreteClassLevel].push_back(instanceIndex);
							}
							if (hasGenotypes) {
								newInst->LoadInstanceFromVector(attributeVector);
							}
							if (hasNumerics) {
								for (int i = 0; i < (int) numericsVector.size(); ++i) {
									newInst->AddNumeric(numericsVector[i]);
								}
							}
							instances.push_back(newInst);
							instanceIds.push_back(ID);
							// instanceIdsToLoad.push_back(ID);
							instancesMask[ID] = instanceIndex;
						} else {
							cerr << "ERROR: loading ARFF @data section. "
									<< "Could not create dataset instance for line number "
									<< lineNumber << endl;
							return false;
						}
						++instanceIndex;
					} // make new instance

					// happy lights
					if ((lineNumber - 1) && ((lineNumber % 100) == 0)) {
						cout << lineNumber << " ";
						cout.flush();
					}
					if ((lineNumber - 1) && ((lineNumber % 1000) == 0)) {
						cout << endl << Timestamp();
					}
				} // while reading file lines

			} // keyword = data
			break;
		} // end switch

	} // end while
	cout << endl;

	dataStream.close();

	cout << Timestamp() << "There are " << NumInstances()
			<< " instances in the data set" << endl;
	cout << Timestamp() << "There are " << instancesMask.size()
			<< " instances in the instance mask" << endl;
	if (instancesMask.size() == 0) {
		cerr << "ERROR: no instances in the instance mask" << endl;
		return false;
	}

	if (hasContinuousPhenotypes) {
		continuousPhenotypeMinMax = make_pair(minPheno, maxPheno);
		cout << Timestamp() << "Continuous phenotypes." << endl;
	} else {
		cout << Timestamp() << "There are " << classIndexes.size()
				<< " classes in the data set" << endl;
	}

	if (hasNumerics) {
		// find the min and max values for each numeric attribute
		// used in diff/distance calculation metrics
		vector<NumericLevel> numericColumn;
		for (unsigned int i = 0; i < NumNumerics(); ++i) {
			GetNumericValues(i, numericColumn);
			double minElement = *numericColumn.begin();
			double maxElement = *numericColumn.begin();
			for (vector<NumericLevel>::const_iterator it = numericColumn.begin();
					it != numericColumn.end(); ++it) {
				if ((*it != MISSING_NUMERIC_VALUE) && (*it < minElement)) {
					minElement = *it;
				}
				if ((*it != MISSING_NUMERIC_VALUE) && (*it > maxElement)) {
					maxElement = *it;
				}
			}
			numericsMinMax.push_back(
					make_pair<double, double>(minElement, maxElement));
		}
	}

	if (hasGenotypes) {
		UpdateAllLevelCounts();
	}

	return true;
}

ArffAttributeType ArffDataset::GetTypeOf(unsigned int columnIndex) {
	if (columnIndex < attributeTypes.size()) {
		return attributeTypes[columnIndex];
	} else
		return (ARFF_ERROR_TYPE);
}

void ArffDataset::PrintNominalsMapping() {
	cout << Timestamp() << "Nominals and their accepted values:" << endl;
	map<string, vector<string> >::const_iterator mit = nominalValues.begin();
	for (; mit != nominalValues.end(); ++mit) {
		cout << (*mit).first << ":";
		vector<string>::const_iterator it = (*mit).second.begin();
		for (; it != (*mit).second.end(); ++it) {
			cout << " " << *it;
		}
		cout << endl;
	}
}
