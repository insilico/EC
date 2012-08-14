/*
 * AttributeRanker.cpp
 *
 *  Created on: Aug 13, 2012
 *      Author: bwhite
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include "AttributeRanker.h"
#include "Dataset.h"
#include "Insilico.h"

using namespace std;

AttributeRanker::AttributeRanker(Dataset* ds) {
	dataset = ds;
	classificationAccuracy = 1.0;
}

AttributeRanker::~AttributeRanker() {
}

AttributeScores AttributeRanker::GetScores() {
	return scores;
}

void AttributeRanker::WriteScores(string baseFilename) {
	string resultsFilename = baseFilename + ".ranks";
	ofstream outFile;
	outFile.open(resultsFilename.c_str());
	if (outFile.bad()) {
		cerr << "ERROR: Could not open scores file " << resultsFilename
				<< "for writing" << endl;
		exit(1);
	}
	PrintScores(outFile);
	outFile.close();
}

void AttributeRanker::PrintScores(ofstream& outStream) {
	for (AttributeScoresCIt scoresIt = scores.begin(); scoresIt != scores.end();
			++scoresIt) {
		// TODO: save and restore output precision
		outStream << fixed << setprecision(8)
				<< scoresIt->first << "\t"
				<< scoresIt->second << endl;
	}
}

double AttributeRanker::GetClassificationError() {
	return classificationAccuracy;
}
