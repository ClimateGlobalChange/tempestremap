///////////////////////////////////////////////////////////////////////////////
///
///	\file    VerticalInterpolate.cpp
///	\author  Paul Ullrich
///	\version November 27, 2019
///
///	<remarks>
///		Copyright 2000-2014 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#include "CommandLine.h"
#include "Exception.h"
#include "Announce.h"
#include "DataArray1D.h"
#include "DataArray2D.h"
#include "NetCDFUtilities.h"

#include "netcdfcpp.h"

#include <cfloat>
#include <string>
#include <iostream>

///////////////////////////////////////////////////////////////////////////////

void ParseDelimiter(
	const std::string & str,
	std::vector<std::string> & vecParsed,
	char szDelimiter
) {
	size_t iLast = 0;
	for (size_t i = 0; i < str.length(); i++) {
		if (str[i] == szDelimiter) {
			vecParsed.push_back(str.substr(iLast, i-iLast));
			iLast = i+1;
		}
	}
	if (iLast != str.length()) {
		vecParsed.push_back(str.substr(iLast, str.length()-iLast));
	}
}

///////////////////////////////////////////////////////////////////////////////

bool IsBetween(
	double dA,
	double dB,
	double dX
) {
	if (dA < dB) {
		if ((dX >= dA) && (dX < dB)) {
			return true;
		} else {
			return false;
		}

	} else if (dA > dB) {
		if ((dX >= dB) && (dX < dA)) {
			return true;
		} else {
			return false;
		}

	} else {
		if (dX == dA) {
			return true;
		} else {
			return false;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	// Turn off fatal errors in NetCDF
	NcError error(NcError::silent_nonfatal);

try {

	// Input filename
	std::string strInputFilenames;

	// Output filename
	std::string strOutputFilename;

	// Variable to interpolate
	std::string strVariables;

	// Old dimension variable
	std::string strOldDimensionVar;

	// New dimension variable
	std::string strNewDimensionVar;

	// Values to interpolate to
	std::string strValues;

	// First level to use for interpolation
	int iLevelBegin;

	// Last level to use for interpolation
	int iLevelLast;

	// For non-monotonic fields match the last instance
	bool fLastMatch;

	// Only interpolate points within range -- don't perform extrapolation
	bool fOnlyInRange;

	// Output the interpolating level instead
	bool fOutputInterpLevel;

	// Output the value of the vertical dimension on the interpolating level instead
	bool fOutputInterpLevelValue;

	// Check weights
	bool fCheckWeights;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFilenames, "in_files", "");
		CommandLineString(strOutputFilename, "out_file", "");
		CommandLineString(strVariables, "vars", "");
		CommandLineString(strOldDimensionVar, "olddim", "");
		CommandLineString(strNewDimensionVar, "newdim", "");
		CommandLineString(strValues, "values", "");
		CommandLineInt(iLevelBegin, "lev_begin", 0);
		CommandLineInt(iLevelLast, "lev_last", -1);
		CommandLineBool(fLastMatch, "use_last_match");
		CommandLineBool(fOnlyInRange, "only_in_range");
		CommandLineBool(fOutputInterpLevel, "output_interp_level");
		CommandLineBool(fOutputInterpLevelValue, "output_interp_level_value");
		CommandLineBool(fCheckWeights, "check_weights");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	// Output banner
	AnnounceBanner();

	// Check arguments
	if (strInputFilenames == "") {
		_EXCEPTIONT("No input filename(s) specified");
	}
	if (strOutputFilename == "") {
		_EXCEPTIONT("No output filename specified");
	}
	if (strVariables == "") {
		_EXCEPTIONT("No variable(s) specified");
	}
	if (strOldDimensionVar == "") {
		_EXCEPTIONT("No old dimension variable specified");
	}
	if (strNewDimensionVar == "") {
		_EXCEPTIONT("No new dimension variable specified");
	}
	if (strValues == "") {
		_EXCEPTIONT("No values specified");
	}
	if (fOutputInterpLevel && fOutputInterpLevelValue) {
		_EXCEPTIONT("Only one of \"--output_interp_level\" and \"--output_interp_level_value\" allowed");
	}

	// Parse input filenames
	std::vector<NcFile *> vecInputNcFiles;
	std::vector<std::string> vecInputFilenames;
	ParseDelimiter(strInputFilenames, vecInputFilenames, ';');

	// Open the input NetCDF files
	for (int f = 0; f < vecInputFilenames.size(); f++) {
		Announce("Opening file \"%s\"", vecInputFilenames[f].c_str());
		NcFile * pncInput = new NcFile(vecInputFilenames[f].c_str());
		if (!pncInput->is_valid()) {
			_EXCEPTION1("Unable to open input file \"%s\"", vecInputFilenames[f].c_str());
		}
		vecInputNcFiles.push_back(pncInput);
	}

	// Open the output NetCDF file
	NcFile ncOutput(strOutputFilename.c_str(), NcFile::Replace);
	if (!ncOutput.is_valid()) {
		_EXCEPTION1("Unable to open output file \"%s\"", strOutputFilename.c_str());
	}

	// Parse input variables
	std::vector<std::string> vecInputVariables;
	vecInputVariables.push_back(strNewDimensionVar);
	ParseDelimiter(strVariables, vecInputVariables, ',');

	// Find input variables
	std::vector<NcVar *> vecInputNcVars;
	for (int v = 0; v < vecInputVariables.size(); v++) {
		NcVar * var = NULL;
		Announce("Opening variable \"%s\"", vecInputVariables[v].c_str());
		for (int f = 0; f < vecInputNcFiles.size(); f++) {
			var = vecInputNcFiles[f]->get_var(vecInputVariables[v].c_str());
			if (var != NULL) {
				break;
			}
		}
		if (var == NULL) {
			_EXCEPTION1("Unable to find input variable \"%s\" in any input files", vecInputVariables[v].c_str());
		}
		vecInputNcVars.push_back(var);
	}

	if (vecInputVariables.size() == 0) {
		_EXCEPTIONT("Logic error");
	}

	_ASSERT(vecInputVariables.size() == vecInputNcVars.size());

	// Verify dimensionality
	Announce("Verifying consistency among variable dimensions");
	std::vector<std::string> vecDimNames;
	std::vector<int> vecDimSizes;
	int iVerticalDimIx = (-1);
	int nVarDims = vecInputNcVars[0]->num_dims();
	for (int d = 0; d < vecInputNcVars[0]->num_dims(); d++) {
		vecDimNames.push_back(vecInputNcVars[0]->get_dim(d)->name());
		vecDimSizes.push_back(vecInputNcVars[0]->get_dim(d)->size());
		if (strNewDimensionVar == vecDimNames[d]) {
			_EXCEPTION1("New dimension variable \"%s\" already exists as a dimension in input files", strNewDimensionVar.c_str());
		}
		if (strOldDimensionVar == vecDimNames[d]) {
			if (iVerticalDimIx == (-1)) {
				iVerticalDimIx = d;
			} else {
				_EXCEPTIONT("Multiple dimensions with the same name?");
			}
		}
	}

	// Number of levels along interpolating dimension
	long nVerticalDimSize = vecDimSizes[iVerticalDimIx];
	if (nVerticalDimSize < 2) {
		_EXCEPTION1("Vertical dimension \"%s\" must have size >= 2",
			strNewDimensionVar.c_str());
	}

	// Range for search along interpolating dimension
	long lLevelBegin = iLevelBegin;
	long lLevelLast = iLevelLast;
	if (lLevelLast == (-1)) {
		lLevelLast = nVerticalDimSize-1;
	}
	if (lLevelBegin < 0) {
		_EXCEPTION1("--lev_begin (%li) must be nonnegative", lLevelBegin);
	}
	if (lLevelLast < 0) {
		_EXCEPTION1("--lev_last (%li) must be nonnegative", lLevelLast);
	}
	if (lLevelLast >= nVerticalDimSize) {
		_EXCEPTION2("--lev_last (%li) out of range (must be less than "
			"interpolating dimension size %li)", lLevelLast, nVerticalDimSize);
	}
	if (lLevelLast <= lLevelBegin) {
		_EXCEPTION2("--lev_last (%li) must be greater than --lev_begin (%li)",
			lLevelLast, lLevelBegin);
	}

	// Validate spatial dimension
	int nSpatialDims = vecInputNcVars[0]->num_dims() - iVerticalDimIx - 1;
	if (nSpatialDims > 2) {
		_EXCEPTION3("No more than two spatial dimensions in variable \"%s\" "
			"permitted beyond the vertical dimension (%i..%i)",
			vecInputVariables[0].c_str(),
			iVerticalDimIx+1,
			vecInputNcVars[0]->num_dims());
	}
	if (nSpatialDims < 1) {
		_EXCEPTION3("At least one spatial dimension in variable \"%s\" "
			"required beyond the vertical dimension (%i..%i)",
			vecInputVariables[0].c_str(),
			iVerticalDimIx+1,
			vecInputNcVars[0]->num_dims());
	}

	for (int v = 1; v < vecInputVariables.size(); v++) {
		if (vecInputNcVars[v]->num_dims() != vecDimSizes.size()) {
			_EXCEPTION2("Variable \"%s\" dimensionality does not match dimension variable \"%s\"", vecInputVariables[v].c_str(), vecInputVariables[0].c_str());
		}
		for (int d = 0; d < vecDimSizes.size(); d++) {
			if (vecInputNcVars[v]->get_dim(d)->size() != vecDimSizes[d]) {
				_EXCEPTION2("Variable \"%s\" dimensionality does not match dimension variable \"%s\"", vecInputVariables[v].c_str(), vecInputVariables[0].c_str());
			}
		}
	}

	_ASSERT(iVerticalDimIx != (-1));
	_ASSERT(vecDimNames.size() == vecDimSizes.size());

	// Copy dimension variables to output file
	Announce("Copying dimension variables to output file");
	std::vector<NcDim *> vecOutputNcDims;
	vecOutputNcDims.resize(vecDimNames.size());
	for (int d = 0; d < vecDimNames.size(); d++) {
		if (vecDimNames[d] == strOldDimensionVar) {
			vecOutputNcDims[d] = NULL;
			continue;
		}
		CopyNcVarIfExists(*vecInputNcFiles[0], ncOutput, vecDimNames[d].c_str());
		NcDim * dimOut = ncOutput.get_dim(vecDimNames[d].c_str());
		if (dimOut == NULL) {
			_EXCEPTION1("Error copying dimension \"%s\" to output file",
				vecDimNames[d].c_str());
		}
		vecOutputNcDims[d] = dimOut;
	}

	// Parse values
	std::vector<std::string> vecValues;
	std::vector<double> vecValuesDouble;
	ParseDelimiter(strValues, vecValues, ',');
	int nValues = vecValues.size();
	vecValuesDouble.resize(vecValues.size());
	for (int i = 0; i < vecValues.size(); i++) {
		vecValuesDouble[i] = std::stod(vecValues[i]);
	}

	// Create new dimension
	Announce("Create new output dimension");
	NcDim * dimNewCoordinate =
		ncOutput.add_dim(
			strNewDimensionVar.c_str(),
			vecValues.size());
	if (dimNewCoordinate == NULL) {
		_EXCEPTION1("Unable to create new dimension \"%s\" in output file",
			strNewDimensionVar.c_str());
	}
	for (int d = 0; d < vecOutputNcDims.size(); d++) {
		if (vecOutputNcDims[d] == NULL) {
			vecOutputNcDims[d] = dimNewCoordinate;
		}
	}

	NcVar * varNewCoordinate =
		ncOutput.add_var(
			strNewDimensionVar.c_str(),
			ncDouble,
			dimNewCoordinate);
	if (varNewCoordinate == NULL) {
		_EXCEPTION1("Unable to create new variable \"%s\" in output file",
			strNewDimensionVar.c_str());
	}

	varNewCoordinate->put(&(vecValuesDouble[0]), vecValues.size());
	CopyNcVarAttributes(vecInputNcVars[0], varNewCoordinate);

	// Create new variables
	Announce("Creating output variables");
	std::vector<NcVar *> vecOutputNcVars;
	vecOutputNcVars.resize(vecInputNcVars.size());
	for (int v = 1; v < vecInputNcVars.size(); v++) {
		vecOutputNcVars[v] =
			ncOutput.add_var(
				vecInputVariables[v].c_str(),
				ncDouble,
				vecOutputNcDims.size(),
				const_cast<const NcDim **>(&(vecOutputNcDims[0])));

		if (vecOutputNcVars[v] == NULL) {
			_EXCEPTION1("Error creating output variable \"%s\" in output file",
				vecInputVariables[v].c_str());
		}

		CopyNcVarAttributes(vecInputNcVars[v], vecOutputNcVars[v]);
	}

	// Perform interpolation
	AnnounceStartBlock("Performing interpolation");
	
	// Single column sizes
	std::vector<long> vecVarColumnSize;
	vecVarColumnSize.resize(nVarDims);
	for (int d = 0; d < nVarDims; d++) {
		vecVarColumnSize[d] = 1;
	}
	vecVarColumnSize[iVerticalDimIx] = nVerticalDimSize;

	std::vector<long> vecVarColumnPos;
	vecVarColumnPos.resize(nVarDims);

	// Single horizontal slice sizes
	std::vector<long> vecVarSliceSize;
	vecVarSliceSize.resize(nVarDims);
	for (int d = 0; d <= iVerticalDimIx; d++) {
		vecVarSliceSize[d] = 1;
	}

	long nSpatialDOFs = 1;
	for (int d = iVerticalDimIx+1 ; d < vecDimSizes.size(); d++) {
		nSpatialDOFs *= vecDimSizes[d];
		vecVarSliceSize[d] = vecDimSizes[d];
	}

	// Allocate data
	DataArray2D<int> nLevelIndexLower(nSpatialDOFs, nValues);
	DataArray2D<double> dWeightLower(nSpatialDOFs, nValues);

	DataArray1D<double> dSliceData1(nSpatialDOFs);
	DataArray1D<double> dSliceData2(nSpatialDOFs);
	DataArray1D<double> * pSliceDataLower = &dSliceData2;
	DataArray1D<double> * pSliceDataUpper = &dSliceData1;

	DataArray2D<double> dOutputData(nValues, nSpatialDOFs);

	// Load old dimension (only used when --output_interp_level_value)
	DataArray1D<double> dOldDimValues;
	if (fOutputInterpLevelValue) {
		_ASSERT(vecInputNcFiles.size() > 0);
		NcVar * varOldDim = vecInputNcFiles[0]->get_var(strOldDimensionVar.c_str());
		if (varOldDim == NULL) {
			_EXCEPTION1("Unable to find dimension variable \"%s\" in input files",
				strOldDimensionVar.c_str());
		}
		if (varOldDim->num_dims() != 1) {
			_EXCEPTION1("Dimension variable \"%s\" has dimensionality greater than 1",
				strOldDimensionVar.c_str());
		}
		dOldDimValues.Allocate(varOldDim->get_dim(0)->size());
		varOldDim->set_cur((long)0);
		varOldDim->get(&(dOldDimValues[0]), dOldDimValues.GetRows());
	}

	// Loop through auxiliary dimensions
	std::vector<long> vecVarSlicePos;
	vecVarSlicePos.resize(nVarDims);

	int nAuxOptions = 1;
	for (int d = 0; d < iVerticalDimIx; d++) {
		nAuxOptions *= vecDimSizes[d];
	}
	for (int a = 0; a < nAuxOptions; a++) {
	//for (int a = 0; a < 1; a++) {

		// Build the auxiliary index to this slice
		std::string strAuxIndex;
		int ax = a;
		for (int d = iVerticalDimIx-1; d >= 0; d--) {
			vecVarSlicePos[d] = ax % vecDimSizes[d];
			ax /= vecDimSizes[d];
			strAuxIndex = std::to_string((long long) vecVarSlicePos[d]) + "," + strAuxIndex;
		}

		_ASSERT((nSpatialDims == 1) || (nSpatialDims == 2));

		std::string strSpatialRange;
		if (nSpatialDims == 1) {
			strSpatialRange = ":";
		} else if (nSpatialDims == 2) {
			strSpatialRange = ":,:";
		}

		// Check for fillvalue
		double dFillValue = DBL_MAX;
		NcAtt * attFillValue = vecInputNcVars[0]->get_att("_FillValue");
		if (attFillValue != NULL) {
			dFillValue = attFillValue->as_float(0);
		}

		///////////////////////////////////////////////////////////
		// Build the level index
		AnnounceStartBlock("Building level index %s(%s:,%s)",
			vecInputVariables[0].c_str(),
			strAuxIndex.c_str(),
			strSpatialRange.c_str());

		nLevelIndexLower.Zero();

		// Number of spatial points where the interpolation point has been found
		long lDone = 0;

		// Default level index to "not found"
		for (int i = 0; i < nSpatialDOFs; i++) {
			for (int k = 0; k < nValues; k++) {
				nLevelIndexLower[i][k] = (-1);
			}
		}

		for (long l = lLevelBegin; l <= lLevelLast; l++) {

			// Load slice data
			vecVarSlicePos[iVerticalDimIx] = l;
			//for (int d = 0; d < nVarDims; d++) {
			//	printf("%li %li;", vecVarSlicePos[d], vecVarSliceSize[d]);
			//}
			//printf("\n");

			vecInputNcVars[0]->set_cur(&(vecVarSlicePos[0]));
			vecInputNcVars[0]->get(
				&((*pSliceDataUpper)[0]),
				&(vecVarSliceSize[0]));

			//printf("%li %1.15e\n", l, (*pSliceDataUpper)[174477]);

			// Is vecValuesDouble[k] below the lowermost level value?
			if ((!fOnlyInRange) && (l == lLevelBegin+1)) {
				for (int i = 0; i < nSpatialDOFs; i++) {
					const double dLower = (*pSliceDataLower)[i];
					const double dUpper = (*pSliceDataUpper)[i];

					if ((dLower == dFillValue) || (dUpper == dFillValue)) {
						continue;
					}

					// Field increasing with level
					if (dUpper > dLower) {
						for (int k = 0; k < nValues; k++) {
							if (vecValuesDouble[k] <= (*pSliceDataLower)[i]) {
								if (!fLastMatch && (nLevelIndexLower[i][k] != (-1))) {
									continue;
								}
								nLevelIndexLower[i][k] = 0;
								dWeightLower[i][k] = 1.0;
								lDone++;
							}
						}

					// Field decreasing with level
					} else {
						for (int k = 0; k < nValues; k++) {
							if (vecValuesDouble[k] >= (*pSliceDataLower)[i]) {
								if (!fLastMatch && (nLevelIndexLower[i][k] != (-1))) {
									continue;
								}
								nLevelIndexLower[i][k] = 0;
								dWeightLower[i][k] = 1.0;
								lDone++;
							}
						}
					}
				}
			}

			// Is vecValuesDouble[k] above the uppermost level value?
			if ((!fOnlyInRange) && (l == lLevelLast)) {
				for (int i = 0; i < nSpatialDOFs; i++) {
					const double dLower = (*pSliceDataLower)[i];
					const double dUpper = (*pSliceDataUpper)[i];

					if ((dLower == dFillValue) || (dUpper == dFillValue)) {
						continue;
					}

					// Field increasing with level
					if (dUpper > dLower) {
						for (int k = 0; k < nValues; k++) {
							if (vecValuesDouble[k] >= (*pSliceDataUpper)[i]) {
								if (!fLastMatch && (nLevelIndexLower[i][k] != (-1))) {
									continue;
								}
								nLevelIndexLower[i][k] = (int)(nVerticalDimSize-2);
								dWeightLower[i][k] = 0.0;
								lDone++;
							}
						}

					// Field decreasing with level
					} else {
						for (int k = 0; k < nValues; k++) {
							if (vecValuesDouble[k] <= (*pSliceDataUpper)[i]) {
								if (!fLastMatch && (nLevelIndexLower[i][k] != (-1))) {
									continue;
								}
								nLevelIndexLower[i][k] = (int)(nVerticalDimSize-2);
								dWeightLower[i][k] = 0.0;
								lDone++;
							}
						}
					}
				}
			}

			// Is vecValuesDouble[k] somewhere in between?
			if (l != lLevelBegin) {
				for (int i = 0; i < nSpatialDOFs; i++) {
					const double dLower = (*pSliceDataLower)[i];
					const double dUpper = (*pSliceDataUpper)[i];

					if ((dLower == dFillValue) || (dUpper == dFillValue)) {
						continue;
					}

					// Field increasing with level
					for (int k = 0; k < nValues; k++) {
						if (IsBetween(dLower, dUpper, vecValuesDouble[k])) {
							if (!fLastMatch && (nLevelIndexLower[i][k] != (-1))) {
								continue;
							}
							//if (fabs(dUpper - dLower) < 1.0e-12) {
							//	continue;
							//}
							nLevelIndexLower[i][k] = (int)(l);
							dWeightLower[i][k] =
								(dUpper - vecValuesDouble[k]) / (dUpper - dLower);
							lDone++;
						}
					}
				}
			}

			Announce("%s(%s%li/%li,%s) (%li/%li)",
				vecInputVariables[0].c_str(),
				strAuxIndex.c_str(),
				l, nVerticalDimSize,
				strSpatialRange.c_str(),
				lDone, nValues * nSpatialDOFs);

			// Swap arrays for next load
			DataArray1D<double> * pSliceDataTemp = pSliceDataUpper;
			pSliceDataUpper = pSliceDataLower;
			pSliceDataLower = pSliceDataTemp;
		}

		if (lDone > nValues * nSpatialDOFs) {
			Announce("WARNING: Coordinate field is not monotonic in the vertical (%li/%li)",
				lDone, nSpatialDOFs);
		}
		AnnounceEndBlock("Done");

		///////////////////////////////////////////////////////////
		// Checking weights
		if (fCheckWeights) {
			for (int i = 0; i < nSpatialDOFs; i++) {
				for (int k = 0; k < nValues; k++) {
					if (nLevelIndexLower[i][k] != -1) {
						if ((dWeightLower[i][k] < 0.0) || (dWeightLower[i][k] > 1.0)) {
							_EXCEPTION2("Weight out of range (%i %i)", i, k);
						}
					}
				}
			}
		}

		///////////////////////////////////////////////////////////
		// Linearly interpolate using level index
		for (int v = 1; v < vecInputNcVars.size(); v++) {
			AnnounceStartBlock("Interpolating %s(%s:,%s)",
				vecInputVariables[v].c_str(),
				strAuxIndex.c_str(),
				strSpatialRange.c_str());

			NcAtt * attFillValue = vecOutputNcVars[v]->get_att("_FillValue");
			if (attFillValue != NULL) {
				double dFillValue = attFillValue->as_double(0);
				for (int i = 0; i < nSpatialDOFs; i++) {
					for (int k = 0; k < nValues; k++) {
						dOutputData[k][i] = dFillValue;
					}
				}
			} else {
				dOutputData.Zero();
			}
			if (fOutputInterpLevel) {
				for (int i = 0; i < nSpatialDOFs; i++) {
					for (int k = 0; k < nValues; k++) {
						if (nLevelIndexLower[i][k] != -1) {
							dOutputData[k][i] = (double) nLevelIndexLower[i][k];
						}
					}
				}

			} else if (fOutputInterpLevelValue) {
				for (int i = 0; i < nSpatialDOFs; i++) {
					for (int k = 0; k < nValues; k++) {
						if (nLevelIndexLower[i][k] != -1) {
							const double dWeight = dWeightLower[i][k];
							dOutputData[k][i] =
								dWeight * dOldDimValues[nLevelIndexLower[i][k]]
								+ (1.0 - dWeight) * dOldDimValues[nLevelIndexLower[i][k]+1];
						}
					}
				}

			} else {

				for (long l = lLevelBegin; l <= lLevelLast; l++) {

					// Load slice data
					vecVarSlicePos[iVerticalDimIx] = l;
					//for (int d = 0; d < nVarDims; d++) {
					//	printf("%li %li;", vecVarSlicePos[d], vecVarSliceSize[d]);
					//}
					//printf("\n");

					Announce("%s(%s%li/%li,%s)",
						vecInputVariables[v].c_str(),
						strAuxIndex.c_str(),
						l, nVerticalDimSize,
						strSpatialRange.c_str());

					vecInputNcVars[v]->set_cur(&(vecVarSlicePos[0]));
					vecInputNcVars[v]->get(
						&((*pSliceDataUpper)[0]),
						&(vecVarSliceSize[0]));

					if (l != lLevelBegin) {
						for (int i = 0; i < nSpatialDOFs; i++) {
							for (int k = 0; k < nValues; k++) {
								if (nLevelIndexLower[i][k] == l-1) {
									const double dWeight = dWeightLower[i][k];
									if ((dWeight < 0.0) || (dWeight > 1.0)) {
										printf("WARNING: %1.5e\n", dWeight);
									}
									dOutputData[k][i] =
										dWeight * (*pSliceDataLower)[i]
										+ (1.0 - dWeight) * (*pSliceDataUpper)[i];
								}
							}
						}
					}

					// Swap arrays for next load
					DataArray1D<double> * pSliceDataTemp = pSliceDataUpper;
					pSliceDataUpper = pSliceDataLower;
					pSliceDataLower = pSliceDataTemp;
				}
			}

			AnnounceEndBlock("Done");

			// Write results
			AnnounceStartBlock("Writing results");
			std::vector<long> vecVarSliceOutputSizes = vecVarSliceSize;
			vecVarSliceOutputSizes[iVerticalDimIx] = nValues;
			vecOutputNcVars[v]->set_cur(&(vecVarSlicePos[0]));
			vecOutputNcVars[v]->put(
				&(dOutputData[0][0]),
				&(vecVarSliceOutputSizes[0]));
			AnnounceEndBlock("Done");

		}

	}

	AnnounceEndBlock(NULL);

	AnnounceBanner();

} catch(Exception & e) {
	std::cout << e.ToString().c_str() << std::endl;
	return (-1);

} catch(...) {
	return (-2);
}
}

