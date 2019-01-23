///////////////////////////////////////////////////////////////////////////////
///
///	\file    GenerateRLLMesh.cpp
///	\author  Paul Ullrich
///	\version March 7, 2014
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
#include "SparseMatrix.h"
#include "Exception.h"
#include "Announce.h"
#include "DataArray2D.h"
#include "NetCDFUtilities.h"

#include <cmath>
#include <iostream>

#include "netcdfcpp.h"

///////////////////////////////////////////////////////////////////////////////

typedef std::pair<int, double> IndexCoefficientPair;

typedef std::vector< std::vector<IndexCoefficientPair> > ConservativeMap1D;

///////////////////////////////////////////////////////////////////////////////

void Generate1DOverlapMap(
	bool fFlipped,
	const DataArray1D<double> dXedge,
	const DataArray1D<double> dXedgeout,
	ConservativeMap1D & matX1D
) {

	const int nX = dXedge.GetRows() - 1;
	const int nXout = dXedgeout.GetRows() - 1;

	static double Tol = 1.0e-12;

	matX1D.clear();
	matX1D.resize(nXout);

	int io = 0;
	int ii = 0;
	for (;;) {

		double dXi0 = dXedge[ii];
		double dXi1 = dXedge[ii+1];

		double dDeltaXi = dXi1 - dXi0;
		double dDeltaXo = dXedgeout[io+1] - dXedgeout[io];

		bool fLeftIn;
		bool fRightIn;

		if (!fFlipped) {

			if ((ii == 0) && (dXedgeout[io+1] < dXi0)) {
				io++;
				continue;
			}
			if ((ii == nX-1) && (dXedgeout[io] > dXi1)) {
				break;
			}

			fLeftIn = (dXedgeout[io] > dXi0 - Tol);
			fRightIn = (dXedgeout[io+1] < dXi1 + Tol);

		} else {

			if ((ii == 0) && (dXedgeout[io+1] > dXi0)) {
				io++;
				continue;
			}
			if ((ii == nX-1) && (dXedgeout[io] < dXi1)) {
				break;
			}

			fLeftIn = (dXedgeout[io] < dXi0 + Tol);
			fRightIn = (dXedgeout[io+1] > dXi1 - Tol);
		}

		// Output volume fully contained in input
		if (fLeftIn && fRightIn) {
			matX1D[io].push_back(
				IndexCoefficientPair(
					ii, (dXedgeout[io+1] - dXedgeout[io]) / dDeltaXi * fabs(dDeltaXi / dDeltaXo)));

			if (fabs(dXedgeout[io+1] - dXedge[ii+1]) < Tol) {
				ii++;
			}
			io++;

		// Right edge of output volume outside input
		} else if ((fLeftIn) && (!fRightIn)) {
			matX1D[io].push_back(
				IndexCoefficientPair(
					ii, (dXedge[ii+1] - dXedgeout[io]) / dDeltaXi * fabs(dDeltaXi / dDeltaXo)));
			ii++;

		// Left edge of output volume outside input
		} else if ((!fLeftIn) && (fRightIn)) {
			matX1D[io].push_back(
				IndexCoefficientPair(
					ii, (dXedgeout[io+1] - dXedge[ii]) / dDeltaXi * fabs(dDeltaXi / dDeltaXo)));

			if (fabs(dXedgeout[io+1] - dXedge[ii+1]) < Tol) {
				ii++;
			}
			io++;

		// Input volume fully contained in output
		} else {
			matX1D[io].push_back(
				IndexCoefficientPair(
					ii, fabs(dDeltaXi / dDeltaXo)));
			ii++;
		} 

		// Done looking at output volumes
		if (io >= nXout) {
			break;
		}

		// Done looking at input volumes
		if (ii >= nX) {
			break;
		}
	}

	// Validate the map
	DataArray1D<double> dTotals(nX);
	for (io = 0; io < nXout; io++) {
		double dDeltaXo = fabs(dXedgeout[io+1] - dXedgeout[io]);
		for (int i = 0; i < matX1D[io].size(); i++) {
			if (matX1D[io][i].second < 0.0) {
				_EXCEPTION1("Out of range (%1.15e)", matX1D[io][i].second);
			}
			if (matX1D[io][i].second > 1.0) {
				_EXCEPTION1("Out of range (%1.15e)", matX1D[io][i].second);
			}
			dTotals[matX1D[io][i].first] +=
				matX1D[io][i].second * dDeltaXo;
		}
	}
	for (int i = 0; i < nX; i++) {
		double dDeltaXi = fabs(dXedge[i+1] - dXedge[i]);
		if (fabs(dTotals[i] - dDeltaXi) > Tol) {
			_EXCEPTION2("Map did not validate (%1.15e : %1.15e)",
				dTotals[i], dDeltaXi);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

template<typename T>
inline void FindReplace(
	DataArray2D<T> & data,
	const T & find,
	const T & replace
) {
	const int nX = data.GetRows();
	const int nY = data.GetColumns();

	for (int i = 0; i < nX; i++) {
	for (int j = 0; j < nY; j++) {
		if (data[i][j] == find) {
			data[i][j] = replace;
		}
	}
	}
}

///////////////////////////////////////////////////////////////////////////////

template<typename T>
inline double TotalMatrix(
	const DataArray2D<T> & data,
	const DataArray1D<double> & dXEdge,
	const DataArray1D<double> & dYEdge
) {
	double dTotal = 0.0;
	for (int i = 0; i < data.GetRows(); i++) {
	for (int j = 0; j < data.GetColumns(); j++) {
		dTotal += static_cast<double>(data[i][j])
			* fabs(dXEdge[i+1] - dXEdge[i])
			* fabs(dYEdge[j+1] - dYEdge[j]);
	}
	}

	return dTotal;
}

///////////////////////////////////////////////////////////////////////////////

template<typename T, typename Tout = T>
inline void ApplyMap(
	const DataArray2D<T> & dataIn,
	const ConservativeMap1D & mapX,
	const ConservativeMap1D & mapY,
	DataArray2D<Tout> & dataOut
) {
	const int nX = dataIn.GetRows();
	const int nY = dataIn.GetColumns();

	const int nXout = dataOut.GetRows();
	const int nYout = dataOut.GetColumns();

	for (int i = 0; i < nXout; i++) {
	for (int j = 0; j < nYout; j++) {

		double dOut = 0.0;

		for (int ix = 0; ix < mapX[i].size(); ix++) {
		for (int iy = 0; iy < mapY[j].size(); iy++) {

			const int ii = mapX[i][ix].first;
			const int jj = mapY[j][iy].first;

			const double dCoeffX = mapX[i][ix].second;
			const double dCoeffY = mapY[j][iy].second;
/*
			if (ii >= dataIn.GetRows()) {
				_EXCEPTION();
			}
			if (jj >= dataIn.GetColumns()) {
				_EXCEPTION();
			}
*/
			dOut += dCoeffX * dCoeffY * static_cast<double>(dataIn[ii][jj]);
		}
		}

		dataOut[i][j] = static_cast<Tout>(dOut);
	}
	}
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	NcError error(NcError::silent_nonfatal);

try {
	// Number of x volumes in output mesh
	int nXout;

	// Number of y volumes in output mesh
	int nYout;

	// First x line on mesh
	double dXBegin = 0.0;

	// Last x line on mesh
	double dXEnd = 0.0;

	// First y line on mesh
	double dYBegin = 0.0;

	// Last y line on mesh
	double dYEnd = 0.0;

	// Output filename
	std::string strInputFile;

	// Output filename
	std::string strOutputFile;

	// Variable to remap
	std::string strVariable;

	// Name of x dimension
	std::string strVarNameX;

	// Name of y dimension
	std::string strVarNameY;

	// Find missing data values with this value
	std::string strFind;

	// Replace missing data values with this value
	std::string strReplace;

	// Validate conservation
	bool fValidate;

	// Output double
	bool fOutputDouble;

	// Parse the command line
	BeginCommandLine()
		CommandLineInt(nXout, "nx", 64);
		CommandLineInt(nYout, "ny", 64);
		//CommandLineDouble(dXBegin, "xbegin", 0.0);
		//CommandLineDouble(dXEnd, "xend", 0.0);
		//CommandLineDouble(dYBegin, "ybegin", 0.0);
		//CommandLineDouble(dYEnd, "yend", 0.0);
		CommandLineString(strInputFile, "in", "");
		CommandLineString(strOutputFile, "out", "");
		CommandLineString(strVariable, "var", "");
		CommandLineString(strVarNameX, "xvarname", "lon");
		CommandLineString(strVarNameY, "yvarname", "lat");
		CommandLineString(strFind, "find", "");
		CommandLineString(strReplace, "replace", "");
		CommandLineBool(fValidate, "validate");
		CommandLineBool(fOutputDouble, "doubleout");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	// Check command line arguments
	if (strInputFile == "") {
		_EXCEPTIONT("No input file specified");
	}
	if (strOutputFile == "") {
		_EXCEPTIONT("No output file specified");
	}
	if (strVariable == "") {
		_EXCEPTIONT("No data variable specified");
	}
	if (nXout <= 0) {
		_EXCEPTIONT("At least one X output volume required");
	}
	if (nYout <= 0) {
		_EXCEPTIONT("At least one Y output volume required");
	}

	// Load NetCDF input file
	NcFile ncfilein(strInputFile.c_str(), NcFile::ReadOnly);
	if (!ncfilein.is_valid()) {
		_EXCEPTIONT("Unable to open input file");
	}

	// Load NetCDF output file
	NcFile ncfileout(strOutputFile.c_str(), NcFile::Replace);
	if (!ncfileout.is_valid()) {
		_EXCEPTIONT("Unable to open output file");
	}

	// Get dimensions of input data
	NcVar * varX = ncfilein.get_var(strVarNameX.c_str());
	if (varX == NULL) {
		_EXCEPTION1("Unable to find variable \"%s\" in input file",
			strVarNameX.c_str());
	}
	if ((varX->type() != ncDouble) && (varX->type() != ncFloat)) {
		_EXCEPTIONT("Only double or float type supported for dimension variables");
	}

	NcVar * varY = ncfilein.get_var(strVarNameY.c_str());
	if (varY == NULL) {
		_EXCEPTION1("Unable to find variable \"%s\" in input file",
			strVarNameY.c_str());
	}
	if ((varY->type() != ncDouble) && (varY->type() != ncFloat)) {
		_EXCEPTIONT("Only double or float type supported for dimension variables");
	}

	// Load variable
	NcVar * varData = ncfilein.get_var(strVariable.c_str());
	if (varData == NULL) {
		_EXCEPTION1("Unable to find variable \"%s\" in input file",
			strVariable.c_str());
	}

	// Identify relevant dimensions
	AnnounceBanner();
	AnnounceStartBlock("Identifying relevant dimensions");

	int iVarDims = varData->num_dims();
	if (iVarDims < 2) {
		_EXCEPTIONT("Input variable must have at least 2 dimensions");
	}
	NcDim * dimX = varData->get_dim(iVarDims-2);
	NcDim * dimY = varData->get_dim(iVarDims-1);

	if (dimX->size() <= 1) {
		_EXCEPTION1("Dimension %s must have size >= 1 (storage order?)",
			dimX->name());
	}
	if (dimY->size() <= 1) {
		_EXCEPTION1("Dimension %s must have size >= 1 (storage order?)",
			dimY->name());
	}

	Announce("Data Dimension X: %s (%lu)", dimX->name(), dimX->size());
	Announce("Data Dimension Y: %s (%lu)", dimY->name(), dimY->size());

	int iDimX = 0;
	int iDimY = 0;

	DataArray1D<long> lSizeX(varX->num_dims());
	DataArray1D<long> lSizeY(varY->num_dims());

	if (varX->num_dims() != 1) {
		Announce("Multiple dimensions detected in X");
		iDimX = (-1);
		for (int d = 0; d < varX->num_dims(); d++) {
			NcDim * dim = varX->get_dim(d);
			if (dim == dimX) {
				Announce("Match found for %s (%i)", dim->name(), d);
				iDimX = d;
				break;
			}
		}
		if (iDimX == (-1)) {
			_EXCEPTION1("No matching dimension found in %s",
				strVarNameX.c_str());
		}
		for (int d = 0; d < varX->num_dims(); d++) {
			lSizeX[d] = 1;
		}
	}
	long nX = dimX->size();
	lSizeX[iDimX] = dimX->size();

	if (varY->num_dims() != 1) {
		Announce("Multiple dimensions detected in Y");
		iDimY = (-1);
		for (int d = 0; d < varY->num_dims(); d++) {
			NcDim * dim = varY->get_dim(d);
			if (dim == dimY) {
				Announce("Match found for %s (%i)", dim->name(), d);
				iDimY = d;
				break;
			}
		}
		if (iDimY == (-1)) {
			_EXCEPTION1("No matching dimension found in %s",
				strVarNameY.c_str());
		}
		for (int d = 0; d < varY->num_dims(); d++) {
			lSizeY[d] = 1;
		}
	}
	long nY = dimY->size();
	lSizeY[iDimY] = dimY->size();

	AnnounceEndBlock("Done");

	// Load in data from dimension
	AnnounceStartBlock("Load dimension data");

	DataArray1D<double> dX(dimX->size());
	DataArray1D<double> dY(dimY->size());

	varX->get(dX, (long*)(lSizeX));
	varY->get(dY, (long*)(lSizeY));

	bool fFlippedX = false;
	bool fFlippedY = false;

	if (dX[1] < dX[0]) {
		fFlippedX = true;
		Announce("Flipped X dimension detected");
	}
	if (dY[1] < dY[0]) {
		fFlippedY = true;
		Announce("Flipped Y dimension detected");
	}

	// Identify volume boundaries
	Announce("Identifying volume boundaries");
	DataArray1D<double> dXedge(nX+1);
	DataArray1D<double> dYedge(nY+1);

	dXedge[0] = dX[0] - 0.5 * (dX[1] - dX[0]);
	for (int i = 1; i < nX; i++) {
		dXedge[i] = 0.5 * (dX[i-1] + dX[i]);
	}
	dXedge[nX] = dX[nX-1] + 0.5 * (dX[nX-1] - dX[nX-2]);

	dYedge[0] = dY[0] - 0.5 * (dY[1] - dY[0]);
	for (int j = 1; j < nY; j++) {
		dYedge[j] = 0.5 * (dY[j-1] + dY[j]);
	}
	dYedge[nY] = dY[nY-1] + 0.5 * (dY[nY-1] - dY[nY-2]);

	Announce("X Bounds: [%1.5f, %1.5f]", dXedge[0], dXedge[nX]);
	Announce("Y Bounds: [%1.5f, %1.5f]", dYedge[0], dYedge[nY]);

	AnnounceEndBlock("Done");

	// Generating output mesh
	AnnounceStartBlock("Generating output mesh");

	DataArray1D<double> dXout(nXout);
	DataArray1D<double> dYout(nYout);

	DataArray1D<double> dXedgeout(nXout+1);
	DataArray1D<double> dYedgeout(nYout+1);

	if (dXEnd == dXBegin) {
		dXBegin = dXedge[0];
		dXEnd = dXedge[nX];
	}

	if (dYEnd == dYBegin) {
		dYBegin = dYedge[0];
		dYEnd = dYedge[nY];
	}

	double dDeltaX = (dXEnd - dXBegin) / static_cast<double>(nXout);
	double dDeltaY = (dYEnd - dYBegin) / static_cast<double>(nYout);

	bool fFlippedXout = (dDeltaX < 0.0);
	bool fFlippedYout = (dDeltaY < 0.0);

	for (int i = 0; i <= nXout; i++) {
		dXedgeout[i] = dXBegin + dDeltaX * static_cast<double>(i);
	}
	for (int i = 0; i < nXout; i++) {
		dXout[i] = 0.5 * (dXedgeout[i] + dXedgeout[i+1]);
	}
	for (int j = 0; j <= nYout; j++) {
		dYedgeout[j] = dYBegin + dDeltaY * static_cast<double>(j);
	}
	for (int j = 0; j < nYout; j++) {
		dYout[j] = 0.5 * (dYedgeout[j] + dYedgeout[j+1]);
	}

	NcDim * dimXout = ncfileout.add_dim(strVarNameX.c_str(), nXout);
	NcDim * dimYout = ncfileout.add_dim(strVarNameY.c_str(), nYout);

	NcVar * varXout = ncfileout.add_var(strVarNameX.c_str(), ncDouble, dimXout);
	varXout->put(dXout, nXout);
	CopyNcVarAttributes(varX, varXout);

	NcVar * varYout = ncfileout.add_var(strVarNameY.c_str(), ncDouble, dimYout);
	varYout->put(dYout, nYout);
	CopyNcVarAttributes(varY, varYout);

	// Data size of one input instance
	DataArray1D<long> lDataSize(varData->num_dims());
	int nData = 1;
	for (int d = 0; d < varData->num_dims()-2; d++) {
		nData *= varData->get_dim(d)->size();
		lDataSize[d] = 1;
	}
	lDataSize[varData->num_dims()-2] = nX;
	lDataSize[varData->num_dims()-1] = nY;

	// Data size of one output instance
	DataArray1D<long> lDataSizeOut = lDataSize;
	lDataSizeOut[varData->num_dims()-2] = nXout;
	lDataSizeOut[varData->num_dims()-1] = nYout;

	Announce("%i data instance(s) found", nData);

	// Pointer to output dimensions
	DataArray1D<NcDim*> vecDimOut(varData->num_dims());
	for (int d = 0; d < varData->num_dims()-2; d++) {
		vecDimOut[d] =
			ncfileout.add_dim(
				varData->get_dim(d)->name(),
				varData->get_dim(d)->size());
	}
	vecDimOut[varData->num_dims()-2] = dimXout;
	vecDimOut[varData->num_dims()-1] = dimYout;

	const NcDim ** pvardims = (const NcDim **)(&(vecDimOut[0]));

	AnnounceEndBlock("Done");

	// Generate 1D maps
	AnnounceStartBlock("Generating 1D maps");

	ConservativeMap1D matX1D;
	ConservativeMap1D matY1D;

	Generate1DOverlapMap(fFlippedX, dXedge, dXedgeout, matX1D);
	Generate1DOverlapMap(fFlippedY, dYedge, dYedgeout, matY1D);

	AnnounceEndBlock("Done");

	// Load data
	AnnounceStartBlock("Begin remapping");

	DataArray2D<char> dByteData;
	DataArray2D<short> dShortData;
	DataArray2D<int> dIntData;
	DataArray2D<float> dFloatData;
	DataArray2D<double> dDoubleData;

	DataArray2D<char> dByteDataOut;
	DataArray2D<short> dShortDataOut;
	DataArray2D<int> dIntDataOut;
	DataArray2D<float> dFloatDataOut;
	DataArray2D<double> dDoubleDataOut;

	char bFindByte = 0;
	char bReplaceByte = 0;

	short sFindShort = 0;
	short sReplaceShort = 0;

	int iFindInt = 0;
	int iReplaceInt = 0;

	float dFindFloat = 0.0f;
	float dReplaceFloat = 0.0f;

	double dFindDouble = 0.0;
	double dReplaceDouble = 0.0;

	NcVar * varDataOut = NULL;

	if (varData->type() == ncByte) {
		dByteData.Allocate(nX, nY);
		dByteDataOut.Allocate(nXout, nYout);

		if (strFind != "") {
			bFindByte = (char)atoi(strFind.c_str());
		}
		if (strReplace != "") {
			bReplaceByte = (char)atoi(strReplace.c_str());
		}

	} else if (varData->type() == ncShort) {
		dShortData.Allocate(nX, nY);
		dShortDataOut.Allocate(nXout, nYout);

		if (strFind != "") {
			sFindShort = (short)atoi(strFind.c_str());
		}
		if (strReplace != "") {
			sReplaceShort = (short)atoi(strReplace.c_str());
		}

	} else if (varData->type() == ncInt) {
		dIntData.Allocate(nX, nY);
		dIntDataOut.Allocate(nXout, nYout);

		if (strFind != "") {
			iFindInt = (int)atoi(strFind.c_str());
		}
		if (strReplace != "") {
			iReplaceInt = (int)atoi(strReplace.c_str());
		}

	} else if (varData->type() == ncFloat) {
		dFloatData.Allocate(nX, nY);
		dFloatDataOut.Allocate(nXout, nYout);

		if (strFind != "") {
			dFindFloat = (float)atof(strFind.c_str());
		}
		if (strReplace != "") {
			dReplaceFloat = (float)atof(strReplace.c_str());
		}

	} else if (varData->type() == ncDouble) {
		dDoubleData.Allocate(nX, nY);
		dDoubleDataOut.Allocate(nXout, nYout);

		if (strFind != "") {
			dFindDouble = (double)atof(strFind.c_str());
		}
		if (strReplace != "") {
			dReplaceDouble = (double)atof(strReplace.c_str());
		}

	} else {
		_EXCEPTIONT("Invalid datatype");
	}

	if (fOutputDouble) {
		dDoubleDataOut.Allocate(nXout, nYout);
	}

	NcType vartype = varData->type();
	if (fOutputDouble) {
		vartype = ncDouble;
	}
	varDataOut = ncfileout.add_var(strVariable.c_str(), vartype, vecDimOut.GetRows(), pvardims);
	if (varDataOut == NULL) {
		_EXCEPTIONT("Unable to create output variable");
	}

	// Copy variable attributes
	CopyNcVarAttributes(varData, varDataOut);

	// Remap
	for (int i = 0; i < nData; i++) {
		char szBuffer[256];
		sprintf(szBuffer, "Remapping data instance %i", i);
		AnnounceStartBlock(szBuffer);

		DataArray1D<long> lDataIx(varData->num_dims());
		lDataIx[varData->num_dims()-1] = 0;
		lDataIx[varData->num_dims()-2] = 0;

		int iData = i;
		for (int d = varData->num_dims()-2; d >= 0; d--) {
			lDataIx[d] = iData % varData->get_dim(d)->size();
			iData /= varData->get_dim(d)->size();
		}

		varData->set_cur((long*)(lDataIx));
		varDataOut->set_cur((long*)(lDataIx));

		double dTotalMassIn = 0.0;
		double dTotalMassOut = 0.0;

		// Read data
		{
			Announce("Reading data");
			if (varData->type() == ncByte) {
				varData->get(&(dByteData[0][0]), lDataSize);
			} else if (varData->type() == ncShort) {
				varData->get(&(dShortData[0][0]), lDataSize);
			} else if (varData->type() == ncInt) {
				varData->get(&(dIntData[0][0]), lDataSize);
			} else if (varData->type() == ncFloat) {
				varData->get(&(dFloatData[0][0]), lDataSize);
			} else if (varData->type() == ncDouble) {
				varData->get(&(dDoubleData[0][0]), lDataSize);
			} else {
				_EXCEPTIONT("Invalid datatype");
			}
		}

		// Apply find and replace
		if (strFind != "") {
			Announce("Applying find and replace");

			if (varData->type() == ncByte) {
				FindReplace<char>(dByteData, bFindByte, bReplaceByte);
			} else if (varData->type() == ncShort) {
				FindReplace<short>(dShortData, sFindShort, sReplaceShort);
			} else if (varData->type() == ncInt) {
				FindReplace<int>(dIntData, iFindInt, iReplaceInt);
			} else if (varData->type() == ncShort) {
				FindReplace<float>(dFloatData, dFindFloat, dReplaceFloat);
			} else if (varData->type() == ncDouble) {
				FindReplace<double>(dDoubleData, dFindDouble, dReplaceDouble);
			}
		}

		// Apply maps
		{
			if (fValidate) {
				if (varData->type() == ncByte) {
					dTotalMassIn = TotalMatrix<char>(dByteData, dXedge, dYedge);
				} else if (varData->type() == ncShort) {
					dTotalMassIn = TotalMatrix<short>(dShortData, dXedge, dYedge);
				} else if (varData->type() == ncInt) {
					dTotalMassIn = TotalMatrix<int>(dIntData, dXedge, dYedge);
				} else if (varData->type() == ncFloat) {
					dTotalMassIn = TotalMatrix<float>(dFloatData, dXedge, dYedge);
				} else if (varData->type() == ncDouble) {
					dTotalMassIn = TotalMatrix<double>(dDoubleData, dXedge, dYedge);
				} else {
					_EXCEPTIONT("Invalid datatype");
				}
				Announce("Input mass: %1.15e", dTotalMassIn);
			}

			Announce("Applying map");

			if (varData->type() == ncByte) {
				if (fOutputDouble) {
					ApplyMap<char,double>(dByteData, matX1D, matY1D, dDoubleDataOut);
				} else {
					ApplyMap<char>(dByteData, matX1D, matY1D, dByteDataOut);
				}

			} else if (varData->type() == ncShort) {
				if (fOutputDouble) {
					ApplyMap<short,double>(dShortData, matX1D, matY1D, dDoubleDataOut);
				} else {
					ApplyMap<short>(dShortData, matX1D, matY1D, dShortDataOut);
				}

			} else if (varData->type() == ncInt) {
				if (fOutputDouble) {
					ApplyMap<int,double>(dIntData, matX1D, matY1D, dDoubleDataOut);
				} else {
					ApplyMap<int>(dIntData, matX1D, matY1D, dIntDataOut);
				}

			} else if (varData->type() == ncShort) {
				if (fOutputDouble) {
					ApplyMap<float,double>(dFloatData, matX1D, matY1D, dDoubleDataOut);
				} else {
					ApplyMap<float>(dFloatData, matX1D, matY1D, dFloatDataOut);
				}

			} else if (varData->type() == ncDouble) {
				ApplyMap<double>(dDoubleData, matX1D, matY1D, dDoubleDataOut);
			}

			if (fValidate) {
				if ((fOutputDouble) || (varData->type() == ncDouble)) {
					dTotalMassOut = TotalMatrix<double>(dDoubleDataOut, dXedgeout, dYedgeout);
				} else if (varData->type() == ncByte) {
					dTotalMassOut = TotalMatrix<char>(dByteDataOut, dXedgeout, dYedgeout);
				} else if (varData->type() == ncShort) {
					dTotalMassOut = TotalMatrix<short>(dShortDataOut, dXedgeout, dYedgeout);
				} else if (varData->type() == ncInt) {
					dTotalMassOut = TotalMatrix<int>(dIntDataOut, dXedgeout, dYedgeout);
				} else if (varData->type() == ncFloat) {
					dTotalMassOut = TotalMatrix<float>(dFloatDataOut, dXedgeout, dYedgeout);
				} else {
					_EXCEPTIONT("Invalid datatype");
				}
				Announce("Output mass: %1.15e", dTotalMassOut);
			}

		}

		// Write
		{
			Announce("Writing data");
			if ((fOutputDouble) || (varData->type() == ncDouble)) {
				varDataOut->put(&(dDoubleDataOut[0][0]), lDataSizeOut);
			} else if (varData->type() == ncByte) {
				varDataOut->put(&(dByteDataOut[0][0]), lDataSizeOut);
			} else if (varData->type() == ncShort) {
				varDataOut->put(&(dShortDataOut[0][0]), lDataSizeOut);
			} else if (varData->type() == ncInt) {
				varDataOut->put(&(dIntDataOut[0][0]), lDataSizeOut);
			} else if (varData->type() == ncFloat) {
				varDataOut->put(&(dFloatDataOut[0][0]), lDataSizeOut);
			} else {
				_EXCEPTIONT("Invalid datatype");
			}
		}
	}

	AnnounceEndBlock("Done");

	return (0);

} catch(Exception & e) {
	Announce(e.ToString().c_str());
	return (-1);

} catch(...) {
	return (-2);
}
}

///////////////////////////////////////////////////////////////////////////////

