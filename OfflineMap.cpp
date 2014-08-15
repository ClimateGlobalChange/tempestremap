///////////////////////////////////////////////////////////////////////////////
///
///	\file    OfflineMap.cpp
///	\author  Paul Ullrich
///	\version August 14, 2014
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

#include "OfflineMap.h"

#include "netcdfcpp.h"
#include "NetCDFUtilities.h"

#include "Announce.h"
#include "Exception.h"
#include "DataVector.h"
#include "DataMatrix.h"

#include <cmath>

///////////////////////////////////////////////////////////////////////////////

void OfflineMap::InitializeOutputDimensionsFromFile(
	const std::string & strOutputMesh
) {
	NcFile ncOutputMesh(strOutputMesh.c_str(), NcFile::ReadOnly);

	// Check for rectilinear attribute
	bool fRectilinear = false;
	for (int a = 0; a < ncOutputMesh.num_atts(); a++) {
		if (strcmp(ncOutputMesh.get_att(a)->name(), "rectilinear")) {
			fRectilinear = true;
			break;
		}
	}

	// Non-rectilinear
	if (!fRectilinear) {
		int nCol = ncOutputMesh.get_dim("num_elem")->size();
		m_vecOutputDimSizes.push_back(nCol);
		m_vecOutputDimNames.push_back("ncol");
		return;
	}

	// Obtain rectilinear attributes
	int nDim0Size = ncOutputMesh.get_att("rectilinear_dim0_size")->as_int(0);
	int nDim1Size = ncOutputMesh.get_att("rectilinear_dim1_size")->as_int(0);

	std::string strDim0Name =
		ncOutputMesh.get_att("rectilinear_dim0_name")->as_string(0);
	std::string strDim1Name =
		ncOutputMesh.get_att("rectilinear_dim1_name")->as_string(0);

	m_vecOutputDimSizes.resize(2);
	m_vecOutputDimSizes[0] = nDim0Size;
	m_vecOutputDimSizes[1] = nDim1Size;

	m_vecOutputDimNames.resize(2);
	m_vecOutputDimNames[0] = strDim0Name;
	m_vecOutputDimNames[1] = strDim1Name;
}

///////////////////////////////////////////////////////////////////////////////

void OfflineMap::Apply(
	const DataVector<double> & vecAreaInput,
	const DataVector<double> & vecAreaOutput,
	const std::string & strInputDataFile,
	const std::string & strOutputDataFile,
	const std::vector<std::string> & vecVariables
) {
	NcFile ncInput(strInputDataFile.c_str(), NcFile::ReadOnly);
	NcFile ncOutput(strOutputDataFile.c_str(), NcFile::Replace);

	// Check for time dimension
	NcDim * dimTime = NULL;
	for (int d = 0; d < ncInput.num_dims(); d++) {
		if (strcmp(ncInput.get_dim(d)->name(), "time") == 0) {
			dimTime = ncInput.get_dim(d);
			break;
		}
	}

	bool fHasTime = false;
	int nTime = 1;
	if (dimTime != NULL) {
		nTime = dimTime->size();
		fHasTime = true;
	}

	// Check for ncol dimension
	NcDim * dimNCol = ncInput.get_dim("ncol");
	int nCol = dimNCol->size();

	if (nCol != m_mapRemap.GetColumns()) {
		_EXCEPTION2("\nMismatch between offline map size (%i) and "
			"data.ncol (%i)", m_mapRemap.GetColumns(), nCol);
	}

	// Output columns
	int nColOut = m_mapRemap.GetRows();

	// Data
	DataVector<float> dataIn;
	dataIn.Initialize(nCol);

	DataMatrix<float> dataOut;
	if (m_vecOutputDimSizes.size() == 1) {
		dataOut.Initialize(m_vecOutputDimSizes[0], 1);
	} else {
		dataOut.Initialize(m_vecOutputDimSizes[0], m_vecOutputDimSizes[1]);
	}

	DataVector<double> dataInDouble;
	dataInDouble.Initialize(nCol);

	DataVector<double> dataOutDouble;
	dataOutDouble.Initialize(nColOut);

	// Output
	CopyNcFileAttributes(&ncInput, &ncOutput);

	NcDim * dimTimeOut = NULL;
	if (fHasTime) {
		dimTimeOut = ncOutput.add_dim("time", dimTime->size());
	}

	std::vector<NcDim *> vecDim;
	for (int d = 0; d < m_vecOutputDimSizes.size(); d++) {
		vecDim.push_back(ncOutput.add_dim(
			m_vecOutputDimNames[d].c_str(),
			m_vecOutputDimSizes[d]));
	}

	// Loop through all variables
	for (int v = 0; v < vecVariables.size(); v++) {
		NcVar * var = ncInput.get_var(vecVariables[v].c_str());

		AnnounceStartBlock(vecVariables[v].c_str());

		// Create new output variable
		NcVar * varOut = NULL;
		if (fHasTime) {
			if (vecDim.size() == 1) {
				varOut = ncOutput.add_var(
					vecVariables[v].c_str(), ncFloat,
						dimTimeOut, vecDim[0]);
			} else {
				varOut = ncOutput.add_var(
					vecVariables[v].c_str(), ncFloat,
						dimTimeOut, vecDim[0], vecDim[1]);
			}

		} else {
			if (vecDim.size() == 1) {
				varOut = ncOutput.add_var(
					vecVariables[v].c_str(), ncFloat, vecDim[0]);
			} else {
				varOut = ncOutput.add_var(
					vecVariables[v].c_str(), ncFloat,
						vecDim[0], vecDim[1]);
			}
		}

		CopyNcVarAttributes(var, varOut);

		// Loop through all times
		for (int t = 0; t < nTime; t++) {

			// Get the data
			if (fHasTime) {
				var->set_cur(t,0);
				var->get(&(dataIn[0]), 1, nCol);
			} else {
				var->set_cur((long)0);
				var->get(&(dataIn[0]), nCol);
			}

			// Cast the data to double
			for (int i = 0; i < nCol; i++) {
				dataInDouble[i] = static_cast<double>(dataIn[i]);
			}

			// Announce input mass
			double dInputMass = 0.0;
			for (int i = 0; i < nCol; i++) {
				dInputMass += dataInDouble[i] * vecAreaInput[i];
			}
			Announce("Input Mass:  %1.15e", dInputMass);

			// Apply the offline map to the data
			m_mapRemap.Apply(dataInDouble, dataOutDouble);

			// Announce output mass
			double dOutputMass = 0.0;
			for (int i = 0; i < nColOut; i++) {
				dOutputMass += dataOutDouble[i] * vecAreaOutput[i];
			}
			Announce("Output Mass: %1.15e", dOutputMass);

			// Cast the data to float
			int ix = 0;
			for (int i = 0; i < dataOut.GetRows(); i++) {
			for (int j = 0; j < dataOut.GetColumns(); j++) {
				dataOut[i][j] = static_cast<float>(dataOutDouble[ix]);
				ix++;
			}
			}

			// Write the data
			if (fHasTime) {
				if (m_vecOutputDimSizes.size() == 1) {
					varOut->set_cur(t,0);
					varOut->put(&(dataOut[0][0]),
						1, dataOut.GetRows());
				} else {
					varOut->set_cur(t,0,0);
					varOut->put(&(dataOut[0][0]),
						1, dataOut.GetRows(), dataOut.GetColumns());
				}
			} else {
				if (m_vecOutputDimSizes.size() == 1) {
					varOut->set_cur((long)0);
					varOut->put(&(dataOut[0][0]), dataOut.GetRows());
				} else {
					varOut->set_cur(0,0);
					varOut->put(&(dataOut[0][0]),
						dataOut.GetRows(), dataOut.GetColumns());
				}
			}
		}
		AnnounceEndBlock(NULL);
	}
}

///////////////////////////////////////////////////////////////////////////////

void OfflineMap::Read(
	const std::string & strInput
) {
	NcFile ncMap(strInput.c_str(), NcFile::ReadOnly);

	NcDim * dimNA = ncMap.get_dim("n_s");

	NcVar * varRow = ncMap.get_var("row");
	NcVar * varCol = ncMap.get_var("col");
	NcVar * varS   = ncMap.get_var("S");

	int nS = dimNA->size();

	DataVector<int> vecRow;
	vecRow.Initialize(nS);

	DataVector<int> vecCol;
	vecCol.Initialize(nS);

	DataVector<double> vecS;
	vecS.Initialize(nS);

	varRow->set_cur((long)0);
	varRow->get(&(vecRow[0]), nS);

	varCol->set_cur((long)0);
	varCol->get(&(vecCol[0]), nS);

	varS->set_cur((long)0);
	varS->get(&(vecS[0]), nS);

	m_mapRemap.SetEntries(vecRow, vecCol, vecS);
}

///////////////////////////////////////////////////////////////////////////////

void OfflineMap::Write(
	const std::string & strOutput
) {
	NcFile ncMap(strOutput.c_str(), NcFile::Replace);

	DataVector<int> vecRow;
	DataVector<int> vecCol;
	DataVector<double> vecS;

	m_mapRemap.GetEntries(vecRow, vecCol, vecS);

	int nS = vecRow.GetRows();

	NcDim * dimNA = ncMap.add_dim("n_s", nS);

	NcVar * varRow = ncMap.add_var("row", ncInt, dimNA);
	NcVar * varCol = ncMap.add_var("col", ncInt, dimNA);
	NcVar * varS = ncMap.add_var("S", ncDouble, dimNA);

	varRow->set_cur((long)0);
	varRow->put(&(vecRow[0]), nS);

	varCol->set_cur((long)0);
	varCol->put(&(vecCol[0]), nS);

	varS->set_cur((long)0);
	varS->put(&(vecS[0]), nS);
}

///////////////////////////////////////////////////////////////////////////////

bool OfflineMap::IsFirstOrder(
	double dTolerance
) {

	// Get map entries
	DataVector<int> dataRows;
	DataVector<int> dataCols;
	DataVector<double> dataEntries;

	m_mapRemap.GetEntries(dataRows, dataCols, dataEntries);

	// Calculate row sums
	DataVector<double> dRowSums;
	dRowSums.Initialize(m_mapRemap.GetRows());

	for (int i = 0; i < dataRows.GetRows(); i++) {
		dRowSums[dataRows[i]] += dataEntries[i];
	}

	// Verify all row sums are equal to 1
	bool fFirstOrder = true;
	for (int i = 0; i < dRowSums.GetRows(); i++) {
		if (fabs(dRowSums[i] - 1.0) > dTolerance) {
			fFirstOrder = false;
			Announce("OfflineMap is not first order in row %i (%1.15e)",
				i, dRowSums[i]);
		}
	}

	return fFirstOrder;
}

///////////////////////////////////////////////////////////////////////////////

bool OfflineMap::IsConservative(
	const DataVector<double> & vecInputAreas,
	const DataVector<double> & vecOutputAreas,
	double dTolerance
) {
	if (vecInputAreas.GetRows() != m_mapRemap.GetColumns()) {
		_EXCEPTIONT("vecInputAreas has not been computed");
	}
	if (vecOutputAreas.GetRows() != m_mapRemap.GetRows()) {
		_EXCEPTIONT("vecOutputAreas has not been computed");
	}

	// Get map entries
	DataVector<int> dataRows;
	DataVector<int> dataCols;
	DataVector<double> dataEntries;

	m_mapRemap.GetEntries(dataRows, dataCols, dataEntries);

	// Calculate column sums
	DataVector<double> dColumnSums;
	dColumnSums.Initialize(m_mapRemap.GetColumns());

	for (int i = 0; i < dataRows.GetRows(); i++) {
		dColumnSums[dataCols[i]] +=
			dataEntries[i] * vecOutputAreas[dataRows[i]];
	}

	// Verify all column sums equal the input Jacobian
	bool fConservative = true;
	for (int i = 0; i < dColumnSums.GetRows(); i++) {
		if (fabs(dColumnSums[i] - vecInputAreas[i]) > dTolerance) {
			fConservative = false;
			Announce("OfflineMap is not conservative in column "
				"%i (%1.15e / %1.15e)",
				i, dColumnSums[i], vecInputAreas[i]);
		}
	}

	return fConservative;
}

///////////////////////////////////////////////////////////////////////////////

