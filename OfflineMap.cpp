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

void OfflineMap::InitializeInputDimensionsFromFile(
	const std::string & strInputMesh
) {
	NcFile ncInputMesh(strInputMesh.c_str(), NcFile::ReadOnly);

	// Check for rectilinear attribute
	bool fRectilinear = false;
	for (int a = 0; a < ncInputMesh.num_atts(); a++) {
		if (strcmp(ncInputMesh.get_att(a)->name(), "rectilinear") == 0) {
			fRectilinear = true;
			break;
		}
	}

	// Non-rectilinear
	if (!fRectilinear) {
		int nElements = ncInputMesh.get_dim("num_elem")->size();
		m_vecInputDimSizes.push_back(nElements);
		m_vecInputDimNames.push_back("num_elem");
		return;
	}

	// Obtain rectilinear attributes
	int nDim0Size = ncInputMesh.get_att("rectilinear_dim0_size")->as_int(0);
	int nDim1Size = ncInputMesh.get_att("rectilinear_dim1_size")->as_int(0);

	std::string strDim0Name =
		ncInputMesh.get_att("rectilinear_dim0_name")->as_string(0);
	std::string strDim1Name =
		ncInputMesh.get_att("rectilinear_dim1_name")->as_string(0);

	m_vecInputDimSizes.resize(2);
	m_vecInputDimSizes[0] = nDim0Size;
	m_vecInputDimSizes[1] = nDim1Size;

	m_vecInputDimNames.resize(2);
	m_vecInputDimNames[0] = strDim0Name;
	m_vecInputDimNames[1] = strDim1Name;
}

///////////////////////////////////////////////////////////////////////////////

void OfflineMap::InitializeOutputDimensionsFromFile(
	const std::string & strOutputMesh
) {
	NcFile ncOutputMesh(strOutputMesh.c_str(), NcFile::ReadOnly);

	// Check for rectilinear attribute
	bool fRectilinear = false;
	for (int a = 0; a < ncOutputMesh.num_atts(); a++) {
		if (strcmp(ncOutputMesh.get_att(a)->name(), "rectilinear") == 0) {
			fRectilinear = true;
			break;
		}
	}

	// Non-rectilinear
	if (!fRectilinear) {
		int nElements = ncOutputMesh.get_dim("num_elem")->size();
		m_vecOutputDimSizes.push_back(nElements);
		m_vecOutputDimNames.push_back("num_elem");
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

NcDim * NcFile_GetDimIfExists(
	NcFile & ncInput,
	const std::string & strDimName,
	int nSize
) {
	for (int d = 0; d < ncInput.num_dims(); d++) {
		NcDim * dim = ncInput.get_dim(d);
		if (strcmp(dim->name(), strDimName.c_str()) == 0) {
			if (dim->size() != nSize) {
				_EXCEPTION3("NetCDF file has dimension \"%s\" with mismatched"
					" size %i != %i", strDimName.c_str(), dim->size(), nSize);
			}
			return ncInput.get_dim(d);
		}
	}
	return ncInput.add_dim(strDimName.c_str(), nSize);
}

///////////////////////////////////////////////////////////////////////////////

void OfflineMap::Apply(
	const std::string & strInputDataFile,
	const std::string & strOutputDataFile,
	const std::vector<std::string> & vecVariables,
	const std::string & strNColName,
	bool fOutputDouble,
	bool fAppend
) {
	// Open input file
	NcFile ncInput(strInputDataFile.c_str(), NcFile::ReadOnly);

	// Open output file
	NcFile::FileMode eOpenMode = NcFile::Replace;
	if (fAppend) {
		eOpenMode = NcFile::Write;
	}

	NcFile ncOutput(strOutputDataFile.c_str(), eOpenMode);

	// Number of source and target regions
	int nSourceCount = m_dInputAreas.GetRows();
	int nTargetCount = m_dOutputAreas.GetRows();

	// Check for rectilinear data
	bool fInputRectilinear;
	if (m_vecInputDimSizes.size() == 1) {
		fInputRectilinear = false;
	} else if (m_vecInputDimSizes.size() == 2) {
		fInputRectilinear = true;
	} else {
		_EXCEPTIONT("m_vecInputDimSizes undefined");
	}

	if (fInputRectilinear) {
		if (nSourceCount != m_vecInputDimSizes[0] * m_vecInputDimSizes[1]) {
			_EXCEPTIONT("Rectilinear input expected in finite volume form");
		}
	}

	bool fOutputRectilinear;
	if (m_vecOutputDimSizes.size() == 1) {
		fOutputRectilinear = false;
	} else if (m_vecOutputDimSizes.size() == 2) {
		fOutputRectilinear = true;
	} else {
		_EXCEPTIONT("m_vecOutputDimSizes undefined");
	}

	if (fOutputRectilinear) {
		if (nTargetCount != m_vecOutputDimSizes[0] * m_vecOutputDimSizes[1]) {
			printf("%i %i\n", nTargetCount, m_vecOutputDimSizes[0] * m_vecOutputDimSizes[1]);
			_EXCEPTIONT("Rectilinear output expected in finite volume form");
		}
	}

	DataVector<float> dataIn;
	if (m_vecInputDimSizes.size() == 1) {
		dataIn.Initialize(nSourceCount);
	} else {
		dataIn.Initialize(m_vecInputDimSizes[0] * m_vecInputDimSizes[1]);
	}

	DataVector<float> dataOut;
	if (m_vecOutputDimSizes.size() == 1) {
		dataOut.Initialize(nTargetCount);
	} else {
		dataOut.Initialize(m_vecOutputDimSizes[0] * m_vecOutputDimSizes[1]);
	}

	DataVector<double> dataInDouble;
	dataInDouble.Initialize(nSourceCount);

	DataVector<double> dataOutDouble;
	dataOutDouble.Initialize(nTargetCount);

	// Output
	if (!fAppend) {
		CopyNcFileAttributes(&ncInput, &ncOutput);
	}

	NcDim * dim0;
	NcDim * dim1;

	dim0 = NcFile_GetDimIfExists(
		ncOutput,
		m_vecOutputDimNames[0].c_str(),
		m_vecOutputDimSizes[0]);

	if (fOutputRectilinear) {
		dim1 = NcFile_GetDimIfExists(
			ncOutput,
			m_vecOutputDimNames[1].c_str(),
			m_vecOutputDimSizes[1]);
	}

	// Loop through all variables
	for (int v = 0; v < vecVariables.size(); v++) {
		NcVar * var = ncInput.get_var(vecVariables[v].c_str());

		AnnounceStartBlock(vecVariables[v].c_str());

		// Loop through all dimensions and add missing dimensions to output
		int nVarTotalEntries = 1;

		DataVector<NcDim *> vecDims;
		vecDims.Initialize(var->num_dims());

		DataVector<NcDim *> vecDimsOut;
		if (fOutputRectilinear) {
			if (fInputRectilinear) {
				vecDimsOut.Initialize(var->num_dims());
			} else {
				vecDimsOut.Initialize(var->num_dims()+1);
			}

			vecDimsOut[vecDimsOut.GetRows()-2] = dim0;
			vecDimsOut[vecDimsOut.GetRows()-1] = dim1;
		} else {
			if (fInputRectilinear) {
				vecDimsOut.Initialize(var->num_dims()-1);
			} else {
				vecDimsOut.Initialize(var->num_dims());
			}

			vecDimsOut[vecDimsOut.GetRows()-1] = dim0;
		}

		int nFreeDims = var->num_dims() - m_vecInputDimSizes.size();

		DataVector<long> vecDimSizes;
		vecDimSizes.Initialize(nFreeDims);

		for (int d = 0; d < nFreeDims; d++) {
			vecDims[d] = var->get_dim(d);

			long nDimSize = vecDims[d]->size();
			std::string strDimName = vecDims[d]->name();

			vecDimSizes[d] = nDimSize;
			nVarTotalEntries *= nDimSize;

			vecDimsOut[d] =
				NcFile_GetDimIfExists(
					ncOutput,
					strDimName.c_str(),
					nDimSize);
		}

		// Create new output variable
		NcVar * varOut;
		if (fOutputDouble) {
			varOut =
				ncOutput.add_var(
					vecVariables[v].c_str(),
					ncDouble,
					vecDimsOut.GetRows(),
					(const NcDim**)&(vecDimsOut[0]));

		} else {
			varOut =
				ncOutput.add_var(
					vecVariables[v].c_str(),
					ncFloat,
					vecDimsOut.GetRows(),
					(const NcDim**)&(vecDimsOut[0]));
		}

		CopyNcVarAttributes(var, varOut);

		// Input and output counts
		DataVector<long> nCountsIn;
		nCountsIn.Initialize(vecDims.GetRows());

		DataVector<long> nCountsOut;
		nCountsOut.Initialize(vecDimsOut.GetRows());

		// Get size
		DataVector<long> nGet;
		nGet.Initialize(nCountsIn.GetRows());
		for (int d = 0; d < nGet.GetRows()-1; d++) {
			nGet[d] = 1;
		}
		if (fInputRectilinear) {
			nGet[nGet.GetRows()-2] = m_vecInputDimSizes[0];
			nGet[nGet.GetRows()-1] = m_vecInputDimSizes[1];
		} else {
			nGet[nGet.GetRows()-1] = nSourceCount;
		}

		// Put size
		DataVector<long> nPut;
		nPut.Initialize(nCountsOut.GetRows());
		for (int d = 0; d < nPut.GetRows()-1; d++) {
			nPut[d] = 1;
		}
		if (fOutputRectilinear) {
			nPut[nPut.GetRows()-2] = m_vecOutputDimSizes[0];
			nPut[nPut.GetRows()-1] = m_vecOutputDimSizes[1];
		} else {
			nPut[nPut.GetRows()-1] = nTargetCount;
		}

		// Loop through all entries
		for (int t = 0; t < nVarTotalEntries; t++) {

			long tt = static_cast<long>(t);
			for (int d = 0; d < vecDimSizes.GetRows(); d++) {
				nCountsIn[d]  = tt % vecDimSizes[d];
				nCountsOut[d] = tt % vecDimSizes[d];
				tt -= nCountsOut[d] * vecDimSizes[d];
			}

			for (int d = vecDimSizes.GetRows(); d < nCountsIn.GetRows(); d++) {
				nCountsIn[d] = 0;
			}
			for (int d = vecDimSizes.GetRows(); d < nCountsOut.GetRows(); d++) {
				nCountsOut[d] = 0;
			}

			// Get the data
			var->set_cur(&(nCountsIn[0]));

			// Load data as Float, cast to Double
			if (var->type() == ncFloat) {
				var->get(&(dataIn[0]), &(nGet[0]));

				for (int i = 0; i < nSourceCount; i++) {
					dataInDouble[i] = static_cast<double>(dataIn[i]);
				}

			// Load data as Double
			} else if (var->type() == ncDouble) {
				var->get(&(dataInDouble[0]), &(nGet[0]));

			} else {
				_EXCEPTIONT("Invalid variable type");
			}

			// Announce input mass
			double dInputMass = 0.0;
			double dInputMin  = dataInDouble[0];
			double dInputMax  = dataInDouble[0];
			for (int i = 0; i < nSourceCount; i++) {
				dInputMass += dataInDouble[i] * m_dInputAreas[i];
				if (dataInDouble[i] < dInputMin) {
					dInputMin = dataInDouble[i];
				}
				if (dataInDouble[i] > dInputMax) {
					dInputMax = dataInDouble[i];
				}
			}
			Announce(" Input Mass: %1.15e Min %1.10e Max %1.10e",
				dInputMass, dInputMin, dInputMax);

			// Apply the offline map to the data
			m_mapRemap.Apply(dataInDouble, dataOutDouble);

			// Announce output mass
			double dOutputMass = 0.0;
			double dOutputMin  = dataOutDouble[0];
			double dOutputMax  = dataOutDouble[0];
			for (int i = 0; i < nTargetCount; i++) {
				dOutputMass += dataOutDouble[i] * m_dOutputAreas[i];
				if (dataOutDouble[i] < dOutputMin) {
					dOutputMin = dataOutDouble[i];
				}
				if (dataOutDouble[i] > dOutputMax) {
					dOutputMax = dataOutDouble[i];
				}
			}
			Announce("Output Mass: %1.15e Min %1.10e Max %1.10e",
				dOutputMass, dOutputMin, dOutputMax);

			// Write the data
			if (fOutputDouble) {
				varOut->set_cur(&(nCountsOut[0]));
				varOut->put(&(dataOutDouble[0]), &(nPut[0]));

			} else {
				// Cast the data to float
				for (int i = 0; i < dataOut.GetRows(); i++) {
					dataOut[i] = static_cast<float>(dataOutDouble[i]);
				}

				// Write the data as float
				varOut->set_cur(&(nCountsOut[0]));
				varOut->put(&(dataOut[0]), &(nPut[0]));
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

	// Read input dimensions entries
	NcDim * dimSrcGridRank = ncMap.get_dim("src_grid_rank");
	NcDim * dimDstGridRank = ncMap.get_dim("dst_grid_rank");

	int nSrcGridDims = (int)(dimSrcGridRank->size());
	int nDstGridDims = (int)(dimDstGridRank->size());

	NcVar * varSrcGridDims = ncMap.get_var("src_grid_dims");
	NcVar * varDstGridDims = ncMap.get_var("dst_grid_dims");

	m_vecInputDimSizes.resize(nSrcGridDims);
	m_vecInputDimNames.resize(nSrcGridDims);

	m_vecOutputDimSizes.resize(nDstGridDims);
	m_vecOutputDimNames.resize(nDstGridDims);

	varSrcGridDims->get(&(m_vecInputDimSizes[0]), nSrcGridDims);
	varDstGridDims->get(&(m_vecOutputDimSizes[0]), nDstGridDims);

	for (int i = 0; i < nSrcGridDims; i++) {
		char szDim[64];
		sprintf(szDim, "name%i", i);
		m_vecInputDimNames[i] = varSrcGridDims->get_att(szDim)->as_string(0);
	}

	for (int i = 0; i < nDstGridDims; i++) {
		char szDim[64];
		sprintf(szDim, "name%i", i);
		m_vecOutputDimNames[i] = varDstGridDims->get_att(szDim)->as_string(0);
	}

	// Input and Output mesh resolutions
	NcDim * dimNA = ncMap.get_dim("n_a");
	NcDim * dimNB = ncMap.get_dim("n_b");

	int nA = dimNA->size();
	int nB = dimNB->size();

	m_dInputAreas.Initialize(nA);
	m_dOutputAreas.Initialize(nB);

	// Read areas
	NcVar * varAreaA = ncMap.get_var("area_a");
	varAreaA->get(&(m_dInputAreas[0]), nA);

	NcVar * varAreaB = ncMap.get_var("area_b");
	varAreaB->get(&(m_dOutputAreas[0]), nB);

	// Read SparseMatrix entries
	NcDim * dimNS = ncMap.get_dim("n_s");

	NcVar * varRow = ncMap.get_var("row");
	NcVar * varCol = ncMap.get_var("col");
	NcVar * varS   = ncMap.get_var("S");

	int nS = dimNS->size();

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

	// Attributes
	ncMap.add_att("Title", "TempestRemap Offline Regridding Weight Generator");

	// Write output dimensions entries
	int nSrcGridDims = (int)(m_vecInputDimSizes.size());
	int nDstGridDims = (int)(m_vecOutputDimSizes.size());

	NcDim * dimSrcGridRank = ncMap.add_dim("src_grid_rank", nSrcGridDims);
	NcDim * dimDstGridRank = ncMap.add_dim("dst_grid_rank", nDstGridDims);

	NcVar * varSrcGridDims =
		ncMap.add_var("src_grid_dims", ncInt, dimSrcGridRank);
	NcVar * varDstGridDims =
		ncMap.add_var("dst_grid_dims", ncInt, dimDstGridRank);

	varSrcGridDims->put(&(m_vecInputDimSizes[0]), nSrcGridDims);
	varDstGridDims->put(&(m_vecOutputDimSizes[0]), nDstGridDims);

	int nA = (int)(m_dInputAreas.GetRows());
	int nB = (int)(m_dOutputAreas.GetRows());

	for (int i = 0; i < m_vecInputDimSizes.size(); i++) {
		char szDim[64];
		sprintf(szDim, "name%i", i);
		varSrcGridDims->add_att(szDim, m_vecInputDimNames[i].c_str());
	}

	for (int i = 0; i < m_vecOutputDimSizes.size(); i++) {
		char szDim[64];
		sprintf(szDim, "name%i", i);
		varDstGridDims->add_att(szDim, m_vecOutputDimNames[i].c_str());
	}

	// Input and Output mesh resolutions
	NcDim * dimNA = ncMap.add_dim("n_a", nA);
	NcDim * dimNB = ncMap.add_dim("n_b", nB);

	// Write areas
	NcVar * varAreaA = ncMap.add_var("area_a", ncDouble, dimNA);
	varAreaA->put(&(m_dInputAreas[0]), nA);

	NcVar * varAreaB = ncMap.add_var("area_b", ncDouble, dimNB);
	varAreaB->put(&(m_dOutputAreas[0]), nB);

	// Write frac
	DataVector<double> dFrac;

	dFrac.Initialize(nA);
	for (int i = 0; i < nA; i++) {
		dFrac[i] = 1.0;
	}
	NcVar * varFracA = ncMap.add_var("frac_a", ncDouble, dimNA);
	varFracA->put(&(dFrac[0]), nA);

	dFrac.Initialize(nB);
	for (int i = 0; i < nB; i++) {
		dFrac[i] = 1.0;
	}
	NcVar * varFracB = ncMap.add_var("frac_b", ncDouble, dimNB);
	varFracB->put(&(dFrac[0]), nB);

	// Write SparseMatrix entries
	DataVector<int> vecRow;
	DataVector<int> vecCol;
	DataVector<double> vecS;

	m_mapRemap.GetEntries(vecRow, vecCol, vecS);

	int nS = vecRow.GetRows();
	NcDim * dimNS = ncMap.add_dim("n_s", nS);

	NcVar * varRow = ncMap.add_var("row", ncInt, dimNS);
	NcVar * varCol = ncMap.add_var("col", ncInt, dimNS);
	NcVar * varS = ncMap.add_var("S", ncDouble, dimNS);

	varRow->set_cur((long)0);
	varRow->put(&(vecRow[0]), nS);

	varCol->set_cur((long)0);
	varCol->put(&(vecCol[0]), nS);

	varS->set_cur((long)0);
	varS->put(&(vecS[0]), nS);
}

///////////////////////////////////////////////////////////////////////////////

bool OfflineMap::IsConsistent(
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
	bool fConsistent = true;
	for (int i = 0; i < dRowSums.GetRows(); i++) {
		if (fabs(dRowSums[i] - 1.0) > dTolerance) {
			fConsistent = false;
			Announce("OfflineMap is not consistent in row %i (%1.15e)",
				i, dRowSums[i]);
		}
	}

	return fConsistent;
}

///////////////////////////////////////////////////////////////////////////////

bool OfflineMap::IsConservative(
	double dTolerance
) {
/*
	if (vecInputAreas.GetRows() != m_mapRemap.GetColumns()) {
		_EXCEPTIONT("vecInputAreas / mapRemap dimension mismatch");
	}
	if (vecOutputAreas.GetRows() != m_mapRemap.GetRows()) {
		_EXCEPTIONT("vecOutputAreas / mapRemap dimension mismatch");
	}
*/
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
			dataEntries[i] * m_dOutputAreas[dataRows[i]];
	}

	// Verify all column sums equal the input Jacobian
	bool fConservative = true;
	for (int i = 0; i < dColumnSums.GetRows(); i++) {
		if (fabs(dColumnSums[i] - m_dInputAreas[i]) > dTolerance) {
			fConservative = false;
			Announce("OfflineMap is not conservative in column "
				"%i (%1.15e / %1.15e)",
				i, dColumnSums[i], m_dInputAreas[i]);
		}
	}

	return fConservative;
}

///////////////////////////////////////////////////////////////////////////////

bool OfflineMap::IsMonotone(
	double dTolerance
) {

	// Get map entries
	DataVector<int> dataRows;
	DataVector<int> dataCols;
	DataVector<double> dataEntries;

	m_mapRemap.GetEntries(dataRows, dataCols, dataEntries);

	// Verify all entries are in the range [0,1]
	bool fMonotone = true;
	for (int i = 0; i < dataRows.GetRows(); i++) {
		if ((dataEntries[i] < -dTolerance) ||
			(dataEntries[i] > 1.0 + dTolerance)
		) {
			fMonotone = false;

			Announce("OfflineMap is not monotone in entry (%i): %1.15e",
				i, dataEntries[i]);
		}
	}

	return fMonotone;
}

///////////////////////////////////////////////////////////////////////////////

