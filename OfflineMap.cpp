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
#include "GridElements.h"
#include "FiniteElementTools.h"

#include "Announce.h"
#include "Exception.h"
#include "DataVector.h"
#include "DataMatrix.h"

#include <cmath>

///////////////////////////////////////////////////////////////////////////////

void OfflineMap::InitializeSourceDimensionsFromFile(
	const std::string & strSourceMesh
) {
	// Open the source mesh
	NcFile ncSourceMesh(strSourceMesh.c_str(), NcFile::ReadOnly);
	if (!ncSourceMesh.is_valid()) {
		_EXCEPTION1("Unable to open mesh \"%s\"", strSourceMesh.c_str());
	}

	// Check for grid dimensions (SCRIP format grid)
	NcVar * varGridDims = ncSourceMesh.get_var("grid_dims");
	if (varGridDims != NULL) {
		NcDim * dimGridRank = varGridDims->get_dim(0);

		m_vecSourceDimSizes.resize(dimGridRank->size());
		varGridDims->get(&(m_vecSourceDimSizes[0]), dimGridRank->size());

		if (dimGridRank->size() == 1) {
			m_vecSourceDimNames.push_back("num_elem");
		} else if (dimGridRank->size() == 2) {
			m_vecSourceDimNames.push_back("lon");
			m_vecSourceDimNames.push_back("lat");
		} else {
			_EXCEPTIONT("Source grid grid_rank must be < 3");
		}
		return;
	}

	// Check for rectilinear attribute
	NcAtt * attRectilinear = ncSourceMesh.get_att("rectilinear");

	// No rectilinear attribute
	if (attRectilinear == NULL) {
		int nElements = ncSourceMesh.get_dim("num_elem")->size();
		m_vecSourceDimSizes.push_back(nElements);
		m_vecSourceDimNames.push_back("num_elem");
		return;
	}

	// Obtain rectilinear attributes (dimension sizes)
	NcAtt * attRectilinearDim0Size =
		ncSourceMesh.get_att("rectilinear_dim0_size");
	NcAtt * attRectilinearDim1Size =
		ncSourceMesh.get_att("rectilinear_dim1_size");

	if (attRectilinearDim0Size == NULL) {
		_EXCEPTIONT("Missing attribute \"rectilinear_dim0_size\"");
	}
	if (attRectilinearDim1Size == NULL) {
		_EXCEPTIONT("Missing attribute \"rectilinear_dim1_size\"");
	}

	int nDim0Size = attRectilinearDim0Size->as_int(0);
	int nDim1Size = attRectilinearDim1Size->as_int(0);

	// Obtain rectilinear attributes (dimension names)
	NcAtt * attRectilinearDim0Name =
		ncSourceMesh.get_att("rectilinear_dim0_name");
	NcAtt * attRectilinearDim1Name =
		ncSourceMesh.get_att("rectilinear_dim1_name");

	if (attRectilinearDim0Name == NULL) {
		_EXCEPTIONT("Missing attribute \"rectilinear_dim0_name\"");
	}
	if (attRectilinearDim1Name == NULL) {
		_EXCEPTIONT("Missing attribute \"rectilinear_dim1_name\"");
	}

	std::string strDim0Name = attRectilinearDim0Name->as_string(0);
	std::string strDim1Name = attRectilinearDim1Name->as_string(0);

	// Push rectilinear attributes into array
	m_vecSourceDimSizes.resize(2);
	m_vecSourceDimSizes[0] = nDim0Size;
	m_vecSourceDimSizes[1] = nDim1Size;

	m_vecSourceDimNames.resize(2);
	m_vecSourceDimNames[0] = strDim0Name;
	m_vecSourceDimNames[1] = strDim1Name;
}

///////////////////////////////////////////////////////////////////////////////

void OfflineMap::InitializeTargetDimensionsFromFile(
	const std::string & strTargetMesh
) {
	// Open the source mesh
	NcFile ncTargetMesh(strTargetMesh.c_str(), NcFile::ReadOnly);
	if (!ncTargetMesh.is_valid()) {
		_EXCEPTION1("Unable to open mesh \"%s\"", strTargetMesh.c_str());
	}

	// Check for grid dimensions (SCRIP format grid)
	NcVar * varGridDims = ncTargetMesh.get_var("grid_dims");
	if (varGridDims != NULL) {
		NcDim * dimGridRank = varGridDims->get_dim(0);

		m_vecTargetDimSizes.resize(dimGridRank->size());
		varGridDims->get(&(m_vecTargetDimSizes[0]), dimGridRank->size());

		if (dimGridRank->size() == 1) {
			m_vecTargetDimNames.push_back("num_elem");
		} else if (dimGridRank->size() == 2) {
			m_vecTargetDimNames.push_back("lon");
			m_vecTargetDimNames.push_back("lat");
		} else {
			_EXCEPTIONT("Target grid grid_rank must be < 3");
		}
		return;
	}

	// Check for rectilinear attribute
	NcAtt * attRectilinear = ncTargetMesh.get_att("rectilinear");

	// No rectilinear attribute
	if (attRectilinear == NULL) {
		int nElements = ncTargetMesh.get_dim("num_elem")->size();
		m_vecTargetDimSizes.push_back(nElements);
		m_vecTargetDimNames.push_back("num_elem");
		return;
	}

	// Obtain rectilinear attributes (dimension sizes)
	NcAtt * attRectilinearDim0Size =
		ncTargetMesh.get_att("rectilinear_dim0_size");
	NcAtt * attRectilinearDim1Size =
		ncTargetMesh.get_att("rectilinear_dim1_size");

	if (attRectilinearDim0Size == NULL) {
		_EXCEPTIONT("Missing attribute \"rectilinear_dim0_size\"");
	}
	if (attRectilinearDim1Size == NULL) {
		_EXCEPTIONT("Missing attribute \"rectilinear_dim1_size\"");
	}

	int nDim0Size = attRectilinearDim0Size->as_int(0);
	int nDim1Size = attRectilinearDim1Size->as_int(0);

	// Obtain rectilinear attributes (dimension names)
	NcAtt * attRectilinearDim0Name =
		ncTargetMesh.get_att("rectilinear_dim0_name");
	NcAtt * attRectilinearDim1Name =
		ncTargetMesh.get_att("rectilinear_dim1_name");

	if (attRectilinearDim0Name == NULL) {
		_EXCEPTIONT("Missing attribute \"rectilinear_dim0_name\"");
	}
	if (attRectilinearDim1Name == NULL) {
		_EXCEPTIONT("Missing attribute \"rectilinear_dim1_name\"");
	}

	std::string strDim0Name = attRectilinearDim0Name->as_string(0);
	std::string strDim1Name = attRectilinearDim1Name->as_string(0);

	// Push rectilinear attributes into array
	m_vecTargetDimSizes.resize(2);
	m_vecTargetDimSizes[0] = nDim1Size;
	m_vecTargetDimSizes[1] = nDim0Size;

	m_vecTargetDimNames.resize(2);
	m_vecTargetDimNames[0] = strDim1Name;
	m_vecTargetDimNames[1] = strDim0Name;
}

///////////////////////////////////////////////////////////////////////////////

void OfflineMap::InitializeCoordinatesFromMeshFV(
	const Mesh & mesh,
	DataVector<double> & dCenterLon,
	DataVector<double> & dCenterLat,
	DataMatrix<double> & dVertexLon,
	DataMatrix<double> & dVertexLat
) {
	int nFaces = mesh.faces.size();

	dCenterLon.Initialize(nFaces);
	dCenterLat.Initialize(nFaces);

	// Count maximum number of Nodes per Face
	int nNodesPerFace = 0;
	for (int i = 0; i < nFaces; i++) {
		if (mesh.faces[i].edges.size() > nNodesPerFace) {
			nNodesPerFace = mesh.faces[i].edges.size();
		}
	}

	dVertexLon.Initialize(nFaces, nNodesPerFace);
	dVertexLat.Initialize(nFaces, nNodesPerFace);

	// Store coordinates of each Node and Face centerpoint
	for (int i = 0; i < nFaces; i++) {

		const Face & face = mesh.faces[i];

		int nNodes = face.edges.size();

		double dXc = 0.0;
		double dYc = 0.0;
		double dZc = 0.0;

		for (int j = 0; j < nNodes; j++) {
			const Node & node = mesh.nodes[face[j]];

			double dX = node.x;
			double dY = node.y;
			double dZ = node.z;

			dXc += dX;
			dYc += dY;
			dZc += dZ;

			double dLonV = atan2(dY, dX);
			double dLatV = acos(dZ);

			dVertexLon[i][j] = dLonV / M_PI * 180.0;
			dVertexLat[i][j] = dLatV / M_PI * 180.0;
		}

		dXc /= static_cast<double>(nNodes);
		dYc /= static_cast<double>(nNodes);
		dZc /= static_cast<double>(nNodes);

		double dMag = sqrt(dXc * dXc + dYc * dYc + dZc * dZc);

		dXc /= dMag;
		dYc /= dMag;
		dZc /= dMag;

		double dLonC = atan2(dYc, dXc);
		double dLatC = asin(dZc);

		dCenterLon[i] = dLonC / M_PI * 180.0;
		dCenterLat[i] = dLatC / M_PI * 180.0;
	}
}

///////////////////////////////////////////////////////////////////////////////

void OfflineMap::InitializeCoordinatesFromMeshFE(
	const Mesh & mesh,
	int nP,
	const DataMatrix3D<int> & dataGLLnodes,
	DataVector<double> & dCenterLon,
	DataVector<double> & dCenterLat,
	DataMatrix<double> & dVertexLon,
	DataMatrix<double> & dVertexLat
) {
	int nFaces = mesh.faces.size();

	if (nFaces != dataGLLnodes.GetSubColumns()) {
		_EXCEPTIONT("Mismatch between mesh and dataGLLnodes");
	}

	// Determine maximum index in dataGLLnodes
	int iMaxNodeIx = dataGLLnodes[0][0][0];
	for (int i = 0; i < dataGLLnodes.GetRows(); i++) {
	for (int j = 0; j < dataGLLnodes.GetColumns(); j++) {
	for (int k = 0; k < dataGLLnodes.GetSubColumns(); k++) {
		if (dataGLLnodes[i][j][k] > iMaxNodeIx) {
			iMaxNodeIx = dataGLLnodes[i][j][k];
		}
	}
	}
	}

	dCenterLon.Initialize(iMaxNodeIx);
	dCenterLat.Initialize(iMaxNodeIx);

	dVertexLon.Initialize(iMaxNodeIx, 1);
	dVertexLat.Initialize(iMaxNodeIx, 1);

	DataVector<double> dG;
	GetDefaultNodalLocations(nP, dG);

	for (int i = 0; i < dataGLLnodes.GetRows(); i++) {
	for (int j = 0; j < dataGLLnodes.GetColumns(); j++) {
	for (int k = 0; k < dataGLLnodes.GetSubColumns(); k++) {
		const Face & face = mesh.faces[k];

		Node node;

		ApplyLocalMap(
			face,
			mesh.nodes,
			dG[i],
			dG[j],
			node);

		int iNode = dataGLLnodes[j][i][k] - 1;

		double dLon = atan2(node.y, node.x);
		double dLat = asin(node.z);

		dCenterLon[iNode] = dLon / M_PI * 180.0;
		dCenterLat[iNode] = dLat / M_PI * 180.0;
	}
	}
	}
}

///////////////////////////////////////////////////////////////////////////////

void OfflineMap::InitializeSourceCoordinatesFromMeshFV(
	const Mesh & meshSource
) {
	InitializeCoordinatesFromMeshFV(
		meshSource,
		m_dSourceCenterLon,
		m_dSourceCenterLat,
		m_dSourceVertexLon,
		m_dSourceVertexLat);
}

///////////////////////////////////////////////////////////////////////////////

void OfflineMap::InitializeTargetCoordinatesFromMeshFV(
	const Mesh & meshTarget
) {
	InitializeCoordinatesFromMeshFV(
		meshTarget,
		m_dTargetCenterLon,
		m_dTargetCenterLat,
		m_dTargetVertexLon,
		m_dTargetVertexLat);
}

///////////////////////////////////////////////////////////////////////////////

void OfflineMap::InitializeSourceCoordinatesFromMeshFE(
	const Mesh & meshSource,
	int nP,
	const DataMatrix3D<int> & dataGLLnodesSource
) {
	InitializeCoordinatesFromMeshFE(
		meshSource,
		nP,
		dataGLLnodesSource,
		m_dSourceCenterLon,
		m_dSourceCenterLat,
		m_dSourceVertexLon,
		m_dSourceVertexLat);
}

///////////////////////////////////////////////////////////////////////////////

void OfflineMap::InitializeTargetCoordinatesFromMeshFE(
	const Mesh & meshTarget,
	int nP,
	const DataMatrix3D<int> & dataGLLnodesTarget
) {
	InitializeCoordinatesFromMeshFE(
		meshTarget,
		nP,
		dataGLLnodesTarget,
		m_dTargetCenterLon,
		m_dTargetCenterLat,
		m_dTargetVertexLon,
		m_dTargetVertexLat);
}

///////////////////////////////////////////////////////////////////////////////

NcDim * NcFile_GetDimIfExists(
	NcFile & ncSource,
	const std::string & strDimName,
	int nSize
) {
	NcDim * dim = ncSource.get_dim(strDimName.c_str());
	if (dim == NULL) {
		return ncSource.add_dim(strDimName.c_str(), nSize);
	}
	
	if (dim->size() != nSize) {
		_EXCEPTION3("NetCDF file has dimension \"%s\" with mismatched"
			" size %i != %i", strDimName.c_str(), dim->size(), nSize);
	}
	return dim;
}

///////////////////////////////////////////////////////////////////////////////

void OfflineMap::Apply(
	const std::string & strSourceDataFile,
	const std::string & strTargetDataFile,
	const std::vector<std::string> & vecVariables,
	const std::string & strNColName,
	bool fTargetDouble,
	bool fAppend
) {
	// Open source data file
	NcFile ncSource(strSourceDataFile.c_str(), NcFile::ReadOnly);

	// Open target data file
	NcFile::FileMode eOpenMode = NcFile::Replace;
	if (fAppend) {
		eOpenMode = NcFile::Write;
	}

	NcFile ncTarget(strTargetDataFile.c_str(), eOpenMode);
	if (!ncTarget.is_valid()) {
		_EXCEPTION1("Cannot open target data file \"%s\"",
			strTargetDataFile.c_str());
	}

	// Number of source and target regions
	int nSourceCount = m_dSourceAreas.GetRows();
	int nTargetCount = m_dTargetAreas.GetRows();

	// Check for rectilinear data
	bool fSourceRectilinear;
	if (m_vecSourceDimSizes.size() == 1) {
		fSourceRectilinear = false;
	} else if (m_vecSourceDimSizes.size() == 2) {
		fSourceRectilinear = true;
	} else {
		_EXCEPTIONT("m_vecSourceDimSizes undefined");
	}

	if (fSourceRectilinear) {
		if (nSourceCount != m_vecSourceDimSizes[0] * m_vecSourceDimSizes[1]) {
			_EXCEPTIONT("Rectilinear input expected in finite volume form");
		}
	}

	bool fTargetRectilinear;
	if (m_vecTargetDimSizes.size() == 1) {
		fTargetRectilinear = false;
	} else if (m_vecTargetDimSizes.size() == 2) {
		fTargetRectilinear = true;
	} else {
		_EXCEPTIONT("m_vecTargetDimSizes undefined");
	}

	if (fTargetRectilinear) {
		if (nTargetCount != m_vecTargetDimSizes[0] * m_vecTargetDimSizes[1]) {
			printf("%i %i\n", nTargetCount,
				m_vecTargetDimSizes[0] * m_vecTargetDimSizes[1]);

			_EXCEPTIONT("Rectilinear output expected in finite volume form");
		}
	}

	DataVector<float> dataIn;
	if (m_vecSourceDimSizes.size() == 1) {
		dataIn.Initialize(nSourceCount);
	} else {
		dataIn.Initialize(m_vecSourceDimSizes[0] * m_vecSourceDimSizes[1]);
	}

	DataVector<float> dataOut;
	if (m_vecTargetDimSizes.size() == 1) {
		dataOut.Initialize(nTargetCount);
	} else {
		dataOut.Initialize(m_vecTargetDimSizes[0] * m_vecTargetDimSizes[1]);
	}

	DataVector<double> dataInDouble;
	dataInDouble.Initialize(nSourceCount);

	DataVector<double> dataOutDouble;
	dataOutDouble.Initialize(nTargetCount);

	// Target
	if (!fAppend) {
		CopyNcFileAttributes(&ncSource, &ncTarget);
	}

	NcDim * dim0;
	NcDim * dim1;

	if (!fTargetRectilinear) {
		dim0 = NcFile_GetDimIfExists(
			ncTarget,
			strNColName.c_str(),
			nTargetCount);

	} else {
		dim0 = NcFile_GetDimIfExists(
			ncTarget,
			m_vecTargetDimNames[0].c_str(),
			m_vecTargetDimSizes[0]);

		dim1 = NcFile_GetDimIfExists(
			ncTarget,
			m_vecTargetDimNames[1].c_str(),
			m_vecTargetDimSizes[1]);
	}

	// Loop through all variables
	for (int v = 0; v < vecVariables.size(); v++) {
		NcVar * var = ncSource.get_var(vecVariables[v].c_str());
		if (var == NULL) {
			_EXCEPTION1("Variable \"%s\" does not exist in source file",
				vecVariables[v].c_str());
		}

		AnnounceStartBlock(vecVariables[v].c_str());

		// Check for _FillValue
		float flFillValue = 0.0f;
		double dFillValue = 0.0;
		for (int a = 0; a < var->num_atts(); a++) {
			NcAtt * att = var->get_att(a);
			if (strcmp(att->name(), "_FillValue") == 0) {
				if (att->type() == ncDouble) {
					dFillValue = att->as_double(0);
				} else if (att->type() == ncFloat) {
					flFillValue = att->as_float(0);
				} else {
					_EXCEPTIONT("Invalid type for _FillValue");
				}
			}
		}

		// Loop through all dimensions and add missing dimensions to output
		int nVarTotalEntries = 1;

		DataVector<NcDim *> vecDims;
		vecDims.Initialize(var->num_dims());

		DataVector<NcDim *> vecDimsOut;
		if (fTargetRectilinear) {
			if (fSourceRectilinear) {
				vecDimsOut.Initialize(var->num_dims());
			} else {
				vecDimsOut.Initialize(var->num_dims()+1);
			}

			vecDimsOut[vecDimsOut.GetRows()-2] = dim0;
			vecDimsOut[vecDimsOut.GetRows()-1] = dim1;
		} else {
			if (fSourceRectilinear) {
				vecDimsOut.Initialize(var->num_dims()-1);
			} else {
				vecDimsOut.Initialize(var->num_dims());
			}

			vecDimsOut[vecDimsOut.GetRows()-1] = dim0;
		}

		int nFreeDims = var->num_dims() - m_vecSourceDimSizes.size();

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
					ncTarget,
					strDimName.c_str(),
					nDimSize);
		}

		// Create new output variable
		NcVar * varOut;
		if (fTargetDouble) {
			varOut =
				ncTarget.add_var(
					vecVariables[v].c_str(),
					ncDouble,
					vecDimsOut.GetRows(),
					(const NcDim**)&(vecDimsOut[0]));

		} else {
			varOut =
				ncTarget.add_var(
					vecVariables[v].c_str(),
					ncFloat,
					vecDimsOut.GetRows(),
					(const NcDim**)&(vecDimsOut[0]));
		}

		CopyNcVarAttributes(var, varOut);

		// Source and output counts
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
		if (fSourceRectilinear) {
			nGet[nGet.GetRows()-2] = m_vecSourceDimSizes[0];
			nGet[nGet.GetRows()-1] = m_vecSourceDimSizes[1];
		} else {
			nGet[nGet.GetRows()-1] = nSourceCount;
		}

		// Put size
		DataVector<long> nPut;
		nPut.Initialize(nCountsOut.GetRows());
		for (int d = 0; d < nPut.GetRows()-1; d++) {
			nPut[d] = 1;
		}
		if (fTargetRectilinear) {
			nPut[nPut.GetRows()-2] = m_vecTargetDimSizes[0];
			nPut[nPut.GetRows()-1] = m_vecTargetDimSizes[1];
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

				if (flFillValue != 0.0f) {
					for (int i = 0; i < nSourceCount; i++) {
						if (dataIn[i] == flFillValue) {
							dataInDouble[i] = 0.0;
						} else {
							dataInDouble[i] = static_cast<double>(dataIn[i]);
						}
					}

				} else {
					for (int i = 0; i < nSourceCount; i++) {
						dataInDouble[i] = static_cast<double>(dataIn[i]);
					}
				}

			// Load data as Double
			} else if (var->type() == ncDouble) {
				var->get(&(dataInDouble[0]), &(nGet[0]));

				if (dFillValue != 0.0) {
					for (int i = 0; i < nSourceCount; i++) {
						if (dataInDouble[i] == dFillValue) {
							dataInDouble[i] = 0.0;
						}
					}
				}

			} else {
				_EXCEPTIONT("Invalid variable type");
			}

			// Announce input mass
			double dSourceMass = 0.0;
			double dSourceMin  = dataInDouble[0];
			double dSourceMax  = dataInDouble[0];
			for (int i = 0; i < nSourceCount; i++) {
				dSourceMass += dataInDouble[i] * m_dSourceAreas[i];
				if (dataInDouble[i] < dSourceMin) {
					dSourceMin = dataInDouble[i];
				}
				if (dataInDouble[i] > dSourceMax) {
					dSourceMax = dataInDouble[i];
				}
			}
			Announce("Source Mass: %1.15e Min %1.10e Max %1.10e",
				dSourceMass, dSourceMin, dSourceMax);

			// Apply the offline map to the data
			m_mapRemap.Apply(dataInDouble, dataOutDouble);

			// Announce output mass
			double dTargetMass = 0.0;
			double dTargetMin  = dataOutDouble[0];
			double dTargetMax  = dataOutDouble[0];
			for (int i = 0; i < nTargetCount; i++) {
				dTargetMass += dataOutDouble[i] * m_dTargetAreas[i];
				if (dataOutDouble[i] < dTargetMin) {
					dTargetMin = dataOutDouble[i];
				}
				if (dataOutDouble[i] > dTargetMax) {
					dTargetMax = dataOutDouble[i];
				}
			}
			Announce("Target Mass: %1.15e Min %1.10e Max %1.10e",
				dTargetMass, dTargetMin, dTargetMax);

			// Write the data
			if (fTargetDouble) {
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
	const std::string & strSource
) {
	NcFile ncMap(strSource.c_str(), NcFile::ReadOnly);
	if (!ncMap.is_valid()) {
		_EXCEPTION1("Unable to open input map file \"%s\"",
			strSource.c_str());
	}

	// Read input dimensions entries
	NcDim * dimSrcGridRank = ncMap.get_dim("src_grid_rank");
	if (dimSrcGridRank == NULL) {
		_EXCEPTION1("Map file \"%s\" does not contain variable"
			"\"src_grid_rank\"", strSource.c_str());
	}

	NcDim * dimDstGridRank = ncMap.get_dim("dst_grid_rank");
	if (dimDstGridRank == NULL) {
		_EXCEPTION1("Map file \"%s\" does not contain variable"
			"\"dst_grid_rank\"", strSource.c_str());
	}

	int nSrcGridDims = (int)(dimSrcGridRank->size());
	int nDstGridDims = (int)(dimDstGridRank->size());

	NcVar * varSrcGridDims = ncMap.get_var("src_grid_dims");
	if (varSrcGridDims == NULL) {
		_EXCEPTION1("Map file \"%s\" does not contain variable"
			"\"src_grid_dims\"", strSource.c_str());
	}

	NcVar * varDstGridDims = ncMap.get_var("dst_grid_dims");
	if (varSrcGridDims == NULL) {
		_EXCEPTION1("Map file \"%s\" does not contain variable"
			"\"dst_grid_dims\"", strSource.c_str());
	}

	m_vecSourceDimSizes.resize(nSrcGridDims);
	m_vecSourceDimNames.resize(nSrcGridDims);

	m_vecTargetDimSizes.resize(nDstGridDims);
	m_vecTargetDimNames.resize(nDstGridDims);

	varSrcGridDims->get(&(m_vecSourceDimSizes[0]), nSrcGridDims);
	varDstGridDims->get(&(m_vecTargetDimSizes[0]), nDstGridDims);

	for (int i = 0; i < nSrcGridDims; i++) {
		char szDim[64];
		sprintf(szDim, "name%i", i);
		m_vecSourceDimNames[i] = varSrcGridDims->get_att(szDim)->as_string(0);
	}

	for (int i = 0; i < nDstGridDims; i++) {
		char szDim[64];
		sprintf(szDim, "name%i", i);
		m_vecTargetDimNames[i] = varDstGridDims->get_att(szDim)->as_string(0);
	}

	// Source and Target mesh resolutions
	NcDim * dimNA = ncMap.get_dim("n_a");
	if (dimNA == NULL) {
		_EXCEPTIONT("Input map missing dimension \"n_a\"");
	}

	NcDim * dimNB = ncMap.get_dim("n_b");
	if (dimNB == NULL) {
		_EXCEPTIONT("Input map missing dimension \"n_b\"");
	}

	int nA = dimNA->size();
	int nB = dimNB->size();

	NcDim * dimNVA = ncMap.get_dim("nv_a");
	if (dimNA == NULL) {
		_EXCEPTIONT("Input map missing dimension \"nv_a\"");
	}

	NcDim * dimNVB = ncMap.get_dim("nv_b");
	if (dimNB == NULL) {
		_EXCEPTIONT("Input map missing dimension \"nv_b\"");
	}

	int nVA = dimNVA->size();
	int nVB = dimNVB->size();

	// Read coordinates
	NcVar * varYCA = ncMap.get_var("yc_a");
	NcVar * varYCB = ncMap.get_var("yc_b");

	NcVar * varXCA = ncMap.get_var("xc_a");
	NcVar * varXCB = ncMap.get_var("xc_b");

	NcVar * varYVA = ncMap.get_var("yv_a");
	NcVar * varYVB = ncMap.get_var("yv_b");

	NcVar * varXVA = ncMap.get_var("xv_a");
	NcVar * varXVB = ncMap.get_var("xv_b");

	if (varYCA == NULL) {
		_EXCEPTION1("Map file \"%s\" does not contain variable \"yc_a\":"
			"\nPossibly an earlier version of map",
			strSource.c_str());
	}
	if (varYCB == NULL) {
		_EXCEPTION1("Map file \"%s\" does not contain variable \"yc_b\":"
			"\nPossibly an earlier version of map",
			strSource.c_str());
	}
	if (varXCA == NULL) {
		_EXCEPTION1("Map file \"%s\" does not contain variable \"xc_a\":"
			"\nPossibly an earlier version of map",
			strSource.c_str());
	}
	if (varXCB == NULL) {
		_EXCEPTION1("Map file \"%s\" does not contain variable \"xc_b\":"
			"\nPossibly an earlier version of map",
			strSource.c_str());
	}
	if (varYVA == NULL) {
		_EXCEPTION1("Map file \"%s\" does not contain variable \"yv_a\":"
			"\nPossibly an earlier version of map",
			strSource.c_str());
	}
	if (varYVB == NULL) {
		_EXCEPTION1("Map file \"%s\" does not contain variable \"yv_b\":"
			"\nPossibly an earlier version of map",
			strSource.c_str());
	}
	if (varXVA == NULL) {
		_EXCEPTION1("Map file \"%s\" does not contain variable \"xv_a\":"
			"\nPossibly an earlier version of map",
			strSource.c_str());
	}
	if (varXVB == NULL) {
		_EXCEPTION1("Map file \"%s\" does not contain variable \"xv_b\":"
			"\nPossibly an earlier version of map",
			strSource.c_str());
	}

	m_dSourceCenterLat.Initialize(nA);
	m_dTargetCenterLat.Initialize(nB);

	m_dSourceCenterLon.Initialize(nA);
	m_dTargetCenterLon.Initialize(nB);

	m_dSourceVertexLat.Initialize(nA, nVA);
	m_dTargetVertexLat.Initialize(nB, nVB);

	m_dSourceVertexLon.Initialize(nA, nVA);
	m_dTargetVertexLon.Initialize(nB, nVB);

	varYCA->get(&(m_dSourceCenterLat[0]), nA);
	varYCB->get(&(m_dTargetCenterLat[0]), nB);

	varXCA->get(&(m_dSourceCenterLon[0]), nA);
	varXCB->get(&(m_dTargetCenterLon[0]), nB);

	varYVA->get(&(m_dSourceVertexLat[0][0]), nA, nVA);
	varYVB->get(&(m_dTargetVertexLat[0][0]), nB, nVB);

	varXVA->get(&(m_dSourceVertexLon[0][0]), nA, nVA);
	varXVB->get(&(m_dTargetVertexLon[0][0]), nB, nVB);

	// Read areas
	m_dSourceAreas.Initialize(nA);
	m_dTargetAreas.Initialize(nB);

	NcVar * varAreaA = ncMap.get_var("area_a");
	if (varAreaA == NULL) {
		_EXCEPTION1("Map file \"%s\" does not contain variable \"area_a\"",
			strSource.c_str());
	}
	varAreaA->get(&(m_dSourceAreas[0]), nA);

	NcVar * varAreaB = ncMap.get_var("area_b");
	if (varAreaB == NULL) {
		_EXCEPTION1("Map file \"%s\" does not contain variable \"area_b\"",
			strSource.c_str());
	}
	varAreaB->get(&(m_dTargetAreas[0]), nB);

	// Read SparseMatrix entries
	NcDim * dimNS = ncMap.get_dim("n_s");
	if (dimNS == NULL) {
		_EXCEPTION1("Map file \"%s\" does not contain dimension \"n_s\"",
			strSource.c_str());
	}

	NcVar * varRow = ncMap.get_var("row");
	if (varRow == NULL) {
		_EXCEPTION1("Map file \"%s\" does not contain variable \"row\"",
			strSource.c_str());
	}

	NcVar * varCol = ncMap.get_var("col");
	if (varRow == NULL) {
		_EXCEPTION1("Map file \"%s\" does not contain variable \"col\"",
			strSource.c_str());
	}

	NcVar * varS = ncMap.get_var("S");
	if (varRow == NULL) {
		_EXCEPTION1("Map file \"%s\" does not contain variable \"S\"",
			strSource.c_str());
	}

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

	// Decrement vecRow and vecCol
	for (int i = 0; i < vecRow.GetRows(); i++) {
		vecRow[i]--;
		vecCol[i]--;
	}

	// Set the entries of the map
	m_mapRemap.SetEntries(vecRow, vecCol, vecS);
}

///////////////////////////////////////////////////////////////////////////////

void OfflineMap::Write(
	const std::string & strTarget
) {
	NcFile ncMap(strTarget.c_str(), NcFile::Replace);
	if (!ncMap.is_valid()) {
		_EXCEPTION1("Unable to open output map file \"%s\"",
			strTarget.c_str());
	}

	// Attributes
	ncMap.add_att("Title", "TempestRemap Offline Regridding Weight Generator");

	// Write output dimensions entries
	int nSrcGridDims = (int)(m_vecSourceDimSizes.size());
	int nDstGridDims = (int)(m_vecTargetDimSizes.size());

	NcDim * dimSrcGridRank = ncMap.add_dim("src_grid_rank", nSrcGridDims);
	NcDim * dimDstGridRank = ncMap.add_dim("dst_grid_rank", nDstGridDims);

	NcVar * varSrcGridDims =
		ncMap.add_var("src_grid_dims", ncInt, dimSrcGridRank);
	NcVar * varDstGridDims =
		ncMap.add_var("dst_grid_dims", ncInt, dimDstGridRank);

	varSrcGridDims->put(&(m_vecSourceDimSizes[0]), nSrcGridDims);
	varDstGridDims->put(&(m_vecTargetDimSizes[0]), nDstGridDims);

	int nA = (int)(m_dSourceAreas.GetRows());
	int nB = (int)(m_dTargetAreas.GetRows());

	for (int i = 0; i < m_vecSourceDimSizes.size(); i++) {
		char szDim[64];
		sprintf(szDim, "name%i", i);
		varSrcGridDims->add_att(szDim, m_vecSourceDimNames[i].c_str());
	}

	for (int i = 0; i < m_vecTargetDimSizes.size(); i++) {
		char szDim[64];
		sprintf(szDim, "name%i", i);
		varDstGridDims->add_att(szDim, m_vecTargetDimNames[i].c_str());
	}

	// Source and Target mesh resolutions
	NcDim * dimNA = ncMap.add_dim("n_a", nA);
	NcDim * dimNB = ncMap.add_dim("n_b", nB);

	// Number of nodes per Face
	int nSourceNodesPerFace = m_dSourceVertexLon.GetColumns();
	int nTargetNodesPerFace = m_dTargetVertexLon.GetColumns();

	NcDim * dimNVA = ncMap.add_dim("nv_a", nSourceNodesPerFace);
	NcDim * dimNVB = ncMap.add_dim("nv_b", nTargetNodesPerFace);

	// Write coordinates
	NcVar * varYCA = ncMap.add_var("yc_a", ncDouble, dimNA);
	NcVar * varYCB = ncMap.add_var("yc_b", ncDouble, dimNB);

	NcVar * varXCA = ncMap.add_var("xc_a", ncDouble, dimNA);
	NcVar * varXCB = ncMap.add_var("xc_b", ncDouble, dimNB);

	NcVar * varYVA = ncMap.add_var("yv_a", ncDouble, dimNA, dimNVA);
	NcVar * varYVB = ncMap.add_var("yv_b", ncDouble, dimNB, dimNVB);

	NcVar * varXVA = ncMap.add_var("xv_a", ncDouble, dimNA, dimNVA);
	NcVar * varXVB = ncMap.add_var("xv_b", ncDouble, dimNB, dimNVB);

	varYCA->add_att("units", "degrees");
	varYCB->add_att("units", "degrees");

	varXCA->add_att("units", "degrees");
	varXCB->add_att("units", "degrees");

	varYVA->add_att("units", "degrees");
	varYVB->add_att("units", "degrees");

	varXVA->add_att("units", "degrees");
	varXVB->add_att("units", "degrees");

	// Verify dimensionality
	if (m_dSourceCenterLon.GetRows() != nA) {
		_EXCEPTIONT("Mismatch between m_dSourceCenterLon and nA");
	}
	if (m_dSourceCenterLat.GetRows() != nA) {
		_EXCEPTIONT("Mismatch between m_dSourceCenterLat and nA");
	}
	if (m_dTargetCenterLon.GetRows() != nB) {
		_EXCEPTIONT("Mismatch between m_dTargetCenterLon and nB");
	}
	if (m_dTargetCenterLat.GetRows() != nB) {
		_EXCEPTIONT("Mismatch between m_dTargetCenterLat and nB");
	}
	if (m_dSourceVertexLon.GetRows() != nA) {
		_EXCEPTIONT("Mismatch between m_dSourceVertexLon and nA");
	}
	if (m_dSourceVertexLat.GetRows() != nA) {
		_EXCEPTIONT("Mismatch between m_dSourceVertexLat and nA");
	}
	if (m_dTargetVertexLon.GetRows() != nB) {
		_EXCEPTIONT("Mismatch between m_dTargetVertexLon and nB");
	}
	if (m_dTargetVertexLat.GetRows() != nB) {
		_EXCEPTIONT("Mismatch between m_dTargetVertexLat and nB");
	}

	varYCA->put(&(m_dSourceCenterLat[0]), nA);
	varYCB->put(&(m_dTargetCenterLat[0]), nB);

	varXCA->put(&(m_dSourceCenterLon[0]), nA);
	varXCB->put(&(m_dTargetCenterLon[0]), nB);

	varYVA->put(&(m_dSourceVertexLat[0][0]), nA, nSourceNodesPerFace);
	varYVB->put(&(m_dTargetVertexLat[0][0]), nB, nTargetNodesPerFace);

	varXVA->put(&(m_dSourceVertexLon[0][0]), nA, nSourceNodesPerFace);
	varXVB->put(&(m_dTargetVertexLon[0][0]), nB, nTargetNodesPerFace);

	// Write areas
	NcVar * varAreaA = ncMap.add_var("area_a", ncDouble, dimNA);
	varAreaA->put(&(m_dSourceAreas[0]), nA);

	NcVar * varAreaB = ncMap.add_var("area_b", ncDouble, dimNB);
	varAreaB->put(&(m_dTargetAreas[0]), nB);

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

	// Increment vecRow and vecCol
	for (int i = 0; i < vecRow.GetRows(); i++) {
		vecRow[i]++;
		vecCol[i]++;
	}

	// Load in data
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
	if (vecSourceAreas.GetRows() != m_mapRemap.GetColumns()) {
		_EXCEPTIONT("vecSourceAreas / mapRemap dimension mismatch");
	}
	if (vecTargetAreas.GetRows() != m_mapRemap.GetRows()) {
		_EXCEPTIONT("vecTargetAreas / mapRemap dimension mismatch");
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
			dataEntries[i] * m_dTargetAreas[dataRows[i]];
	}

	// Verify all column sums equal the input Jacobian
	bool fConservative = true;
	for (int i = 0; i < dColumnSums.GetRows(); i++) {
		if (fabs(dColumnSums[i] - m_dSourceAreas[i]) > dTolerance) {
			fConservative = false;
			Announce("OfflineMap is not conservative in column "
				"%i (%1.15e / %1.15e)",
				i, dColumnSums[i], m_dSourceAreas[i]);
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

