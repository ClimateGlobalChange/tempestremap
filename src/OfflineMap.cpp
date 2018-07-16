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
			m_vecSourceDimNames.push_back("lat");
			m_vecSourceDimNames.push_back("lon");

			int iTemp = m_vecSourceDimSizes[0];
			m_vecSourceDimSizes[0] = m_vecSourceDimSizes[1];
			m_vecSourceDimSizes[1] = iTemp;
		} else {
			_EXCEPTIONT("Source grid grid_rank must be < 3");
		}

		// Number of faces
		NcDim * dimGridSize = ncSourceMesh.get_dim("grid_size");
		if (dimGridSize == NULL) {
			_EXCEPTIONT("Missing \"grid_size\" dimension in grid file");
		}

		// Number of grid corners
		NcDim * dimGridCorners = ncSourceMesh.get_dim("grid_corners");
		if (dimGridCorners == NULL) {
			_EXCEPTIONT("Missing \"grid_corners\" dimension in grid file");
		}

		// Pull grid center information from file
		NcVar * varGridCenterLon = ncSourceMesh.get_var("grid_center_lon");
		if (varGridCenterLon == NULL) {
			_EXCEPTIONT("Missing \"grid_center_lon\" variable in grid file");
		}

		m_dSourceCenterLon.Initialize(dimGridSize->size());

		varGridCenterLon->get(
			&(m_dSourceCenterLon[0]),
			dimGridSize->size());

		NcVar * varGridCenterLat = ncSourceMesh.get_var("grid_center_lat");
		if (varGridCenterLat == NULL) {
			_EXCEPTIONT("Missing \"grid_center_lat\" variable in grid file");
		}

		m_dSourceCenterLat.Initialize(dimGridSize->size());

		varGridCenterLat->get(
			&(m_dSourceCenterLat[0]),
			dimGridSize->size());

		// Pull grid vertex information from file
		NcVar * varGridVertexLon = ncSourceMesh.get_var("grid_corner_lon");
		if (varGridVertexLon == NULL) {
			_EXCEPTIONT("Missing \"grid_corner_lon\" variable in grid file");
		}

		m_dSourceVertexLon.Initialize(
			dimGridSize->size(),
			dimGridCorners->size());

		varGridVertexLon->get(
			&(m_dSourceVertexLon[0][0]),
			dimGridSize->size(),
			dimGridCorners->size());

		NcVar * varGridVertexLat = ncSourceMesh.get_var("grid_corner_lat");
		if (varGridVertexLat == NULL) {
			_EXCEPTIONT("Missing \"grid_corner_lat\" variable in grid file");
		}

		m_dSourceVertexLat.Initialize(
			dimGridSize->size(),
			dimGridCorners->size());

		varGridVertexLat->get(
			&(m_dSourceVertexLat[0][0]),
			dimGridSize->size(),
			dimGridCorners->size());

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

void OfflineMap::InitializeSourceDimensions(
  const std::vector<std::string>& p_srcDimNames,
  const std::vector<int>& p_srcDimSizes
) {
  m_vecSourceDimNames.clear();
  m_vecSourceDimNames.resize(p_srcDimNames.size());
  std::copy(p_srcDimNames.begin(), p_srcDimNames.end(), m_vecSourceDimNames.begin());
  
  m_vecSourceDimSizes.clear();
  m_vecSourceDimSizes.resize(p_srcDimSizes.size());
  std::copy(p_srcDimSizes.begin(), p_srcDimSizes.end(), m_vecSourceDimSizes.begin());
  
  return;
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
		if (dimGridRank == NULL) {
			_EXCEPTIONT("Variable \"grid_dims\" has no dimensions");
		}

		m_vecTargetDimSizes.resize(dimGridRank->size());
		varGridDims->get(&(m_vecTargetDimSizes[0]), dimGridRank->size());

		//std::cout << "TEST " << dimGridRank->size() << std::endl;
		if (dimGridRank->size() == 1) {
			m_vecTargetDimNames.push_back("num_elem");
		} else if (dimGridRank->size() == 2) {
			m_vecTargetDimNames.push_back("lat");
			m_vecTargetDimNames.push_back("lon");

			int iTemp = m_vecTargetDimSizes[0];
			m_vecTargetDimSizes[0] = m_vecTargetDimSizes[1];
			m_vecTargetDimSizes[1] = iTemp;

		} else {
			_EXCEPTIONT("Target grid grid_rank must be < 3");
		}

		std::cout << "TEST " << m_vecTargetDimSizes.size() << std::endl;

		// Number of faces
		NcDim * dimGridSize = ncTargetMesh.get_dim("grid_size");
		if (dimGridSize == NULL) {
			_EXCEPTIONT("Missing \"grid_size\" dimension in grid file");
		}

		// Number of grid corners
		NcDim * dimGridCorners = ncTargetMesh.get_dim("grid_corners");
		if (dimGridCorners == NULL) {
			_EXCEPTIONT("Missing \"grid_corners\" dimension in grid file");
		}

		// Pull grid center information from file
		NcVar * varGridCenterLon = ncTargetMesh.get_var("grid_center_lon");
		if (varGridCenterLon == NULL) {
			_EXCEPTIONT("Missing \"grid_center_lon\" variable in grid file");
		}

		m_dTargetCenterLon.Initialize(dimGridSize->size());

		varGridCenterLon->get(
			&(m_dTargetCenterLon[0]),
			dimGridSize->size());

		NcVar * varGridCenterLat = ncTargetMesh.get_var("grid_center_lat");
		if (varGridCenterLat == NULL) {
			_EXCEPTIONT("Missing \"grid_center_lat\" variable in grid file");
		}

		m_dTargetCenterLat.Initialize(dimGridSize->size());

		varGridCenterLat->get(
			&(m_dTargetCenterLat[0]),
			dimGridSize->size());

		// Pull grid vertex information from file
		NcVar * varGridVertexLon = ncTargetMesh.get_var("grid_corner_lon");
		if (varGridVertexLon == NULL) {
			_EXCEPTIONT("Missing \"grid_corner_lon\" variable in grid file");
		}

		m_dTargetVertexLon.Initialize(
			dimGridSize->size(),
			dimGridCorners->size());

		varGridVertexLon->get(
			&(m_dTargetVertexLon[0][0]),
			dimGridSize->size(),
			dimGridCorners->size());

		NcVar * varGridVertexLat = ncTargetMesh.get_var("grid_corner_lat");
		if (varGridVertexLat == NULL) {
			_EXCEPTIONT("Missing \"grid_corner_lat\" variable in grid file");
		}

		m_dTargetVertexLat.Initialize(
			dimGridSize->size(),
			dimGridCorners->size());

		varGridVertexLat->get(
			&(m_dTargetVertexLat[0][0]),
			dimGridSize->size(),
			dimGridCorners->size());

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
	m_vecTargetDimSizes[0] = nDim0Size;
	m_vecTargetDimSizes[1] = nDim1Size;

	m_vecTargetDimNames.resize(2);
	m_vecTargetDimNames[0] = strDim0Name;
	m_vecTargetDimNames[1] = strDim1Name;
}

///////////////////////////////////////////////////////////////////////////////

void OfflineMap::InitializeTargetDimensions(
  const std::vector<std::string>& p_tgtDimNames,
  const std::vector<int>& p_tgtDimSizes
) {
  m_vecTargetDimNames.clear();
  m_vecTargetDimNames.resize(p_tgtDimNames.size());
  std::copy(p_tgtDimNames.begin(), p_tgtDimNames.end(), m_vecTargetDimNames.begin());
  
  m_vecTargetDimSizes.clear();
  m_vecTargetDimSizes.resize(p_tgtDimSizes.size());
  std::copy(p_tgtDimSizes.begin(), p_tgtDimSizes.end(), m_vecTargetDimSizes.begin());

  return;
}

///////////////////////////////////////////////////////////////////////////////

void OfflineMap::InitializeCoordinatesFromMeshFV(
	const Mesh & mesh,
	DataVector<double> & dCenterLon,
	DataVector<double> & dCenterLat,
	DataMatrix<double> & dVertexLon,
	DataMatrix<double> & dVertexLat,
	bool fLatLon
) {
	int nFaces = mesh.faces.size();

	// Count maximum number of Nodes per Face
	int nNodesPerFace = 0;
	for (int i = 0; i < nFaces; i++) {
		if (mesh.faces[i].edges.size() > nNodesPerFace) {
			nNodesPerFace = mesh.faces[i].edges.size();
		}
	}

	// Check if already initialized
	if (dCenterLon.GetRows() != 0) {
		return;
	}

	dVertexLon.Initialize(nFaces, nNodesPerFace);
	dVertexLat.Initialize(nFaces, nNodesPerFace);

	if ((fLatLon) && (nNodesPerFace != 4)) {
		_EXCEPTIONT("Logic error");
	}

	dCenterLon.Initialize(nFaces);
	dCenterLat.Initialize(nFaces);

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

			XYZtoRLL_Deg(
				dX, dY, dZ,
				dVertexLon[i][j],
				dVertexLat[i][j]);
		}

		if ((fLatLon) && (nNodes == 3)) {
			dVertexLon[i][3] = dVertexLon[i][0];
			dVertexLat[i][3] = dVertexLat[i][0];
		}

		dXc /= static_cast<double>(nNodes);
		dYc /= static_cast<double>(nNodes);
		dZc /= static_cast<double>(nNodes);

		double dMag = sqrt(dXc * dXc + dYc * dYc + dZc * dZc);

		dXc /= dMag;
		dYc /= dMag;
		dZc /= dMag;

		XYZtoRLL_Deg(
			dXc,
			dYc,
			dZc,
			dCenterLon[i],
			dCenterLat[i]);

		// Modify vertex coordinates of polar volumes on latlon grid
		if (fLatLon) {

			// Change longitudes of polar volumes
			int nNodesMod = dVertexLat.GetColumns();
			for (int j = 0; j < nNodesMod; j++) {
				if (fabs(fabs(dVertexLat[i][j]) - 90.0) < 1.0e-12) {
					int jn = (j + 1) % nNodesMod;
					int jp = (j + nNodesMod - 1) % nNodesMod;

					if (fabs(fabs(dVertexLat[i][jn]) - 90.0) > 1.0e-12) {
						dVertexLon[i][j] = dVertexLon[i][jn];
					} else if (fabs(fabs(dVertexLat[i][jp]) - 90.0) > 1.0e-12) {
						dVertexLon[i][j] = dVertexLon[i][jp];
					} else {
						_EXCEPTIONT("Logic error");
					}
				}
			}
		}
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

		XYZtoRLL_Deg(
			node.x,
			node.y,
			node.z,
			dCenterLon[iNode],
			dCenterLat[iNode]);
	}
	}
	}
}

///////////////////////////////////////////////////////////////////////////////

void OfflineMap::InitializeRectilinearCoordinateVector(
	int nLon,
	int nLat,
	const DataMatrix<double> & dVertexLon,
	const DataMatrix<double> & dVertexLat,
	bool fLonFirst,
	DataVector<double> & dCenterLon,
	DataVector<double> & dCenterLat,
	DataVector<double> & dVectorCenterLon,
	DataVector<double> & dVectorCenterLat,
	DataMatrix<double> & dVectorBoundsLon,
	DataMatrix<double> & dVectorBoundsLat
) {
	if (dCenterLon.GetRows() != nLon * nLat) {
		_EXCEPTION3("CenterLon has incorrect dimension:\n"
			"\t[%i x %i] expected, [%i] found",
			nLon, nLat, dCenterLon.GetRows());
	}
	if (dCenterLat.GetRows() != nLon * nLat) {
		_EXCEPTION3("CenterLat has incorrect dimension:\n"
			"\t[%i x %i] expected, [%i] found",
			nLon, nLat, dCenterLat.GetRows());
	}

	const int nFaces = dVertexLon.GetRows();

	if (nFaces != nLon * nLat) {
		_EXCEPTIONT("Number of faces must be tensor product of nLon and nLat");
	}

	dCenterLon.Initialize(nFaces);
	dCenterLat.Initialize(nFaces);

	dVectorCenterLon.Initialize(nLon);
	dVectorCenterLat.Initialize(nLat);

	dVectorBoundsLon.Initialize(nLon, 2);
	dVectorBoundsLat.Initialize(nLat, 2);

	if (!fLonFirst) {
		for (int i = 0; i < nLon; i++) {
			dVectorBoundsLon[i][0] = 720.0;
			dVectorBoundsLon[i][1] = -720.0;
			for (int k = 0; k < dVertexLon.GetColumns(); k++) {
				if (fabs(fabs(dVertexLat[i][k]) - 90.0) < 1.0e-12) {
					continue;
				}
				if (dVertexLon[i][k] < dVectorBoundsLon[i][0]) {
					dVectorBoundsLon[i][0] = dVertexLon[i][k];
				}
				if (dVertexLon[i][k] > dVectorBoundsLon[i][1]) {
					dVectorBoundsLon[i][1] = dVertexLon[i][k];
				}
			}

			if ((dVectorBoundsLon[i][0] < 90.0) &&
			    (dVectorBoundsLon[i][1] > 270.0)
			) {
				if (i == 0) {
					double dTemp = dVectorBoundsLon[i][0];
					dVectorBoundsLon[i][0] = dVectorBoundsLon[i][1] - 360.0;
					dVectorBoundsLon[i][1] = dTemp;
				} else {
					double dTemp = dVectorBoundsLon[i][1];
					dVectorBoundsLon[i][1] = dVectorBoundsLon[i][0] + 360.0;
					dVectorBoundsLon[i][0] = dTemp;
				}
			}

			dVectorCenterLon[i] = 0.5 * (
				  dVectorBoundsLon[i][0]
				+ dVectorBoundsLon[i][1]);
		}
		for (int j = 0; j < nLat; j++) {
			dVectorBoundsLat[j][0] = dVertexLat[j * nLon][0];
			dVectorBoundsLat[j][1] = dVertexLat[j * nLon][0];
			for (int k = 0; k < dVertexLat.GetColumns(); k++) {
				if (dVertexLat[j * nLon][k] < dVectorBoundsLat[j][0]) {
					dVectorBoundsLat[j][0] = dVertexLat[j * nLon][k];
				}
				if (dVertexLat[j * nLon][k] > dVectorBoundsLat[j][1]) {
					dVectorBoundsLat[j][1] = dVertexLat[j * nLon][k];
				}
			}
			dVectorCenterLat[j] = 0.5 * (
				  dVectorBoundsLat[j][0]
				+ dVectorBoundsLat[j][1]);
		}

		for (int i = 0; i < nFaces; i++) {
			dCenterLon[i] = dVectorCenterLon[i%nLon];
			dCenterLat[i] = dVectorCenterLat[i/nLon];
		}

	} else {
		for (int i = 0; i < nLon; i++) {
			//dVectorCenterLon[i] = dCenterLon[i * nLat];
			dVectorBoundsLon[i][0] = 720.0;
			dVectorBoundsLon[i][1] = -720.0;
			for (int k = 0; k < dVertexLon.GetColumns(); k++) {
				if (fabs(fabs(dVertexLat[i * nLat][k]) - 90.0) < 1.0e-12) {
					continue;
				}
				if (dVertexLon[i * nLat][k] < dVectorBoundsLon[i][0]) {
					dVectorBoundsLon[i][0] = dVertexLon[i * nLat][k];
				}
				if (dVertexLon[i * nLat][k] > dVectorBoundsLon[i][1]) {
					dVectorBoundsLon[i][1] = dVertexLon[i * nLat][k];
				}
			}

			if ((dVectorBoundsLon[i][0] < 90.0) &&
			    (dVectorBoundsLon[i][1] > 270.0)
			) {
				if (i == 0) {
					double dTemp = dVectorBoundsLon[i][0];
					dVectorBoundsLon[i][0] = dVectorBoundsLon[i][1] - 360.0;
					dVectorBoundsLon[i][1] = dTemp;
				} else {
					double dTemp = dVectorBoundsLon[i][1];
					dVectorBoundsLon[i][1] = dVectorBoundsLon[i][0] + 360.0;
					dVectorBoundsLon[i][0] = dTemp;
				}
			}
			dVectorCenterLon[i] = 0.5 * (
				  dVectorBoundsLon[i][0]
				+ dVectorBoundsLon[i][1]);
		}
		for (int j = 0; j < nLat; j++) {
			//dVectorCenterLat[j] = dCenterLat[j];
			dVectorBoundsLat[j][0] = dVertexLat[j][0];
			dVectorBoundsLat[j][1] = dVertexLat[j][0];
			for (int k = 0; k < dVertexLat.GetColumns(); k++) {
				if (dVertexLat[j][k] < dVectorBoundsLat[j][0]) {
					dVectorBoundsLat[j][0] = dVertexLat[j][k];
				}
				if (dVertexLat[j][k] > dVectorBoundsLat[j][1]) {
					dVectorBoundsLat[j][1] = dVertexLat[j][k];
				}
			}
			dVectorCenterLat[j] = 0.5 * (
				  dVectorBoundsLat[j][0]
				+ dVectorBoundsLat[j][1]);
		}

		for (int i = 0; i < nFaces; i++) {
			dCenterLon[i] = dVectorCenterLon[i/nLat];
			dCenterLat[i] = dVectorCenterLat[i%nLat];
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void OfflineMap::InitializeSourceCoordinatesFromMeshFV(
	const Mesh & meshSource
) {
	// Check if these arrays have been read from file
	if ((m_dSourceVertexLon.IsInitialized()) ||
		(m_dSourceVertexLat.IsInitialized()) ||
		(m_dSourceCenterLon.IsInitialized()) ||
		(m_dSourceCenterLat.IsInitialized())
	) {
		if ((m_dSourceVertexLon.IsInitialized()) &&
			(m_dSourceVertexLat.IsInitialized()) &&
			(m_dSourceCenterLon.IsInitialized()) &&
			(m_dSourceCenterLat.IsInitialized())
		) {
			return;
		}

		_EXCEPTIONT("Logic error");
	}

	// Generate arrays
	bool fLatLon = false;
	if ((m_vecSourceDimNames[0] == "lat") &&
	    (m_vecSourceDimNames[1] == "lon")
	) {
		fLatLon = true;
	}
	if ((m_vecSourceDimNames[0] == "lon") &&
	    (m_vecSourceDimNames[1] == "lat")
	) {
		fLatLon = true;
	}

	InitializeCoordinatesFromMeshFV(
		meshSource,
		m_dSourceCenterLon,
		m_dSourceCenterLat,
		m_dSourceVertexLon,
		m_dSourceVertexLat,
		fLatLon);

	// Initialize coordinate arrays
	if (fLatLon) {
		if (m_vecSourceDimNames[0] == "lon") {
			InitializeRectilinearCoordinateVector(
				m_vecSourceDimSizes[0],
				m_vecSourceDimSizes[1],
				m_dSourceVertexLon,
				m_dSourceVertexLat,
				true,
				m_dSourceCenterLon,
				m_dSourceCenterLat,
				m_dVectorSourceCenterLon,
				m_dVectorSourceCenterLat,
				m_dVectorSourceBoundsLon,
				m_dVectorSourceBoundsLat
			);

		} else if (m_vecSourceDimNames[1] == "lon") {
			InitializeRectilinearCoordinateVector(
				m_vecSourceDimSizes[1],
				m_vecSourceDimSizes[0],
				m_dSourceVertexLon,
				m_dSourceVertexLat,
				false,
				m_dSourceCenterLon,
				m_dSourceCenterLat,
				m_dVectorSourceCenterLon,
				m_dVectorSourceCenterLat,
				m_dVectorSourceBoundsLon,
				m_dVectorSourceBoundsLat
			);

		} else {
			_EXCEPTIONT("LatLon specified but no dimensions have name \"lon\"");
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void OfflineMap::InitializeTargetCoordinatesFromMeshFV(
	const Mesh & meshTarget
) {
	// Check if these arrays have been read from file
	if ((m_dTargetVertexLon.IsInitialized()) ||
		(m_dTargetVertexLat.IsInitialized()) ||
		(m_dTargetCenterLon.IsInitialized()) ||
		(m_dTargetCenterLat.IsInitialized())
	) {
		if ((m_dTargetVertexLon.IsInitialized()) &&
			(m_dTargetVertexLat.IsInitialized()) &&
			(m_dTargetCenterLon.IsInitialized()) &&
			(m_dTargetCenterLat.IsInitialized())
		) {
			return;
		}

		_EXCEPTIONT("Logic error");
	}

	// Initialize coordinate arrays
	bool fLatLon = false;
	if ((m_vecTargetDimNames[0] == "lat") &&
	    (m_vecTargetDimNames[1] == "lon")
	) {
		fLatLon = true;
	}
	if ((m_vecTargetDimNames[0] == "lon") &&
	    (m_vecTargetDimNames[1] == "lat")
	) {
		fLatLon = true;
	}

	InitializeCoordinatesFromMeshFV(
		meshTarget,
		m_dTargetCenterLon,
		m_dTargetCenterLat,
		m_dTargetVertexLon,
		m_dTargetVertexLat,
		fLatLon);

	// Initialize vector coordinate
	if (fLatLon) {
		if (m_vecTargetDimNames[0] == "lon") {
			InitializeRectilinearCoordinateVector(
				m_vecTargetDimSizes[0],
				m_vecTargetDimSizes[1],
				m_dTargetVertexLon,
				m_dTargetVertexLat,
				true,
				m_dTargetCenterLon,
				m_dTargetCenterLat,
				m_dVectorTargetCenterLon,
				m_dVectorTargetCenterLat,
				m_dVectorTargetBoundsLon,
				m_dVectorTargetBoundsLat
			);

		} else if (m_vecTargetDimNames[1] == "lon") {
			InitializeRectilinearCoordinateVector(
				m_vecTargetDimSizes[1],
				m_vecTargetDimSizes[0],
				m_dTargetVertexLon,
				m_dTargetVertexLat,
				false,
				m_dTargetCenterLon,
				m_dTargetCenterLat,
				m_dVectorTargetCenterLon,
				m_dVectorTargetCenterLat,
				m_dVectorTargetBoundsLon,
				m_dVectorTargetBoundsLat
			);

		} else {
			_EXCEPTIONT("LatLon specified but no dimensions have name \"lon\"");
		}
/*
		// Bounds on vertex array
		m_dVectorTargetBoundsLon.Initialize(2);
		m_dVectorTargetBoundsLat.Initialize(2);
		m_dVectorTargetBoundsLon[0] = m_dTargetVertexLon[0][0];
		m_dVectorTargetBoundsLon[1] = m_dTargetVertexLon[0][0];
		m_dVectorTargetBoundsLat[0] = m_dTargetVertexLat[0][0];
		m_dVectorTargetBoundsLat[1] = m_dTargetVertexLat[0][0];

		if ((m_dTargetVertexLon.GetRows() != m_dTargetVertexLat.GetRows()) ||
		    (m_dTargetVertexLon.GetRows() != m_dTargetVertexLat.GetRows())
		) {
			_EXCEPTIONT("Catastrophic vertex array size mismatch");
		}

		for (int i = 0; i < m_dTargetVertexLon.GetRows(); i++) {
		for (int j = 0; j < m_dTargetVertexLon.GetColumns(); j++) {
			if (m_dTargetVertexLon[i][j] < m_dVectorTargetBoundsLon[0]) {
				m_dVectorTargetBoundsLon[0] = m_dTargetVertexLon[i][j];
			}
			if (m_dTargetVertexLon[i][j] > m_dVectorTargetBoundsLon[1]) {
				m_dVectorTargetBoundsLon[1] = m_dTargetVertexLon[i][j];
			}
			if (m_dTargetVertexLat[i][j] < m_dVectorTargetBoundsLat[0]) {
				m_dVectorTargetBoundsLat[0] = m_dTargetVertexLat[i][j];
			}
			if (m_dTargetVertexLat[i][j] > m_dVectorTargetBoundsLat[1]) {
				m_dVectorTargetBoundsLat[1] = m_dTargetVertexLat[i][j];
			}
		}
		}
*/

	}
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
	NcFile & ncFile,
	const std::string & strDimName,
	int nSize
) {
	NcDim * dim = ncFile.get_dim(strDimName.c_str());
	if (dim == NULL) {
		//std::cout << strDimName.c_str() << ", " << (long)nSize << std::endl;
		dim = ncFile.add_dim(strDimName.c_str(), (long)nSize);
		if (dim == NULL) {
			_EXCEPTION2("Failed to add dimension \"%s\" (%i) to file",
				strDimName.c_str(), nSize);
		}
	}
	
	if (dim->size() != nSize) {
		_EXCEPTION3("NetCDF file has dimension \"%s\" with mismatched"
			" size %i != %i", strDimName.c_str(), dim->size(), nSize);
	}
	return dim;
}

///////////////////////////////////////////////////////////////////////////////

void OfflineMap::PreserveVariables(
	const std::string & strSourceDataFile,
	const std::string & strTargetDataFile,
	const std::vector<std::string> & vecPreserveVariables
) {
	// Open source data file
	NcFile ncSource(strSourceDataFile.c_str(), NcFile::ReadOnly);
	if (!ncSource.is_valid()) {
		_EXCEPTION1("Cannot open source data file \"%s\" for reading",
			strSourceDataFile.c_str());
	}

	// Open target data file
	NcFile ncTarget(strTargetDataFile.c_str(), NcFile::Write);
	if (!ncTarget.is_valid()) {
		_EXCEPTION1("Cannot open target data file \"%s\" for writing",
			strTargetDataFile.c_str());
	}

	// Copy over dimensions
	for (int v = 0; v < vecPreserveVariables.size(); v++) {

		if (ncTarget.get_var(vecPreserveVariables[v].c_str()) != NULL) {
			Announce("%s (already exists, skipping)",
				vecPreserveVariables[v].c_str());

		} else {
			Announce("%s", vecPreserveVariables[v].c_str());

			CopyNcVar(ncSource, ncTarget, vecPreserveVariables[v]);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void OfflineMap::PreserveAllVariables(
	const std::string & strSourceDataFile,
	const std::string & strTargetDataFile
) {
	// Open source data file
	NcFile ncSource(strSourceDataFile.c_str(), NcFile::ReadOnly);
	if (!ncSource.is_valid()) {
		_EXCEPTION1("Cannot open source data file \"%s\"",
			strSourceDataFile.c_str());
	}

	// Check for rectilinear data
	bool fSourceRectilinear;
	if (m_vecSourceDimSizes.size() == 1) {
		fSourceRectilinear = false;
		if (m_vecSourceDimSizes.size() < 1) {
			_EXCEPTIONT("vecSourceDimSizes has not been initialized");
		}

	} else if (m_vecSourceDimSizes.size() == 2) {
		fSourceRectilinear = true;
		if (m_vecSourceDimSizes.size() < 2) {
			_EXCEPTIONT("vecSourceDimSizes has not been initialized");
		}

	} else {
		_EXCEPTIONT("m_vecSourceDimSizes undefined");
	}

	// Generate variable list
	std::vector<std::string> vecPreserveVariables;

	for (int v = 0; v < ncSource.num_vars(); v++) {
		NcVar * var = ncSource.get_var(v);
		if (var == NULL) {
			_EXCEPTION1("Error reading variable %i in source file", v);
		}

		if (fSourceRectilinear) {
			if (var->num_dims() >= 2) {
				NcDim * dimA = var->get_dim(var->num_dims()-2);
				NcDim * dimB = var->get_dim(var->num_dims()-1);

				if (dimA->size() == m_vecSourceDimSizes[0]) {
					continue;
				}
				if (dimB->size() == m_vecSourceDimSizes[1]) {
					continue;
				}
				if (strcmp(dimA->name(), m_vecSourceDimNames[0].c_str()) == 0) {
					continue;
				}
				if (strcmp(dimB->name(), m_vecSourceDimNames[1].c_str()) == 0) {
					continue;
				}
			}

		} else {
			int nSourceCount = m_dSourceAreas.GetRows();

			if (var->num_dims() >= 1) {
				NcDim * dim = var->get_dim(var->num_dims()-1);

				if (dim->size() == nSourceCount) {
					continue;
				}
				if (strcmp(dim->name(), m_vecSourceDimNames[0].c_str()) == 0) {
					continue;
				}
			}
		}

		vecPreserveVariables.push_back(var->name());
	}

	// Preserve all variables in the list
	PreserveVariables(
		strSourceDataFile,
		strTargetDataFile,
		vecPreserveVariables);
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

	// Check variable list for "lat" and "lon"
	for (int v = 0; v < vecVariables.size(); v++) {
		if (vecVariables[v] == "lat") {
			_EXCEPTIONT("Latitude variable \"lat\" in variable list will be overwritten on output");
		}
		if (vecVariables[v] == "lon") {
			_EXCEPTIONT("Longitude variable \"lon\" in variable list will be overwritten on output");
		}
	}

	// Open source data file
	NcFile ncSource(strSourceDataFile.c_str(), NcFile::ReadOnly);
	if (!ncSource.is_valid()) {
		_EXCEPTION1("Cannot open source data file \"%s\"",
			strSourceDataFile.c_str());
	}

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

	// Generate variable list
	std::vector<std::string> vecVariableList = vecVariables;

	if (vecVariables.size() == 0) {
		for (int v = 0; v < ncSource.num_vars(); v++) {
			NcVar * var = ncSource.get_var(v);

			if (fSourceRectilinear) {
				if (var->num_dims() < 2) {
					continue;

				} else {
					NcDim * dimA = var->get_dim(var->num_dims()-2);
					NcDim * dimB = var->get_dim(var->num_dims()-1);

					if (dimA->size() != m_vecSourceDimSizes[0]) {
						continue;
					}
					if (dimB->size() != m_vecSourceDimSizes[1]) {
						continue;
					}
				}

			} else {
				if (var->num_dims() < 1) {
					continue;

				} else {
					NcDim * dim = var->get_dim(var->num_dims()-1);

					if (dim->size() != nSourceCount) {
						continue;
					}

					bool fDimensionName = false;
					for (int d = 0; d < m_vecTargetDimNames.size(); d++) {
						const char * szDimName = m_vecTargetDimNames[d].c_str();
						if (strcmp(var->name(), szDimName) == 0) {
							fDimensionName = true;
							break;
						}
					}
					if (fDimensionName) {
						continue;
					}
				}
			}

			vecVariableList.push_back(var->name());
		}
	}

	// Add lat/lon vector dimensions
	if ((m_dVectorTargetCenterLon.GetRows() != 0) &&
		(m_dVectorTargetCenterLat.GetRows() != 0)
	) {
		NcDim * dimLon;
		if (m_vecTargetDimNames[0] == "lon") {
			dimLon = dim0;
		}
		if (m_vecTargetDimNames[1] == "lon") {
			dimLon = dim1;
		}
		NcVar * varLon = ncTarget.get_var("lon");
		if (varLon == NULL) {
			varLon = ncTarget.add_var("lon", ncDouble, dimLon);
			if (varLon == NULL) {
				_EXCEPTIONT("Cannot create variable \"lon\" in target file");
			}

			varLon->put(
				&(m_dVectorTargetCenterLon[0]),
				m_dVectorTargetCenterLon.GetRows());

			varLon->add_att("bounds", "lon_bnds");
			varLon->add_att("units", "degrees_east");
			varLon->add_att("axis", "X");
			varLon->add_att("long_name", "longitude");
			varLon->add_att("standard_name", "longitude");

		} else {
			if (varLon->get_dim(0)->size() != dimLon->size()) {
				_EXCEPTIONT("\"lon\" variable mismatch");
			}
		}

		NcDim * dimLat;
		if (m_vecTargetDimNames[0] == "lat") {
			dimLat = dim0;
		}
		if (m_vecTargetDimNames[1] == "lat") {
			dimLat = dim1;
		}
		NcVar * varLat = ncTarget.get_var("lat");
		if (varLat == NULL) {
			varLat = ncTarget.add_var("lat", ncDouble, dimLat);
			if (varLat == NULL) {
				_EXCEPTIONT("Cannot create variable \"lat\" in target file");
			}

			varLat->put(
				&(m_dVectorTargetCenterLat[0]),
				m_dVectorTargetCenterLat.GetRows());

			varLat->add_att("bounds", "lat_bnds");
			varLat->add_att("units", "degrees_north");
			varLat->add_att("axis", "Y");
			varLat->add_att("long_name", "latitude");
			varLat->add_att("standard_name", "latitude");

		} else {
			if (varLat->get_dim(0)->size() != dimLat->size()) {
				_EXCEPTIONT("\"lat\" variable mismatch");
			}
		}

		// Output bounds variables
		if ((m_dVectorTargetBoundsLon.GetRows() != 0) &&
		    (m_dVectorTargetBoundsLat.GetRows() != 0) &&
			(ncTarget.get_dim("bnds") != NULL)
		) {
			NcDim * dimBounds = ncTarget.add_dim("bnds", 2);

			NcVar * varLonBounds =
				ncTarget.add_var("lon_bnds", ncDouble, dimLon, dimBounds);

			if (varLonBounds == NULL) {
				_EXCEPTIONT("Cannot create variable \"lon_bnds\""
					" in target file");
			}

			varLonBounds->put(&(m_dVectorTargetBoundsLon[0][0]),
				m_dVectorTargetBoundsLon.GetRows(), 2);

			NcVar * varLatBounds =
				ncTarget.add_var("lat_bnds", ncDouble, dimLat, dimBounds);

			if (varLatBounds == NULL) {
				_EXCEPTIONT("Cannot create variable \"lat_bnds\""
					" in target file");
			}

			varLatBounds->put(&(m_dVectorTargetBoundsLat[0][0]),
				m_dVectorTargetBoundsLat.GetRows(), 2);
		}

	// Add lat/lon to unstructured data
	} else if (!fTargetRectilinear) {
		NcVar * varLon = ncTarget.get_var("lon");
		if (varLon == NULL) {
			varLon = ncTarget.add_var("lon", ncDouble, dim0);
			if (m_dTargetCenterLon.GetRows() != dim0->size()) {
				_EXCEPTION2("TargetCenterLon / NCol dimension size mismatch (%i, %i)",
					m_dTargetCenterLon.GetRows(), dim0->size());
			}
			if (varLon == NULL) {
				_EXCEPTIONT("Cannot create variable \"lon\" in target file");
			}

			varLon->put(
				&(m_dTargetCenterLon[0]),
				m_dTargetCenterLon.GetRows());

			//varLon->add_att("bounds", "lon_bnds");
			varLon->add_att("units", "degrees_east");
			varLon->add_att("axis", "X");
			varLon->add_att("long_name", "longitude");
			varLon->add_att("standard_name", "longitude");

		} else {
			if (varLon->get_dim(0)->size() != dim0->size()) {
				_EXCEPTIONT("\"lon\" variable mismatch");
			}
		}

		NcVar * varLat = ncTarget.get_var("lat");
		if (varLat == NULL) {
			varLat = ncTarget.add_var("lat", ncDouble, dim0);
			if (m_dTargetCenterLat.GetRows() != dim0->size()) {
				_EXCEPTION2("TargetCenterLat / NCol dimension size mismatch (%i, %i)",
					m_dTargetCenterLat.GetRows(), dim0->size());
			}
			if (varLat == NULL) {
				_EXCEPTIONT("Cannot create variable \"lat\" in target file");
			}

			varLat->put(
				&(m_dTargetCenterLat[0]),
				m_dTargetCenterLat.GetRows());

			//varLat->add_att("bounds", "lat_bnds");
			varLat->add_att("units", "degrees_north");
			varLat->add_att("axis", "Y");
			varLat->add_att("long_name", "latitude");
			varLat->add_att("standard_name", "latitude");

		} else {
			if (varLat->get_dim(0)->size() != dim0->size()) {
				_EXCEPTIONT("\"lat\" variable mismatch");
			}
		}

	// Non-latitude-longitude rectilinear mesh target
	} else {
		NcVar * varLon = ncTarget.get_var("lon");
		if (varLon == NULL) {
			varLon = ncTarget.add_var("lon", ncDouble, dim0, dim1);
			if (m_dTargetCenterLon.GetRows() != dim0->size() * dim1->size()) {
				_EXCEPTION3("TargetCenterLon / NCol dimension size mismatch (%i, %i x %i)",
					m_dTargetCenterLon.GetRows(), dim0->size(), dim1->size());
			}
			if (varLon == NULL) {
				_EXCEPTIONT("Cannot create variable \"lon\" in target file");
			}

			varLon->put(
				&(m_dTargetCenterLon[0]),
				dim0->size(),
				dim1->size());

			//varLon->add_att("bounds", "lon_bnds");
			varLon->add_att("units", "degrees_east");
			varLon->add_att("axis", "X");
			varLon->add_att("long_name", "longitude");
			varLon->add_att("standard_name", "longitude");

		} else {
			if (varLon->num_dims() != 2) {
				_EXCEPTIONT("\"lon\" variable mismatch");
			}
			if (varLon->get_dim(0)->size() != dim0->size()) {
				_EXCEPTIONT("\"lon\" variable mismatch");
			}
			if (varLon->get_dim(1)->size() != dim1->size()) {
				_EXCEPTIONT("\"lon\" variable mismatch");
			}
		}

		NcVar * varLat = ncTarget.get_var("lat");
		if (varLat == NULL) {
			varLat = ncTarget.add_var("lat", ncDouble, dim0, dim1);
			if (m_dTargetCenterLat.GetRows() != dim0->size() * dim1->size()) {
				_EXCEPTION3("TargetCenterLat / NCol dimension size mismatch (%i, %i x %i)",
					m_dTargetCenterLat.GetRows(), dim0->size(), dim1->size());
			}
			if (varLat == NULL) {
				_EXCEPTIONT("Cannot create variable \"lat\" in target file");
			}

			varLat->put(
				&(m_dTargetCenterLat[0]),
				dim0->size(),
				dim1->size());

			//varLat->add_att("bounds", "lat_bnds");
			varLat->add_att("units", "degrees_north");
			varLat->add_att("axis", "Y");
			varLat->add_att("long_name", "latitude");
			varLat->add_att("standard_name", "latitude");

		} else {
			if (varLat->num_dims() != 2) {
				_EXCEPTIONT("\"lat\" variable mismatch");
			}
			if (varLat->get_dim(0)->size() != dim0->size()) {
				_EXCEPTIONT("\"lat\" variable mismatch");
			}
			if (varLat->get_dim(1)->size() != dim1->size()) {
				_EXCEPTIONT("\"lat\" variable mismatch");
			}
		}
	}

	// Loop through all variables
	for (int v = 0; v < vecVariableList.size(); v++) {
		NcVar * var = ncSource.get_var(vecVariableList[v].c_str());
		if (var == NULL) {
			_EXCEPTION1("Variable \"%s\" does not exist in source file",
				vecVariableList[v].c_str());
		}

		AnnounceStartBlock(vecVariableList[v].c_str());

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

		// Construct an array of dimensions for this variable
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
				if (var->num_dims() == 1) {
					_EXCEPTIONT("Rectilinear source data stored in 1D format");
				}

				vecDimsOut.Initialize(var->num_dims()-1);
			} else {
				vecDimsOut.Initialize(var->num_dims());
			}

			vecDimsOut[vecDimsOut.GetRows()-1] = dim0;
		}

		int nFreeDims = var->num_dims() - m_vecSourceDimSizes.size();

		// Add any missing dimension variables to target file
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

			// Copy over associated variable, if it exists
			NcVar * varAssocDimSource =
				ncSource.get_var(strDimName.c_str());
			if (varAssocDimSource != NULL) {
				NcVar * varAssocDimTarget =
					ncTarget.get_var(strDimName.c_str());
				if (varAssocDimTarget == NULL) {
					CopyNcVar(ncSource, ncTarget, strDimName);
				}

				NcAtt * attAssocDimBounds =
					varAssocDimSource->get_att("bounds");
				if (attAssocDimBounds != NULL) {
					NcVar * varAssocDimBoundsSource =
						ncSource.get_var(attAssocDimBounds->as_string(0));
					NcVar * varAssocDimBoundsTarget =
						ncTarget.get_var(attAssocDimBounds->as_string(0));
					if ((varAssocDimBoundsSource != NULL) &&
					    (varAssocDimBoundsTarget == NULL)
					) {
						CopyNcVar(
							ncSource, ncTarget,
							attAssocDimBounds->as_string(0));
					}
				}
			}
		}

		// Create new output variable
		NcVar * varOut;
		if (fTargetDouble) {
			varOut =
				ncTarget.add_var(
					vecVariableList[v].c_str(),
					ncDouble,
					vecDimsOut.GetRows(),
					(const NcDim**)&(vecDimsOut[0]));

		} else {
			varOut =
				ncTarget.add_var(
					vecVariableList[v].c_str(),
					ncFloat,
					vecDimsOut.GetRows(),
					(const NcDim**)&(vecDimsOut[0]));
		}

		if (varOut == NULL) {
			_EXCEPTION1("Cannot create variable \"%s\" in output file",
				vecVariableList[v].c_str());
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
			for (int d = vecDimSizes.GetRows()-1; d >= 0; d--) {
				nCountsIn[d]  = tt % vecDimSizes[d];
				nCountsOut[d] = tt % vecDimSizes[d];
				tt /= vecDimSizes[d];
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
/*
			for (int d = 0; d < nCountsIn.GetRows(); d++) {
				printf("%li ", nCountsIn[d]);
			}
			printf(" : ");
			for (int d = 0; d < nCountsOut.GetRows(); d++) {
				printf("%li ", nCountsOut[d]);
			}
			printf("\n");
*/
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

	for (int i = 0; i < nSrcGridDims/2; i++) {
		int iTemp = m_vecSourceDimSizes[i];
		m_vecSourceDimSizes[i] = m_vecSourceDimSizes[nSrcGridDims - i - 1];
		m_vecSourceDimSizes[nSrcGridDims - i - 1] = iTemp;
	}

	for (int i = 0; i < nSrcGridDims; i++) {
		char szDim[64];
		sprintf(szDim, "name%i", nSrcGridDims - i - 1);
		NcAtt * attDim = varSrcGridDims->get_att(szDim);
		if (attDim != NULL) {
			m_vecSourceDimNames[i] = attDim->as_string(0);
		} else {
			sprintf(szDim, "dim%i", nSrcGridDims - i - 1);
			m_vecSourceDimNames[i] = szDim;
		}
	}

	for (int i = 0; i < nDstGridDims/2; i++) {
		int iTemp = m_vecTargetDimSizes[i];
		m_vecTargetDimSizes[i] = m_vecTargetDimSizes[nDstGridDims - i - 1];
		m_vecTargetDimSizes[nDstGridDims - i - 1] = iTemp;
	}

	for (int i = 0; i < nDstGridDims; i++) {
		char szDim[64];
		sprintf(szDim, "name%i", nDstGridDims - i - 1);
		NcAtt * attDim = varDstGridDims->get_att(szDim);
		if (attDim != NULL) {
			m_vecTargetDimNames[i] = attDim->as_string(0);
		} else {
			sprintf(szDim, "dim%i", nDstGridDims - i - 1);
			m_vecTargetDimNames[i] = szDim;
		}
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

	// Read vector centers and bounds
	NcDim * dimLatB = ncMap.get_dim("lat_b");
	NcDim * dimLonB = ncMap.get_dim("lon_b");

	if ((dimLatB != NULL) && (dimLonB != NULL)) {
		NcVar * varLatCB = ncMap.get_var("latc_b");
		NcVar * varLonCB = ncMap.get_var("lonc_b");

		m_dVectorTargetCenterLat.Initialize(dimLatB->size());
		m_dVectorTargetCenterLon.Initialize(dimLonB->size());

		varLatCB->get(&(m_dVectorTargetCenterLat[0]), dimLatB->size());
		varLonCB->get(&(m_dVectorTargetCenterLon[0]), dimLonB->size());

		NcVar * varLatBounds = ncMap.get_var("lat_bnds");
		NcVar * varLonBounds = ncMap.get_var("lon_bnds");

		m_dVectorTargetBoundsLat.Initialize(dimLatB->size(), 2);
		m_dVectorTargetBoundsLon.Initialize(dimLonB->size(), 2);

		varLatBounds->get(&(m_dVectorTargetBoundsLat[0][0]),
			dimLatB->size(), 2);
		varLonBounds->get(&(m_dVectorTargetBoundsLon[0][0]),
			dimLonB->size(), 2);
	}

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
	const std::string & strTarget,
	const std::map<std::string, std::string> & mapAttributes
) {
	NcFile ncMap(strTarget.c_str(), NcFile::Replace);
	if (!ncMap.is_valid()) {
		_EXCEPTION1("Unable to open output map file \"%s\"",
			strTarget.c_str());
	}

	// Attributes
	ncMap.add_att("Title", "TempestRemap Offline Regridding Weight Generator");

	// Map dimensions
	int nA = (int)(m_dSourceAreas.GetRows());
	int nB = (int)(m_dTargetAreas.GetRows());

	// Write output dimensions entries
	int nSrcGridDims = (int)(m_vecSourceDimSizes.size());
	int nDstGridDims = (int)(m_vecTargetDimSizes.size());

	NcDim * dimSrcGridRank = ncMap.add_dim("src_grid_rank", nSrcGridDims);
	NcDim * dimDstGridRank = ncMap.add_dim("dst_grid_rank", nDstGridDims);

	NcVar * varSrcGridDims =
		ncMap.add_var("src_grid_dims", ncInt, dimSrcGridRank);
	NcVar * varDstGridDims =
		ncMap.add_var("dst_grid_dims", ncInt, dimDstGridRank);

	char szDim[64];
	if ((nSrcGridDims == 1) && (m_vecSourceDimSizes[0] != nA)) {
		varSrcGridDims->put(&nA, 1);
		varSrcGridDims->add_att("name0", "num_dof");

	} else {
		for (int i = 0; i < m_vecSourceDimSizes.size(); i++) {
			varSrcGridDims->set_cur(nSrcGridDims - i - 1);
			varSrcGridDims->put(&(m_vecSourceDimSizes[i]), 1);
		}

		for (int i = 0; i < m_vecSourceDimSizes.size(); i++) {
			sprintf(szDim, "name%i", i);
			varSrcGridDims->add_att(szDim,
				m_vecSourceDimNames[nSrcGridDims - i - 1].c_str());
		}
	}

	if ((nDstGridDims == 1) && (m_vecTargetDimSizes[0] != nB)) {
		varDstGridDims->put(&nB, 1);
		varDstGridDims->add_att("name0", "num_dof");

	} else {
		for (int i = 0; i < m_vecTargetDimSizes.size(); i++) {
			varDstGridDims->set_cur(nDstGridDims - i - 1);
			varDstGridDims->put(&(m_vecTargetDimSizes[i]), 1);
		}

		for (int i = 0; i < m_vecTargetDimSizes.size(); i++) {
			sprintf(szDim, "name%i", i);
			varDstGridDims->add_att(szDim,
				m_vecTargetDimNames[nDstGridDims - i - 1].c_str());
		}
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

	// Write vector centers
	if ((m_dVectorTargetCenterLat.GetRows() != 0) &&
		(m_dVectorTargetCenterLon.GetRows() != 0)
	) {
		NcDim * dimLatB =
			ncMap.add_dim("lat_b", m_dVectorTargetCenterLat.GetRows());
		NcDim * dimLonB =
			ncMap.add_dim("lon_b", m_dVectorTargetCenterLon.GetRows());

		NcVar * varLatCB = ncMap.add_var("latc_b", ncDouble, dimLatB);
		NcVar * varLonCB = ncMap.add_var("lonc_b", ncDouble, dimLonB);

		varLatCB->put(&(m_dVectorTargetCenterLat[0]), dimLatB->size());
		varLonCB->put(&(m_dVectorTargetCenterLon[0]), dimLonB->size());

		NcDim * dimBounds = ncMap.add_dim("bnds", 2);
		NcVar * varLatBounds =
			ncMap.add_var("lat_bnds", ncDouble, dimLatB, dimBounds);
		NcVar * varLonBounds =
			ncMap.add_var("lon_bnds", ncDouble, dimLonB, dimBounds);

		varLatBounds->put(&(m_dVectorTargetBoundsLat[0][0]),
			m_dVectorTargetBoundsLat.GetRows(), 2);
		varLonBounds->put(&(m_dVectorTargetBoundsLon[0][0]),
			m_dVectorTargetBoundsLon.GetRows(), 2);
	}

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

	// Add global attributes
	std::map<std::string, std::string>::const_iterator iterAttributes =
		mapAttributes.begin();
	for (; iterAttributes != mapAttributes.end(); iterAttributes++) {
		ncMap.add_att(
			iterAttributes->first.c_str(),
			iterAttributes->second.c_str());
	}
}

///////////////////////////////////////////////////////////////////////////////

void OfflineMap::Write(
	const std::string & strTarget
) {
	std::map<std::string, std::string> mapNoAttributes;
	return Write(strTarget, mapNoAttributes);
}

///////////////////////////////////////////////////////////////////////////////

void OfflineMap::SetTranspose(
	const OfflineMap & mapIn
) {
	m_dSourceAreas = mapIn.m_dTargetAreas;
	m_dTargetAreas = mapIn.m_dSourceAreas;

	m_dSourceCenterLon = mapIn.m_dTargetCenterLon;
	m_dSourceCenterLat = mapIn.m_dTargetCenterLat;
	m_dTargetCenterLon = mapIn.m_dSourceCenterLon;
	m_dTargetCenterLat = mapIn.m_dSourceCenterLat;

	m_dSourceVertexLon = mapIn.m_dTargetVertexLon;
	m_dSourceVertexLat = mapIn.m_dTargetVertexLat;
	m_dVectorSourceCenterLon = mapIn.m_dVectorTargetCenterLon;
	m_dVectorSourceCenterLat = mapIn.m_dVectorTargetCenterLat;
	m_dVectorSourceBoundsLon = mapIn.m_dVectorTargetBoundsLon;
	m_dVectorSourceBoundsLat = mapIn.m_dVectorTargetBoundsLat;

	m_dTargetVertexLon = mapIn.m_dSourceVertexLon;
	m_dTargetVertexLat = mapIn.m_dSourceVertexLat;
	m_dVectorTargetCenterLon = mapIn.m_dVectorSourceCenterLon;
	m_dVectorTargetCenterLat = mapIn.m_dVectorSourceCenterLat;
	m_dVectorTargetBoundsLon = mapIn.m_dVectorSourceBoundsLon;
	m_dVectorTargetBoundsLat = mapIn.m_dVectorSourceBoundsLat;

	m_vecSourceDimSizes = mapIn.m_vecTargetDimSizes;
	m_vecSourceDimNames = mapIn.m_vecTargetDimNames;

	m_vecTargetDimSizes = mapIn.m_vecSourceDimSizes;
	m_vecTargetDimNames = mapIn.m_vecSourceDimNames;

	DataVector<int> dataRows;
	DataVector<int> dataCols;
	DataVector<double> dataEntries;

	mapIn.m_mapRemap.GetEntries(dataRows, dataCols, dataEntries);
/*
	DataVector<double> dTargetCoverage(mapIn.m_dTargetAreas.GetRows());
	for (size_t s = 0; s < dataEntries.GetRows(); s++) {
		dTargetCoverage[dataRows[s]] +=
			mapIn.m_dTargetAreas[dataRows[s]];
	}
*/
	for (size_t s = 0; s < dataEntries.GetRows(); s++) {
		dataEntries[s] *=
			mapIn.m_dTargetAreas[dataRows[s]]
			/ mapIn.m_dSourceAreas[dataCols[s]];
	}

	m_mapRemap.SetEntries(dataCols, dataRows, dataEntries);
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

