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
#include "STLStringHelper.h"

#include "Announce.h"
#include "Exception.h"
#include "DataArray1D.h"
#include "DataArray2D.h"

#include <cmath>

///////////////////////////////////////////////////////////////////////////////

void ParseEnforceBounds(
	const std::string & strEnforceBounds,
	EnforceBoundsVector & vecEnforceBounds
) {
	enum ParseMode {
		ParseMode_Variable,
		ParseMode_LowerBound,
		ParseMode_UpperBound
	} eParseMode = ParseMode_Variable;

	if (strEnforceBounds.length() == 0) {
		return;
	}

	int iCommaCount = 0;
	int iLastComma = 0;
	for (int i = 0; i < strEnforceBounds.length(); i++) {
		if (strEnforceBounds[i] == ',') {
			iCommaCount++;
			iLastComma = i;
		}
	}

	EnforceBounds enfbnds;

	// Parse string of the form "<lb>,<ub>"
	if (iCommaCount == 1) {
		enfbnds.strVariable = "";
		enfbnds.strLowerBound = strEnforceBounds.substr(0,iLastComma);
		enfbnds.strUpperBound = strEnforceBounds.substr(iLastComma+1);
		vecEnforceBounds.push_back(enfbnds);

	// Parse string of the form "<var>,<lb>,<ub>[;...]"
	} else {
		int iLast = 0;
		for (int i = 0; i <= strEnforceBounds.length(); i++) {
			if ((i == strEnforceBounds.length()) || (strEnforceBounds[i] == ';')) {
				if (eParseMode != ParseMode_UpperBound) {
					_EXCEPTION1("Malformed bounds string \"%s\", expected \"<var>,<lb>,<ub>[;...]\"",
						strEnforceBounds.c_str());
				}
				enfbnds.strUpperBound = strEnforceBounds.substr(iLast, i-iLast);
				vecEnforceBounds.push_back(enfbnds);
				if (i == strEnforceBounds.length()) {
					break;
				}
				eParseMode = ParseMode_Variable;
				iLast = i+1;
			}
			if (strEnforceBounds[i] == ',') {
				if (eParseMode == ParseMode_Variable) {
					enfbnds.strVariable = strEnforceBounds.substr(iLast, i-iLast);
					eParseMode = ParseMode_LowerBound;
					iLast = i+1;

				} else if (eParseMode == ParseMode_LowerBound) {
					enfbnds.strLowerBound = strEnforceBounds.substr(iLast, i-iLast);
					eParseMode = ParseMode_UpperBound;
					iLast = i+1;

				} else {
					_EXCEPTION1("Malformed bounds string \"%s\", expected \"<var>,<lb>,<ub>[;...]\"",
						strEnforceBounds.c_str());
				}
			}
		}
	}

	// Validate
	for (auto it : vecEnforceBounds) {
		if ((it.strVariable == "") && (vecEnforceBounds.size() != 1)) {
			_EXCEPTIONT("No variable specified in bounds string \"\"");
		}
		if (it.strLowerBound == "") {
			it.strLowerBound = "n";
		}
		if (it.strUpperBound == "") {
			it.strUpperBound = "n";
		}

		bool fLowerBoundIsFloat = STLStringHelper::IsFloat(it.strLowerBound);
		bool fUpperBoundIsFloat = STLStringHelper::IsFloat(it.strUpperBound);

		if ((it.strLowerBound != "n") && (it.strLowerBound != "l") && (it.strLowerBound != "g") && (!fLowerBoundIsFloat)) {
			_EXCEPTION1("Invalid lower bound in bounds string \"%s\", expected \"n\", \"l\", \"g\", or floating point value",
				it.strLowerBound.c_str());
		}
		if ((it.strUpperBound != "n") && (it.strUpperBound != "l") && (it.strUpperBound != "g") && (!fUpperBoundIsFloat)) {
			_EXCEPTION1("Invalid upper bound in bounds string \"%s\", expected \"n\", \"l\", \"g\", or floating point value",
				it.strUpperBound.c_str());
		}
		if (fLowerBoundIsFloat && fUpperBoundIsFloat) {
			double dLowerBound = std::stof(it.strLowerBound);
			double dUpperBound = std::stof(it.strUpperBound);

			if (dLowerBound > dUpperBound) {
				_EXCEPTION2("Lower bound \"%s\" must be less than upper bound \"%s\" in bounds string",
					it.strLowerBound.c_str(), it.strUpperBound.c_str());
			}
		}
	}
	//for (auto it : vecEnforceBounds) {
	//	printf("%s %s %s\n", it.strVariable.c_str(), it.strLowerBound.c_str(), it.strUpperBound.c_str());
	//}

}

///////////////////////////////////////////////////////////////////////////////

void OfflineMap::InitializeDimensionsFromMeshFile(
	const std::string & strMeshFile,
	std::vector<std::string> & vecDimNames,
	std::vector<int> & vecDimSizes,
	DataArray1D<double> & dCenterLon,
	DataArray1D<double> & dCenterLat,
	DataArray2D<double> & dVertexLon,
	DataArray2D<double> & dVertexLat
) {
	// Open the mesh
	NcFile ncMesh(strMeshFile.c_str(), NcFile::ReadOnly);
	if (!ncMesh.is_valid()) {
		_EXCEPTION1("Unable to open mesh file \"%s\"", strMeshFile.c_str());
	}

	// Check for grid dimensions (SCRIP format grid)
	NcVar * varGridDims = ncMesh.get_var("grid_dims");
	if (varGridDims != NULL) {
		NcDim * dimGridRank = varGridDims->get_dim(0);
		if (dimGridRank == NULL) {
			_EXCEPTIONT("Variable \"grid_dims\" has no dimensions");
		}

		vecDimSizes.resize(dimGridRank->size());
		varGridDims->get(&(vecDimSizes[0]), dimGridRank->size());

		if (dimGridRank->size() == 1) {
			vecDimNames.push_back("num_elem");
		} else if (dimGridRank->size() == 2) {
			vecDimNames.push_back("lat");
			vecDimNames.push_back("lon");

			int iTemp = vecDimSizes[0];
			vecDimSizes[0] = vecDimSizes[1];
			vecDimSizes[1] = iTemp;

		} else {
			_EXCEPTION1("Mesh file \"%s\" grid_rank must be < 3", strMeshFile.c_str());
		}

		// Number of faces
		NcDim * dimGridSize = ncMesh.get_dim("grid_size");
		if (dimGridSize == NULL) {
			_EXCEPTION1("Missing \"grid_size\" dimension in mesh file \"%s\"", strMeshFile.c_str());
		}

		// Number of grid corners
		NcDim * dimGridCorners = ncMesh.get_dim("grid_corners");
		if (dimGridCorners == NULL) {
			_EXCEPTION1("Missing \"grid_corners\" dimension in mesh file \"%s\"", strMeshFile.c_str());
		}

		// Pull grid center longitude information from file
		NcVar * varGridCenterLon = ncMesh.get_var("grid_center_lon");
		if (varGridCenterLon == NULL) {
			_EXCEPTION1("Missing \"grid_center_lon\" variable in mesh file \"%s\"", strMeshFile.c_str());
		}

		dCenterLon.Allocate(dimGridSize->size());

		varGridCenterLon->get(
			&(dCenterLon[0]),
			dimGridSize->size());

		// Convert radians to degrees
		NcAtt * attGridCenterLonUnits = varGridCenterLon->get_att("units");
		if (attGridCenterLonUnits != NULL) {
			std::string strGridCenterLonUnits = attGridCenterLonUnits->as_string(0);
			if (strGridCenterLonUnits == "degrees") {

			} else if (strGridCenterLonUnits == "radians") {
				for (int i = 0; i < dCenterLon.GetRows(); i++) {
					dCenterLon[i] *= 180.0 / M_PI;
				}

			} else {
				_EXCEPTION1("Invalid \"units\" attribute for \"grid_center_lon\" variable: Expected \"degrees\" or \"radians\" in mesh file \"%s\"", strMeshFile.c_str());
			}
		}

		// Pull grid center latitude information from file
		NcVar * varGridCenterLat = ncMesh.get_var("grid_center_lat");
		if (varGridCenterLat == NULL) {
			_EXCEPTION1("Missing \"grid_center_lat\" variable in mesh file \"%s\"", strMeshFile.c_str());
		}

		dCenterLat.Allocate(dimGridSize->size());

		varGridCenterLat->get(
			&(dCenterLat[0]),
			dimGridSize->size());

		// Convert radians to degrees
		NcAtt * attGridCenterLatUnits = varGridCenterLat->get_att("units");
		if (attGridCenterLatUnits != NULL) {
			std::string strGridCenterLatUnits = attGridCenterLatUnits->as_string(0);
			if (strGridCenterLatUnits == "degrees") {

			} else if (strGridCenterLatUnits == "radians") {
				for (int i = 0; i < dCenterLat.GetRows(); i++) {
					dCenterLat[i] *= 180.0 / M_PI;
				}

			} else {
				_EXCEPTION1("Invalid \"units\" attribute for \"grid_center_lat\" variable: Expected \"degrees\" or \"radians\" in mesh file \"%s\"", strMeshFile.c_str());
			}
		}

		// Pull longitude grid vertex information from file
		NcVar * varGridVertexLon = ncMesh.get_var("grid_corner_lon");
		if (varGridVertexLon == NULL) {
			_EXCEPTION1("Missing \"grid_corner_lon\" variable in mesh file \"%s\"", strMeshFile.c_str());
		}

		dVertexLon.Allocate(
			dimGridSize->size(),
			dimGridCorners->size());

		varGridVertexLon->get(
			&(dVertexLon[0][0]),
			dimGridSize->size(),
			dimGridCorners->size());

		// Convert radians to degrees
		NcAtt * attGridVertexLonUnits = varGridVertexLon->get_att("units");
		if (attGridVertexLonUnits != NULL) {
			std::string strGridVertexLonUnits = attGridVertexLonUnits->as_string(0);
			if (strGridVertexLonUnits == "degrees") {

			} else if (strGridVertexLonUnits == "radians") {
				for (int i = 0; i < dVertexLon.GetRows(); i++) {
				for (int j = 0; j < dVertexLon.GetColumns(); j++) {
					dVertexLon[i][j] *= 180.0 / M_PI;
				}
				}

			} else {
				_EXCEPTION1("Invalid \"units\" attribute for \"grid_corner_lon\" variable: Expected \"degrees\" or \"radians\" in mesh file \"%s\"", strMeshFile.c_str());
			}
		}

		// Pull latitude grid vertex information from file
		NcVar * varGridVertexLat = ncMesh.get_var("grid_corner_lat");
		if (varGridVertexLat == NULL) {
			_EXCEPTION1("Missing \"grid_corner_lat\" variable in mesh file \"%s\"", strMeshFile.c_str());
		}

		dVertexLat.Allocate(
			dimGridSize->size(),
			dimGridCorners->size());

		varGridVertexLat->get(
			&(dVertexLat[0][0]),
			dimGridSize->size(),
			dimGridCorners->size());

		// Convert radians to degrees
		NcAtt * attGridVertexLatUnits = varGridVertexLat->get_att("units");
		if (attGridVertexLatUnits != NULL) {
			std::string strGridVertexLatUnits = attGridVertexLatUnits->as_string(0);
			if (strGridVertexLatUnits == "degrees") {

			} else if (strGridVertexLatUnits == "radians") {
				for (int i = 0; i < dVertexLat.GetRows(); i++) {
				for (int j = 0; j < dVertexLat.GetColumns(); j++) {
					dVertexLat[i][j] *= 180.0 / M_PI;
				}
				}

			} else {
				_EXCEPTION1("Invalid \"units\" attribute for \"grid_corner_lat\" variable: Expected \"degrees\" or \"radians\" in mesh file \"%s\"", strMeshFile.c_str());
			}
		}

		return;
	}

	// Check for rectilinear attribute
	NcAtt * attRectilinear = ncMesh.get_att("rectilinear");

	// No rectilinear attribute
	if (attRectilinear == NULL) {
		NcDim * dimNumElem = ncMesh.get_dim("num_elem");
		if (dimNumElem == NULL) {
			_EXCEPTION1("Missing dimension \"num_elem\" in mesh file \"%s\"", strMeshFile.c_str());
		}
		int nElements = dimNumElem->size();
		vecDimSizes.push_back(nElements);
		vecDimNames.push_back("num_elem");
		return;
	}

	// Obtain rectilinear attributes (dimension sizes)
	NcAtt * attRectilinearDim0Size =
		ncMesh.get_att("rectilinear_dim0_size");
	NcAtt * attRectilinearDim1Size =
		ncMesh.get_att("rectilinear_dim1_size");

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
		ncMesh.get_att("rectilinear_dim0_name");
	NcAtt * attRectilinearDim1Name =
		ncMesh.get_att("rectilinear_dim1_name");

	if (attRectilinearDim0Name == NULL) {
		_EXCEPTIONT("Missing attribute \"rectilinear_dim0_name\"");
	}
	if (attRectilinearDim1Name == NULL) {
		_EXCEPTIONT("Missing attribute \"rectilinear_dim1_name\"");
	}

	std::string strDim0Name = attRectilinearDim0Name->as_string(0);
	std::string strDim1Name = attRectilinearDim1Name->as_string(0);

	// Push rectilinear attributes into array
	vecDimSizes.resize(2);
	vecDimSizes[0] = nDim0Size;
	vecDimSizes[1] = nDim1Size;

	vecDimNames.resize(2);
	vecDimNames[0] = strDim0Name;
	vecDimNames[1] = strDim1Name;
}

///////////////////////////////////////////////////////////////////////////////

void OfflineMap::InitializeSourceDimensionsFromFile(
	const std::string & strSourceMesh
) {
	InitializeDimensionsFromMeshFile(
		strSourceMesh,
		m_vecSourceDimNames,
		m_vecSourceDimSizes,
		m_dSourceCenterLon,
		m_dSourceCenterLat,
		m_dSourceVertexLon,
		m_dSourceVertexLat);
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
	InitializeDimensionsFromMeshFile(
		strTargetMesh,
		m_vecTargetDimNames,
		m_vecTargetDimSizes,
		m_dTargetCenterLon,
		m_dTargetCenterLat,
		m_dTargetVertexLon,
		m_dTargetVertexLat);
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
	DataArray1D<double> & dCenterLon,
	DataArray1D<double> & dCenterLat,
	DataArray2D<double> & dVertexLon,
	DataArray2D<double> & dVertexLat,
	bool fLatLon,
	int nNodesPerFace
) {
	// Check if already initialized
	if (dCenterLon.GetRows() != 0) {
		return;
	}

	int nFaces = mesh.faces.size();

	// Count maximum number of Nodes per Face
	if (nNodesPerFace == 0) {
    for (int i = 0; i < nFaces; i++) {
      if (mesh.faces[i].edges.size() > nNodesPerFace) {
        nNodesPerFace = mesh.faces[i].edges.size();
      }
    }
  }

	dVertexLon.Allocate(nFaces, nNodesPerFace);
	dVertexLat.Allocate(nFaces, nNodesPerFace);

	if ((fLatLon) && (nNodesPerFace != 4)) {
		_EXCEPTIONT("Logic error");
	}

	dCenterLon.Allocate(nFaces);
	dCenterLat.Allocate(nFaces);

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
	const DataArray3D<int> & dataGLLnodes,
	DataArray1D<double> & dCenterLon,
	DataArray1D<double> & dCenterLat,
	DataArray2D<double> & dVertexLon,
	DataArray2D<double> & dVertexLat
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

	dCenterLon.Allocate(iMaxNodeIx);
	dCenterLat.Allocate(iMaxNodeIx);

	dVertexLon.Allocate(iMaxNodeIx, 1);
	dVertexLat.Allocate(iMaxNodeIx, 1);

	DataArray1D<double> dG;
	GetDefaultNodalLocations(nP, dG);

	for (int i = 0; i < dataGLLnodes.GetRows(); i++) {
	for (int j = 0; j < dataGLLnodes.GetColumns(); j++) {
	for (int k = 0; k < dataGLLnodes.GetSubColumns(); k++) {
		const Face & face = mesh.faces[k];

		Node node;

		ApplyLocalMap(
			face,
			mesh.nodes,
			dG[j],
			dG[i],
			node);

		int iNode = dataGLLnodes[i][j][k] - 1;

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
	const DataArray2D<double> & dVertexLon,
	const DataArray2D<double> & dVertexLat,
	bool fLonFirst,
	DataArray1D<double> & dCenterLon,
	DataArray1D<double> & dCenterLat,
	DataArray1D<double> & dVectorCenterLon,
	DataArray1D<double> & dVectorCenterLat,
	DataArray2D<double> & dVectorBoundsLon,
	DataArray2D<double> & dVectorBoundsLat
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

	dCenterLon.Allocate(nFaces);
	dCenterLat.Allocate(nFaces);

	dVectorCenterLon.Allocate(nLon);
	dVectorCenterLat.Allocate(nLat);

	dVectorBoundsLon.Allocate(nLon, 2);
	dVectorBoundsLat.Allocate(nLat, 2);

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

	// Monotonize
	if ((dVectorCenterLon[1] > dVectorCenterLon[0]) &&
		(dVectorCenterLon[0] > dVectorCenterLon[nLon-1])
	) {
		for (int i = 0; i < nLon; i++) {
			if (dVectorCenterLon[i] > 180.0) {
				dVectorCenterLon[i] -= 360.0;
			}
			if (dVectorBoundsLon[i][0] > 180.0) {
				dVectorBoundsLon[i][0] -= 360.0;
			}
			if (dVectorBoundsLon[i][1] > 180.0) {
				dVectorBoundsLon[i][1] -= 360.0;
			}
		}
		for (int i = 0; i < nFaces; i++) {
			if (dCenterLon[i] > 180.0) {
				dCenterLon[i] -= 360.0;
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void OfflineMap::InitializeSourceCoordinatesFromMeshFV(
	const Mesh & meshSource
) {
	// Check if these arrays have been read from file
	if ((m_dSourceVertexLon.IsAttached()) ||
		(m_dSourceVertexLat.IsAttached()) ||
		(m_dSourceCenterLon.IsAttached()) ||
		(m_dSourceCenterLat.IsAttached())
	) {
		if ((m_dSourceVertexLon.IsAttached()) &&
			(m_dSourceVertexLat.IsAttached()) &&
			(m_dSourceCenterLon.IsAttached()) &&
			(m_dSourceCenterLat.IsAttached())
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
	if ((m_dTargetVertexLon.IsAttached()) ||
		(m_dTargetVertexLat.IsAttached()) ||
		(m_dTargetCenterLon.IsAttached()) ||
		(m_dTargetCenterLat.IsAttached())
	) {
		if ((m_dTargetVertexLon.IsAttached()) &&
			(m_dTargetVertexLat.IsAttached()) &&
			(m_dTargetCenterLon.IsAttached()) &&
			(m_dTargetCenterLat.IsAttached())
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
		m_dVectorTargetBoundsLon.Allocate(2);
		m_dVectorTargetBoundsLat.Allocate(2);
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
	const DataArray3D<int> & dataGLLnodesSource
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
	const DataArray3D<int> & dataGLLnodesTarget
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

	DataArray1D<float> dataIn;
	if (m_vecSourceDimSizes.size() == 1) {
		dataIn.Allocate(nSourceCount);
	} else {
		dataIn.Allocate(m_vecSourceDimSizes[0] * m_vecSourceDimSizes[1]);
	}

	DataArray1D<float> dataOut;
	if (m_vecTargetDimSizes.size() == 1) {
		dataOut.Allocate(nTargetCount);
	} else {
		dataOut.Allocate(m_vecTargetDimSizes[0] * m_vecTargetDimSizes[1]);
	}

	DataArray1D<double> dataInDouble(nSourceCount);
	DataArray1D<double> dataOutDouble(nTargetCount);

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

			NcBool fNoErr = varLon->put(
				&(m_dVectorTargetCenterLon[0]),
				m_dVectorTargetCenterLon.GetRows());

			if (!fNoErr) {
				_EXCEPTION1("Error writing \"lon\" to NetCDF file (%i)", NcError::get_err());
			}

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

			NcBool fNoErr = varLat->put(
				&(m_dVectorTargetCenterLat[0]),
				m_dVectorTargetCenterLat.GetRows());

			if (!fNoErr) {
				_EXCEPTION1("Error writing \"lon\" to NetCDF file (%i)", NcError::get_err());
			}

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
		float flFillValue = m_flFillValueOverride;
		double dFillValue = m_dFillValueOverride;
		for (int a = 0; a < var->num_atts(); a++) {
			NcAtt * att = var->get_att(a);
			if ((strcmp(att->name(), "_FillValue") == 0) ||
			    (strcmp(att->name(), "missing_value") == 0)
			) {
				if (att->type() == ncDouble) {
					dFillValue = att->as_double(0);
				} else if (att->type() == ncFloat) {
					flFillValue = att->as_float(0);
				} else {
					_EXCEPTION1("Invalid type for attribute \"%s\"", att->name());
				}
			}
		}

		// Construct an array of dimensions for this variable
		int nVarTotalEntries = 1;

		DataArray1D<NcDim *> vecDims(var->num_dims());

		DataArray1D<NcDim *> vecDimsOut;
		if (fTargetRectilinear) {

			if (fSourceRectilinear) {
				vecDimsOut.Allocate(var->num_dims());
			} else {
				vecDimsOut.Allocate(var->num_dims()+1);
			}

			vecDimsOut[vecDimsOut.GetRows()-2] = dim0;
			vecDimsOut[vecDimsOut.GetRows()-1] = dim1;

		} else {
			if (fSourceRectilinear) {
				if (var->num_dims() == 1) {
					_EXCEPTIONT("Expected rectilinear source data to be stored in 2D");
				}

				vecDimsOut.Allocate(var->num_dims()-1);
			} else {
				vecDimsOut.Allocate(var->num_dims());
			}

			vecDimsOut[vecDimsOut.GetRows()-1] = dim0;
		}

		int nFreeDims = var->num_dims() - m_vecSourceDimSizes.size();

		// Verify map is compatible with source variable
		if (fSourceRectilinear) {
			if (var->get_dim(var->num_dims()-2)->size() != m_vecSourceDimSizes[0]) {
				_EXCEPTION4("Error: Source variable \"%s\" has inconsistent dimension size "
					"for this map in dimension \"%s\".  Expected %i, found %i.",
					var->name(),
					var->get_dim(var->num_dims()-2)->name(),
					m_vecSourceDimSizes[0],
					var->get_dim(var->num_dims()-2)->size());
			}
			if (var->get_dim(var->num_dims()-1)->size() != m_vecSourceDimSizes[1]) {
				_EXCEPTION4("Error: Source variable \"%s\" has inconsistent dimension size "
					"for this map in dimension \"%s\".  Expected %i, found %i.",
					var->name(),
					var->get_dim(var->num_dims()-1)->name(),
					m_vecSourceDimSizes[1],
					var->get_dim(var->num_dims()-1)->size());
			}

		} else {
			if (var->get_dim(var->num_dims()-1)->size() != nSourceCount) {
				_EXCEPTION4("Error: Source variable \"%s\" has inconsistent dimension size "
					"for this map in dimension \"%s\".  Expected %i, found %i.",
					var->name(),
					var->get_dim(var->num_dims()-1)->name(),
					nSourceCount,
					var->get_dim(var->num_dims()-1)->size());
			}
		}

		// Add any missing dimension variables to target file
		DataArray1D<long> vecDimSizes(nFreeDims);

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

		// Fix type of _FillValue
		{
			NcAtt * attFillValue = varOut->get_att("_FillValue");
			if (attFillValue != NULL) {
				if ((attFillValue->type() == ncFloat) &&
				    (varOut->type() == ncDouble)
				) {
					double dFillValue = attFillValue->as_double(0);
					attFillValue->remove();
					varOut->add_att("_FillValue", dFillValue);

				} else if (
				    (attFillValue->type() == ncDouble) &&
				    (varOut->type() == ncFloat)
				) {
					float flFillValue = attFillValue->as_float(0);
					attFillValue->remove();
					varOut->add_att("_FillValue", flFillValue);

				} else if (attFillValue->type() != varOut->type()) {
					_EXCEPTION1("Type incompatibility between inherited attribute"
						" \"_FillValue\" and variable \"%s\"", varOut->name());
				}
			}
		}

		// Source and output counts
		DataArray1D<long> nCountsIn(vecDims.GetRows());
		DataArray1D<long> nCountsOut(vecDimsOut.GetRows());

		// Get size
		DataArray1D<long> nGet(nCountsIn.GetRows());
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
		DataArray1D<long> nPut(nCountsOut.GetRows());
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
			} else {
				var->get(&(dataInDouble[0]), &(nGet[0]));

				if (dFillValue != 0.0) {
					for (int i = 0; i < nSourceCount; i++) {
						if (dataInDouble[i] == dFillValue) {
							dataInDouble[i] = 0.0;
						}
					}
				}
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
				NcBool fNoErr = varOut->put(&(dataOutDouble[0]), &(nPut[0]));
				if (!fNoErr) {
					_EXCEPTION1("Error writing to NetCDF file (%i)", NcError::get_err());
				}

			} else {
				// Cast the data to float
				for (int i = 0; i < dataOut.GetRows(); i++) {
					dataOut[i] = static_cast<float>(dataOutDouble[i]);
				}

				// Write the data as float
				varOut->set_cur(&(nCountsOut[0]));
				NcBool fNoErr = varOut->put(&(dataOut[0]), &(nPut[0]));
				if (!fNoErr) {
					_EXCEPTION1("Error writing to NetCDF file (%i)", NcError::get_err());
				}

				//}
			}
		}
		AnnounceEndBlock(NULL);
	}
}

///////////////////////////////////////////////////////////////////////////////

void OfflineMap::Read(
	const std::string & strSource,
	std::map<std::string, std::string> * pmapAttributes,
	NcFile::FileFormat * peFileFormat
) {
	NcFile ncMap(strSource.c_str(), NcFile::ReadOnly);
	if (!ncMap.is_valid()) {
		_EXCEPTION1("Unable to open input map file \"%s\"",
			strSource.c_str());
	}

	// Read netcdf file format
	NcFile::FileFormat eFileFormat = ncMap.get_format();

	if (peFileFormat != NULL) {
		*peFileFormat = eFileFormat;
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

	m_dSourceCenterLat.Allocate(nA);
	m_dTargetCenterLat.Allocate(nB);

	m_dSourceCenterLon.Allocate(nA);
	m_dTargetCenterLon.Allocate(nB);

	m_dSourceVertexLat.Allocate(nA, nVA);
	m_dTargetVertexLat.Allocate(nB, nVB);

	m_dSourceVertexLon.Allocate(nA, nVA);
	m_dTargetVertexLon.Allocate(nB, nVB);

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

		m_dVectorTargetCenterLat.Allocate(dimLatB->size());
		m_dVectorTargetCenterLon.Allocate(dimLonB->size());

		varLatCB->get(&(m_dVectorTargetCenterLat[0]), dimLatB->size());
		varLonCB->get(&(m_dVectorTargetCenterLon[0]), dimLonB->size());

		NcVar * varLatBounds = ncMap.get_var("lat_bnds");
		NcVar * varLonBounds = ncMap.get_var("lon_bnds");

		m_dVectorTargetBoundsLat.Allocate(dimLatB->size(), 2);
		m_dVectorTargetBoundsLon.Allocate(dimLonB->size(), 2);

		varLatBounds->get(&(m_dVectorTargetBoundsLat[0][0]),
			dimLatB->size(), 2);
		varLonBounds->get(&(m_dVectorTargetBoundsLon[0][0]),
			dimLonB->size(), 2);
	}

	// Read areas
	m_dSourceAreas.Allocate(nA);
	m_dTargetAreas.Allocate(nB);

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

	// Read masks
	NcVar * varMaskA = ncMap.get_var("mask_a");
	if (varMaskA != NULL) {
		m_iSourceMask.Allocate(nA);
		varMaskA->get(&(m_iSourceMask[0]), nA);
	}

	NcVar * varMaskB = ncMap.get_var("mask_b");
	if (varMaskB != NULL) {
		m_iTargetMask.Allocate(nB);
		varMaskB->get(&(m_iTargetMask[0]), nB);
	}

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

	DataArray1D<int> vecRow(nS);
	DataArray1D<int> vecCol(nS);
	DataArray1D<double> vecS(nS);

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

	// Load file attributes
	if (pmapAttributes != NULL) {
		for (int a = 0; a < ncMap.num_atts(); a++) {
			NcAtt * attGlobal = ncMap.get_att(a);
			if (attGlobal == NULL) {
				_EXCEPTIONT("Logic error");
			}

			std::pair<std::string, std::string> prAttribute(
				attGlobal->name(),
				attGlobal->as_string(0));

			pmapAttributes->insert(prAttribute);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void OfflineMap::Write(
	const std::string & strTarget,
	const std::map<std::string, std::string> & mapAttributes,
	NcFile::FileFormat eOutputFormat
) {
	// Temporarily change error reporting
	NcError error_temp(NcError::verbose_fatal);

	// Open an output file
	NcFile ncMap(strTarget.c_str(), NcFile::Replace, NULL, 0, eOutputFormat);
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
	varAreaA->add_att("units", "steradians");

	NcVar * varAreaB = ncMap.add_var("area_b", ncDouble, dimNB);
	varAreaB->put(&(m_dTargetAreas[0]), nB);
	varAreaB->add_att("units", "steradians");

	// Write masks
	if (m_iSourceMask.IsAttached()) {
		NcVar * varMaskA = ncMap.add_var("mask_a", ncInt, dimNA);
		varMaskA->put(&(m_iSourceMask[0]), nA);
		varMaskA->add_att("units", "unitless");

		if (!m_iTargetMask.IsAttached()) {
			NcVar * varMaskB = ncMap.add_var("mask_b", ncInt, dimNB);
			DataArray1D<int> iTargetMaskTemp(nB);
			for (int i = 0; i < nB; i++) {
				iTargetMaskTemp[i] = 1;
			}
			varMaskB->put(&(iTargetMaskTemp[0]), nB);
			varMaskB->add_att("units", "unitless");
		}
	}

	if (m_iTargetMask.IsAttached()) {
		if (!m_iSourceMask.IsAttached()) {
			NcVar * varMaskA = ncMap.add_var("mask_a", ncInt, dimNA);
			DataArray1D<int> iSourceMaskTemp(nA);
			for (int i = 0; i < nA; i++) {
				iSourceMaskTemp[i] = 1;
			}
			varMaskA->put(&(iSourceMaskTemp[0]), nA);
			varMaskA->add_att("units", "unitless");
		}

		NcVar * varMaskB = ncMap.add_var("mask_b", ncInt, dimNB);
		varMaskB->put(&(m_iTargetMask[0]), nB);
		varMaskB->add_att("units", "unitless");
	}

	// Write SparseMatrix entries
	DataArray1D<int> vecRow;
	DataArray1D<int> vecCol;
	DataArray1D<double> vecS;

	m_mapRemap.GetEntries(vecRow, vecCol, vecS);

	// Calculate and write fractional coverage arrays
	{
		DataArray1D<double> dFracA(nA);
		DataArray1D<double> dFracB(nB);

		for (int i = 0; i < vecS.GetRows(); i++) {
			dFracA[vecCol[i]] += vecS[i] / m_dSourceAreas[vecCol[i]] * m_dTargetAreas[vecRow[i]];
			dFracB[vecRow[i]] += vecS[i];
		}

		NcVar * varFracA = ncMap.add_var("frac_a", ncDouble, dimNA);
		varFracA->put(&(dFracA[0]), nA);
		varFracA->add_att("name", "fraction of target coverage of source dof");
		varFracA->add_att("units", "unitless");

		NcVar * varFracB = ncMap.add_var("frac_b", ncDouble, dimNB);
		varFracB->put(&(dFracB[0]), nB);
		varFracB->add_att("name", "fraction of source coverage of target dof");
		varFracB->add_att("units", "unitless");
	}

	// Increment vecRow and vecCol
	for (int i = 0; i < vecRow.GetRows(); i++) {
		vecRow[i]++;
		vecCol[i]++;
	}

	// Write out data
	int nS = vecRow.GetRows();
	NcDim * dimNS = ncMap.add_dim("n_s", nS);

	NcVar * varRow = ncMap.add_var("row", ncInt, dimNS);
	varRow->add_att("name", "sparse matrix target dof index");
	varRow->add_att("first_index", "1");

	NcVar * varCol = ncMap.add_var("col", ncInt, dimNS);
	varCol->add_att("name", "sparse matrix source dof index");
	varCol->add_att("first_index", "1");

	NcVar * varS = ncMap.add_var("S", ncDouble, dimNS);
	varS->add_att("name", "sparse matrix coefficient");

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

void OfflineMap::SetTranspose(
	const OfflineMap & mapIn
) {
	m_dSourceAreas = mapIn.m_dTargetAreas;
	m_dTargetAreas = mapIn.m_dSourceAreas;

	m_iSourceMask = mapIn.m_iTargetMask;
	m_iTargetMask = mapIn.m_iSourceMask;

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

	DataArray1D<int> dataRows;
	DataArray1D<int> dataCols;
	DataArray1D<double> dataEntries;

	mapIn.m_mapRemap.GetEntries(dataRows, dataCols, dataEntries);
/*
	DataArray1D<double> dTargetCoverage(mapIn.m_dTargetAreas.GetRows());
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

int OfflineMap::IsConsistent(
	double dTolerance,
	const DataArray1D<int> & dataRows,
	const DataArray1D<int> & dataCols,
	const DataArray1D<double> & dataEntries,
	DataArray1D<double> * pdRowSums
) {
	int nCount = 0;

	// Calculate row sums
	if (m_mapRemap.GetRows() < 1) {
		_EXCEPTIONT("IsConservative() called on map with no rows");
	}

	bool fDeleteRowSums = false;
	if (pdRowSums == NULL) {
		pdRowSums = new DataArray1D<double>(m_mapRemap.GetRows());
		fDeleteRowSums = true;
	}
	DataArray1D<double> & dRowSums = (*pdRowSums);

	for (int i = 0; i < dataRows.GetRows(); i++) {
		dRowSums[dataRows[i]] += dataEntries[i];
	}

	// Verify all row sums are equal to 1
	for (int i = 0; i < dRowSums.GetRows(); i++) {
		if (fabs(dRowSums[i] - 1.0) > dTolerance) {
			nCount++;
			if (nCount <= OfflineMapWarningMessageCount) {
				Announce("OfflineMap is not consistent (row %i) [%1.15e != 1.0]",
					i+1, dRowSums[i]);
			}
		}
	}
	if (nCount > OfflineMapWarningMessageCount) {
		if (OfflineMapWarningMessageCount == 0) {
			Announce("OfflineMap is not consistent in %i dofs",
				nCount - OfflineMapWarningMessageCount);
		} else if (OfflineMapWarningMessageCount > 0) {
			Announce("OfflineMap is not consistent in %i more dofs",
				nCount - OfflineMapWarningMessageCount);
		}
	}

	if (fDeleteRowSums) {
		delete pdRowSums;
	}

	return nCount;
}

///////////////////////////////////////////////////////////////////////////////

int OfflineMap::IsConservative(
	double dTolerance,
	const DataArray1D<int> & dataRows,
	const DataArray1D<int> & dataCols,
	const DataArray1D<double> & dataEntries,
	DataArray1D<double> * pdColumnSums

) {
	int nCount = 0;

	// Calculate column sums
	if (m_mapRemap.GetColumns() < 1) {
		_EXCEPTIONT("IsConservative() called on map with no columns");
	}

	bool fDeleteColumnSums = false;
	if (pdColumnSums == NULL) {
		pdColumnSums = new DataArray1D<double>(m_dSourceAreas.GetRows());
		fDeleteColumnSums = true;
	}
	DataArray1D<double> & dColumnSums = (*pdColumnSums);

	if (dColumnSums.GetRows() != m_dSourceAreas.GetRows()) {
		_EXCEPTION2("Assertion failure: dColumnSums.GetRows() (%i) != m_dSourceAreas.GetRows() (%i)", dColumnSums.GetRows(), m_dSourceAreas.GetRows());
	}
	for (int i = 0; i < dataRows.GetRows(); i++) {
		dColumnSums[dataCols[i]] +=
			dataEntries[i] * m_dTargetAreas[dataRows[i]];
	}
	for (int i = 0; i < m_dSourceAreas.GetRows(); i++) {
		dColumnSums[i] /= m_dSourceAreas[i];
	}

	// Verify all column sums equal the input Jacobian
	for (int i = 0; i < dColumnSums.GetRows(); i++) {
		if (fabs(dColumnSums[i] - 1.0) > dTolerance) {
			nCount++;
			if (nCount <= OfflineMapWarningMessageCount) {
				Announce("OfflineMap is not conservative (col %i) [%1.15e != 1.0]",
					i+1, dColumnSums[i]);
			}
		}
	}
	if (nCount > OfflineMapWarningMessageCount) {
		if (OfflineMapWarningMessageCount == 0) {
			Announce("OfflineMap is not conservative in %i dofs",
				nCount - OfflineMapWarningMessageCount);
		} else if (OfflineMapWarningMessageCount > 0) {
			Announce("OfflineMap is not conservative in %i more dofs",
				nCount - OfflineMapWarningMessageCount);
		}
	}

	if (fDeleteColumnSums) {
		delete pdColumnSums;
	}

	return nCount;
}

///////////////////////////////////////////////////////////////////////////////

int OfflineMap::IsMonotone(
	double dTolerance,
	const DataArray1D<int> & dataRows,
	const DataArray1D<int> & dataCols,
	const DataArray1D<double> & dataEntries
) {
	int nCount = 0;

	// Verify all entries are in the range [0,1]
	for (int i = 0; i < dataRows.GetRows(); i++) {
		if (std::isnan(dataEntries[i])) {
			Announce("OfflineMap has NaN (s%i -> t%i)",
				dataCols[i]+1, dataRows[i]+1);
		}
		if ((dataEntries[i] < -dTolerance) ||
			(dataEntries[i] > 1.0 + dTolerance)
		) {
			nCount++;
			if (nCount <= OfflineMapWarningMessageCount) {
				Announce("OfflineMap is not monotone (s%i -> t%i) %1.15e",
					dataCols[i]+1, dataRows[i]+1, dataEntries[i]);
			}
		}
	}
	if (nCount > OfflineMapWarningMessageCount) {
		if (OfflineMapWarningMessageCount == 0) {
			Announce("OfflineMap is not monotone in %i dofs",
				nCount - OfflineMapWarningMessageCount);
		} else if (OfflineMapWarningMessageCount > 0) {
			Announce("OfflineMap is not monotone in %i more dofs",
				nCount - OfflineMapWarningMessageCount);
		}
	}

	return nCount;
}

///////////////////////////////////////////////////////////////////////////////

int OfflineMap::IsConsistent(
	double dTolerance
) {

	// Get map entries
	DataArray1D<int> dataRows;
	DataArray1D<int> dataCols;
	DataArray1D<double> dataEntries;

	m_mapRemap.GetEntries(dataRows, dataCols, dataEntries);

	return IsConsistent(dTolerance, dataRows, dataCols, dataEntries);
}

///////////////////////////////////////////////////////////////////////////////

int OfflineMap::IsConservative(
	double dTolerance
) {
	// Get map entries
	DataArray1D<int> dataRows;
	DataArray1D<int> dataCols;
	DataArray1D<double> dataEntries;

	m_mapRemap.GetEntries(dataRows, dataCols, dataEntries);

	return IsConservative(dTolerance, dataRows, dataCols, dataEntries);
}

///////////////////////////////////////////////////////////////////////////////

int OfflineMap::IsMonotone(
	double dTolerance
) {

	// Get map entries
	DataArray1D<int> dataRows;
	DataArray1D<int> dataCols;
	DataArray1D<double> dataEntries;

	m_mapRemap.GetEntries(dataRows, dataCols, dataEntries);

	return IsMonotone(dTolerance, dataRows, dataCols, dataEntries);
}

///////////////////////////////////////////////////////////////////////////////

bool OfflineMap::CheckMap(
	bool fCheckConsistency,
	bool fCheckConservation,
	bool fCheckMonotonicity,
	double dNormalTolerance,
	double dStrictTolerance,
	double dTotalOverlapArea
) {
	// Get map entries
	DataArray1D<int> dataRows;
	DataArray1D<int> dataCols;
	DataArray1D<double> dataEntries;

	m_mapRemap.GetEntries(dataRows, dataCols, dataEntries);

	// Verify at least one entry
	if (dataEntries.GetRows() == 0) {
		Announce("No entries found in map; aborting");
		return true;
	}
	if (m_mapRemap.GetRows() < 1) {
		_EXCEPTIONT("Assertion failure: m_mapRemap.GetRows() < 1");
	}
	if (m_mapRemap.GetColumns() < 1) {
		_EXCEPTIONT("Assertion failure: m_mapRemap.GetColumns() < 1");
	}
	if (m_dSourceAreas.GetRows() < m_mapRemap.GetColumns()) {
		_EXCEPTIONT("Assertion failure: m_dSourceAreas.GetRows() < m_mapRemap.GetColumns()");
	}
	if (m_dTargetAreas.GetRows() < m_mapRemap.GetRows()) {
		_EXCEPTIONT("Assertion failure: m_dTargetAreas.GetRows() < m_mapRemap.GetRows()");
	}

	DataArray1D<double> dRowSums(m_dTargetAreas.GetRows());
	DataArray1D<double> dColSums(m_dSourceAreas.GetRows());

	// Announce
	AnnounceBanner();
   	AnnounceStartBlock("Analyzing map");

	// Check consistency in individual cells
	int nConsistentFail = 0;
	if (fCheckConsistency) {
		AnnounceStartBlock("Per-dof consistency  (tol %1.5e)", dNormalTolerance);
		nConsistentFail =
			IsConsistent(dNormalTolerance, dataRows, dataCols, dataEntries, &dRowSums);
		if (nConsistentFail == 0) {
			AnnounceEndBlock("PASS");
		} else {
			AnnounceEndBlock(NULL);
		}
	} else {
		for (int i = 0; i < dataRows.GetRows(); i++) {
			dRowSums[dataRows[i]] += dataEntries[i];
		}
	}

	// Check conservation
	int nConservativeFail = 0;
	if (fCheckConservation) {
		AnnounceStartBlock("Per-dof conservation (tol %1.5e)", dNormalTolerance);
		nConservativeFail =
			IsConservative(dNormalTolerance, dataRows, dataCols, dataEntries, &dColSums);
		if (nConservativeFail == 0) {
			AnnounceEndBlock("PASS");
		} else {
			AnnounceEndBlock(NULL);
		}
	} else {
		if (m_dSourceAreas.GetRows() != dColSums.GetRows()) {
			_EXCEPTIONT("Assertion failure: m_dSourceAreas.GetRows() != dColSums.GetRows()");
		}
		for (int i = 0; i < dataRows.GetRows(); i++) {
			dColSums[dataCols[i]] +=
				dataEntries[i] * m_dTargetAreas[dataRows[i]];
		}
	}

	// Check monotonicity
	int nMonotoneFail = 0;
	if (fCheckMonotonicity) {
		AnnounceStartBlock("Per-dof monotonicity (tol %1.5e)", dStrictTolerance);
		nMonotoneFail =
			IsMonotone(dStrictTolerance, dataRows, dataCols, dataEntries);
		if (nMonotoneFail == 0) {
			AnnounceEndBlock("PASS");
		} else {
			AnnounceEndBlock(NULL);
		}

	// Check nominal range of entries
	} else {
		AnnounceStartBlock("Weights within range [-10,+10]");
		for (int i = 0; i < dataRows.GetRows(); i++) {
			if (std::isnan(dataEntries[i])) {
				Announce("OfflineMap has NaN (s%i -> t%i)",
					dataCols[i]+1, dataRows[i]+1);

			} else if ((dataEntries[i] < -10.0) || (dataEntries[i] > 10.0)) {
				Announce("OfflineMap has unusually large weight (s%i -> t%i) %1.15e",
					dataCols[i]+1, dataRows[i]+1, dataEntries[i]);
			}
		}
		AnnounceEndBlock("Done");
	}

	// Basic information
	{
		DataArray1D<int> nHistogramWeights(7);
		DataArray1D<int> nHistogramRows(32);
		DataArray1D<int> nHistogramCols(32);
		DataArray1D<int> nNonzeroRowCount(m_mapRemap.GetRows());
		DataArray1D<int> nNonzeroColCount(m_mapRemap.GetColumns());

		int iMinCol = dataCols[0];
		int iMaxCol = dataCols[0];

		int iMinRow = dataRows[0];
		int iMaxRow = dataRows[0];

		double dMinWeight = dataEntries[0];
		double dMaxWeight = dataEntries[0];

		for (int i = 0; i < dataRows.GetRows(); i++) {
			if (dataCols[i] > iMaxCol) {
				iMaxCol = dataCols[i];
			}
			if (dataCols[i] < iMinCol) {
				iMinCol = dataCols[i];
			}
			if (dataRows[i] > iMaxRow) {
				iMaxRow = dataRows[i];
			}
			if (dataRows[i] < iMinRow) {
				iMinRow = dataRows[i];
			}
			if (!std::isnan(dataEntries[i])) {
				if (dataEntries[i] < dMinWeight) {
					dMinWeight = dataEntries[i];
				}
				if (dataEntries[i] > dMaxWeight) {
					dMaxWeight = dataEntries[i];
				}
			}

			nNonzeroRowCount[dataRows[i]]++;
			nNonzeroColCount[dataCols[i]]++;

			if (dataEntries[i] < -10.0) {
				nHistogramWeights[0]++;
			} else if (dataEntries[i] < -1.0) {
				nHistogramWeights[1]++;
			} else if (dataEntries[i] < - dStrictTolerance) {
				nHistogramWeights[2]++;
			} else if (dataEntries[i] <= 1.0 + dStrictTolerance) {
				nHistogramWeights[3]++;
			} else if (dataEntries[i] < 2.0) {
				nHistogramWeights[4]++;
			} else if (dataEntries[i] < 10.0) {
				nHistogramWeights[5]++;
			} else {
				nHistogramWeights[6]++;
			}
		}

		for (int i = 0; i < nNonzeroRowCount.GetRows(); i++) {
			if (nNonzeroRowCount[i] < 31) {
				nHistogramRows[ nNonzeroRowCount[i] ]++;
			} else {
				nHistogramRows[31]++;
			}
		}

		for (int i = 0; i < nNonzeroColCount.GetRows(); i++) {
			if (nNonzeroColCount[i] < 31) {
				nHistogramCols[ nNonzeroColCount[i] ]++;
			} else {
				nHistogramCols[31]++;
			}
		}

		double dSourceMinArea = DBL_MAX;
		double dSourceMaxArea = -DBL_MAX;
		double dSourceArea = 0.0;
		for (int i = 0; i < m_dSourceAreas.GetRows(); i++) {
			dSourceArea += m_dSourceAreas[i];
			if (m_dSourceAreas[i] < dSourceMinArea) {
				dSourceMinArea = m_dSourceAreas[i];
			}
			if (m_dSourceAreas[i] > dSourceMaxArea) {
				dSourceMaxArea = m_dSourceAreas[i];
			}
		}

		double dTargetMinArea = DBL_MAX;
		double dTargetMaxArea = -DBL_MAX;
		double dTargetArea = 0.0;
		for (int i = 0; i < m_dTargetAreas.GetRows(); i++) {
			dTargetArea += m_dTargetAreas[i];
			if (m_dTargetAreas[i] < dTargetMinArea) {
				dTargetMinArea = m_dTargetAreas[i];
			}
			if (m_dTargetAreas[i] > dTargetMaxArea) {
				dTargetMaxArea = m_dTargetAreas[i];
			}
		}

		double dSourceMinFrac = DBL_MAX;
		double dSourceMaxFrac = -DBL_MAX;
		double dTargetMinFrac = DBL_MAX;
		double dTargetMaxFrac = -DBL_MAX;
		{
			DataArray1D<double> dSourceFrac(m_dSourceAreas.GetRows());
			DataArray1D<double> dTargetFrac(m_dTargetAreas.GetRows());

			for (int i = 0; i < dataEntries.GetRows(); i++) {
				dSourceFrac[dataCols[i]] +=
					dataEntries[i]
					/ m_dSourceAreas[dataCols[i]]
					* m_dTargetAreas[dataRows[i]];
				dTargetFrac[dataRows[i]] += dataEntries[i];
			}

			for (int i = 0; i < m_dSourceAreas.GetRows(); i++) {
				if (dSourceFrac[i] < dSourceMinFrac) {
					dSourceMinFrac = dSourceFrac[i];
				}
				if (dSourceFrac[i] > dSourceMaxFrac) {
					dSourceMaxFrac = dSourceFrac[i];
				}
			}

			for (int i = 0; i < m_dTargetAreas.GetRows(); i++) {
				if (dTargetFrac[i] < dTargetMinFrac) {
					dTargetMinFrac = dTargetFrac[i];
				}
				if (dTargetFrac[i] > dTargetMaxFrac) {
					dTargetMaxFrac = dTargetFrac[i];
				}
			}

		}

		int iSourceMinMask = INT_MAX;
		int iSourceMaxMask = INT_MIN;
		for (int i = 0; i < m_iSourceMask.GetRows(); i++) {
			if (m_iSourceMask[i] < iSourceMinMask) {
				iSourceMinMask = m_iSourceMask[i];
			}
			if (m_iSourceMask[i] > iSourceMaxMask) {
				iSourceMaxMask = m_iSourceMask[i];
			}
		}

		int iTargetMinMask = INT_MAX;
		int iTargetMaxMask = INT_MIN;
		for (int i = 0; i < m_iTargetMask.GetRows(); i++) {
			if (m_iTargetMask[i] < iTargetMinMask) {
				iTargetMinMask = m_iTargetMask[i];
			}
			if (m_iTargetMask[i] > iTargetMaxMask) {
				iTargetMaxMask = m_iTargetMask[i];
			}
		}

		double dRowSumMin = dRowSums[0];
		double dRowSumMax = dRowSums[0];
		for (int i = 1; i < dRowSums.GetRows(); i++) {
			if (dRowSums[i] < dRowSumMin) {
				dRowSumMin = dRowSums[i];
			}
			if (dRowSums[i] > dRowSumMax) {
				dRowSumMax = dRowSums[i];
			}
		}

		double dColSumMin = dColSums[0];
		double dColSumMax = dColSums[0];
		for (int i = 1; i < dColSums.GetRows(); i++) {
			if (dColSums[i] < dColSumMin) {
				dColSumMin = dColSums[i];
			}
			if (dColSums[i] > dColSumMax) {
				dColSumMax = dColSums[i];
			}
		}

		Announce("");
		Announce("  Total nonzero entries: %i", dataEntries.GetRows());
		Announce("   Column index min/max: %i / %i (%i source dofs)",
			iMinCol+1, iMaxCol+1, m_mapRemap.GetColumns());
		Announce("      Row index min/max: %i / %i (%i target dofs)",
			iMinRow+1, iMaxRow+1, m_mapRemap.GetRows());
		Announce("      Source area / 4pi: %1.15e", dSourceArea / (4.0 * M_PI));
		Announce("    Source area min/max: %1.15e / %1.15e",
			dSourceMinArea, dSourceMaxArea);
		if (dSourceMinArea < 0.0) {
			Announce("ERROR: Negative source area detected");
		}
		Announce("      Target area / 4pi: %1.15e", dTargetArea / (4.0 * M_PI));
		Announce("    Target area min/max: %1.15e / %1.15e",
			dTargetMinArea, dTargetMaxArea);
		if (dTargetMinArea < 0.0) {
			Announce("ERROR: Negative source area detected");
		}
		Announce("    Source frac min/max: %1.15e / %1.15e",
			dSourceMinFrac, dSourceMaxFrac);
		Announce("    Target frac min/max: %1.15e / %1.15e",
			dTargetMinFrac, dTargetMaxFrac);
		if (m_iSourceMask.GetRows() != 0) {
			Announce("    Source mask min/max: %i / %i",
				iSourceMinMask, iSourceMaxMask);
		}
		if (m_iTargetMask.GetRows() != 0) {
			Announce("    Target mask min/max: %i / %i",
				iTargetMinMask, iTargetMaxMask);
		}

		Announce("    Map weights min/max: %1.15e / %1.15e",
			dMinWeight, dMaxWeight);
		Announce("       Row sums min/max: %1.15e / %1.15e",
			dRowSumMin, dRowSumMax);
		Announce("   Consist. err min/max: %1.15e / %1.15e",
			dRowSumMin - 1.0, dRowSumMax - 1.0);
		Announce("    Col wt.sums min/max: %1.15e / %1.15e",
			dColSumMin, dColSumMax);
		Announce("   Conserv. err min/max: %1.15e / %1.15e",
			dColSumMin - 1.0, dColSumMax - 1.0);

		if (dTotalOverlapArea != 0.0) {
			double dGlobalConservationError = 0.0;
			for (int i = 0; i < dColSums.GetRows(); i++) {
				dGlobalConservationError += dColSums[i] * m_dSourceAreas[i];
			}
			Announce("    Global conserv. err: %1.15e",
				dGlobalConservationError - dTotalOverlapArea);
		}

		char szBuffer[128];

		Announce("");
		Announce("Histogram of nonzero entries in sparse matrix");
		Announce("..Column 1: Number of nonzero entries (bin minimum)");
		Announce("..Column 2: Number of columns with that many nonzero values");
		Announce("..Column 3: Number of rows with that many nonzero values");
		std::string strHistogram("[");
		for (int i = 0; i < nHistogramRows.GetRows(); i++) {
			if ((nHistogramCols[i] != 0) || (nHistogramRows[i] != 0)) {
				if (strHistogram.length() != 1) {
					strHistogram += ",";
				}
				if (i == nHistogramRows.GetRows()-1) {
					sprintf(szBuffer, "[%i, %i, %i]",
						i, nHistogramCols[i], nHistogramRows[i]);
				} else {
					sprintf(szBuffer, "[%i, %i, %i]",
						i, nHistogramCols[i], nHistogramRows[i]);
				}
				strHistogram += szBuffer;
			}
		}
		strHistogram += "]";
		Announce(strHistogram.c_str());

		Announce("");
		Announce("Histogram of weights");
		Announce("..Column 1: Lower bound on weights");
		Announce("..Column 2: Upper bound on weights");
		Announce("..Column 3: # of weights in that bin");
		strHistogram = "[";
		if (nHistogramWeights[0] != 0) {
			sprintf(szBuffer, "[-inf, -10, %i]", nHistogramWeights[0]);
			strHistogram += szBuffer;
		}
		if (nHistogramWeights[1] != 0) {
			if (strHistogram.length() != 1) {
				strHistogram += ",";
			}
			sprintf(szBuffer, "[-10, -1, %i]", nHistogramWeights[1]);
			strHistogram += szBuffer;
		}
		if (nHistogramWeights[2] != 0) {
			if (strHistogram.length() != 1) {
				strHistogram += ",";
			}
			sprintf(szBuffer, "[-1, 0, %i]", nHistogramWeights[2]);
			strHistogram += szBuffer;
		}
		if (nHistogramWeights[3] != 0) {
			if (strHistogram.length() != 1) {
				strHistogram += ",";
			}
			sprintf(szBuffer, "[0, 1, %i]", nHistogramWeights[3]);
			strHistogram += szBuffer;
		}
		if (nHistogramWeights[4] != 0) {
			if (strHistogram.length() != 1) {
				strHistogram += ",";
			}
			sprintf(szBuffer, "[1, 2, %i]", nHistogramWeights[4]);
			strHistogram += szBuffer;
		}
		if (nHistogramWeights[5] != 0) {
			if (strHistogram.length() != 1) {
				strHistogram += ",";
			}
			sprintf(szBuffer, "[2, 10, %i]", nHistogramWeights[5]);
			strHistogram += szBuffer;
		}
		if (nHistogramWeights[6] != 0) {
			if (strHistogram.length() != 1) {
				strHistogram += ",";
			}
			sprintf(szBuffer, "[10, inf, %i]", nHistogramWeights[6]);
			strHistogram += szBuffer;
		}
		strHistogram += "]";

		Announce("%s", strHistogram.c_str());

		bool fError = false;
		if (nHistogramCols[0] != 0) {
			if (!fError) Announce("");
			fError = true;
			Announce("NOTE: Some source mesh dofs are not captured by this map");
		}
		if (nHistogramRows[0] != 0) {
			if (!fError) Announce("");
			fError = true;
			Announce("NOTE: Some target mesh dofs are not captured by this map");
		}
		if (fabs(dSourceArea - 4.0 * M_PI) > dStrictTolerance) {
			if (!fError) Announce("");
			fError = true;
			Announce("NOTE: Source weights do not agree with sphere area (tol %1.5e)",
				dStrictTolerance);
		}
		if (fabs(dTargetArea - 4.0 * M_PI) > dStrictTolerance) {
			if (!fError) Announce("");
			fError = true;
			Announce("NOTE: Target weights do not agree with sphere area (tol %1.5e)",
				dStrictTolerance);
		}

	}

    AnnounceEndBlock(NULL);
	AnnounceBanner();

	return (
		(nConsistentFail == 0) &&
		(nConservativeFail == 0) &&
		(nMonotoneFail == 0));
}

///////////////////////////////////////////////////////////////////////////////

