///////////////////////////////////////////////////////////////////////////////
///
///	\file    GenerateRectilinearMeshFromFile.cpp
///	\author  Paul Ullrich
///	\version March 31, 2021
///
///	<remarks>
///		Copyright 2021 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#include "GridElements.h"
#include "DataArray2D.h"
#include "Exception.h"
#include "Announce.h"
#include "STLStringHelper.h"
#include "NetCDFUtilities.h"
#include "PolynomialInterp.h"

#include <cmath>
#include <cfloat>
#include <iostream>

#include "netcdfcpp.h"

///////////////////////////////////////////////////////////////////////////////

extern "C" 
int GenerateRectilinearMeshFromFile(
	Mesh & mesh, 
	std::string strInputFile,
	std::string strInputFileLonName,
	std::string strInputFileLatName,
	std::string strOutputFile, 
	std::string strOutputFormat,
	bool fVerbose
) {

	NcError error(NcError::silent_nonfatal);

try {

    // Output format
    STLStringHelper::ToLower(strOutputFormat);

	NcFile::FileFormat eOutputFormat =
		GetNcFileFormatFromString(strOutputFormat);
	if (eOutputFormat == NcFile::BadFormat) {
		_EXCEPTION1("Invalid \"out_format\" value (%s), "
			"expected [Classic|Offset64Bits|Netcdf4|Netcdf4Classic]",
			strOutputFormat.c_str());
	}

	// Load input file
	NcFile ncinfile(strInputFile.c_str(), NcFile::ReadOnly);
	if (!ncinfile.is_valid()) {
		_EXCEPTION1("Unable to load input file \"%s\"", strInputFile.c_str());
	}

	NcVar * varLon = ncinfile.get_var(strInputFileLonName.c_str());
	if (varLon == NULL) {
		_EXCEPTION1("Unable to find variable \"%s\" in input file", strInputFileLonName.c_str());
	}

	NcVar * varLat = ncinfile.get_var(strInputFileLatName.c_str());
	if (varLat == NULL) {
		_EXCEPTION1("Unable to find variable \"%s\" in input file", strInputFileLatName.c_str());
	}

	// The mesh
	mesh.Clear();

	// 2D latitude and longitude arrays
	if (varLon->num_dims() == 2) {
		if (varLat->num_dims() != 2) {
			_EXCEPTION2("Longitude variable \"%s\" and latitude variable \"%s\" must have same dimension count",
				strInputFileLonName.c_str(),
				strInputFileLatName.c_str());
		}

		AnnounceStartBlock("Using 2D latitude-longitude arrays to build mesh");

		if ((std::string(varLon->get_dim(0)->name()) != std::string(varLat->get_dim(0)->name())) ||
		    (std::string(varLon->get_dim(1)->name()) != std::string(varLat->get_dim(1)->name()))
		) {
			_EXCEPTION2("Longitude variable \"%s\" and latitude variable \"%s\" must have same dimensions",
				strInputFileLonName.c_str(),
				strInputFileLatName.c_str());
		}

		long lDim0Size = varLon->get_dim(0)->size();
		long lDim1Size = varLon->get_dim(1)->size();

		if (lDim0Size == 1) {
			_EXCEPTIONT("Longitude variable dimension 0 must have size greater than 1");
		}
		if (lDim1Size == 1) {
			_EXCEPTIONT("Longitude variable dimension 1 must have size greater than 1");
		}

		double dFillValue = DBL_MAX;
		NcAtt * attFillValue = varLon->get_att("_FillValue");
		if (attFillValue != NULL) {
			dFillValue = attFillValue->as_double(0);
		}

		DataArray2D<double> dLon(lDim0Size, lDim1Size);
		DataArray2D<double> dLat(lDim0Size, lDim1Size);

		varLon->get(&(dLon(0,0)), lDim0Size, lDim1Size);
		varLat->get(&(dLat(0,0)), lDim0Size, lDim1Size);

		std::vector<double> dFitX1;
		std::vector<double> dFitLon1;
		std::vector<double> dFitLat1;
		std::vector<double> dFitX2;
		std::vector<double> dFitLon2;
		std::vector<double> dFitLat2;

		dFitX1.reserve(4);
		dFitLon1.reserve(4);
		dFitLat1.reserve(4);
		dFitX2.reserve(4);
		dFitLon2.reserve(4);
		dFitLat2.reserve(4);

		// Verify longitude and latitude are consistent
		AnnounceStartBlock("Counting faces");
		int nFaceCount = 0;
		for (long j = 0; j < lDim0Size; j++) {
		for (long i = 0; i < lDim1Size; i++) {
			if ((dLon(j,i) != dFillValue) && (dLat(j,i) == dFillValue)) {
				_EXCEPTION2("At index (%li,%li) only one of latitude and longitude array is defined", i, j);
			}
			if ((dLon(j,i) == dFillValue) && (dLat(j,i) != dFillValue)) {
				_EXCEPTION2("At index (%li,%li) only one of latitude and longitude array is defined", i, j);
			}
			nFaceCount++;
		}
		}
		Announce("%i faces found", nFaceCount);
		AnnounceEndBlock(NULL);

		// Mesh
		mesh.faces.reserve(nFaceCount);

		// Use linear, quadratic or cubic interpolation to estimate longitude at vertex
		AnnounceStartBlock("Estimating vertex locations");

		DataArray2D<int> ixV(lDim0Size+1, lDim1Size+1);
		//DataArray2D<double> dLonV(lDim0Size+1, lDim1Size+1);
		//DataArray2D<double> dLatV(lDim0Size+1, lDim1Size+1);

		int ixVertexCount = 0;
		for (long j = 0; j < lDim0Size+1; j++) {
		for (long i = 0; i < lDim1Size+1; i++) {

			bool fVertexNeeded = false;
			if ((i != lDim1Size) && (j != lDim0Size) && (dLon(j,i) != dFillValue)) {
				fVertexNeeded = true;
			} else if ((i != lDim1Size) && (j != 0) && (dLon(j-1,i) != dFillValue)) {
				fVertexNeeded = true;
			} else if ((i != 0) && (j != lDim0Size) && (dLon(j,i-1) != dFillValue)) {
				fVertexNeeded = true;
			} else if ((i != 0) && (j != 0) && (dLon(j-1,i-1) != dFillValue)) {
				fVertexNeeded = true;
			}
			if (!fVertexNeeded) {
				ixV(j,i) = InvalidNode;
				continue;
			}

			long ixbegin = (i<2)?(-i):(-2);
			long ixend   = (i>=lDim1Size-2)?(lDim1Size-i):(2);

			long jxbegin = (j<2)?(-j):(-2);
			long jxend   = (j>=lDim0Size-2)?(lDim0Size-j):(2);

			dFitX2.clear();
			dFitLon2.clear();
			dFitLat2.clear();
			for (long ix = ixbegin; ix < ixend; ix++) {
				dFitX1.clear();
				dFitLon1.clear();
				dFitLat1.clear();
				for (long jx = jxbegin; jx < jxend; jx++) {
					if (dLon(j+jx,i+ix) != dFillValue) {
						dFitX1.push_back(static_cast<double>(jx)+0.5);
						dFitLon1.push_back(dLon(j+jx,i+ix));
						dFitLat1.push_back(dLat(j+jx,i+ix));
					}
				}
				if (dFitX1.size() < 2) {
					continue;
				}

				dFitX2.push_back(static_cast<double>(ix)+0.5);
				dFitLon2.push_back(PolynomialInterp::Interpolate(dFitX1.size(), &(dFitX1[0]), &(dFitLon1[0]), 0.0));
				dFitLat2.push_back(PolynomialInterp::Interpolate(dFitX1.size(), &(dFitX1[0]), &(dFitLat1[0]), 0.0));
			}

			if (dFitX2.size() < 2) {
				_EXCEPTION2("Insufficient data to determine vertex coordinate at vertex (%li,%li)", i, j);
			}

			//dLonV(j,i) = PolynomialInterp::Interpolate(dFitX2.size(), &(dFitX2[0]), &(dFitLon2[0]), 0.0);
			//dLatV(j,i) = PolynomialInterp::Interpolate(dFitX2.size(), &(dFitX2[0]), &(dFitLat2[0]), 0.0);
			double dLonVDeg = PolynomialInterp::Interpolate(dFitX2.size(), &(dFitX2[0]), &(dFitLon2[0]), 0.0);
			double dLatVDeg = PolynomialInterp::Interpolate(dFitX2.size(), &(dFitX2[0]), &(dFitLat2[0]), 0.0);

			double dLonVRad = dLonVDeg * M_PI / 180.0;
			double dLatVRad = dLatVDeg * M_PI / 180.0;

			double dX = cos(dLatVRad) * cos(dLonVRad);
			double dY = cos(dLatVRad) * sin(dLonVRad);
			double dZ = sin(dLatVRad);

			ixV(j,i) = ixVertexCount;
			mesh.nodes.push_back(Node(dX,dY,dZ));

			ixVertexCount++;
		}
		}

		AnnounceEndBlock("Done");

		// Build the mesh
		AnnounceStartBlock("Building mesh");

		long lTotalFaces = 0;

		for (long j = 0; j < lDim0Size; j++) {
		for (long i = 0; i < lDim1Size; i++) {

			if (dLon(j,i) == dFillValue) {
				continue;
			}

			lTotalFaces++;

			int iGivenValues =
				  ((ixV(j  ,i  ) != InvalidNode)?(1):(0))
				+ ((ixV(j+1,i  ) != InvalidNode)?(1):(0))
				+ ((ixV(j  ,i+1) != InvalidNode)?(1):(0))
				+ ((ixV(j+1,i+1) != InvalidNode)?(1):(0));

			if (iGivenValues != 4) {
				_EXCEPTION2("Error estimating longitude and latitude at vertex of face (%li,%li)", i, j);
			}

			// Check for orientation
			bool fCounterclockwise = true;

			int ixVf[4];
			ixVf[0] = ixV(j,i);
			ixVf[1] = ixV(j,i+1);
			ixVf[2] = ixV(j+1,i+1);
			ixVf[3] = ixV(j+1,i);
			for (int k = 0; k < 4; k++) {
				int i0 = ixVf[(k+3)%4];
				int i1 = ixVf[k];
				int i2 = ixVf[(k+1)%4];

				double dDX0 = mesh.nodes[i0].x - mesh.nodes[i1].x;
				double dDY0 = mesh.nodes[i0].y - mesh.nodes[i1].y;
				double dDZ0 = mesh.nodes[i0].z - mesh.nodes[i1].z;

				double dDX1 = mesh.nodes[i2].x - mesh.nodes[i1].x;
				double dDY1 = mesh.nodes[i2].y - mesh.nodes[i1].y;
				double dDZ1 = mesh.nodes[i2].z - mesh.nodes[i1].z;

				double dCrossX = dDY0 * dDZ1 - dDZ0 * dDY1;
				double dCrossY = dDZ0 * dDX1 - dDX0 * dDZ1;
				double dCrossZ = dDX0 * dDY1 - dDY0 * dDX1;

				double dDot = dCrossX * mesh.nodes[i1].x + dCrossY * mesh.nodes[i1].y + dCrossZ * mesh.nodes[i1].z;

				if (fabs(dDot) < ReferenceTolerance) {
					printf("\n\nNodes (i0,i1,i2): (%i,%i,%i)\n", i0, i1, i2);
					mesh.nodes[i0].Print("i0");
					mesh.nodes[i1].Print("i1");
					mesh.nodes[i2].Print("i2");
					_EXCEPTION3("Nearly parallel grid lines detected in face (%li,%li) edge %li", i, j, k);
				}

				if (k == 0) {
					if (dDot > 0.0) {
						fCounterclockwise = false;
					}
				} else {
					if (((dDot < 0.0) && !fCounterclockwise) || ((dDot > 0.0) && fCounterclockwise)) {
						_EXCEPTION3("Orientation error in face (%li,%li) edge %li; possible concave face", i, j, k);
					}
				}
			}

			Face face(4);
			if (fCounterclockwise) {
				face.SetNode(0, ixV(j  ,i  ));
				face.SetNode(1, ixV(j  ,i+1));
				face.SetNode(2, ixV(j+1,i+1));
				face.SetNode(3, ixV(j+1,i  ));
			} else {
				face.SetNode(0, ixV(j  ,i  ));
				face.SetNode(1, ixV(j+1,i  ));
				face.SetNode(2, ixV(j+1,i+1));
				face.SetNode(3, ixV(j  ,i+1));
			}
			mesh.faces.push_back(face);
		}
		}

		mesh.vecGridDimSize.resize(2);
		mesh.vecGridDimSize[0] = lDim0Size;
		mesh.vecGridDimSize[1] = lDim1Size;

		mesh.vecGridDimName.resize(2);
		mesh.vecGridDimName[0] = varLon->get_dim(0)->name();
		mesh.vecGridDimName[1] = varLon->get_dim(1)->name();

		AnnounceEndBlock("Done");

		// Store rectilinear coordinates if all points are available
		if (lTotalFaces == lDim0Size * lDim1Size) {
			mesh.vecGridDimSize.resize(2);
			mesh.vecGridDimSize[0] = lDim0Size;
			mesh.vecGridDimSize[1] = lDim1Size;

			mesh.vecGridDimName.resize(2);
			mesh.vecGridDimName[0] = "dim0";
			mesh.vecGridDimName[1] = "dim1";
		}

		AnnounceEndBlock("Done");

	} else {
		_EXCEPTIONT("Not implemented:  At present longitude and latitude variable must have 2 dimensions");
	}

	// Write mesh
	{
		AnnounceStartBlock("Writing mesh");

		mesh.Write(strOutputFile, eOutputFormat);

		AnnounceEndBlock("Done");
	}

	AnnounceBanner();

	return 0;

} catch(Exception & e) {
	Announce(e.ToString().c_str());
	return (0);

} catch(...) {
	return (0);
}
}

///////////////////////////////////////////////////////////////////////////////

