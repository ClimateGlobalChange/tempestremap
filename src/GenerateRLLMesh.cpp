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
#include "GridElements.h"
#include "Exception.h"
#include "Announce.h"

#include <cmath>
#include <iostream>

#include "netcdfcpp.h"

///////////////////////////////////////////////////////////////////////////////

#define ONLY_GREAT_CIRCLES


///////////////////////////////////////////////////////////////////////////////
// 
// Input Parameters:
// Number of longitudes in mesh: int nLongitudes;
// Number of latitudes in mesh: int nLatitudes;
// First longitude line on mesh: double dLonBegin;
// Last longitude line on mesh: double dLonEnd;
// First latitude line on mesh: double dLatBegin;
// Last latitude line on mesh: double dLatEnd;
// Flip latitude and longitude dimension in FaceVector ordering: bool fFlipLatLon;
// Output filename:  std::string strOutputFile;
// 
// Output Parameters: Mesh*
// 
extern "C" Mesh* GenerateRLLMesh(int nLongitudes, int nLatitudes, double dLonBegin, double dLonEnd, double dLatBegin, double dLatEnd, bool fFlipLatLon, std::string strOutputFile) {

	NcError error(NcError::silent_nonfatal);

	// Generate the mesh
	Mesh *mesh = new Mesh();

try {

	// Verify latitude box is increasing
	if (dLatBegin >= dLatEnd) {
		_EXCEPTIONT("--lat_begin and --lat_end must specify a positive interval");
	}
	if (dLonBegin >= dLonEnd) {
		_EXCEPTIONT("--lon_begin and --lon_end must specify a positive interval");
	}

	// Announce
	std::cout << "=========================================================";
	std::cout << std::endl;
	std::cout << "..Generating mesh with resolution [";
	std::cout << nLongitudes << ", " << nLatitudes << "]" << std::endl;
	std::cout << "..Longitudes in range [";
	std::cout << dLonBegin << ", " << dLonEnd << "]" << std::endl;
	std::cout << "..Latitudes in range [";
	std::cout << dLatBegin << ", " << dLatEnd << "]" << std::endl;
	std::cout << std::endl;

	// Convert latitude and longitude interval to radians
	dLonBegin *= M_PI / 180.0;
	dLonEnd   *= M_PI / 180.0;
	dLatBegin *= M_PI / 180.0;
	dLatEnd   *= M_PI / 180.0;

	// Check parameters
	if (nLatitudes < 2) {
		std::cout << "Error: At least 2 latitudes are required." << std::endl;
		return NULL;
	}
	if (nLongitudes < 2) {
		std::cout << "Error: At least 2 longitudes are required." << std::endl;
		return NULL;
	}

	NodeVector & nodes = mesh->nodes;
	FaceVector & faces = mesh->faces;
    mesh->type = Mesh::MeshType_RLL;

	// Change in longitude
	double dDeltaLon = dLonEnd - dLonBegin;
	double dDeltaLat = dLatEnd - dLatBegin;

	// Check if longitudes wrap
	bool fWrapLongitudes = false;
	if (fmod(dDeltaLon, 2.0 * M_PI) < 1.0e-12) {
		fWrapLongitudes = true;
	}
	bool fIncludeSouthPole = (fabs(dLatBegin + 0.5 * M_PI) < 1.0e-12);
	bool fIncludeNorthPole = (fabs(dLatEnd   - 0.5 * M_PI) < 1.0e-12);

	int iSouthPoleOffset = (fIncludeSouthPole)?(1):(0);

	// Increase number of latitudes if south pole is not included
	int iInteriorLatBegin = (fIncludeSouthPole)?(1):(0);
	int iInteriorLatEnd   = (fIncludeNorthPole)?(nLatitudes-1):(nLatitudes);

	// Number of longitude nodes
	int nLongitudeNodes = nLongitudes;
	if (!fWrapLongitudes) {
		nLongitudeNodes++;
	}

	// Generate nodes
	if (fIncludeSouthPole) {
		nodes.push_back(Node(0.0, 0.0, -1.0));
	}
	for (int j = iInteriorLatBegin; j < iInteriorLatEnd+1; j++) {
		for (int i = 0; i < nLongitudeNodes; i++) {
			Real dPhiFrac =
				  static_cast<Real>(j)
				/ static_cast<Real>(nLatitudes);

			Real dLambdaFrac =
				  static_cast<Real>(i)
				/ static_cast<Real>(nLongitudes);

			Real dPhi = dDeltaLat * dPhiFrac + dLatBegin;
			Real dLambda = dDeltaLon * dLambdaFrac + dLonBegin;

			Real dX = cos(dPhi) * cos(dLambda);
			Real dY = cos(dPhi) * sin(dLambda);
			Real dZ = sin(dPhi);

			nodes.push_back(Node(dX, dY, dZ));
		}
	}
	if (fIncludeNorthPole) {
		nodes.push_back(Node(0.0, 0.0, +1.0));
	}

	// Generate south polar faces
	if (fIncludeSouthPole) {
		for (int i = 0; i < nLongitudes; i++) {
			Face face(4);
			face.SetNode(0, 0);
			face.SetNode(1, (i+1) % nLongitudeNodes + 1);
			face.SetNode(2, i + 1);
			face.SetNode(3, 0);

#ifndef ONLY_GREAT_CIRCLES
			face.edges[1].type = Edge::Type_ConstantLatitude;
			face.edges[3].type = Edge::Type_ConstantLatitude;
#endif

			faces.push_back(face);
		}
	}

	// Generate interior faces
	for (int j = iInteriorLatBegin; j < iInteriorLatEnd; j++) {
		int jx = j - iInteriorLatBegin;

		int iThisLatNodeIx =  jx    * nLongitudeNodes + iSouthPoleOffset;
		int iNextLatNodeIx = (jx+1) * nLongitudeNodes + iSouthPoleOffset;

		for (int i = 0; i < nLongitudes; i++) {
			Face face(4);
			face.SetNode(0, iThisLatNodeIx + (i + 1) % nLongitudeNodes);
			face.SetNode(1, iNextLatNodeIx + (i + 1) % nLongitudeNodes);
			face.SetNode(2, iNextLatNodeIx + i);
			face.SetNode(3, iThisLatNodeIx + i);

#ifndef ONLY_GREAT_CIRCLES
			face.edges[1].type = Edge::Type_ConstantLatitude;
			face.edges[3].type = Edge::Type_ConstantLatitude;
#endif

			faces.push_back(face);
		}
	}

	// Generate north polar faces
	if (fIncludeNorthPole) {
		int jx = nLatitudes - iInteriorLatBegin - 1;

		int iThisLatNodeIx =  jx    * nLongitudeNodes + iSouthPoleOffset;
		int iNextLatNodeIx = (jx+1) * nLongitudeNodes + iSouthPoleOffset;

		int iNorthPolarNodeIx = static_cast<int>(nodes.size()-1);
		for (int i = 0; i < nLongitudes; i++) {
			Face face(4);
			face.SetNode(0, iNorthPolarNodeIx);
			face.SetNode(1, iThisLatNodeIx + i);
			face.SetNode(2, iThisLatNodeIx + (i + 1) % nLongitudeNodes);
			face.SetNode(3, iNorthPolarNodeIx);

#ifndef ONLY_GREAT_CIRCLES
			face.edges[1].type = Edge::Type_ConstantLatitude;
			face.edges[3].type = Edge::Type_ConstantLatitude;
#endif

			faces.push_back(face);
		}
	}

	// Reorder the faces
	if (fFlipLatLon) {
		FaceVector faceTemp = mesh->faces;
		mesh->faces.clear();
		for (int i = 0; i < nLongitudes; i++) {
		for (int j = 0; j < nLatitudes; j++) {
			mesh->faces.push_back(faceTemp[j * nLongitudes + i]);
		}
		}
	}

	// Output the mesh
	if (strOutputFile.size()) {

		// Announce
		std::cout << "..Writing mesh to file [" << strOutputFile.c_str() << "] ";
		std::cout << std::endl;

		mesh->Write(strOutputFile);

		// Add rectilinear properties
		if (!fIncludeSouthPole) {
			nLatitudes--;
		}

		NcFile ncOutput(strOutputFile.c_str(), NcFile::Write);
		ncOutput.add_att("rectilinear", "true");

		if (fFlipLatLon) {
			ncOutput.add_att("rectilinear_dim0_size", nLongitudes);
			ncOutput.add_att("rectilinear_dim1_size", nLatitudes);
			ncOutput.add_att("rectilinear_dim0_name", "lon");
			ncOutput.add_att("rectilinear_dim1_name", "lat");
		} else {
			ncOutput.add_att("rectilinear_dim0_size", nLatitudes);
			ncOutput.add_att("rectilinear_dim1_size", nLongitudes);
			ncOutput.add_att("rectilinear_dim0_name", "lat");
			ncOutput.add_att("rectilinear_dim1_name", "lon");
		}
		ncOutput.close();
	}

	// Announce
	std::cout << "..Mesh generator exited successfully" << std::endl;
	std::cout << "=========================================================";
	std::cout << std::endl;

} catch(Exception & e) {
	Announce(e.ToString().c_str());
	return (0);

} catch(...) {
	return (0);
}
	return mesh;
}

///////////////////////////////////////////////////////////////////////////////

#ifdef TEMPEST_DRIVER_MODE

int main(int argc, char** argv) {
    // Number of longitudes in mesh
	int nLongitudes;

	// Number of latitudes in mesh
	int nLatitudes;

	// First longitude line on mesh
	double dLonBegin;

	// Last longitude line on mesh
	double dLonEnd;

	// First latitude line on mesh
	double dLatBegin;

	// Last latitude line on mesh
	double dLatEnd;

	// Flip latitude and longitude dimension in FaceVector ordering
	bool fFlipLatLon;

	// Output filename
	std::string strOutputFile;

	// Parse the command line
	BeginCommandLine()
		CommandLineInt(nLongitudes, "lon", 128);
		CommandLineInt(nLatitudes, "lat", 64);
		CommandLineDouble(dLonBegin, "lon_begin", 0.0);
		CommandLineDouble(dLonEnd, "lon_end", 360.0);
		CommandLineDouble(dLatBegin, "lat_begin", -90.0);
		CommandLineDouble(dLatEnd, "lat_end", 90.0);
		CommandLineBool(fFlipLatLon, "flip");
		CommandLineString(strOutputFile, "file", "outRLLMesh.g");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	// Call the actual mesh generator
	Mesh* mesh = GenerateRLLMesh(nLongitudes, nLatitudes, dLonBegin, dLonEnd, dLatBegin, dLatEnd, fFlipLatLon, strOutputFile);
	if (mesh) delete mesh;
	else return (-1);

	return 0;
}

#endif

///////////////////////////////////////////////////////////////////////////////
