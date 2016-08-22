///////////////////////////////////////////////////////////////////////////////
///
///	\file    GenerateLambertConfConicMesh.cpp
///	\author  Paul Ullrich
///	\version November 17, 2014
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

int GenerateLambertConfConicMesh(int argc, char** argv) {

	NcError error(NcError::silent_nonfatal);

try {
	// Number of columns in mesh
	int nNCol;

	// Number of rows in mesh
	int nNRow;

	// Reference longitude
	double dLon0;

	// Reference latitude
	double dLat0;

	// First standard parallel
	double dLat1;

	// Second standard parallel
	double dLat2;

	// Meters to bottom-left X position
	double dXLL;

	// Meters to bottom-left Y position
	double dYLL;

	// Cell size
	double dDX;

	// Output filename
	std::string strOutputFile;

	// Parse the command line
	BeginCommandLine()
		CommandLineInt(nNCol, "ncol", 5268);
		CommandLineInt(nNRow, "nrow", 4823);
		CommandLineDouble(dLon0, "lon0", -100.0);
		CommandLineDouble(dLat0, "lat0", 42.5);
		CommandLineDouble(dLat1, "lat1", 25.0);
		CommandLineDouble(dLat2, "lat2", 60.0);
		CommandLineDoubleD(dXLL,  "xll", -2015000.0, "(meters)");
		CommandLineDoubleD(dYLL,  "yll", 1785000.0, "(meters)");
		CommandLineDoubleD(dDX,   "dx", 1000.0, "(meters)");
		CommandLineString(strOutputFile, "file", "outLCCMesh.g");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	// Verify latitude box is increasing
	if (dLat1 >= dLat2) {
		_EXCEPTIONT("--lat1 must be less than --lat2");
	}
	if (dLat0 <= dLat1) {
		_EXCEPTIONT("--lat0 must be larger than --lat1");
	}
	if (dLat0 >= dLat2) {
		_EXCEPTIONT("--lat0 must be less than --lat2");
	}

	// Announce
	AnnounceBanner();

	// Convert latitude and longitude to radians
	dLon0 *= M_PI / 180.0;
	dLat0 *= M_PI / 180.0;
	dLat1 *= M_PI / 180.0;
	dLat2 *= M_PI / 180.0;

	// Convert XLL, YLL and DX to Earth radii
	dXLL /= 6.37122e6;
	dYLL /= 6.37122e6;
	dDX  /= 6.37122e6;

	// Calculate N, F and Rho0
	double dN = log(cos(dLat1) / cos(dLat2))
		/ log(tan(0.25 * M_PI + 0.5 * dLat2) / tan(0.25 * M_PI + 0.5 * dLat1));

	double dF = cos(dLat1) * pow(tan(0.25 * M_PI + 0.5 * dLat1), dN) / dN;

	double dRho0 = dF * pow(1.0 / tan(0.25 * M_PI + 0.5 * dLat0), dN);

	// Generate the mesh
	Mesh mesh;

	NodeVector & nodes = mesh.nodes;
	FaceVector & faces = mesh.faces;

	// Announce
	AnnounceStartBlock("Distributing nodes");

	// Add all nodal locations
	for (int i = 0; i <= nNRow; i++) {
	for (int j = 0; j <= nNCol; j++) {

		double dXX = dXLL + dDX * static_cast<double>(i);
		double dYY = dYLL + dDX * static_cast<double>(j);

		double dTheta = atan2(dXX, dRho0 - dYY);
		double dRho = dN / fabs(dN)
			* sqrt(dXX * dXX + (dRho0 - dYY) * (dRho0 - dYY));

		double dLambda = dLon0 + dTheta / dN;
		double dPhi = 2.0 * atan( pow(dF / dRho, 1.0/dN)) - 0.5 * M_PI;

		if ((i == 0) || (i == nNRow)) {
			if ((j == 0) || (j == nNCol)) {
				Announce("Corner: %3.3f %3.3f",
					dLambda * 180.0 / M_PI, dPhi * 180.0 / M_PI);
/*
				double dRho2 = dF * pow(1/tan(0.25 * M_PI + 0.5 * dPhi), dN);
				double dXX2 = dRho2 * sin(dN * (dLambda - dLon0));
				double dYY2 = dRho0 - dRho2 * cos(dN * (dLambda - dLon0));

				printf("%1.5e %1.5e : %1.5e %1.5e\n",
					dXX, dYY, dXX2, dYY2);
*/
			}
		}

		double dX = cos(dPhi) * cos(dLambda);
		double dY = cos(dPhi) * sin(dLambda);
		double dZ = sin(dPhi);

		nodes.push_back(Node(dX, dY, dZ));
	}
	}
	return (-1);

	// Announce
	AnnounceEndBlock("Done");
	AnnounceStartBlock("Assigning faces");

	// Add all faces 
	for (int j = 0; j < nNRow; j++) {
		int iThisYNodeIx =  j    * (nNCol + 1);
		int iNextYNodeIx = (j+1) * (nNCol + 1);

		for (int i = 0; i < nNCol; i++) {
			Face face(4);
			face.SetNode(0, iThisYNodeIx + i);
			face.SetNode(1, iThisYNodeIx + (i + 1));
			face.SetNode(2, iNextYNodeIx + (i + 1));
			face.SetNode(3, iNextYNodeIx + i);

			faces.push_back(face);
		}
	}

	AnnounceEndBlock("Done");

/*
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
	if (!fIncludeSouthPole) {
		nLatitudes++;
	}

	int iInteriorLatBegin = (fIncludeSouthPole)?(1):(1);
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
	for (int j = 1; j < iInteriorLatEnd; j++) {
		int iThisLatNodeIx = (j-1) * nLongitudeNodes + iSouthPoleOffset;
		int iNextLatNodeIx =  j    * nLongitudeNodes + iSouthPoleOffset;

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
		int iThisLatNodeIx = (nLatitudes - 2) * nLongitudeNodes + iSouthPoleOffset;
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
*/
	// Announce
	Announce("Writing mesh to file [%s]", strOutputFile.c_str());

	// Output the mesh
	mesh.Write(strOutputFile);

	// Add rectilinear properties
	NcFile ncOutput(strOutputFile.c_str(), NcFile::Write);
	ncOutput.add_att("rectilinear", "true");
	ncOutput.add_att("rectilinear_dim0_size", nNRow);
	ncOutput.add_att("rectilinear_dim1_size", nNCol);
	ncOutput.add_att("rectilinear_dim0_name", "y");
	ncOutput.add_att("rectilinear_dim1_name", "x");
	ncOutput.close();

	// Announce
	Announce("Mesh generator exited successfully");
	AnnounceBanner();

	return (0);

} catch(Exception & e) {
	Announce(e.ToString().c_str());
	return (-1);

} catch(...) {
	return (-2);
}
	return 0;
}

///////////////////////////////////////////////////////////////////////////////

