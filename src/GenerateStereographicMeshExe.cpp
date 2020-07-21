///////////////////////////////////////////////////////////////////////////////
///
///	\file    GenerateStereographicMeshExe.cpp
///	\author  Paul Ullrich
///	\version February 13, 2020
///
///	<remarks>
///		Copyright 2020 Paul Ullrich
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

#include "TempestRemapAPI.h"

#include <cmath>
#include <iostream>

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	// Longitude of initial point
	double dLonDeg0;

	// Latitude of initial point
	double dLatDeg0;

	// Longitude of final point
	double dLonDeg1;

	// Latitude of final point
	double dLatDeg1;

	// X resolution
	int nResX;

	// Y resolution
	int nResY;

	// Cornerpoints are actually centroids
	bool fCentroids;

	// Longitude of the pole
	double dLonDegP;

	// Latitude of the pole
	double dLatDegP;

	// Output filename
	std::string strOutputFile;

	// NetCDF format
	std::string strOutputFormat;

	// Parse the command line
	BeginCommandLine()
		CommandLineDouble(dLonDegP, "lonp", 0.0);
		CommandLineDouble(dLatDegP, "latp", 90.0);
		CommandLineDouble(dLonDeg0, "lon0", 0.0);
		CommandLineDouble(dLatDeg0, "lat0", 0.0);
		CommandLineDouble(dLonDeg1, "lon1", 0.0);
		CommandLineDouble(dLatDeg1, "lat1", 90.0);
		CommandLineInt(nResX, "resx", 10);
		CommandLineInt(nResY, "resy", 10);
		CommandLineBool(fCentroids, "centroids");
		CommandLineString(strOutputFile, "file", "outStereographicMesh.g");
		CommandLineString(strOutputFormat, "out_format", "Netcdf4");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	// Call the actual mesh generator
    Mesh mesh;
	int err = GenerateStereographicMesh(
		mesh,
		dLonDegP, dLatDegP,
		dLonDeg0, dLatDeg0,
		dLonDeg1, dLatDeg1,
		nResX, nResY,
		fCentroids,
		strOutputFile,
		strOutputFormat);

	if (err) exit(err);
	else return 0;
}

///////////////////////////////////////////////////////////////////////////////
