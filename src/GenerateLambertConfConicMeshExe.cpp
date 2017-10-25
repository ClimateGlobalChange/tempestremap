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

#include "TempestRemapAPI.h"

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

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

	// Announce
	AnnounceBanner();

	// Calculate metadata
    Mesh mesh;
    int err = GenerateLambertConfConicMesh(mesh, nNCol, nNRow, dLon0, dLat0, dLat1, dLat2, dXLL, dYLL, dDX, strOutputFile);
	if (err) exit(err);

	// Done
	AnnounceBanner();

	return 0;
}

///////////////////////////////////////////////////////////////////////////////
