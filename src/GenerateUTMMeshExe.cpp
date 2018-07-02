///////////////////////////////////////////////////////////////////////////////
///
///	\file    GenerateUTMMesh.cpp
///	\author  Paul Ullrich
///	\version July 2, 2018
///
///	<remarks>
///		Copyright 2000-2018 Paul Ullrich
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

	// UTM zone
	int nZone;

	// Number of columns in mesh
	int nCols;

	// Number of rows in mesh
	int nRows;

	// XLL Corner of the mesh
	double dXLLCorner;

	// YLL Corner of the mesh
	double dYLLCorner;

	// Cell size (in meters)
	double dCellSize;

	// Output filename
	std::string strOutputFile;

	// Verbose output
	bool fVerbose;

	// Parse the command line
	BeginCommandLine()
	CommandLineInt(nZone, "zone", 0);
	CommandLineInt(nCols, "cols", 128);
	CommandLineInt(nRows, "rows", 64);
	CommandLineDouble(dXLLCorner, "xllcorner", 0.0);
	CommandLineDouble(dYLLCorner, "yllcorner", 0.0);
	CommandLineDouble(dCellSize, "cellsize", 1000.0);
	CommandLineString(strOutputFile, "file", "outUTMMesh.g");
	CommandLineBool(fVerbose, "verbose");

	ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	// Verify input arguments
	if (nCols <= 0) {
		_EXCEPTIONT("--cols must be positive");
	}
	if (nRows <= 0) {
		_EXCEPTIONT("--rows must be positive");
	}
	if (dCellSize <= 0.0) {
		_EXCEPTIONT("--cellsize must be positive");
	}

	std::cout << "=========================================================";
	std::cout << std::endl;

	// Call the actual mesh generator
	Mesh mesh;
	int err = GenerateUTMMesh(mesh, nZone, nCols, nRows, dXLLCorner, dYLLCorner, dCellSize, strOutputFile, fVerbose);
	if (err) exit(err);

	return 0;
}

///////////////////////////////////////////////////////////////////////////////
