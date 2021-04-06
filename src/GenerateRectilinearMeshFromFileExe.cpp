///////////////////////////////////////////////////////////////////////////////
///
///	\file	GenerateRectilinearMeshFromFileExe.cpp
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

#include "CommandLine.h"
#include "GridElements.h"
#include "Exception.h"
#include "Announce.h"

#include <cmath>
#include <iostream>

#include "TempestRemapAPI.h"

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	// Input file
	std::string strInputFile;

	// Input file longitude name
	std::string strInputFileLonName;

	// Input file latitude name
	std::string strInputFileLatName;

	// Output filename
	std::string strOutputFile;

	// Output format
	std::string strOutputFormat;

	// Verbose output
	bool fVerbose;

	// Parse the command line
	BeginCommandLine()
	CommandLineString(strInputFile, "in_file", "");
	CommandLineString(strInputFileLonName, "in_file_lon", "lon");
	CommandLineString(strInputFileLatName, "in_file_lat", "lat");
	CommandLineString(strOutputFile, "out_file", "outMesh.g");
	CommandLineString(strOutputFormat, "out_format", "Netcdf4");
	CommandLineBool(fVerbose, "verbose");

	ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Call the actual mesh generator
	Mesh mesh;
	int err = GenerateRectilinearMeshFromFile(
		mesh,
		strInputFile,
		strInputFileLonName,
		strInputFileLatName,
		strOutputFile,
		strOutputFormat,
		fVerbose);
	if (err) exit(err);

	return 0;
}

///////////////////////////////////////////////////////////////////////////////
