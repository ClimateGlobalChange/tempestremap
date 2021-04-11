///////////////////////////////////////////////////////////////////////////////
///
///	\file	RestructureDataExe.cpp
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

	// Variable name to restructure
	std::string strVariable;

	// FillValue
	std::string strFillValue;

	// Reference file (containing lat/lon information)
	std::string strRefFile;

	// Reference file longitude name
	std::string strRefFileLonName;

	// Reference file latitude name
	std::string strRefFileLatName;

	// Output filename
	std::string strOutputFile;

	// Output format
	std::string strOutputFormat;

	// Verbose output
	bool fVerbose;

	// Parse the command line
	BeginCommandLine()
	CommandLineString(strInputFile, "in_file", "");
	CommandLineString(strVariable, "var", "");
	CommandLineString(strFillValue, "fillvalue", "");
	CommandLineString(strRefFile, "ref_file", "");
	CommandLineString(strRefFileLonName, "ref_file_lon", "lon");
	CommandLineString(strRefFileLatName, "ref_file_lat", "lat");
	CommandLineString(strOutputFile, "out_file", "");
	CommandLineString(strOutputFormat, "out_format", "Netcdf4");
	CommandLineBool(fVerbose, "verbose");
	ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Call the actual mesh generator
	Mesh mesh;
	int err = RestructureData(
		strInputFile,
		strVariable,
		strFillValue,
		strRefFile,
		strRefFileLonName,
		strRefFileLatName,
		strOutputFile, 
		strOutputFormat,
		fVerbose);
	if (err) exit(err);

	return 0;
}

///////////////////////////////////////////////////////////////////////////////

