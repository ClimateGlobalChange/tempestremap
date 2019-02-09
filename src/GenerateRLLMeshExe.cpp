///////////////////////////////////////////////////////////////////////////////
///
///	\file    GenerateRLLMeshExe.cpp
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

#include "TempestRemapAPI.h"

///////////////////////////////////////////////////////////////////////////////

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

	// Change the arrangement of latitudes to have half a latitude at the poles
	bool fGlobalCap;

	// Flip latitude and longitude dimension in FaceVector ordering
	bool fFlipLatLon;

    // Input filename
    std::string strInputFile;

    // Input mesh is global
    bool fForceGlobal;

    // Verbose output
    bool fVerbose;

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
	CommandLineBool(fGlobalCap, "global_cap");
    CommandLineBool(fFlipLatLon, "flip");
    CommandLineString(strInputFile, "in_file", "");
    CommandLineBool(fForceGlobal, "in_global");
    CommandLineBool(fVerbose, "verbose");
    CommandLineString(strOutputFile, "file", "outRLLMesh.g");

    ParseCommandLine(argc, argv);
    EndCommandLine(argv)
/*
    // Verify latitude box is increasing
    if (dLatBegin >= dLatEnd) {
        _EXCEPTIONT("--lat_begin and --lat_end must specify a positive interval");
    }
    if (dLonBegin >= dLonEnd) {
        _EXCEPTIONT("--lon_begin and --lon_end must specify a positive interval");
    }
 */
    std::cout << "=========================================================";
    std::cout << std::endl;

	// Call the actual mesh generator
    Mesh mesh;
    int err = GenerateRLLMesh(
		mesh,
		nLongitudes, nLatitudes,
		dLonBegin, dLonEnd,
		dLatBegin, dLatEnd,
		fGlobalCap,
		fFlipLatLon,
		fForceGlobal,
		strInputFile, strOutputFile,
		fVerbose);
	if (err) exit(err);

	return 0;
}

///////////////////////////////////////////////////////////////////////////////
