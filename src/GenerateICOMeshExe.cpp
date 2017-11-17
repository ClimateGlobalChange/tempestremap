///////////////////////////////////////////////////////////////////////////////
///
///	\file    GenerateICOMeshExe.cpp
///	\author  Paul Ullrich
///	\version September 9, 2014
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

#include "TempestRemapAPI.h"

#include <cmath>
#include <iostream>

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {
    
	// Resolution
	int nResolution;

	// Dual mesh
	bool fDual;

	// Output filename
	std::string strOutputFile;

	// Parse the command line
	BeginCommandLine()
		CommandLineInt(nResolution, "res", 10);
		CommandLineBool(fDual, "dual");
		CommandLineString(strOutputFile, "file", "outICOMesh.g");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	// Call the actual mesh generator
    Mesh mesh;
	int err = GenerateICOMesh(mesh, nResolution, fDual, strOutputFile);
	if (err) exit(err);
	else return 0;
}

///////////////////////////////////////////////////////////////////////////////
