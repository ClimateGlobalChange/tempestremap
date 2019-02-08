///////////////////////////////////////////////////////////////////////////////
///
///	\file    GenerateOverlapMesh.cpp
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

#include "Announce.h"
#include "CommandLine.h"
#include "GridElements.h"

#include "TempestRemapAPI.h"

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	// Input mesh A
	std::string strMeshA;

	// Input mesh B
	std::string strMeshB;

	// Output mesh file
	std::string strOverlapMesh;

	// Output format
	std::string strOutputFormat;

	// Overlap grid generation method
	std::string strMethod;

	// No validation of the meshes
	bool fNoValidate;

	// Concave elements may be present in mesh A
	bool fHasConcaveFacesA;

	// Concave elements may be present in mesh B
	bool fHasConcaveFacesB;

	// Allow for the case of no overlap element found (enabling this may
	// cause the mesh generator to generate an incomplete overlap mesh)
	bool fAllowNoOverlap;

	// Verbose
	bool fVerbose;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strMeshA, "a", "");
		CommandLineString(strMeshB, "b", "");
		CommandLineString(strOverlapMesh, "out", "overlap.g");
		CommandLineString(strOutputFormat, "out_format", "netcdf4");
		CommandLineStringD(strMethod, "method", "fuzzy", "(fuzzy|exact|mixed)");
		CommandLineBool(fNoValidate, "novalidate");
		CommandLineBool(fHasConcaveFacesA, "concavea");
		CommandLineBool(fHasConcaveFacesB, "concaveb");
		CommandLineBool(fAllowNoOverlap, "allow_no_overlap");
		CommandLineBool(fVerbose, "verbose");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Call the actual mesh generator
    Mesh meshOverlap;
    int err =
		GenerateOverlapMesh(
			strMeshA, strMeshB,
 			meshOverlap, strOverlapMesh, strOutputFormat,
			strMethod, fNoValidate,
			fHasConcaveFacesA, fHasConcaveFacesB,
			fAllowNoOverlap,
			fVerbose);

	if (err) exit(err);

	AnnounceBanner();

	return 0;
}

///////////////////////////////////////////////////////////////////////////////
