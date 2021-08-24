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
#include "STLStringHelper.h"

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

	// Overlap grid generation algorithm
	std::string strAlgorithm;

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
		CommandLineStringD(strAlgorithm, "alg", "kdx", "(edge|kdx|lint)");
		CommandLineBool(fNoValidate, "novalidate");
		CommandLineBool(fHasConcaveFacesA, "concavea");
		CommandLineBool(fHasConcaveFacesB, "concaveb");
		CommandLineBool(fAllowNoOverlap, "allow_no_overlap");
		CommandLineBool(fVerbose, "verbose");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Change algorithm to lowercase
	STLStringHelper::ToLower(strAlgorithm);

	// Call the actual mesh generator
	Mesh meshOverlap;

	// Edge algorithm
	if (strAlgorithm == "edge") {
		if (fHasConcaveFacesA) {
			Announce("WARNING: --alg \"edge\" does not support concave faces.  Argument --concavea ignored.");
		}
		if (fHasConcaveFacesB) {
			Announce("WARNING: --alg \"edge\" does not support concave faces.  Argument --concaveb ignored.");
		}
		if (fAllowNoOverlap) {
			Announce("WARNING: Argument --allow_no_overlap has no effect with --alg \"edge\"");
		}

		int err =
			GenerateOverlapMeshEdge(
				strMeshA, strMeshB,
 				meshOverlap, strOverlapMesh, strOutputFormat,
				strMethod, fNoValidate);

		if (err) exit(err);

	// KD-tree search
	} else if (strAlgorithm == "kdx" ) {
		int err =
			GenerateOverlapMeshKdx(
				strMeshA, strMeshB,
 				meshOverlap, strOverlapMesh, strOutputFormat,
				strMethod, fNoValidate,
				fHasConcaveFacesA, fHasConcaveFacesB,
				fAllowNoOverlap,
				fVerbose);

		if (err) exit(err);

	// Line search with interval tree
	} else if (strAlgorithm == "lint" ) {
		if (fAllowNoOverlap) {
			Announce("WARNING: Argument --allow_no_overlap has no effect with --alg \"lint\"");
		}

		int err =
			GenerateOverlapMeshLint(
				strMeshA, strMeshB,
 				meshOverlap, strOverlapMesh, strOutputFormat,
				strMethod, fNoValidate,
				fHasConcaveFacesA, fHasConcaveFacesB,
				fAllowNoOverlap,
				fVerbose);

		if (err) exit(err);

	} else {
		_EXCEPTIONT("Invalid value of --alg, expected \"edge\", \"kdx\" or \"lint\"");
	}

	AnnounceBanner();

	return 0;
}

///////////////////////////////////////////////////////////////////////////////
