///////////////////////////////////////////////////////////////////////////////
///
///	\file    GenerateOverlapMeshParallelExe.cpp
///	\author  Paul Ullrich
///	\version August 23, 2021
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

#include "Announce.h"
#include "CommandLine.h"
#include "GridElements.h"
#include "STLStringHelper.h"

#include "TempestRemapAPI.h"

#if defined(TEMPEST_MPIOMP)
#include <mpi.h>
#endif

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

#if defined(TEMPEST_MPIOMP)
	MPI_Init(&argc, &argv);
	AnnounceOnlyOutputOnRankZero();
#endif

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

	// Run in parallel
	bool fParallel;

	// Temporary file directory
	std::string strTempDir;

	// Verbose
	bool fVerbose;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strMeshA, "a", "");
		CommandLineString(strMeshB, "b", "");
		CommandLineString(strOverlapMesh, "out", "overlap.g");
		CommandLineString(strOutputFormat, "out_format", "netcdf4");
		CommandLineStringD(strMethod, "method", "fuzzy", "(fuzzy|exact|mixed)");
		CommandLineStringD(strAlgorithm, "alg", "lint", "(lint)");
		CommandLineBool(fNoValidate, "novalidate");
		CommandLineBool(fHasConcaveFacesA, "concavea");
		CommandLineBool(fHasConcaveFacesB, "concaveb");
		CommandLineBool(fAllowNoOverlap, "allow_no_overlap");
		CommandLineString(strTempDir, "tmpdir", "/tmp");
		CommandLineBool(fVerbose, "verbose");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

#if !defined(TEMPEST_MPIOMP)
	_EXCEPTIONT("TempestRemap was not compiled with MPI: Parallel operation not supported");
#endif

	// Change algorithm to lowercase
	STLStringHelper::ToLower(strAlgorithm);

	// Call the actual mesh generator
	Mesh meshOverlap;

	// Edge algorithm
	if (strAlgorithm == "lint" ) {
		if (fAllowNoOverlap) {
			Announce("WARNING: Argument --allow_no_overlap has no effect with --alg \"lint\"");
		}

		int err =
			GenerateOverlapMeshLint(
				strMeshA, strMeshB,
 				meshOverlap, strOverlapMesh, strOutputFormat,
				strMethod, fNoValidate,
				fHasConcaveFacesA, fHasConcaveFacesB,
				true, strTempDir,
				fVerbose);

		if (err) exit(err);

	} else {
		_EXCEPTIONT("Invalid value of --alg, expected \"lint\"");
	}

	AnnounceBanner();

#if defined(TEMPEST_MPIOMP)
	MPI_Finalize();
#endif

	return 0;
}

///////////////////////////////////////////////////////////////////////////////
