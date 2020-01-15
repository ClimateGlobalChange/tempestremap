///////////////////////////////////////////////////////////////////////////////
///
///	\file    GenerateGLLMetaDataExe.cpp
///	\author  Paul Ullrich
///	\version August 12, 2014
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
#include "Exception.h"
#include "GridElements.h"
#include "DataArray3D.h"

#include "TempestRemapAPI.h"

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	// Input mesh
	std::string strMesh;

	// Polynomial order
	int nP;

	// Use of bubble to adjust areas
	bool fNoBubble;

	// Output metadata file
	std::string strOutput;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strMesh, "mesh", "");
		CommandLineInt(nP, "np", 4);
		CommandLineBool(fNoBubble, "no_bubble");
		CommandLineString(strOutput, "out", "gllmeta.nc");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Calculate metadata
	DataArray3D<int> dataGLLnodes;
	DataArray3D<double> dataGLLJacobian;
    Mesh mesh;

	int err =
		GenerateGLLMetaData(
			strMesh,
			mesh,
			nP,
			fNoBubble,
			strOutput,
			dataGLLnodes,
			dataGLLJacobian);

	if (err) exit(err);

	// Done
	AnnounceBanner();

	return 0;
}

///////////////////////////////////////////////////////////////////////////////
