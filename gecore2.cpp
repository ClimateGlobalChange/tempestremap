///////////////////////////////////////////////////////////////////////////////
///
///	\file    gecore2.cpp
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
#include "Exception.h"
#include "GridElements.h"
#include "OverlapMesh.h"

#include <cmath>

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

try {

	// Input mesh file
	std::string strInputMesh;

	// Output mesh file
	std::string strOutputMesh;

	// Input data file
	std::string strInputData;

	// Output data file
	std::string strOutputData;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputMesh, "in_mesh", "");
		CommandLineString(strOutputMesh, "out_mesh", "");
		CommandLineString(strInputData, "in_data", "");
		CommandLineString(strOutputData, "out_data", "");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Load input mesh
	AnnounceStartBlock("Loading input mesh");
	Mesh meshInput(strInputMesh);
	meshInput.RemoveZeroEdges();
	AnnounceEndBlock(NULL);

	// Validate mesh
	AnnounceStartBlock("Validate input mesh");
	meshInput.Validate();
	AnnounceEndBlock(NULL);

	// Load output mesh
	AnnounceStartBlock("Loading output mesh");
	Mesh meshOutput(strOutputMesh);
	meshOutput.RemoveZeroEdges();
	AnnounceEndBlock(NULL);

	// Validate mesh
	AnnounceStartBlock("Validate output mesh");
	meshOutput.Validate();
	AnnounceEndBlock(NULL);

	// Construct the edge map on both meshes
	AnnounceStartBlock("Constructing edge map on input mesh");
	meshInput.ConstructEdgeMap();
	AnnounceEndBlock(NULL);

	AnnounceStartBlock("Constructing edge map on output mesh");
	meshOutput.ConstructEdgeMap();
	AnnounceEndBlock(NULL);

	// Construct the reverse node array on both meshes
	AnnounceStartBlock("Constructing reverse node array on input mesh");
	meshInput.ConstructReverseNodeArray();
	AnnounceEndBlock(NULL);

	AnnounceStartBlock("Constructing reverse node array on output mesh");
	meshOutput.ConstructReverseNodeArray();
	AnnounceEndBlock(NULL);

#pragma message "Need to deal with coincident nodes on First and Second mesh"
	// Construct the overlap mesh
	Mesh meshOverlap;

	AnnounceStartBlock("Construct overlap mesh");
	GenerateOverlapMesh(meshInput, meshOutput, meshOverlap);
	AnnounceEndBlock(NULL);

	// Write the overlap mesh
	AnnounceStartBlock("Writing overlap mesh");
	meshOverlap.Write("overlap.g");
	AnnounceEndBlock(NULL);

	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}
}

///////////////////////////////////////////////////////////////////////////////

