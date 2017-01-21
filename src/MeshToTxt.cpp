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
#include "netcdfcpp.h"

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	NcError error(NcError::silent_nonfatal);

try {

	// Input mesh file
	std::string strInputMesh;

	// Output node file
	std::string strOutputNodes;

	// Output face file
	std::string strOutputFaces;

	// Pad faces
	bool fPad;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputMesh, "in", "");
		CommandLineString(strOutputNodes, "out_nodes", "nodes.dat");
		CommandLineString(strOutputFaces, "out_faces", "faces.dat");
		CommandLineBool(fPad, "pad");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	if (strInputMesh == "") {
		return (-1);
	}

	AnnounceBanner();

	// Load input mesh
	AnnounceStartBlock("Loading input mesh");
	Mesh meshInput(strInputMesh);
	AnnounceEndBlock(NULL);

	// Output nodes
	AnnounceStartBlock("Writing nodes");
	FILE * fpNodes = fopen(strOutputNodes.c_str(), "w");
	for (int i = 0; i < meshInput.nodes.size(); i++) {
		fprintf(fpNodes, "%1.10e %1.10e %1.10e\n",
			static_cast<double>(meshInput.nodes[i].x),
			static_cast<double>(meshInput.nodes[i].y),
			static_cast<double>(meshInput.nodes[i].z));
	}
	fclose(fpNodes);
	AnnounceEndBlock("Done!");

	// Maximum number of nodes per face
	int nMaximumNodes = 0;
	if (fPad) {
		for (int i = 0; i < meshInput.faces.size(); i++) {
			if (meshInput.faces[i].edges.size() > nMaximumNodes) {
				nMaximumNodes = meshInput.faces[i].edges.size();
			}
		}
	}

	// Output faces
	AnnounceStartBlock("Writing faces");
	FILE * fpFaces = fopen(strOutputFaces.c_str(), "w");
	for (int i = 0; i < meshInput.faces.size(); i++) {
		for (int j = 0; j < meshInput.faces[i].edges.size(); j++) {
			fprintf(fpFaces, "%i", meshInput.faces[i][j] + 1);
			if (j != meshInput.faces[i].edges.size()-1) {
				fprintf(fpFaces, " ");
			}
		}
		if (fPad) {
			int j = meshInput.faces[i].edges.size();
			for (; j < nMaximumNodes; j++) {
				fprintf(fpFaces, " %i", meshInput.faces[i][0] + 1);
			}
		}
		fprintf(fpFaces, "\n");
	}
	fclose(fpFaces);
	AnnounceEndBlock("Done!");

	AnnounceBanner();

	return (0);

} catch(Exception & e) {
	Announce(e.ToString().c_str());
	return (-1);

} catch(...) {
	return (-2);
}
}

///////////////////////////////////////////////////////////////////////////////

