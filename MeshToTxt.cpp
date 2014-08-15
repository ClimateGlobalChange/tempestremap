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

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

try {

	// Input mesh file
	std::string strInputMesh;

	// Output node file
	std::string strOutputNodes;

	// Output face file
	std::string strOutputFaces;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputMesh, "in", "");
		CommandLineString(strOutputNodes, "out_nodes", "nodes.dat");
		CommandLineString(strOutputFaces, "out_faces", "faces.dat");

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
		fprintf(fpFaces, "\n");
	}
	fclose(fpFaces);
	AnnounceEndBlock("Done!");

	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}
}

///////////////////////////////////////////////////////////////////////////////

