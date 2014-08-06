///////////////////////////////////////////////////////////////////////////////
///
///	\file    GenerateRLLMesh.cpp
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

#include <cmath>
#include <iostream>

///////////////////////////////////////////////////////////////////////////////

#define ONLY_GREAT_CIRCLES

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

try {
	// Number of longitudes in mesh
	int nLongitudes;

	// Number of latitudes in mesh
	int nLatitudes;

	// Output filename
	std::string strOutputFile;

	// Parse the command line
	BeginCommandLine()
		CommandLineInt(nLongitudes, "lon", 128);
		CommandLineInt(nLatitudes, "lat", 64);
		CommandLineString(strOutputFile, "file", "outRLLMesh.g");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	// Announce
	std::cout << "=========================================================";
	std::cout << std::endl;
	std::cout << "..Generating mesh with resolution [";
	std::cout << nLongitudes << "," << nLatitudes << "]";
	std::cout << std::endl;

	// Check parameters
	if (nLatitudes < 2) {
		std::cout << "Error: At least 2 latitudes are required." << std::endl;
		return (-1);
	}
	if (nLongitudes < 2) {
		std::cout << "Error: At least 2 longitudes are required." << std::endl;
		return (-1);
	}

	// Generate the mesh
	Mesh mesh;

	NodeVector & nodes = mesh.nodes;
	FaceVector & faces = mesh.faces;

	// Generate nodes
	nodes.push_back(Node(0.0, 0.0, -1.0));
	for (int j = 1; j < nLatitudes; j++) {
		for (int i = 0; i < nLongitudes; i++) {
			Real dPhiFrac =
				  static_cast<Real>(j)
				/ static_cast<Real>(nLatitudes);

			Real dLambdaFrac =
				  static_cast<Real>(i)
				/ static_cast<Real>(nLongitudes);

			Real dPhi = M_PI * (dPhiFrac - 0.5);
			Real dLambda = 2.0 * M_PI * dLambdaFrac;

			Real dX = cos(dPhi) * cos(dLambda);
			Real dY = cos(dPhi) * sin(dLambda);
			Real dZ = sin(dPhi);

			nodes.push_back(Node(dX, dY, dZ));
		}
	}
	nodes.push_back(Node(0.0, 0.0, +1.0));

	// Generate south polar faces
	for (int i = 0; i < nLongitudes; i++) {
		Face face(4);
		face.SetNode(0, 0);
		face.SetNode(1, (i+1) % nLongitudes + 1);
		face.SetNode(2, i + 1);
		face.SetNode(3, 0);

#ifndef ONLY_GREAT_CIRCLES
		face.edges[1].type = Edge::Type_ConstantLatitude;
		face.edges[3].type = Edge::Type_ConstantLatitude;
#endif

		faces.push_back(face);
	}

	// Generate interior faces
	for (int j = 1; j < nLatitudes-1; j++) {
		int iThisLatNodeIx = (j-1) * nLongitudes + 1;
		int iNextLatNodeIx =  j    * nLongitudes + 1;

		for (int i = 0; i < nLongitudes; i++) {
			Face face(4);
			face.SetNode(0, iThisLatNodeIx + (i + 1) % nLongitudes);
			face.SetNode(1, iNextLatNodeIx + (i + 1) % nLongitudes);
			face.SetNode(2, iNextLatNodeIx + i);
			face.SetNode(3, iThisLatNodeIx + i);

#ifndef ONLY_GREAT_CIRCLES
			face.edges[1].type = Edge::Type_ConstantLatitude;
			face.edges[3].type = Edge::Type_ConstantLatitude;
#endif

			faces.push_back(face);
		}
	}

	// Generate north polar faces
	{
		int iThisLatNodeIx = (nLatitudes - 2) * nLongitudes + 1;
		int iNorthPolarNodeIx = static_cast<int>(nodes.size()-1);
		for (int i = 0; i < nLongitudes; i++) {
			Face face(4);
			face.SetNode(0, iNorthPolarNodeIx);
			face.SetNode(1, iThisLatNodeIx + i);
			face.SetNode(2, iThisLatNodeIx + (i + 1) % nLongitudes);
			face.SetNode(3, iNorthPolarNodeIx);

#ifndef ONLY_GREAT_CIRCLES
			face.edges[1].type = Edge::Type_ConstantLatitude;
			face.edges[3].type = Edge::Type_ConstantLatitude;
#endif

			faces.push_back(face);
		}
	}

	// Announce
	std::cout << "..Writing mesh to file [" << strOutputFile.c_str() << "] ";
	std::cout << std::endl;

	// Output the mesh
	mesh.Write(strOutputFile);

	// Announce
	std::cout << "..Mesh generator exited successfully" << std::endl;
	std::cout << "=========================================================";
	std::cout << std::endl;

} catch(Exception & e) {
	std::cout << e.ToString() << std::endl;
}
}

///////////////////////////////////////////////////////////////////////////////

