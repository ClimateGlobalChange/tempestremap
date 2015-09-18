///////////////////////////////////////////////////////////////////////////////
///
///	\file    GenerateVolumetricMesh.cpp
///	\author  Paul Ullrich
///	\version September 18, 2015
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
#include "FiniteElementTools.h"
#include "GaussLobattoQuadrature.h"
#include "Exception.h"

#include "netcdfcpp.h"

#include <cmath>
#include <iostream>

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	NcError error(NcError::silent_nonfatal);

try {
	// Input filename
	std::string strInputFile;

	// Output mesh filename
	std::string strOutputMesh;

	// Output connectivity filename
	std::string strOutputConnectivity;

	// Number of elements in mesh
	int nP;

	// Merge faces
	bool fCGLL = true;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in", "");
		CommandLineString(strOutputMesh, "out_mesh", "");
		CommandLineString(strOutputConnectivity, "out_connect", "");
		CommandLineInt(nP, "np", 2);
		//CommandLineBool(fCGLL, "cgll");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	// Check file names
	if (strInputFile == "") {
		std::cout << "ERROR: No input file specified" << std::endl;
		return (-1);
	}
	if (nP < 1) {
		std::cout << "ERROR: --np must be >= 2" << std::endl;
		return (-1);
	}

	std::cout << "=========================================================";

	// Load input mesh
	std::cout << std::endl;
	std::cout << "..Loading input mesh" << std::endl;

	Mesh meshIn(strInputFile);
	meshIn.RemoveZeroEdges();

	// Number of elements
	int nElements = meshIn.faces.size();

	// Gauss-Lobatto quadrature nodes and weights
	std::cout << "..Computing sub-volume boundaries" << std::endl;
	DataVector<double> dG;
	DataVector<double> dW;

	GaussLobattoQuadrature::GetPoints(nP, 0.0, 1.0, dG, dW);

	// Accumulated weight vector
	DataVector<double> dAccumW(nP+1);
	dAccumW[0] = 0.0;
	for (int i = 1; i < nP+1; i++) {
		dAccumW[i] = dAccumW[i-1] + dW[i-1];
	}
	if (fabs(dAccumW[dAccumW.GetRows()-1] - 1.0) > 1.0e-14) {
		_EXCEPTIONT("Logic error in accumulated weight");
	}

	// Merge face map
	DataMatrix3D<int> dataGLLnodes(nP, nP, nElements);
	std::vector<Node> vecNodes;
	std::map<Node, int> mapFaces;

	// Generate new mesh
	std::cout << "..Generating sub-volumes" << std::endl;
	Mesh meshOut;

	for (size_t f = 0; f < nElements; f++) {

		const Face & face = meshIn.faces[f];

		if (face.edges.size() != 4) {
			_EXCEPTIONT("Input mesh must only contain quadrilaterals");
		}

		const Node & node0 = meshIn.nodes[face[0]];
		const Node & node1 = meshIn.nodes[face[1]];
		const Node & node2 = meshIn.nodes[face[2]];
		const Node & node3 = meshIn.nodes[face[3]];

		for (int q = 0; q < nP; q++) {
		for (int p = 0; p < nP; p++) {

			// Build unique node array if CGLL
			if (fCGLL) {

				// Get local nodal location
				Node nodeGLL;
				Node dDx1G;
				Node dDx2G;

				ApplyLocalMap(
					face,
					meshIn.nodes,
					dG[p],
					dG[q],
					nodeGLL,
					dDx1G,
					dDx2G);

				// Determine if this is a unique Node
				std::map<Node, int>::const_iterator iter =
					mapFaces.find(nodeGLL);

				if (iter == mapFaces.end()) {

					// Insert new unique node into map
					int ixNode = static_cast<int>(mapFaces.size());
					mapFaces.insert(std::pair<Node, int>(nodeGLL, ixNode));
					dataGLLnodes[q][p][f] = ixNode + 1;
					vecNodes.push_back(nodeGLL);

				} else {
					dataGLLnodes[q][p][f] = iter->second + 1;
				}

			// Non-unique node array if DGLL
			} else {
				dataGLLnodes[q][p][f] = nP * nP * f + q * nP + p;
			}

			// Get volumetric region
			Node nodeOut0 =
				InterpolateQuadrilateralNode(
					node0, node1, node2, node3,
					dAccumW[p], dAccumW[q]);

			Node nodeOut1 =
				InterpolateQuadrilateralNode(
					node0, node1, node2, node3,
					dAccumW[p+1], dAccumW[q]);

			Node nodeOut2 =
				InterpolateQuadrilateralNode(
					node0, node1, node2, node3,
					dAccumW[p+1], dAccumW[q+1]);

			Node nodeOut3 =
				InterpolateQuadrilateralNode(
					node0, node1, node2, node3,
					dAccumW[p], dAccumW[q+1]);

			int nNodeStart = meshOut.nodes.size();
			meshOut.nodes.push_back(nodeOut0);
			meshOut.nodes.push_back(nodeOut1);
			meshOut.nodes.push_back(nodeOut2);
			meshOut.nodes.push_back(nodeOut3);

			Face faceNew(4);
			faceNew.SetNode(0, nNodeStart);
			faceNew.SetNode(1, nNodeStart+1);
			faceNew.SetNode(2, nNodeStart+2);
			faceNew.SetNode(3, nNodeStart+3);

			meshOut.faces.push_back(faceNew);
		}
		}
	}

	// Build connectivity and write to file
	if (strOutputConnectivity != "") {

		std::cout << "..Constructing connectivity file" << std::endl;

		std::vector< std::set<int> > vecConnectivity;
		vecConnectivity.resize(mapFaces.size());

		for (size_t f = 0; f < nElements; f++) {

			for (int q = 0; q < nP; q++) {
			for (int p = 0; p < nP; p++) {

				std::set<int> & setLocalConnectivity =
					vecConnectivity[dataGLLnodes[q][p][f]-1];

				// Connect in all directions
				if (p != 0) {
					setLocalConnectivity.insert(
						dataGLLnodes[q][p-1][f]);
				}
				if (p != (nP-1)) {
					setLocalConnectivity.insert(
						dataGLLnodes[q][p+1][f]);
				}
				if (q != 0) {
					setLocalConnectivity.insert(
						dataGLLnodes[q-1][p][f]);
				}
				if (q != (nP-1)) {
					setLocalConnectivity.insert(
						dataGLLnodes[q+1][p][f]);
				}
			}
			}
		}

		// Open output file
		FILE * fp = fopen(strOutputConnectivity.c_str(), "w");
		for (size_t f = 0; f < vecConnectivity.size(); f++) {
			const Node & node = vecNodes[f];

			double dLon = atan2(node.y, node.x);
			double dLat = asin(node.z);

			if (dLon < 0.0) {
				dLon += 2.0 * M_PI;
			}

			fprintf(fp, "%1.14f,", dLon / M_PI * 180.0);
			fprintf(fp, "%1.14f", dLat / M_PI * 180.0);

			std::set<int>::const_iterator iter = vecConnectivity[f].begin();
			for (; iter != vecConnectivity[f].end(); iter++) {
				fprintf(fp, ",%i", *iter);
			}
			if (f != vecConnectivity.size()-1) {
				fprintf(fp,"\n");
			}
		}
		fclose(fp);
	}

	// Equalize coincident nodes
	std::cout << "..Equalizing coincident nodes" << std::endl;
	EqualizeCoincidentNodes(meshOut);

	// Write the mesh
	if (strOutputMesh != "") {
		std::cout << "..Writing mesh" << std::endl;
		meshOut.Write(strOutputMesh);
	}

	// Announce
	std::cout << "..Mesh generator exited successfully" << std::endl;
	std::cout << "=========================================================";
	std::cout << std::endl;

} catch(Exception & e) {
	std::cout << e.ToString() << std::endl;
}
}

///////////////////////////////////////////////////////////////////////////////


