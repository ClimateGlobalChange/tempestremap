///////////////////////////////////////////////////////////////////////////////
///
///	\file    GenerateICOMesh.cpp
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

#include <cmath>
#include <iostream>

///////////////////////////////////////////////////////////////////////////////

struct LonLatNode {
	double lon;
	double lat;

	///	<summary>
	///		Constructor
	///	</summary>
	LonLatNode(
		double _lon,
		double _lat
	) :
		lon(_lon),
		lat(_lat)
	{ }
};

typedef std::vector<LonLatNode> LonLatNodeVector;

///////////////////////////////////////////////////////////////////////////////

int InsertSubNode(
	int ix0,
	int ix1,
	double alpha,
	NodeVector & vecNodes
) {
	double dDeltaX = (vecNodes[ix1].x - vecNodes[ix0].x);
	double dDeltaY = (vecNodes[ix1].y - vecNodes[ix0].y);
	double dDeltaZ = (vecNodes[ix1].z - vecNodes[ix0].z);
	double dCartLength =
		sqrt(dDeltaX*dDeltaX + dDeltaY*dDeltaY + dDeltaZ*dDeltaZ);

	double dGamma = acos(0.5 * dCartLength);
	double dTheta = acos(1.0 - 0.5 * dCartLength * dCartLength);
	double dAlphaTheta = alpha * dTheta;
	double dBeta = M_PI - dGamma - dAlphaTheta;

	alpha = sin(dAlphaTheta) / sin(dBeta) / dCartLength;

	double dX = vecNodes[ix0].x + (vecNodes[ix1].x - vecNodes[ix0].x) * alpha;
	double dY = vecNodes[ix0].y + (vecNodes[ix1].y - vecNodes[ix0].y) * alpha;
	double dZ = vecNodes[ix0].z + (vecNodes[ix1].z - vecNodes[ix0].z) * alpha;

	// Project to sphere
	double dRadius = sqrt(dX*dX + dY*dY + dZ*dZ);

	dX /= dRadius;
	dY /= dRadius;
	dZ /= dRadius;

	// Index
	int ix = vecNodes.size();

	// Insert node
	vecNodes.push_back(Node(dX, dY, dZ));

	return ix;
}

///////////////////////////////////////////////////////////////////////////////

int InsertTriFaceCentroidNode(
	int ix0,
	int ix1,
	int ix2,
	NodeVector & vecNodes
) {
	double dX = (vecNodes[ix0].x + vecNodes[ix1].x + vecNodes[ix2].x) / 3.0;
	double dY = (vecNodes[ix0].y + vecNodes[ix1].y + vecNodes[ix2].y) / 3.0;
	double dZ = (vecNodes[ix0].z + vecNodes[ix1].z + vecNodes[ix2].z) / 3.0;

	// Project to sphere
	double dRadius = sqrt(dX*dX + dY*dY + dZ*dZ);

	dX /= dRadius;
	dY /= dRadius;
	dZ /= dRadius;

	// Index
	int ix = vecNodes.size();

	// Insert node
	vecNodes.push_back(Node(dX, dY, dZ));

	return ix;
}

///////////////////////////////////////////////////////////////////////////////

void GenerateEdgeVertices(
	int nRefineLevel,
	int ix0,
	int ix1,
	NodeVector & vecNodes,
	MultiEdge & edge
) {
	edge.clear();
	edge.push_back(ix0);

	for (int i = 1; i < nRefineLevel; i++) {

		// Nodes along line in Cartesian geometry
		double alpha =
			static_cast<double>(i) / static_cast<double>(nRefineLevel);

		// Insert node along edge
		int ixNode = InsertSubNode(ix0, ix1, alpha, vecNodes);

		// Add node to edge
		edge.push_back(ixNode);
	}

	edge.push_back(ix1);
}

///////////////////////////////////////////////////////////////////////////////

void GenerateFacesFromTriangle(
	int nRefineLevel,
	const MultiEdge & edge0,
	const MultiEdge & edge1,
	const MultiEdge & edge2,
	NodeVector & vecNodes,
	FaceVector & vecFaces
) {
	int i;
	int j;
	int k;

	int ixEndNode;

	int ixInt;

	// MultiEdges
	MultiEdge edgeBot;
	MultiEdge edgeMid;
	MultiEdge edgeTop;

	// Initial bottom edge
	edgeBot.push_back(edge0[0]);

	// Loop over all refined faces
	for (j = 0; j < nRefineLevel; j++) {

		// Generate top level vertices
		if (j == nRefineLevel-1) {
			edgeTop = edge2;
		} else {
			GenerateEdgeVertices(
				j+1, edge0[j+1], edge1[j+1], vecNodes, edgeTop);
		}

		// Generate faces
		for (i = 0; i < 2*j+1; i++) {
			// Downward pointing faces
			if (i % 2 == 0) {
				int ix = i/2;

				Face face(3);
				face.SetNode(0, edgeBot[ix]);
				face.SetNode(1, edgeTop[ix]);
				face.SetNode(2, edgeTop[ix+1]);

				vecFaces.push_back(face);

			// Upward pointing faces
			} else {
				int ix = (i-1)/2;

				Face face(3);
				face.SetNode(0, edgeTop[ix+1]);
				face.SetNode(1, edgeBot[ix+1]);
				face.SetNode(2, edgeBot[ix]);

				vecFaces.push_back(face);
			}
		}

		// New bottom edge
		edgeBot = edgeTop;
	}
}

///////////////////////////////////////////////////////////////////////////////

void ConvertFromLonLatToCartesian(
	const LonLatNodeVector & vecLonLatNodes,
	NodeVector & vecNodes
) {
	vecNodes.resize(vecLonLatNodes.size());

	// Loop over all nodes
	int i;
	for (i = 0; i < vecLonLatNodes.size(); i++) {
		vecNodes[i].x =
			sin(vecLonLatNodes[i].lon) * cos(vecLonLatNodes[i].lat);
		vecNodes[i].y =
			cos(vecLonLatNodes[i].lon) * cos(vecLonLatNodes[i].lat);
		vecNodes[i].z =
			sin(vecLonLatNodes[i].lat);
	}
}

///////////////////////////////////////////////////////////////////////////////

void GenerateIcosahedralQuadGrid(
	int nRefineLevel,
	NodeVector & vecNodes,
	FaceVector & vecFaces
) {
	// Latitude of nodes (Northern Hemisphere)
	const double NodeLat = atan(0.5);

	// Store all icosahedral nodes
	LonLatNodeVector vecLonLatNodes;

	vecLonLatNodes.push_back(LonLatNode(0.0,          -0.5*M_PI));
	vecLonLatNodes.push_back(LonLatNode(0.0,          -NodeLat));
	vecLonLatNodes.push_back(LonLatNode(2.0*M_PI*0.2, -NodeLat));
	vecLonLatNodes.push_back(LonLatNode(2.0*M_PI*0.4, -NodeLat));
	vecLonLatNodes.push_back(LonLatNode(2.0*M_PI*0.6, -NodeLat));
	vecLonLatNodes.push_back(LonLatNode(2.0*M_PI*0.8, -NodeLat));
	vecLonLatNodes.push_back(LonLatNode(2.0*M_PI*0.1, +NodeLat));
	vecLonLatNodes.push_back(LonLatNode(2.0*M_PI*0.3, +NodeLat));
	vecLonLatNodes.push_back(LonLatNode(2.0*M_PI*0.5, +NodeLat));
	vecLonLatNodes.push_back(LonLatNode(2.0*M_PI*0.7, +NodeLat));
	vecLonLatNodes.push_back(LonLatNode(2.0*M_PI*0.9, +NodeLat));
	vecLonLatNodes.push_back(LonLatNode(0.0,          +0.5*M_PI));

	// Convert icosahedral nodes to Cartesian geometry
	ConvertFromLonLatToCartesian(vecLonLatNodes, vecNodes);

	// Vector of edges
	MultiEdgeVector vecEdges;
	vecEdges.resize(30);

	// Generate vertices along edges
	for (int i = 0; i < 5; i++) {
		GenerateEdgeVertices(
			nRefineLevel, 0, i+1, vecNodes, vecEdges[i]);
	}

	for (int i = 0; i < 5; i++) {
		GenerateEdgeVertices(
			nRefineLevel, i+1, ((i+1)%5)+1, vecNodes, vecEdges[i+5]);
	}

	GenerateEdgeVertices(nRefineLevel, 1, 6, vecNodes, vecEdges[10]);
	GenerateEdgeVertices(nRefineLevel, 6, 2, vecNodes, vecEdges[11]);
	GenerateEdgeVertices(nRefineLevel, 2, 7, vecNodes, vecEdges[12]);
	GenerateEdgeVertices(nRefineLevel, 7, 3, vecNodes, vecEdges[13]);
	GenerateEdgeVertices(nRefineLevel, 3, 8, vecNodes, vecEdges[14]);
	GenerateEdgeVertices(nRefineLevel, 8, 4, vecNodes, vecEdges[15]);
	GenerateEdgeVertices(nRefineLevel, 4, 9, vecNodes, vecEdges[16]);
	GenerateEdgeVertices(nRefineLevel, 9, 5, vecNodes, vecEdges[17]);
	GenerateEdgeVertices(nRefineLevel, 5, 10, vecNodes, vecEdges[18]);
	GenerateEdgeVertices(nRefineLevel, 10, 1, vecNodes, vecEdges[19]);

	for (int i = 0; i < 5; i++) {
		GenerateEdgeVertices(
			nRefineLevel, i+6, ((i+1)%5)+6, vecNodes, vecEdges[i+20]);
	}

	for (int i = 0; i < 5; i++) {
		GenerateEdgeVertices(
			nRefineLevel, i+6, 11, vecNodes, vecEdges[i+25]);
	}

	// Generate south polar faces
	for (int i = 0; i < 5; i++) {
		GenerateFacesFromTriangle(
			nRefineLevel,
			vecEdges[i],
			vecEdges[(i+1)%5],
			vecEdges[i+5],
			vecNodes,
			vecFaces
		);
	}

	// Generate south equatorial faces
	for (int i = 0; i < 5; i++) {
		GenerateFacesFromTriangle(
			nRefineLevel,
			vecEdges[2*i+10],
			vecEdges[i+5],
			vecEdges[2*i+11],
			vecNodes,
			vecFaces
		);
	}

	// Generate north equatorial faces
	for (int i = 0; i < 5; i++) {
		GenerateFacesFromTriangle(
			nRefineLevel,
			vecEdges[i+20],
			vecEdges[2*i+11],
			vecEdges[2*((i+1)%5)+10].Flip(),
			vecNodes,
			vecFaces
		);
	}

	// Generate north polar faces
	for (int i = 0; i < 5; i++) {
		GenerateFacesFromTriangle(
			nRefineLevel,
			vecEdges[i+25],
			vecEdges[i+20],
			vecEdges[((i+1)%5)+25].Flip(),
			vecNodes,
			vecFaces
		);
	}
}

///////////////////////////////////////////////////////////////////////////////

void Dual(
	Mesh & mesh
) {
	const int EdgeCountHexagon = 6;

	// Generate ReverseNodeArray
	mesh.ConstructReverseNodeArray();

	// Backup Nodes and Faces
	NodeVector nodesOld = mesh.nodes;
	FaceVector facesOld = mesh.faces;

	mesh.nodes.clear();
	mesh.faces.clear();

	// Generate new Node array
	for (int i = 0; i < facesOld.size(); i++) {
		Node node;
		for (int j = 0; j < facesOld[i].edges.size(); j++) {
			node.x += nodesOld[facesOld[i][j]].x;
			node.y += nodesOld[facesOld[i][j]].y;
			node.z += nodesOld[facesOld[i][j]].z;
		}

		node.x /= static_cast<double>(facesOld[i].edges.size());
		node.y /= static_cast<double>(facesOld[i].edges.size());
		node.z /= static_cast<double>(facesOld[i].edges.size());

		double dMag = node.Magnitude();

		node.x /= dMag;
		node.y /= dMag;
		node.z /= dMag;

		mesh.nodes.push_back(node);
	}

	// Generate new Face array
	for (int i = 0; i < nodesOld.size(); i++) {
		const int nEdges = mesh.revnodearray[i].size();

		Face face(EdgeCountHexagon);
		Face faceTemp(EdgeCountHexagon);

		int ixNode = 0;
		std::set<int>::const_iterator iter = mesh.revnodearray[i].begin();
		for (; iter != mesh.revnodearray[i].end(); iter++) {
			faceTemp.SetNode(ixNode, *iter);
			ixNode++;
		}

		// Reorient Faces
		Node nodeCentral = nodesOld[i];
		Node node0 = mesh.nodes[faceTemp[0]] - nodeCentral;

		Node nodeCross = CrossProduct(mesh.nodes[faceTemp[0]], nodeCentral);
		
		double dNode0Mag = node0.Magnitude();

		// Determine the angles about the central Node of each Face Node
		std::vector<double> dAngles;
		dAngles.resize(faceTemp.edges.size());
		dAngles[0] = 0.0;

		for (int j = 1; j < nEdges; j++) {
			Node nodeDiff = mesh.nodes[faceTemp[j]] - nodeCentral;
			double dNodeDiffMag = nodeDiff.Magnitude();

			double dSide = DotProduct(nodeCross, nodeDiff);

			double dDotNorm =
				DotProduct(node0, nodeDiff) / (dNode0Mag * dNodeDiffMag);

			double dAngle;
			if (dDotNorm > 1.0) {
				dDotNorm = 1.0;
			}

			dAngles[j] = acos(dDotNorm);

			if (dSide > 0.0) {
				dAngles[j] = - dAngles[j] + 2.0 * M_PI;
			}
		}

		// Orient each Face by putting Nodes in order of increasing angle
		double dCurrentAngle = 0.0;
		face.SetNode(0, faceTemp[0]);
		for (int j = 1; j < nEdges; j++) {
			int ixNextNode = 1;
			double dNextAngle = 2.0 * M_PI;

			for (int k = 1; k < nEdges; k++) {
				if ((dAngles[k] > dCurrentAngle) && (dAngles[k] < dNextAngle)) {
					ixNextNode = k;
					dNextAngle = dAngles[k];
				}
			}

			face.SetNode(j, faceTemp[ixNextNode]);
			dCurrentAngle = dNextAngle;
		}

		// Fill in missing edges
		for (int j = nEdges; j < EdgeCountHexagon; j++) {
			face.SetNode(j, face[nEdges-1]);
		}

		mesh.faces.push_back(face);
	}
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

try {

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
		CommandLineString(strOutputFile, "out", "outICOMesh.g");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	// Generate Mesh
	AnnounceBanner();
	AnnounceStartBlock("Generating Mesh");
	Mesh mesh;
	GenerateIcosahedralQuadGrid(nResolution, mesh.nodes, mesh.faces);
	AnnounceEndBlock("Done");

	// Generate the dual grid
	if (fDual) {
		Dual(mesh);
	}

	// Output the mesh
	AnnounceStartBlock("Writing Mesh to file");
	Announce("Mesh size: Nodes [%i] Elements [%i]",
		mesh.nodes.size(), mesh.faces.size());
	mesh.Write(strOutputFile);
	AnnounceEndBlock("Done");

} catch(Exception & e) {
	std::cout << e.ToString() << std::endl;
}
}

///////////////////////////////////////////////////////////////////////////////


