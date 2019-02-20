///////////////////////////////////////////////////////////////////////////////
///
///	\file    GenerateCSMesh.cpp
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

#include "GridElements.h"
#include "Exception.h"
#include "Announce.h"
#include "STLStringHelper.h"
#include "NetCDFUtilities.h"

#include "netcdfcpp.h"

#include <cmath>
#include <iostream>

///////////////////////////////////////////////////////////////////////////////

int InsertCSSubNode(
	int ix0,
	int ix1,
	Real alpha,
	NodeVector & nodes
) {
	Real dX = nodes[ix0].x + (nodes[ix1].x - nodes[ix0].x) * alpha;
	Real dY = nodes[ix0].y + (nodes[ix1].y - nodes[ix0].y) * alpha;
	Real dZ = nodes[ix0].z + (nodes[ix1].z - nodes[ix0].z) * alpha;

	// Project to sphere
	Real dRadius = sqrt(dX*dX + dY*dY + dZ*dZ);

	dX /= dRadius;
	dY /= dRadius;
	dZ /= dRadius;

	// Index
	int ix = nodes.size();

	// Insert node
	nodes.push_back(Node(dX, dY, dZ));

	return ix;
}

///////////////////////////////////////////////////////////////////////////////

void GenerateCSMultiEdgeVertices(
	int nRefineLevel,
	int ix0,
	int ix1,
	NodeVector & nodes,
	MultiEdge & edge
) {
	edge.clear();
	edge.push_back(ix0);

	int i;
	for (i = 1; i < nRefineLevel; i++) {

		// Nodes along line in Cartesian geometry
		Real alpha =
			static_cast<Real>(i) / static_cast<Real>(nRefineLevel);

		alpha = 0.5 * (tan(0.25 * M_PI * (2.0 * alpha - 1.0)) + 1.0);

		// Insert node along edge
		int ixNode = InsertCSSubNode(ix0, ix1, alpha, nodes);

		// Add node to edge
		edge.push_back(ixNode);
	}

	edge.push_back(ix1);
}

///////////////////////////////////////////////////////////////////////////////

void GenerateFacesFromQuad(
	int nResolution,
	int iPanel,
	const MultiEdge & edge0,
	const MultiEdge & edge1,
	const MultiEdge & edge2,
	const MultiEdge & edge3,
	NodeVector & nodes,
	FaceVector & vecFaces
) {
	MultiEdge edgeTop;
	MultiEdge edgeBot = edge0;

	for (int j = 0; j < nResolution; j++) {

		// Generate top level edge
		if (j != nResolution-1) {
			int ix0 = edge1[j+1];
			int ix1 = edge2[j+1];

			GenerateCSMultiEdgeVertices(nResolution, ix0, ix1, nodes, edgeTop);

		} else {
			edgeTop = edge3;
		}

		// Generate face
		for (int i = 0; i < nResolution; i++) {
			Face face(4);
			face.SetNode(0, edgeBot[i+1]);
			face.SetNode(1, edgeTop[i+1]);
			face.SetNode(2, edgeTop[i]);
			face.SetNode(3, edgeBot[i]);

			vecFaces.push_back(face);
		}

		// Increment row
		edgeBot = edgeTop;
	}
}


///////////////////////////////////////////////////////////////////////////////
// 
// Input Parameters:
// Number of elements in mesh: int nResolution;
// Alternate arrangement: bool fAlt;
// Output filename:  std::string strOutputFile;
// 
// Output Parameters: Mesh*
// 
extern "C" 
int GenerateCSMesh(
	Mesh & mesh,
	int nResolution,
	bool fAlt,
	std::string strOutputFile,
	std::string strOutputFormat
) {

	NcError error(NcError::silent_nonfatal);

try {

    // Check command line parameters (data type arguments)
    STLStringHelper::ToLower(strOutputFormat);

	NcFile::FileFormat eOutputFormat =
		GetNcFileFormatFromString(strOutputFormat);
	if (eOutputFormat == NcFile::BadFormat) {
		_EXCEPTION1("Invalid \"out_format\" value (%s), "
			"expected [Classic|Offset64Bits|Netcdf4|Netcdf4Classic]",
			strOutputFormat.c_str());
	}

	// Announce
	std::cout << "=========================================================";
	std::cout << std::endl;
	std::cout << "..Generating mesh with resolution [" << nResolution << "]";
	std::cout << std::endl;

	NodeVector & nodes = mesh.nodes;
	FaceVector & faces = mesh.faces;
    mesh.type = Mesh::MeshType_CubedSphere;

	// Generate corner points
	Real dInvDeltaX = 1.0 / sqrt(3.0);

	nodes.push_back(Node(+dInvDeltaX, -dInvDeltaX, -dInvDeltaX));
	nodes.push_back(Node(+dInvDeltaX, +dInvDeltaX, -dInvDeltaX));
	nodes.push_back(Node(-dInvDeltaX, +dInvDeltaX, -dInvDeltaX));
	nodes.push_back(Node(-dInvDeltaX, -dInvDeltaX, -dInvDeltaX));
	nodes.push_back(Node(+dInvDeltaX, -dInvDeltaX, +dInvDeltaX));
	nodes.push_back(Node(+dInvDeltaX, +dInvDeltaX, +dInvDeltaX));
	nodes.push_back(Node(-dInvDeltaX, +dInvDeltaX, +dInvDeltaX));
	nodes.push_back(Node(-dInvDeltaX, -dInvDeltaX, +dInvDeltaX));

	// Generate edges
	MultiEdgeVector vecMultiEdges;
	vecMultiEdges.resize(12);

	GenerateCSMultiEdgeVertices(nResolution, 0, 1, nodes, vecMultiEdges[0]);
	GenerateCSMultiEdgeVertices(nResolution, 1, 2, nodes, vecMultiEdges[1]);
	GenerateCSMultiEdgeVertices(nResolution, 2, 3, nodes, vecMultiEdges[2]);
	GenerateCSMultiEdgeVertices(nResolution, 3, 0, nodes, vecMultiEdges[3]);

	GenerateCSMultiEdgeVertices(nResolution, 0, 4, nodes, vecMultiEdges[4]);
	GenerateCSMultiEdgeVertices(nResolution, 1, 5, nodes, vecMultiEdges[5]);
	GenerateCSMultiEdgeVertices(nResolution, 2, 6, nodes, vecMultiEdges[6]);
	GenerateCSMultiEdgeVertices(nResolution, 3, 7, nodes, vecMultiEdges[7]);

	GenerateCSMultiEdgeVertices(nResolution, 4, 5, nodes, vecMultiEdges[8]);
	GenerateCSMultiEdgeVertices(nResolution, 5, 6, nodes, vecMultiEdges[9]);
	GenerateCSMultiEdgeVertices(nResolution, 6, 7, nodes, vecMultiEdges[10]);
	GenerateCSMultiEdgeVertices(nResolution, 7, 4, nodes, vecMultiEdges[11]);

	// Generate equatorial faces
	GenerateFacesFromQuad(
		nResolution,
		0,
		vecMultiEdges[0],
		vecMultiEdges[4],
		vecMultiEdges[5],
		vecMultiEdges[8],
		nodes,
		faces);

	GenerateFacesFromQuad(
		nResolution,
		1,
		vecMultiEdges[1],
		vecMultiEdges[5],
		vecMultiEdges[6],
		vecMultiEdges[9],
		nodes,
		faces);

	GenerateFacesFromQuad(
		nResolution,
		2,
		vecMultiEdges[2],
		vecMultiEdges[6],
		vecMultiEdges[7],
		vecMultiEdges[10],
		nodes,
		faces);

	GenerateFacesFromQuad(
		nResolution,
		3,
		vecMultiEdges[3],
		vecMultiEdges[7],
		vecMultiEdges[4],
		vecMultiEdges[11],
		nodes,
		faces);

	// Generate south polar face
	GenerateFacesFromQuad(
		nResolution,
		5,
		vecMultiEdges[2].Flip(),
		vecMultiEdges[3],
		vecMultiEdges[1].Flip(),
		vecMultiEdges[0],
		nodes,
		faces);

	// Generate north polar face
	GenerateFacesFromQuad(
		nResolution,
		4,
		vecMultiEdges[8],
		vecMultiEdges[11].Flip(),
		vecMultiEdges[9],
		vecMultiEdges[10].Flip(),
		nodes,
		faces);

	// Alternative arrangement of nodes on Faces
	if (fAlt) {
		for (int i = 0; i < faces.size(); i++) {
			int ix[4];
			for (int j = 0; j < 4; j++) {
				ix[j] = faces[i][j];
			}
			for (int j = 0; j < 4; j++) {
				faces[i].SetNode((j+1)%4, ix[j]);
			}
		}
	}

	// Output the mesh
	if (strOutputFile.size()) {

		// Announce
		std::cout << "..Writing mesh to file [" << strOutputFile.c_str() << "] ";
		std::cout << std::endl;

		mesh.Write(strOutputFile, eOutputFormat);
	}

	// Announce
	std::cout << "..Mesh generator exited successfully" << std::endl;
	std::cout << "=========================================================";
	std::cout << std::endl;

} catch(Exception & e) {
	Announce(e.ToString().c_str());
	return (0);

} catch(...) {
	return (0);
}
    return 0;
}

///////////////////////////////////////////////////////////////////////////////
