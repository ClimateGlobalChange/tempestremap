///////////////////////////////////////////////////////////////////////////////
///
///	\file    OverlapMesh.cpp
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

#include "OverlapMesh.h"

#include <iostream>

///////////////////////////////////////////////////////////////////////////////

void GenerateOverlapMesh(
	const Mesh & meshFirst,
	const Mesh & meshSecond,
	Mesh & meshOverlap
) {
	const NodeVector & nodevecFirst = meshFirst.nodes;
	const NodeVector & nodevecSecond = meshSecond.nodes;

	// Loop through all faces of the first mesh
	for (int k = 0; k < meshFirst.faces.size(); k++) {

		//std::vector<int> ixFaceQueue;
		const Node & nodeBegin = nodevecFirst[meshFirst.faces[k][0]];
/*
		const Node & node0 = nodevecFirst[meshFirst.faces[k][0]];
		const Node & node1 = nodevecFirst[meshFirst.faces[k][1]];
		const Node & node2 = nodevecFirst[meshFirst.faces[k][2]];
		const Node & node3 = nodevecFirst[meshFirst.faces[k][3]];

		Node newnode(
			0.25 * (node0.x + node1.x + node2.x + node3.x),
			0.25 * (node0.y + node1.y + node2.y + node3.y),
			0.25 * (node0.z + node1.z + node2.z + node3.z));
*/
		//printf("%1.10e %1.10e %1.10e\n", newnode.x, newnode.y, newnode.z);

		// Find the face(s) that contains this node
		// Note: This algorithm can likely be dramatically improved
		std::set<int> setStartingFaces;

		for (int l = 0; l < meshSecond.faces.size(); l++) {
			Face::NodeLocation loc;
			int ixLocation;

			meshSecond.faces[l].ContainsNode(
				nodevecSecond,
				nodeBegin,
				loc,
				ixLocation);

			if (loc == Face::NodeLocation_Interior) {
				setStartingFaces.insert(l);
				break;
			}
			if (loc == Face::NodeLocation_Edge) {
				setStartingFaces.insert(l);
			}
			if (loc == Face::NodeLocation_Corner) {
				setStartingFaces.insert(l);
			}
		}

		if (setStartingFaces.size() == 0) {
			_EXCEPTIONT("Cannot find starting face");
		}

		// Find the actual starting face
		if (setStartingFaces.size() > 1) {

			// Do something special
		}

		// Current face on second mesh
		int ixSecondFace = *(setStartingFaces.begin());

		// Trace along edge
		for (;;) {
			// Equation for edge: For a great circle arc

			// Find all intersections between this edge and the second face
			const Face & faceCurrent = meshSecond.faces[ixSecondFace];
			for (int i = 0; i < faceCurrent.edges.size(); i++) {

			}
			break;
		}

		// Look at neighboring faces and check for overlap
	}
}

///////////////////////////////////////////////////////////////////////////////

