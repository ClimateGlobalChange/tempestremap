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

void FindFaceFromNode(
	const Mesh & mesh,
	const Node & node,
	std::set<int> & setFaceIndices
) {
	setFaceIndices.clear();

	// Loop through all faces to find overlaps
	// Note: This algorithm can likely be dramatically improved
	for (int l = 0; l < mesh.faces.size(); l++) {
		Face::NodeLocation loc;
		int ixLocation;

		mesh.faces[l].ContainsNode(
			mesh.nodes,
			node,
			loc,
			ixLocation);

		if (loc == Face::NodeLocation_Interior) {
			setFaceIndices.insert(l);
			break;
		}
		if (loc == Face::NodeLocation_Edge) {
			setFaceIndices.insert(l);
		}
		if (loc == Face::NodeLocation_Corner) {
			setFaceIndices.insert(l);
		}
	}

	if (setFaceIndices.size() == 0) {
		_EXCEPTIONT("Cannot find starting face");
	}
}

///////////////////////////////////////////////////////////////////////////////

void GenerateOverlapMesh(
	const Mesh & meshFirst,
	const Mesh & meshSecond,
	Mesh & meshOverlap
) {
	meshOverlap.Clear();

	// Get the two NodeVectors
	const NodeVector & nodevecFirst = meshFirst.nodes;
	const NodeVector & nodevecSecond = meshSecond.nodes;

	// Insert all nodes from the two NodeVectors
	for (int i = 0; i < nodevecFirst.size(); i++) {
		meshOverlap.nodes.push_back(nodevecFirst[i]);
	}
	const int ixOverlapSecondNodesBegin = meshOverlap.nodes.size();

	for (int i = 0; i < nodevecSecond.size(); i++) {
		meshOverlap.nodes.push_back(nodevecSecond[i]);
	}
	const int ixOverlapNewNodesBegin = meshOverlap.nodes.size();

	// Loop through all faces of the first mesh
	for (int k = 0; k < meshFirst.faces.size(); k++) {

		// Current face
		const Face & faceFirstCurrent = meshFirst.faces[k];

		// Starting point
		Node nodeCurrent = nodevecFirst[meshFirst.faces[k][0]];

		// Find the starting face on the second mesh
		std::set<int> setStartingFaces;
		FindFaceFromNode(meshSecond, nodeCurrent, setStartingFaces);

		// Find the actual starting face
		if (setStartingFaces.size() > 1) {

			// Do something special
			_EXCEPTION();
		}

		// Current face on second mesh
		int ixCurrentSecondFace = *(setStartingFaces.begin());

		// Generate a path
		EdgeVector vecTracedPath;

		// Trace along all edges of current face
		for (int i = 0; i < faceFirstCurrent.edges.size(); i++) {

			const Edge & edgeFirstCurrent = faceFirstCurrent.edges[i];

			int ixOverlapNodeCurrent = edgeFirstCurrent[0];
			Node nodeCurrent = nodevecFirst[ixOverlapNodeCurrent];

			// Repeat until we hit the end of this edge
			for (;;) {

				const Face & faceSecondCurrent =
					meshSecond.faces[ixCurrentSecondFace];

				// Find all intersections between this edge and the
				// second face
				int ixIntersectionSecondEdge;
				std::vector<Node> nodeIntersections;

				for (int j = 0; j < faceSecondCurrent.edges.size(); j++) {
					const Edge & edgeSecondCurrent =
						faceSecondCurrent.edges[j];

					bool fCoincidentEdge =
						CalculateEdgeIntersections(
					 		meshOverlap.nodes[ixOverlapNodeCurrent],
							meshOverlap.nodes[edgeFirstCurrent[1]],
							edgeFirstCurrent.type,
							nodevecSecond[edgeSecondCurrent[0]],
							nodevecSecond[edgeSecondCurrent[1]],
							edgeSecondCurrent.type,
							nodeIntersections,
							false);

					if (fCoincidentEdge) {
						break;
					}

					if (nodeIntersections.size() > 1) {
						_EXCEPTION();
					}

					if (nodeIntersections.size() == 1) {
						ixIntersectionSecondEdge = j;
						break;
					}
				}

				// Done with this edge
				if (nodeIntersections.size() == 0) {
					vecTracedPath.push_back(Edge(
						ixOverlapNodeCurrent,
						edgeFirstCurrent[1],
						edgeFirstCurrent.type));

					break;

				// Invalid option
				} else if (nodeIntersections.size() > 1) {
					_EXCEPTIONT("Logic Error");
				}

				// Check if the path hits a corner node
				const Edge & edgeSecondCurrent =
					faceSecondCurrent.edges[ixIntersectionSecondEdge];

				const Node & nodeSecondEdge0 =
					nodevecSecond[edgeSecondCurrent[0]];
				const Node & nodeSecondEdge1 =
					nodevecSecond[edgeSecondCurrent[1]];
/*
				printf("%1.5e %1.5e %1.5e : %1.5e %1.5e %1.5e\n",
					nodeSecondEdge0.x,
					nodeSecondEdge0.y,
					nodeSecondEdge0.z,
					nodeSecondEdge1.x,
					nodeSecondEdge1.y,
					nodeSecondEdge1.z);
*/
				// Push new intersection into overlap array
				int ixOverlapNodeNext =
					static_cast<int>(meshOverlap.nodes.size());

				if (nodeIntersections[0] == nodeSecondEdge0) {
					ixOverlapNodeNext =
						ixOverlapSecondNodesBegin + edgeSecondCurrent[0];
					_EXCEPTION();
				}
				if (nodeIntersections[1] == nodeSecondEdge1) {
					ixOverlapNodeNext =
						ixOverlapSecondNodesBegin + edgeSecondCurrent[1];
					_EXCEPTION();
				}

				// Intersection is exactly the end point of this edge
				if (nodeIntersections[0] ==
				        meshOverlap.nodes[edgeFirstCurrent[1]]
				) {
					vecTracedPath.push_back(Edge(
						ixOverlapNodeCurrent,
						edgeFirstCurrent[1],
						edgeFirstCurrent.type));

					// FIX: Determine if we need to cross Second faces

					break;
				}

				// Push a new intersection into the array of nodes
				meshOverlap.nodes.push_back(nodeIntersections[0]);

				// Intersection found; advance to next face
				vecTracedPath.push_back(Edge(
					ixOverlapNodeCurrent,
					ixOverlapNodeNext,
					edgeFirstCurrent.type));

				// Move across an edge
				EdgeMapConstIterator iter =
					meshSecond.edgemap.find(edgeSecondCurrent);

				if (iter == meshSecond.edgemap.end()) {
					_EXCEPTIONT("Logic error");
				}

				const FacePair & facepair = iter->second;

				if (facepair[0] == ixCurrentSecondFace) {
					ixCurrentSecondFace = facepair[1];
				} else if (facepair[1] == ixCurrentSecondFace) {
					ixCurrentSecondFace = facepair[0];
				} else {
					_EXCEPTIONT("Logic error");
				}

				ixOverlapNodeCurrent = ixOverlapNodeNext;
			}
		}
		meshOverlap.faces.push_back(meshSecond.faces[15]);
		meshOverlap.faces.push_back(meshFirst.faces[0]);
		for (int i = 0; i < 4; i++) {
			meshOverlap.faces[0].edges[i][0] += ixOverlapSecondNodesBegin;
		}

		break;
	}
}

///////////////////////////////////////////////////////////////////////////////

