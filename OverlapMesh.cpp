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

///	<summary>
///		A segment connecting two nodes that also has an associated face
///		on the First and Second mesh.
///	</summary>
class PathSegment : public Edge {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	PathSegment(
		int a_node0,
		int a_node1,
		Edge::Type a_type,
		int a_ixFirstFace,
		int a_ixSecondFace
	) :
		Edge(a_node0, a_node1, a_type),
		ixFirstFace(a_ixFirstFace),
		ixSecondFace(a_ixSecondFace)
	{ }

public:
	///	<summary>
	///		Origin face on first mesh.
	///	</summary>
	int ixFirstFace;

	///	<summary>
	///		Origin face on second mesh.
	///	</summary>
	int ixSecondFace;
};

///	<summary>
///		A vector of PathSegments.
///	</summary>
typedef std::vector<PathSegment> PathSegmentVector;

///////////////////////////////////////////////////////////////////////////////

void FindFaceFromNode(
	const Mesh & mesh,
	const Node & node,
	std::vector<int> & vecFaceIndices
) {
	vecFaceIndices.clear();

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
			vecFaceIndices.push_back(l);
			break;
		}
		if (loc == Face::NodeLocation_Edge) {
			vecFaceIndices.push_back(l);
		}
		if (loc == Face::NodeLocation_Corner) {
			vecFaceIndices.push_back(l);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void FindFaceNearNode(
	const Mesh & mesh,
	const Node & nodeBegin,
	const Node & nodeEnd,
	const Edge::Type edgetype,
	const std::vector<int> & vecPossibleFaces,
	std::vector<int> & vecFaceIndices
) {
	// Clear vector of nearby faces
	vecFaceIndices.clear();

	// Find the nudged node towards nodeEnd with type edgetype
	Node nodeNudged;
	NudgeAlongEdge(nodeBegin, nodeEnd, edgetype, nodeNudged);

	// Loop through all faces
	for (int i = 0; i < vecPossibleFaces.size(); i++) {

		Face::NodeLocation loc;
		int ixLocation;

		int ixPossibleFace = vecPossibleFaces[i];

		mesh.faces[ixPossibleFace].ContainsNode(
			mesh.nodes,
			nodeNudged,
			loc,
			ixLocation);

		if (loc == Face::NodeLocation_Exterior) {
			continue;
		} else {
			vecFaceIndices.push_back(ixPossibleFace);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void FindFaceNearNode(
	const Mesh & mesh,
	int ixNode,
	const Node & nodeEnd,
	const Edge::Type edgetype,
	std::vector<int> & vecFaceIndices
) {
	const Node & nodeBegin = mesh.nodes[ixNode];

	// Clear vector of nearby faces
	vecFaceIndices.clear();

	// Find the nudged node towards nodeEnd with type edgetype
	Node nodeNudged;
	NudgeAlongEdge(nodeBegin, nodeEnd, edgetype, nodeNudged);

	//printf("%1.14e %1.14e %1.14e\n", nodeBegin.x, nodeBegin.y, nodeBegin.z);
	//printf("%1.14e %1.14e %1.14e\n", nodeEnd.x, nodeEnd.y, nodeEnd.z);
	//printf("%1.14e %1.14e %1.14e\n", nodeNudged.x, nodeNudged.y, nodeNudged.z);

	// Find this nudged node in neighboring faces
	const std::set<int> & vecNearbyFaces = mesh.revnodearray[ixNode];

	// Loop through all faces
	std::set<int>::const_iterator iter = vecNearbyFaces.begin();
	for (; iter != vecNearbyFaces.end(); iter++) {

		Face::NodeLocation loc;
		int ixLocation;

		mesh.faces[*iter].ContainsNode(
			mesh.nodes,
			nodeNudged,
			loc,
			ixLocation);

		if (loc == Face::NodeLocation_Exterior) {
			continue;
		} else {
			vecFaceIndices.push_back(*iter);
		}
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
	int ixCurrentFirstFace = 0;
	for (; ixCurrentFirstFace < meshFirst.faces.size(); ixCurrentFirstFace++) {
	//for (int k = 36; k < 37; k++) {

		// Current face
		const Face & faceFirstCurrent = meshFirst.faces[ixCurrentFirstFace];

		// Starting point
		Node nodeCurrent =
			nodevecFirst[meshFirst.faces[ixCurrentFirstFace][0]];

		// Find the starting face on the second mesh
		std::vector<int> vecStartingFaces;
		FindFaceFromNode(meshSecond, nodeCurrent, vecStartingFaces);

		// No faces found
		if (vecStartingFaces.size() == 0) {
			_EXCEPTIONT("No initial face found");
		}

		// Current face on second mesh
		int ixCurrentSecondFace = vecStartingFaces[0];

		printf("Faces: %i %i\n", ixCurrentFirstFace, ixCurrentSecondFace);

		// This node lies on the boundary between faces; find the
		// actual starting face by nudging along FirstEdge.
		if (vecStartingFaces.size() > 1) {

			std::vector<int> vecFaceIndices;
			FindFaceNearNode(
				meshSecond,
				nodeCurrent,
				nodevecFirst[faceFirstCurrent[1]],
				faceFirstCurrent.edges[0].type,
				vecStartingFaces,
				vecFaceIndices);

			if (vecFaceIndices.size() == 0) {
				_EXCEPTIONT("No nearby face found");
			} else if (vecFaceIndices.size() == 1) {
				ixCurrentSecondFace = vecFaceIndices[0];
			} else {
				_EXCEPTIONT("Coincident edges away from node!");
			}
		}

		// Generate a path
		PathSegmentVector vecTracedPath;

		// Trace along all edges of current face
		for (int i = 0; i < faceFirstCurrent.edges.size(); i++) {

			// Equal node indices indicate a non-edge
			if (faceFirstCurrent.edges[i][0] == faceFirstCurrent.edges[i][1]) {
				continue;
			}

			// Initialize the trace
			const Edge & edgeFirstCurrent = faceFirstCurrent.edges[i];

			const Node & nodeFirstEnd = nodevecFirst[edgeFirstCurrent[1]];

			int ixOverlapNodeCurrent = edgeFirstCurrent[0];
			Node nodeCurrent = nodevecFirst[ixOverlapNodeCurrent];

			// Repeat until we hit the end of this edge
			for (;;) {

				const Face & faceSecondCurrent =
					meshSecond.faces[ixCurrentSecondFace];

				// Find all intersections between this edge and the
				// second face
				bool fCoincidentEdge = false;

				int ixIntersectionSecondEdge;
				std::vector<Node> nodeIntersections;

				for (int j = 0; j < faceSecondCurrent.edges.size(); j++) {
					const Edge & edgeSecondCurrent =
						faceSecondCurrent.edges[j];

					// Equal node indices indicating a non-edge
					if (edgeSecondCurrent[0] == edgeSecondCurrent[1]) {
						continue;
					}

					fCoincidentEdge =
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
					vecTracedPath.push_back(PathSegment(
						ixOverlapNodeCurrent,
						edgeFirstCurrent[1],
						edgeFirstCurrent.type,
						ixCurrentFirstFace,
						ixCurrentSecondFace));

					break;
				}

				// Coincident edges
				if (fCoincidentEdge) {
					_EXCEPTIONT("Not implemented");
				}

				// Invalid option
				if (nodeIntersections.size() > 1) {
					_EXCEPTIONT("Logic Error");
				}

				// Check if the path hits a corner node
				const Edge & edgeSecondCurrent =
					faceSecondCurrent.edges[ixIntersectionSecondEdge];

				const Node & nodeSecondEdge0 =
					nodevecSecond[edgeSecondCurrent[0]];
				const Node & nodeSecondEdge1 =
					nodevecSecond[edgeSecondCurrent[1]];

				// Push new intersection into overlap array
				int ixOverlapNodeNext =
					static_cast<int>(meshOverlap.nodes.size());

				if (nodeIntersections[0] == nodeSecondEdge0) {
					ixOverlapNodeNext =
						ixOverlapSecondNodesBegin + edgeSecondCurrent[0];

#pragma message "What if nodeIntersections[0] is also the FirstEnd?"
					std::vector<int> vecNearbyFaces;
					FindFaceNearNode(
						meshSecond,
						edgeSecondCurrent[0],
						nodeFirstEnd,
						edgeFirstCurrent.type,
						vecNearbyFaces);

					if (vecNearbyFaces.size() == 0) {
						_EXCEPTIONT("No nearby faces found!");
					} else if (vecNearbyFaces.size() == 1) {
						ixCurrentSecondFace = vecNearbyFaces[0];
					} else if (vecNearbyFaces.size() > 1) {
						_EXCEPTIONT("Not implemented");
					}

					vecTracedPath.push_back(PathSegment(
						ixOverlapNodeCurrent,
						ixOverlapNodeNext,
						edgeFirstCurrent.type,
						ixCurrentFirstFace,
						ixCurrentSecondFace));

					ixOverlapNodeCurrent = ixOverlapNodeNext;

					continue;
				}
				if (nodeIntersections[0] == nodeSecondEdge1) {
					ixOverlapNodeNext =
						ixOverlapSecondNodesBegin + edgeSecondCurrent[1];

#pragma message "What if nodeIntersections[0] is also the FirstEnd?"
					std::vector<int> vecNearbyFaces;
					FindFaceNearNode(
						meshSecond,
						edgeSecondCurrent[1],
						nodeFirstEnd,
						edgeFirstCurrent.type,
						vecNearbyFaces);

					if (vecNearbyFaces.size() == 0) {
						_EXCEPTIONT("No nearby faces found!");
					} else if (vecNearbyFaces.size() == 1) {
						ixCurrentSecondFace = vecNearbyFaces[0];
					} else if (vecNearbyFaces.size() > 1) {
						_EXCEPTIONT("Not implemented");
					}

					vecTracedPath.push_back(PathSegment(
						ixOverlapNodeCurrent,
						ixOverlapNodeNext,
						edgeFirstCurrent.type,
						ixCurrentFirstFace,
						ixCurrentSecondFace));

					ixOverlapNodeCurrent = ixOverlapNodeNext;

					continue;
				}

				// Intersection is exactly the end point of this edge
				if (nodeIntersections[0] ==
				        meshOverlap.nodes[edgeFirstCurrent[1]]
				) {
					vecTracedPath.push_back(PathSegment(
						ixOverlapNodeCurrent,
						edgeFirstCurrent[1],
						edgeFirstCurrent.type,
						ixCurrentFirstFace,
						ixCurrentSecondFace));

					// FIX: Determine if we need to cross Second faces

					break;
				}

				// Push a new intersection into the array of nodes
				meshOverlap.nodes.push_back(nodeIntersections[0]);

				// Intersection found; advance to next face
				vecTracedPath.push_back(PathSegment(
					ixOverlapNodeCurrent,
					ixOverlapNodeNext,
					edgeFirstCurrent.type,
					ixCurrentFirstFace,
					ixCurrentSecondFace));

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
	}

	meshOverlap.faces.push_back(meshSecond.faces[8]);
	meshOverlap.faces.push_back(meshFirst.faces[2]);
	for (int i = 0; i < 4; i++) {
		meshOverlap.faces[0].edges[i][0] += ixOverlapSecondNodesBegin;
	}

}

///////////////////////////////////////////////////////////////////////////////

