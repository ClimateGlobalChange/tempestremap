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
///		Enumerator defining a type of intersection.
///	</summary>
enum IntersectType {
	IntersectType_None,
	IntersectType_Edge,
	IntersectType_Node
};

///	<summary>
///		A segment connecting two nodes that also has an associated face
///		on the First and Second mesh.
///	</summary>
class PathSegment : public Edge {

public:
	///	<summary>
	///		Constructor for node intersections.
	///	</summary>
	PathSegment(
		int a_node0,
		int a_node1,
		Edge::Type a_type,
		int a_ixFirstFace,
		int a_ixSecondFace,
		IntersectType a_inttype,
		int a_ixIntersect
	) :
		Edge(a_node0, a_node1, a_type),
		ixFirstFace(a_ixFirstFace),
		ixSecondFace(a_ixSecondFace),
		inttype(a_inttype),
		ixIntersect(a_ixIntersect)
	{
	}

	///	<summary>
	///		Constructor for edge intersections.
	///	</summary>
	PathSegment(
		int a_node0,
		int a_node1,
		Edge::Type a_type,
		int a_ixFirstFace,
		int a_ixSecondFace,
		int a_ixIntersect,
		const Edge & a_edgeIntersect
	) :
		Edge(a_node0, a_node1, a_type),
		ixFirstFace(a_ixFirstFace),
		ixSecondFace(a_ixSecondFace),
		ixIntersect(a_ixIntersect),
		edgeIntersect(a_edgeIntersect)
	{
		inttype = IntersectType_Edge;
	}

public:
	///	<summary>
	///		Origin face on first mesh.
	///	</summary>
	int ixFirstFace;

	///	<summary>
	///		Origin face on second mesh.
	///	</summary>
	int ixSecondFace;

	///	<summary>
	///		Type of intersection that occurs to end this PathSegment.
	///	</summary>
	IntersectType inttype;

	///	<summary>
	///		This is the local index of the edge which one will hit if moving
	///		counter-clockwise around ixSecondFace.
	///	</summary>
	int ixIntersect;

	///	<summary>
	///		If inttype is IntersectType_Edge, this is the index of the edge
	///		on the SecondMesh which has been intersected.
	///	</summary>
	Edge edgeIntersect;
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

int FindFaceNearNode(
	const Mesh & mesh,
	const Node & nodeBegin,
	const Node & nodeEnd,
	const Edge::Type edgetype,
	const std::vector<int> & vecPossibleFaces
) {
	// Find the nudged node towards nodeEnd with type edgetype
	Node nodeNudged;
	NudgeAlongEdge(nodeBegin, nodeEnd, edgetype, nodeNudged);

	// Loop through all faces
	for (int i = 0; i < vecPossibleFaces.size(); i++) {

		Face::NodeLocation loc;
		int ixLocation;

		int ixPossibleFace = vecPossibleFaces[i];

		const Face & face = mesh.faces[ixPossibleFace];

		face.ContainsNode(
			mesh.nodes,
			nodeNudged,
			loc,
			ixLocation);

		if (loc == Face::NodeLocation_Exterior) {
			continue;

		} else if (loc == Face::NodeLocation_Edge) {
			const Edge & edge = face.edges[ixLocation];
			const Node & node0 = mesh.nodes[edge[0]];
			const Node & node1 = mesh.nodes[edge[1]];

			Node nodeDeltaA(
				nodeNudged.x - nodeBegin.x,
				nodeNudged.y - nodeBegin.y,
				nodeNudged.z - nodeBegin.z);

			Node nodeDeltaB(
				node1.x - node0.x,
				node1.y - node0.y,
				node1.z - node0.z);

			printf("%1.5e %1.5e %1.5e : %1.5e %1.5e %1.5e\n",
				nodeDeltaA.x, nodeDeltaA.y, nodeDeltaA.z,
				nodeDeltaB.x, nodeDeltaB.y, nodeDeltaB.z);

			_EXCEPTION();

		} else if (loc == Face::NodeLocation_Corner) {
			_EXCEPTIONT("Coincident nodes away from fixed node!");

		} else {
			return ixPossibleFace;
		}
	}
	return 0;
}

///////////////////////////////////////////////////////////////////////////////

int FindFaceNearNode(
	const Mesh & mesh,
	int ixNode,
	const Node & nodeEnd,
	const Edge::Type edgetype
) {
	const Node & nodeBegin = mesh.nodes[ixNode];

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

		const Face & face = mesh.faces[*iter];

		face.ContainsNode(
			mesh.nodes,
			nodeNudged,
			loc,
			ixLocation);

		if (loc == Face::NodeLocation_Exterior) {
			continue;

		} else if (loc == Face::NodeLocation_Edge) {
			const Edge & edge = face.edges[ixLocation];
			const Node & node0 = mesh.nodes[edge[0]];
			const Node & node1 = mesh.nodes[edge[1]];

			Node nodeDeltaA(
				nodeEnd.x - nodeBegin.x,
				nodeEnd.y - nodeBegin.y,
				nodeEnd.z - nodeBegin.z);

			Node nodeDeltaB(
				node1.x - node0.x,
				node1.y - node0.y,
				node1.z - node0.z);

			double dDot =
				+ nodeDeltaA.x * nodeDeltaB.x
				+ nodeDeltaA.y * nodeDeltaB.y
				+ nodeDeltaA.z * nodeDeltaB.z;

			if (dDot > 0.0) {
				return (*iter);
			}

		} else if (loc == Face::NodeLocation_Corner) {
			_EXCEPTIONT("Coincident nodes away from fixed node!");

		} else {
			return (*iter);
		}
	}

	_EXCEPTIONT("No face found!");
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
	const NodeVector & nodevecOverlap = meshOverlap.nodes;

	// Insert all nodes from the two NodeVectors
	for (int i = 0; i < nodevecFirst.size(); i++) {
		meshOverlap.nodes.push_back(nodevecFirst[i]);
	}
	const int ixOverlapSecondNodesBegin = meshOverlap.nodes.size();

	for (int i = 0; i < nodevecSecond.size(); i++) {
		meshOverlap.nodes.push_back(nodevecSecond[i]);
	}
	const int ixOverlapNewNodesBegin = meshOverlap.nodes.size();

	int ixCurrentFirstFace = 0;
	//for (; ixCurrentFirstFace < meshFirst.faces.size(); ixCurrentFirstFace++) {
	for (int ixCurrentFirstFace = 0; ixCurrentFirstFace < 9; ixCurrentFirstFace++) {

		///////////////////////////////////////////////////////////////////////
		// Build the PathSegmentVector around this face

		// Generate a path
		PathSegmentVector vecTracedPath;

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
			_EXCEPTIONT("No initial face found!");
		}

		// Current face on second mesh
		int ixCurrentSecondFace = vecStartingFaces[0];

		// This node lies on the boundary between faces; find the
		// actual starting face by nudging along FirstEdge.
		if (vecStartingFaces.size() > 1) {

			std::vector<int> vecFaceIndices;
			ixCurrentSecondFace =
				FindFaceNearNode(
					meshSecond,
					nodeCurrent,
					nodevecFirst[faceFirstCurrent[1]],
					faceFirstCurrent.edges[0].type,
					vecStartingFaces);
		}

		printf("\nFaces: %i %i\n", ixCurrentFirstFace, ixCurrentSecondFace);

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

				// Index within faceSecondCurrent of the intersection
				// Node or Edge.
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
						nodeIntersections.clear();
						continue;
					}

					if (nodeIntersections.size() > 1) {
						_EXCEPTIONT("Not implemented: Non-convex intersections");
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
						ixCurrentSecondFace,
						IntersectType_None,
						0));

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
					vecTracedPath.push_back(PathSegment(
						ixOverlapNodeCurrent,
						ixOverlapNodeNext,
						edgeFirstCurrent.type,
						ixCurrentFirstFace,
						ixCurrentSecondFace,
						IntersectType_Node,
						ixIntersectionSecondEdge));

					ixCurrentSecondFace =
						FindFaceNearNode(
							meshSecond,
							edgeSecondCurrent[0],
							nodeFirstEnd,
							edgeFirstCurrent.type);

					ixOverlapNodeCurrent = ixOverlapNodeNext;

					continue;
				}
				if (nodeIntersections[0] == nodeSecondEdge1) {
					ixOverlapNodeNext =
						ixOverlapSecondNodesBegin + edgeSecondCurrent[1];

#pragma message "What if nodeIntersections[0] is also the FirstEnd?"

					vecTracedPath.push_back(PathSegment(
						ixOverlapNodeCurrent,
						ixOverlapNodeNext,
						edgeFirstCurrent.type,
						ixCurrentFirstFace,
						ixCurrentSecondFace,
						IntersectType_Node,
						(ixIntersectionSecondEdge + 1)
							% faceSecondCurrent.edges.size()));

					ixCurrentSecondFace =
						FindFaceNearNode(
							meshSecond,
							edgeSecondCurrent[1],
							nodeFirstEnd,
							edgeFirstCurrent.type);

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
						ixCurrentSecondFace,
						IntersectType_None,
						0));

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
					ixCurrentSecondFace,
					ixIntersectionSecondEdge,
					edgeSecondCurrent));

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

		///////////////////////////////////////////////////////////////////////
		// Build the SubMesh associated with this face
/*
		// Build a set of all nodes on the boundary of faceFirstCurrent
		std::map<int, int> mapExitNodes;
		std::map<int, int> mapExitEdges;
		for (int j = 0; j < vecTracedPath.size(); j++) {
			if (vecTracedPath[j].inttype == IntersectType_Node) {
				mapExitNodes.insert(
					std::pair<int,int>(
						vecTracedPath[j][0],
						(j+1) % vecTracedPath.size()));
			}
			if (vecTracedPath[j].inttype == IntersectType_Edge) {
				mapExitEdges.insert(
					std::pair<int,int>(
						vecTracedPath[j][0],
						(j+1) % vecTracedPath.size()));
			}
		}
*/
		// Array indicating which elements of vecTracedPath have been used
		std::vector<bool> vecTracedPathUsed;
		vecTracedPathUsed.resize(vecTracedPath.size(), false);

		// The set of interior faces (from the Second mesh)
		std::set<int> setInteriorFaces;
		for (int j = 0; j < vecTracedPath.size(); j++) {
			setInteriorFaces.insert(vecTracedPath[j].ixSecondFace);
		}

		// Loop through all possible starting PathSegments
		for (int j = 0; j < vecTracedPath.size(); j++) {

			// Ignore starting PathSegments that have already been used
			if (vecTracedPathUsed[j]) {
				continue;
			}

			// Build the new face
			Face faceOverlap(0);

			// Origin node of this face
			int ixOverlapOriginNode = vecTracedPath[j][0];

			// Current second face
			int ixCurrentSecondFace =
				vecTracedPath[j].ixSecondFace;

			const Face & faceSecondCurrent =
				meshSecond.faces[ixCurrentSecondFace];

			// Search may require multiple trips along FirstMesh and SecondMesh
			for (;;) {

				// Loop along edge of FirstFace until we find a branch
				// into SecondMesh or hit our origin.
				for (;; j++) {
					faceOverlap.edges.push_back(vecTracedPath[j]);

					// Mark traced path edge as used
					if (vecTracedPathUsed[j]) {
						_EXCEPTIONT("Trying to reuse traced path edge");
					}

					vecTracedPathUsed[j] = true;

					printf("P%i: %i %i\n",
						j, vecTracedPath[j][0], vecTracedPath[j][1]);

					// Found a branch into the interior
					if (vecTracedPath[j].inttype != IntersectType_None) {
						break;
					}

					// Hit origin; finish this face
					if (vecTracedPath[j][1] == ixOverlapOriginNode) {
						goto ContinueToNextFace;
					}
				}

				// Determine the index of intersection
				int ixCurrentSecondEdge = vecTracedPath[j].ixIntersect;

				int ixCurrentOverlapNode = vecTracedPath[j][1];

				// Loop around the interior of faceSecondCurrent
				for (;;) {

					const Edge & edgeSecondCurrent =
						faceSecondCurrent.edges[ixCurrentSecondEdge];

					// Determine if this edge exits onto FirstMesh
					int ixExitNode = InvalidNode;

					int k;
					for (
						k = (j + 1) % vecTracedPath.size();
						k != j;
						k = (k + 1) % vecTracedPath.size()
					) {
						int ixIntersectNode =
							faceSecondCurrent[vecTracedPath[k].ixIntersect];

						IntersectType inttype = vecTracedPath[k].inttype;
/*
						printf("%i %i %i\n",
							vecTracedPath[k][1],
							ixOverlapSecondNodesBegin + edgeSecondCurrent[0],
							ixOverlapSecondNodesBegin + edgeSecondCurrent[1]);
*/
						// Check for node intersections
						if (inttype == IntersectType_Node) {
							if (vecTracedPath[k][1] ==
								ixOverlapSecondNodesBegin
									+ edgeSecondCurrent[0]
							) {
								ixExitNode = vecTracedPath[k][1];
								break;
							}

							if (vecTracedPath[k][1] ==
								ixOverlapSecondNodesBegin
									+ edgeSecondCurrent[1]
							) {
								ixExitNode = vecTracedPath[k][1];
								break;
							}
						}

						// Check for edge intersections
						if ((inttype == IntersectType_Edge) &&
							(edgeSecondCurrent == vecTracedPath[k].edgeIntersect)
						) {
							ixExitNode = vecTracedPath[k][1];
							break;
						}
					}

					// This edge exits; push the edge into the Face and loop
					if (ixExitNode != InvalidNode) {

						printf("S: %i %i\n",
							ixCurrentOverlapNode, ixExitNode);

						faceOverlap.edges.push_back(Edge(
							ixCurrentOverlapNode,
							ixExitNode,
							edgeSecondCurrent.type));

#pragma message "FIX: What if edge only touches meshFirst then reverses back into meshSecond?"
						j = (k + 1) % vecTracedPath.size();

						if (ixExitNode == ixOverlapOriginNode) {
							goto ContinueToNextFace;
						} else {
							break;
						}
					}

					printf("T: %i (%i) %i\n",
						ixCurrentOverlapNode,
						ixOverlapSecondNodesBegin + edgeSecondCurrent[0],
						ixOverlapSecondNodesBegin + edgeSecondCurrent[1]);

					// Push this edge into the overlap face
					faceOverlap.edges.push_back(Edge(
						ixCurrentOverlapNode,
						ixOverlapSecondNodesBegin + edgeSecondCurrent[1],
						edgeSecondCurrent.type));

					// Advance the edge
					ixCurrentSecondEdge =
						(ixCurrentSecondEdge + 1)
							% faceSecondCurrent.edges.size();

					ixCurrentOverlapNode =
						ixOverlapSecondNodesBegin + edgeSecondCurrent[1];

					if (ixCurrentOverlapNode == ixOverlapOriginNode) {
						goto ContinueToNextFace;
					}
				}
			}

ContinueToNextFace:
			printf("PUSH %lu\n", faceOverlap.edges.size());
			// Push this Face into the overlap Mesh
			meshOverlap.faces.push_back(faceOverlap);

			// Reset the search counter
			j = 0;

			// Check if all segments have been used
			int nUsedSegments = 0;
			for (int k = 0; k < vecTracedPath.size(); k++) {
				if (vecTracedPathUsed[k]) {
					nUsedSegments++;
				}
			}
			if (nUsedSegments == static_cast<int>(vecTracedPath.size())) {
				break;
			}
		}

		// FIX:REMOVE
		//break;
	}
/*
	meshOverlap.faces.clear();
	meshOverlap.faces.push_back(meshSecond.faces[8]);
	meshOverlap.faces.push_back(meshFirst.faces[1]);
	for (int i = 0; i < 4; i++) {
		meshOverlap.faces[0].edges[i][0] += ixOverlapSecondNodesBegin;
	}
*/
}

///////////////////////////////////////////////////////////////////////////////

