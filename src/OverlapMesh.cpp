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

#include "Defines.h"
#include "OverlapMesh.h"
#include "MeshUtilitiesFuzzy.h"
#include "MeshUtilitiesExact.h"

#include "Announce.h"

#include "kdtree.h"

#include <unistd.h>
#include <iostream>
#include <queue>

///////////////////////////////////////////////////////////////////////////////

#define VERBOSE

//#define CHECK_AREAS

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
///		on the Source and Target mesh.
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
		int a_ixSourceFace,
		int a_ixTargetFace,
		IntersectType a_inttype,
		int a_ixIntersect
	) :
		Edge(a_node0, a_node1, a_type),
		ixSourceFace(a_ixSourceFace),
		ixTargetFace(a_ixTargetFace),
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
		int a_ixSourceFace,
		int a_ixTargetFace,
		int a_ixIntersect,
		const Edge & a_edgeIntersect
	) :
		Edge(a_node0, a_node1, a_type),
		ixSourceFace(a_ixSourceFace),
		ixTargetFace(a_ixTargetFace),
		ixIntersect(a_ixIntersect),
		edgeIntersect(a_edgeIntersect)
	{
		inttype = IntersectType_Edge;
	}

public:
	///	<summary>
	///		Origin face on source mesh.
	///	</summary>
	int ixSourceFace;

	///	<summary>
	///		Origin face on target mesh.
	///	</summary>
	int ixTargetFace;

	///	<summary>
	///		Type of intersection that occurs to end this PathSegment.
	///	</summary>
	IntersectType inttype;

	///	<summary>
	///		This is the local index of the edge which one will hit if moving
	///		counter-clockwise around ixTargetFace.
	///	</summary>
	int ixIntersect;

	///	<summary>
	///		If inttype is IntersectType_Edge, this is the index of the edge
	///		on the TargetMesh which has been intersected.
	///	</summary>
	Edge edgeIntersect;
};

///	<summary>
///		A vector of PathSegments.
///	</summary>
typedef std::vector<PathSegment> PathSegmentVector;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Generate a PathSegmentVector describing the path around the face
///		ixCurrentSourceFace.
///	</summary>
template <
	class MeshUtilities,
	class NodeIntersectType
>
void GeneratePath(
	const Mesh & meshSource,
	const Mesh & meshTarget,
	const std::vector<int> & vecTargetNodeMap,
	int ixCurrentSourceFace,
	PathSegmentVector & vecTracedPath,
	Mesh & meshOverlap
) {

	// MeshUtilities object
	MeshUtilities utils;

	// Get the NodeVectors
	const NodeVector & nodevecSource = meshSource.nodes;
	const NodeVector & nodevecTarget = meshTarget.nodes;

	NodeVector & nodevecOverlap = meshOverlap.nodes;

	// Current face
	const Face & faceSourceCurrent = meshSource.faces[ixCurrentSourceFace];

	// Starting point
	Node nodeCurrent =
		nodevecSource[faceSourceCurrent[0]];

	// Find the starting face on the target mesh
	FindFaceStruct aFindFaceStruct;
	utils.FindFaceFromNode(
		meshTarget,
		nodeCurrent,
		aFindFaceStruct);

	// No faces found
	if (aFindFaceStruct.vecFaceIndices.size() == 0) {
		nodeCurrent.Print("node");
		_EXCEPTIONT("No initial face found:\n"
			"    Mesh B must be a submesh of mesh A");
	}

	// Current face on target mesh
	int ixCurrentTargetFace = aFindFaceStruct.vecFaceIndices[0];

	// Verify starting Node is not on the Exterior
	if (aFindFaceStruct.loc == Face::NodeLocation_Exterior) {
		_EXCEPTIONT("Invalid Node location returned from FindFaceFromNode");

	// Starting Node is on the interior of a Face on meshTarget
	} else if (aFindFaceStruct.loc == Face::NodeLocation_Interior) {
		ixCurrentTargetFace = aFindFaceStruct.vecFaceIndices[0];

	// Starting Node is on the Node of a Face on meshTarget
	} else {
/*
		nodevecSource[faceSourceCurrent[0]].PrintX("0");
		nodevecSource[faceSourceCurrent[1]].PrintX("1");

		for (int i = 0; i < aFindFaceStruct.vecFaceIndices.size(); i++) {
			printf("%i\n", aFindFaceStruct.vecFaceIndices[i]);
		}
*/
		ixCurrentTargetFace =
			utils.FindFaceNearNode(
				meshTarget,
				//nodeCurrent,
				nodevecSource[faceSourceCurrent[0]],
				nodevecSource[faceSourceCurrent[1]],
				faceSourceCurrent.edges[0].type,
				aFindFaceStruct);

		const Face & face0 = meshTarget.faces[ixCurrentTargetFace];

	}

	// Starting information
	printf("\nFaces: %i %i\n", ixCurrentSourceFace, ixCurrentTargetFace);
#ifdef VERBOSE
	printf("Starting Node: %i\n", meshSource.faces[ixCurrentSourceFace][0]);
	printf("Next Node: %i\n", meshSource.faces[ixCurrentSourceFace][1]);
#endif

	// Trace along all edges of current face
	for (int i = 0; i < faceSourceCurrent.edges.size(); i++) {

		// Equal node indices indicate a non-edge
		if (faceSourceCurrent.edges[i][0] == faceSourceCurrent.edges[i][1]) {
			continue;
		}

		// Initialize the trace
		const Edge & edgeSourceCurrent = faceSourceCurrent.edges[i];

		const Node & nodeSourceBegin = nodevecSource[edgeSourceCurrent[0]];
		const Node & nodeSourceEnd = nodevecSource[edgeSourceCurrent[1]];

		int ixOverlapNodeCurrent = edgeSourceCurrent[0];
		Node nodeCurrent = nodevecSource[ixOverlapNodeCurrent];

		NodeIntersectType nodeLastIntersection = nodeSourceBegin;

		// Repeat until we hit the end of this edge
		for (;;) {

			const Face & faceTargetCurrent =
				meshTarget.faces[ixCurrentTargetFace];

#ifdef VERBOSE
			printf("--- (%i) ---\n", ixCurrentTargetFace);
#endif

/*
			meshOverlap.nodes[edgeSourceCurrent[0]].PrintMX();
			meshOverlap.nodes[edgeSourceCurrent[1]].PrintMX();
*/
			// Find all intersections between this edge and the
			// second face
			bool fCoincidentEdge = false;

			// Index within faceTargetCurrent of the intersection
			// Node or Edge.
			int ixIntersectionTargetEdge;

			std::vector<NodeIntersectType> nodeIntersections;

			for (int j = 0; j < faceTargetCurrent.edges.size(); j++) {
				const Edge & edgeTargetCurrent =
					faceTargetCurrent.edges[j];

				// Equal node indices indicating a zero edge
				if (edgeTargetCurrent[0] == edgeTargetCurrent[1]) {
					_EXCEPTIONT("Zero Edge detected");
				}
/*
				fCoincidentEdge =
					utilsFuzzy.CalculateEdgeIntersections(
				 		//meshOverlap.nodes[ixOverlapNodeCurrent],
						meshOverlap.nodes[edgeSourceCurrent[0]],
						meshOverlap.nodes[edgeSourceCurrent[1]],
						edgeSourceCurrent.type,
						nodevecTarget[edgeTargetCurrent[0]],
						nodevecTarget[edgeTargetCurrent[1]],
						edgeTargetCurrent.type,
						nodeIntersections);

				printf(" - - - \n");
*/
				if (edgeTargetCurrent[0] < edgeTargetCurrent[1]) {
					fCoincidentEdge =
						utils.CalculateEdgeIntersections(
							nodevecOverlap[edgeSourceCurrent[0]],
							nodevecOverlap[edgeSourceCurrent[1]],
							edgeSourceCurrent.type,
							nodevecTarget[edgeTargetCurrent[0]],
							nodevecTarget[edgeTargetCurrent[1]],
							edgeTargetCurrent.type,
							nodeIntersections);
				} else {
					fCoincidentEdge =
						utils.CalculateEdgeIntersections(
							nodevecOverlap[edgeSourceCurrent[0]],
							nodevecOverlap[edgeSourceCurrent[1]],
							edgeSourceCurrent.type,
							nodevecTarget[edgeTargetCurrent[1]],
							nodevecTarget[edgeTargetCurrent[0]],
							edgeTargetCurrent.type,
							nodeIntersections);
				}

				for (int i = 0; i < nodeIntersections.size(); i++) {
					bool fEqualNodes =
						utils.AreNodesEqual(
							nodeIntersections[i],
							nodeLastIntersection);

/*
					meshOverlap.nodes[edgeSourceCurrent[1]].Print("f1");
					nodeIntersections[i].Print("n0");
					nodeLastIntersection.Print("l0");
*/
					if (fEqualNodes) {
						nodeIntersections.erase(nodeIntersections.begin()+i);
						i--;
					}

/*
						bool fBeginPoint =
							utils.AreNodesEqual(
								nodeIntersections[i],
								meshOverlap.nodes[edgeSourceCurrent[0]]);

						if (!fBeginPoint) {
							nodeIntersections.erase(nodeIntersections.begin()+i);
							i--;
						}
					}
*/
				}

				if (fCoincidentEdge) {
					nodeIntersections.clear();
					continue;
				}

				if (nodeIntersections.size() > 1) {
					_EXCEPTIONT("Not implemented: Non-convex intersections");
				}

				if (nodeIntersections.size() == 1) {
					ixIntersectionTargetEdge = j;
					break;
				}
			}

			// Done with this edge
			if (nodeIntersections.size() == 0) {

				// DEBUG: Verify that endpoint is in this element
				Face::NodeLocation loc;
				int ixLocation;

				utils.ContainsNode(
					meshTarget.faces[ixCurrentTargetFace],
					meshTarget.nodes,
					nodevecOverlap[edgeSourceCurrent[1]],
					loc,
					ixLocation);

				if (loc == Face::NodeLocation_Exterior) {
					_EXCEPTIONT("Edge endpoint failure");
				}

				printf("Done with Edge (%i)\n", loc);

				// Add the final PathSegment for this Edge
				vecTracedPath.push_back(PathSegment(
					ixOverlapNodeCurrent,
					edgeSourceCurrent[1],
					edgeSourceCurrent.type,
					ixCurrentSourceFace,
					ixCurrentTargetFace,
					IntersectType_None,
					0));

				break;
			}

/*
			meshOverlap.nodes[edgeSourceCurrent[0]].PrintX("Begin");
			nodeIntersections[0].PrintX("Int");
			meshOverlap.nodes[edgeSourceCurrent[1]].PrintX("End");
*/
			// Coincident edges
			if (fCoincidentEdge) {
				_EXCEPTIONT("Not implemented");
			}

			// Invalid option
			if (nodeIntersections.size() > 1) {
				_EXCEPTIONT("Logic Error");
			}

			// Set last intersection
			nodeLastIntersection = nodeIntersections[0];

			// Find next face on meshTarget
			const Edge & edgeTargetCurrent =
				faceTargetCurrent.edges[ixIntersectionTargetEdge];

			const Node & nodeTargetEdge0 =
				nodevecTarget[edgeTargetCurrent[0]];
			const Node & nodeTargetEdge1 =
				nodevecTarget[edgeTargetCurrent[1]];

			int ixOverlapNodeNext;

			// Special case:  Intersection with Edge is exactly an
			// beginpoint / endpoint of edgeSourceCurrent
/*
			nodeIntersections[0].Print("n");
			meshOverlap.nodes[edgeSourceCurrent[1]].Print("1");
*/
/*
			nodeIntersections[0].PrintMX();
			nodeIntersections[0].PrintNorm();
			nodeTargetEdge0.PrintMX();
			nodeTargetEdge1.PrintMX();
			meshOverlap.nodes[edgeSourceCurrent[0]].PrintMX();
			meshOverlap.nodes[edgeSourceCurrent[1]].PrintMX();
*/
/*
			bool fIntersectEndpoint =
				utils.AreNodesEqual(
					nodeIntersections[0],
					meshOverlap.nodes[edgeSourceCurrent[1]]);
*/
			// Determine if this Face contains the endpoint of the Edge
			Face::NodeLocation loc0;
			int ixLocation0;

			utils.ContainsNode(
				meshTarget.faces[ixCurrentTargetFace],
				meshTarget.nodes,
				nodevecOverlap[edgeSourceCurrent[1]],
				loc0,
				ixLocation0);

			bool fIntersectEndpoint = false;
			if ((loc0 == Face::NodeLocation_Edge) ||
				(loc0 == Face::NodeLocation_Corner)
			) {
				fIntersectEndpoint = true;
			}

/*
			if (dynamic_cast<MeshUtilitiesExact*>(&utils) != NULL) {
				FixedPoint fp1 =
					DotProductX(nodeTargetEdge0,
						CrossProductX(
							nodeIntersections[0],
							meshOverlap.nodes[edgeSourceCurrent[1]]));

				if (fp1.IsZero()) {
					fIntersectEndpoint = true;
				}
			}
*/
			if (fIntersectEndpoint) {

				// Next edge
				int iNext = (i + 1) % faceSourceCurrent.edges.size();

				const Edge & edgeSourceNext =
					faceSourceCurrent.edges[iNext];

				// Update TargetMesh face
				EdgeMapConstIterator iter =
					meshTarget.edgemap.find(edgeTargetCurrent);

				if (iter == meshTarget.edgemap.end()) {
					_EXCEPTIONT("Logic error");
				}

				int ixNextTargetFace;

/*
				nodeIntersections[0].Print("i0");
				nodeSourceEnd.Print("i0");
				nodeTargetEdge0.Print("n0");
				nodeTargetEdge1.Print("n1");
*/

				// Path hits the beginpoint of the Edge
				if (utils.AreNodesEqual(nodeIntersections[0], nodeTargetEdge0)) {
					ixNextTargetFace =
						utils.FindFaceNearNode(
							meshTarget,
							edgeTargetCurrent[0],
							nodevecSource[edgeSourceNext[1]],
							edgeSourceNext.type);
#ifdef VERBOSE
					if (ixNextTargetFace == ixCurrentTargetFace) {
						printf("WARNING: Face does not change across Edge (1)\n");
					}
#endif
					// If face changes insert a Node bifurcation
					if (ixNextTargetFace != ixCurrentTargetFace) {
						vecTracedPath.push_back(PathSegment(
							ixOverlapNodeCurrent,
							edgeSourceCurrent[1],
							edgeSourceCurrent.type,
							ixCurrentSourceFace,
							ixCurrentTargetFace,
							IntersectType_Node,
							ixIntersectionTargetEdge));
					}

				// Path hits the endpoint of the Edge
				} else if (
					utils.AreNodesEqual(nodeIntersections[0], nodeTargetEdge1)
				) {
					ixNextTargetFace =
						utils.FindFaceNearNode(
							meshTarget,
							edgeTargetCurrent[1],
							nodevecSource[edgeSourceNext[1]],
							edgeSourceNext.type);
#ifdef VERBOSE
					if (ixNextTargetFace == ixCurrentTargetFace) {
						printf("WARNING: Face does not change across Edge (2)\n");
					}
#endif
					// If face changes insert a Node bifuraction
					if (ixNextTargetFace != ixCurrentTargetFace) {
						vecTracedPath.push_back(PathSegment(
							ixOverlapNodeCurrent,
							edgeSourceCurrent[1],
							edgeSourceCurrent.type,
							ixCurrentSourceFace,
							ixCurrentTargetFace,
							IntersectType_Node,
							(ixIntersectionTargetEdge + 1)
								% faceTargetCurrent.edges.size()));
					}

				// Path hits the Edge directly
				} else {
					const FacePair & facepair = iter->second;

					// Find the set of possible Faces
					std::vector<int> vecPossibleFaces;
					vecPossibleFaces.push_back(facepair[0]);
					vecPossibleFaces.push_back(facepair[1]);

					// Build the FindFaceStruct
					FindFaceStruct aNextFindFaceStruct;

					aNextFindFaceStruct.vecFaceIndices.push_back(facepair[0]);
					aNextFindFaceStruct.vecFaceIndices.push_back(facepair[1]);

					const Face & face0 = meshTarget.faces[facepair[0]];
					const Face & face1 = meshTarget.faces[facepair[1]];

					int iEdge0 = face0.GetEdgeIndex(edgeTargetCurrent);
					int iEdge1 = face1.GetEdgeIndex(edgeTargetCurrent);

					aNextFindFaceStruct.vecFaceLocations.push_back(iEdge0);
					aNextFindFaceStruct.vecFaceLocations.push_back(iEdge1);

					aNextFindFaceStruct.loc = Face::NodeLocation_Edge;

					// Find the next Face index
					ixNextTargetFace =
						utils.FindFaceNearNode(
							meshTarget,
							nodevecSource[edgeSourceNext[0]],
							nodevecSource[edgeSourceNext[1]],
							edgeSourceNext.type,
							aNextFindFaceStruct);

#ifdef VERBOSE
					if (ixNextTargetFace == ixCurrentTargetFace) {
						printf("WARNING: Face does not change across Edge (3)\n");
					}
#endif
					// If face changes insert an Edge bifurcation
					if (ixNextTargetFace != ixCurrentTargetFace) {
						vecTracedPath.push_back(PathSegment(
							ixOverlapNodeCurrent,
							edgeSourceCurrent[1],
							edgeSourceCurrent.type,
							ixCurrentSourceFace,
							ixCurrentTargetFace,
							ixIntersectionTargetEdge,
							edgeTargetCurrent));
					}
				}

				// Remain on the same face
				if (ixNextTargetFace == ixCurrentTargetFace) {
					vecTracedPath.push_back(PathSegment(
						ixOverlapNodeCurrent,
						edgeSourceCurrent[1],
						edgeSourceCurrent.type,
						ixCurrentSourceFace,
						ixCurrentTargetFace,
						IntersectType_None,
						ixIntersectionTargetEdge));

				}

				// Update OverlapNodeCurrent
				ixOverlapNodeCurrent = edgeSourceCurrent[1];

				ixCurrentTargetFace = ixNextTargetFace;

				break;
			}

			// SourceEdge hits nodeTargetEdge0
			if (utils.AreNodesEqual(nodeIntersections[0], nodeTargetEdge0)) {
				ixOverlapNodeNext =
					vecTargetNodeMap[edgeTargetCurrent[0]];

				vecTracedPath.push_back(PathSegment(
					ixOverlapNodeCurrent,
					ixOverlapNodeNext,
					edgeSourceCurrent.type,
					ixCurrentSourceFace,
					ixCurrentTargetFace,
					IntersectType_Node,
					ixIntersectionTargetEdge));

				int ixPrevTargetFace = ixCurrentTargetFace;

				ixCurrentTargetFace =
					utils.FindFaceNearNode(
						meshTarget,
						edgeTargetCurrent[0],
						nodeSourceEnd,
						edgeSourceCurrent.type);
#ifdef VERBOSE
				if (ixPrevTargetFace == ixCurrentTargetFace) {
					printf("WARNING: Face does not change across Edge (4)\n");
				}
#endif
				ixOverlapNodeCurrent = ixOverlapNodeNext;

				if (ixOverlapNodeNext == edgeSourceCurrent[1]) {
					break;
				}

				continue;

			// SourceEdge hits nodeTargetEdge1
			} else if (
				utils.AreNodesEqual(nodeIntersections[0], nodeTargetEdge1)
			) {
				ixOverlapNodeNext =
					vecTargetNodeMap[edgeTargetCurrent[1]];

				vecTracedPath.push_back(PathSegment(
					ixOverlapNodeCurrent,
					ixOverlapNodeNext,
					edgeSourceCurrent.type,
					ixCurrentSourceFace,
					ixCurrentTargetFace,
					IntersectType_Node,
					(ixIntersectionTargetEdge + 1)
						% faceTargetCurrent.edges.size()));

				int ixPrevTargetFace = ixCurrentTargetFace;

				ixCurrentTargetFace =
					utils.FindFaceNearNode(
						meshTarget,
						edgeTargetCurrent[1],
						nodeSourceEnd,
						edgeSourceCurrent.type);

#ifdef VERBOSE
				if (ixPrevTargetFace == ixCurrentTargetFace) {
					printf("WARNING: Face does not change across Edge (5)\n");
				}
#endif
				ixOverlapNodeCurrent = ixOverlapNodeNext;

				if (ixOverlapNodeNext == edgeSourceCurrent[1]) {
					break;
				}

				continue;

			// General intersection between edgeSourceCurrent and
			// edgeTargetCurrent.
			} else {
				// Push a new intersection into the array of nodes
				ixOverlapNodeNext =
					static_cast<int>(meshOverlap.nodes.size());

				nodevecOverlap.push_back((Node)(nodeIntersections[0]));
/*
				NodeMap::const_iterator iterNodeMap =
					mapOverlapNodes.find((Node)(nodeIntersections[0]));

				if (iterNodeMap == mapOverlapNodes.end()) {
					meshOverlap.nodes.push_back((Node)(nodeIntersections[0]));
					mapOverlapNodes.insert(
						NodeMap::value_type(
							nodeIntersections[0], ixOverlapNodeNext));
				} else {
					ixOverlapNodeNext = iterNodeMap->second;
				}
*/
				// Intersection found with edge
				vecTracedPath.push_back(PathSegment(
					ixOverlapNodeCurrent,
					ixOverlapNodeNext,
					edgeSourceCurrent.type,
					ixCurrentSourceFace,
					ixCurrentTargetFace,
					ixIntersectionTargetEdge,
					edgeTargetCurrent));

				// Update OverlapNodeCurrent
				ixOverlapNodeCurrent = ixOverlapNodeNext;

				// Update TargetMesh face
				EdgeMapConstIterator iter =
					meshTarget.edgemap.find(edgeTargetCurrent);

				if (iter == meshTarget.edgemap.end()) {
					_EXCEPTIONT("Logic error");
				}

				const FacePair & facepair = iter->second;

				int ixPrevTargetFace = ixCurrentTargetFace;

				std::vector<int> vecPossibleFaces;
				vecPossibleFaces.push_back(facepair[0]);
				vecPossibleFaces.push_back(facepair[1]);

				// Face always changes across GreatCircleArcs
				if ((edgeSourceCurrent.type == Edge::Type_GreatCircleArc) &&
					(edgeTargetCurrent.type == Edge::Type_GreatCircleArc)
				) {
					if (ixPrevTargetFace == facepair[0]) {
						ixCurrentTargetFace = facepair[1];
					} else if (ixPrevTargetFace == facepair[1]) {
						ixCurrentTargetFace = facepair[0];
					} else {
						_EXCEPTIONT("Logic error");
					}

				} else {
					_EXCEPTIONT("Not implemented");
				}
/*
				const Face & face0 = meshTarget.faces[ixCurrentTargetFace];

				Face::NodeLocation loc0;
				int ixLocation0;
				utils.ContainsNode(
					face0,
					nodevecTarget,
					nodeIntersections[0],
					loc0,
					ixLocation0);
				
				printf("%i %i %i\n", ixCurrentTargetFace, loc0, ixLocation0);
*/
/*
				// Build the FindFaceStruct
				FindFaceStruct aNextFindFaceStruct;

				aNextFindFaceStruct.vecFaceIndices.push_back(facepair[0]);
				aNextFindFaceStruct.vecFaceIndices.push_back(facepair[1]);

				const Face & face0 = meshTarget.faces[facepair[0]];
				const Face & face1 = meshTarget.faces[facepair[1]];

				int iEdge0 = face0.GetEdgeIndex(edgeTargetCurrent);
				int iEdge1 = face1.GetEdgeIndex(edgeTargetCurrent);

				aNextFindFaceStruct.vecFaceLocations.push_back(iEdge0);
				aNextFindFaceStruct.vecFaceLocations.push_back(iEdge1);

				aNextFindFaceStruct.loc = Face::NodeLocation_Edge;

				// Find the next Face index near this Node
				ixCurrentTargetFace = 
					utils.FindFaceNearNode(
						meshTarget,
						nodeIntersections[0],
						meshOverlap.nodes[edgeSourceCurrent[1]],
						edgeSourceCurrent.type,
						aNextFindFaceStruct);
*/
#ifdef VERBOSE
				if (ixPrevTargetFace == ixCurrentTargetFace) {
					printf("WARNING: Face does not change across Edge (3)\n");
				}
#endif
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

bool GenerateOverlapFaces(
	const Mesh & meshTarget,
	const std::vector<int> & vecTargetNodeMap,
	const PathSegmentVector & vecTracedPath,
	int ixCurrentSourceFace,
	Mesh & meshOverlap
) {

	// Array indicating which elements of vecTracedPath have been used
	std::vector<bool> vecTracedPathUsed;
	vecTracedPathUsed.resize(vecTracedPath.size(), false);

	// The set of interior faces (from the Target mesh)
	std::set<int> setTargetFacesAdded;
	for (int j = 0; j < vecTracedPath.size(); j++) {
		setTargetFacesAdded.insert(vecTracedPath[j].ixTargetFace);

#ifdef VERBOSE
		printf("%i %i : %i\n", vecTracedPath[j][0], vecTracedPath[j][1], vecTracedPath[j].ixTargetFace);
#endif
	}

	// Set of faces from meshTarget that should be added
	std::set<int> setTargetFacesToAdd;

	// Loop through all possible starting PathSegments
	for (int j = 0; j < vecTracedPath.size(); j++) {

		// Ignore starting PathSegments that have already been used
		if (vecTracedPathUsed[j]) {
			continue;
		}

		// Build the new face
		Face faceOverlap(0);

		// Origin node of this fac
		int ixOverlapOriginNode = vecTracedPath[j][0];

		// Current second face
		int ixCurrentTargetFace =
			vecTracedPath[j].ixTargetFace;

		const Face & faceTargetCurrent =
			meshTarget.faces[ixCurrentTargetFace];

		// Search may require multiple trips along SourceMesh and TargetMesh
		for (;;) {

			// Loop along edge of SourceFace until we find a branch
			// into TargetMesh or hit our origin.
			for (;; j++) {
				if (j >= vecTracedPath.size()) {
					printf("ERROR: Invalid TracedPathEdge index\n");
					return false;
				}

				faceOverlap.edges.push_back(vecTracedPath[j]);

				// Mark traced path edge as used
				if (vecTracedPathUsed[j]) {
					_EXCEPTIONT("Trying to reuse traced path edge");
				}

				vecTracedPathUsed[j] = true;
#ifdef VERBOSE
				printf("P%i: %i %i\n",
					j, vecTracedPath[j][0], vecTracedPath[j][1]);
#endif
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
			int ixCurrentTargetEdge = vecTracedPath[j].ixIntersect;

			int ixCurrentOverlapNode = vecTracedPath[j][1];

			// Loop around the interior of faceTargetCurrent
			int nEdgesCompleted = 0;
			for (;;) {

				const Edge & edgeTargetCurrent =
					faceTargetCurrent.edges[ixCurrentTargetEdge];

				// Check for infinite loop
				if (nEdgesCompleted > faceTargetCurrent.edges.size()) {
					printf("Possible infinite loop - aborting\n");
					return false;
				}

				nEdgesCompleted++;

				// Identical endpoints; advance the edge
				if (edgeTargetCurrent[0] == edgeTargetCurrent[1]) {
					ixCurrentTargetEdge =
						(ixCurrentTargetEdge + 1)
							% faceTargetCurrent.edges.size();

					ixCurrentOverlapNode = //edgeTargetCurrent[1];
						vecTargetNodeMap[edgeTargetCurrent[1]];

					continue;
				}

				// Determine if this edge exits onto SourceMesh
				int ixExitNode = InvalidNode;

				int k;
				for (
					k = (j + 1) % vecTracedPath.size();
					k != j;
					k = (k + 1) % vecTracedPath.size()
				) {
					int ixIntersectNode =
						faceTargetCurrent[vecTracedPath[k].ixIntersect];

					IntersectType inttype = vecTracedPath[k].inttype;

					if (ixCurrentOverlapNode == vecTracedPath[k][1]) {
						continue;
					}

					// Check for node intersections
					if (inttype == IntersectType_Node) {
						if (vecTracedPath[k][1] == //edgeTargetCurrent[0]
							vecTargetNodeMap[edgeTargetCurrent[0]]
						) {
							ixExitNode = vecTracedPath[k][1];
							break;
						}

						if (vecTracedPath[k][1] == //edgeTargetCurrent[1]
							vecTargetNodeMap[edgeTargetCurrent[1]]
						) {
							ixExitNode = vecTracedPath[k][1];
							break;
						}
					}

					// Check for edge intersections
					if ((inttype == IntersectType_Edge) &&
						(edgeTargetCurrent == vecTracedPath[k].edgeIntersect)
					) {
						ixExitNode = vecTracedPath[k][1];
						break;
					}
				}

				// Add the interior face to the list of faces to be added
				EdgeMapConstIterator iter =
					meshTarget.edgemap.find(edgeTargetCurrent);

				if (iter == meshTarget.edgemap.end()) {
					_EXCEPTIONT("Logic error");
				}

				const FacePair & facepair = iter->second;

				if (facepair[0] == ixCurrentTargetFace) {
					setTargetFacesToAdd.insert(facepair[1]);
				} else if (facepair[1] == ixCurrentTargetFace) {
					setTargetFacesToAdd.insert(facepair[0]);
				} else {
					_EXCEPTIONT("Logic error");
				}

				// An exit node is found
				if (ixExitNode != InvalidNode) {

					// Check if the traced path becomes active
					int jNext = (k + 1) % vecTracedPath.size();
					if (vecTracedPath[jNext].ixTargetFace ==
						ixCurrentTargetFace
					) {
#ifdef VERBOSE
						printf("S: %i %i\n",
							ixCurrentOverlapNode, ixExitNode);
#endif
						faceOverlap.edges.push_back(Edge(
							ixCurrentOverlapNode,
							ixExitNode,
							edgeTargetCurrent.type));

						j = jNext;
						if (ixExitNode == ixOverlapOriginNode) {
							goto ContinueToNextFace;
						} else {
							break;
						}
					}
				}

#ifdef VERBOSE
				printf("T: %i (%i) %i\n",
					ixCurrentOverlapNode,
					//edgeTargetCurrent[0],
					//edgeTargetCurrent[1]);
					vecTargetNodeMap[edgeTargetCurrent[0]],
					vecTargetNodeMap[edgeTargetCurrent[1]]);
#endif
				// Push this edge into the overlap mesh
				faceOverlap.edges.push_back(Edge(
					ixCurrentOverlapNode,
					//edgeTargetCurrent[1],
					vecTargetNodeMap[edgeTargetCurrent[1]],
					edgeTargetCurrent.type));

				// Advance the edge
				ixCurrentTargetEdge =
					(ixCurrentTargetEdge + 1)
						% faceTargetCurrent.edges.size();

				ixCurrentOverlapNode =
					//edgeTargetCurrent[1];
					vecTargetNodeMap[edgeTargetCurrent[1]];

				if (ixCurrentOverlapNode == ixOverlapOriginNode) {
					goto ContinueToNextFace;
				}
			}
		}

ContinueToNextFace:
		if (faceOverlap.edges.size() < 3) {
			printf("ERROR: Overlap Face containing only two Nodes detected\n");
			return false;
		}

#ifdef VERBOSE
		printf("PUSH %lu\n", faceOverlap.edges.size());
/*
		printf("AREA %1.15e\n", CalculateFaceArea(faceOverlap, meshOverlap.nodes));

		for (int k = 0; k < faceOverlap.edges.size(); k++) {
			meshOverlap.nodes[faceOverlap[k]].Print("n");
		}
*/
#endif

		// Push this Face into the overlap Mesh
		if (CalculateFaceArea(faceOverlap, meshOverlap.nodes) >= 1.0e-12) {
			meshOverlap.faces.push_back(faceOverlap);
			meshOverlap.vecSourceFaceIx.push_back(ixCurrentSourceFace);
			meshOverlap.vecTargetFaceIx.push_back(ixCurrentTargetFace);
		}

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

	///////////////////////////////////////////////////////////////////////
	// Find interior faces from meshTarget to add to meshOverlap

	// Remove added faces from set of faces to add
	std::set<int>::iterator iterToAdd;
	std::set<int>::iterator iterAdded = setTargetFacesAdded.begin();
	for (; iterAdded != setTargetFacesAdded.end(); iterAdded++) {
		iterToAdd = setTargetFacesToAdd.find(*iterAdded);

		if (iterToAdd != setTargetFacesToAdd.end()) {
			setTargetFacesToAdd.erase(iterToAdd);
		}
	}

	for (;;) {

		// Take the next Face from the ToAdd set
		if (setTargetFacesToAdd.size() == 0) {
			break;
		}

		iterToAdd = setTargetFacesToAdd.begin();

		int ixFaceToAdd = *iterToAdd;

		setTargetFacesToAdd.erase(iterToAdd);
		setTargetFacesAdded.insert(ixFaceToAdd);

		const Face & faceTargetCurrent = meshTarget.faces[ixFaceToAdd];

		// Add this face to meshOverlap
		Face faceOverlapCurrent(faceTargetCurrent.edges.size());

		for (int i = 0; i < faceOverlapCurrent.edges.size(); i++) {
			faceOverlapCurrent.edges[i][0] =
				//faceTargetCurrent.edges[i][0];
				vecTargetNodeMap[faceTargetCurrent.edges[i][0]];
			faceOverlapCurrent.edges[i][1] =
				//faceTargetCurrent.edges[i][1];
				vecTargetNodeMap[faceTargetCurrent.edges[i][1]];
			faceOverlapCurrent.edges[i].type =
				faceTargetCurrent.edges[i].type;
		}

		if (CalculateFaceArea(faceOverlapCurrent, meshOverlap.nodes) >= 1.0e-12) {
			meshOverlap.faces.push_back(faceOverlapCurrent);
			meshOverlap.vecSourceFaceIx.push_back(ixCurrentSourceFace);
			meshOverlap.vecTargetFaceIx.push_back(ixFaceToAdd);
		}

		// Add further interior faces that are Edge-neighbors of this Face
		bool fMoreFacesToAdd = false;

		for (int i = 0; i < faceTargetCurrent.edges.size(); i++) {

			if (faceTargetCurrent.edges[i][0] ==
				faceTargetCurrent.edges[i][1]
			) {
				continue;
			}

			EdgeMapConstIterator iterEdge =
				meshTarget.edgemap.find(faceTargetCurrent.edges[i]);

			if (iterEdge == meshTarget.edgemap.end()) {
				_EXCEPTIONT("Edge not found in EdgeMap");
			}

			int ixOtherFace;
			if (iterEdge->second[0] == ixFaceToAdd) {
				ixOtherFace = iterEdge->second[1];
			} else if (iterEdge->second[1] == ixFaceToAdd) {
				ixOtherFace = iterEdge->second[0];
			} else {
				_EXCEPTIONT("EdgeMap consistency error");
			}

			if (setTargetFacesAdded.find(ixOtherFace) ==
				setTargetFacesAdded.end()
			) {
				setTargetFacesToAdd.insert(ixOtherFace);
				fMoreFacesToAdd = true;
			}
		}
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////

void GenerateOverlapMesh_v1(
	const Mesh & meshSource,
	const Mesh & meshTarget,
	Mesh & meshOverlap,
	OverlapMeshMethod method
) {
	meshOverlap.Clear();

	// Get the two NodeVectors
	const NodeVector & nodevecSource = meshSource.nodes;
	const NodeVector & nodevecTarget = meshTarget.nodes;

	// Check for coincident Nodes between meshSource and meshTarget
	std::vector<int> vecTargetNodeMap;

	int nCoincidentNodes =
		BuildCoincidentNodeVector(
			meshSource, meshTarget, vecTargetNodeMap);

	Announce("Number of coincident nodes between mesh A and B [%i]",
		nCoincidentNodes);

	// Insert all nodes from the two NodeVectors
	for (int i = 0; i < nodevecSource.size(); i++) {
		meshOverlap.nodes.push_back(nodevecSource[i]);
		//mapOverlapNodes.insert(NodeMap::value_type(nodevecSource[i], i));
	}

	for (int i = 0; i < nodevecTarget.size(); i++) {
		if (vecTargetNodeMap[i] == InvalidNode) {
			int ix = static_cast<int>(meshOverlap.nodes.size());
			meshOverlap.nodes.push_back(nodevecTarget[i]);
			//mapOverlapNodes.insert(NodeMap::value_type(nodevecTarget[i], ix));
			vecTargetNodeMap[i] = ix;
		}
	}

/*
	// Estimate meshOverlap size
	int nMaximumFaceCount;
	if (meshSource.faces.size() > meshTarget.faces.size()) {
		nMaximumFaceCount = static_cast<int>(meshSource.faces.size());
	} else {
		nMaximumFaceCount = static_cast<int>(meshTarget.faces.size());
	}
	meshOverlap.faces.reserve(2 * nMaximumFaceCount);
*/
	// Loop through all Faces on the first Mesh
	int ixCurrentSourceFace = 0;

#pragma message "OpenMP here"
	for (; ixCurrentSourceFace < meshSource.faces.size(); ixCurrentSourceFace++) {
	//for (int ixCurrentSourceFace = 853; ixCurrentSourceFace < 854; ixCurrentSourceFace++) {

#ifdef CHECK_AREAS
		Real dSourceFaceArea = meshSource.CalculateFaceArea(ixCurrentSourceFace);
#endif

		// Generate the path
		PathSegmentVector vecTracedPath;

		// Fuzzy arithmetic (standard floating point operations)
		if (method == OverlapMeshMethod_Fuzzy) {
			GeneratePath<MeshUtilitiesFuzzy, Node>(
				meshSource,
				meshTarget,
				vecTargetNodeMap,
				ixCurrentSourceFace,
				vecTracedPath,
				meshOverlap
			);

			GenerateOverlapFaces(
				meshTarget,
				vecTargetNodeMap,
				vecTracedPath,
				ixCurrentSourceFace,
				meshOverlap
			);
		}

		// Exact arithmetic
		if (method == OverlapMeshMethod_Exact) {
			GeneratePath<MeshUtilitiesExact, NodeExact>(
				meshSource,
				meshTarget,
				vecTargetNodeMap,
				ixCurrentSourceFace,
				vecTracedPath,
				meshOverlap
			);

			GenerateOverlapFaces(
				meshTarget,
				vecTargetNodeMap,
				vecTracedPath,
				ixCurrentSourceFace,
				meshOverlap
			);
		}

		// Mixed method; try Fuzzy arithmetic first
		if (method == OverlapMeshMethod_Mixed) {
			int nInitialOverlapNodes = meshOverlap.nodes.size();
			int nInitialOverlapFaces = meshOverlap.faces.size();

			try {
				GeneratePath<MeshUtilitiesFuzzy, Node>(
					meshSource,
					meshTarget,
					vecTargetNodeMap,
					ixCurrentSourceFace,
					vecTracedPath,
					meshOverlap
				);

				GenerateOverlapFaces(
					meshTarget,
					vecTargetNodeMap,
					vecTracedPath,
					ixCurrentSourceFace,
					meshOverlap
				);

			} catch(Exception & e) {
				printf("WARNING: Fuzzy arithmetic operations failed "
					"with message:\n  \"%s\"\n  Trying exact arithmetic",
					e.ToString().c_str());

				vecTracedPath.clear();

				meshOverlap.nodes.resize(nInitialOverlapNodes);
				meshOverlap.faces.resize(nInitialOverlapFaces);
				meshOverlap.vecSourceFaceIx.resize(nInitialOverlapFaces);
				meshOverlap.vecTargetFaceIx.resize(nInitialOverlapFaces);

				GeneratePath<MeshUtilitiesExact, NodeExact>(
					meshSource,
					meshTarget,
					vecTargetNodeMap,
					ixCurrentSourceFace,
					vecTracedPath,
					meshOverlap
				);

				GenerateOverlapFaces(
					meshTarget,
					vecTargetNodeMap,
					vecTracedPath,
					ixCurrentSourceFace,
					meshOverlap
				);
			}
		}

#ifdef CHECK_AREAS
		int nMeshOverlapPrevFaces = meshOverlap.faces.size();
#endif

#ifdef CHECK_AREAS
		int nMeshOverlapCurrentFaces = meshOverlap.faces.size();

		Real dOverlapAreas = 0.0;
		for (int i = nMeshOverlapPrevFaces; i < nMeshOverlapCurrentFaces; i++) {
			dOverlapAreas += meshOverlap.CalculateFaceArea(i);
		}

		if (fabs(dOverlapAreas - dSourceFaceArea) > ReferenceTolerance) {
			printf("Area inconsistency (%i : %1.15e %1.15e)\n",
				ixCurrentSourceFace,
				dSourceFaceArea,
				dOverlapAreas);
			_EXCEPTION();
		}
#endif

/*
		if (!fSuccess) {
			_EXCEPTIONT("OverlapMesh generation failed");
		}
*/
	}
}

///////////////////////////////////////////////////////////////////////////////

template <
	class MeshUtilities,
	class NodeIntersectType
>
void GenerateOverlapFace(
	const Mesh & meshSource,
	const Mesh & meshTarget,
	int iSourceFace,
	int iTargetFace,
	NodeVector & nodevecOutput
) {
/*
  // Sutherland–Hodgman algorithm (pseudocode)
  List outputList = subjectPolygon;
  for (Edge clipEdge in clipPolygon) do
     List inputList = outputList;
     outputList.clear();
     Point S = inputList.last;
     for (Point E in inputList) do
        if (E inside clipEdge) then
           if (S not inside clipEdge) then
              outputList.add(ComputeIntersection(S,E,clipEdge));
           end if
           outputList.add(E);
        else if (S inside clipEdge) then
           outputList.add(ComputeIntersection(S,E,clipEdge));
        end if
        S = E;
     done
  done
*/
	MeshUtilities utils;

	// Use Sutherland–Hodgman algorithm to clip polygons
	const NodeVector & nodesTarget = meshTarget.nodes;
	const NodeVector & nodesSource = meshSource.nodes;

	const Face & faceTarget = meshTarget.faces[iTargetFace];
	const Face & faceSource = meshSource.faces[iSourceFace];

	const EdgeVector & evecTarget = faceTarget.edges;
	const EdgeVector & evecSource = faceSource.edges;

	// List outputList = subjectPolygon
	for (int i = 0; i < evecTarget.size(); i++) {
		nodevecOutput.push_back(nodesTarget[evecTarget[i][0]]);
	}

	//std::cout << nodevecOutput.size() << std::endl;

	// for (Edge clipEdge in clipPolygon) do
	for (int i = 0; i < evecSource.size(); i++) {

		// All points outside of source polygon
		if (nodevecOutput.size() == 0) {
			break;
		}

		// List inputList = outputList;
		NodeVector nodevecInput = nodevecOutput;

		// outputList.clear();
		nodevecOutput.clear();

		// Point S = inputList.last;
		Node nodeS = nodevecInput[nodevecInput.size()-1];

		int iNodeEdgeSideS =
			utils.FindNodeEdgeSide(
				nodesSource[evecSource[i][0]],
				nodesSource[evecSource[i][1]],
				evecSource[i].type,
				nodeS);

		//printf("===================\n");
		//printf("S Side: %i\n", iNodeEdgeSideS);

		// for (Point E in inputList) do
		for (int iNodeE = 0; iNodeE < nodevecInput.size(); iNodeE++) {
			const Node & nodeE = nodevecInput[iNodeE];

			// if (E inside clipEdge) then
			int iNodeEdgeSideE =
				utils.FindNodeEdgeSide(
					nodesSource[evecSource[i][0]],
					nodesSource[evecSource[i][1]],
					evecSource[i].type,
					nodeE);

			//printf("E Side: %i\n", iNodeEdgeSideE);

			if (iNodeEdgeSideE >= 0) {

				// if (S not inside clipEdge) then
				if (iNodeEdgeSideS < 0) {

					// outputList.add(ComputeIntersection(S,E,clipEdge));
					std::vector<Node> vecIntersections;
					bool fCoincident =
						utils.CalculateEdgeIntersectionsSemiClip(
							nodeS,
							nodeE,
							Edge::Type_GreatCircleArc,
							nodesSource[evecSource[i][0]],
							nodesSource[evecSource[i][1]],
							evecSource[i].type,
							vecIntersections);

/*
					vecIntersections[0].Print("X");

					if (vecIntersections.size() != 1) {
						printf("%i %i\n", iNodeEdgeSideE, iNodeEdgeSideS);
						std::cout << vecIntersections.size() << std::endl;

						nodesSource[evecSource[i][0]].Print("S0");
						nodesSource[evecSource[i][1]].Print("S1");
						nodeS.Print("S");
						nodeE.Print("E");
						_EXCEPTIONT("Logic error");
					}
*/
					if (vecIntersections.size() != 1) {
						_EXCEPTIONT("Logic error");
					}

					//vecIntersections[0].Print("Add");
					nodevecOutput.push_back(vecIntersections[0]);
				}

				// outputList.add(E);
				nodevecOutput.push_back(nodeE);

			// else if (S inside clipEdge) then
			} else if (iNodeEdgeSideS >= 0) {

				// outputList.add(ComputeIntersection(S,E,clipEdge));
				std::vector<Node> vecIntersections;
				bool fCoincident =
					utils.CalculateEdgeIntersectionsSemiClip(
						nodeS,
						nodeE,
						Edge::Type_GreatCircleArc,
						nodesSource[evecSource[i][0]],
						nodesSource[evecSource[i][1]],
						evecSource[i].type,
						vecIntersections);
/*
				nodesSource[evecSource[i][0]].Print("S0");
				nodesSource[evecSource[i][1]].Print("S1");
				nodeS.Print("S");
				nodeE.Print("E");
				vecIntersections[0].Print("X");

				if (vecIntersections.size() != 1) {
					_EXCEPTIONT("Logic error");
				}

				nodevecOutput.push_back(vecIntersections[0]);
*/
				if (fCoincident) {
					//nodevecOutput.push_back(nodeE);
				} else if (vecIntersections.size() == 1) {
					nodevecOutput.push_back(vecIntersections[0]);

				} else {
					_EXCEPTIONT("Logic error");
				}

				//vecIntersections[0].Print("Add");
			}

			// S = E;
			nodeS = nodeE;
			iNodeEdgeSideS = iNodeEdgeSideE;
		}
/*
		std::cout << nodevecOutput.size() << std::endl;
		for (int i = 0; i < nodevecOutput.size(); i++) {
			nodevecOutput[i].Print("");
		}
*/
	}

	//std::cout << nodevecOutput.size() << std::endl;
/*
	// Check for repeats
	for (int i = 0; i < nodevecOutput.size(); i++) {
		if (nodevecOutput[i] == nodevecOutput[(i+1)%(nodevecOutput.size())]) {
			nodevecOutput.erase(nodevecOutput.begin()+i);
			i--;
		}
 	}
*/
	// If the overlap consists of fewer than three nodes, ignore
	if (nodevecOutput.size() < 3) {
		nodevecOutput.clear();
	}
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Find the Face in meshTarget containing node, searching near
///		meshTarget Face with index ixTargetFaceSeed.
///	</summary>
template <
	class MeshUtilities,
	class NodeIntersectType
>
int FindFaceContainingNode(
	const Mesh & meshSource,
	const Mesh & meshTarget,
	int ixSourceFaceSeed,
	int ixTargetFaceSeed
) {
	if (ixSourceFaceSeed > meshSource.faces.size()) {
		_EXCEPTIONT("SourceFaceSeed greater than Source mesh size");
	}
	if (ixTargetFaceSeed > meshTarget.faces.size()) {
		_EXCEPTIONT("TargetFaceSeed greater than Target mesh size");
	}

	const Node & node =
		meshSource.nodes[meshSource.faces[ixSourceFaceSeed][0]];

	const EdgeMap & edgemapTarget = meshTarget.edgemap;

	MeshUtilities utils;

	Face::NodeLocation loc;
	int ixLocation;

	std::set<int> setExaminedTargetFaces;
	std::queue<int> queueTargetFaces;

	queueTargetFaces.push(ixTargetFaceSeed);
	setExaminedTargetFaces.insert(ixTargetFaceSeed);

	while (!queueTargetFaces.empty()) {

		int ixCurrentTargetFace = queueTargetFaces.front();
		queueTargetFaces.pop();

		const Face & faceTarget = meshTarget.faces[ixCurrentTargetFace];

		utils.ContainsNode(
			faceTarget,
			meshTarget.nodes,
			node,
			loc,
			ixLocation);

		// Node found within interior of target face; overlap detected
		if (loc == Face::NodeLocation_Interior) {
			return ixCurrentTargetFace;
		}

		// Node on boundary of target face; check for overlap
		if (loc != Face::NodeLocation_Exterior) {
			NodeVector nodevecOverlap;

			GenerateOverlapFace<MeshUtilitiesFuzzy, Node>(
				meshSource,
				meshTarget,
				ixSourceFaceSeed,
				ixCurrentTargetFace,
				nodevecOverlap
			);

			if (nodevecOverlap.size() > 2) {
				return ixCurrentTargetFace;
			}
		}

		// Add all neighboring Faces into the queue of Target Faces
		for (int i = 0; i < faceTarget.edges.size(); i++) {
			EdgeMapConstIterator iter =
				edgemapTarget.find(faceTarget.edges[i]);

			if (iter == edgemapTarget.end()) {
				_EXCEPTIONT("Missing Edge in Target EdgeMap");
			}

			const FacePair & facepair = iter->second;

			int iPushFace;
			if (facepair[0] == ixCurrentTargetFace) {
				iPushFace = facepair[1];

			} else if (facepair[1] == ixCurrentTargetFace) {
				iPushFace = facepair[0];

			} else {
				_EXCEPTIONT("EdgeMap error");
			}

			std::set<int>::const_iterator iterFace =
				setExaminedTargetFaces.find(iPushFace);

			if (iterFace == setExaminedTargetFaces.end()) {
				queueTargetFaces.push(iPushFace);
				setExaminedTargetFaces.insert(iPushFace);
			}
		}
	}

	_EXCEPTIONT("Unable to find Target Face");
}

///////////////////////////////////////////////////////////////////////////////

void GenerateOverlapMeshFromFace(
	const Mesh & meshSource,
	const Mesh & meshTarget,
	int ixSourceFace,
	Mesh & meshOverlap,
	NodeMap & nodemapOverlap,
	OverlapMeshMethod method,
	int ixTargetFaceSeed = 0
) {
	// Verify the EdgeMap exists in both meshSource and meshTarget
	if (meshSource.edgemap.size() == 0) {
		_EXCEPTIONT("EdgeMap in meshSource must be constructed prior"
			" to GenerateOverlapFace");
	}
	if (meshTarget.edgemap.size() == 0) {
		_EXCEPTIONT("EdgeMap in meshTarget must be constructed prior"
			" to GenerateOverlapFace");
	}

	// Set of Faces on the Target Mesh that overlap ixSourceFace
	std::set<int> setExaminedTargetFaces;

	// Get the two NodeVectors
	const NodeVector & nodevecSource = meshSource.nodes;
	const NodeVector & nodevecTarget = meshTarget.nodes;

	// EdgeMap
	const EdgeMap & edgemapTarget = meshTarget.edgemap;

	// Find the starting face on the target mesh
	MeshUtilitiesFuzzy utils;
/*
	FindFaceStruct aFindFaceStruct;
	utils.FindFaceFromNode(
		meshTarget,
		meshSource.nodes[meshSource.faces[ixSourceFace][0]],
		aFindFaceStruct);

	if (aFindFaceStruct.vecFaceIndices.size() < 1) {
		_EXCEPTIONT("No target face found");
	}
*/
	int ixCurrentTargetFace =
		FindFaceContainingNode<MeshUtilitiesFuzzy, Node>(
			meshSource,
			meshTarget,
			ixSourceFace,
			ixTargetFaceSeed);

	std::cout << "\tFirst overlap match " << ixCurrentTargetFace << std::endl;
/*
	// Verify starting Node is not on the Exterior
	if (aFindFaceStruct.loc == Face::NodeLocation_Exterior) {
		_EXCEPTIONT("Invalid Node location returned from FindFaceFromNode");
	}

	// Current target Face
	int ixCurrentTargetFace = aFindFaceStruct.vecFaceIndices[0];
*/
	std::queue<int> queueTargetFaces;
	queueTargetFaces.push(ixCurrentTargetFace);
	setExaminedTargetFaces.insert(ixCurrentTargetFace);

	while (!queueTargetFaces.empty()) {
/*
	meshSource.nodes[meshSource.faces[ixSourceFace][0]].Print("S");
	meshTarget.nodes[meshTarget.faces[ixCurrentTargetFace][0]].Print("T");

	for (int i = 0; i < meshSource.faces[ixSourceFace].edges.size(); i++) {
		meshSource.nodes[meshSource.faces[ixSourceFace][i]].Print("");
	}
	printf("\n");
	for (int i = 0; i < meshTarget.faces[ixCurrentTargetFace].edges.size(); i++) {
		meshTarget.nodes[meshTarget.faces[ixCurrentTargetFace][i]].Print("");
	}
*/
		// Get the next target face
		ixCurrentTargetFace = queueTargetFaces.front();
		queueTargetFaces.pop();

		const Face & faceTarget = meshTarget.faces[ixCurrentTargetFace];

		// Find the overlap polygon
		NodeVector nodevecOutput;

		GenerateOverlapFace<MeshUtilitiesFuzzy, Node>(
			meshSource,
			meshTarget,
			ixSourceFace,
			ixCurrentTargetFace,
			nodevecOutput
		);

		if (nodevecOutput.size() == 0) {

		} else if (nodevecOutput.size() < 3) {
			_EXCEPTIONT("Overlap polygon consists of fewer than 3 nodes");

		} else {

			// Add all neighboring Faces into the queue of Target Faces
			for (int i = 0; i < faceTarget.edges.size(); i++) {
				EdgeMapConstIterator iter =
					edgemapTarget.find(faceTarget.edges[i]);

				if (iter == edgemapTarget.end()) {
					_EXCEPTIONT("Missing Edge in Target EdgeMap");
				}

				const FacePair & facepair = iter->second;

				int iPushFace;
				if (facepair[0] == ixCurrentTargetFace) {
					iPushFace = facepair[1];

				} else if (facepair[1] == ixCurrentTargetFace) {
					iPushFace = facepair[0];

				} else {
					_EXCEPTIONT("EdgeMap error");
				}

				std::set<int>::const_iterator iterFace =
					setExaminedTargetFaces.find(iPushFace);

				if (iterFace == setExaminedTargetFaces.end()) {
					queueTargetFaces.push(iPushFace);
					setExaminedTargetFaces.insert(iPushFace);
				}
			}

			std::cout << "\tOverlap with Face " << ixCurrentTargetFace
				<< std::endl;

			// Calculate face area
			Face faceTemp(nodevecOutput.size());
			for (int i = 0; i < nodevecOutput.size(); i++) {
				faceTemp.SetNode(i, i);
			}
			Real dArea = CalculateFaceArea(faceTemp, nodevecOutput);

			if (dArea < 1.0e-13) {
				continue;
/*
				printf("%i %i %1.15e %i\n", ixSourceFace, ixCurrentTargetFace, dArea, nodevecOutput.size());
				for (int i = 0; i < nodevecOutput.size(); i++) {
					nodevecOutput[i].Print("");
				}
				//_EXCEPTION();
*/
			}

			// Insert Face into meshOverlap
			Face faceNew(nodevecOutput.size());
			for (int i = 0; i < nodevecOutput.size(); i++) {
				NodeMapConstIterator iter
					= nodemapOverlap.find(nodevecOutput[i]);

				if (iter != nodemapOverlap.end()) {
					faceNew.SetNode(i, iter->second);
				} else {
					int iNextNodeMapOverlapIx = nodemapOverlap.size();
					faceNew.SetNode(i, iNextNodeMapOverlapIx);
					nodemapOverlap.insert(
						NodeMapPair(nodevecOutput[i], iNextNodeMapOverlapIx));
				}
			}
			meshOverlap.faces.push_back(faceNew);

			meshOverlap.vecSourceFaceIx.push_back(ixSourceFace);
			meshOverlap.vecTargetFaceIx.push_back(ixCurrentTargetFace);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void GenerateOverlapMesh_v2(
	const Mesh & meshSource,
	const Mesh & meshTarget,
	Mesh & meshOverlap,
	OverlapMeshMethod method
) {
	NodeMap nodemapOverlap;

	// Create a KD tree for fast searching of the target mesh
	kdtree * kdTarget = kd_create(3);

	for (int i = 0; i < meshTarget.faces.size(); i++) {
		int ixNodeCorner = meshTarget.faces[i][0];
		kd_insert3(
			kdTarget,
			meshTarget.nodes[ixNodeCorner].x,
			meshTarget.nodes[ixNodeCorner].y,
			meshTarget.nodes[ixNodeCorner].z,
			(void*)(((char*)(nullptr)) + i));
	}

	// Generate Overlap mesh for each Face
	for (int i = 0; i < meshSource.faces.size(); i++) {
		std::cout << "Source Face " << i << std::endl;

		// Find a Target face near this source face
		int ixNodeCorner = meshSource.faces[i][0];

		kdres * kdresTarget =
			kd_nearest3(
				kdTarget,
				meshSource.nodes[ixNodeCorner].x,
				meshSource.nodes[ixNodeCorner].y,
				meshSource.nodes[ixNodeCorner].z);

		void * pdata = kd_res_item_data(kdresTarget);

		int iTargetFaceSeed = ((char*)(pdata)) - nullptr;

		std::cout << "\tNearest target face " << iTargetFaceSeed << std::endl;

		// Generate the overlap mesh associated with this source face
		GenerateOverlapMeshFromFace(
			meshSource,
			meshTarget,
			i,
			meshOverlap,
			nodemapOverlap,
			method,
			iTargetFaceSeed);
	}

	// Destroy the KD tree
	kd_free(kdTarget);

	// Insert all Nodes from nodemapOverlap into meshOverlap.nodes
	meshOverlap.nodes.resize(nodemapOverlap.size());

	NodeMapConstIterator iter = nodemapOverlap.begin();
	for (; iter != nodemapOverlap.end(); iter++) {
		meshOverlap.nodes[iter->second] = iter->first;
	}

	// Calculate Face areas
	double dTotalAreaOverlap = meshOverlap.CalculateFaceAreas();
	Announce("Overlap Mesh Geometric Area: %1.15e (%1.15e)", dTotalAreaOverlap, 4.0 * M_PI);

}

///////////////////////////////////////////////////////////////////////////////
