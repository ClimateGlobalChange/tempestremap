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

#include "Announce.h"

#include <unistd.h>
#include <iostream>

///////////////////////////////////////////////////////////////////////////////

#define VERBOSE

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

///	<summary>
///		Location data returned from FindFaceFromNode()
///	</summary>
struct FindFaceStruct {
	
	///	<summary>
	///		A vector of face indices indicating possible Faces.
	///	</summary>
	std::vector<int> vecFaceIndices;

	///	<summary>
	///		A vector of locations on each Face.  If loc is NodeLocation_Corner,
	///		this corresponds to the associated corner of the Face.  If loc
	///		is NodeLocation_Edge, this corresponds to the associated Edge of
	///		the Face.  If loc is NodeLocation_Interior, this value is
	///		undefined.
	///	</summary>
	std::vector<int> vecFaceLocations;

	///	<summary>
	///		The NodeLocation where this Node lies.
	///	</summary>
	Face::NodeLocation loc;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Find all face indices that contain this node.
///	</summary>
void FindFaceFromNode(
	const Mesh & mesh,
	const Node & node,
	FindFaceStruct & aFindFaceStruct
) {
#ifdef VERBOSE
	printf("---\n");
#endif
	// Reset the FaceStruct
	aFindFaceStruct.vecFaceIndices.clear();
	aFindFaceStruct.vecFaceLocations.clear();
	aFindFaceStruct.loc = Face::NodeLocation_Undefined;

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

		if (loc == Face::NodeLocation_Exterior) {
			continue;
		}

#ifdef VERBOSE
		printf("%i\n", l);
		printf("n: %1.5e %1.5e %1.5e\n", node.x, node.y, node.z);
		printf("n0: %1.5e %1.5e %1.5e\n",
			mesh.nodes[mesh.faces[l][0]].x,
			mesh.nodes[mesh.faces[l][0]].y,
			mesh.nodes[mesh.faces[l][0]].z);
		printf("n1: %1.5e %1.5e %1.5e\n",
			mesh.nodes[mesh.faces[l][1]].x,
			mesh.nodes[mesh.faces[l][1]].y,
			mesh.nodes[mesh.faces[l][1]].z);
		printf("n2: %1.5e %1.5e %1.5e\n",
			mesh.nodes[mesh.faces[l][2]].x,
			mesh.nodes[mesh.faces[l][2]].y,
			mesh.nodes[mesh.faces[l][2]].z);
		printf("n3: %1.5e %1.5e %1.5e\n",
			mesh.nodes[mesh.faces[l][3]].x,
			mesh.nodes[mesh.faces[l][3]].y,
			mesh.nodes[mesh.faces[l][3]].z);
#endif

		if (aFindFaceStruct.loc == Face::NodeLocation_Undefined) {
			aFindFaceStruct.loc = loc;
		}

		// Node is in the interior of this face
		if (loc == Face::NodeLocation_Interior) {
			if (loc != aFindFaceStruct.loc) {
				_EXCEPTIONT("No consensus on location of Node");
			}

			aFindFaceStruct.vecFaceIndices.push_back(l);
			aFindFaceStruct.vecFaceLocations.push_back(ixLocation);
			break;
		}

		// Node is on the edge of this face
		if (loc == Face::NodeLocation_Edge) {
			if (loc != aFindFaceStruct.loc) {
				_EXCEPTIONT("No consensus on location of Node");
			}

			aFindFaceStruct.vecFaceIndices.push_back(l);
			aFindFaceStruct.vecFaceLocations.push_back(ixLocation);
		}

		// Node is at the corner of this face
		if (loc == Face::NodeLocation_Corner) {
			if (loc != aFindFaceStruct.loc) {
				_EXCEPTIONT("No consensus on location of Node");
			}

			aFindFaceStruct.vecFaceIndices.push_back(l);
			aFindFaceStruct.vecFaceLocations.push_back(ixLocation);
		}
	}

	// Edges can only have two adjacent Faces
	if (aFindFaceStruct.loc == Face::NodeLocation_Edge) {
		if (aFindFaceStruct.vecFaceIndices.size() != 2) {
			printf("n: %1.5e %1.5e %1.5e\n", node.x, node.y, node.z);
			_EXCEPTION1("Multiple co-located edges detected (%i)",
				(int)(aFindFaceStruct.vecFaceIndices.size()));
		}
	}

	// Corners must have at least three adjacent Faces
	if (aFindFaceStruct.loc == Face::NodeLocation_Corner) {
		if (aFindFaceStruct.vecFaceIndices.size() < 3) {
			printf("n: %1.5e %1.5e %1.5e\n", node.x, node.y, node.z);
			_EXCEPTION1("Two Faced corner detected (%i)",
				(int)(aFindFaceStruct.vecFaceIndices.size()));
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Find all face indices that contain this node from a list of
///		possible Faces.
///	</summary>
void FindFaceFromNode(
	const Mesh & mesh,
	const Node & node,
	const std::vector<int> & vecPossibleFaces,
	FindFaceStruct & aFindFaceStruct
) {
#ifdef VERBOSE
	printf("---\n");
#endif
	// Reset the FaceStruct
	aFindFaceStruct.vecFaceIndices.clear();
	aFindFaceStruct.vecFaceLocations.clear();
	aFindFaceStruct.loc = Face::NodeLocation_Undefined;

	// Loop through all faces to find overlaps
	// Note: This algorithm can likely be dramatically improved
	for (int l = 0; l < vecPossibleFaces.size(); l++) {
		Face::NodeLocation loc;
		int ixLocation;

		mesh.faces[vecPossibleFaces[l]].ContainsNode(
			mesh.nodes,
			node,
			loc,
			ixLocation);

		if (loc == Face::NodeLocation_Exterior) {
			continue;
		}

		if (aFindFaceStruct.loc == Face::NodeLocation_Undefined) {
			aFindFaceStruct.loc = loc;
		}

		// Node is in the interior of this face
		if (loc == Face::NodeLocation_Interior) {
			if (loc != aFindFaceStruct.loc) {
				_EXCEPTIONT("No consensus on location of Node");
			}

			aFindFaceStruct.vecFaceIndices.push_back(vecPossibleFaces[l]);
			aFindFaceStruct.vecFaceLocations.push_back(ixLocation);
			break;
		}

		// Node is on the edge of this face
		if (loc == Face::NodeLocation_Edge) {
			if (loc != aFindFaceStruct.loc) {
				_EXCEPTIONT("No consensus on location of Node");
			}

			aFindFaceStruct.vecFaceIndices.push_back(vecPossibleFaces[l]);
			aFindFaceStruct.vecFaceLocations.push_back(ixLocation);
		}

		// Node is at the corner of this face
		if (loc == Face::NodeLocation_Corner) {
			if (loc != aFindFaceStruct.loc) {
				_EXCEPTIONT("No consensus on location of Node");
			}

			aFindFaceStruct.vecFaceIndices.push_back(vecPossibleFaces[l]);
			aFindFaceStruct.vecFaceLocations.push_back(ixLocation);
		}
	}

	// Edges can only have two adjacent Faces
	if (aFindFaceStruct.loc == Face::NodeLocation_Edge) {
		if (aFindFaceStruct.vecFaceIndices.size() != 2) {
			_EXCEPTIONT("Multiple co-located edges detected");
		}
	}

	// Corners must have at least three faces
	if (aFindFaceStruct.loc == Face::NodeLocation_Corner) {
		if (aFindFaceStruct.vecFaceIndices.size() < 3) {
			_EXCEPTIONT("At least three Faces needed at each corner");
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

int FindFaceNearNode(
	const Mesh & mesh,
	int ixNode,
	const Node & nodeEnd,
	const Edge::Type edgetype
) {
#ifdef VERBOSE
	printf("---\n");
#endif
	static const double Tolerance = 1.0e-12;

	// Get the begin node
	const Node & nodeBegin = mesh.nodes[ixNode];

	// Get the set of faces adjacent this node
	const std::set<int> & setNearbyFaces = mesh.revnodearray[ixNode];

	if (setNearbyFaces.size() < 3) {
		_EXCEPTIONT("Insufficient Faces at Corner; at least three Faces expected");
	}

	// Direction along second curve projected onto
	// surface of the sphere at nodeBegin
	Node nodeDirA;

	GetLocalDirection(nodeBegin, nodeEnd, edgetype, nodeDirA);

	// Loop through all faces
	std::set<int>::const_iterator iter = setNearbyFaces.begin();
	for (; iter != setNearbyFaces.end(); iter++) {

		const Face & face = mesh.faces[*iter];

		int nEdges = face.edges.size();

		// Find the edge attached to this node
		int ixThisEdge = (-1);
		for (int j = 0; j < nEdges; j++) {
			int jNext = (j + 1) % nEdges;
			if (face[j] == face[jNext]) {
				continue;
			}
			if (ixNode == face[j]) {
				ixThisEdge = j;
				break;
			}
		}
		if (ixThisEdge == (-1)) {
			_EXCEPTIONT("Logic error");
		}

		int ixPrevEdge = (ixThisEdge + nEdges - 1) % nEdges;
		for (int j = ixPrevEdge; j != ixThisEdge; j--) {
			if (face[j] != face[ixThisEdge]) {
				break;
			}
			if (j == (-1)) {
				j = nEdges - 1;
			}
		}

		if (face[ixPrevEdge] == face[ixThisEdge]) {
			_EXCEPTIONT("Invalid Face: Only one distinct Edge found");
		}

		// Bounding edges on Node
		const Edge & edgeThis = face.edges[ixThisEdge];
		const Edge & edgePrev = face.edges[ixPrevEdge];

		if (edgeThis[0] != edgePrev[1]) {
			_EXCEPTIONT("Invalid Node indices on Edges");
		}

		const Node & node0 = mesh.nodes[edgePrev[0]];
		const Node & node1 = mesh.nodes[edgePrev[1]];
		const Node & node2 = mesh.nodes[edgeThis[1]];

		if (node1 != nodeBegin) {
			_EXCEPTIONT("Logic error");
		}

		// Direction vectors towards each of the nodes
		Node nodeDirL;
		Node nodeDirR;

		GetLocalDirection(nodeBegin, node0, edgePrev.type, nodeDirL);
		GetLocalDirection(nodeBegin, node2, edgeThis.type, nodeDirR);
/*
		// Verify coplanarity
		Node nodeX;
		nodeX.x = nodeDirA.y * nodeDirL.z - nodeDirA.z * nodeDirL.y;
		nodeX.y = nodeDirA.z * nodeDirL.x - nodeDirA.x * nodeDirL.z;
		nodeX.z = nodeDirA.x * nodeDirL.y - nodeDirA.y * nodeDirL.x;

		double dDotPlanar =
			  nodeX.x * nodeDirR.x
			+ nodeX.y * nodeDirR.y
			+ nodeX.z * nodeDirR.z;

		printf("Planarity: %1.10e\n", dDotPlanar);
*/
		double dMagA = nodeDirA.Magnitude();
		double dMagL = nodeDirL.Magnitude();
		double dMagR = nodeDirR.Magnitude();

		if (fabs(dMagA) < Tolerance) {
			_EXCEPTIONT("Zero magnitude of path Edge (possible zero Edge)");
		}
		if (fabs(dMagL) < Tolerance) {
			_EXCEPTIONT("Zero magnitude of incoming Edge");
		}
		if (fabs(dMagR) < Tolerance) {
			_EXCEPTIONT("Zero magnitude of outgoing Edge");
		}

		double dDotLR =
			  nodeDirL.x * nodeDirR.x
			+ nodeDirL.y * nodeDirR.y
			+ nodeDirL.z * nodeDirR.z;

		double dDotLA =
			  nodeDirL.x * nodeDirA.x
			+ nodeDirL.y * nodeDirA.y
			+ nodeDirL.z * nodeDirA.z;

		double dDotRA =
			  nodeDirR.x * nodeDirA.x
			+ nodeDirR.y * nodeDirA.y
			+ nodeDirR.z * nodeDirA.z;

		double dNormDotLR = dDotLR / dMagL / dMagR;
		double dNormDotLA = dDotLA / dMagL / dMagA;
		double dNormDotRA = dDotRA / dMagR / dMagA;

		// These values mimic the monotone structure of acos()
		double dAngleLR = 1.0 - dNormDotLR;
		double dAngleLA = 1.0 - dNormDotLA;
		double dAngleRA = 1.0 - dNormDotRA;

#ifdef VERBOSE
		printf("Face: %i\n", (*iter));
		printf("Type: %i %i %i\n", edgePrev.type, edgeThis.type, edgetype);
		nodeDirA.Print("DirA");
		nodeDirL.Print("DirL");
		nodeDirR.Print("DirR");
		printf("Angles: %1.5e %1.5e %1.5e\n", dAngleLR, dAngleLA, dAngleRA);
#endif

		// Sanity check for angles
		if (dAngleLR < - Tolerance) {
			_EXCEPTIONT("Logic error");
		}
		if (dAngleLA < - Tolerance) {
			_EXCEPTIONT("Logic error");
		}
		if (dAngleRA < - Tolerance) {
			_EXCEPTIONT("Logic error");
		}
		if (fabs(dAngleLR) < Tolerance) {
			_EXCEPTIONT("Zero angle detected, cannot continue");
		}
		if (fabs(dAngleLR - 2.0) < Tolerance) {
			_EXCEPTIONT("180 degree angle detected, cannot continue");
		}

		// Direction vector inside face
		if ((dAngleLA < dAngleLR - Tolerance) &&
			(dAngleRA < dAngleLR - Tolerance)
		) {
			return (*iter);
		}

		// This line makes zero angle with the L edge (special case)
		if (dAngleLA <= Tolerance) {

			// Coincident edges
			if (edgetype == edgePrev.type) {
				continue;
			}

			// Great circle arc emerging along a line of constant latitude
			if ((edgetype == Edge::Type_GreatCircleArc) &&
				(edgePrev.type == Edge::Type_ConstantLatitude)
			) {
				if ((node2.z > node1.z) && (nodeBegin.z < - Tolerance)) {
					return (*iter);
				} else if ((node2.z < node1.z) && (nodeBegin.z > Tolerance)) {
					return (*iter);
				} else {
					continue;
				}
			}

			// Line of constant latitude emerging along a great circle arc
			if ((edgetype == Edge::Type_ConstantLatitude) &&
				(edgePrev.type == Edge::Type_GreatCircleArc)
			) {
				if ((node2.z < node1.z) && (nodeBegin.z < - Tolerance)) {
					return (*iter);
				} else if ((node2.z > node1.z) && (nodeBegin.z > Tolerance)) {
					return (*iter);
				} else {
					continue;
				}
			}
		}

		// This line makes zero angle with the R edge (special case)
		if (dAngleRA <= Tolerance) {

			// Coincident edges
			if (edgetype == edgeThis.type) {
				return (*iter);
			}

			// Great circle arc emerging along a line of constant latitude
			if ((edgetype == Edge::Type_GreatCircleArc) &&
				(edgeThis.type == Edge::Type_ConstantLatitude)
			) {
				if ((node0.z > node1.z) && (nodeBegin.z < Tolerance)) {
					return (*iter);
				} else if ((node0.z < node1.z) && (nodeBegin.z > Tolerance)) {
					return (*iter);
				} else if (fabs(nodeBegin.z) <= Tolerance) {
					return (*iter);
				} else {
					continue;
				}
			}

			// Line of constant latitude emerging along a great circle arc
			if ((edgetype == Edge::Type_ConstantLatitude) &&
				(edgeThis.type == Edge::Type_GreatCircleArc)
			) {
				if ((node0.z < node1.z) && (nodeBegin.z < - Tolerance)) {
					return (*iter);
				} else if ((node0.z > node1.z) && (nodeBegin.z > Tolerance)) {
					return (*iter);
				} else if (fabs(nodeBegin.z) <= Tolerance) {
					return (*iter);
				} else {
					continue;
				}
			}

		}
	}

	_EXCEPTIONT("Logic Error: No exit Face found from Node");
}

///////////////////////////////////////////////////////////////////////////////

int FindFaceNearNodeOld(
	const Mesh & mesh,
	int ixNode,
	const Node & nodeEnd,
	const Edge::Type edgetype
) {
	static const double Tolerance = 1.0e-12;

	// Beginning node
	const Node & nodeBegin = mesh.nodes[ixNode];

	// Set of possible nodes
	std::vector< std::pair<int,int> > vecPossibleFaces;

	// Get the set of faces adjacent this node
	const std::set<int> & setNearbyFaces = mesh.revnodearray[ixNode];

	std::set<int>::const_iterator iter = setNearbyFaces.begin();
	for (; iter != setNearbyFaces.end(); iter++) {

		const Face & face = mesh.faces[*iter];

		// Edge id, if the point lies directly on an edge
		int iEdgeId = (-1);

		// Loop through all edges attached to this node
		int nAdjacentEdges = 0;

		int nPairedEdges = 0;

		for (int i = 0; i < face.edges.size(); i++) {

			const Edge & edge = face.edges[i];

			// Zero edge
			if (edge[0] == edge[1]) {
				continue;
			}

			// Node does not appear on edge
			if ((edge[0] != ixNode) && (edge[1] != ixNode)) {
				continue;
			}

			// Increment number of adjacent edges
			nAdjacentEdges++;

			// Get nodes
			const Node & na = mesh.nodes[edge[0]];
			const Node & nb = mesh.nodes[edge[1]];

			// Edge is a great circle arc
			if (edge.type == Edge::Type_GreatCircleArc) {
				double dDotNorm =
					  (na.y * nb.z - nb.y * na.z) * nodeEnd.x
					+ (nb.x * na.z - na.x * nb.z) * nodeEnd.y
					+ (na.x * nb.y - nb.x * na.y) * nodeEnd.z;

				if (fabs(dDotNorm) < Tolerance) {
					if (iEdgeId != (-1)) {
						_EXCEPTIONT("Already found on an edge");
					}
					iEdgeId = i;
				}
				if (dDotNorm > - Tolerance) {
					nPairedEdges++;
				}

			// Edge is a line of constant latitude
			} else if (edge.type == Edge::Type_ConstantLatitude) {
				double dAlignment = (na.x * nb.y - nb.x * na.y);
				double dDotNorm = dAlignment / fabs(dAlignment) * (nodeBegin.z - na.z);

				if (fabs(dDotNorm) < Tolerance) {
					if (iEdgeId != (-1)) {
						_EXCEPTIONT("Already found on an edge");
					}
					iEdgeId = i;
				}
				if (dDotNorm > - Tolerance) {
					nPairedEdges++;
				}

			} else {
				_EXCEPTIONT("Invalid EdgeType");
			}
		}

		if (nAdjacentEdges != 2) {
			_EXCEPTIONT("Logic error");
		}
		if (nPairedEdges == 2) {
			vecPossibleFaces.push_back(std::pair<int,int>(*iter, iEdgeId));
		}
	}

	// Only one possible face
	if (vecPossibleFaces.size() == 1) {
		return (vecPossibleFaces[0].first);
	}

	// Two faces share an edge and node lies close to this edge
	if (vecPossibleFaces.size() == 2) {
		const Face & face0 = mesh.faces[vecPossibleFaces[0].first];
		const Face & face1 = mesh.faces[vecPossibleFaces[1].first];

		const Edge & edge  = face0.edges[vecPossibleFaces[0].second];
		const Edge & edge1 = face1.edges[vecPossibleFaces[1].second];

		if (!(edge == edge1)) {
			_EXCEPTIONT("Non-coincident edge");
		}

		const Node & na = mesh.nodes[edge[0]];
		const Node & nb = mesh.nodes[edge[1]];

		if ((edgetype == Edge::Type_GreatCircleArc) &&
			(edge.type == Edge::Type_ConstantLatitude) &&
			(fabs(nodeEnd.z) > Tolerance)
		) {
			if ((nodeBegin == na) || (nodeBegin == nb)) {
				double dDotNorm =
					  (na.y * nb.z - nb.y * na.z) * nodeEnd.x
					+ (nb.x * na.z - na.x * nb.z) * nodeEnd.y
					+ (na.x * nb.y - nb.x * na.y) * nodeEnd.z;

				if (dDotNorm < 0.0) {
					return (vecPossibleFaces[0].first);
				} else if (dDotNorm > 0.0) {
					return (vecPossibleFaces[1].first);
				} else {
					_EXCEPTIONT("Unable to identify faces");
				}
			}
		}

		Node nodeDeltaA(
			nodeEnd.x - nodeBegin.x,
			nodeEnd.y - nodeBegin.y,
			nodeEnd.z - nodeBegin.z);

		Node nodeDeltaB(
			nb.x - na.x,
			nb.y - na.y,
			nb.z - na.z);

		double dDot =
			+ nodeDeltaA.x * nodeDeltaB.x
			+ nodeDeltaA.y * nodeDeltaB.y
			+ nodeDeltaA.z * nodeDeltaB.z;

		if (dDot > 0.0) {
			return (vecPossibleFaces[0].first);
		} else {
			return (vecPossibleFaces[1].first);
		}
	}

	_EXCEPTIONT("No face found matching criteria");
}

///////////////////////////////////////////////////////////////////////////////

int FindFaceNearNode(
	const Mesh & mesh,
	const Node & nodeBegin,
	const Node & nodeEnd,
	const Edge::Type edgetype,
	const FindFaceStruct & aFindFaceStruct
) {
	static const double Tolerance = 1.0e-12;

	const std::vector<int> & vecPossibleFaces =
		aFindFaceStruct.vecFaceIndices;

	// Check that there is actually a decision to be made
	if (vecPossibleFaces.size() < 2) {
		_EXCEPTIONT("vecPossibleFaces must contain at least two Faces.");

	// Node is on an edge
	} else if (aFindFaceStruct.loc == Face::NodeLocation_Edge) {

		// Get Faces and Edges
		if (aFindFaceStruct.vecFaceIndices.size() != 2) {
			_EXCEPTIONT("Logic failure");
		}

		const Face & face0 = mesh.faces[aFindFaceStruct.vecFaceIndices[0]];
		const Face & face1 = mesh.faces[aFindFaceStruct.vecFaceIndices[1]];

		const Edge & edge0 = face0.edges[aFindFaceStruct.vecFaceLocations[0]];
		const Edge & edge1 = face1.edges[aFindFaceStruct.vecFaceLocations[1]];

		const Node & node0 = mesh.nodes[edge0[0]];
		const Node & node1 = mesh.nodes[edge0[1]];

		if (edge0 != edge1) {
			_EXCEPTIONT("Logic failure");
		}

		// Calculate all intersections between edges
		std::vector<Node> nodeIntersections;

		bool fCoincident =
			CalculateEdgeIntersections(
				node0, node1, edge0.type,
				nodeBegin, nodeEnd, edgetype,
				nodeIntersections);

		// If edges are coincident check orientation to determine face 
		if (fCoincident) {

			// Calculate orientation
			double dOrientation =
				  (node1.x - node0.x) * (nodeEnd.x - nodeBegin.x)
				+ (node1.y - node0.y) * (nodeEnd.y - nodeBegin.y)
				+ (node1.z - node0.z) * (nodeEnd.z - nodeBegin.z);

			//if (edge0.type != edgetype) {
			//	_EXCEPTIONT("Logic error");
			//}

			if (dOrientation == 0.0) {
				_EXCEPTIONT("Logic error");

			// Oriented along with face0
			} else if (dOrientation > 0.0) {
				return aFindFaceStruct.vecFaceIndices[0];

			// Oriented along with face1
			} else {
				return aFindFaceStruct.vecFaceIndices[1];
			}
		}

#ifdef VERBOSE
		printf("Nodes: %i %i\n", edge0[0], edge0[1]);
		printf("Intersects: %lu\n", nodeIntersections.size());
		printf("type: %i %i\n", edge0.type, edgetype);
		node0.Print("node0");
		node1.Print("node1");
		nodeBegin.Print("nodeB");
		nodeEnd.Print("nodeE");
#endif

		if (nodeIntersections.size() == 0) {
			_EXCEPTIONT("Logic failure");
		}

		// Both Edges are lines of constant latitude
		if ((edge0.type == Edge::Type_ConstantLatitude) &&
			(edgetype == Edge::Type_ConstantLatitude)
		) {
			_EXCEPTIONT("Non-coincident lines of constant latitude found.");
		}

		// Great circle arc originating on a great circle arc
		if ((edge0.type == Edge::Type_GreatCircleArc) &&
			(edgetype == Edge::Type_GreatCircleArc)
		) {
			// Outward-pointing normal to first great circle arc
			Node nodeNormal;
			nodeNormal.x = (node0.y * node1.z) - (node0.z * node1.y);
			nodeNormal.y = (node0.z * node1.x) - (node0.x * node1.z);
			nodeNormal.z = (node0.x * node1.y) - (node0.y * node1.x);

			// Direction along second great circle arc projected onto
			// surface of the sphere at nodeBegin
			Node nodeDir;

			GetLocalDirection(nodeBegin, nodeEnd, edgetype, nodeDir);

			// Dot product to determine direction of motion
			double dDot =
				  nodeNormal.x * nodeDir.x
				+ nodeNormal.y * nodeDir.y
				+ nodeNormal.z * nodeDir.z;
/*
			printf("Norm: %1.5e %1.5e %1.5e\n", nodeNormal.x, nodeNormal.y, nodeNormal.z);
			printf("Dir: %1.5e %1.5e %1.5e\n", nodeDir.x, nodeDir.y, nodeDir.z);
			printf("%1.5e\n", dDot);
*/
			if (dDot > Tolerance) {
				return aFindFaceStruct.vecFaceIndices[0];

			} else if (dDot < -Tolerance) {
				return aFindFaceStruct.vecFaceIndices[1];

			} else {
				_EXCEPTIONT("Logic error");
			}
		}

		// Great circle arc originating on a line of constant latitude
		if ((edge0.type == Edge::Type_ConstantLatitude) &&
			(edgetype == Edge::Type_GreatCircleArc)
		) {
			// Check for equatorial failure
			if ((fabs(nodeBegin.z) < Tolerance) &&
				(nodeIntersections.size() > 1)
			) {
				_EXCEPTIONT("Equatorial failure in intersection");
			}

			// Find top Face and bottom Face
			int ixTopFace;
			int ixBottomFace;

			// Orientation of edge indicates which Face is on top
			if (IsPositivelyOrientedEdge(node0, node1)) {
				ixTopFace    = aFindFaceStruct.vecFaceIndices[0];
				ixBottomFace = aFindFaceStruct.vecFaceIndices[1];
			} else {
				ixTopFace    = aFindFaceStruct.vecFaceIndices[1];
				ixBottomFace = aFindFaceStruct.vecFaceIndices[0];
			}
/*
			printf("Pos: %i\n" ,IsPositivelyOrientedEdge(node0, node1));
			printf("Face: %i %i\n", aFindFaceStruct.vecFaceIndices[0], aFindFaceStruct.vecFaceIndices[1]);
			printf("nodeInt: %i\n", nodeIntersections.size());
*/
			// Only begin node intersects
			if (nodeIntersections.size() == 1) {
				if (nodeEnd.z > node0.z) {
					return ixTopFace;
				} else if (nodeEnd.z < node0.z) {
					return ixBottomFace;
				} else if (nodeBegin.z >= Tolerance) {
					return ixTopFace;
				} else if (nodeBegin.z <= Tolerance) {
					return ixBottomFace;
				} else {
					_EXCEPTIONT("Logic error");
				}

			// Great circle arc intersects twice
			} else if (nodeIntersections.size() == 2) {
				if (nodeBegin.z >= Tolerance) {
					return ixTopFace;
				} else if (nodeBegin.z <= Tolerance) {
					return ixBottomFace;
				} else {
					_EXCEPTIONT("Logic error");
				}

			} else {
				_EXCEPTIONT("Too many intersections detected");
			}
		}

		// Line of constant latitude originating on a great circle arc
		if ((edgetype == Edge::Type_ConstantLatitude) &&
			(edge0.type == Edge::Type_GreatCircleArc)
		) {
			// Exterior pointing normal to great circle arc
			Node nodeNormal;
			nodeNormal.x = (node1.y * node0.z) - (node1.z * node0.y);
			nodeNormal.y = (node1.z * node0.x) - (node1.x * node0.z);
			//nodeNormal.z = (node1.x * node0.y) - (node1.y * node0.x);

			// Direction of line of constant latitude
			Node nodeDir;

			GetLocalDirection(nodeBegin, nodeEnd, edgetype, nodeDir);

			// Dot product
			double dDot = nodeDir.x * nodeNormal.x + nodeDir.y * nodeNormal.y;
/*
			printf("F: %i %i\n",
				aFindFaceStruct.vecFaceIndices[0],
				aFindFaceStruct.vecFaceIndices[1]);
			printf("E0: %i %i\n", edge0[0], edge0[1]);
			printf("E1: %i %i\n", edge1[0], edge1[1]);
			nodeNormal.Print("Norm");
			nodeDir.Print("Dir ");
			printf("Dot: %1.5e\n", dDot);

			printf("%i %i\n",
				aFindFaceStruct.vecFaceIndices[0],
				aFindFaceStruct.vecFaceIndices[1]);
*/
			if (dDot > Tolerance) {
				return aFindFaceStruct.vecFaceIndices[1];
			} else if (dDot < - Tolerance) {
				return aFindFaceStruct.vecFaceIndices[0];
			}

			// Zero dot product (quadratic intersection)
			int ixTopFace;
			int ixBottomFace;

			// Orientation of edge indicates which Face is on top
			if (IsPositivelyOrientedEdge(node0, node1)) {
				ixTopFace    = aFindFaceStruct.vecFaceIndices[0];
				ixBottomFace = aFindFaceStruct.vecFaceIndices[1];
			} else {
				ixTopFace    = aFindFaceStruct.vecFaceIndices[1];
				ixBottomFace = aFindFaceStruct.vecFaceIndices[0];
			}

			if (nodeBegin.z > Tolerance) {
				return ixTopFace;

			} else if (nodeBegin.z < Tolerance) {
				return ixBottomFace;

			} else {
				_EXCEPTIONT("Logic error");
			}
		}

	// Node is at a corner
	} else if (aFindFaceStruct.loc == Face::NodeLocation_Corner) {

		if (aFindFaceStruct.vecFaceLocations.size() < 3) {
			_EXCEPTIONT("Logic error");
		}

		int ixLocation = aFindFaceStruct.vecFaceLocations[0];

		int ixFace = aFindFaceStruct.vecFaceIndices[0];

		int ixNode = mesh.faces[ixFace][ixLocation];

		return FindFaceNearNode(
			mesh,
			ixNode,
			nodeEnd,
			edgetype);

	// This function does not handle Interior or Exterior nodes
	} else {
		_EXCEPTIONT("Invalid Node location");
	}

/*
	// Locations
	Face::NodeLocation locConsensus;

	std::vector< std::pair<int,int> > vecPossibleFaces2;

	// Determine if the node is shared between edges or nodes
	for (int i = 0; i < vecPossibleFaces.size(); i++) {

		int ixPossibleFace = vecPossibleFaces[i];

		const Face & face = mesh.faces[ixPossibleFace];

		Face::NodeLocation loc;
		int ixLocation;

		face.ContainsNode(
			mesh.nodes,
			nodeBegin,
			loc,
			ixLocation);

		if (loc == Face::NodeLocation_Exterior) {
			_EXCEPTIONT("Node detected outside of Face");
		}

		if (i == 0) {
			locConsensus = loc;
		} else if (loc != locConsensus) {
			_EXCEPTIONT("Multiple consensus locations");
		}

		vecPossibleFaces2.push_back(
			std::pair<int,int>(ixPossibleFace, ixLocation));
	}

	// Node is coincident with multiple edges
	if (locConsensus == Face::NodeLocation_Corner) {

		// Find the nudged node towards nodeEnd with type edgetype
		Node nodeNudged;
		NudgeAlongEdge(nodeBegin, nodeEnd, edgetype, nodeNudged);

		//printf("A: %1.5e %1.5e %1.5e\n", nodeBegin.x, nodeBegin.y, nodeBegin.z);

		int ixFaceSuggestion = (-1);
		for (int i = 0; i < vecPossibleFaces.size(); i++) {

			// Determine if this edge lies between angle subtended by
			// the two adjacent edges.
			int ixPossibleFace = vecPossibleFaces[i];
			const Face & face = mesh.faces[ixPossibleFace];
			int ixCurrentEdge = vecPossibleFaces2[i].second;

			int nFaceEdgeCount = face.edges.size();

			int ixPrevEdge =
				(ixCurrentEdge + nFaceEdgeCount - 1) % nFaceEdgeCount;

			const Edge & edgeL = face.edges[ixPrevEdge];
			const Edge & edgeR = face.edges[ixCurrentEdge];

			int iEdgeSideL = EdgeSide(mesh, edgeL, nodeNudged);
			int iEdgeSideR = EdgeSide(mesh, edgeR, nodeNudged);

			if ((iEdgeSideL == (-1)) && (iEdgeSideR >= 0)) {
				if (ixFaceSuggestion != (-1)) {
					_EXCEPTIONT("Logic error");
				}
				ixFaceSuggestion = ixPossibleFace;
			}
		}
		return ixFaceSuggestion;
	}
*/
//////////////////////////////////////////////////////////////////////////////
/*
	// Find the nudged node towards nodeEnd with type edgetype
	double dNudgeAmt = 1.0e-6;

	Node nodeNudged;
	NudgeAlongEdge(nodeBegin, nodeEnd, edgetype, nodeNudged, dNudgeAmt);


	// Flag indicating that this node is on the edge between two elements
	bool fNodeOnEdge = false;

	// Set of possible nodes
	std::vector< std::pair<int,int> > vecPossibleFaces2;

	// Loop through all faces
	while (!fNodeOnEdge) {
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
				fNodeOnEdge = true;
				vecPossibleFaces2.push_back(
					std::pair<int,int>(ixPossibleFace, ixLocation));

			} else if (loc == Face::NodeLocation_Corner) {
				_EXCEPTIONT("Coincident nodes away from fixed node!");

			} else {
				return ixPossibleFace;
			}
		}

		// Reduce nudge size and try again
		dNudgeAmt /= 2.0;
		NudgeAlongEdge(nodeBegin, nodeEnd, edgetype, nodeNudged, dNudgeAmt);

		if (dNudgeAmt < Tolerance) {
			_EXCEPTIONT("No face found matching criteria");
		}
	}

	// No obvious face found; attempt to identify from available options
	//if (!fNodeOnEdge) {
	//	_EXCEPTIONT("No face found matching criteria");
	//}

	// Two faces share an edge and nudged node lies close to this edge
	if (vecPossibleFaces2.size() == 2) {
		const Face & face0 = mesh.faces[vecPossibleFaces2[0].first];
		const Face & face1 = mesh.faces[vecPossibleFaces2[1].first];

		const Edge & edge  = face0.edges[vecPossibleFaces2[0].second];
		const Edge & edge1 = face1.edges[vecPossibleFaces2[1].second];

		if (!(edge == edge1)) {
			_EXCEPTIONT("Non-coincident edge");
		}

		const Node & na = mesh.nodes[edge[0]];
		const Node & nb = mesh.nodes[edge[1]];

		if ((edgetype == Edge::Type_GreatCircleArc) &&
			(edge.type == Edge::Type_ConstantLatitude) &&
			(fabs(nodeEnd.z) > Tolerance)
		) {
			if ((nodeBegin == na) || (nodeBegin == nb)) {
				double dDotNorm =
					  (na.y * nb.z - nb.y * na.z) * nodeEnd.x
					+ (nb.x * na.z - na.x * nb.z) * nodeEnd.y
					+ (na.x * nb.y - nb.x * na.y) * nodeEnd.z;

				if (dDotNorm < 0.0) {
					return (vecPossibleFaces2[0].first);
				} else if (dDotNorm > 0.0) {
					return (vecPossibleFaces2[1].first);
				} else {
					_EXCEPTIONT("Unable to identify faces");
				}
			}
		}

		Node nodeDeltaA(
			nodeEnd.x - nodeBegin.x,
			nodeEnd.y - nodeBegin.y,
			nodeEnd.z - nodeBegin.z);

		Node nodeDeltaB(
			nb.x - na.x,
			nb.y - na.y,
			nb.z - na.z);

		double dDot =
			+ nodeDeltaA.x * nodeDeltaB.x
			+ nodeDeltaA.y * nodeDeltaB.y
			+ nodeDeltaA.z * nodeDeltaB.z;

		if (dDot > 0.0) {
			return (vecPossibleFaces2[0].first);
		} else {
			return (vecPossibleFaces2[1].first);
		}
	}

	_EXCEPTION1("Invalid number of possible detections on edge (%i)",
		vecPossibleFaces2.size());
*/
}

///////////////////////////////////////////////////////////////////////////////

int FindFaceNearNodeOld(
	const Mesh & mesh,
	const Node & nodeBegin,
	const Node & nodeEnd,
	const Edge::Type edgetype,
	const FindFaceStruct & aFindFaceStruct
) {
	const double Tolerance = 1.0e-12;

	const std::vector<int> & vecPossibleFaces =
		aFindFaceStruct.vecFaceIndices;

	// Find the nudged node towards nodeEnd with type edgetype
	double dNudgeAmt = 1.0e-6;

	Node nodeNudged;
	NudgeAlongEdge(nodeBegin, nodeEnd, edgetype, nodeNudged, dNudgeAmt);

	// Flag indicating that this node is on the edge between two elements
	bool fNodeOnEdge = false;

	// Set of possible nodes
	std::vector< std::pair<int,int> > vecPossibleFaces2;

	// Loop through all faces
	while (!fNodeOnEdge) {
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
				fNodeOnEdge = true;
				vecPossibleFaces2.push_back(
					std::pair<int,int>(ixPossibleFace, ixLocation));

			} else if (loc == Face::NodeLocation_Corner) {
				_EXCEPTIONT("Coincident nodes away from fixed node!");

			} else {
				return ixPossibleFace;
			}
		}

		// Reduce nudge size and try again
		dNudgeAmt /= 2.0;
		NudgeAlongEdge(nodeBegin, nodeEnd, edgetype, nodeNudged, dNudgeAmt);

		if (dNudgeAmt < Tolerance) {
			_EXCEPTIONT("No face found matching criteria");
		}
	}

	// No obvious face found; attempt to identify from available options
	//if (!fNodeOnEdge) {
	//	_EXCEPTIONT("No face found matching criteria");
	//}

	// Two faces share an edge and nudged node lies close to this edge
	if (vecPossibleFaces2.size() == 2) {
		const Face & face0 = mesh.faces[vecPossibleFaces2[0].first];
		const Face & face1 = mesh.faces[vecPossibleFaces2[1].first];

		const Edge & edge  = face0.edges[vecPossibleFaces2[0].second];
		const Edge & edge1 = face1.edges[vecPossibleFaces2[1].second];

		if (!(edge == edge1)) {
			_EXCEPTIONT("Non-coincident edge");
		}

		const Node & na = mesh.nodes[edge[0]];
		const Node & nb = mesh.nodes[edge[1]];

		if ((edgetype == Edge::Type_GreatCircleArc) &&
			(edge.type == Edge::Type_ConstantLatitude) &&
			(fabs(nodeEnd.z) > Tolerance)
		) {
			if ((nodeBegin == na) || (nodeBegin == nb)) {
				double dDotNorm =
					  (na.y * nb.z - nb.y * na.z) * nodeEnd.x
					+ (nb.x * na.z - na.x * nb.z) * nodeEnd.y
					+ (na.x * nb.y - nb.x * na.y) * nodeEnd.z;

				if (dDotNorm < 0.0) {
					return (vecPossibleFaces2[0].first);
				} else if (dDotNorm > 0.0) {
					return (vecPossibleFaces2[1].first);
				} else {
					_EXCEPTIONT("Unable to identify faces");
				}
			}
		}

		Node nodeDeltaA(
			nodeEnd.x - nodeBegin.x,
			nodeEnd.y - nodeBegin.y,
			nodeEnd.z - nodeBegin.z);

		Node nodeDeltaB(
			nb.x - na.x,
			nb.y - na.y,
			nb.z - na.z);

		double dDot =
			+ nodeDeltaA.x * nodeDeltaB.x
			+ nodeDeltaA.y * nodeDeltaB.y
			+ nodeDeltaA.z * nodeDeltaB.z;

		if (dDot > 0.0) {
			return (vecPossibleFaces2[0].first);
		} else {
			return (vecPossibleFaces2[1].first);
		}
	}

	_EXCEPTION1("Invalid number of possible detections on edge (%i)",
		vecPossibleFaces2.size());
}

/*
///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Returns true if two elements overlap.
///	</summary>
bool AreElementsOverlapped(
	const Mesh & meshFirst,
	const Mesh & meshSecond,
	int iElementFirst,
	int iElementSecond
) {

}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Find all elements on meshSecond that overlap elements on meshFirst.
///	</summary>
void FindAllOverlapElements(
	const Mesh & meshFirst,
	const Mesh & meshSecond,
	int iFirstElement
	std::vector<int> & vecOverlapElements
) {

	std::set<int> setOverlapElements;

	// Find the starting face on the second mesh
	FindFaceFromNode(meshSecond, nodeCurrent, vecOverlapElements);

	for (int i = 0; i < vecOverlapElements.size(); i++) {
		setOverlapElements.insert(vecOverlapElements[i]);
	}

	// No faces found
	if (vecOverlapElements.size() == 0) {
		_EXCEPTIONT("No initial face found!");
	}

	// Multiple faces found.
	if (vecOverlapElements.size() > 1) {
	}
}
*/
///////////////////////////////////////////////////////////////////////////////

void GenerateOverlapMesh(
	const Mesh & meshFirst,
	const Mesh & meshSecond,
	Mesh & meshOverlap,
	int nVerbosity
) {
	meshOverlap.Clear();

	// Get the two NodeVectors
	const NodeVector & nodevecFirst = meshFirst.nodes;
	const NodeVector & nodevecSecond = meshSecond.nodes;
	const NodeVector & nodevecOverlap = meshOverlap.nodes;

	// Construct the coincident node vector
	std::vector<int> vecSecondNodeMap;
	int nCoincidentNodes =
		BuildCoincidentNodeVector(
			meshFirst, meshSecond, vecSecondNodeMap);

	Announce("Number of coincident nodes [%i]", nCoincidentNodes);

	// Insert all nodes from the two NodeVectors
	for (int i = 0; i < nodevecFirst.size(); i++) {
		meshOverlap.nodes.push_back(nodevecFirst[i]);
	}
	const int ixOverlapSecondNodesBegin = meshOverlap.nodes.size();

	for (int i = 0; i < nodevecSecond.size(); i++) {
		meshOverlap.nodes.push_back(nodevecSecond[i]);
		if (vecSecondNodeMap[i] != InvalidNode) {
			continue;
		}
	}
	const int ixOverlapNewNodesBegin = meshOverlap.nodes.size();

	// Build the meshSecond node map
	for (int i = 0; i < vecSecondNodeMap.size(); i++) {
		if (vecSecondNodeMap[i] == InvalidNode) {
			vecSecondNodeMap[i] = ixOverlapSecondNodesBegin + i;
		}
	}

	int ixCurrentFirstFace = 0;

#pragma message "OpenMP here"
	//for (; ixCurrentFirstFace < meshFirst.faces.size(); ixCurrentFirstFace++) {
	for (int ixCurrentFirstFace = 389; ixCurrentFirstFace < 390; ixCurrentFirstFace++) {

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
		FindFaceStruct aFindFaceStruct;
		FindFaceFromNode(meshSecond, nodeCurrent, aFindFaceStruct);

		// No faces found
		if (aFindFaceStruct.vecFaceIndices.size() == 0) {
			_EXCEPTIONT("No initial face found!");
		}

		// Current face on second mesh
		int ixCurrentSecondFace = aFindFaceStruct.vecFaceIndices[0];

		// This node lies on the boundary between faces; find the
		// actual starting face by nudging along FirstEdge.
		if (aFindFaceStruct.vecFaceIndices.size() > 1) {

			ixCurrentSecondFace =
				FindFaceNearNode(
					meshSecond,
					nodeCurrent,
					nodevecFirst[faceFirstCurrent[1]],
					faceFirstCurrent.edges[0].type,
					aFindFaceStruct);
/*
			int ixOldSecondFace =
				FindFaceNearNodeOld(
					meshSecond,
					nodeCurrent,
					nodevecFirst[faceFirstCurrent[1]],
					faceFirstCurrent.edges[0].type,
					aFindFaceStruct);

			printf("A: %i %i\n", ixCurrentSecondFace, ixOldSecondFace);

			if (ixCurrentSecondFace != ixOldSecondFace) {
				_EXCEPTION();
			}
*/
		}

#ifdef VERBOSE
		printf("\nStarting Node: %i", meshFirst.faces[ixCurrentFirstFace][0]);
		printf("\nNext Node: %i", meshFirst.faces[ixCurrentFirstFace][1]);
#endif
		printf("\nFaces: %i %i\n", ixCurrentFirstFace, ixCurrentSecondFace);

		// Trace along all edges of current face
		for (int i = 0; i < faceFirstCurrent.edges.size(); i++) {

			std::cout << ixCurrentSecondFace << std::endl;

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
							//meshOverlap.nodes[edgeFirstCurrent[0]],
							meshOverlap.nodes[edgeFirstCurrent[1]],
							edgeFirstCurrent.type,
							nodevecSecond[edgeSecondCurrent[0]],
							nodevecSecond[edgeSecondCurrent[1]],
							edgeSecondCurrent.type,
							nodeIntersections,
							false);
/*
					for (int i = 0; i < nodeIntersections.size(); i++) {
						if (nodeIntersections[i] == meshOverlap.nodes[ixOverlapNodeCurrent]) {
							nodeIntersections.erase(nodeIntersections.begin()+i);
							i--;
						}
					}
*/
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
/*
					Face::NodeLocation loc;
					int ixLocation;

					meshSecond.faces[ixCurrentSecondFace].ContainsNode(
						meshSecond.nodes,
						meshOverlap.nodes[edgeFirstCurrent[1]],
						loc,
						ixLocation);

					if (loc == Face::NodeLocation_Edge) {
						_EXCEPTIONT("Intersection failure");
					}
*/
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

				// Find next face on meshSecond
				const Edge & edgeSecondCurrent =
					faceSecondCurrent.edges[ixIntersectionSecondEdge];

				const Node & nodeSecondEdge0 =
					nodevecSecond[edgeSecondCurrent[0]];
				const Node & nodeSecondEdge1 =
					nodevecSecond[edgeSecondCurrent[1]];

				int ixOverlapNodeNext;

				// Special case:  Intersection with Edge is exactly an
				// beginpoint / endpoint of edgeFirstCurrent
				if (nodeIntersections[0] ==
				        meshOverlap.nodes[edgeFirstCurrent[1]]
				) {
					// Next edge
					int iNext = (i + 1) % faceFirstCurrent.edges.size();

					const Edge & edgeFirstNext =
						faceFirstCurrent.edges[iNext];

					// Update SecondMesh face
					EdgeMapConstIterator iter =
						meshSecond.edgemap.find(edgeSecondCurrent);

					if (iter == meshSecond.edgemap.end()) {
						_EXCEPTIONT("Logic error");
					}

					int ixNextSecondFace;

					// Path hits the beginpoint of the Edge
					if (nodeIntersections[0] == nodeSecondEdge0) {
						ixNextSecondFace =
							FindFaceNearNode(
								meshSecond,
								edgeSecondCurrent[0],
								nodevecFirst[edgeFirstNext[1]],
								edgeFirstNext.type);
/*
						if (ixNextSecondFace == ixCurrentSecondFace) {
							printf("WARNING: Face does not change across Edge (1)\n");
						}
*/
						// If face changes insert a Node bifurcation
						if (ixNextSecondFace != ixCurrentSecondFace) {
							vecTracedPath.push_back(PathSegment(
								ixOverlapNodeCurrent,
								edgeFirstCurrent[1],
								edgeFirstCurrent.type,
								ixCurrentFirstFace,
								ixCurrentSecondFace,
								IntersectType_Node,
								ixIntersectionSecondEdge));
						}

					// Path hits the endpoint of the Edge
					} else if (nodeIntersections[0] == nodeSecondEdge1) {
						ixNextSecondFace =
							FindFaceNearNode(
								meshSecond,
								edgeSecondCurrent[1],
								nodevecFirst[edgeFirstNext[1]],
								edgeFirstNext.type);
/*
						if (ixNextSecondFace == ixCurrentSecondFace) {
							printf("WARNING: Face does not change across Edge (2)\n");
						}
*/
						// If face changes insert a Node bifuraction
						if (ixNextSecondFace != ixCurrentSecondFace) {
							vecTracedPath.push_back(PathSegment(
								ixOverlapNodeCurrent,
								edgeFirstCurrent[1],
								edgeFirstCurrent.type,
								ixCurrentFirstFace,
								ixCurrentSecondFace,
								IntersectType_Node,
								(ixIntersectionSecondEdge + 1)
									% faceSecondCurrent.edges.size()));
						}

					// Path hits the Edge directly
					} else {
						const FacePair & facepair = iter->second;

						std::vector<int> vecPossibleFaces;
						vecPossibleFaces.push_back(facepair[0]);
						vecPossibleFaces.push_back(facepair[1]);

						FindFaceStruct aNextFindFaceStruct;

						FindFaceFromNode(
							meshSecond,
							nodevecFirst[edgeFirstNext[0]],
							vecPossibleFaces,
							aNextFindFaceStruct);

						ixNextSecondFace =
							FindFaceNearNode(
								meshSecond,
								nodevecFirst[edgeFirstNext[0]],
								nodevecFirst[edgeFirstNext[1]],
								edgeFirstNext.type,
								aNextFindFaceStruct);
/*
						if (ixNextSecondFace == ixCurrentSecondFace) {
							printf("WARNING: Face does not change across Edge (3)\n");
						}
*/
/*
						int ixOldSecondFace =
							FindFaceNearNodeOld(
								meshSecond,
								nodevecFirst[edgeFirstNext[0]],
								nodevecFirst[edgeFirstNext[1]],
								edgeFirstNext.type,
								aNextFindFaceStruct);

						printf("B: %i %i\n", ixNextSecondFace, ixOldSecondFace);

						if (ixCurrentSecondFace != ixOldSecondFace) {
							_EXCEPTION();
						}
*/
						// If face changes insert an Edge bifurcation
						if (ixNextSecondFace != ixCurrentSecondFace) {
							vecTracedPath.push_back(PathSegment(
								ixOverlapNodeCurrent,
								edgeFirstCurrent[1],
								edgeFirstCurrent.type,
								ixCurrentFirstFace,
								ixCurrentSecondFace,
								ixIntersectionSecondEdge,
								edgeSecondCurrent));
						}
					}

					// Remain on the same face
					if (ixNextSecondFace == ixCurrentSecondFace) {
						vecTracedPath.push_back(PathSegment(
							ixOverlapNodeCurrent,
							edgeFirstCurrent[1],
							edgeFirstCurrent.type,
							ixCurrentFirstFace,
							ixCurrentSecondFace,
							IntersectType_None,
							ixIntersectionSecondEdge));

					}

					// Update OverlapNodeCurrent
					ixOverlapNodeCurrent = edgeFirstCurrent[1];

					ixCurrentSecondFace = ixNextSecondFace;

					break;
				}

				// FirstEdge hits nodeSecondEdge0
				if (nodeIntersections[0] == nodeSecondEdge0) {
					ixOverlapNodeNext =
						vecSecondNodeMap[edgeSecondCurrent[0]];

					vecTracedPath.push_back(PathSegment(
						ixOverlapNodeCurrent,
						ixOverlapNodeNext,
						edgeFirstCurrent.type,
						ixCurrentFirstFace,
						ixCurrentSecondFace,
						IntersectType_Node,
						ixIntersectionSecondEdge));

					int ixPrevSecondFace = ixCurrentSecondFace;

					ixCurrentSecondFace =
						FindFaceNearNode(
							meshSecond,
							edgeSecondCurrent[0],
							nodeFirstEnd,
							edgeFirstCurrent.type);

					if (ixPrevSecondFace == ixCurrentSecondFace) {
						printf("WARNING: Face does not change across Edge (4)\n");
					}

					ixOverlapNodeCurrent = ixOverlapNodeNext;

					if (ixOverlapNodeNext == edgeFirstCurrent[1]) {
						break;
					}

					continue;
				
				// FirstEdge hits nodeSecondEdge1
				} else if (nodeIntersections[0] == nodeSecondEdge1) {
					ixOverlapNodeNext =
						vecSecondNodeMap[edgeSecondCurrent[1]];

					vecTracedPath.push_back(PathSegment(
						ixOverlapNodeCurrent,
						ixOverlapNodeNext,
						edgeFirstCurrent.type,
						ixCurrentFirstFace,
						ixCurrentSecondFace,
						IntersectType_Node,
						(ixIntersectionSecondEdge + 1)
							% faceSecondCurrent.edges.size()));

					int ixPrevSecondFace = ixCurrentSecondFace;

					ixCurrentSecondFace =
						FindFaceNearNode(
							meshSecond,
							edgeSecondCurrent[1],
							nodeFirstEnd,
							edgeFirstCurrent.type);

					if (ixPrevSecondFace == ixCurrentSecondFace) {
						printf("WARNING: Face does not change across Edge (5)\n");
					}

					ixOverlapNodeCurrent = ixOverlapNodeNext;

					if (ixOverlapNodeNext == edgeFirstCurrent[1]) {
						break;
					}

					continue;

				// General intersection between edgeFirstCurrent and
				// edgeSecondCurrent.
				} else {

					// Push a new intersection into the array of nodes
					ixOverlapNodeNext =
						static_cast<int>(meshOverlap.nodes.size());

					meshOverlap.nodes.push_back(nodeIntersections[0]);

					// Intersection found with edge
					vecTracedPath.push_back(PathSegment(
						ixOverlapNodeCurrent,
						ixOverlapNodeNext,
						edgeFirstCurrent.type,
						ixCurrentFirstFace,
						ixCurrentSecondFace,
						ixIntersectionSecondEdge,
						edgeSecondCurrent));

					// Update OverlapNodeCurrent
					ixOverlapNodeCurrent = ixOverlapNodeNext;

					// Update SecondMesh face
					EdgeMapConstIterator iter =
						meshSecond.edgemap.find(edgeSecondCurrent);

					if (iter == meshSecond.edgemap.end()) {
						_EXCEPTIONT("Logic error");
					}

					const FacePair & facepair = iter->second;

					std::vector<int> vecPossibleFaces;
					vecPossibleFaces.push_back(facepair[0]);
					vecPossibleFaces.push_back(facepair[1]);

					FindFaceStruct aNextFindFaceStruct;

					FindFaceFromNode(
						meshSecond,
						nodeIntersections[0],
						vecPossibleFaces,
						aNextFindFaceStruct);

					int ixPrevSecondFace = ixCurrentSecondFace;

					ixCurrentSecondFace = 
						FindFaceNearNode(
							meshSecond,
							nodeIntersections[0],
							meshOverlap.nodes[edgeFirstCurrent[1]],
							edgeFirstCurrent.type,
							aNextFindFaceStruct);

					if (ixPrevSecondFace == ixCurrentSecondFace) {
						printf("WARNING: Face does not change across Edge (3)\n");
					}
/*
					int ixOldSecondFace =
						FindFaceNearNodeOld(
							meshSecond,
							nodeIntersections[0],
							meshOverlap.nodes[edgeFirstCurrent[1]],
							edgeFirstCurrent.type,
							aNextFindFaceStruct);

					printf("C: %i %i\n", ixCurrentSecondFace, ixOldSecondFace);

					if (ixCurrentSecondFace != ixOldSecondFace) {
						printf("Type: %i\n", edgeFirstCurrent.type);
						printf("Loc: %i\n", aNextFindFaceStruct.loc);
						_EXCEPTION();
					}
*/
				}
			}
		}

		///////////////////////////////////////////////////////////////////////
		// Find faces associated with the TracedPath

		// Array indicating which elements of vecTracedPath have been used
		std::vector<bool> vecTracedPathUsed;
		vecTracedPathUsed.resize(vecTracedPath.size(), false);

		// The set of interior faces (from the Second mesh)
		std::set<int> setSecondFacesAdded;
		for (int j = 0; j < vecTracedPath.size(); j++) {
			setSecondFacesAdded.insert(vecTracedPath[j].ixSecondFace);

#ifdef VERBOSE
			printf("%i %i : %i\n", vecTracedPath[j][0], vecTracedPath[j][1], vecTracedPath[j].ixSecondFace);
#endif
		}

		// Set of faces from meshSecond that should be added
		std::set<int> setSecondFacesToAdd;

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
				int nEdgesCompleted = 0;
				for (;;) {

					const Edge & edgeSecondCurrent =
						faceSecondCurrent.edges[ixCurrentSecondEdge];

					// Check for infinite loop
					if (nEdgesCompleted > faceSecondCurrent.edges.size()) {
						_EXCEPTIONT("Possible infinite loop - aborting");
						//printf("Possible infinite loop - aborting\n");
						//goto Done;
					}

					nEdgesCompleted++;

					// Identical endpoints; advance the edge
					if (edgeSecondCurrent[0] == edgeSecondCurrent[1]) {
						ixCurrentSecondEdge =
							(ixCurrentSecondEdge + 1)
								% faceSecondCurrent.edges.size();

						ixCurrentOverlapNode =
							vecSecondNodeMap[edgeSecondCurrent[1]];

						continue;
					}

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

						if (ixCurrentOverlapNode == vecTracedPath[k][1]) {
							continue;
						}

						// Check for node intersections
						if (inttype == IntersectType_Node) {
							//printf("%i %i: %i %i %i\n", j, k, vecTracedPath[k][1], vecSecondNodeMap[edgeSecondCurrent[0]], vecSecondNodeMap[edgeSecondCurrent[1]]);
							if (vecTracedPath[k][1] ==
								vecSecondNodeMap[edgeSecondCurrent[0]]
							) {
								ixExitNode = vecTracedPath[k][1];
								break;
							}

							if (vecTracedPath[k][1] ==
								vecSecondNodeMap[edgeSecondCurrent[1]]
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

					// Add the interior face to the list of faces to be added
					EdgeMapConstIterator iter =
						meshSecond.edgemap.find(edgeSecondCurrent);

					if (iter == meshSecond.edgemap.end()) {
						_EXCEPTIONT("Logic error");
					}

					const FacePair & facepair = iter->second;

					if (facepair[0] == ixCurrentSecondFace) {
						setSecondFacesToAdd.insert(facepair[1]);
					} else if (facepair[1] == ixCurrentSecondFace) {
						setSecondFacesToAdd.insert(facepair[0]);
					} else {
						_EXCEPTIONT("Logic error");
					}

					// An exit node is found
					if (ixExitNode != InvalidNode) {

						// Check if the traced path becomes active
						int jNext = (k + 1) % vecTracedPath.size();
						if (vecTracedPath[jNext].ixSecondFace ==
							ixCurrentSecondFace
						) {

							printf("S: %i %i\n",
								ixCurrentOverlapNode, ixExitNode);

							faceOverlap.edges.push_back(Edge(
								ixCurrentOverlapNode,
								ixExitNode,
								edgeSecondCurrent.type));

							j = jNext;
							if (ixExitNode == ixOverlapOriginNode) {
								goto ContinueToNextFace;
							} else {
								break;
							}
						}
					}

					printf("T: %i (%i) %i\n",
						ixCurrentOverlapNode,
						vecSecondNodeMap[edgeSecondCurrent[0]],
						vecSecondNodeMap[edgeSecondCurrent[1]]);

					// Push this edge into the overlap mesh
					faceOverlap.edges.push_back(Edge(
						ixCurrentOverlapNode,
						vecSecondNodeMap[edgeSecondCurrent[1]],
						edgeSecondCurrent.type));

					// Advance the edge
					ixCurrentSecondEdge =
						(ixCurrentSecondEdge + 1)
							% faceSecondCurrent.edges.size();

					ixCurrentOverlapNode =
						vecSecondNodeMap[edgeSecondCurrent[1]];

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

		///////////////////////////////////////////////////////////////////////
		// Find interior faces from meshSecond to add to meshOverlap

		// Remove added faces from set of faces to add
		std::set<int>::iterator iterToAdd;
		std::set<int>::iterator iterAdded = setSecondFacesAdded.begin();
		for (; iterAdded != setSecondFacesAdded.end(); iterAdded++) {
			iterToAdd = setSecondFacesToAdd.find(*iterAdded);

			if (iterToAdd != setSecondFacesToAdd.end()) {
				setSecondFacesToAdd.erase(iterToAdd);
			}
		}

		iterToAdd = setSecondFacesToAdd.begin();
		while (iterToAdd != setSecondFacesToAdd.end()) {

			const Face & faceSecondCurrent =
				meshSecond.faces[*iterToAdd];

			// Add this face to meshOverlap
			Face faceOverlapCurrent(faceSecondCurrent.edges.size());

			for (int i = 0; i < faceOverlapCurrent.edges.size(); i++) {
				faceOverlapCurrent.edges[i][0] =
					vecSecondNodeMap[faceSecondCurrent.edges[i][0]];
				faceOverlapCurrent.edges[i][1] =
					vecSecondNodeMap[faceSecondCurrent.edges[i][1]];
				faceOverlapCurrent.edges[i].type =
					faceSecondCurrent.edges[i].type;
			}

			meshOverlap.faces.push_back(faceOverlapCurrent);

			setSecondFacesAdded.insert(*iterToAdd);

			// Add further interior faces
			bool fMoreFacesToAdd = false;

			for (int i = 0; i < faceSecondCurrent.edges.size(); i++) {

				if (faceSecondCurrent.edges[i][0] ==
					faceSecondCurrent.edges[i][1]
				) {
					continue;
				}

				EdgeMapConstIterator iterEdge =
					meshSecond.edgemap.find(faceSecondCurrent.edges[i]);

				if (iterEdge == meshSecond.edgemap.end()) {
					_EXCEPTIONT("Edge not found in EdgeMap");
				}

				int iOtherFace;
				if (iterEdge->second[0] == *iterToAdd) {
					iOtherFace = iterEdge->second[1];
				} else if (iterEdge->second[1] == *iterToAdd) {
					iOtherFace = iterEdge->second[0];
				} else {
					_EXCEPTIONT("EdgeMap consistency error");
				}

				if (setSecondFacesAdded.find(iOtherFace) ==
					setSecondFacesAdded.end()
				) {
					setSecondFacesToAdd.insert(iOtherFace);
					fMoreFacesToAdd = true;
				}
			}

			if (fMoreFacesToAdd) {
				setSecondFacesToAdd.erase(iterToAdd);
				iterToAdd = setSecondFacesToAdd.begin();
			} else {
				iterToAdd++;
			}
		}
	}
/*
Done:
	meshOverlap.faces.clear();
	meshOverlap.faces.push_back(meshFirst.faces[4003]);
	meshOverlap.faces.push_back(meshSecond.faces[62613]);
	for (int i = 0; i < meshOverlap.faces[1].edges.size(); i++) {
		meshOverlap.faces[1].edges[i][0] =
			vecSecondNodeMap[meshOverlap.faces[1].edges[i][0]];
	}
*/
}

///////////////////////////////////////////////////////////////////////////////

