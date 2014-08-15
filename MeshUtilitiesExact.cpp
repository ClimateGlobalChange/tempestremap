///////////////////////////////////////////////////////////////////////////////
///
///	\file    MeshUtilitiesExact.cpp
///	\author  Paul Ullrich
///	\version August 7, 2014
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

#include "MeshUtilitiesExact.h"

#include "Exception.h"

#ifdef USE_EXACT_ARITHMETIC

///////////////////////////////////////////////////////////////////////////////

bool MeshUtilitiesExact::AreNodesEqual(
	const Node & node0,
	const Node & node1
) {
	Node nodeCross = CrossProductX(node0, node1);

	if (nodeCross.fx.IsZero() &&
		nodeCross.fy.IsZero() &&
		nodeCross.fz.IsZero()
	) {
		return true;
	}

	return false;
}

///////////////////////////////////////////////////////////////////////////////

void MeshUtilitiesExact::ContainsNode(
	const Face & face,
	const NodeVector & nodevec,
	const Node & node,
	Face::NodeLocation & loc,
	int & ixLocation
) const {

	// Set of edges which "contain" this node
	std::set<int> setContainedEdgeIx;

	// Loop through all Edges of this face
	for (int i = 0; i < face.edges.size(); i++) {

		// Paired Edges means Edge is invalid
		if (face.edges[i][0] == face.edges[i][1]) {
			_EXCEPTIONT("Zero Edge detected");
		}

		// Check which side of the Face this Edge is on
		const Node & na = nodevec[face.edges[i][0]];
		const Node & nb = nodevec[face.edges[i][1]];

		if (face.edges[i].type == Edge::Type_GreatCircleArc) {

			FixedPoint fpDotNorm = DotProductX(CrossProductX(na, nb), node);
/*
			Node nx = CrossProductX(na, nb);
			nx.PrintMX();
			FixedPoint fpX = nx.fx * node.fx;
			FixedPoint fpY = nx.fy * node.fy;
			FixedPoint fpZ = nx.fz * node.fz;

			printf("NOW:\n");
			FixedPoint fpDot = fpX + fpZ;

			node.PrintMX(); printf("\n");
			fpX.Print(); printf("\n");
			fpY.Print(); printf("\n");
			fpZ.Print(); printf("\n");
*/
/*
			printf("a := Vector("); na.PrintMX(); printf(");");
			printf("b := Vector("); na.PrintMX(); printf(");");
			printf("N: "); fpDotNorm.Print(); printf("\n");
*/
/*
			Node nodeCross = CrossProduct(na, nb);
			Node nodeCrossX = CrossProductX(na, nb);
			long double dDotNorm = DotProduct(CrossProduct(na, nb), node);

			FixedPoint fp0 = nodeCrossX.fx * node.fx;
			FixedPoint fp1 = nodeCrossX.fy * node.fy;
			FixedPoint fp2 = nodeCrossX.fz * node.fz;

			printf("(X) %1.15e : ", nodeCross.x); nodeCrossX.fx.Print(); printf("\n");
			printf("(Y) %1.15e : ", nodeCross.y); nodeCrossX.fy.Print(); printf("\n");
			printf("(Z) %1.15e : ", nodeCross.z); nodeCrossX.fz.Print(); printf("\n");

			printf("(X) %1.15e : ", nodeCross.x * node.x); fp0.Print(); printf("\n");
			printf("(Y) %1.15e : ", nodeCross.y * node.y); fp1.Print(); printf("\n");
			printf("(Z) %1.15e : ", nodeCross.z * node.z); fp2.Print(); printf("\n");

			printf("%1.15e : ", dDotNorm); fpDotNorm.Print(); printf("\n");
*/
			if (fpDotNorm.IsNegative()) {
				loc = Face::NodeLocation_Exterior;
				ixLocation = 0;
				return;
			}
			if (fpDotNorm.IsZero()) {
				setContainedEdgeIx.insert(i);
			}

		} else if (face.edges[i].type == Edge::Type_ConstantLatitude) {
			_EXCEPTIONT("Unimplemented");

		} else {
			_EXCEPTIONT("Invalid EdgeType");
		}
	}

	// Check if the node is contained on an edge
	if (setContainedEdgeIx.size() == 1) {
		loc = Face::NodeLocation_Edge;
		ixLocation = *(setContainedEdgeIx.begin());
		return;
	}

	// Node is coincident with a corner of this face
	if (setContainedEdgeIx.size() == 2) {

		std::set<int>::iterator iter;

		iter = setContainedEdgeIx.begin();
		int ix0 = *(iter);

		iter++;
		int ix1 = *(iter);

		if ((ix0 == 0) && (ix1 != 1)) {
			ixLocation = 0;
		} else {
			ixLocation = ix1;
		}

		loc = Face::NodeLocation_Corner;
		return;
	}

	// Node occurs in more than two edges; error.
	if (setContainedEdgeIx.size() > 2) {
		_EXCEPTIONT("Logic error: Node occurs in more than two edges");
	}

	// Default; node occurs in the interior of the face
	loc = Face::NodeLocation_Interior;
	ixLocation = 0;
	return;
}

///////////////////////////////////////////////////////////////////////////////

bool MeshUtilitiesExact::CalculateEdgeIntersections(
	const Node & nodeFirstBegin,
	const Node & nodeFirstEnd,
	Edge::Type typeFirst,
	const Node & nodeSecondBegin,
	const Node & nodeSecondEnd,
	Edge::Type typeSecond,
	std::vector<Node> & nodeIntersections,
	bool fIncludeFirstBeginNode
) {
	
	// Make a locally modifyable version of the Nodes
	Node node11;
	Node node12;
	Node node21;
	Node node22;

	// Second edge is a line of constant latitude; first is a great circle arc
	if ((typeFirst  == Edge::Type_ConstantLatitude) &&
		(typeSecond == Edge::Type_GreatCircleArc)
	) {
		node11 = nodeSecondBegin;
		node12 = nodeSecondEnd;
		node21 = nodeFirstBegin;
		node22 = nodeFirstEnd;

		typeFirst = Edge::Type_GreatCircleArc;
		typeSecond = Edge::Type_ConstantLatitude;

	} else {
		node11 = nodeFirstBegin;
		node12 = nodeFirstEnd;
		node21 = nodeSecondBegin;
		node22 = nodeSecondEnd;
	}

	// Check for coincident nodes
	if (AreNodesEqual(node11, node12)) {
		_EXCEPTIONT("Coincident nodes used to define edge");
	}
	if (AreNodesEqual(node21, node22)) {
		_EXCEPTIONT("Coincident nodes used to define edge");
	}

	// Clear the intersection vector
	nodeIntersections.clear();

	// Both edges are great circle arcs
	if ((typeFirst  == Edge::Type_GreatCircleArc) &&
		(typeSecond == Edge::Type_GreatCircleArc)
	) {

		// Cross products
		Node nodeN11xN12(CrossProductX(node11, node12));
		Node nodeN21xN22(CrossProductX(node21, node22));
/*
		FixedPoint fp1 = node21.fy * node22.fz;
		FixedPoint fp2 = node21.fz * node22.fy;

		printf("NOW:\n");
		FixedPoint fp3 = fp1 - fp2;

		nodeN11xN12.PrintMX();
		node21.PrintMX();
		node22.PrintMX();
		fp1.Print(); printf("\n");
		fp2.Print(); printf("\n");
		fp3.Print(); printf("\n");
		nodeN21xN22.PrintMX();
*/
		// Check for coincident lines
		FixedPoint fpDot11 = DotProductX(nodeN21xN22, node11);
		FixedPoint fpDot12 = DotProductX(nodeN21xN22, node12);
		FixedPoint fpDot21 = DotProductX(nodeN11xN12, node21);
		FixedPoint fpDot22 = DotProductX(nodeN11xN12, node22);

		// A line which is coincident with both planes
		Node nodeLine;

		// Determine if either fpDot1 or fpDot2 are zero
		bool fp11_isZero = fpDot11.IsZero();
		bool fp12_isZero = fpDot12.IsZero();
		bool fp21_isZero = fpDot21.IsZero();
		bool fp22_isZero = fpDot22.IsZero();
/*
		printf("B: %i %i %i %i\n",
			fp11_isZero,
			fp12_isZero,
			fp21_isZero,
			fp22_isZero);
*/
		// Coincident planes
		if (fp11_isZero && fp21_isZero) {
			return true;

		// node11 is coplanar with the second arc
		} else if (fp11_isZero) {
			nodeLine = node11;

		// node12 is coplanar with the second arc
		} else if (fp12_isZero) {
			nodeLine = node12;

		// node21 is coplanar with the first arc
		} else if (fp21_isZero) {
			nodeLine = node21;

		// node22 is coplanar with the first arc
		} else if (fp22_isZero) {
			nodeLine = node22;

		// Line of intersection is the cross product of cross products
		} else {
			nodeLine = CrossProductX(nodeN11xN12, nodeN21xN22);
/*
			printf("X: "); nodeLine.fx.Print(); printf("\n");
			printf("Y: "); nodeLine.fy.Print(); printf("\n");
			printf("Z: "); nodeLine.fz.Print(); printf("\n");
*/
/*
			// Verify coplanarity
			FixedPoint fpDotDebug1 = DotProductX(nodeLine, nodeN11xN12);
			FixedPoint fpDotDebug2 = DotProductX(nodeLine, nodeN21xN22);
			FixedPoint fpDotDebug3 = DotProductX(nodeN11xN12, nodeLine);
			FixedPoint fpDotDebug4 = DotProductX(nodeN21xN22, nodeLine);

			printf("DD1: "); fpDotDebug1.Print(); printf("\n");
			printf("DD2: "); fpDotDebug2.Print(); printf("\n");
			printf("DD3: "); fpDotDebug3.Print(); printf("\n");
			printf("DD4: "); fpDotDebug4.Print(); printf("\n");

			if (!fpDotDebug1.IsZero() || !fpDotDebug2.IsZero()) {
				_EXCEPTIONT("Logic error");
			}
			if (!fpDotDebug3.IsZero() || !fpDotDebug4.IsZero()) {
				_EXCEPTIONT("Logic error");
			}
*/
		}

		// True if nodeLine is in the fan of node11 and node12.  False if
		// -nodeLine is in the fan of node11 and node12.
		bool fFirstPosModeInRange = true;
		bool fSecondPosModeInRange = true;

		FixedPoint fpDenom;
		FixedPoint fpC;
		FixedPoint fpD;

		// Determine if nodeLine is in the fan of [node11, node12]
		if (! nodeN11xN12.fx.IsZero()) {
			fpDenom = nodeN11xN12.fx;

			fpC = nodeLine.fy * node12.fz - nodeLine.fz * node12.fy;
			fpD = nodeLine.fz * node11.fy - nodeLine.fy * node11.fz;

		} else if (! nodeN11xN12.fy.IsZero()) {
			fpDenom = nodeN11xN12.fy;

			fpC = nodeLine.fz * node12.fx - nodeLine.fx * node12.fz;
			fpD = nodeLine.fx * node11.fz - nodeLine.fz * node11.fx;

		} else if (! nodeN11xN12.fz.IsZero()) {
			fpDenom = nodeN11xN12.fz;

			fpC = nodeLine.fx * node12.fy - nodeLine.fy * node12.fx;
			fpD = nodeLine.fy * node11.fx - nodeLine.fx * node11.fy;

		} else {
			_EXCEPTIONT("Zero Cross product detected");
		}

		if (fpC.IsNegative() && fpD.IsPositive()) {
			return false;
		} else if (fpC.IsPositive() && fpD.IsNegative()) {
			return false;
		} else if (fpC.IsNonNegative() && fpD.IsNonNegative()) {
			fFirstPosModeInRange = fpDenom.IsPositive();
		} else {
			fFirstPosModeInRange = fpDenom.IsNegative();
		}
/*
		fpC.Print(); printf("\n");
		fpD.Print(); printf("\n");
		fpDenom.Print(); printf("\n");
*/
		// If we have reached this point then nodeLine is in the fan
		// of [node11, node12].  Now determine if nodeLine is in the
		// fan of [node21, node22].
		if (! nodeN21xN22.fx.IsZero()) {
			fpDenom = nodeN21xN22.fx;

			fpC = nodeLine.fy * node22.fz - nodeLine.fz * node22.fy;
			fpD = nodeLine.fz * node21.fy - nodeLine.fy * node21.fz;

		} else if (! nodeN21xN22.fy.IsZero()) {
			fpDenom = nodeN21xN22.fy;

			fpC = nodeLine.fz * node22.fx - nodeLine.fx * node22.fz;
			fpD = nodeLine.fx * node21.fz - nodeLine.fz * node21.fx;

		} else if (! nodeN21xN22.fz.IsZero()) {
			fpDenom = nodeN21xN22.fz;

			fpC = nodeLine.fx * node22.fy - nodeLine.fy * node22.fx;
			fpD = nodeLine.fy * node21.fx - nodeLine.fx * node21.fy;

		} else {
			_EXCEPTIONT("Zero Cross product detected");
		}
/*
		fpC.Print(); printf("\n");
		fpD.Print(); printf("\n");
		fpDenom.Print(); printf("\n");
*/
		if (fpC.IsNegative() && fpD.IsPositive()) {
			return false;
		} else if (fpC.IsPositive() && fpD.IsNegative()) {
			return false;
		} else if (fpC.IsNonNegative() && fpD.IsNonNegative()) {
			fSecondPosModeInRange = fpDenom.IsPositive();
		} else {
			fSecondPosModeInRange = fpDenom.IsNegative();
		}

		// Verify fans are not antipodal
		if (fFirstPosModeInRange != fSecondPosModeInRange) {
			return false;
		}

		// Solution exists
		if (!fFirstPosModeInRange) {
			nodeLine.fx.Negate();
			nodeLine.fy.Negate();
			nodeLine.fz.Negate();
		}

		// Verify coplanarity
		FixedPoint fpDotDebug1 = DotProductX(nodeLine, nodeN11xN12);
		FixedPoint fpDotDebug2 = DotProductX(nodeLine, nodeN21xN22);
		FixedPoint fpDotDebug3 = DotProductX(nodeN11xN12, nodeLine);
		FixedPoint fpDotDebug4 = DotProductX(nodeN21xN22, nodeLine);
/*
		printf("DD1: "); fpDotDebug1.Print(); printf("\n");
		printf("DD2: "); fpDotDebug2.Print(); printf("\n");
		printf("DD3: "); fpDotDebug3.Print(); printf("\n");
		printf("DD4: "); fpDotDebug4.Print(); printf("\n");
*/
		if (!fpDotDebug1.IsZero() || !fpDotDebug2.IsZero()) {
			_EXCEPTIONT("Logic error");
		}
		if (!fpDotDebug3.IsZero() || !fpDotDebug4.IsZero()) {
			_EXCEPTIONT("Logic error");
		}

		// Insert solution
		nodeIntersections.push_back(nodeLine);
		return false;

/*
		printf("---\n");
		FixedPoint fpa = nodeN21xN22.fx * node11.fx;
		FixedPoint fpb = nodeN21xN22.fy * node11.fy;
	!	FixedPoint fpc = nodeN21xN22.fz * node11.fz;

		printf("X: %1.15e : ", nodeN11xN12.x); nodeN11xN12.fx.Print(); printf("\n");
		printf("Y: %1.15e : ", nodeN11xN12.y); nodeN11xN12.fy.Print(); printf("\n");
		printf("Z: %1.15e : ", nodeN11xN12.z); nodeN11xN12.fz.Print(); printf("\n");

		printf("X: %1.15e : ", nodeN21xN22.x); nodeN21xN22.fx.Print(); printf("\n");
		printf("Y: %1.15e : ", nodeN21xN22.y); nodeN21xN22.fy.Print(); printf("\n");
		printf("Z: %1.15e : ", nodeN21xN22.z); nodeN21xN22.fz.Print(); printf("\n");

		printf("X: %1.15e : ", node11.x); node11.fx.Print(); printf("\n");
		printf("Y: %1.15e : ", node11.y); node11.fy.Print(); printf("\n");
		printf("Z: %1.15e : ", node11.z); node11.fz.Print(); printf("\n");

		printf("I1: %1.15e : ", dDot1); fp1.Print(); printf("\n");
		printf("I2: %1.15e : ", dDot2); fp2.Print(); printf("\n");
*/

	}

	_EXCEPTIONT("Not implemented");
}

///////////////////////////////////////////////////////////////////////////////

int MeshUtilitiesExact::FindFaceNearNode(
	const Mesh & mesh,
	int ixNode,
	const Node & nodeEnd,
	const Edge::Type edgetype
) {

	// Get the reference point
	const Node & nodeBegin = mesh.nodes[ixNode];

	// Get the set of faces adjacent this node
	const std::set<int> & setNearbyFaces = mesh.revnodearray[ixNode];

	if (setNearbyFaces.size() < 3) {
		_EXCEPTIONT("Insufficient Faces at Corner; at least three Faces expected");
	}

	// Direction along second curve projected onto
	// surface of the sphere at nodeBegin (up to a scaling factor)
	FixedPoint fpDotNbNb = DotProductX(nodeBegin, nodeBegin);
	FixedPoint fpDotNeNb = DotProductX(nodeEnd, nodeBegin);

	Node nodeLocalE =
		ScalarProductX(fpDotNbNb, nodeEnd)
		- ScalarProductX(fpDotNeNb, nodeBegin);
/*
	printf("B := Vector("); nodeBegin.PrintMX(); printf(");\n");
	printf("E := Vector("); nodeEnd.PrintMX(); printf(");\n");
	printf("P := Vector("); nodeLocalE.PrintMX(); printf(");\n");

	printf("BBNorm: "); fpDotNbNb.Print(); printf("\n");
	printf("BENorm: "); fpDotNeNb.Print(); printf("\n");
*/
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

		if (!AreNodesEqual(node1, nodeBegin)) {
			_EXCEPTIONT("Logic error");
		}

		// Local left and right vectors (up to a scaling factor)
		FixedPoint fpDotN0Nb = DotProductX(node0, nodeBegin);
		FixedPoint fpDotN2Nb = DotProductX(node2, nodeBegin);
/*
		printf("n0 := Vector("); node0.PrintMX(); printf(");\n");
		printf("n1 := Vector("); node1.PrintMX(); printf(");\n");
		printf("n2 := Vector("); node2.PrintMX(); printf(");\n");

		printf("N0Nb: "); fpDotN0Nb.Print(); printf("\n");
		printf("N2Nb: "); fpDotN2Nb.Print(); printf("\n");
*/
		Node nodeLocal0 =
			ScalarProductX(fpDotNbNb, node0)
			- ScalarProductX(fpDotN0Nb, nodeBegin);

		Node nodeLocal2 =
			ScalarProductX(fpDotNbNb, node2)
			- ScalarProductX(fpDotN2Nb, nodeBegin);

		Node nodeLocalCross = CrossProductX(nodeLocal0, nodeLocal2);

/*
		printf("l0 := Vector("); nodeLocal0.PrintMX(); printf(");\n");
		printf("l2 := Vector("); nodeLocal2.PrintMX(); printf(");\n");
		printf("nx := Vector("); nodeLocalCross.PrintMX(); printf(");\n");
*/
/*
		// Verify coplanarity
		FixedPoint fpCoplanar = DotProductX(nodeLocalCross, nodeLocalE);

		if (! fpCoplanar.IsZero()) {
			printf("Coplanar: "); fpCoplanar.Print(); printf("\n");
			_EXCEPTIONT("Coplanarity failure");
		}
*/
		// Determine if nodeLocalE is in the fan of [nodeLocal0, nodeLocal2]
		FixedPoint fpDenom;
		FixedPoint fp2;
		FixedPoint fp0;
		if (! nodeLocalCross.fx.IsZero()) {
			fpDenom = nodeLocalCross.fx;

			fp2 = nodeLocalE.fy * nodeLocal2.fz - nodeLocalE.fz * nodeLocal2.fy;
			fp0 = nodeLocalE.fz * nodeLocal0.fy - nodeLocalE.fy * nodeLocal0.fz;

		} else if (! nodeLocalCross.fy.IsZero()) {
			fpDenom = nodeLocalCross.fy;

			fp2 = nodeLocalE.fz * nodeLocal2.fx - nodeLocalE.fx * nodeLocal2.fz;
			fp0 = nodeLocalE.fx * nodeLocal0.fz - nodeLocalE.fz * nodeLocal0.fx;

		} else if (! nodeLocalCross.fz.IsZero()) {
			fpDenom = nodeLocalCross.fz;

			fp2 = nodeLocalE.fx * nodeLocal2.fy - nodeLocalE.fy * nodeLocal2.fx;
			fp0 = nodeLocalE.fy * nodeLocal0.fx - nodeLocalE.fx * nodeLocal0.fy;

		} else {
			_EXCEPTIONT("Zero Cross product detected");
		}

		if (fpDenom.IsPositive()) {
			if (fp2.IsNonNegative() && fp0.IsPositive()) {
				return (*iter);
			}

		} else {
			if (fp2.IsNonPositive() && fp0.IsNegative()) {
				return (*iter);
			}
		}
	}

	_EXCEPTIONT("Logic Error: No exit Face found from Node");
}

///////////////////////////////////////////////////////////////////////////////

int MeshUtilitiesExact::FindFaceNearNode(
	const Mesh & mesh,
	const Node & nodeBegin,
	const Node & nodeEnd,
	const Edge::Type edgetype,
	const FindFaceStruct & aFindFaceStruct
) {

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

			FixedPoint fpOrientation =
				DotProductX(node1 - node0, nodeEnd - nodeBegin);

			// Perpendicular or zero-length edges
			if (fpOrientation.IsZero()) {
				_EXCEPTIONT("Logic error");

			// Oriented along with face0
			} else if (fpOrientation.IsPositive()) {
				return aFindFaceStruct.vecFaceIndices[0];

			// Oriented along with face1
			} else {
				return aFindFaceStruct.vecFaceIndices[1];
			}
		}

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
			Node nodeNormal(CrossProductX(node0, node1));

			// Calculate alignment of outward normal and direction vector
			// at nodeBegin.
			FixedPoint fpDotNbNn = DotProductX(nodeBegin, nodeNormal);
			FixedPoint fpDotNeNn = DotProductX(nodeEnd, nodeNormal);
			FixedPoint fpDotNbNb = DotProductX(nodeBegin, nodeBegin);
			FixedPoint fpDotNeNb = DotProductX(nodeEnd, nodeBegin);

			FixedPoint fpAlign = fpDotNeNn * fpDotNbNb - fpDotNeNb * fpDotNbNn;

			if (fpAlign.IsPositive()) {
				return aFindFaceStruct.vecFaceIndices[0];
			} else if (fpAlign.IsNegative()) {
				return aFindFaceStruct.vecFaceIndices[1];
			} else {
				_EXCEPTIONT("Logic error");
			}
		}

		// Great circle arc originating on a line of constant latitude
		if ((edge0.type == Edge::Type_ConstantLatitude) &&
			(edgetype == Edge::Type_GreatCircleArc)
		) {
			_EXCEPTIONT("Unimplemented");
		}

		// Line of constant latitude originating on a great circle arc
		if ((edgetype == Edge::Type_ConstantLatitude) &&
			(edge0.type == Edge::Type_GreatCircleArc)
		) {
			_EXCEPTIONT("Unimplemented");
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
	}

	// This function does not handle Interior or Exterior nodes
	_EXCEPTIONT("Invalid Node location");
}

///////////////////////////////////////////////////////////////////////////////

#endif

