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

///////////////////////////////////////////////////////////////////////////////

bool MeshUtilitiesExact::AreNodesEqual(
	const Node & node0,
	const Node & node1
) {
	if ((node0.x == node1.x) &&
		(node0.y == node1.y) &&
		(node0.z == node1.z)
	) {
		return true;
	}
	return false;

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

		// Check for coincident lines
		FixedPoint fpDot1 = DotProductX(nodeN11xN12, node21);
		FixedPoint fpDot2 = DotProductX(nodeN21xN22, node11);

		// A line which is coincident with both planes
		Node nodeLine;
		FixedPoint fp1 = DotProductX(nodeN11xN12, node21);
		FixedPoint fp2 = DotProductX(nodeN21xN22, node11);

		bool fp1_isZero = fp1.IsZero();
		bool fp2_isZero = fp2.IsZero();

		if (fp1_isZero && fp2_isZero) {
			return true;

		} else if (fp1_isZero) {
			nodeLine = node21;

		} else if (fp2_isZero) {
			nodeLine = node11;

		} else {
			nodeLine = CrossProduct(nodeN11xN12, nodeN21xN22);
/*
			printf("X: "); nodeLine.fx.Print(); printf("\n");
			printf("Y: "); nodeLine.fy.Print(); printf("\n");
			printf("Z: "); nodeLine.fz.Print(); printf("\n");

			// Verify coplanarity
			FixedPoint fpDotDebug1 = DotProduct(nodeLine, nodeN11xN12);
			FixedPoint fpDotDebug2 = DotProduct(nodeLine, nodeN21xN22);

			printf("DD1: "); fpDotDebug1.Print(); printf("\n");
			printf("DD2: "); fpDotDebug2.Print(); printf("\n");

			if (!fpDotDebug1.IsZero() || !fpDotDebug2.IsZero()) {
				_EXCEPTIONT("Logic error");
			}
*/
		}

		// True if nodeLine is in the fan of node11 and node12.  False if
		// -nodeLine is in the fan of node11 and node12.
		bool fPosModeInRange = true;

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
		} else if (fpC.IsPositive() && fpD.IsPositive()) {
			fPosModeInRange = fpDenom.IsPositive();
		} else {
			fPosModeInRange = !fpDenom.IsPositive();
		}

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

		if (fpC.IsNegative() && fpD.IsPositive()) {
			return false;
		} else if (fpC.IsPositive() && fpD.IsNegative()) {
			return false;
		} else if (fpC.IsPositive() && fpD.IsPositive()) {
			if (!fPosModeInRange) {
				return false;
			}
		} else {
			if (fPosModeInRange) {
				return false;
			}
		}

		// Solution exists
		nodeIntersections.push_back(nodeLine);

/*
		printf("---\n");
		FixedPoint fpa = nodeN21xN22.fx * node11.fx;
		FixedPoint fpb = nodeN21xN22.fy * node11.fy;
	!	FixedPoint fpc = nodeN21xN22.fz * node11.fz;

		printf("X: %1.15Le : ", nodeN11xN12.x); nodeN11xN12.fx.Print(); printf("\n");
		printf("Y: %1.15Le : ", nodeN11xN12.y); nodeN11xN12.fy.Print(); printf("\n");
		printf("Z: %1.15Le : ", nodeN11xN12.z); nodeN11xN12.fz.Print(); printf("\n");

		printf("X: %1.15Le : ", nodeN21xN22.x); nodeN21xN22.fx.Print(); printf("\n");
		printf("Y: %1.15Le : ", nodeN21xN22.y); nodeN21xN22.fy.Print(); printf("\n");
		printf("Z: %1.15Le : ", nodeN21xN22.z); nodeN21xN22.fz.Print(); printf("\n");

		printf("X: %1.15Le : ", node11.x); node11.fx.Print(); printf("\n");
		printf("Y: %1.15Le : ", node11.y); node11.fy.Print(); printf("\n");
		printf("Z: %1.15Le : ", node11.z); node11.fz.Print(); printf("\n");

		printf("I1: %1.15Le : ", dDot1); fp1.Print(); printf("\n");
		printf("I2: %1.15Le : ", dDot2); fp2.Print(); printf("\n");
*/

	}

	_EXCEPTIONT("Not implemented");
}

///////////////////////////////////////////////////////////////////////////////

void MeshUtilitiesExact::FindFaceFromNode(
	const Mesh & mesh,
	const Node & node,
	FindFaceStruct & aFindFaceStruct
) {
	_EXCEPTION();
}

///////////////////////////////////////////////////////////////////////////////

int MeshUtilitiesExact::FindFaceNearNode(
	const Mesh & mesh,
	int ixNode,
	const Node & nodeEnd,
	const Edge::Type edgetype
) {
	_EXCEPTION();
	return 0;
}

///////////////////////////////////////////////////////////////////////////////

int MeshUtilitiesExact::FindFaceNearNode(
	const Mesh & mesh,
	const Node & nodeBegin,
	const Node & nodeEnd,
	const Edge::Type edgetype,
	const FindFaceStruct & aFindFaceStruct
) {
	_EXCEPTION();
	return 0;
}

///////////////////////////////////////////////////////////////////////////////

