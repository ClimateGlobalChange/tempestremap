///////////////////////////////////////////////////////////////////////////////
///
///	\file    MeshUtilitiesExact.h
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

#ifndef _MESHUTILITIESEXACT_H_
#define _MESHUTILITIESEXACT_H_

#include "Defines.h"
#include "GridElements.h"
#include "GridElementsExact.h"
#include "MeshUtilities.h"

#include <vector>

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Various implementations of methods for determining Faces from Nodes.
///	</summary>
class MeshUtilitiesExact : public MeshUtilities {

public:
	///	<summary>
	///		Overwrite the fuzzy coordinate values with exact values.
	///	</summary>
	inline Node ToRealCoords(
		NodeExact & node
	) {
		Node nodeReal;

		nodeReal.x = node.fx.ToReal();
		nodeReal.y = node.fy.ToReal();
		nodeReal.z = node.fz.ToReal();

		Real dMag = nodeReal.Magnitude();

		if (dMag == 0.0) {
			node.fx.Print(); printf("\n");
			node.fy.Print(); printf("\n");
			node.fz.Print(); printf("\n");
			_EXCEPTIONT("Zero magnitude Node detected");
		}

		nodeReal.x /= dMag;
		nodeReal.y /= dMag;
		nodeReal.z /= dMag;

		return nodeReal;
	}

	///	<summary>
	///		Determine if two Nodes are equal.
	///	</summary>
	bool AreNodesEqual(
		const NodeExact & node0,
		const NodeExact & node1
	);

	///	<summary>
	///		Determine if face contains node, and whether
	///		the Node is along an edge or at a corner.
	///	</summary>
	virtual void ContainsNode(
		const Face & face,
		const NodeVector & nodevec,
		const Node & node,
		Face::NodeLocation & loc,
		int & ixLocation
	) const;

	///	<summary>
	///		Calculate all intersections between the Edge connecting
	///		nodeFirstBegin and nodeFirstEnd with type typeFirst and the Edge
	///		connecting nodeSecondBegin and nodeSecondEnd with type typeSecond.
	///		Intersections are recorded in nodeIntersections.
	///	</summary>
	///	<returns>
	///		Returns true if lines are coincident, false otherwise.
	///
	///		If lines are coincident, intersections includes any nodes of Second
	///		that are contained in First, ordered from FirstBegin to FirstEnd.
	///	</returns>
	bool CalculateEdgeIntersections(
		const Node & nodeFirstBegin,
		const Node & nodeFirstEnd,
		Edge::Type typeFirst,
		const Node & nodeSecondBegin,
		const Node & nodeSecondEnd,
		Edge::Type typeSecond,
		std::vector<NodeExact> & nodeIntersections,
		bool fIncludeFirstBeginNode = false
	);

	///	<summary>
	///		Find the Face that is near ixNode in the direction of nodeEnd.
	///	</summary>
	int FindFaceNearNode(
		const Mesh & mesh,
		int ixNode,
		const NodeExact & nodeEnd,
		const Edge::Type edgetype
	);

	///	<summary>
	///		Find the Face that is near nodeBegin in the direction of nodeEnd.
	///	</summary>
	int FindFaceNearNode(
		const Mesh & mesh,
		const NodeExact & nodeBegin,
		const NodeExact & nodeEnd,
		const Edge::Type edgetype,
		const FindFaceStruct & aFindFaceStruct
	);
};

///////////////////////////////////////////////////////////////////////////////

#endif

