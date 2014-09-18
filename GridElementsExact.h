///////////////////////////////////////////////////////////////////////////////
///
///	\file    GridElementsExact.h
///	\author  Paul Ullrich
///	\version September 17, 2014
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

#ifndef _GRIDELEMENTSEXACT_H_
#define _GRIDELEMENTSEXACT_H_

///////////////////////////////////////////////////////////////////////////////

#include "Defines.h"

#include "GridElements.h"
#include "FixedPoint.h"

#include <vector>
#include <set>
#include <map>
#include <string>
#include <cmath>

#include "Exception.h"
#include "DataVector.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A single point in 3D Cartesian geometry.
///	</summary>
class NodeExact {

public:
	///	<summary>
	///		Fixed point Cartesian coordinates (x,y,z) of this Node.
	///	</summary>
	FixedPoint fx;
	FixedPoint fy;
	FixedPoint fz;

public:
	///	<summary>
	///		Default constructor.
	///	</summary>
	NodeExact() {
		fx.Set(0.0);
		fy.Set(0.0);
		fz.Set(0.0);
	}

	///	<summary>
	///		Constructor.
	///	</summary>
	NodeExact(
		Real _x,
		Real _y,
		Real _z
	) {
		fx.Set(_x);
		fy.Set(_y);
		fz.Set(_z);
	}

	///	<summary>
	///		Copy constructor.
	///	</summary>
	NodeExact(const NodeExact & node) {
		fx = node.fx;
		fy = node.fy;
		fz = node.fz;
	}

	///	<summary>
	///		Assignment operator.
	///	</summary>
	const NodeExact & operator=(const NodeExact & node) {
		fx = node.fx;
		fy = node.fy;
		fz = node.fz;

		return (*this);
	}

	///	<summary>
	///		Copy constructor from Node.
	///	</summary>
	NodeExact(const Node & node) {
		fx.Set(node.x);
		fy.Set(node.y);
		fz.Set(node.z);
	}

	///	<summary>
	///		Assignment operator.
	///	</summary>
	const NodeExact & operator=(const Node & node) {
		fx.Set(node.x);
		fy.Set(node.y);
		fz.Set(node.z);

		return (*this);
	}

	///	<summary>
	///		Difference between two nodes.
	///	</summary>
	NodeExact operator-(const NodeExact & node) const {
		NodeExact nodeDiff;

		nodeDiff.fx = fx - node.fx;
		nodeDiff.fy = fy - node.fy;
		nodeDiff.fz = fz - node.fz;

		return nodeDiff;
	}

	///	<summary>
	///		Casting operator to a Node.
	///	</summary>
	operator Node() const {
		Node node;
		node.x = fx.ToReal();
		node.y = fy.ToReal();
		node.z = fz.ToReal();
		return node;
	}

public:
	///	<summary>
	///		Output node to stdout.
	///	</summary>
	void Print(const char * szName) const {
		printf("%s:\n", szName);
		printf("  X: "); fx.Print(); printf("\n");
		printf("  Y: "); fy.Print(); printf("\n");
		printf("  Z: "); fz.Print(); printf("\n");
	}

	///	<summary>
	///		Output node to stdout.
	///	</summary>
	void PrintMX() const {
		printf("[");
		fx.Print(); printf(", ");
		fy.Print(); printf(", ");
		fz.Print(); printf("]\n");
	}

	///	<summary>
	///		Output node to stdout.
	///	</summary>
	void PrintNorm() const {
		Real dx = fx.ToReal();
		Real dy = fy.ToReal();
		Real dz = fz.ToReal();

		Real mag = sqrt(dx * dx + dy * dy + dz * dz);

		printf("%1.15e %1.15e %1.15e\n", dx / mag, dy / mag, dz / mag);
	}
};

///	<summary>
///		A vector for the storage of Nodes.
///	</summary>
typedef std::vector<NodeExact> NodeExactVector;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Calculate the exact dot product between two Nodes.
///	</summary>
inline FixedPoint DotProductX(
	const NodeExact & node1,
	const NodeExact & node2
) {
	return (node1.fx * node2.fx + node1.fy * node2.fy + node1.fz * node2.fz); 
}

///	<summary>
///		Calculate the exact cross product between two Nodes.
///	</summary>
inline NodeExact CrossProductX(
	const NodeExact & node1,
	const NodeExact & node2
) {
	NodeExact nodeCross;
	nodeCross.fx = node1.fy * node2.fz - node1.fz * node2.fy;
	nodeCross.fy = node1.fz * node2.fx - node1.fx * node2.fz;
	nodeCross.fz = node1.fx * node2.fy - node1.fy * node2.fx;
	return nodeCross;
}

///	<summary>
///		Calculate the product of a Node with a scalar.
///	</summary>
inline NodeExact ScalarProductX(
	const FixedPoint & fp,
	const NodeExact & node
) {
	NodeExact nodeProduct(node);
	nodeProduct.fx *= fp;
	nodeProduct.fy *= fp;
	nodeProduct.fz *= fp;
	return nodeProduct;
}

///////////////////////////////////////////////////////////////////////////////

#endif
