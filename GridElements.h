///////////////////////////////////////////////////////////////////////////////
///
///	\file    GridElements.h
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

#ifndef _GRIDELEMENTS_H_
#define _GRIDELEMENTS_H_

///////////////////////////////////////////////////////////////////////////////

#include <vector>
#include <set>
#include <map>
#include <string>
#include <cmath>

#include "Exception.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A single point in 3D Cartesian geometry.
///	</summary>
class Node {

public:
	///	<summary>
	///		Cartesian coordinates (x,y,z) of this node.
	///	</summary>
	double x;
	double y;
	double z;

public:
	///	<summary>
	///		Default constructor.
	///	</summary>
	Node() :
		x(0.0),
		y(0.0),
		z(0.0)
	{ }

	///	<summary>
	///		Constructor.
	///	</summary>
	Node(
		double _x,
		double _y,
		double _z
	) :
		x(_x),
		y(_y),
		z(_z)
	{ }

	///	<summary>
	///		Equality operator using floating point tolerance.
	///	</summary>
	bool operator== (const Node & node) const {
		static const double Tolerance = 1.0e-12;

		if ((fabs(x - node.x) < Tolerance) &&
			(fabs(y - node.y) < Tolerance) &&
			(fabs(z - node.z) < Tolerance)
		) {
			return true;
		}
		return false;
	}

	///	<summary>
	///		Magnitude of this node.
	///	</summary>
	double Magnitude() const {
		return sqrt(x * x + y * y + z * z);
	}
};

///	<summary>
///		A vector for the storage of Nodes.
///	</summary>
typedef std::vector<Node> NodeVector;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A node index.
///	</summary>
typedef int NodeIndex;

///	<summary>
///		A vector for the storage of Node indices.
///	</summary>
typedef std::vector<NodeIndex> NodeIndexVector;

///	<summary>
///		An index indicating this Node is invalid.
///	</summary>
static const NodeIndex InvalidNode = (-1);

///	<summary>
///		An index indicating this Face is invalid.
///	</summary>
static const NodeIndex InvalidFace = (-1);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		An edge connects two nodes.
///	</summary>
class Edge {

public:
	///	<summary>
	///		Type of edge.
	///	</summary>
	enum Type {
		Type_GreatCircleArc = 0,
		Type_Default = Type_GreatCircleArc,
		Type_ConstantLatitude = 1
	};

public:
	///	<summary>
	///		Node indices representing the endpoints of this edge.
	///	</summary>
	int node[2];

	///	<summary>
	///		The type of this edge.
	///	</summary>
	Type type;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	Edge(
		int node0 = InvalidNode,
		int node1 = InvalidNode,
		Type _type = Type_Default
	) {
		node[0] = node0;
		node[1] = node1;
		type = _type;
	}

	///	<summary>
	///		Virtual destructor.
	///	</summary>
	virtual ~Edge()
	{ }

	///	<summary>
	///		Flip the order of the nodes stored in the segment.  Note that this
	///		does not affect the comparator properties of the segment, and so
	///		this method can be treated as const.
	///	</summary>
	void Flip() const {
		int ixTemp = node[0];
		const_cast<int&>(node[0]) = node[1];
		const_cast<int&>(node[1]) = ixTemp;
	}

	///	<summary>
	///		Accessor.
	///	</summary>
	int operator[](int i) const {
		return node[i];
	}

	int & operator[](int i) {
		return node[i];
	}

	///	<summary>
	///		Get the nodes as an ordered pair.
	///	</summary>
	void GetOrderedNodes(
		int & ixNodeSmall,
		int & ixNodeBig
	) const {
		if (node[0] < node[1]) {
			ixNodeSmall = node[0];
			ixNodeBig   = node[1];
		} else {
			ixNodeSmall = node[1];
			ixNodeBig   = node[0];
		}
	}

	///	<summary>
	///		Comparator.
	///	</summary>
	bool operator<(const Edge & edge) const {

		// Order the local nodes
		int ixNodeSmall;
		int ixNodeBig;
		GetOrderedNodes(ixNodeSmall, ixNodeBig);

		// Order the nodes in edge
		int ixEdgeNodeSmall;
		int ixEdgeNodeBig;
		edge.GetOrderedNodes(ixEdgeNodeSmall, ixEdgeNodeBig);

		// Compare
		if (ixNodeSmall < ixEdgeNodeSmall) {
			return true;
		} else if (ixNodeSmall > ixEdgeNodeSmall) {
			return false;
		} else if (ixNodeBig < ixEdgeNodeBig) {
			return true;
		} else {
			return false;
		}
	}

	///	<summary>
	///		Equality operator.
	///	</summary>
	bool operator==(const Edge & edge) const {
		if ((node[0] == edge.node[0]) &&
			(node[1] == edge.node[1])
		) {
			return true;

		} else if (
			(node[0] == edge.node[1]) &&
			(node[1] == edge.node[0])
		) {
			return true;
		}

		return false;
	}

	///	<summary>
	///		Return the node that is shared between segments.
	///	</summary>
	int CommonNode(
		const Edge & edge
	) const {
		if (edge[0] == node[0]) {
			return node[0];
		} else if (edge[0] == node[1]) {
			return node[1];
		} else if (edge[1] == node[0]) {
			return node[0];
		} else if (edge[1] == node[1]) {
			return node[1];
		} else {
			return InvalidNode;
		}
	}
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		An edge connects two nodes with a sub-array of interior nodes.
///	</summary>
class MultiEdge : public std::vector<int> {

public:
	///	<summary>
	///		Flip the edge.
	///	</summary>
	MultiEdge Flip() const {
		MultiEdge edgeFlip;
		for (int i = size()-1; i >= 0; i--) {
			edgeFlip.push_back((*this)[i]);
		}
		return edgeFlip;
	}

};

typedef std::vector<MultiEdge> MultiEdgeVector;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A pair of face indices, typically on opposite sides of an Edge.
///	</summary>
class FacePair {

public:
	///	<summary>
	///		Indices of the Faces in this pair.
	///	</summary>
	int face[2];

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	FacePair() {
		face[0] = InvalidFace;
		face[1] = InvalidFace;
	}

	///	<summary>
	///		Add a face to this FacePair.
	///	</summary>
	void AddFace(int ixFace) {
		if (face[0] == InvalidFace) {
			face[0] = ixFace;

		} else if (face[1] == InvalidFace) {
			face[1] = ixFace;

		} else {
			_EXCEPTIONT("FacePair already has a full set of Faces.");
		}
	}

	///	<summary>
	///		Does this FacePair have a complete set of Faces?
	///	</summary>
	bool IsComplete() const {
		return ((face[0] != InvalidFace) && (face[1] != InvalidFace));
	}

	///	<summary>
	///		Accessor.
	///	</summary>
	int operator[](int i) const {
		return face[i];
	}
};

///////////////////////////////////////////////////////////////////////////////

typedef std::vector<Edge> EdgeVector;

typedef std::map<Edge, FacePair> EdgeMap;

typedef EdgeMap::value_type EdgeMapPair;

typedef EdgeMap::iterator EdgeMapIterator;

typedef EdgeMap::const_iterator EdgeMapConstIterator;

typedef std::vector<EdgeMap::iterator> EdgeMapIteratorVector;

typedef std::set<Edge> EdgeSet;

typedef std::pair<Edge, FacePair> EdgePair;

typedef std::vector<EdgePair> EdgeMapVector;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A face.
///	</summary>
class Face {

public:
	///	<summary>
	///		Vector of node indices bounding this face, stored in
	///		counter-clockwise order.
	///	</summary>
	EdgeVector edges;

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	Face(
		int edge_count
	) {
		edges.resize(edge_count);
	}

	///	<summary>
	///		Accessor.
	///	</summary>
	inline int operator[](int ix) const {
		return edges[ix][0];
	}

	///	<summary>
	///		Set a node.
	///	</summary>
	void SetNode(int ixLocal, int ixNode) {
		int nEdges = static_cast<int>(edges.size());
		edges[ixLocal][0] = ixNode;

		int ixPrev = (ixLocal + nEdges - 1) % nEdges;
		edges[ixPrev][1] = ixNode;
	}

public:
	///	<summary>
	///		Possible locations of nodes.
	///	</summary>
	enum NodeLocation {
		NodeLocation_Exterior = 0,
		NodeLocation_Default = NodeLocation_Exterior,
		NodeLocation_Interior = 1,
		NodeLocation_Edge = 2,
		NodeLocation_Corner = 3
	};

	///	<summary>
	///		Determine if this face contains the specified point, and whether
	///		the point is along an edge or at a corner.
	///	</summary>
	void ContainsNode(
		const NodeVector & nodevec,
		const Node & node,
		NodeLocation & loc,
		int & ixLocation
	) const;
};

///	<summary>
///		A vector of Faces.
///	</summary>
typedef std::vector<Face> FaceVector;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A reverse node array stores all faces associated with a given node.
///	</summary>
typedef std::vector< std::set<int> > ReverseNodeArray;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A mesh.
///	</summary>
class Mesh {

public:
	///	<summary>
	///		Vector of Nodes for this mesh.
	///	</summary>
	NodeVector nodes;

	///	<summary>
	///		Vector of Faces for this mesh.
	///	<summary>
	FaceVector faces;

	///	<summary>
	///		EdgeMap for this mesh.
	///	</summary>
	EdgeMap edgemap;

	///	<summary>
	///		ReverseNodeArray for this mesh.
	///	</summary>
	ReverseNodeArray revnodearray;

public:
	///	<summary>
	///		Default constructor.
	///	</summary>
	Mesh() {
	}

	///	<summary>
	///		Constructor with input mesh parameter.
	///	</summary>
	Mesh(const std::string & strFile) {
		Read(strFile);
	}

public:
	///	<summary>
	///		Clear the contents of the mesh.
	///	</summary>
	void Clear();

	///	<summary>
	///		Construct the EdgeMap from the NodeVector and FaceVector.
	///	</summary>
	void ConstructEdgeMap();

	///	<summary>
	///		Construct the ReverseNodeArray from the NodeVector and FaceVector.
	///	</summary>
	void ConstructReverseNodeArray();

	///	<summary>
	///		Write the mesh to a NetCDF file.
	///	</summary>
	void Write(const std::string & strFile) const;

	///	<summary>
	///		Read the mesh to a NetCDF file.
	///	</summary>
	void Read(const std::string & strFile);

};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Calculate all intersections between the two edges.
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
	const Edge::Type typeFirst,
	const Node & nodeSecondBegin,
	const Node & nodeSecondEnd,
	const Edge::Type typeSecond,
	std::vector<Node> & nodeIntersections,
	bool fIncludeFirstBeginNode = true
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Calculate the node which is a small increment closer to nodeEnd
///		from nodeBegin.
///	</summary>
void NudgeAlongEdge(
	const Node & nodeBegin,
	const Node & nodeEnd,
	const Edge::Type type,
	Node & nodeNudged
);

///////////////////////////////////////////////////////////////////////////////

#endif

