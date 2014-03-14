///////////////////////////////////////////////////////////////////////////////
///
///	\file    OverlapMesh.h
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

#ifndef _OVERLAPMESH_H_
#define _OVERLAPMESH_H_

#include "GridElements.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Find the face(s) in mesh that contain node.
///	</summary>
void FindFaceFromNode(
	const Mesh & mesh,
	const Node & node,
	std::vector<int> & vecFaceIndices
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Find the face index in mesh that is near nodeBegin and is found by
///		traveling towards nodeEnd along an edge of type edgetype.  Only faces
///		from vecPossibleFaces are considered.
///	</summary>
void FindFaceNearNode(
	const Mesh & mesh,
	const Node & nodeBegin,
	const Node & nodeEnd,
	const Edge::Type edgetype,
	const std::vector<int> & vecPossibleFaces,
	std::vector<int> & vecFaceIndices
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Find the face index in mesh that is near the node with index ixNode
///		and is found by traveling towards nodeEnd along an edge of type
///		edgetype.
///	</summary>
void FindFaceNearNode(
	const Mesh & mesh,
	int ixNode,
	const Node & nodeEnd,
	const Edge::Type edgetype,
	std::vector<int> & vecFaceIndices
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Generate the mesh obtained by overlapping meshes meshFirst and
///		meshSecond.
///	</summary>
void GenerateOverlapMesh(
	const Mesh & meshFirst,
	const Mesh & meshSecond,
	Mesh & meshOverlap
);

///////////////////////////////////////////////////////////////////////////////

#endif

