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

#include "Defines.h"

#include "GridElements.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Method to use to generate overlap mesh
///	</summary>
enum OverlapMeshMethod {
	OverlapMeshMethod_Fuzzy,
	OverlapMeshMethod_Exact,
	OverlapMeshMethod_Mixed,
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Generate the mesh obtained by overlapping meshes meshSource and
///		meshTarget.
///	</summary>
void GenerateOverlapMesh_v1(
	const Mesh & meshSource,
	const Mesh & meshTarget,
	Mesh & meshOverlap,
    OverlapMeshMethod method,
    const bool verbose=true
);

///////////////////////////////////////////////////////////////////////////////
/*
///	<summary>
///		Generate the mesh obtained by overlapping Face iSourceFace in
///		meshSource with meshTarget.
///	</summary>
void GenerateOverlapMeshFromFace(
	const Mesh & meshSource,
	const Mesh & meshTarget,
	int iSourceFace,
	Mesh & meshOverlap,
	OverlapMeshMethod method
);
*/

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Compute the overlap polygon between two Faces.
///	</summary>
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
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Generate the mesh obtained by overlapping meshes meshSource and
///		meshTarget.
///	</summary>
void GenerateOverlapMesh_v2(
	const Mesh & meshSource,
	const Mesh & meshTarget,
	Mesh & meshOverlap,
    OverlapMeshMethod method,
    const bool verbose=true
);

///////////////////////////////////////////////////////////////////////////////

#endif

