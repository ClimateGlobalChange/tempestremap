///////////////////////////////////////////////////////////////////////////////
///
///	\file    OverlapMesh.h
///	\author  Paul Ullrich
///	\version August 19, 2021
///
///	<remarks>
///		Copyright 2021 Paul Ullrich
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
	OverlapMeshMethod_Mixed
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Generate the mesh obtained by overlapping meshes meshSource and
///		meshTarget.
///	</summary>
void GenerateOverlapMeshEdge(
	const Mesh & meshSource,
	const Mesh & meshTarget,
	Mesh & meshOverlap,
    OverlapMeshMethod method,
    const bool fVerbose = true
);

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
void GenerateOverlapMeshKdx(
	const Mesh & meshSource,
	const Mesh & meshTarget,
	Mesh & meshOverlap,
    OverlapMeshMethod method,
	const bool fAllowNoOverlap,
    const bool fVerbose = true
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Generate the mesh obtained by overlapping meshes meshSource and
///		meshTarget.
///	</summary>
void GenerateOverlapMeshLint(
	const Mesh & meshSource,
	const Mesh & meshTarget,
	Mesh & meshOverlap,
    OverlapMeshMethod method,
	const bool fParallel,
	std::string strTempDir,
    const bool fVerbose = true
);

///////////////////////////////////////////////////////////////////////////////

#endif

