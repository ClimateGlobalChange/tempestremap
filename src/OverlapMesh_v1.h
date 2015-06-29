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
///		Generate the mesh obtained by overlapping meshes meshFirst and
///		meshSecond.
///	</summary>
void GenerateOverlapMesh(
	const Mesh & meshFirst,
	const Mesh & meshSecond,
	Mesh & meshOverlap,
	OverlapMeshMethod method
);

///////////////////////////////////////////////////////////////////////////////

#endif

