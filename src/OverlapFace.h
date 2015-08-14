///////////////////////////////////////////////////////////////////////////////
///
///	\file    OverlapFace.h
///	\author  Paul Ullrich
///	\version June 29, 2015
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

#ifndef _OVERLAPFACE_H_
#define _OVERLAPFACE_H_

#include "Defines.h"

#include "GridElements.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Given a Face index in the source mesh, construct the Mesh that
///		overlaps that Face.
///	</summary>
void OverlapFace(
	const Mesh & meshSource,
	const Mesh & meshTarget,
	int iFace,
	Mesh & meshOverlap,
	OverlapMeshMethod method
) {
}

///////////////////////////////////////////////////////////////////////////////


