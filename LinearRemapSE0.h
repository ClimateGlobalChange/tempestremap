///////////////////////////////////////////////////////////////////////////////
///
///	\file    LinearRemapSE0.h
///	\author  Paul Ullrich
///	\version August 14, 2014
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

#ifndef _LINEARREMAPSE0_H_
#define _LINEARREMAPSE0_H_

///////////////////////////////////////////////////////////////////////////////

#include "GridElements.h"
#include "DataMatrix3D.h"

class OfflineMap;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Generate the OfflineMap for linear conserative element-average
///		spectral element to element average remapping.
///	</summary>
void LinearRemapSE0(
	const Mesh & meshInput,
	const Mesh & meshOutput,
	const Mesh & meshOverlap,
	const DataMatrix3D<int> & dataGLLNodes,
	const DataMatrix3D<double> & dataGLLJacobian,
	OfflineMap & mapRemap
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Generate the OfflineMap for cubic conserative element-average
///		spectral element to element average remapping.
///	</summary>
void LinearRemapSE4(
	const Mesh & meshInput,
	const Mesh & meshOutput,
	const Mesh & meshOverlap,
	const DataMatrix3D<int> & dataGLLNodes,
	const DataMatrix3D<double> & dataGLLJacobian,
	bool fMonotone,
	OfflineMap & mapRemap
);

///////////////////////////////////////////////////////////////////////////////

#endif

