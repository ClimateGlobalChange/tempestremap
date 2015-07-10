///////////////////////////////////////////////////////////////////////////////
///
///	\file    LinearRemapFV.h
///	\author  Paul Ullrich
///	\version September 1, 2014
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

#ifndef _LINEARREMAPFV_H_
#define _LINEARREMAPFV_H_

#include "DataMatrix3D.h"

///////////////////////////////////////////////////////////////////////////////

class Mesh;
class OfflineMap;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Generate the OfflineMap for remapping from finite volumes to finite
///		volumes.
///	</summary>
void LinearRemapFVtoFV(
	const Mesh & meshInput,
	const Mesh & meshOutput,
	const Mesh & meshOverlap,
	int nOrder,
	OfflineMap & mapRemap
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Generate the OfflineMap for remapping from finite volumes to finite
///		elements using simple sampling of the FV reconstruction.
///	</summary>
void LinearRemapFVtoGLL_Simple(
	const Mesh & meshInput,
	const Mesh & meshOutput,
	const Mesh & meshOverlap,
	const DataMatrix3D<int> & dataGLLNodes,
	const DataMatrix3D<double> & dataGLLJacobian,
	const DataVector<double> & dataGLLNodalArea,
	int nOrder,
	OfflineMap & mapRemap,
	int nMonotoneType,
	bool fContinuous
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Generate the OfflineMap for remapping from finite volumes to finite
///		elements using a new experimental method.
///	</summary>
void LinearRemapFVtoGLL_Volumetric(
	const Mesh & meshInput,
	const Mesh & meshOutput,
	const Mesh & meshOverlap,
	const DataMatrix3D<int> & dataGLLNodes,
	const DataMatrix3D<double> & dataGLLJacobian,
	const DataVector<double> & dataGLLNodalArea,
	int nOrder,
	OfflineMap & mapRemap,
	int nMonotoneType,
	bool fContinuous
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Generate the OfflineMap for remapping from finite volumes to finite
///		elements.
///	</summary>
void LinearRemapFVtoGLL(
	const Mesh & meshInput,
	const Mesh & meshOutput,
	const Mesh & meshOverlap,
	const DataMatrix3D<int> & dataGLLNodes,
	const DataMatrix3D<double> & dataGLLJacobian,
	const DataVector<double> & dataGLLNodalArea,
	int nOrder,
	OfflineMap & mapRemap,
	int nMonotoneType,
	bool fContinuous
);

///////////////////////////////////////////////////////////////////////////////
/*
///	<summary>
///		Generate the OfflineMap for remapping from finite elements to finite
///		elements.
///	</summary>
void LinearRemapGLLtoGLL(
	const Mesh & meshInput,
	const Mesh & meshOutput,
	const Mesh & meshOverlap,
	const DataMatrix3D<int> & dataGLLNodesIn,
	const DataMatrix3D<double> & dataGLLJacobianIn,
	const DataMatrix3D<int> & dataGLLNodesOut,
	const DataMatrix3D<double> & dataGLLJacobianOut,
	const DataVector<double> & dataNodalAreaOut,
	int nPin,
	int nPout,
	OfflineMap & mapRemap,
	int nMonotoneType,
	bool fContinuousIn,
	bool fContinuousOut
);
*/
///////////////////////////////////////////////////////////////////////////////

#endif

