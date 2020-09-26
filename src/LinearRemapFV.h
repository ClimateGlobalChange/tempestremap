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

#include "DataArray1D.h"
#include "DataArray2D.h"
#include "DataArray3D.h"

///////////////////////////////////////////////////////////////////////////////

class Mesh;
class OfflineMap;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Generate the OfflineMap for remapping from finite volumes to finite
///		volumes using a constant reconstruction.
///	</summary>
void LinearRemapFVtoFV_np1(
	const Mesh & meshInput,
	const Mesh & meshOutput,
	const Mesh & meshOverlap,
	OfflineMap & mapRemap,
	bool verbose = false
);

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
	OfflineMap & mapRemap,
	bool verbose = false
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
	const DataArray3D<int> & dataGLLNodes,
	const DataArray3D<double> & dataGLLJacobian,
	const DataArray1D<double> & dataGLLNodalArea,
	int nOrder,
	OfflineMap & mapRemap,
	int nMonotoneType,
	bool fContinuous,
	bool fNoConservation
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
	const DataArray3D<int> & dataGLLNodes,
	const DataArray3D<double> & dataGLLJacobian,
	const DataArray1D<double> & dataGLLNodalArea,
	int nOrder,
	OfflineMap & mapRemap,
	int nMonotoneType,
	bool fContinuous,
	bool fNoConservation
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
	const DataArray3D<int> & dataGLLNodes,
	const DataArray3D<double> & dataGLLJacobian,
	const DataArray1D<double> & dataGLLNodalArea,
	int nOrder,
	OfflineMap & mapRemap,
	int nMonotoneType,
	bool fContinuous,
	bool fNoConservation
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Generate the OfflineMap for remapping from finite elements to finite
///		elements.
///	</summary>
void LinearRemapGLLtoGLL2(
	const Mesh & meshInput,
	const Mesh & meshOutput,
	const Mesh & meshOverlap,
	const DataArray3D<int> & dataGLLNodesIn,
	const DataArray3D<double> & dataGLLJacobianIn,
	const DataArray3D<int> & dataGLLNodesOut,
	const DataArray3D<double> & dataGLLJacobianOut,
	const DataArray1D<double> & dataNodalAreaOut,
	int nPin,
	int nPout,
	int nMonotoneType,
	bool fContinuousIn,
	bool fContinuousOut,
	bool fNoConservation,
	OfflineMap & mapRemap
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Generate the OfflineMap for remapping from finite elements to finite
///		elements.
///	</summary>
void LinearRemapGLLtoGLL2_Pointwise(
	const Mesh & meshInput,
	const Mesh & meshOutput,
	const Mesh & meshOverlap,
	const DataArray3D<int> & dataGLLNodesIn,
	const DataArray3D<double> & dataGLLJacobianIn,
	const DataArray3D<int> & dataGLLNodesOut,
	const DataArray3D<double> & dataGLLJacobianOut,
	const DataArray1D<double> & dataNodalAreaOut,
	int nPin,
	int nPout,
	int nMonotoneType,
	bool fContinuousIn,
	bool fContinuousOut,
	OfflineMap & mapRemap
);

///////////////////////////////////////////////////////////////////////////////

#endif

