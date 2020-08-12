///////////////////////////////////////////////////////////////////////////////
///
///	\file    FiniteVolumeTools.h
///	\author  Paul Ullrich
///	\version August 12, 2020
///
///	<remarks>
///		Copyright 2020 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#ifndef _FINITEVOLUMETOOLS_H_
#define _FINITEVOLUMETOOLS_H_

#include "Defines.h"
#include "GridElements.h"
#include "DataArray1D.h"
#include "DataArray2D.h"
#include "TriangularQuadrature.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Face index and distance metric pair.
///	</summary>
typedef std::pair<int, int> FaceDistancePair;

///	<summary>
///		Vector storing adjacent Faces.
///	</summary>
typedef std::vector<FaceDistancePair> AdjacentFaceVector;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Get the centroid of a Face by averaging vertices in Cartesian space.
///	</summary>
Node GetFaceCentroid(
	const Face & face,
	const NodeVector & nodes
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Get a vector of adjacent faces using edge information to determine
///		adjacencies.
///	</summary>
void GetAdjacentFaceVectorByEdge(
	const Mesh & mesh,
	int iFaceInitial,
	int nRequiredFaceSetSize,
	AdjacentFaceVector & vecFaces
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Get a vector of adjacent faces using vertex information to determine
///		adjacencies.
///	</summary>
void GetAdjacentFaceVectorByNode(
	const Mesh & mesh,
	int iFaceInitial,
	int nRequiredFaceSetSize,
	AdjacentFaceVector & vecFaces
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Build the integration array, an operator that integrates a polynomial
///		reconstruction over all overlap faces.
///	</summary>
void BuildIntegrationArray(
	const Mesh & meshInput,
	const Mesh & meshOverlap,
	const TriangularQuadratureRule & triquadrule,
	int ixFirstFace,
	int ixOverlapBegin,
	int ixOverlapEnd,
	int nOrder,
	DataArray2D<double> & dIntArray
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Build the fit array, which maps area averages in adjacent faces
///		to the coefficients of the polynomial reconstruction.
///	</summary>
void BuildFitArray(
	const Mesh & mesh,
	const TriangularQuadratureRule & triquadrule,
	int ixFirst,
	const AdjacentFaceVector & vecAdjFaces,
	int nOrder,
	int nFitWeightsExponent,
	const DataArray1D<double> & dConstraint,
	DataArray2D<double> & dFitArray,
	DataArray1D<double> & dFitWeights
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Invert the fit array using a pseudo-inverse.
///	</summary>
bool InvertFitArray_Corrected(
	const DataArray1D<double> & dConstraint,
	DataArray2D<double> & dFitArray,
	DataArray1D<double> & dFitWeights,
	DataArray2D<double> & dFitArrayPlus
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Invert the fit array using LAPACK's built-in least squares solver.
///	</summary>
void InvertFitArray_LeastSquares(
	const DataArray1D<double> & dConstraint,
	DataArray2D<double> & dFitArray,
	DataArray1D<double> & dFitWeights,
	DataArray2D<double> & dFitArrayPlus
);

///////////////////////////////////////////////////////////////////////////////

#endif // _FINITEVOLUMETOOLS_H_

