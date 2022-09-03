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
inline Node GetFaceCentroid(
	const Face & face,
	const NodeVector & nodes
) {
	Node nodeRef;

	nodeRef.x = 0.0;
	nodeRef.y = 0.0;
	nodeRef.z = 0.0;

	for (int i = 0; i < face.edges.size(); i++) {
		nodeRef.x += nodes[face[i]].x;
		nodeRef.y += nodes[face[i]].y;
		nodeRef.z += nodes[face[i]].z;
	}
	nodeRef.x /= static_cast<double>(face.edges.size());
	nodeRef.y /= static_cast<double>(face.edges.size());
	nodeRef.z /= static_cast<double>(face.edges.size());

	return nodeRef;
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Get a perpendicular tangent basis to the given unit-length node.
///	</summary>
inline void GetTangentBasis(
	const Node & nodeRef,
	Node & nodeA,
	Node & nodeB
) {
	// One basis vector has zero y component
	if (fabs(nodeRef.z) > 0.5) {
		double dMagXZ = sqrt(nodeRef.x * nodeRef.x + nodeRef.z * nodeRef.z);
		double dInvMagXZ = 1.0 / dMagXZ;

		nodeA.x = - nodeRef.z * dInvMagXZ;
		nodeA.y = 0.0;
		nodeA.z = nodeRef.x * dInvMagXZ;

		nodeB.x = nodeRef.y * nodeRef.x * dInvMagXZ;
		nodeB.y = - dMagXZ;
		nodeB.z = nodeRef.y * nodeRef.z * dInvMagXZ;

	// One basis vector has zero z component
	} else {
		double dMagXY = sqrt(nodeRef.x * nodeRef.x + nodeRef.y * nodeRef.y);
		double dInvMagXY = 1.0 / dMagXY;

		nodeA.x = - nodeRef.y * dInvMagXY;
		nodeA.y = nodeRef.x * dInvMagXY;
		nodeA.z = 0.0;

		nodeB.x = - nodeRef.z * nodeRef.x * dInvMagXY;
		nodeB.y = - nodeRef.z * nodeRef.y * dInvMagXY;
		nodeB.z = dMagXY;
	}
}

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
///		Find the triangle that contains a two-dimensional point
///	</summary>

void GetTriangleThatContainsPoint(
	const Mesh & mesh,
	int iFaceInitial,
	int & iFaceFinal,
	double dX,
	double dY
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Find the triangle that contains a two-dimensional point
///	</summary>

void GetFaceThatContainsPoint(
	const Mesh & mesh,
	int iFaceInitial,
	int & iFaceFinal,
	double dX,
	double dY,
	double dZ
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Returns true if a face contains a given point
///	</summary>

bool fFaceContainsPoint(
	const Mesh & mesh,
	int iFace,
	double dX,
	double dY,
	double dZ

);

///	<summary>
///		Returns true if a face contains a given point
///	</summary>

///////////////////////////////////////////////////////////////////////////////

bool fFaceContainsPoint(
	const NodeVector & nodesP,
	double dX,
	double dY,
	double dZ

);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Barycentric coordinates of a two-dimensional point
///	</summary>

void BarycentricCoordinates(
	const Mesh & mesh,
	int iFace,
	double dX,
	double dY,
	double & dA,
	double & dB
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Determine if a triangle contains a point
///	</summary>

bool fTriangleContainsPoint(
	const Mesh & mesh,
	int iFace,
	double dX,
	double dY
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
///		Multiply a 3-dimensional vector by a 3x3 matrix.
///	</summary>

void MatVectorMult(
	const DataArray2D<double> & dMat,
	DataArray1D<double> & dRHS,
	DataArray1D<double> & dOutput
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Intersection of a point and triangle.
///	</summary>

void TriangleLineIntersection(
	Node & nodeQ,
	NodeVector & nodesP,
	DataArray1D<double> & dCoeffs,
	double & dCond
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Find the coordinates of a point in the basis of a quadrilateral.
///	</summary>

void NewtonQuadrilateral(
	Node & nodeQ,
	NodeVector & nodesP,
	DataArray1D<double> & dCoeffs,
	bool & fConverged
);

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Build the integration array, an operator that integrates a polynomial
///		reconstruction (specified as a vector of polynomial coefficients) to
///		the integral over all overlap faces.
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
///		Build the fit array, which maps the coefficients of the polynomial
///		reconstruction to the area averages over the AdjacentFaceVector.
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

///	<summary>
///		Construct the dual mesh.
///	</summary>

void Dual(
	Mesh & mesh
);

///////////////////////////////////////////////////////////////////////////////

#endif // _FINITEVOLUMETOOLS_H_

