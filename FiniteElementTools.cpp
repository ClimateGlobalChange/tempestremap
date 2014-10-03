///////////////////////////////////////////////////////////////////////////////
///
///	\file    FiniteElementTools.cpp
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

#include "FiniteElementTools.h"
#include "PolynomialInterp.h"
#include "GridElements.h"
#include "GaussLobattoQuadrature.h"

#include <map>

///////////////////////////////////////////////////////////////////////////////

void ApplyLocalMap(
	const Face & face,
	const NodeVector & nodes,
	double dAlpha,
	double dBeta,
	Node & nodeG,
	Node & dDx1G,
	Node & dDx2G
) {
	// Calculate nodal locations on the plane
	double dXc =
		  nodes[face[0]].x * (1.0 - dAlpha) * (1.0 - dBeta)
		+ nodes[face[1]].x *        dAlpha  * (1.0 - dBeta)
		+ nodes[face[2]].x *        dAlpha  *        dBeta
		+ nodes[face[3]].x * (1.0 - dAlpha) *        dBeta;

	double dYc =
		  nodes[face[0]].y * (1.0 - dAlpha) * (1.0 - dBeta)
		+ nodes[face[1]].y *        dAlpha  * (1.0 - dBeta)
		+ nodes[face[2]].y *        dAlpha  *        dBeta
		+ nodes[face[3]].y * (1.0 - dAlpha) *        dBeta;

	double dZc =
		  nodes[face[0]].z * (1.0 - dAlpha) * (1.0 - dBeta)
		+ nodes[face[1]].z *        dAlpha  * (1.0 - dBeta)
		+ nodes[face[2]].z *        dAlpha  *        dBeta
		+ nodes[face[3]].z * (1.0 - dAlpha) *        dBeta;

	double dR = sqrt(dXc * dXc + dYc * dYc + dZc * dZc);

	// Check if this Node exists in the NodeMap
	nodeG.x = dXc / dR;
	nodeG.y = dYc / dR;
	nodeG.z = dZc / dR;

	// Pointwise basis vectors in Cartesian geometry
	Node dDx1F(
		(1.0 - dBeta) * (nodes[face[1]].x - nodes[face[0]].x)
		+      dBeta  * (nodes[face[2]].x - nodes[face[3]].x),
		(1.0 - dBeta) * (nodes[face[1]].y - nodes[face[0]].y)
		+      dBeta  * (nodes[face[2]].y - nodes[face[3]].y),
		(1.0 - dBeta) * (nodes[face[1]].z - nodes[face[0]].z)
		+      dBeta  * (nodes[face[2]].z - nodes[face[3]].z));

	Node dDx2F(
		(1.0 - dAlpha) * (nodes[face[3]].x - nodes[face[0]].x)
		+      dAlpha  * (nodes[face[2]].x - nodes[face[1]].x),
		(1.0 - dAlpha) * (nodes[face[3]].y - nodes[face[0]].y)
		+      dAlpha  * (nodes[face[2]].y - nodes[face[1]].y),
		(1.0 - dAlpha) * (nodes[face[3]].z - nodes[face[0]].z)
		+      dAlpha  * (nodes[face[2]].z - nodes[face[1]].z));

	// Pointwise basis vectors in spherical geometry
	double dDenomTerm = 1.0 / (dR * dR * dR);

	dDx1G = Node(
		- dXc * (dYc * dDx1F.y + dZc * dDx1F.z)
			+ (dYc * dYc + dZc * dZc) * dDx1F.x,
		- dYc * (dXc * dDx1F.x + dZc * dDx1F.z)
			+ (dXc * dXc + dZc * dZc) * dDx1F.y,
		- dZc * (dXc * dDx1F.x + dYc * dDx1F.y)
			+ (dXc * dXc + dYc * dYc) * dDx1F.z);

	dDx2G = Node(
		- dXc * (dYc * dDx2F.y + dZc * dDx2F.z)
			+ (dYc * dYc + dZc * dZc) * dDx2F.x,
		- dYc * (dXc * dDx2F.x + dZc * dDx2F.z)
			+ (dXc * dXc + dZc * dZc) * dDx2F.y,
		- dZc * (dXc * dDx2F.x + dYc * dDx2F.y)
			+ (dXc * dXc + dYc * dYc) * dDx2F.z);

	dDx1G.x *= dDenomTerm;
	dDx1G.y *= dDenomTerm;
	dDx1G.z *= dDenomTerm;

	dDx2G.x *= dDenomTerm;
	dDx2G.y *= dDenomTerm;
	dDx2G.z *= dDenomTerm;
}

///////////////////////////////////////////////////////////////////////////////

void ApplyInverseMap(
	const Face & face,
	const NodeVector & nodes,
	const Node & node,
	double & dAlpha,
	double & dBeta
) {
	// Forward map quantities
	Node nodeG;
	Node nodeDx1G;
	Node nodeDx2G;

	// First guess
	dAlpha = 0.5;
	dBeta = 0.5;

	// Map matrix
	double dMap[2][2];

	double dF[2];

	// Loop until convergence
	double dDelta = 1.0;

	// Try all pairs of 
	for (int i = 0; i < 10; i++) {
/*
		nodes[face[0]].Print("f0");
		nodes[face[1]].Print("f1");
		nodes[face[2]].Print("f2");
		nodes[face[3]].Print("f3");
*/
		// Apply forward map
		ApplyLocalMap(
			face, nodes,
			dAlpha, dBeta,
			nodeG, nodeDx1G, nodeDx2G);

		// Pick the two Cartesian components with greatest chance of success
		const Node & nodeRef = nodes[face[0]];

		if ((fabs(nodeRef.x) >= fabs(nodeRef.y)) &&
			(fabs(nodeRef.x) >= fabs(nodeRef.z))
		) {
			dMap[0][0] = nodeDx1G.y;
			dMap[0][1] = nodeDx2G.y;
			dMap[1][0] = nodeDx1G.z;
			dMap[1][1] = nodeDx2G.z;

			dF[0] = (nodeG.y - node.y);
			dF[1] = (nodeG.z - node.z);

		} else if (
			(fabs(nodeRef.y) >= fabs(nodeRef.x)) &&
			(fabs(nodeRef.y) >= fabs(nodeRef.z))
		) {
			dMap[0][0] = nodeDx1G.x;
			dMap[0][1] = nodeDx2G.x;
			dMap[1][0] = nodeDx1G.z;
			dMap[1][1] = nodeDx2G.z;

			dF[0] = (nodeG.x - node.x);
			dF[1] = (nodeG.z - node.z);

		} else {
			dMap[0][0] = nodeDx1G.x;
			dMap[0][1] = nodeDx2G.x;
			dMap[1][0] = nodeDx1G.y;
			dMap[1][1] = nodeDx2G.y;

			dF[0] = (nodeG.x - node.x);
			dF[1] = (nodeG.y - node.y);
		}

		double dDet = dMap[0][0] * dMap[1][1] - dMap[0][1] * dMap[1][0];

		if (fabs(dDet) < ReferenceTolerance) {
			_EXCEPTIONT("Zero determinant in map inverse");
		}

		// Apply Newton's method
		double dDeltaAlpha =
			1.0 / dDet * (  dMap[1][1] * dF[0] - dMap[0][1] * dF[1]);

		double dDeltaBeta =
			1.0 / dDet * (- dMap[1][0] * dF[0] + dMap[0][0] * dF[1]);

		dAlpha -= dDeltaAlpha;
		dBeta  -= dDeltaBeta;

		double dDeltaNorm = fabs(dDeltaAlpha) + fabs(dDeltaBeta);

		if (dDeltaNorm < ReferenceTolerance) {
			break;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

double GenerateMetaData(
	const Mesh & mesh,
	int nP,
	bool fBubble,
	DataMatrix3D<int> & dataGLLnodes,
	DataMatrix3D<double> & dataGLLJacobian
) {

	// Number of Faces
	int nElements = static_cast<int>(mesh.faces.size());

	// Initialize data structures
	dataGLLnodes.Initialize(nP, nP, nElements);
	dataGLLJacobian.Initialize(nP, nP, nElements);

	std::map<Node, int> mapNodes;

	// GLL Quadrature nodes
	DataVector<double> dG;
	DataVector<double> dW;
	GaussLobattoQuadrature::GetPoints(nP, 0.0, 1.0, dG, dW);

	// Accumulated Jacobian
	double dAccumulatedJacobian = 0.0;

	// Growing array of Jacobian values
	std::vector<double> vecGLLJacobian;

	// Verify face areas are available
	if (fBubble) {
		if (mesh.vecFaceArea.GetRows() != nElements) {
			_EXCEPTIONT("Face area information unavailable or incorrect");
		}
	}

	// Write metadata
	for (int k = 0; k < nElements; k++) {
		const Face & face = mesh.faces[k];
		const NodeVector & nodevec = mesh.nodes;

		if (face.edges.size() != 4) {
			_EXCEPTIONT("Input mesh must only contain quadrilateral elements");
		}

		double dFaceNumericalArea = 0.0;

		for (int j = 0; j < nP; j++) {
		for (int i = 0; i < nP; i++) {
/*
			double dXc =
				  nodevec[face[0]].x * (1.0 - dG[i]) * (1.0 - dG[j])
				+ nodevec[face[1]].x *        dG[i]  * (1.0 - dG[j])
				+ nodevec[face[2]].x *        dG[i]  *        dG[j]
				+ nodevec[face[3]].x * (1.0 - dG[i]) *        dG[j];

			double dYc =
				  nodevec[face[0]].y * (1.0 - dG[i]) * (1.0 - dG[j])
				+ nodevec[face[1]].y *        dG[i]  * (1.0 - dG[j])
				+ nodevec[face[2]].y *        dG[i]  *        dG[j]
				+ nodevec[face[3]].y * (1.0 - dG[i]) *        dG[j];

			double dZc =
				  nodevec[face[0]].z * (1.0 - dG[i]) * (1.0 - dG[j])
				+ nodevec[face[1]].z *        dG[i]  * (1.0 - dG[j])
				+ nodevec[face[2]].z *        dG[i]  *        dG[j]
				+ nodevec[face[3]].z * (1.0 - dG[i]) *        dG[j];

			double dR = sqrt(dXc * dXc + dYc * dYc + dZc * dZc);

			// Check if this Node exists in the NodeMap
			Node nodeGLL;
			nodeGLL.x = dXc / dR;
			nodeGLL.y = dYc / dR;
			nodeGLL.z = dZc / dR;

			// Calculate Jacobian

			// Pointwise basis vectors in Cartesian geometry
			Node dDx1F(
				(1.0 - dG[j]) * (nodevec[face[1]].x - nodevec[face[0]].x)
				+      dG[j]  * (nodevec[face[2]].x - nodevec[face[3]].x),
				(1.0 - dG[j]) * (nodevec[face[1]].y - nodevec[face[0]].y)
				+      dG[j]  * (nodevec[face[2]].y - nodevec[face[3]].y),
				(1.0 - dG[j]) * (nodevec[face[1]].z - nodevec[face[0]].z)
				+      dG[j]  * (nodevec[face[2]].z - nodevec[face[3]].z));

			Node dDx2F(
				(1.0 - dG[i]) * (nodevec[face[3]].x - nodevec[face[0]].x)
				+      dG[i]  * (nodevec[face[2]].x - nodevec[face[1]].x),
				(1.0 - dG[i]) * (nodevec[face[3]].y - nodevec[face[0]].y)
				+      dG[i]  * (nodevec[face[2]].y - nodevec[face[1]].y),
				(1.0 - dG[i]) * (nodevec[face[3]].z - nodevec[face[0]].z)
				+      dG[i]  * (nodevec[face[2]].z - nodevec[face[1]].z));

			// Pointwise basis vectors in spherical geometry
			double dDenomTerm = 1.0 / (dR * dR * dR);

			Node dDx1G(
				- dXc * (dYc * dDx1F.y + dZc * dDx1F.z)
					+ (dYc * dYc + dZc * dZc) * dDx1F.x,
				- dYc * (dXc * dDx1F.x + dZc * dDx1F.z)
					+ (dXc * dXc + dZc * dZc) * dDx1F.y,
				- dZc * (dXc * dDx1F.x + dYc * dDx1F.y)
					+ (dXc * dXc + dYc * dYc) * dDx1F.z);

			Node dDx2G(
				- dXc * (dYc * dDx2F.y + dZc * dDx2F.z)
					+ (dYc * dYc + dZc * dZc) * dDx2F.x,
				- dYc * (dXc * dDx2F.x + dZc * dDx2F.z)
					+ (dXc * dXc + dZc * dZc) * dDx2F.y,
				- dZc * (dXc * dDx2F.x + dYc * dDx2F.y)
					+ (dXc * dXc + dYc * dYc) * dDx2F.z);

			dDx1G.x *= dDenomTerm;
			dDx1G.y *= dDenomTerm;
			dDx1G.z *= dDenomTerm;

			dDx2G.x *= dDenomTerm;
			dDx2G.y *= dDenomTerm;
			dDx2G.z *= dDenomTerm;
*/
			// Get local map vectors
			Node nodeGLL;
			Node dDx1G;
			Node dDx2G;

			ApplyLocalMap(
				face,
				nodevec,
				dG[i],
				dG[j],
				nodeGLL,
				dDx1G,
				dDx2G);

			// Determine if this is a unique Node
			std::map<Node, int>::const_iterator iter = mapNodes.find(nodeGLL);
			if (iter == mapNodes.end()) {

				// Insert new unique node into map
				int ixNode = static_cast<int>(mapNodes.size());
				mapNodes.insert(std::pair<Node, int>(nodeGLL, ixNode));
				dataGLLnodes[j][i][k] = ixNode + 1;

			} else {
				dataGLLnodes[j][i][k] = iter->second + 1;
			}

			// Cross product gives local Jacobian
			Node nodeCross = CrossProduct(dDx1G, dDx2G);

			double dJacobian = sqrt(
				  nodeCross.x * nodeCross.x
				+ nodeCross.y * nodeCross.y
				+ nodeCross.z * nodeCross.z);

			// Element area weighted by local GLL weights
			dJacobian *= dW[i] * dW[j];

			if (dJacobian <= 0.0) {
				_EXCEPTIONT("Nonpositive Jacobian detected");
			}

			dFaceNumericalArea += dJacobian;

			dataGLLJacobian[j][i][k] = dJacobian;
		}
		}

		// Apply bubble adjustment to area
		if (fBubble && (dFaceNumericalArea != mesh.vecFaceArea[k])) {
			double dMassDifference = mesh.vecFaceArea[k] - dFaceNumericalArea;
			for (int j = 0; j < nP; j++) {
			for (int i = 0; i < nP; i++) {
				dataGLLJacobian[j][i][k] +=
					dMassDifference * dW[i] * dW[j];
 			}
			}

			dFaceNumericalArea += dMassDifference;
		}

		// Accumulate area from element
		dAccumulatedJacobian += dFaceNumericalArea;
	}

	return dAccumulatedJacobian;
}

///////////////////////////////////////////////////////////////////////////////

void GenerateUniqueJacobian(
	const DataMatrix3D<int> & dataGLLnodes,
	const DataMatrix3D<double> & dataGLLJacobian,
	DataVector<double> & dataUniqueJacobian
) {
	// Verify correct array sizes
	if ((dataGLLnodes.GetRows() != dataGLLJacobian.GetRows()) ||
		(dataGLLnodes.GetColumns() != dataGLLJacobian.GetColumns()) ||
		(dataGLLnodes.GetSubColumns() != dataGLLJacobian.GetSubColumns())
	) {
		_EXCEPTIONT("Dimension mismatch in dataGLLnodes / dataGLLJacobian");
	}

	// Find the maximum index of GLLnodes
	int iMaximumIndex = 0;

	for (int i = 0; i < dataGLLnodes.GetRows(); i++) {
	for (int j = 0; j < dataGLLnodes.GetColumns(); j++) {
	for (int k = 0; k < dataGLLnodes.GetSubColumns(); k++) {
		if (dataGLLnodes[i][j][k] > iMaximumIndex) {
			iMaximumIndex = dataGLLnodes[i][j][k];
		}
	}
	}
	}

	// Resize unique Jacobian array
	dataUniqueJacobian.Initialize(iMaximumIndex);

	for (int i = 0; i < dataGLLnodes.GetRows(); i++) {
	for (int j = 0; j < dataGLLnodes.GetColumns(); j++) {
	for (int k = 0; k < dataGLLnodes.GetSubColumns(); k++) {
		dataUniqueJacobian[dataGLLnodes[i][j][k]-1] +=
			dataGLLJacobian[i][j][k];
	}
	}
	}
}

///////////////////////////////////////////////////////////////////////////////

void GenerateDiscontinuousJacobian(
	const DataMatrix3D<double> & dataGLLJacobian,
	DataVector<double> & dataDiscontinuousJacobian
) {

	// Resize unique Jacobian array
	dataDiscontinuousJacobian.Initialize(
		  dataGLLJacobian.GetRows()
		* dataGLLJacobian.GetColumns()
		* dataGLLJacobian.GetSubColumns());

	int nP = dataGLLJacobian.GetRows();

	for (int i = 0; i < dataGLLJacobian.GetRows(); i++) {
	for (int j = 0; j < dataGLLJacobian.GetColumns(); j++) {
	for (int k = 0; k < dataGLLJacobian.GetSubColumns(); k++) {
		dataDiscontinuousJacobian[k * nP * nP + i * nP + j] =
			dataGLLJacobian[i][j][k];
	}
	}
	}

}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Global coefficient arrays (kept global to minimize reallocation)
///	</summary>
DataVector<double> g_dCoeffAlpha;
DataVector<double> g_dCoeffBeta;

void SampleGLLFiniteElement(
	bool fMonotone,
	int nP,
	double dAlpha,
	double dBeta,
	DataMatrix<double> & dCoeff
) {
	// Interpolation coefficients
	g_dCoeffAlpha.Initialize(nP);
	g_dCoeffBeta.Initialize(nP);

	// Non-monotone interpolation
	if (!fMonotone) {

		if (nP > 4) {
			// GLL Quadrature nodes on [0,1]
			DataVector<double> dG;
			DataVector<double> dW;
			GaussLobattoQuadrature::GetPoints(nP, 0.0, 1.0, dG, dW);

			// Get interpolation coefficients in each direction
			PolynomialInterp::LagrangianPolynomialCoeffs(
				nP, dG, g_dCoeffAlpha, dAlpha);

			PolynomialInterp::LagrangianPolynomialCoeffs(
				nP, dG, g_dCoeffBeta, dBeta);
		}

		// Map dAlpha and dBeta to [-1,1]
		dAlpha = 2.0 * dAlpha - 1.0;
		dBeta  = 2.0 * dBeta  - 1.0;

		// Second order monotone interpolation
		if (nP == 2) {
			g_dCoeffAlpha[0] = 0.5 * (1.0 - dAlpha);
			g_dCoeffAlpha[1] = 0.5 * (1.0 + dAlpha);

			g_dCoeffBeta[0] = 0.5 * (1.0 - dBeta);
			g_dCoeffBeta[1] = 0.5 * (1.0 + dBeta);

		// Third order interpolation
		} else if (nP == 3) {
			g_dCoeffAlpha[0] = 0.5 * (dAlpha * dAlpha - dAlpha);
			g_dCoeffAlpha[1] = 1.0 - dAlpha * dAlpha;
			g_dCoeffAlpha[2] = 0.5 * (dAlpha * dAlpha + dAlpha);

			g_dCoeffBeta[0] = 0.5 * (dBeta * dBeta - dBeta);
			g_dCoeffBeta[1] = 1.0 - dBeta * dBeta;
			g_dCoeffBeta[2] = 0.5 * (dBeta * dBeta + dBeta);

		// Fourth order interpolation
		} else if (nP == 4) {
			g_dCoeffAlpha[0] = -1.0/8.0
				* (dAlpha - 1.0) * (5.0 * dAlpha * dAlpha - 1.0);
			g_dCoeffAlpha[1] = - sqrt(5.0)/8.0
				* (sqrt(5.0) - 5.0 * dAlpha)
				* (dAlpha * dAlpha - 1.0);
			g_dCoeffAlpha[2] = - sqrt(5.0)/8.0
				* (sqrt(5.0) + 5.0 * dAlpha)
				* (dAlpha * dAlpha - 1.0);
			g_dCoeffAlpha[3] =  1.0/8.0
				* (dAlpha + 1.0) * (5.0 * dAlpha * dAlpha - 1.0);

			g_dCoeffBeta[0] = -1.0/8.0
				* (dBeta - 1.0) * (5.0 * dBeta * dBeta - 1.0);
			g_dCoeffBeta[1] = - sqrt(5.0)/8.0
				* (sqrt(5.0) - 5.0 * dBeta)
				* (dBeta * dBeta - 1.0);
			g_dCoeffBeta[2] = - sqrt(5.0)/8.0
				* (sqrt(5.0) + 5.0 * dBeta)
				* (dBeta * dBeta - 1.0);
			g_dCoeffBeta[3] =  1.0/8.0
				* (dBeta + 1.0) * (5.0 * dBeta * dBeta - 1.0);
		}

	// Monotone interpolation
	} else {

		// Map dAlpha and dBeta to [-1,1]
		dAlpha = 2.0 * dAlpha - 1.0;
		dBeta  = 2.0 * dBeta  - 1.0;

		// Second order monotone interpolation
		if (nP == 2) {
			g_dCoeffAlpha[0] = 0.5 * (1.0 - dAlpha);
			g_dCoeffAlpha[1] = 0.5 * (1.0 + dAlpha);

			g_dCoeffBeta[0] = 0.5 * (1.0 - dBeta);
			g_dCoeffBeta[1] = 0.5 * (1.0 + dBeta);

		// Third order monotone interpolation
		} else if (nP == 3) {

			if (dAlpha < 0.0) {
				g_dCoeffAlpha[0] = dAlpha * dAlpha;
				g_dCoeffAlpha[1] = 1.0 - dAlpha * dAlpha;
			} else {
				g_dCoeffAlpha[1] = 1.0 - dAlpha * dAlpha;
				g_dCoeffAlpha[2] = dAlpha * dAlpha;
			}
			if (dBeta < 0.0) {
				g_dCoeffBeta[0] = dBeta * dBeta;
				g_dCoeffBeta[1] = 1.0 - dBeta * dBeta;
			} else {
				g_dCoeffBeta[1] = 1.0 - dBeta * dBeta;
				g_dCoeffBeta[2] = dBeta * dBeta;
			}
/*
			if (dAlpha < 0.0) {
				g_dCoeffAlpha[0] = - dAlpha;
				g_dCoeffAlpha[1] = 1.0 + dAlpha;
			} else {
				g_dCoeffAlpha[1] = 1.0 - dAlpha;
				g_dCoeffAlpha[2] = dAlpha;
			}
			if (dBeta < 0.0) {
				g_dCoeffBeta[0] = - dBeta;
				g_dCoeffBeta[1] = 1.0 + dBeta;
			} else {
				g_dCoeffBeta[1] = 1.0 - dBeta;
				g_dCoeffBeta[2] = dBeta;
			}
*/
		// Fourth order monotone interpolation
		} else if (nP == 4) {

			const double dGLL1 = 1.0/sqrt(5.0);

			const double dA0 = (1.0 + sqrt(5.0)) / 16.0;
			const double dB0 = (5.0 + sqrt(5.0)) / 16.0;
			const double dC0 = (-5.0 - 5.0 * sqrt(5.0)) / 16.0;
			const double dD0 = (-25.0 - 5.0 * sqrt(5.0)) / 16.0;

			const double dA1 = 0.5;
			const double dB1 = -(3.0 / 4.0) * sqrt(5.0);
			const double dC1 = 0.0;
			const double dD1 = (5.0 / 4.0) * sqrt(5.0);

			if ((dAlpha >= -dGLL1) && (dAlpha <= dGLL1)) {
				g_dCoeffAlpha[1] =
					dA1 + dAlpha * (dB1 + dAlpha * (dC1 + dAlpha * dD1));
				g_dCoeffAlpha[2] =
					1.0 - g_dCoeffAlpha[1];
			} else if (dAlpha < -dGLL1) {
				g_dCoeffAlpha[0] =
					dA0 + dAlpha * (dB0 + dAlpha * (dC0 + dAlpha * dD0));
				g_dCoeffAlpha[1] =
					1.0 - g_dCoeffAlpha[0];
			} else {
				g_dCoeffAlpha[3] =
					dA0 - dAlpha * (dB0 - dAlpha * (dC0 - dAlpha * dD0));
				g_dCoeffAlpha[2] =
					1.0 - g_dCoeffAlpha[3];
			}

			if ((dBeta >= -dGLL1) && (dBeta <= dGLL1)) {
				g_dCoeffBeta[1] =
					dA1 + dBeta * (dB1 + dBeta * (dC1 + dBeta * dD1));
				g_dCoeffBeta[2] =
					1.0 - g_dCoeffBeta[1];
			} else if (dBeta < -dGLL1) {
				g_dCoeffBeta[0] =
					dA0 + dBeta * (dB0 + dBeta * (dC0 + dBeta * dD0));
				g_dCoeffBeta[1] =
					1.0 - g_dCoeffBeta[0];
			} else {
				g_dCoeffBeta[3] =
					dA0 - dBeta * (dB0 - dBeta * (dC0 - dBeta * dD0));
				g_dCoeffBeta[2] =
					1.0 - g_dCoeffBeta[3];
			}

		} else {
			_EXCEPTIONT("Not implemented");
		}
	}

	// Combine coefficients
	dCoeff.Initialize(nP, nP);

	for (int i = 0; i < nP; i++) {
	for (int j = 0; j < nP; j++) {
		dCoeff[j][i] = g_dCoeffAlpha[i] * g_dCoeffBeta[j];
	}
	}
/*
	// DEBUG: Override
	if (fMonotone && (nP == 4)) {

		const double dGLL1 = 1.0/sqrt(5.0);

		if ((dAlpha < -dGLL1) && (dBeta < -dGLL1)) {
			dCoeff[0][0] = 0.090903957942705;
			dCoeff[1][0] = 0.251678634003417;
			dCoeff[0][1] = 0.251678634003417;
			dCoeff[1][1] = 0.405738774050462;
			return;
		}
		if ((dAlpha < -dGLL1) && (dBeta < dGLL1)) {
			dCoeff[1][0] = 0.125362728546156;
			dCoeff[2][0] = 0.125362728546156;
			dCoeff[1][1] = 0.374637271453844;
			dCoeff[2][1] = 0.374637271453844;
			return;
		}
		if (dAlpha < -dGLL1) {
			dCoeff[2][0] = 0.251678634003417;
			dCoeff[3][0] = 0.090903957942705;
			dCoeff[3][1] = 0.251678634003417;
			dCoeff[2][1] = 0.405738774050462;
			return;
		}
		if ((dAlpha < dGLL1) && (dBeta < -dGLL1)) {
			dCoeff[0][1] = 0.125362728546156;
			dCoeff[0][2] = 0.125362728546156;
			dCoeff[1][1] = 0.374637271453844;
			dCoeff[1][2] = 0.374637271453844;
			return;
		}
		if ((dAlpha < dGLL1) && (dBeta < dGLL1)) {
			dCoeff[1][1] = 0.25;
			dCoeff[1][2] = 0.25;
			dCoeff[2][1] = 0.25;
			dCoeff[2][2] = 0.25;
			return;
		}
		if (dAlpha < dGLL1) {
			dCoeff[3][1] = 0.125362728546156;
			dCoeff[3][2] = 0.125362728546156;
			dCoeff[2][1] = 0.374637271453844;
			dCoeff[2][2] = 0.374637271453844;
			return;
		}
		if (dBeta < -dGLL1) {
			dCoeff[0][2] = 0.251678634003417;
			dCoeff[0][3] = 0.090903957942705;
			dCoeff[1][3] = 0.251678634003417;
			dCoeff[1][2] = 0.405738774050462;
			return;
		}
		if (dBeta < dGLL1) {
			dCoeff[1][3] = 0.125362728546156;
			dCoeff[2][3] = 0.125362728546156;
			dCoeff[1][2] = 0.374637271453844;
			dCoeff[2][2] = 0.374637271453844;
			return;
		}

		dCoeff[3][2] = 0.251678634003417;
		dCoeff[3][3] = 0.090903957942705;
		dCoeff[2][3] = 0.251678634003417;
		dCoeff[2][2] = 0.405738774050462;
		return;
	}
*/
}

///////////////////////////////////////////////////////////////////////////////

