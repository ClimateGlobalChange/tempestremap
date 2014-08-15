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
#include "GridElements.h"
#include "GaussLobattoQuadrature.h"

#include <map>

///////////////////////////////////////////////////////////////////////////////

double GenerateMetaData(
	const Mesh & mesh,
	int nP,
	DataMatrix3D<int> & dataGLLnodes,
	DataMatrix3D<double> & dataGLLJacobian
) {

	// Number of Faces
	int nElements = static_cast<int>(mesh.faces.size());

	// Initialize data structures
	dataGLLnodes.Initialize(nP, nP, nElements);
	dataGLLJacobian.Initialize(nP, nP, nElements);

	// Generate GLL metadata

	std::map<Node, int> mapNodes;

	// GLL Quadrature nodes
	DataVector<double> dG;
	DataVector<double> dW;
	GaussLobattoQuadrature::GetPoints(nP, 0.0, 1.0, dG, dW);

	// Accumulated Jacobian
	double dAccumulatedJacobian = 0.0;

	// Growing array of Jacobian values
	std::vector<double> vecGLLJacobian;

	// Write metadata
	for (int k = 0; k < nElements; k++) {
		const Face & face = mesh.faces[k];
		const NodeVector & nodevec = mesh.nodes;

		if (face.edges.size() != 4) {
			_EXCEPTIONT("Input mesh must only contain quadrilateral elements");
		}

		for (int j = 0; j < nP; j++) {
		for (int i = 0; i < nP; i++) {

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

			std::map<Node, int>::const_iterator iter = mapNodes.find(nodeGLL);
			if (iter == mapNodes.end()) {

				// Insert new unique node into map
				int ixNode = static_cast<int>(mapNodes.size());
				mapNodes.insert(std::pair<Node, int>(nodeGLL, ixNode));
				dataGLLnodes[j][i][k] = ixNode + 1;

			} else {
				dataGLLnodes[j][i][k] = iter->second + 1;
			}

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

			dAccumulatedJacobian += dJacobian;

			dataGLLJacobian[j][i][k] = dJacobian;
		}
		}
	}

	return dAccumulatedJacobian;
}

///////////////////////////////////////////////////////////////////////////////

