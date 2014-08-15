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

#include "LinearRemapSE0.h"
#include "OfflineMap.h"

///////////////////////////////////////////////////////////////////////////////

void LinearRemapSE0(
	const Mesh & meshInput,
	const Mesh & meshOutput,
	const Mesh & meshOverlap,
	const DataMatrix3D<int> & dataGLLNodes,
	const DataMatrix3D<double> & dataGLLJacobian,
	OfflineMap & mapRemap
) {
	// Order of the polynomial interpolant
	int nP = dataGLLNodes.GetRows();

	// Total cell Jacobian
	double dTotalJacobian;

	// Get SparseMatrix represntation of the OfflineMap
	SparseMatrix<double> & smatMap = mapRemap.GetSparseMatrix();

	// Loop through all elements in the overlap mesh
	int ixLastFirstMeshFace = (-1);
	int ixCurrentFirstMeshFace;
	int ixCurrentSecondMeshFace;

	for (int i = 0; i < meshOverlap.faces.size(); i++) {
		ixCurrentFirstMeshFace = meshOverlap.vecFirstFaceIx[i];
		ixCurrentSecondMeshFace = meshOverlap.vecSecondFaceIx[i];

		double dSecondFaceArea =
			meshOutput.vecFaceArea[ixCurrentSecondMeshFace];

		// Recalculate total element Jacobian for this first Mesh element
		if (ixLastFirstMeshFace != ixCurrentFirstMeshFace) {

			// Calculate total element Jacobian
			dTotalJacobian = 0.0;
			for (int p = 0; p < nP; p++) {
			for (int q = 0; q < nP; q++) {
				dTotalJacobian += dataGLLJacobian[p][q][ixCurrentFirstMeshFace];
			}
			}
		}

		// Determine remap coefficients
		for (int p = 0; p < nP; p++) {
		for (int q = 0; q < nP; q++) {
			int ixFirstGlobal = dataGLLNodes[p][q][ixCurrentFirstMeshFace] - 1;

			smatMap(ixCurrentSecondMeshFace, ixFirstGlobal) +=
				dataGLLJacobian[p][q][ixCurrentFirstMeshFace]
				/ dTotalJacobian
				* meshOverlap.vecFaceArea[i]
				/ dSecondFaceArea;
		}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

