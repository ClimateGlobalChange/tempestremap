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
#include "FiniteElementTools.h"
#include "GaussLobattoQuadrature.h"

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

			ixLastFirstMeshFace = ixCurrentFirstMeshFace;
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

void LinearRemapSE4(
	const Mesh & meshInput,
	const Mesh & meshOutput,
	const Mesh & meshOverlap,
	const DataMatrix3D<int> & dataGLLNodes,
	const DataMatrix3D<double> & dataGLLJacobian,
	bool fMonotone,
	OfflineMap & mapRemap
) {
	// Triangular quadrature to use for integration
	// Dunavant, D.A. "High Degree Efficient Symmetrical Gaussian Quadrature
	// Rules for the Triangle."  J. Numer. Meth. Eng., 21, pp 1129-1148.
/*
	const int TriQuadraturePoints = 6;

	const double TriQuadratureG[6][3] = {
		{0.108103018168070, 0.445948490915965, 0.445948490915965},
		{0.445948490915965, 0.108103018168070, 0.445948490915965},
		{0.445948490915965, 0.445948490915965, 0.108103018168070},
		{0.816847572980458, 0.091576213509771, 0.091576213509771},
		{0.091576213509771, 0.816847572980458, 0.091576213509771},
		{0.091576213509771, 0.091576213509771, 0.816847572980458}};

	const double TriQuadratureW[6] =
		{0.223381589678011, 0.223381589678011, 0.223381589678011,
		 0.109951743655322, 0.109951743655322, 0.109951743655322};
*/
	const double TriQuadraturePoints = 1;

	const double TriQuadratureG[1][3] = {
		{0.333333333333333, 0.333333333333333, 0.333333333333333}};

	const double TriQuadratureW[1] =
		{1.000000000000000};

	// Order of the polynomial interpolant
	int nP = dataGLLNodes.GetRows();

	// GLL Quadrature nodes on quadrilateral elements
	DataVector<double> dG;
	DataVector<double> dW;
	GaussLobattoQuadrature::GetPoints(nP, 0.0, 1.0, dG, dW);

	// Total cell Jacobian
	double dTotalJacobian;

	// Sample coefficients
	DataMatrix<double> dSampleCoeff;

	// Get SparseMatrix represntation of the OfflineMap
	SparseMatrix<double> & smatMap = mapRemap.GetSparseMatrix();

	// NodeVector from meshOverlap
	const NodeVector & nodesOverlap = meshOverlap.nodes;
	const NodeVector & nodesFirst   = meshInput.nodes;

	// Loop through all elements in the overlap mesh
	int ixLastFirstMeshFace = (-1);
	int ixCurrentFirstMeshFace;
	int ixCurrentSecondMeshFace;

	// Triangle area and quadrature nodes within triangles
	DataVector<double> dTriArea;
	DataMatrix3D<double> dTriNodes;

	// Current sub-triangle
	int iCurrentTri = 0;

	// Construct set of quadrature nodes on each Overlap Face
	std::set<int> setFirstFacesVisited;

	for (int i = 0; i < meshOverlap.faces.size(); i++) {
		ixCurrentFirstMeshFace = meshOverlap.vecFirstFaceIx[i];
		ixCurrentSecondMeshFace = meshOverlap.vecSecondFaceIx[i];

		const Face & faceOverlap = meshOverlap.faces[i];
		const Face & faceFirst   = meshInput.faces[ixCurrentFirstMeshFace];

		if (faceFirst.edges.size() != 4) {
			_EXCEPTIONT("Only quadrilateral elements allowed for SE remapping");
		}

		double dFirstFaceArea =
			meshInput.vecFaceArea[ixCurrentFirstMeshFace];
		double dSecondFaceArea =
			meshOutput.vecFaceArea[ixCurrentSecondMeshFace];

		// Calculate quantities relevant to this FirstMesh Face
		if (ixLastFirstMeshFace != ixCurrentFirstMeshFace) {

			if (setFirstFacesVisited.find(ixCurrentFirstMeshFace)
				!= setFirstFacesVisited.end()
			) {
				_EXCEPTIONT("Overmap Mesh must be sorted on First Face index");
			}
			setFirstFacesVisited.insert(ixCurrentFirstMeshFace);

			// Determine how many sub-Faces are here
			int nSubFaces = 0;
			int nTotalSubTriangles = 0;
			for (;;) {
				if (i + nSubFaces == meshOverlap.vecFirstFaceIx.size()) {
					break;
				}

				int iSubFirstMeshFace =
					meshOverlap.vecFirstFaceIx[i + nSubFaces];

				if (iSubFirstMeshFace != ixCurrentFirstMeshFace) {
					break;
				}

				nTotalSubTriangles +=
					meshOverlap.faces[i + nSubFaces].edges.size() - 2;

				nSubFaces++;
			}

			// Allocate quadrature point array
			dTriArea.Initialize(nTotalSubTriangles);
			dTriNodes.Initialize(nTotalSubTriangles, TriQuadraturePoints, 2);

			// Calculate triangle areas and quadrature points
			iCurrentTri = 0;

			for (int j = 0; j < nSubFaces; j++) {
				const Face & faceSub = meshOverlap.faces[i + j];

				int nSubTriangles = faceSub.edges.size() - 2;

				for (int k = 0; k < nSubTriangles; k++) {

					// Cornerpoints of triangle
					const Node & node0 = nodesOverlap[faceSub[0]];
					const Node & node1 = nodesOverlap[faceSub[k+1]];
					const Node & node2 = nodesOverlap[faceSub[k+2]];

					// Calculate the area of the modified Face
					Face faceTri(3);
					faceTri.SetNode(0, faceSub[0]);
					faceTri.SetNode(1, faceSub[k+1]);
					faceTri.SetNode(2, faceSub[k+2]);

					dTriArea[iCurrentTri] =
						CalculateFaceArea(faceTri, nodesOverlap);

					// Coordinates of quadrature Node
					for (int l = 0; l < TriQuadraturePoints; l++) {
						Node nodeQuadrature;
						nodeQuadrature.x =
							  TriQuadratureG[l][0] * node0.x
							+ TriQuadratureG[l][1] * node1.x
							+ TriQuadratureG[l][2] * node2.x;

						nodeQuadrature.y =
							  TriQuadratureG[l][0] * node0.y
							+ TriQuadratureG[l][1] * node1.y
							+ TriQuadratureG[l][2] * node2.y;

						nodeQuadrature.z =
							  TriQuadratureG[l][0] * node0.z
							+ TriQuadratureG[l][1] * node1.z
							+ TriQuadratureG[l][2] * node2.z;

						double dMag = sqrt(
							  nodeQuadrature.x * nodeQuadrature.x
							+ nodeQuadrature.y * nodeQuadrature.y
							+ nodeQuadrature.z * nodeQuadrature.z);

						nodeQuadrature.x /= dMag;
						nodeQuadrature.y /= dMag;
						nodeQuadrature.z /= dMag;

						// Find components of quadrature point in basis
						// of the first Face
						ApplyInverseMap(
							faceFirst,
							nodesFirst,
							nodeQuadrature,
							dTriNodes[iCurrentTri][l][0],
							dTriNodes[iCurrentTri][l][1]);

						// Check inverse map value
						if ((dTriNodes[iCurrentTri][l][0] < 0.0) ||
							(dTriNodes[iCurrentTri][l][0] > 1.0) ||
							(dTriNodes[iCurrentTri][l][1] < 0.0) ||
							(dTriNodes[iCurrentTri][l][1] > 1.0)
						) {
							_EXCEPTIONT("Inverse Map out of range");
						}
					}

					// Increment current triangle
					iCurrentTri++;
				}
			}

			// Calculate total element Jacobian
			dTotalJacobian = 0.0;
			for (int j = 0; j < nTotalSubTriangles; j++) {
			for (int k = 0; k < TriQuadraturePoints; k++) {

				// Interpolate Jacobian to Node
				SampleGLLFiniteElement(
					fMonotone, nP,
					dTriNodes[j][k][0],
					dTriNodes[j][k][1],
					dSampleCoeff);

				for (int p = 0; p < nP; p++) {
				for (int q = 0; q < nP; q++) {
					dTotalJacobian +=
						TriQuadratureW[k]
						* dTriArea[j]
						* dSampleCoeff[p][q]
						* dataGLLJacobian[p][q][ixCurrentFirstMeshFace]
						/ (dW[p] * dW[q]);
				}
				}
			}
			}

			ixLastFirstMeshFace = ixCurrentFirstMeshFace;

			// Reset current triangle index
			iCurrentTri = 0;
		}

		// Loop through all sub-triangles for this OverlapMesh Face
		for (int j = 0; j < faceOverlap.edges.size()-2; j++) {

			// Get contribution from each quadrature point
			for (int k = 0; k < TriQuadraturePoints; k++) {

				// Sample the finite element at this point
				SampleGLLFiniteElement(
					fMonotone, nP,
					dTriNodes[iCurrentTri][k][0],
					dTriNodes[iCurrentTri][k][1],
					dSampleCoeff);

				// Add sample coefficients to the map
				for (int p = 0; p < nP; p++) {
				for (int q = 0; q < nP; q++) {
					int ixNodeFirst =
						dataGLLNodes[p][q][ixCurrentFirstMeshFace] - 1;

					double dJacobian =
						dataGLLJacobian[p][q][ixCurrentFirstMeshFace];

					smatMap(ixCurrentSecondMeshFace, ixNodeFirst) +=
						TriQuadratureW[k]
						* dTriArea[iCurrentTri]
						* dSampleCoeff[p][q]
						//* dJacobian / (dW[p] * dW[q])
						// / dTotalJacobian
						/ dSecondFaceArea;
				}
				}
			}

			// Increment current triangle
			iCurrentTri++;
		}

		// Compute triangular quadrature rule over each triangle
/*
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
*/
	}
}

///////////////////////////////////////////////////////////////////////////////

