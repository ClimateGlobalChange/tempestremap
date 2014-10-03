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

#include "LinearRemapFV.h"
#include "GridElements.h"
#include "OfflineMap.h"
#include "FiniteElementTools.h"
#include "GaussLobattoQuadrature.h"
#include "TriangularQuadrature.h"

#include "Announce.h"
#include "MathHelper.h"

#include <map>

///////////////////////////////////////////////////////////////////////////////

extern "C" {
	///	General matrix solver from CLAPACK
	int dgesv_(
		int * n,
		int * nrhs,
		double * a,
		int * lda,
		int * ipiv,
		double * b,
		int * ldb,
		int * info);

	///	Compute the QR decomposition from CLAPACK
	int dgeqrf_(
		int * m,
		int * n,
		double * a,
		int * lda,
		double * tau,
		double * work,
		int * lwork,
		int * info);

	// Compute the matrix Q from DGEQRF from CLAPACK
	int dorgqr_(
		int * m,
		int * n,
		int * k,
		double * a,
		int * lda,
		double * tau,
		double * work,
		int * lwork,
		int * info);

	/// Compute an LU Factorization from CLAPACK
	int dgetrf_(
		int * m,
		int * n,
		double * a,
		int * lda,
		int * ipiv,
		int * info);

	///	Compute the matrix inverse from the LU factorization from CLAPACK
	int dgetri_(
		int * n,
		double * a,
		int * lda,
		int * ipiv,
		double * work,
		int * lwork,
		int * info);

	/// Matrix vector multiply
	int dgemv_(
		char * trans,
		int * m,
		int * n,
		double * alpha,
		double * a,
		int * lda,
		double * x,
		int * incx,
		double * beta,
		double * y,
		int * incy);

	/// Symmetric positive definite matrix solver from CLAPACK
	int dposv_(
		char * uplo,
		int * n,
		int * nrhs,
		double * a,
		int * lda,
		double * b,
		int * ldb,
		int * info);

	/// Symmetric matrix solver from CLAPACK
	int dsysv_(
		char * uplo,
		int * n,
		int * nrhs,
		double * a,
		int * lda,
		int * ipiv,
		double * b,
		int * ldb,
		double * work,
		int * lwork,
		int * info);

	/// Constrained least squares solver
	int dgglse_(
		int * m,
		int * n,
		int * p,
		double * a,
		int * lda,
		double * b,
		int * ldb,
		double * c,
		double * d,
		double * x,
		double * work,
		int * lwork,
		int * info);

};

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

void GetAdjacentFaceVector(
	const Mesh & mesh,
	int iFaceInitial,
	int nRequiredFaceSetSize,
	AdjacentFaceVector & vecFaces
) {
	// Ensure the ReverseNodeArray has been constructed
	if (mesh.revnodearray.size() == 0) {
		_EXCEPTIONT("ReverseNodeArray is required");
	}

	// Insert the initial Face
	vecFaces.push_back(FaceDistancePair(iFaceInitial,1));

	// Initial Face
	const Face & faceInitial = mesh.faces[iFaceInitial];

	// Set of nodes at the "perimeter" of faceInitial
	std::set<int> setPerimeterNodes;

	for (int j = 0; j < faceInitial.edges.size(); j++) {
		setPerimeterNodes.insert(faceInitial[j]);
	}

	// Set of Faces already added to vecFaces
	std::set<int> setFaces;
	setFaces.insert(iFaceInitial);

	// Generate the set of faces
	int iDistance = 1;
	for (;;) {
		if (vecFaces.size() >= nRequiredFaceSetSize) {
			break;
		}

		// Increment distance metric
		iDistance++;

		// Store old perimeter nodes
		std::set<int> setPerimeterNodesOld = setPerimeterNodes;

		setPerimeterNodes.clear();

		// Loop through all perimeter nodes and add corresponding Faces
		std::set<int>::const_iterator iterNode = setPerimeterNodesOld.begin();

		for (; iterNode != setPerimeterNodesOld.end(); iterNode++) {
			const std::set<int> & setAdjFaces =
				mesh.revnodearray[*iterNode];

			std::set<int>::const_iterator iterFace = setAdjFaces.begin();
			for (; iterFace != setAdjFaces.end(); iterFace++) {

				// Verify this Face has not already been added
				if (setFaces.find(*iterFace) != setFaces.end()) {
					continue;
				}

				// Add this Face to vecFaces and setFaces
				vecFaces.push_back(FaceDistancePair(*iterFace, iDistance));
				setFaces.insert(*iterFace);

				// Add new nodes to perimeter node set
				const Face & faceAdj = mesh.faces[*iterFace];

				for (int j = 0; j < faceAdj.edges.size(); j++) {
					if (setPerimeterNodesOld.find(faceAdj[j]) ==
						setPerimeterNodesOld.end()
					) {
						setPerimeterNodes.insert(faceAdj[j]);
					}
				}
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void BuildFitArray(
	const Mesh & mesh,
	int ixFirst,
	const AdjacentFaceVector & vecAdjFaces,
	int nOrder,
	DataVector<double> & dConstraint,
	DataMatrix<double> & dFitArray,
	DataVector<double> & dFitWeights,
	DataMatrix<double> & dFitArrayPlus
) {

	// Reference to active Face
	const Face & faceFirst = mesh.faces[ixFirst];

	// Number of coefficients
	int nCoefficients = nOrder * (nOrder + 1) / 2;

	// Number of adjacent Faces
	int nAdjFaces = vecAdjFaces.size();

	// Initialize arrays,
	dFitArray.Initialize(nCoefficients, nAdjFaces);
	dFitWeights.Initialize(nAdjFaces);
	dFitArrayPlus.Initialize(nAdjFaces, nCoefficients);

	// Triangular quadrature rule
	TriangularQuadratureRule triquadrule(4);

	const DataMatrix<double> & dG = triquadrule.GetG();
	const DataVector<double> & dW = triquadrule.GetW();

	// Coordinate axes
	const Node & nodeRef = mesh.nodes[faceFirst[0]];

	Node nodeA1 = mesh.nodes[faceFirst[1]] - nodeRef;
	Node nodeA2 = mesh.nodes[faceFirst[2]] - nodeRef;

	Node nodeC = CrossProduct(nodeA1, nodeA2);

	// Fit matrix
	DataMatrix<double> dFit;
	dFit.Initialize(3,3);

	dFit[0][0] = nodeA1.x; dFit[0][1] = nodeA1.y; dFit[0][2] = nodeA1.z;
	dFit[1][0] = nodeA2.x; dFit[1][1] = nodeA2.y; dFit[1][2] = nodeA2.z;
	dFit[2][0] = nodeC.x;  dFit[2][1] = nodeC.y;  dFit[2][2] = nodeC.z;

	// Loop through all adjacent Faces
	for (int iAdjFace = 0; iAdjFace < vecAdjFaces.size(); iAdjFace++) {

		const FaceDistancePair & fdp = vecAdjFaces[iAdjFace];

		const Face & faceAdj = mesh.faces[fdp.first];

		// Adjacent face area
		double dAdjFaceArea = mesh.vecFaceArea[fdp.first];

		// Loop through all sub-triangles
		for (int j = 0; j < faceAdj.edges.size()-2; j++) {

			const Node & node0 = mesh.nodes[faceAdj[0]];
			const Node & node1 = mesh.nodes[faceAdj[j+1]];
			const Node & node2 = mesh.nodes[faceAdj[j+2]];

			// Calculate sub-triangular area
			Face faceTri(3);
			faceTri.SetNode(0, faceAdj[0]);
			faceTri.SetNode(1, faceAdj[j+1]);
			faceTri.SetNode(2, faceAdj[j+2]);

			double dTriArea = CalculateFaceArea(faceTri, mesh.nodes);

			// Loop through all triangular quadrature nodes
			for (int k = 0; k < triquadrule.GetPoints(); k++) {

				double * dGL = dG[k];

				// Get the nodal location of this point
				double dX[3];

				dX[0] = dGL[0] * node0.x + dGL[1] * node1.x + dGL[2] * node2.x;
				dX[1] = dGL[0] * node0.y + dGL[1] * node1.y + dGL[2] * node2.y;
				dX[2] = dGL[0] * node0.z + dGL[1] * node1.z + dGL[2] * node2.z;

				double dMag =
					sqrt(dX[0] * dX[0] + dX[1] * dX[1] + dX[2] * dX[2]);

				dX[0] /= dMag;
				dX[1] /= dMag;
				dX[2] /= dMag;

				dX[0] -= nodeRef.x;
				dX[1] -= nodeRef.y;
				dX[2] -= nodeRef.z;

				// Find the coefficients for this point
				int n = 3;
				int nrhs = 1;
				int lda = 3;
				int ipiv[9];
				int ldb = 3;
				int info;

#pragma message "Pre-triangularize this matrix for efficiency"
				DataMatrix<double> dFitTemp;
				dFitTemp = dFit;
				dgesv_(
					&n, &nrhs, &(dFitTemp[0][0]), &lda, ipiv, dX, &ldb, &info);

				if (info != 0) {
					_EXCEPTIONT("Solve failure in dgesv");
				}

				// Loop through all coefficients
				int ixp = 0;
				for (int p = 0; p < nOrder; p++) {
				for (int q = 0; q < nOrder - p; q++) {

					dFitArray[ixp][iAdjFace] +=
						  IPow(dX[0], p)
						* IPow(dX[1], q)
						* dW[k]
						* dTriArea
						/ dAdjFaceArea;

					ixp++;
				}
				}
			}
		}
/*
		// Reweight the fit array
		dFitWeights[iAdjFace] =
			pow(static_cast<double>(fdp.second),
			    - static_cast<double>(nOrder+1));

		for (int j = 0; j < dFitArray.GetRows(); j++) {
			dFitArray[j][iAdjFace] *= dFitWeights[iAdjFace];
		}
*/
	}

	// First order
	if (nOrder == 1) {

		// Invert the fit array using Moore-Penrose pseudoinverse
		DataMatrix<double> dFit2;
		dFit2.Initialize(nCoefficients, nCoefficients);

		for (int j = 0; j < nCoefficients; j++) {
		for (int k = 0; k < nCoefficients; k++) {
		for (int l = 0; l < nAdjFaces; l++) {
			dFit2[j][k] += dFitArray[j][l] * dFitArray[k][l];
		}
		}
		}

		// Calculate pseudoinverse of FitArray
		int m = nCoefficients;
		int n = nCoefficients;
		int lda = nCoefficients;
		int info;

		DataVector<int> iPIV;
		iPIV.Initialize(nCoefficients);

		DataVector<double> dWork;
		dWork.Initialize(nCoefficients);

		int lWork = nCoefficients;

		dgetrf_(&m, &n, &(dFit2[0][0]), &lda, &(iPIV[0]), &info);

		dgetri_(&n, &(dFit2[0][0]), &lda, &(iPIV[0]), &(dWork[0]), &lWork, &info);

		// Calculate pseudoinverse
		for (int j = 0; j < nAdjFaces; j++) {
		for (int k = 0; k < nCoefficients; k++) {
		for (int l = 0; l < nCoefficients; l++) {
			dFitArrayPlus[j][k] += dFitArray[l][j] * dFit2[l][k];
		}
		}
		}
/*
		for (int j = 0; j < nAdjFaces; j++) {
		for (int k = 0; k < nCoefficients; k++) {
			dFitArrayPlus[j][k] *= dFitWeights[j];
		}
		}
*/
		return;
	}
/*
	for (int i = 0; i < dConstraint.GetRows(); i++) {
		printf("%1.15e\n", dConstraint[i]);
	}

	for (int i = 0; i < dFitArray.GetRows(); i++) {
		for (int j = 0; j < dFitArray.GetColumns(); j++) {
			printf("%1.15e  ", dFitArray[i][j]);
		}
		printf(";\n");
	}
*/
	// Compute QR factorization of the constraint
	DataMatrix<double> dQ;
	dQ.Initialize(nCoefficients, nCoefficients);

	double dR;

	{
		int m = nCoefficients;
		int n = 1;
		int lda = m;

		double tau;

		DataVector<double> dWork;
		int lwork = nCoefficients;
		dWork.Initialize(nCoefficients);

		int info;

		memcpy(&(dQ[0][0]), &(dConstraint[0]), nCoefficients * sizeof(double));

		dgeqrf_(&m, &n, &(dQ[0][0]), &lda, &tau, &(dWork[0]), &lwork, &info);
		if (info != 0) {
			_EXCEPTION1("Error in dgeqrf: %i", info);
		}

		dR = dQ[0][0];

		int k = 1;
		n = nCoefficients;
		dorgqr_(&m, &n, &k, &(dQ[0][0]), &lda, &tau, &(dWork[0]), &lwork, &info);
		if (info != 0) {
			_EXCEPTION1("Error in dorgqr: %i", info);
		}
	}
/*
	printf("R: %1.15e\n", dR);
	printf("Q:\n");
	for (int i = 0; i < nCoefficients; i++) {
		for (int j = 0; j < nCoefficients; j++) {
			printf("%1.15e  ", dQ[j][i]);
		}
		printf("\n");
	}
*/
	// Calculate G = F * Q 
	DataMatrix<double> dGG;
	dGG.Initialize(nCoefficients, nAdjFaces);

	for (int i = 0; i < nCoefficients; i++) {
	for (int j = 0; j < nAdjFaces; j++) {
		for (int k = 0; k < nCoefficients; k++) {
			dGG[i][j] += dFitArray[k][j] * dQ[i][k];
		}
	}
	}
/*
	printf("G:\n");
	for (int i = 0; i < nAdjFaces; i++) {
		for (int j = 0; j < nCoefficients; j++) {
			printf("%1.15e  ", dGG[j][i]);
		}
		printf("\n");
	}
*/
	// Calculate Moore-Penrose pseudoinverse of G(:,2:p)
	DataMatrix<double> dGxPlus;
	dGxPlus.Initialize(nAdjFaces, nCoefficients-1);

	{
		// Gx2 = G(:, 2:p)^T * G(:,2:p)
		DataMatrix<double> dGx2;
		dGx2.Initialize(nCoefficients-1, nCoefficients-1);

		for (int i = 0; i < nCoefficients-1; i++) {
		for (int j = 0; j < nCoefficients-1; j++) {
			for (int k = 0; k < nAdjFaces; k++) {
				dGx2[i][j] += dGG[i+1][k] * dGG[j+1][k];
			}
		}
		}

		int m = nCoefficients-1;
		int n = nCoefficients-1;
		int lda = nCoefficients-1;
		int info;

		DataVector<int> iPIV;
		iPIV.Initialize(nCoefficients-1);

		DataVector<double> dWork;
		dWork.Initialize(nCoefficients-1);

		int lWork = nCoefficients-1;

		dgetrf_(&m, &n, &(dGx2[0][0]), &lda, &(iPIV[0]), &info);

		dgetri_(&n, &(dGx2[0][0]), &lda, &(iPIV[0]), &(dWork[0]), &lWork, &info);

		// Calculate pseudoinverse
		for (int i = 0; i < nAdjFaces; i++) {
		for (int j = 0; j < nCoefficients-1; j++) {
			for (int k = 0; k < nCoefficients-1; k++) {
				dGxPlus[i][j] += dGG[k+1][i] * dGx2[k][j];
			}
		}
		}
	}
/*
	printf("Gp:\n");
	for (int i = 0; i < nCoefficients-1; i++) {
		for (int j = 0; j < nAdjFaces; j++) {
			printf("%1.15e  ", dGxPlus[j][i]);
		}
		printf("\n");
	}
*/
	// Z = G+ * (I - G(:,1) * R^{-1} * e0T)
	DataMatrix<double> dZ;
	dZ.Initialize(nAdjFaces, nCoefficients-1);

	{
		DataMatrix<double> dSubZ;
		dSubZ.Initialize(nAdjFaces, nAdjFaces);

		for (int i = 0; i < nAdjFaces; i++) {
			dSubZ[i][i] = 1.0;
		}

		for (int i = 0; i < nAdjFaces; i++) {
			dSubZ[0][i] -= dGG[0][i] / dR;
		}

		for (int i = 0; i < nAdjFaces; i++) {
		for (int j = 0; j < nCoefficients-1; j++) {
			for (int k = 0; k < nAdjFaces; k++) {
				dZ[i][j] += dGxPlus[k][j] * dSubZ[i][k];
			}
		}
		}
	}
/*
	printf("Z:\n");
	for (int i = 0; i < nCoefficients-1; i++) {
		for (int j = 0; j < nAdjFaces; j++) {
			printf("%1.15e  ", dZ[j][i]);
		}
		printf("\n");
	}
*/
	// Fhat = Q(:,1) * R^{-1} * e0T + Q(:,2:p) * Z
	for (int i = 0; i < nCoefficients; i++) {
		dFitArrayPlus[0][i] += dQ[0][i] / dR;
	}
	for (int i = 0; i < nAdjFaces; i++) {
	for (int j = 0; j < nCoefficients; j++) {
		for (int k = 0; k < nCoefficients-1; k++) {
			dFitArrayPlus[i][j] += dQ[k+1][j] * dZ[i][k];
		}
	}
	}
/*
	printf("Fp:\n");
	for (int i = 0; i < nCoefficients; i++) {
		for (int j = 0; j < nAdjFaces; j++) {
			printf("%1.15e  ", dFitArrayPlus[j][i]);
		}
		printf("\n");
	}

	_EXCEPTION();
*/
/*
	// Invert the fit array using Moore-Penrose pseudoinverse
	DataMatrix<double> dFit2;
	dFit2.Initialize(nCoefficients, nCoefficients);

	for (int j = 0; j < nCoefficients; j++) {
	for (int k = 0; k < nCoefficients; k++) {
	for (int l = 0; l < nAdjFaces; l++) {
		dFit2[j][k] += dFitArray[j][l] * dFitArray[k][l];
	}
	}
	}

	// Calculate pseudoinverse of FitArray
	int m = nCoefficients;
	int n = nCoefficients;
	int lda = nCoefficients;
	int info;

	DataVector<int> iPIV;
	iPIV.Initialize(nCoefficients);

	DataVector<double> dWork;
	dWork.Initialize(nCoefficients);

	int lWork = nCoefficients;

	dgetrf_(&m, &n, &(dFit2[0][0]), &lda, &(iPIV[0]), &info);

	dgetri_(&n, &(dFit2[0][0]), &lda, &(iPIV[0]), &(dWork[0]), &lWork, &info);

	// Calculate pseudoinverse
	for (int j = 0; j < nAdjFaces; j++) {
	for (int k = 0; k < nCoefficients; k++) {
	for (int l = 0; l < nCoefficients; l++) {
		dFitArrayPlus[j][k] += dFitArray[l][j] * dFit2[l][k];
	}
	}
	}

	for (int j = 0; j < nAdjFaces; j++) {
	for (int k = 0; k < nCoefficients; k++) {
		dFitArrayPlus[j][k] *= dFitWeights[j];
	}
	}
*/
}

///////////////////////////////////////////////////////////////////////////////

void LinearRemapFVtoFV(
	const Mesh & meshInput,
	const Mesh & meshOutput,
	const Mesh & meshOverlap,
	int nOrder,
	OfflineMap & mapRemap
) {
	// Verify ReverseNodeArray has been calculated
	if (meshInput.revnodearray.size() == 0) {
		_EXCEPTIONT("ReverseNodeArray has not been calculated for meshInput");
	}

	// Triangular quadrature rule
	TriangularQuadratureRule triquadrule(4);

	const DataMatrix<double> & dG = triquadrule.GetG();
	const DataVector<double> & dW = triquadrule.GetW();

	// Get SparseMatrix represntation of the OfflineMap
	SparseMatrix<double> & smatMap = mapRemap.GetSparseMatrix();

	// Number of elements needed
	int nCoefficients = nOrder * (nOrder + 1) / 2;

	int nRequiredFaceSetSize = nCoefficients;

	// Current overlap face
	int ixOverlap = 0;

	// Loop through all faces on meshInput
	for (int ixFirst = 0; ixFirst < meshInput.faces.size(); ixFirst++) {

		// Output every 100 elements
		if (ixFirst % 100 == 0) {
			Announce("Element %i", ixFirst);
		}

		// This Face
		const Face & faceFirst = meshInput.faces[ixFirst];

		// Coordinate axes
		const Node & nodeRef = meshInput.nodes[faceFirst[0]];

		Node nodeA1 = meshInput.nodes[faceFirst[1]] - nodeRef;
		Node nodeA2 = meshInput.nodes[faceFirst[2]] - nodeRef;

		Node nodeC = CrossProduct(nodeA1, nodeA2);

		// Fit matrix
		DataMatrix<double> dFit;
		dFit.Initialize(3,3);

		dFit[0][0] = nodeA1.x; dFit[0][1] = nodeA1.y; dFit[0][2] = nodeA1.z;
		dFit[1][0] = nodeA2.x; dFit[1][1] = nodeA2.y; dFit[1][2] = nodeA2.z;
		dFit[2][0] = nodeC.x;  dFit[2][1] = nodeC.y;  dFit[2][2] = nodeC.z;

		// Number of overlapping Faces and triangles
		int nOverlapFaces = 0;
		int nTotalOverlapTriangles = 0;

		// Determine how many overlap Faces and triangles are present
		int ixOverlapTemp = ixOverlap;
		for (; ixOverlapTemp < meshOverlap.faces.size(); ixOverlapTemp++) {

			const Face & faceOverlap = meshOverlap.faces[ixOverlapTemp];

			if (meshOverlap.vecFirstFaceIx[ixOverlapTemp] != ixFirst) {
				break;
			}

			nOverlapFaces++;
			nTotalOverlapTriangles += faceOverlap.edges.size() - 2;
		}

		// Build integration array
		DataMatrix<double> dIntArray;
		dIntArray.Initialize(nCoefficients, nOverlapFaces);

		// Loop through all overlap Faces
		for (int i = 0; i < nOverlapFaces; i++) {
			const Face & faceOverlap = meshOverlap.faces[ixOverlap + i];

			const NodeVector & nodesOverlap = meshOverlap.nodes;

			int nOverlapTriangles = faceOverlap.edges.size() - 2;

			// Loop over all sub-triangles of this Overlap Face
			for (int j = 0; j < nOverlapTriangles; j++) {

				// Cornerpoints of triangle
				const Node & node0 = nodesOverlap[faceOverlap[0]];
				const Node & node1 = nodesOverlap[faceOverlap[j+1]];
				const Node & node2 = nodesOverlap[faceOverlap[j+2]];

				// Calculate the area of the modified Face
				Face faceTri(3);
				faceTri.SetNode(0, faceOverlap[0]);
				faceTri.SetNode(1, faceOverlap[j+1]);
				faceTri.SetNode(2, faceOverlap[j+2]);

				double dTriArea =
					CalculateFaceArea(faceTri, nodesOverlap);

				for (int k = 0; k < triquadrule.GetPoints(); k++) {
					double * dGL = dG[k];

					// Get the nodal location of this point
					double dX[3];

					dX[0] = dGL[0] * node0.x + dGL[1] * node1.x + dGL[2] * node2.x;
					dX[1] = dGL[0] * node0.y + dGL[1] * node1.y + dGL[2] * node2.y;
					dX[2] = dGL[0] * node0.z + dGL[1] * node1.z + dGL[2] * node2.z;

					double dMag =
						sqrt(dX[0] * dX[0] + dX[1] * dX[1] + dX[2] * dX[2]);

					dX[0] /= dMag;
					dX[1] /= dMag;
					dX[2] /= dMag;

					dX[0] -= nodeRef.x;
					dX[1] -= nodeRef.y;
					dX[2] -= nodeRef.z;

					// Find the coefficients for this point
					int n = 3;
					int nrhs = 1;
					int lda = 3;
					int ipiv[3];
					int ldb = 3;
					int info;

					DataMatrix<double> dFitTemp;
					dFitTemp = dFit;
					dgesv_(
						&n, &nrhs, &(dFitTemp[0][0]), &lda, ipiv, dX, &ldb, &info);

					// Sample this point
					int ixp = 0;
					for (int p = 0; p < nOrder; p++) {
					for (int q = 0; q < nOrder - p; q++) {
						dIntArray[ixp][i] +=
							  IPow(dX[0], p)
							* IPow(dX[1], q)
							* dW[k]
							* dTriArea;

						ixp++;
					}
					}
				}
			}
		}

		// Determine the conservative constraint equation
		DataVector<double> dConstraint;
		dConstraint.Initialize(nCoefficients);

		double dFirstArea = meshInput.vecFaceArea[ixFirst];

		for (int p = 0; p < nCoefficients; p++) {
		for (int j = 0; j < nOverlapFaces; j++) {
			dConstraint[p] += dIntArray[p][j] / dFirstArea;
		}
		}

		// Set of Faces to use in building the reconstruction and associated
		// distance metric.
		AdjacentFaceVector vecAdjFaces;

		GetAdjacentFaceVector(
			meshInput,
			ixFirst,
			nRequiredFaceSetSize,
			vecAdjFaces);

		// Number of adjacent Faces
		int nAdjFaces = vecAdjFaces.size();

		// Least squares arrays
		DataMatrix<double> dFitArray;
		DataVector<double> dFitWeights;
		DataMatrix<double> dFitArrayPlus;

		BuildFitArray(
			meshInput,
			ixFirst,
			vecAdjFaces,
			nOrder,
			dConstraint,
			dFitArray,
			dFitWeights,
			dFitArrayPlus
		);
/*
		printf("\n");
		for (int i = 0; i < nAdjFaces; i++) {
		for (int p = 0; p < nCoefficients; p++) {
			printf("%1.3e  ", dFitArray[p][i]);
		}
			printf("\n");
		}
*/
		// Test conservative constraint equation
		DataVector<double> dColumnSums;
		dColumnSums.Initialize(nAdjFaces);
/*
		for (int p = 0; p < nCoefficients; p++) {
			printf("%1.10e\n", dConstraint[p]);
		}
*/
/*
		for (int p = 0; p < nCoefficients; p++) {
		for (int i = 0; i < nAdjFaces; i++) {
			dColumnSums[i] += dConstraint[p] * dFitArrayPlus[i][p];
			printf("%1.3e  ", dFitArrayPlus[i][p]);
		}
		printf(";\n");
		}

		for (int i = 0; i < nAdjFaces; i++) {
			printf("%1.10e\n", dColumnSums[i]);
		}
		_EXCEPTION();
*/

		// Multiply integration array and fit array
		DataMatrix<double> dComposedArray;
		dComposedArray.Initialize(nAdjFaces, nOverlapFaces);

		for (int i = 0; i < nAdjFaces; i++) {
		for (int j = 0; j < nOverlapFaces; j++) {
		for (int k = 0; k < nCoefficients; k++) {
			dComposedArray[i][j] += dIntArray[k][j] * dFitArrayPlus[i][k];
		}
		}
		}
/*
		DataVector<double> dRowSums;
		dRowSums.Initialize(nOverlapFaces);

		DataVector<double> dColSums;
		dColSums.Initialize(nAdjFaces);

		for (int i = 0; i < nAdjFaces; i++) {
		for (int j = 0; j < nOverlapFaces; j++) {
			dRowSums[j] += dComposedArray[i][j] / meshOverlap.vecFaceArea[ixOverlap + j];
			dColSums[i] += dComposedArray[i][j];
		}
		}

		printf("\n");
		for (int j = 0; j < nOverlapFaces; j++) {
			printf("%1.15e\n", dRowSums[j]);
		}
		printf("\n");
		printf("%1.15e %1.15e\n", dColSums[0], meshInput.vecFaceArea[ixFirst]);
		_EXCEPTION();
*/
		// Put composed array into map
		for (int i = 0; i < vecAdjFaces.size(); i++) {
		for (int j = 0; j < nOverlapFaces; j++) {
			int ixFirstFace = vecAdjFaces[i].first;
			int ixSecondFace = meshOverlap.vecSecondFaceIx[ixOverlap + j];

			smatMap(ixSecondFace, ixFirstFace) +=
				dComposedArray[i][j]
				/ meshOutput.vecFaceArea[ixSecondFace];
		}
		}

		// Increment the current overlap index
		ixOverlap += nOverlapFaces;

		//_EXCEPTION();
	}
}

///////////////////////////////////////////////////////////////////////////////

void ForceIntArrayConsistencyConservation(
	const DataVector<double> & vecSourceArea,
	const DataVector<double> & vecTargetArea,
	DataMatrix<double> & dCoeff,
	bool fMonotone
) {

	// Number of free coefficients
	int nCoeff = dCoeff.GetRows() * dCoeff.GetColumns();

	// Number of conditions
	int nCondConservation = dCoeff.GetColumns();
	int nCondConsistency  = dCoeff.GetRows();

	// One condition is dropped due to linear dependence
	int nCond = nCondConservation + nCondConsistency - 1;

	// Product matrix
	DataMatrix<double> dCCt;
	dCCt.Initialize(nCond, nCond);

	DataMatrix<double> dC;
	dC.Initialize(nCoeff, nCond);

	DataVector<double> dRHS;
	dRHS.Initialize(nCoeff + nCond);

	// RHS
	int ix = 0;
	for (int i = 0; i < dCoeff.GetRows(); i++) {
	for (int j = 0; j < dCoeff.GetColumns(); j++) {
		dRHS[ix] = dCoeff[i][j];
		ix++;
	}
	}

	// Consistency
	ix = 0;
	for (int i = 0; i < dCoeff.GetRows(); i++) {
		for (int j = 0; j < dCoeff.GetColumns(); j++) {
			dC[i * dCoeff.GetColumns() + j][ix] = 1.0;
		}
		dRHS[nCoeff + ix] = 1.0;
		ix++;
	}

	// Conservation
	for (int j = 0; j < dCoeff.GetColumns()-1; j++) {
		for (int i = 0; i < dCoeff.GetRows(); i++) {
			dC[i * dCoeff.GetColumns() + j][ix] = vecTargetArea[i];
		}
		dRHS[nCoeff + ix] = vecSourceArea[j];
		ix++;
	}

	// Calculate CCt
	double dP = 0.0;
	for (int i = 0; i < dCoeff.GetRows(); i++) {
		dP += vecTargetArea[i] * vecTargetArea[i];
	}

	for (int i = 0; i < dCoeff.GetRows(); i++) {
		dCCt[i][i] = static_cast<double>(dCoeff.GetColumns());
		for (int j = 0; j < dCoeff.GetColumns()-1; j++) {
			dCCt[i][dCoeff.GetRows() + j] = vecTargetArea[i];
			dCCt[dCoeff.GetRows() + j][i] = vecTargetArea[i];
		}
	}

	for (int i = 0; i < dCoeff.GetColumns()-1; i++) {
		int ix = dCoeff.GetRows() + i;
		dCCt[ix][ix] = dP;
	}
/*
	for (int i = 0; i < nCond; i++) {
	for (int j = 0; j < nCond; j++) {
		for (int k = 0; k < nCoeff; k++) {
			dCCt[i][j] += dC[k][i] * dC[k][j];
		}
	}
	}
*/
/*
	FILE * fp = fopen("cct.dat", "w");
	for (int i = 0; i < nCond; i++) {
		for (int j = 0; j < nCond; j++) {
			fprintf(fp, "%1.15e\t", dCCt[i][j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	_EXCEPTION();
*/
	// Calculate C*r1 - r2
	char trans = 'n';
	int m = nCond;
	int n = nCoeff;
	int lda = m;
	int incx = 1;
	int incy = 1;
	double posone = 1.0;
	double negone = -1.0;
	dgemv_(
		&trans,
		&m,
		&n,
		&posone,
		&(dC[0][0]),
		&lda,
		&(dRHS[0]),
		&incx,
		&negone,
		&(dRHS[nCoeff]),
		&incy);

	// Solve the general system
	int nrhs = 1;
	int ldb = nCond;

	int nInfo;
/*
	DataVector<int> iPIV;
	iPIV.Initialize(nCond);

	dgesv_(
		&m,
		&nrhs,
		&(dCCt[0][0]),
		&lda,
		&(iPIV[0]),
		&(dRHS[nCoeff]),
		&ldb,
		&nInfo);
*/
	char uplo = 'U';
	dposv_(
		&uplo,
		&m,
		&nrhs,
		&(dCCt[0][0]),
		&lda,
		&(dRHS[nCoeff]),
		&ldb,
		&nInfo);

	if (nInfo != 0) {
		_EXCEPTION1("Unable to solve SPD Schur system: %i\n", nInfo);
	}

	// Obtain coefficients
	trans = 't';
	dgemv_(
		&trans,
		&m,
		&n,
		&negone,
		&(dC[0][0]),
		&lda,
		&(dRHS[nCoeff]),
		&incx,
		&posone,
		&(dRHS[0]),
		&incy);

	// Store coefficients in array
	ix = 0;
	for (int i = 0; i < dCoeff.GetRows(); i++) {
	for (int j = 0; j < dCoeff.GetColumns(); j++) {
		dCoeff[i][j] = dRHS[ix];
		ix++;
	}
	}

	// Force monotonicity
	if (fMonotone) {

		// Calculate total element Jacobian
		double dTotalJacobian = 0.0;
		for (int i = 0; i < vecSourceArea.GetRows(); i++) {
			dTotalJacobian += vecSourceArea[i];
		}

		// Determine low-order remap coefficients
		DataMatrix<double> dMonoCoeff;
		dMonoCoeff.Initialize(dCoeff.GetRows(), dCoeff.GetColumns());

		for (int i = 0; i < dCoeff.GetRows(); i++) {
		for (int j = 0; j < dCoeff.GetColumns(); j++) {
			dMonoCoeff[i][j] =
				vecSourceArea[j]
				/ dTotalJacobian;
		}
		}

		// Compute scaling factor
		double dA = 0.0;
		for (int i = 0; i < dCoeff.GetRows(); i++) {
		for (int j = 0; j < dCoeff.GetColumns(); j++) {
			if (dCoeff[i][j] < 0.0) {
				double dNewA =
					- dCoeff[i][j] / fabs(dMonoCoeff[i][j] - dCoeff[i][j]);

				if (dNewA > dA) {
					dA = dNewA;
				}
			}
		}
		}

		for (int i = 0; i < dCoeff.GetRows(); i++) {
		for (int j = 0; j < dCoeff.GetColumns(); j++) {
			dCoeff[i][j] = (1.0 - dA) * dCoeff[i][j] + dA * dMonoCoeff[i][j];
		}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void LinearRemapFVtoGLL(
	const Mesh & meshInput,
	const Mesh & meshOutput,
	const Mesh & meshOverlap,
	const DataMatrix3D<int> & dataGLLNodes,
	const DataMatrix3D<double> & dataGLLJacobian,
	int nOrder,
	OfflineMap & mapRemap,
	bool fMonotone
) {
	// Verify ReverseNodeArray has been calculated
	if (meshInput.revnodearray.size() == 0) {
		_EXCEPTIONT("ReverseNodeArray has not been calculated for meshInput");
	}

	// Triangular quadrature rule
	TriangularQuadratureRule triquadrule(4);

	const DataMatrix<double> & dG = triquadrule.GetG();
	const DataVector<double> & dW = triquadrule.GetW();

	// Get SparseMatrix represntation of the OfflineMap
	SparseMatrix<double> & smatMap = mapRemap.GetSparseMatrix();

	// Order of the finite element method
	int nP = dataGLLNodes.GetRows();

	// Sample coefficients
	DataMatrix<double> dSampleCoeff;
	dSampleCoeff.Initialize(nP, nP);

	// Number of elements needed
	int nCoefficients = nOrder * (nOrder + 1) / 2;

	int nRequiredFaceSetSize = nCoefficients;

	// Current overlap face
	int ixOverlap = 0;
/*
	// Generate the unique Jacobian for each point
	DataVector<double> dataUniqueJacobian;
	GenerateUniqueJacobian(dataGLLNodes, dataGLLJacobian, dataUniqueJacobian);
*/
	// Build the integration array for each element on meshOutput
	DataMatrix3D<double> dGlobalIntArray;
	dGlobalIntArray.Initialize(
		nCoefficients,
		meshOverlap.faces.size(),
		nP * nP);

	// Number of overlap Faces per source Face
	DataVector<int> nAllOverlapFaces;
	nAllOverlapFaces.Initialize(meshInput.faces.size());

	DataVector<int> nAllTotalOverlapTriangles;
	nAllTotalOverlapTriangles.Initialize(meshInput.faces.size());

	for (int ixFirst = 0; ixFirst < meshInput.faces.size(); ixFirst++) {

		int ixOverlapTemp = ixOverlap;
		for (; ixOverlapTemp < meshOverlap.faces.size(); ixOverlapTemp++) {

			const Face & faceOverlap = meshOverlap.faces[ixOverlapTemp];

			if (meshOverlap.vecFirstFaceIx[ixOverlapTemp] != ixFirst) {
				break;
			}

			nAllOverlapFaces[ixFirst]++;
			nAllTotalOverlapTriangles[ixFirst] += faceOverlap.edges.size() - 2;
		}

		// Increment the current overlap index
		ixOverlap += nAllOverlapFaces[ixFirst];
	}

	// Loop through all faces on meshInput
	ixOverlap = 0;

	for (int ixFirst = 0; ixFirst < meshInput.faces.size(); ixFirst++) {

		// Output every 100 elements
		if (ixFirst % 100 == 0) {
			Announce("Element %i", ixFirst);
		}

		// This Face
		const Face & faceFirst = meshInput.faces[ixFirst];

		// Area of the First Face
		double dFirstArea = meshInput.vecFaceArea[ixFirst];

		// Coordinate axes
		const Node & nodeRef = meshInput.nodes[faceFirst[0]];

		Node nodeA1 = meshInput.nodes[faceFirst[1]] - nodeRef;
		Node nodeA2 = meshInput.nodes[faceFirst[2]] - nodeRef;

		Node nodeC = CrossProduct(nodeA1, nodeA2);

		// Fit matrix
		DataMatrix<double> dFit;
		dFit.Initialize(3,3);

		dFit[0][0] = nodeA1.x; dFit[0][1] = nodeA1.y; dFit[0][2] = nodeA1.z;
		dFit[1][0] = nodeA2.x; dFit[1][1] = nodeA2.y; dFit[1][2] = nodeA2.z;
		dFit[2][0] = nodeC.x;  dFit[2][1] = nodeC.y;  dFit[2][2] = nodeC.z;

		// Number of overlapping Faces and triangles
		int nOverlapFaces = nAllOverlapFaces[ixFirst];
		int nTotalOverlapTriangles = nAllTotalOverlapTriangles[ixFirst];

		// Loop through all Overlap Faces
		for (int i = 0; i < nOverlapFaces; i++) {

			// Quantities from the overlap Mesh
			const Face & faceOverlap = meshOverlap.faces[ixOverlap + i];

			const NodeVector & nodesOverlap = meshOverlap.nodes;

			int nOverlapTriangles = faceOverlap.edges.size() - 2;

			// Quantities from the Second Mesh
			int ixSecond = meshOverlap.vecSecondFaceIx[ixOverlap + i];

			const NodeVector & nodesSecond = meshOutput.nodes;

			const Face & faceSecond = meshOutput.faces[ixSecond];

			// Loop over all sub-triangles of this Overlap Face
			for (int j = 0; j < nOverlapTriangles; j++) {

				// Cornerpoints of triangle
				const Node & node0 = nodesOverlap[faceOverlap[0]];
				const Node & node1 = nodesOverlap[faceOverlap[j+1]];
				const Node & node2 = nodesOverlap[faceOverlap[j+2]];

				// Calculate the area of the modified Face
				Face faceTri(3);
				faceTri.SetNode(0, faceOverlap[0]);
				faceTri.SetNode(1, faceOverlap[j+1]);
				faceTri.SetNode(2, faceOverlap[j+2]);

				double dTriArea =
					CalculateFaceArea(faceTri, nodesOverlap);

				for (int k = 0; k < triquadrule.GetPoints(); k++) {
					double * dGL = dG[k];

					// Get the nodal location of this point
					double dX[3];

					dX[0] = dGL[0] * node0.x + dGL[1] * node1.x + dGL[2] * node2.x;
					dX[1] = dGL[0] * node0.y + dGL[1] * node1.y + dGL[2] * node2.y;
					dX[2] = dGL[0] * node0.z + dGL[1] * node1.z + dGL[2] * node2.z;

					double dMag =
						sqrt(dX[0] * dX[0] + dX[1] * dX[1] + dX[2] * dX[2]);

					dX[0] /= dMag;
					dX[1] /= dMag;
					dX[2] /= dMag;

					Node nodeQuadrature(dX[0], dX[1], dX[2]);

					// Find the coefficients for this point
					int n = 3;
					int nrhs = 1;
					int lda = 3;
					int ipiv[3];
					int ldb = 3;
					int info;

					DataMatrix<double> dFitTemp;
					dFitTemp = dFit;
					dgesv_(
						&n, &nrhs, &(dFitTemp[0][0]), &lda, ipiv, dX, &ldb, &info);

					// Find the components of this quadrature point in the basis
					// of the finite element.
					double dAlpha;
					double dBeta;

					ApplyInverseMap(
						faceSecond,
						nodesSecond,
						nodeQuadrature,
						dAlpha,
						dBeta);

					// Check inverse map value
					if ((dAlpha < 0.0) || (dAlpha > 1.0) ||
						(dBeta  < 0.0) || (dBeta  > 1.0)
					) {
						_EXCEPTION2("Inverse Map out of range (%1.5e %1.5e)",
							dAlpha, dBeta);
					}

					// Sample this point in the GLL element
					int ixs = 0;
					for (int s = 0; s < nP; s++) {
					for (int t = 0; t < nP; t++) {

						// Sample the finite element at this point
						SampleGLLFiniteElement(
							fMonotone, nP,
							dAlpha,
							dBeta,
							dSampleCoeff);

						int ixp = 0;
						for (int p = 0; p < nOrder; p++) {
						for (int q = 0; q < nOrder - p; q++) {

							double dIntUpdate =
								  IPow(dX[0], p)
								* IPow(dX[1], q)
								* dW[k]
								* dSampleCoeff[s][t]
								* dTriArea;

							dGlobalIntArray[ixp][ixOverlap + i][ixs] +=
								dIntUpdate / dataGLLJacobian[s][t][ixSecond];

							ixp++;
						}
						}
						ixs++;
					}
					}
				}
			}
		}

		// Increment the current overlap index
		ixOverlap += nOverlapFaces;
	}

	// Reverse map
	std::vector< std::vector<int> > vecReverseFaceIx;
	vecReverseFaceIx.resize(meshOutput.faces.size());
	for (int i = 0; i < meshOverlap.faces.size(); i++) {
		int ixSecond = meshOverlap.vecSecondFaceIx[i];

		vecReverseFaceIx[ixSecond].push_back(i);
	}

	// Force consistency and conservation
	for (int ixSecond = 0; ixSecond < meshOutput.faces.size(); ixSecond++) {

		DataMatrix<double> dCoeff;
		dCoeff.Initialize(
			nP * nP,
			vecReverseFaceIx[ixSecond].size());

		for (int i = 0; i < vecReverseFaceIx[ixSecond].size(); i++) {
			int ixOverlap = vecReverseFaceIx[ixSecond][i];

			for (int s = 0; s < nP * nP; s++) {
				dCoeff[s][i] = dGlobalIntArray[0][ixOverlap][s];
			}
		}
/*
		for (int s = 0; s < nP * nP; s++) {
			double dConsistency = 0.0;
			for (int i = 0; i < dCoeff.GetRows(); i++) {
				dConsistency += dCoeff[i][s];
			}
			printf("%1.15e\n", dConsistency);
		}
*/

		// Target areas
		DataVector<double> vecTargetArea;
		vecTargetArea.Initialize(nP * nP);

		for (int i = 0; i < dCoeff.GetRows(); i++) {
			int ixOverlap = vecReverseFaceIx[ixSecond][i];

			for (int s = 0; s < nP * nP; s++) {
				vecTargetArea[s] =
					dataGLLJacobian[s/nP][s%nP][ixSecond];
					// meshOverlap.vecFaceArea[ixOverlap];
			}
		}

		// Source areas
		DataVector<double> vecSourceArea;
		vecSourceArea.Initialize(vecReverseFaceIx[ixSecond].size());

		for (int i = 0; i < vecReverseFaceIx[ixSecond].size(); i++) {
			int ixOverlap = vecReverseFaceIx[ixSecond][i];
			vecSourceArea[i] = meshOverlap.vecFaceArea[ixOverlap];
		}

		ForceIntArrayConsistencyConservation(
			vecSourceArea,
			vecTargetArea,
			dCoeff,
			fMonotone);

		for (int i = 0; i < dCoeff.GetRows(); i++) {
			double dConsistency = 0.0;
			for (int j = 0; j < dCoeff.GetColumns(); j++) {
				dConsistency += dCoeff[i][j];
			}
			//printf("%1.15e\n", dConsistency);
		}

		for (int i = 0; i < vecReverseFaceIx[ixSecond].size(); i++) {
			int ixOverlap = vecReverseFaceIx[ixSecond][i];

			for (int s = 0; s < nP * nP; s++) {
				//printf("%1.15e %1.15e\n", dGlobalIntArray[0][ixOverlap][s], dCoeff[s][i]);
				dGlobalIntArray[0][ixOverlap][s] = dCoeff[s][i];
			}
		}

/*
		for (int i = 0; i < dCoeff.GetRows(); i++) {
			double dConsistency = 0.0;
			for (int j = 0; j < dCoeff.GetColumns(); j++) {
				dConsistency += dCoeff[i][j];
			}
			printf("%1.15e\n", dConsistency);
		}

		for (int i = 0; i < dCoeff.GetRows(); i++) {
			int ixFirst = vecReverseFaceIx[ixSecond][i];

			for (int s = 0; s < dCoeff.GetColumns(); s++) {
				vecTargetArea[i] += dCoeff[i][s]
					* dataGLLJacobian[s/nP][s%nP][ixSecond]
					 meshInput.vecFaceArea[ixFirst];
			}
			printf("%1.15e\n", vecTargetArea[i]);
		}
*/
/*
		double dConsistency = 0.0;
		double dConservation = 0.0;

		int ixFirst = meshOverlap.vecFirstFaceIx[i];
		int ixSecond = meshOverlap.vecSecondFaceIx[i];

		for (int s = 0; s < nP * nP; s++) {
			//dConsistency += dGlobalIntArray[0][i][s];
			dConservation += dGlobalIntArray[0][i][s]
				* dataGLLJacobian[s/nP][s%nP][ixSecond]
				/ meshInput.vecFaceArea[ixFirst];

			printf("%1.15e\n", dataGLLJacobian[s/nP][s%nP][ixSecond]);
		}

		//printf("Consistency: %1.15e\n", dConsistency);
		printf("Conservation: %1.15e\n", dConservation);

		_EXCEPTION();
*/
	}
/*
	// Check consistency
	DataMatrix<double> dIntSums;
	dIntSums.Initialize(meshOutput.faces.size(), nP * nP);

	for (int i = 0; i < meshOverlap.faces.size(); i++) {
	for (int s = 0; s < nP * nP; s++) {
		int ixSecond = meshOverlap.vecSecondFaceIx[i];

		dIntSums[ixSecond][s] += dGlobalIntArray[0][i][s];
	}
	}
	for (int i = 0; i < dIntSums.GetRows(); i++) {
	for (int s = 0; s < dIntSums.GetColumns(); s++) {
		printf("%1.15e\n", dIntSums[i][s]);
	}
	}
*/
/*
	// Check conservation
	DataVector<double> dMassSums;
	dMassSums.Initialize(meshInput.faces.size());

	for (int i = 0; i < meshOverlap.faces.size(); i++) {
	for (int s = 0; s < nP * nP; s++) {
		int ixFirst = meshOverlap.vecFirstFaceIx[i];
		int ixSecond = meshOverlap.vecSecondFaceIx[i];

		dMassSums[ixFirst] += dGlobalIntArray[0][i][s]
			* dataGLLJacobian[s/nP][s%nP][ixSecond]
			/ meshInput.vecFaceArea[ixFirst];
	}
	}

	for (int i = 0; i < dMassSums.GetRows(); i++) {
		printf("%1.15e\n", dMassSums[i]);
	}
*/
/*
	DataVector<double> dOverlapMass;
	dOverlapMass.Initialize(meshOutput.faces.size());

	for (int i = 0; i < meshOverlap.faces.size(); i++) {
	for (int s = 0; s < nP * nP; s++) {
		dOverlapMass[i] += dGlobalIntArray[0][i][s];
	}
	}
*/
	// Impose conservative and consistent conditions on integration array
	//_EXCEPTION();

	// Construct finite-volume fit matrix and compose with integration operator
	ixOverlap = 0;

	for (int ixFirst = 0; ixFirst < meshInput.faces.size(); ixFirst++) {

		// Output every 100 elements
		if (ixFirst % 100 == 0) {
			Announce("Element %i", ixFirst);
		}

		// This Face
		const Face & faceFirst = meshInput.faces[ixFirst];

		// Area of the First Face
		double dFirstArea = meshInput.vecFaceArea[ixFirst];

		// Coordinate axes
		const Node & nodeRef = meshInput.nodes[faceFirst[0]];

		Node nodeA1 = meshInput.nodes[faceFirst[1]] - nodeRef;
		Node nodeA2 = meshInput.nodes[faceFirst[2]] - nodeRef;

		Node nodeC = CrossProduct(nodeA1, nodeA2);

		// Fit matrix
		DataMatrix<double> dFit;
		dFit.Initialize(3,3);

		dFit[0][0] = nodeA1.x; dFit[0][1] = nodeA1.y; dFit[0][2] = nodeA1.z;
		dFit[1][0] = nodeA2.x; dFit[1][1] = nodeA2.y; dFit[1][2] = nodeA2.z;
		dFit[2][0] = nodeC.x;  dFit[2][1] = nodeC.y;  dFit[2][2] = nodeC.z;

		// Number of overlapping Faces and triangles
		int nOverlapFaces = nAllOverlapFaces[ixFirst];
		int nTotalOverlapTriangles = nAllTotalOverlapTriangles[ixFirst];

		// Verify equal partition of mass in integration array
		double dTotal = 0.0;
		for (int i = 0; i < nOverlapFaces; i++) {
			int ixSecond = meshOverlap.vecSecondFaceIx[ixOverlap + i];

			for (int s = 0; s < nP * nP; s++) {
				dTotal += dGlobalIntArray[0][ixOverlap + i][s]
					* dataGLLJacobian[s/nP][s%nP][ixSecond]
					/ dFirstArea;
			}
		}
		if (fabs(dTotal - 1.0) > 1.0e-8) {
			printf("%1.15e\n", dTotal);
			_EXCEPTION();
		}

		// Determine the conservative constraint equation
		DataVector<double> dConstraint;
		dConstraint.Initialize(nCoefficients);

		for (int p = 0; p < nCoefficients; p++) {
		for (int i = 0; i < nOverlapFaces; i++) {
			int ixSecond = meshOverlap.vecSecondFaceIx[ixOverlap + i];

			for (int s = 0; s < nP * nP; s++) {
				dConstraint[p] += dGlobalIntArray[p][ixOverlap + i][s]
					* dataGLLJacobian[s/nP][s%nP][ixSecond]
					/ dFirstArea;
			}
		}
		}
/*
		for (int p = 0; p < nCoefficients; p++) {
		for (int j = 0; j < nOverlapFaces * nP * nP; j++) {
			dConstraint[p] += dIntArray[p][j] / dFirstArea;
		}
		}
*/
		// Set of Faces to use in building the reconstruction and associated
		// distance metric.
		AdjacentFaceVector vecAdjFaces;

		GetAdjacentFaceVector(
			meshInput,
			ixFirst,
			nRequiredFaceSetSize,
			vecAdjFaces);

		// Number of adjacent Faces
		int nAdjFaces = vecAdjFaces.size();

		// Least squares arrays
		DataMatrix<double> dFitArray;
		DataVector<double> dFitWeights;
		DataMatrix<double> dFitArrayPlus;

		BuildFitArray(
			meshInput,
			ixFirst,
			vecAdjFaces,
			nOrder,
			dConstraint,
			dFitArray,
			dFitWeights,
			dFitArrayPlus
		);

/*
		DataVector<double> dRowSum;
		dRowSum.Initialize(nCoefficients);

		for (int i = 0; i < nAdjFaces; i++) {
		for (int k = 0; k < nCoefficients; k++) {
			dRowSum[k] += dFitArrayPlus[i][k];
		}
		}

		for (int k = 0; k < nCoefficients; k++) {
			printf("%1.15e\n", dRowSum[k]);
		}
		_EXCEPTION();
*/

		// Multiply integration array and fit array
		DataMatrix<double> dComposedArray;
		dComposedArray.Initialize(nAdjFaces, nOverlapFaces * nP * nP);

		for (int j = 0; j < nOverlapFaces; j++) {
			int ixSecond = meshOverlap.vecSecondFaceIx[ixOverlap + j];

			for (int i = 0; i < nAdjFaces; i++) {
			for (int s = 0; s < nP * nP; s++) {
			for (int k = 0; k < nCoefficients; k++) {
				dComposedArray[i][j * nP * nP + s] +=
					dGlobalIntArray[k][ixOverlap + j][s]
					* dFitArrayPlus[i][k];
			}
			}
			}
		}

		// Put composed array into map
		for (int i = 0; i < vecAdjFaces.size(); i++) {
		for (int j = 0; j < nOverlapFaces; j++) {
			int ixFirstFace = vecAdjFaces[i].first;
			int ixSecondFace = meshOverlap.vecSecondFaceIx[ixOverlap + j];

			for (int s = 0; s < nP; s++) {
			for (int t = 0; t < nP; t++) {

				int jx = j * nP * nP + s * nP + t;

				int ixSecondNode = ixSecondFace * nP * nP + s * nP + t;//dataGLLNodes[s][t][ixSecondFace]-1;

				smatMap(ixSecondNode, ixFirstFace) +=
					dComposedArray[i][jx];
			}
			}
		}
		}

		// Increment the current overlap index
		ixOverlap += nOverlapFaces;
	}

/*
	double dTotal = 0.0;
	for (int i = 0; i < dNumericalSecondArea.GetRows(); i++) {
		if (dNumericalSecondArea[i] == 0.0) {
			printf("%i\n", i);
		}
		dTotal += dNumericalSecondArea[i];
	}
	printf("%1.15e\n", dTotal);
	_EXCEPTION();
*/
/*
	// Loop through all entries in the sparse matrix and divide by the
	// numerical area.
	DataVector<int> dataRows;
	DataVector<int> dataCols;
	DataVector<double> dataEntries;

	smatMap.GetEntries(dataRows, dataCols, dataEntries);

	for (int i = 0; i < dataEntries.GetRows(); i++) {
		int ixSecond = dataRows[i];
		dataEntries[i] /= dNumericalSecondArea[ixSecond];
	}

	smatMap.SetEntries(dataRows, dataCols, dataEntries);
*/
/*
	for (int i = 0; i < meshOutput.faces.size(); i++) {
		double dTotalArea = 0.0;

		int ixs = 0;
		for (int s = 0; s < nP; s++) {
		for (int t = 0; t < nP; t++) {
			dTotalArea += dNumericalSecondArea[i * nP * nP + ixs];

			ixs++;
		}
		}

		if (fabs(dTotalArea - meshOutput.vecFaceArea[i]) > 1.0e-12) {
			printf("%1.15e %1.15e\n", dTotalArea, meshOutput.vecFaceArea[i]);
		}
		
	}
	_EXCEPTION();
*/
}

///////////////////////////////////////////////////////////////////////////////

