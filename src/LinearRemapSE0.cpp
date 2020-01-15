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
#include "TriangularQuadrature.h"
#include "MeshUtilitiesFuzzy.h"

#include "Announce.h"

///////////////////////////////////////////////////////////////////////////////


extern "C" {
/*
	/// Routine from DQED.h to solve constrained bounded least squares problems
	void dbocls_(
		double * w,
		int * mdw,
		int * mcon,
		int * mrows,
		int * ncols,
		double * bl,
		double * bu,
		int * ind,
		int * iopt,
		double * x,
		double * rnormc,
		double * rnorm,
		int * mode,
		double * rw,
		int * iw);
*/
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
}

///////////////////////////////////////////////////////////////////////////////

void LinearRemapSE0(
	const Mesh & meshInput,
	const Mesh & meshOutput,
	const Mesh & meshOverlap,
	const DataArray3D<int> & dataGLLNodes,
	const DataArray3D<double> & dataGLLJacobian,
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
		ixCurrentFirstMeshFace = meshOverlap.vecSourceFaceIx[i];
		ixCurrentSecondMeshFace = meshOverlap.vecTargetFaceIx[i];

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

void ForceConsistencyConservation(
	const DataArray1D<double> & vecSourceArea,
	const DataArray1D<double> & vecTargetArea,
	DataArray2D<double> & dCoeff,
	bool fMonotone
) {
	int ncols = dCoeff.GetRows() * dCoeff.GetColumns();
	int mrows = dCoeff.GetRows() * dCoeff.GetColumns();
	//int mcon = dCoeff.GetRows() + dCoeff.GetColumns();
	int mcon = dCoeff.GetColumns();

	int mdw = mcon + mrows;

	// Allocate arrays
	DataArray2D<double> dW(ncols+mcon+1, mdw);

	DataArray1D<double> dBL(ncols + mcon);
	DataArray1D<double> dBU(ncols + mcon);
	DataArray1D<double> dX(2*(ncols+mcon)+2);
	DataArray1D<double> dRW(6 * ncols + 5 * mcon);
	DataArray1D<int> iInd(ncols + mcon);

	DataArray1D<int> iOpt(32);
	iOpt[0] = 99;

	DataArray1D<int> iIW(2*(ncols+mcon));

	double dRNormC = -1.0;
	double dRNorm = -1.0;
	int iMode = -1;

	// Matrix elements
	int ix = 0;
/*
	// Consistency
	for (int i = 0; i < dCoeff.GetRows(); i++) {
		for (int j = 0; j < dCoeff.GetColumns(); j++) {
			dW[i * dCoeff.GetColumns() + j][ix] = 1.0;
		}
		ix++;
	}
*/
	// Conservation
	for (int j = 0; j < dCoeff.GetColumns(); j++) {
		for (int i = 0; i < dCoeff.GetRows(); i++) {
			dW[i * dCoeff.GetColumns() + j][ix] = vecTargetArea[i];
		}
		ix++;
	}

	// Least squares problem:
	// Minimize L2 error between high-order coefficients and new coefficients
	ix = 0;

	for (int i = 0; i < dCoeff.GetRows(); i++) {
	for (int j = 0; j < dCoeff.GetColumns(); j++) {
		dW[ix][ix + mcon] = 1.0;
		dW[ncols][ix + mcon] = dCoeff[i][j];
		ix++;
	}
	}

	FILE * fp = fopen("w.txt", "w");
	for (int i = 0; i < dW.GetRows(); i++) {
		for (int j = 0; j < dW.GetColumns(); j++) {
			fprintf(fp, "%1.10e ", dW[i][j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	// Apply constraints
	ix = 0;

	// Monotonicity
	if (fMonotone) {
		for (int i = 0; i < ncols; i++) {
			dBL[ix] = 0.0;
			dBU[ix] = 1.0;
			iInd[ix] = 3;
			ix++;
		}

	} else {
		for (int i = 0; i < ncols; i++) {
			iInd[ix] = 4;
			ix++;
		}
	}
/*
	// Consistency
	for (int i = 0; i < dCoeff.GetRows(); i++) {
		dBL[ix] = 1.0;
		dBU[ix] = 1.0;
		iInd[ix] = 3;
		ix++;
	}
*/
	// Conservation
	for (int j = 0; j < dCoeff.GetColumns(); j++) {
		dBL[ix] = vecSourceArea[j];
		dBU[ix] = vecSourceArea[j];
		iInd[ix] = 3;
		ix++;
	}
/*
	for (int i = 0; i < dBL.GetRows(); i++) {
		printf("B%i: %1.15e %1.15e %i\n", i, dBL[i], dBU[i], iInd[i]);
	}
*/
	// Solve least squares problem
	DataArray2D<double> dWin = dW;
/*
	dbocls_(
		&(dW[0][0]),
		&mdw,
		&mcon,
		&mrows,
		&ncols,
		dBL,
		dBU,
		iInd,
		iOpt,
		dX,
		&dRNormC,
		&dRNorm,
		&iMode,
		dRW,
		iIW);
*/
	if (iMode >= 0) {
		printf("\nSUCCESS (%i)\n", iMode);
	} else {
		_EXCEPTION1("\nFAIL (%i)\n", iMode);
	}

	printf("\nNorm: %1.15e %1.15e\n", dRNorm, dRNormC);
/*
	ix = 0;
	for (int i = 0; i < dCoeff.GetColumns(); i++) {
		double dProduct = 0.0;
		for (int j = 0; j < dCoeff.GetColumns() * dCoeff.GetRows(); j++) {
			dProduct += dWin[j][i] * dX[j];
		}
		printf("%1.15e %1.15e\n", dProduct, vecSourceArea[i]);
	}

	for (int i = 0; i < dX.GetRows(); i++) {
		printf("X%i: %1.15e\n", i, dX[i]);
	}

	for (int i = 0; i < dBL.GetRows(); i++) {
		printf("B%i: %1.15e %1.15e %i\n", i, dBL[i], dBU[i], iInd[i]);
	}

	for (int i = 0; i < dCoeff.GetRows(); i++) {
	for (int j = 0; j < dCoeff.GetColumns(); j++) {
		//printf("%1.15e %1.15e\n", dCoeff[i][j], dX[ix]);
		printf("%1.15e, ", dX[ix]);
		ix++;
	}
	}
*/
}

///////////////////////////////////////////////////////////////////////////////

void ForceConsistencyConservation2(
	const DataArray1D<double> & vecSourceArea,
	const DataArray1D<double> & vecTargetArea,
	DataArray2D<double> & dCoeff,
	bool fMonotone
) {
	// Number of free coefficients
	int nCoeff = dCoeff.GetRows() * dCoeff.GetColumns();

	// Number of conditions
	int nCondConservation = dCoeff.GetColumns();
	int nCondConsistency  = dCoeff.GetRows();

	// One condition is dropped due to linear dependence
	int nCond = nCondConservation + nCondConsistency - 1;

	// The Lagrangian matrix and RHS vector
	DataArray2D<double> dLagrangian(nCoeff + nCond, nCoeff + nCond);

	DataArray1D<double> dRHS(nCoeff + nCond);

	// Build the Lagrangian matrix
	int ix = 0;

	for (int i = 0; i < dCoeff.GetRows(); i++) {
	for (int j = 0; j < dCoeff.GetColumns(); j++) {

		// Least squares conditions
		dLagrangian[ix][ix] = 1.0;
		dLagrangian[nCoeff + i][ix] = - 1.0;

		if (j != dCoeff.GetColumns()-1) {
			dLagrangian[nCoeff + nCondConsistency + j][ix] = - vecTargetArea[i];
		}

		// Consistency
		dLagrangian[ix][nCoeff + i] = -1.0;

		// Conservation
		if (j != dCoeff.GetColumns()-1) {
			dLagrangian[ix][nCoeff + nCondConsistency + j] = - vecTargetArea[i];
		}

		// RHS vector
		dRHS[ix] = dCoeff[i][j];

		ix++;
	}
	}

	FILE * fp = fopen("lagrangian.dat", "w");
	for (int i = 0; i < dLagrangian.GetRows(); i++) {
		for (int j = 0; j < dLagrangian.GetColumns(); j++) {
			fprintf(fp, "%1.15e\t", dLagrangian[i][j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	// Consistency and conservation of RHS vector
	for (int i = 0; i < dCoeff.GetRows(); i++) {
		dRHS[nCoeff + i] = -1.0;
	}
	for (int i = 0; i < dCoeff.GetColumns()-1; i++) {
		dRHS[nCoeff + nCondConsistency + i] = - vecSourceArea[i];
	}

	// Solve the system
	DataArray1D<int> iPIV(dLagrangian.GetRows());

	int n = dLagrangian.GetRows();
	int nRHS = 1;
	int nLDA = n;
	int nLDB = n;
	int nInfo;
/*
	dgesv_(
		&n,
		&nRHS,
		&(dLagrangian[0][0]),
		&nLDA,
		&(iPIV[0]),
		&(dRHS[0]),
		&nLDB,
		&nInfo);
*/
	char uplo = 'U';

	int lWork = n;
	DataArray1D<double> dWork(lWork);

	dsysv_(
		&uplo,
		&n,
		&nRHS,
		&(dLagrangian[0][0]),
		&nLDA,
		&(iPIV[0]),
		&(dRHS[0]),
		&nLDB,
		&(dWork[0]),
		&lWork,
		&nInfo);

	if (nInfo != 0) {
		_EXCEPTION1("Cannot solve target system: %i", nInfo);
	}

	/*
	printf("Info: %i\n", nInfo);

	// Output the result
	ix = 0;
	for (int i = 0; i < dCoeff.GetRows(); i++) {
	for (int j = 0; j < dCoeff.GetColumns(); j++) {
		printf("%1.15e %1.15e\n", dCoeff[i][j], dRHS[ix]);
		ix++;
	}
	}

	// Verify consistency
	ix = 0;
	for (int i = 0; i < dCoeff.GetRows(); i++) {
		double dRowSum = 0.0;
		for (int j = 0; j < dCoeff.GetColumns(); j++) {
			dRowSum += dRHS[ix];
			ix++;
		}
		printf("Row %i: %1.15e\n", i, dRowSum);
	}

	// Verify conservation
	ix = 0;
	DataArray1D<double> dColumnSums(dCoeff.GetColumns());

	for (int i = 0; i < dCoeff.GetRows(); i++) {
		double dRowSum = 0.0;
		for (int j = 0; j < dCoeff.GetColumns(); j++) {
			dColumnSums[j] += dRHS[ix] * vecTargetArea[i];
			ix++;
		}
	}
	for (int j = 0; j < dCoeff.GetColumns(); j++) {
		printf("Column %i: %1.15e %1.15e\n", j, dColumnSums[j], vecSourceArea[j]);
	}
	*/

	// Store coefficients in array
	ix = 0;
	for (int i = 0; i < dCoeff.GetRows(); i++) {
	for (int j = 0; j < dCoeff.GetColumns(); j++) {
		dCoeff[i][j] = dRHS[ix];
		ix++;
	}
	}
}

///////////////////////////////////////////////////////////////////////////////

void ForceConsistencyConservation3(
	const DataArray1D<double> & vecSourceArea,
	const DataArray1D<double> & vecTargetArea,
	DataArray2D<double> & dCoeff,
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
	DataArray2D<double> dCCt(nCond, nCond);
	DataArray2D<double> dC(nCoeff, nCond);

	DataArray1D<double> dRHS(nCoeff + nCond);

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
	DataArray1D<int> iPIV(nCond);

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
		DataArray2D<double> dMonoCoeff(dCoeff.GetRows(), dCoeff.GetColumns());

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

void LinearRemapSE4(
	const Mesh & meshInput,
	const Mesh & meshOutput,
	const Mesh & meshOverlap,
	const DataArray3D<int> & dataGLLNodes,
	const DataArray3D<double> & dataGLLJacobian,
	int nMonotoneType,
	bool fContinuousIn,
	bool fNoConservation,
	OfflineMap & mapRemap
) {
	// Order of the polynomial interpolant
	int nP = dataGLLNodes.GetRows();

	// Triangular quadrature rule
	TriangularQuadratureRule triquadrule(4);

	int TriQuadraturePoints = triquadrule.GetPoints();

	const DataArray2D<double> & TriQuadratureG = triquadrule.GetG();

	const DataArray1D<double> & TriQuadratureW = triquadrule.GetW();

	// Sample coefficients
	DataArray2D<double> dSampleCoeff(nP, nP);

	// GLL Quadrature nodes on quadrilateral elements
	DataArray1D<double> dG;
	DataArray1D<double> dW;
	GaussLobattoQuadrature::GetPoints(nP, 0.0, 1.0, dG, dW);

	// Get SparseMatrix represntation of the OfflineMap
	SparseMatrix<double> & smatMap = mapRemap.GetSparseMatrix();

	// NodeVector from meshOverlap
	const NodeVector & nodesOverlap = meshOverlap.nodes;
	const NodeVector & nodesFirst   = meshInput.nodes;

	// Vector of source areas
	DataArray1D<double> vecSourceArea(nP * nP);

	DataArray1D<double> vecTargetArea;

	DataArray2D<double> dCoeff;

	// Current Overlap Face
	int ixOverlap = 0;

	// Loop over all input Faces
	for (int ixFirst = 0; ixFirst < meshInput.faces.size(); ixFirst++) {

		const Face & faceFirst = meshInput.faces[ixFirst];

		if (faceFirst.edges.size() != 4) {
			_EXCEPTIONT("Only quadrilateral elements allowed for SE remapping");
		}

		// Output every 100 elements
		if (ixFirst % 1000 == 0) {
			Announce("Element %i/%i", ixFirst,meshInput.faces.size());
		}

		// Number of overlapping Faces and triangles
		int nOverlapFaces = 0;
		int nTotalOverlapTriangles = 0;

		// Determine how many overlap Faces and triangles are present
		int ixOverlapTemp = ixOverlap;
		for (; ixOverlapTemp < meshOverlap.faces.size(); ixOverlapTemp++) {

			const Face & faceOverlap = meshOverlap.faces[ixOverlapTemp];

			if (meshOverlap.vecSourceFaceIx[ixOverlapTemp] != ixFirst) {
				break;
			}

			nOverlapFaces++;
			nTotalOverlapTriangles += faceOverlap.edges.size() - 2;
		}

		// No overlaps
		if (nOverlapFaces == 0) {
			continue;
		}

		// Allocate remap coefficients array for meshFirst Face
		DataArray3D<double> dRemapCoeff(nP, nP, nOverlapFaces);

		// Find the local remap coefficients
		for (int j = 0; j < nOverlapFaces; j++) {
			const Face & faceOverlap = meshOverlap.faces[ixOverlap + j];

			int nOverlapTriangles = faceOverlap.edges.size() - 2;

			// Loop over all sub-triangles of this Overlap Face
			for (int k = 0; k < nOverlapTriangles; k++) {

				// Cornerpoints of triangle
				const Node & node0 = nodesOverlap[faceOverlap[0]];
				const Node & node1 = nodesOverlap[faceOverlap[k+1]];
				const Node & node2 = nodesOverlap[faceOverlap[k+2]];

				// Calculate the area of the modified Face
				Face faceTri(3);
				faceTri.SetNode(0, faceOverlap[0]);
				faceTri.SetNode(1, faceOverlap[k+1]);
				faceTri.SetNode(2, faceOverlap[k+2]);

				double dTriangleArea =
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
					double dAlpha;
					double dBeta;

					ApplyInverseMap(
						faceFirst,
						nodesFirst,
						nodeQuadrature,
						dAlpha,
						dBeta);

					// Check inverse map value
					if ((dAlpha < -InverseMapTolerance)      ||
						(dAlpha > 1.0 + InverseMapTolerance) ||
						(dBeta  < -InverseMapTolerance)      ||
						(dBeta  > 1.0 + InverseMapTolerance)
					) {
						printf("\n==== BEGIN DEBUGGING INFO ====\n");
						printf("WARNING (%s, Line %u) Inverse map out of range\n",
							__FILE__, __LINE__);
						int ixSecond = meshOverlap.vecTargetFaceIx[ixOverlap + j];
						printf("Source face ix %i, Target face ix %i, Overlap face ix %i\n",
							ixFirst, ixSecond, ixOverlap + j);
						printf("Face nodes:\n");
						for (int x = 0; x < faceFirst.edges.size(); x++) {
							nodesFirst[faceFirst[x]].Print("");
						}
						printf("Quadrature node:\n");
						nodeQuadrature.Print("");
						printf("Alpha, Beta: %1.15e %1.15e\n", dAlpha, dBeta);
						printf("==== END DEBUGGING INFO ====\n");
						//_EXCEPTION2("Inverse Map out of range (%1.5e %1.5e)",
						//	dAlpha, dBeta);
					}

					// Sample the finite element at this point
					SampleGLLFiniteElement(
						nMonotoneType,
						nP,
						dAlpha,
						dBeta,
						dSampleCoeff);

					// Add sample coefficients to the map
					for (int p = 0; p < nP; p++) {
					for (int q = 0; q < nP; q++) {

						dRemapCoeff[p][q][j] +=
							TriQuadratureW[l]
							* dTriangleArea
							* dSampleCoeff[p][q]
							/ meshOverlap.vecFaceArea[ixOverlap + j];
					}
					}
				}
			}
		}

		// Force consistency and conservation
		if (!fNoConservation) {
			double dTargetArea = 0.0;
			for (int j = 0; j < nOverlapFaces; j++) {
				dTargetArea += meshOverlap.vecFaceArea[ixOverlap + j];
			}

			for (int p = 0; p < nP; p++) {
			for (int q = 0; q < nP; q++) {
				vecSourceArea[p * nP + q] = dataGLLJacobian[p][q][ixFirst];
			}
			}

			// Source elements are completely covered by target volumes
			if (fabs(meshInput.vecFaceArea[ixFirst] - dTargetArea) <= HighTolerance) {
				vecTargetArea.Allocate(nOverlapFaces);
				for (int j = 0; j < nOverlapFaces; j++) {
					vecTargetArea[j] = meshOverlap.vecFaceArea[ixOverlap + j];
				}

				dCoeff.Allocate(nOverlapFaces, nP * nP);

				for (int j = 0; j < nOverlapFaces; j++) {
				for (int p = 0; p < nP; p++) {
				for (int q = 0; q < nP; q++) {
					dCoeff[j][p * nP + q] = dRemapCoeff[p][q][j];
				}
				}
				}

			// Target volumes only partially cover source elements
			} else if (meshInput.vecFaceArea[ixFirst] - dTargetArea > HighTolerance) {
				double dExtraneousArea = meshInput.vecFaceArea[ixFirst] - dTargetArea;

				vecTargetArea.Allocate(nOverlapFaces+1);
				for (int j = 0; j < nOverlapFaces; j++) {
					vecTargetArea[j] = meshOverlap.vecFaceArea[ixOverlap + j];
				}
				vecTargetArea[nOverlapFaces] = dExtraneousArea;

				Announce("Partial volume: %i (%1.10e / %1.10e)",
					ixFirst, dTargetArea, meshInput.vecFaceArea[ixFirst]);

				if (dTargetArea > meshInput.vecFaceArea[ixFirst]) {
					_EXCEPTIONT("Partial element area exceeds total element area");
				}

				dCoeff.Allocate(nOverlapFaces+1, nP * nP);

				for (int j = 0; j < nOverlapFaces; j++) {
				for (int p = 0; p < nP; p++) {
				for (int q = 0; q < nP; q++) {
					dCoeff[j][p * nP + q] = dRemapCoeff[p][q][j];
				}
				}
				}
				for (int p = 0; p < nP; p++) {
				for (int q = 0; q < nP; q++) {
					dCoeff[nOverlapFaces][p * nP + q] =
						dataGLLJacobian[p][q][ixFirst];
				}
				}
				for (int j = 0; j < nOverlapFaces; j++) {
				for (int p = 0; p < nP; p++) {
				for (int q = 0; q < nP; q++) {
					dCoeff[nOverlapFaces][p * nP + q] -=
						dRemapCoeff[p][q][j]
						* meshOverlap.vecFaceArea[ixOverlap + j];
				}
				}
				}
				for (int p = 0; p < nP; p++) {
				for (int q = 0; q < nP; q++) {
					dCoeff[nOverlapFaces][p * nP + q] /= dExtraneousArea;
				}
				}

			// Source elements only partially cover target volumes
			} else {
				printf("\n==== BEGIN DEBUGGING INFO ====\n");
				printf("EXCEPTION (%s, Line %u) Target grid must be a subset of source grid\n",
					__FILE__, __LINE__);
				printf("Source face ix %i, Overlap face ix [%i,%i]\n",
					ixFirst, ixOverlap, ixOverlap+nOverlapFaces-1);
				printf("Target faces / overlap area:\n");
				for (int j = 0; j < nOverlapFaces; j++) {
					printf("  (%i) (%i) %1.15e\n",
						ixOverlap + j,
						meshOverlap.vecTargetFaceIx[ixOverlap + j],
						meshOverlap.vecFaceArea[ixOverlap + j]);
				}
				printf("Source nodes / source area:\n");
				for (int p = 0; p < nP; p++) {
				for (int q = 0; q < nP; q++) {
					printf("  (%i,%i) %1.15e\n", p, q, dataGLLJacobian[p][q][ixFirst]);
				}
				}
				printf("==== END DEBUGGING INFO ====\n");

				_EXCEPTION2("Target grid must be a subset of source grid:"
					"\nInput mesh area (%1.15e) Target area (%1.15e)",
					meshInput.vecFaceArea[ixFirst],
					dTargetArea);
			}

			ForceConsistencyConservation3(
				vecSourceArea,
				vecTargetArea,
				dCoeff,
				(nMonotoneType != 0));

			for (int j = 0; j < nOverlapFaces; j++) {
			for (int p = 0; p < nP; p++) {
			for (int q = 0; q < nP; q++) {
				dRemapCoeff[p][q][j] = dCoeff[j][p * nP + q];
			}
			}
			}
		}

		// Put these remap coefficients into the SparseMatrix map
		for (int j = 0; j < nOverlapFaces; j++) {
			int ixSecondFace = meshOverlap.vecTargetFaceIx[ixOverlap + j];

			for (int p = 0; p < nP; p++) {
			for (int q = 0; q < nP; q++) {

				if (dRemapCoeff[p][q][j] == 0.0) {
					continue;
				}

				if (fContinuousIn) {
					int ixFirstNode = dataGLLNodes[p][q][ixFirst] - 1;

					smatMap(ixSecondFace, ixFirstNode) +=
						dRemapCoeff[p][q][j]
						* meshOverlap.vecFaceArea[ixOverlap + j]
						/ meshOutput.vecFaceArea[ixSecondFace];

				} else {
					int ixFirstNode = ixFirst * nP * nP + p * nP + q;

					smatMap(ixSecondFace, ixFirstNode) +=
						dRemapCoeff[p][q][j]
						* meshOverlap.vecFaceArea[ixOverlap + j]
						/ meshOutput.vecFaceArea[ixSecondFace];
				}
			}
			}
		}

		// Increment the current overlap index
		ixOverlap += nOverlapFaces;
	}
}

///////////////////////////////////////////////////////////////////////////////

void LinearRemapGLLtoGLL_Pointwise(
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
) {

	// Triangular quadrature rule
	TriangularQuadratureRule triquadrule(4);

	const DataArray2D<double> & dG = triquadrule.GetG();
	const DataArray1D<double> & dW = triquadrule.GetW();

	// Sample coefficients
	DataArray2D<double> dSampleCoeffIn(nPin, nPin);
	DataArray2D<double> dSampleCoeffOut(nPout, nPout);

	// Build the integration array for each element on meshOverlap
	DataArray3D<double> dGlobalIntArray(
		nPout * nPout,
		meshOverlap.faces.size(),
		nPin * nPin);

	DataArray2D<double> dConstantIntArray(
		meshOverlap.faces.size(),
		nPin * nPin);

	// Number of overlap Faces per source Face
	DataArray1D<int> nAllOverlapFaces(meshInput.faces.size());

	int ixOverlap = 0;

	for (int ixFirst = 0; ixFirst < meshInput.faces.size(); ixFirst++) {

		int ixOverlapTemp = ixOverlap;
		for (; ixOverlapTemp < meshOverlap.faces.size(); ixOverlapTemp++) {

			const Face & faceOverlap = meshOverlap.faces[ixOverlapTemp];

			if (meshOverlap.vecSourceFaceIx[ixOverlapTemp] != ixFirst) {
				break;
			}

			nAllOverlapFaces[ixFirst]++;
		}

		// Increment the current overlap index
		ixOverlap += nAllOverlapFaces[ixFirst];
	}

	// Mesh utilities
	MeshUtilitiesFuzzy meshutil;

	// Construct overlap
	ixOverlap = 0;

	int ixTotal = 0;

	std::set< std::pair<int, int> > setFound;

	for (int ixFirst = 0; ixFirst < meshInput.faces.size(); ixFirst++) {

		// Output every 100 elements
		if (ixFirst % 1000 == 0) {
                  Announce("Element %i/%i", ixFirst,meshInput.faces.size());
		}

		// Quantities from the First Mesh
		const Face & faceFirst = meshInput.faces[ixFirst];

		const NodeVector & nodesFirst = meshInput.nodes;

		// Number of overlapping Faces and triangles
		int nOverlapFaces = nAllOverlapFaces[ixFirst];

		// Loop through all Overlap Faces
		for (int i = 0; i < nOverlapFaces; i++) {

			// Quantities from the overlap Mesh
			const Face & faceOverlap = meshOverlap.faces[ixOverlap + i];

			const NodeVector & nodesOverlap = meshOverlap.nodes;

			int nOverlapTriangles = faceOverlap.edges.size() - 2;

			// Quantities from the Second Mesh
			int ixSecond = meshOverlap.vecTargetFaceIx[ixOverlap + i];

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

					// Get the nodal location of this point
					double dX[3];

					dX[0] = dG(k,0) * node0.x + dG(k,1) * node1.x + dG(k,2) * node2.x;
					dX[1] = dG(k,0) * node0.y + dG(k,1) * node1.y + dG(k,2) * node2.y;
					dX[2] = dG(k,0) * node0.z + dG(k,1) * node1.z + dG(k,2) * node2.z;

					double dMag =
						sqrt(dX[0] * dX[0] + dX[1] * dX[1] + dX[2] * dX[2]);

					dX[0] /= dMag;
					dX[1] /= dMag;
					dX[2] /= dMag;

					Node nodeQuadrature(dX[0], dX[1], dX[2]);

					// Find the components of this quadrature point in the basis
					// of the first Face.
					double dAlphaIn;
					double dBetaIn;

					ApplyInverseMap(
						faceFirst,
						nodesFirst,
						nodeQuadrature,
						dAlphaIn,
						dBetaIn);

					// Find the components of this quadrature point in the basis
					// of the second Face.
					double dAlphaOut;
					double dBetaOut;

					ApplyInverseMap(
						faceSecond,
						nodesSecond,
						nodeQuadrature,
						dAlphaOut,
						dBetaOut);

					// Check inverse map value
					if ((dAlphaIn < -InverseMapTolerance)      ||
						(dAlphaIn > 1.0 + InverseMapTolerance) ||
						(dBetaIn  < -InverseMapTolerance)      ||
						(dBetaIn  > 1.0 + InverseMapTolerance)
					) {
						printf("\n==== BEGIN DEBUGGING INFO ====\n");
						printf("WARNING (%s, Line %u) Inverse map out of range\n",
							__FILE__, __LINE__);
						printf("Source face ix %i, Target face ix %i, Overlap face ix %i\n",
							ixFirst, ixSecond, ixOverlap + i);
						printf("Face nodes:\n");
						for (int x = 0; x < faceFirst.edges.size(); x++) {
							nodesFirst[faceFirst[x]].Print("");
						}
						printf("Quadrature node:\n");
						nodeQuadrature.Print("");
						printf("Alpha, Beta: %1.15e %1.15e\n", dAlphaIn, dBetaIn);
						printf("==== END DEBUGGING INFO ====\n");
						//_EXCEPTION2("Inverse Map out of range (%1.5e %1.5e)",
						//	dAlphaIn, dBetaIn);
					}

					// Check inverse map value
					if ((dAlphaOut < -InverseMapTolerance)      ||
						(dAlphaOut > 1.0 + InverseMapTolerance) ||
						(dBetaOut  < -InverseMapTolerance)      ||
						(dBetaOut  > 1.0 + InverseMapTolerance)
					) {
						printf("\n==== BEGIN DEBUGGING INFO ====\n");
						printf("WARNING (%s, Line %u) Inverse map out of range\n",
							__FILE__, __LINE__);
						printf("Source face ix %i, Target face ix %i, Overlap face ix %i\n",
							ixFirst, ixSecond, ixOverlap + i);
						printf("Face nodes:\n");
						for (int x = 0; x < faceSecond.edges.size(); x++) {
							nodesSecond[faceSecond[x]].Print("");
						}
						printf("Quadrature node:\n");
						nodeQuadrature.Print("");
						printf("Alpha, Beta: %1.15e %1.15e\n", dAlphaOut, dBetaOut);
						printf("==== END DEBUGGING INFO ====\n");
						//_EXCEPTION2("Inverse Map out of range (%1.5e %1.5e)",
						//	dAlphaOut, dBetaOut);
					}

					// Sample the First finite element at this point
					SampleGLLFiniteElement(
						nMonotoneType,
						nPin,
						dAlphaIn,
						dBetaIn,
						dSampleCoeffIn);

					// Sample the Second finite element at this point
					SampleGLLFiniteElement(
						nMonotoneType,
						nPout,
						dAlphaOut,
						dBetaOut,
						dSampleCoeffOut);

					// Compute overlap integral
					int ixs = 0;
					for (int s = 0; s < nPin; s++) {
					for (int t = 0; t < nPin; t++) {
						dConstantIntArray[ixOverlap + i][ixs] +=
							  dSampleCoeffIn[s][t]
							* dW[k]
							* dTriArea
							/ meshOverlap.vecFaceArea[ixOverlap + i];
					}
					}
				}
			}

			{

				// Sample pointwise
				DataArray1D<double> dGL;
				DataArray1D<double> dWL;

				GaussLobattoQuadrature::GetPoints(nPout, 0.0, 1.0, dGL, dWL);

				int ixp = 0;
				for (int p = 0; p < nPout; p++) {
				for (int q = 0; q < nPout; q++) {

					Node node;
					Node dDx1G;
					Node dDx2G;

					ApplyLocalMap(
						faceSecond,
						nodesSecond,
						dGL[p],
						dGL[q],
						node,
						dDx1G,
						dDx2G);

					Face::NodeLocation loc;
					int ixLocation;

					meshutil.ContainsNode(
						faceFirst,
						nodesFirst,
						node,
						loc,
						ixLocation);

					if (loc == Face::NodeLocation_Exterior) {
						ixp++;
						continue;
					}

					std::pair<int,int> pairFound(ixSecond, p * nPout + q);

					if (setFound.find(pairFound) != setFound.end()) {
						ixp++;
						continue;
					} else {
						setFound.insert(pairFound);
					}

					ixTotal++;

					// Find the components of this quadrature point in the
					// basis of the input Face.
					double dAlphaIn;
					double dBetaIn;

					ApplyInverseMap(
						faceFirst,
						nodesFirst,
						node,
						dAlphaIn,
						dBetaIn);

					// Check inverse map value
					if ((dAlphaIn < -InverseMapTolerance)      ||
						(dAlphaIn > 1.0 + InverseMapTolerance) ||
						(dBetaIn  < -InverseMapTolerance)      ||
						(dBetaIn  > 1.0 + InverseMapTolerance)
					) {
						printf("\n==== BEGIN DEBUGGING INFO ====\n");
						printf("WARNING (%s, Line %u) Inverse map out of range\n",
							__FILE__, __LINE__);
						printf("Source face ix %i, Target face ix %i, Overlap face ix %i\n",
							ixFirst, ixSecond, ixOverlap + i);
						printf("Face nodes:\n");
						for (int x = 0; x < faceFirst.edges.size(); x++) {
							nodesFirst[faceFirst[x]].Print("");
						}
						printf("Quadrature node:\n");
						node.Print("");
						printf("Alpha, Beta: %1.15e %1.15e\n", dAlphaIn, dBetaIn);
						printf("==== END DEBUGGING INFO ====\n");
						//_EXCEPTION2("Inverse Map out of range (%1.5e %1.5e)",
						//	dAlphaIn, dBetaIn);
					}

					// Sample the First finite element at this point
					SampleGLLFiniteElement(
						nMonotoneType,
						nPin,
						dAlphaIn,
						dBetaIn,
						dSampleCoeffIn);

					int ixs = 0;
					for (int s = 0; s < nPin; s++) {
					for (int t = 0; t < nPin; t++) {

						if (dGlobalIntArray[ixp][ixOverlap + i][ixs] != 0.0) {
							_EXCEPTION();
						}

						dGlobalIntArray[ixp][ixOverlap + i][ixs] +=
							dSampleCoeffIn[s][t];

						ixs++;
					}
					}

					ixp++;
				}
				}

			}

		}

		// Force consistency and conservation of ConstantIntArray
		DataArray1D<double> dSourceArea(nPin * nPin);

		for (int s = 0; s < nPin; s++) {
		for (int t = 0; t < nPin; t++) {
			dSourceArea[s * nPin + t] = dataGLLJacobianIn[s][t][ixFirst];
		}
		}

		DataArray1D<double> dTargetArea(nOverlapFaces);
		for (int i = 0; i < nOverlapFaces; i++) {
			dTargetArea[i] = meshOverlap.vecFaceArea[ixOverlap + i];
		}

		DataArray2D<double> dCoeff(nOverlapFaces, nPin * nPin);

		for (int i = 0; i < nOverlapFaces; i++) {
		for (int s = 0; s < nPin; s++) {
		for (int t = 0; t < nPin; t++) {
			dCoeff[i][s * nPin + t] =
				dConstantIntArray[ixOverlap + i][s * nPin + t];
		}
		}
		}

		ForceConsistencyConservation3(
			dSourceArea,
			dTargetArea,
			dCoeff,
			(nMonotoneType != 0));

		for (int i = 0; i < nOverlapFaces; i++) {
			int ixSecondFace = meshOverlap.vecTargetFaceIx[ixOverlap + i];

			for (int s = 0; s < nPin; s++) {
			for (int t = 0; t < nPin; t++) {
				dConstantIntArray[ixOverlap + i][s * nPin + t] =
					dCoeff[i][s * nPin + t]
					* meshOverlap.vecFaceArea[ixOverlap + i]
					/ meshOutput.vecFaceArea[ixSecondFace];

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
		int ixSecond = meshOverlap.vecTargetFaceIx[i];

		vecReverseFaceIx[ixSecond].push_back(i);
	}

	// Force consistency and conservation on linear sub-map
	int ixDistBegin = 0;

	for (int ixSecond = 0; ixSecond < meshOutput.faces.size(); ixSecond++) {

		// Coefficients
		DataArray2D<double> dCoeff(
			nPout * nPout,
			vecReverseFaceIx[ixSecond].size() * nPin * nPin);

		for (int i = 0; i < vecReverseFaceIx[ixSecond].size(); i++) {

			int ixOverlap = vecReverseFaceIx[ixSecond][i];

			int ixs = 0;
			for (int s = 0; s < nPin; s++) {
			for (int t = 0; t < nPin; t++) {

				int ixp = 0;
				for (int p = 0; p < nPout; p++) {
				for (int q = 0; q < nPout; q++) {
					dCoeff[ixp][i * nPin * nPin + ixs] =
						dGlobalIntArray[ixp][ixOverlap][ixs];

					ixp++;
				}
				}

				ixs++;
			}
			}
		}

		// Source areas
		DataArray1D<double> dSourceArea(
			vecReverseFaceIx[ixSecond].size() * nPin * nPin);

		for (int i = 0; i < vecReverseFaceIx[ixSecond].size(); i++) {

			int ixOverlap = vecReverseFaceIx[ixSecond][i];

			int ixs = 0;
			for (int s = 0; s < nPin; s++) {
			for (int t = 0; t < nPin; t++) {
				dSourceArea[i * nPin * nPin + ixs] =
					dConstantIntArray[ixOverlap][ixs];

				ixs++;
			}
			}
		}

		// Target areas
		double dTotalTargetArea = 0.0;

		DataArray1D<double> dTargetArea(nPout * nPout);

		{
			int ixp = 0;
			for (int p = 0; p < nPout; p++) {
			for (int q = 0; q < nPout; q++) {

				dTargetArea[ixp] =
					dataGLLJacobianOut[p][q][ixSecond]
					/ meshOutput.vecFaceArea[ixSecond];

				ixp++;
			}
			}
		}

		// Force consistency and conservation
		ForceConsistencyConservation3(
			dSourceArea,
			dTargetArea,
			dCoeff,
			(nMonotoneType != 0));

/*
		// Check column sums (conservation)
		for (int i = 0; i < dCoeff.GetColumns(); i++) {
			double dColSum = 0.0;
			for (int j = 0; j < dCoeff.GetRows(); j++) {
				dColSum += dCoeff[j][i] * dTargetArea[j];
			}
			printf("Col %i: %1.15e\n", i, dColSum / dSourceArea[i]);
		}

		// Check row sums (consistency)
		for (int j = 0; j < dCoeff.GetRows(); j++) {
			double dRowSum = 0.0;
			for (int i = 0; i < dCoeff.GetColumns(); i++) {
				dRowSum += dCoeff[j][i];
			}
			printf("Row %i: %1.15e\n", j, dRowSum);
		}
		_EXCEPTION();
*/
/*
		// Coefficients
		for (int i = 0; i < vecReverseFaceIx[ixSecond].size(); i++) {

			ixOverlap = vecReverseFaceIx[ixSecond][i];

			if ((ixOverlap < 0) || (ixOverlap > dGlobalIntArray.GetColumns())) {
				_EXCEPTION();
			}

			int ixs = 0;
			for (int s = 0; s < nPin; s++) {
			for (int t = 0; t < nPin; t++) {

				int ixp = 0;
				for (int p = 0; p < nPout; p++) {
				for (int q = 0; q < nPout; q++) {
					dGlobalIntArray[ixp][ixOverlap][ixs] =
						dCoeff[ixp][i * nPin * nPin + ixs];

					ixp++;
				}
				}

				ixs++;
			}
			}
		}
*/
	}

	// Get SparseMatrix represntation of the OfflineMap
	SparseMatrix<double> & smatMap = mapRemap.GetSparseMatrix();

	// Compose the integration operator with the output map
	ixOverlap = 0;

	for (int ixFirst = 0; ixFirst < meshInput.faces.size(); ixFirst++) {

		// Output every 100 elements
		if (ixFirst % 1000 == 0) {
                  Announce("Element %i/%i", ixFirst,meshInput.faces.size());
		}

		// This Face
		const Face & faceFirst = meshInput.faces[ixFirst];

		// Number of overlapping Faces and triangles
		int nOverlapFaces = nAllOverlapFaces[ixFirst];

		// Put composed array into map
		for (int i = 0; i < nOverlapFaces; i++) {
			int ixSecond = meshOverlap.vecTargetFaceIx[ixOverlap + i];

			int ixs = 0;
			for (int s = 0; s < nPin; s++) {
			for (int t = 0; t < nPin; t++) {

				int ixFirstNode;
				if (fContinuousIn) {
					ixFirstNode = dataGLLNodesIn[s][t][ixFirst] - 1;
				} else {
					ixFirstNode = ixFirst * nPin * nPin + s * nPin + t;
				}

				int ixp = 0;
				for (int p = 0; p < nPout; p++) {
				for (int q = 0; q < nPout; q++) {

					int ixSecondNode;
					if (fContinuousOut) {
						ixSecondNode = dataGLLNodesOut[p][q][ixSecond] - 1;

						smatMap(ixSecondNode, ixFirstNode) +=
							dGlobalIntArray[ixp][ixOverlap + i][ixs]
							* dataGLLJacobianOut[p][q][ixSecond]
							/ dataNodalAreaOut[ixSecondNode];
							// dataIntAreaOut[ixSecondNode];

					} else {
						ixSecondNode = ixSecond * nPout * nPout + p * nPout + q;
					
						smatMap(ixSecondNode, ixFirstNode) +=
							dGlobalIntArray[ixp][ixOverlap + i][ixs];
							// dataIntAreaOut[ixSecondNode];
					}

					ixp++;
				}
				}

				ixs++;
			}
			}
		}

		// Increment the current overlap index
		ixOverlap += nOverlapFaces;
	}
}

///////////////////////////////////////////////////////////////////////////////

void LinearRemapGLLtoGLL_Integrated(
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
) {
	// Triangular quadrature rule
	TriangularQuadratureRule triquadrule(4);

	const DataArray2D<double> & dG = triquadrule.GetG();
	const DataArray1D<double> & dW = triquadrule.GetW();

	// Sample coefficients
	DataArray2D<double> dSampleCoeffIn(nPin, nPin);
	DataArray2D<double> dSampleCoeffOut(nPout, nPout);

	// Build the integration array for each element on meshOverlap
	DataArray3D<double> dGlobalIntArray(
		nPout * nPout,
		meshOverlap.faces.size(),
		nPin * nPin);

	DataArray2D<double> dConstantIntArray(
		meshOverlap.faces.size(),
		nPin * nPin);

	// Number of overlap Faces per source Face
	DataArray1D<int> nAllOverlapFaces(meshInput.faces.size());

	int ixOverlap = 0;

	for (int ixFirst = 0; ixFirst < meshInput.faces.size(); ixFirst++) {

		int ixOverlapTemp = ixOverlap;
		for (; ixOverlapTemp < meshOverlap.faces.size(); ixOverlapTemp++) {

			const Face & faceOverlap = meshOverlap.faces[ixOverlapTemp];

			if (meshOverlap.vecSourceFaceIx[ixOverlapTemp] != ixFirst) {
				break;
			}

			nAllOverlapFaces[ixFirst]++;
		}

		// Increment the current overlap index
		ixOverlap += nAllOverlapFaces[ixFirst];
	}

	// Mesh utilities
	MeshUtilitiesFuzzy meshutil;

	// Construct overlap
	ixOverlap = 0;

	int ixTotal = 0;

	std::set< std::pair<int, int> > setFound;

	for (int ixFirst = 0; ixFirst < meshInput.faces.size(); ixFirst++) {

		// Output every 100 elements
		if (ixFirst % 1000 == 0) {
                  Announce("Element %i/%i", ixFirst,meshInput.faces.size());
		}

		// Quantities from the First Mesh
		const Face & faceFirst = meshInput.faces[ixFirst];

		const NodeVector & nodesFirst = meshInput.nodes;

		// Number of overlapping Faces and triangles
		int nOverlapFaces = nAllOverlapFaces[ixFirst];

		// Loop through all Overlap Faces
		for (int i = 0; i < nOverlapFaces; i++) {

			// Quantities from the overlap Mesh
			const Face & faceOverlap = meshOverlap.faces[ixOverlap + i];

			const NodeVector & nodesOverlap = meshOverlap.nodes;

			int nOverlapTriangles = faceOverlap.edges.size() - 2;

			// Quantities from the Second Mesh
			int ixSecond = meshOverlap.vecTargetFaceIx[ixOverlap + i];

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

					// Get the nodal location of this point
					double dX[3];

					dX[0] = dG(k,0) * node0.x + dG(k,1) * node1.x + dG(k,2) * node2.x;
					dX[1] = dG(k,0) * node0.y + dG(k,1) * node1.y + dG(k,2) * node2.y;
					dX[2] = dG(k,0) * node0.z + dG(k,1) * node1.z + dG(k,2) * node2.z;

					double dMag =
						sqrt(dX[0] * dX[0] + dX[1] * dX[1] + dX[2] * dX[2]);

					dX[0] /= dMag;
					dX[1] /= dMag;
					dX[2] /= dMag;

					Node nodeQuadrature(dX[0], dX[1], dX[2]);

					// Find the components of this quadrature point in the basis
					// of the first Face.
					double dAlphaIn;
					double dBetaIn;

					ApplyInverseMap(
						faceFirst,
						nodesFirst,
						nodeQuadrature,
						dAlphaIn,
						dBetaIn);

					// Find the components of this quadrature point in the basis
					// of the second Face.
					double dAlphaOut;
					double dBetaOut;

					ApplyInverseMap(
						faceSecond,
						nodesSecond,
						nodeQuadrature,
						dAlphaOut,
						dBetaOut);

					// Check inverse map value
					if ((dAlphaIn < -InverseMapTolerance)      ||
						(dAlphaIn > 1.0 + InverseMapTolerance) ||
						(dBetaIn  < -InverseMapTolerance)      ||
						(dBetaIn  > 1.0 + InverseMapTolerance)
					) {
						printf("\n==== BEGIN DEBUGGING INFO ====\n");
						printf("WARNING (%s, Line %u) Inverse map out of range\n",
							__FILE__, __LINE__);
						printf("Source face ix %i, Target face ix %i, Overlap face ix %i\n",
							ixFirst, ixSecond, ixOverlap + i);
						printf("Face nodes:\n");
						for (int x = 0; x < faceFirst.edges.size(); x++) {
							nodesFirst[faceFirst[x]].Print("");
						}
						printf("Quadrature node:\n");
						nodeQuadrature.Print("");
						printf("Alpha, Beta: %1.15e %1.15e\n", dAlphaIn, dBetaIn);
						printf("==== END DEBUGGING INFO ====\n");
						//_EXCEPTION2("Inverse Map out of range (%1.5e %1.5e)",
						//	dAlphaIn, dBetaIn);
					}

					// Check inverse map value
					if ((dAlphaOut < -InverseMapTolerance)      ||
						(dAlphaOut > 1.0 + InverseMapTolerance) ||
						(dBetaOut  < -InverseMapTolerance)      ||
						(dBetaOut  > 1.0 + InverseMapTolerance)
					) {
						printf("\n==== BEGIN DEBUGGING INFO ====\n");
						printf("WARNING (%s, Line %u) Inverse map out of range\n",
							__FILE__, __LINE__);
						printf("Source face ix %i, Target face ix %i, Overlap face ix %i\n",
							ixFirst, ixSecond, ixOverlap + i);
						printf("Face nodes:\n");
						for (int x = 0; x < faceSecond.edges.size(); x++) {
							nodesSecond[faceSecond[x]].Print("");
						}
						printf("Quadrature node:\n");
						nodeQuadrature.Print("");
						printf("Alpha, Beta: %1.15e %1.15e\n", dAlphaOut, dBetaOut);
						printf("==== END DEBUGGING INFO ====\n");
						//_EXCEPTION2("Inverse Map out of range (%1.5e %1.5e)",
						//	dAlphaOut, dBetaOut);
					}

					// Sample the First finite element at this point
					SampleGLLFiniteElement(
						nMonotoneType,
						nPin,
						dAlphaIn,
						dBetaIn,
						dSampleCoeffIn);

					// Sample the Second finite element at this point
					SampleGLLFiniteElement(
						nMonotoneType,
						nPout,
						dAlphaOut,
						dBetaOut,
						dSampleCoeffOut);

					// Compute overlap integral
					int ixs = 0;
					for (int s = 0; s < nPin; s++) {
					for (int t = 0; t < nPin; t++) {

						dConstantIntArray[ixOverlap + i][ixs] +=
							  dSampleCoeffIn[s][t]
							* dW[k]
							* dTriArea
							/ meshOverlap.vecFaceArea[ixOverlap + i];

						int ixp = 0;
						for (int p = 0; p < nPout; p++) {
						for (int q = 0; q < nPout; q++) {

							dGlobalIntArray[ixp][ixOverlap + i][ixs] +=
								  dSampleCoeffIn[s][t]
								* dSampleCoeffOut[p][q]
								* dW[k]
								* dTriArea
								/ dataGLLJacobianOut[p][q][ixSecond];

							ixp++;
						}
						}
						ixs++;

					}
					}
				}
			}
		}

		// Force consistency and conservation of ConstantIntArray
		DataArray1D<double> dSourceArea(nPin * nPin);

		for (int s = 0; s < nPin; s++) {
		for (int t = 0; t < nPin; t++) {
			dSourceArea[s * nPin + t] = dataGLLJacobianIn[s][t][ixFirst];
		}
		}

		DataArray1D<double> dTargetArea(nOverlapFaces);
		for (int i = 0; i < nOverlapFaces; i++) {
			dTargetArea[i] = meshOverlap.vecFaceArea[ixOverlap + i];
		}

		DataArray2D<double> dCoeff(nOverlapFaces, nPin * nPin);

		for (int i = 0; i < nOverlapFaces; i++) {
		for (int s = 0; s < nPin; s++) {
		for (int t = 0; t < nPin; t++) {
			dCoeff[i][s * nPin + t] =
				dConstantIntArray[ixOverlap + i][s * nPin + t];
		}
		}
		}

		ForceConsistencyConservation3(
			dSourceArea,
			dTargetArea,
			dCoeff,
			(nMonotoneType != 0));

		for (int i = 0; i < nOverlapFaces; i++) {
			int ixSecondFace = meshOverlap.vecTargetFaceIx[ixOverlap + i];

			for (int s = 0; s < nPin; s++) {
			for (int t = 0; t < nPin; t++) {
				dConstantIntArray[ixOverlap + i][s * nPin + t] =
					dCoeff[i][s * nPin + t]
					* meshOverlap.vecFaceArea[ixOverlap + i]
					/ meshOutput.vecFaceArea[ixSecondFace];

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
		int ixSecond = meshOverlap.vecTargetFaceIx[i];

		vecReverseFaceIx[ixSecond].push_back(i);
	}

	// Force consistency and conservation on linear sub-map
	int ixDistBegin = 0;

	for (int ixSecond = 0; ixSecond < meshOutput.faces.size(); ixSecond++) {

		// Coefficients
		DataArray2D<double> dCoeff(
			nPout * nPout,
			vecReverseFaceIx[ixSecond].size() * nPin * nPin);

		for (int i = 0; i < vecReverseFaceIx[ixSecond].size(); i++) {

			int ixOverlap = vecReverseFaceIx[ixSecond][i];

			int ixs = 0;
			for (int s = 0; s < nPin; s++) {
			for (int t = 0; t < nPin; t++) {

				int ixp = 0;
				for (int p = 0; p < nPout; p++) {
				for (int q = 0; q < nPout; q++) {
					dCoeff[ixp][i * nPin * nPin + ixs] =
						dGlobalIntArray[ixp][ixOverlap][ixs];

					ixp++;
				}
				}

				ixs++;
			}
			}
		}

		// Source areas
		DataArray1D<double> dSourceArea(
			vecReverseFaceIx[ixSecond].size() * nPin * nPin);

		for (int i = 0; i < vecReverseFaceIx[ixSecond].size(); i++) {

			int ixOverlap = vecReverseFaceIx[ixSecond][i];

			int ixs = 0;
			for (int s = 0; s < nPin; s++) {
			for (int t = 0; t < nPin; t++) {
				dSourceArea[i * nPin * nPin + ixs] =
					dConstantIntArray[ixOverlap][ixs];

				ixs++;
			}
			}
		}

		// Target areas
		double dTotalTargetArea = 0.0;

		DataArray1D<double> dTargetArea(nPout * nPout);

		{
			int ixp = 0;
			for (int p = 0; p < nPout; p++) {
			for (int q = 0; q < nPout; q++) {

				dTargetArea[ixp] =
					dataGLLJacobianOut[p][q][ixSecond]
					/ meshOutput.vecFaceArea[ixSecond];

				ixp++;
			}
			}
		}

		// Force consistency and conservation
		ForceConsistencyConservation3(
			dSourceArea,
			dTargetArea,
			dCoeff,
			(nMonotoneType != 0));
/*
		// Check column sums (conservation)
		for (int i = 0; i < dCoeff.GetColumns(); i++) {
			double dColSum = 0.0;
			for (int j = 0; j < dCoeff.GetRows(); j++) {
				dColSum += dCoeff[j][i] * dTargetArea[j];
			}
			printf("Col %i: %1.15e\n", i, dColSum / dSourceArea[i]);
		}

		// Check row sums (consistency)
		for (int j = 0; j < dCoeff.GetRows(); j++) {
			double dRowSum = 0.0;
			for (int i = 0; i < dCoeff.GetColumns(); i++) {
				dRowSum += dCoeff[j][i];
			}
			printf("Row %i: %1.15e\n", j, dRowSum);
		}
		_EXCEPTION();
*/

		// Coefficients
		for (int i = 0; i < vecReverseFaceIx[ixSecond].size(); i++) {

			ixOverlap = vecReverseFaceIx[ixSecond][i];

			if ((ixOverlap < 0) || (ixOverlap > dGlobalIntArray.GetColumns())) {
				_EXCEPTION();
			}

			int ixs = 0;
			for (int s = 0; s < nPin; s++) {
			for (int t = 0; t < nPin; t++) {

				int ixp = 0;
				for (int p = 0; p < nPout; p++) {
				for (int q = 0; q < nPout; q++) {
					dGlobalIntArray[ixp][ixOverlap][ixs] =
						dCoeff[ixp][i * nPin * nPin + ixs];

					ixp++;
				}
				}

				ixs++;
			}
			}
		}

	}

	// Get SparseMatrix represntation of the OfflineMap
	SparseMatrix<double> & smatMap = mapRemap.GetSparseMatrix();

	// Compose the integration operator with the output map
	ixOverlap = 0;

	for (int ixFirst = 0; ixFirst < meshInput.faces.size(); ixFirst++) {

		// Output every 100 elements
		if (ixFirst % 1000 == 0) {
                  Announce("Element %i/%i", ixFirst,meshInput.faces.size());
		}

		// This Face
		const Face & faceFirst = meshInput.faces[ixFirst];

		// Number of overlapping Faces and triangles
		int nOverlapFaces = nAllOverlapFaces[ixFirst];

		// Put composed array into map
		for (int i = 0; i < nOverlapFaces; i++) {
			int ixSecond = meshOverlap.vecTargetFaceIx[ixOverlap + i];

			int ixs = 0;
			for (int s = 0; s < nPin; s++) {
			for (int t = 0; t < nPin; t++) {

				int ixFirstNode;
				if (fContinuousIn) {
					ixFirstNode = dataGLLNodesIn[s][t][ixFirst] - 1;
				} else {
					ixFirstNode = ixFirst * nPin * nPin + s * nPin + t;
				}

				int ixp = 0;
				for (int p = 0; p < nPout; p++) {
				for (int q = 0; q < nPout; q++) {

					int ixSecondNode;
					if (fContinuousOut) {
						ixSecondNode = dataGLLNodesOut[p][q][ixSecond] - 1;

						smatMap(ixSecondNode, ixFirstNode) +=
							dGlobalIntArray[ixp][ixOverlap + i][ixs]
							* dataGLLJacobianOut[p][q][ixSecond]
							/ dataNodalAreaOut[ixSecondNode];
							// dataIntAreaOut[ixSecondNode];

					} else {
						ixSecondNode = ixSecond * nPout * nPout + p * nPout + q;
					
						smatMap(ixSecondNode, ixFirstNode) +=
							dGlobalIntArray[ixp][ixOverlap + i][ixs];
							// dataIntAreaOut[ixSecondNode];
					}

					ixp++;
				}
				}

				ixs++;
			}
			}
		}

		// Increment the current overlap index
		ixOverlap += nOverlapFaces;
	}
}

///////////////////////////////////////////////////////////////////////////////

