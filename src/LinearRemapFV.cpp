///////////////////////////////////////////////////////////////////////////////
///
///	\file    LinearRemapFV.cpp
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

#include "Defines.h"
#include "LinearRemapFV.h"
#include "FiniteVolumeTools.h"
#include "GridElements.h"
#include "OfflineMap.h"
#include "FiniteElementTools.h"
#include "GaussQuadrature.h"
#include "GaussLobattoQuadrature.h"
#include "TriangularQuadrature.h"
#include "MeshUtilitiesFuzzy.h"
#include "OverlapMesh.h"

#include "Announce.h"
#include "MathHelper.h"
#include "kdtree.h"

#include <cstring>
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

	// Unconstrained least squares solve (used for pseudoinverse)
	int dgelss_(
		int * m,
		int * n,
		int * nrhs,
		double * a,
		int * lda,
		double * b,
		int * ldb,
		double * s,
		double * rcond,
		int * rank,
		double * work,
		int * lwork,
		int * info);
};

///////////////////////////////////////////////////////////////////////////////

void LinearRemapFVtoFV_np1(
	const Mesh & meshInput,
	const Mesh & meshOutput,
	const Mesh & meshOverlap,
	OfflineMap & mapRemap
) {
	// Verify ReverseNodeArray has been calculated
	if (meshInput.revnodearray.size() == 0) {
		_EXCEPTIONT("ReverseNodeArray has not been calculated for meshInput");
	}

	// Get SparseMatrix represntation of the OfflineMap
	SparseMatrix<double> & smatMap = mapRemap.GetSparseMatrix();

	// Current overlap face
	int ixOverlap = 0;

	// Loop through all faces on meshInput
	for (int ixFirst = 0; ixFirst < meshInput.faces.size(); ixFirst++) {

		// Output every 1000 elements
		if (ixFirst % 1000 == 0) {
			Announce("Element %i/%i", ixFirst, meshInput.faces.size());
		}

		// This Face
		const Face & faceFirst = meshInput.faces[ixFirst];

		// Find the set of Faces that overlap faceFirst
		int ixOverlapBegin = ixOverlap;
		int ixOverlapEnd = ixOverlapBegin;
	
		for (; ixOverlapEnd < meshOverlap.faces.size(); ixOverlapEnd++) {
			if (meshOverlap.vecSourceFaceIx[ixOverlapEnd] != ixFirst) {
				break;
			}
		}

		int nOverlapFaces = ixOverlapEnd - ixOverlapBegin;

		// Put composed array into map
		for (int j = 0; j < nOverlapFaces; j++) {
			int ixFirstFace = meshOverlap.vecSourceFaceIx[ixOverlap + j];
			int ixSecondFace = meshOverlap.vecTargetFaceIx[ixOverlap + j];

			smatMap(ixSecondFace, ixFirstFace) +=
				meshOverlap.vecFaceArea[ixOverlap + j]
				/ meshOutput.vecFaceArea[ixSecondFace];

			if (smatMap(ixSecondFace, ixFirstFace) > 10.0) {
				printf("%i %i %i\n", ixFirstFace, ixSecondFace, ixOverlap+j);
				printf("Input:\n");
				for (int i = 0; i < meshInput.faces[ixFirstFace].edges.size(); i++) {
					const Node & node = meshInput.nodes[ meshInput.faces[ixFirstFace][i] ];
					printf("%i,%1.15e,%1.15e;\n",
						i, atan2(node.y, node.x), asin(node.z));
				}
				printf("Output:\n");
				for (int i = 0; i < meshOutput.faces[ixSecondFace].edges.size(); i++) {
					const Node & node = meshOutput.nodes[ meshOutput.faces[ixSecondFace][i] ];
					printf("%i,%1.15e,%1.15e;\n",
						i, atan2(node.y, node.x), asin(node.z));
				}
				printf("Overlap:\n");
				for (int i = 0; i < meshOverlap.faces[ixOverlap + j].edges.size(); i++) {
					const Node & node = meshOverlap.nodes[ meshOverlap.faces[ixOverlap + j][i] ];
					printf("%i,%1.15e,%1.15e;\n",
						i, atan2(node.y, node.x), asin(node.z));
				}


				printf("%1.15e\n", meshInput.vecFaceArea[ixFirstFace]);
				printf("%1.15e\n", meshOverlap.vecFaceArea[ixOverlap + j]);
				printf("%1.15e\n", meshOutput.vecFaceArea[ixSecondFace]);
				_EXCEPTIONT("Anomalous map weight detected");
			}
		}

		// Increment the current overlap index
		ixOverlap += nOverlapFaces;

	}
}

///////////////////////////////////////////////////////////////////////////////

void LinearRemapFVtoFV(
	const Mesh & meshInput,
	const Mesh & meshOutput,
	const Mesh & meshOverlap,
	int nOrder,
	OfflineMap & mapRemap
) {
	// Use streamlined helper function for first order
	if (nOrder == 1) {
		return LinearRemapFVtoFV_np1(
			meshInput, meshOutput, meshOverlap, mapRemap);
	}

	// Order of triangular quadrature rule
	const int TriQuadRuleOrder = 4;

	// Verify ReverseNodeArray has been calculated
	if (meshInput.revnodearray.size() == 0) {
		_EXCEPTIONT("ReverseNodeArray has not been calculated for meshInput");
	}

	// Triangular quadrature rule
	TriangularQuadratureRule triquadrule(TriQuadRuleOrder);

	// Get SparseMatrix represntation of the OfflineMap
	SparseMatrix<double> & smatMap = mapRemap.GetSparseMatrix();

	// Number of elements needed
#ifdef RECTANGULAR_TRUNCATION
	const int nCoefficients = nOrder * nOrder;
#endif
#ifdef TRIANGULAR_TRUNCATION 
	const int nCoefficients = nOrder * (nOrder + 1) / 2;
#endif

//#pragma message "This should be a command-line parameter"
	// Number of faces you need
	const int nRequiredFaceSetSize = nCoefficients;

	// Fit weight exponent
	const int nFitWeightsExponent = nOrder + 2;

	// Announcemnets
	Announce("Triangular quadrature rule order %i", TriQuadRuleOrder);
	Announce("Number of coefficients: %i", nCoefficients);
	Announce("Required adjacency set size: %i", nRequiredFaceSetSize);
	Announce("Fit weights exponent: %i", nFitWeightsExponent);

	// Current overlap face
	int ixOverlap = 0;

	// Arrays used in analysis
	DataArray2D<double> dIntArray;
	DataArray1D<double> dConstraint(nCoefficients);

	// Loop through all faces on meshInput
	for (int ixFirst = 0; ixFirst < meshInput.faces.size(); ixFirst++) {

		// Output every 100 elements
		if (ixFirst % 1000 == 0) {
			Announce("Element %i/%i", ixFirst, meshInput.faces.size());
		}

		// This Face
		const Face & faceFirst = meshInput.faces[ixFirst];

		// Find the set of Faces that overlap faceFirst
		int ixOverlapBegin = ixOverlap;
		int ixOverlapEnd = ixOverlapBegin;
	
		for (; ixOverlapEnd < meshOverlap.faces.size(); ixOverlapEnd++) {
			if (meshOverlap.vecSourceFaceIx[ixOverlapEnd] != ixFirst) {
				break;
			}
		}

		int nOverlapFaces = ixOverlapEnd - ixOverlapBegin;

		// Build integration array, which maps polynomial coefficients to
		// area integrals.
		BuildIntegrationArray(
			meshInput,
			meshOverlap,
			triquadrule,
			ixFirst,
			ixOverlapBegin,
			ixOverlapEnd,
			nOrder,
			dIntArray);

		// Set of Faces to use in building the reconstruction and associated
		// distance metric.
		AdjacentFaceVector vecAdjFaces;

//#ifdef RECTANGULAR_TRUNCATION
//		GetAdjacentFaceVectorByNode(
//#endif
//#ifdef TRIANGULAR_TRUNCATION
		GetAdjacentFaceVectorByEdge(
//#endif
			meshInput,
			ixFirst,
			nRequiredFaceSetSize,
			vecAdjFaces);

		// Number of adjacent Faces
		int nAdjFaces = vecAdjFaces.size();

		// Determine the conservative constraint equation
		double dFirstArea = meshInput.vecFaceArea[ixFirst];

		dConstraint.Zero();
		for (int p = 0; p < nCoefficients; p++) {
			for (int j = 0; j < nOverlapFaces; j++) {
				dConstraint[p] += dIntArray(p,j);
			}
			dConstraint[p] /= dFirstArea;
		}

		// Build the fit array from the integration operator
		DataArray2D<double> dFitArray;
		DataArray1D<double> dFitWeights;
		DataArray2D<double> dFitArrayPlus;

		BuildFitArray(
			meshInput,
			triquadrule,
			ixFirst,
			vecAdjFaces,
			nOrder,
			nFitWeightsExponent,
			dConstraint,
			dFitArray,
			dFitWeights
		);

		// Compute the inverse fit array
		bool fSuccess =
			InvertFitArray_Corrected(
				dConstraint,
				dFitArray,
				dFitWeights,
				dFitArrayPlus
			);

		// Build the composition, which maps average values in adjacent cells
		// to the integrated values of the reconstruction in overlap faces.
		DataArray2D<double> dComposedArray(nAdjFaces, nOverlapFaces);
		if (fSuccess) {

			// Multiply integration array and inverse fit array
			for (int i = 0; i < nAdjFaces; i++) {
			for (int j = 0; j < nOverlapFaces; j++) {
			for (int k = 0; k < nCoefficients; k++) {
				dComposedArray(i,j) += dIntArray(k,j) * dFitArrayPlus(i,k);
			}
			}
			}

		// Unable to invert fit array, drop to 1st order.  In this case
		// dFitArrayPlus(0,0) = 1 and all other entries are zero.
		} else {
			dComposedArray.Zero();
			for (int j = 0; j < nOverlapFaces; j++) {
				dComposedArray(0,j) += dIntArray(0,j);
			}
		}


/*
		for (int j = 0; j < nOverlapFaces; j++) {
			dComposedArray(0,j) = meshOverlap.vecFaceArea[ixOverlap + j];
		}

		for (int i = 0; i < nAdjFaces; i++) {
		for (int j = 0; j < nOverlapFaces; j++) {
		for (int k = 1; k < nCoefficients; k++) {
			double dOverlapArea =
				meshOverlap.vecFaceArea[ixOverlap + j];

			dComposedArray(i,j) +=
				(dIntArray(k,j) - dConstraint[k] * dOverlapArea)
					* dFitArrayPlus(i,k);
		}
		}
		}
*/

		// Put composed array into map
		for (int i = 0; i < vecAdjFaces.size(); i++) {
		for (int j = 0; j < nOverlapFaces; j++) {
			int ixFirstFace = vecAdjFaces[i].first;
			int ixSecondFace = meshOverlap.vecTargetFaceIx[ixOverlap + j];

			smatMap(ixSecondFace, ixFirstFace) +=
				dComposedArray(i,j)
				/ meshOutput.vecFaceArea[ixSecondFace];
		}
		}

		// Increment the current overlap index
		ixOverlap += nOverlapFaces;
	}
}

///////////////////////////////////////////////////////////////////////////////

void LinearRemapFVtoFVInvDist(
	const Mesh & meshInput,
	const Mesh & meshOutput,
	const Mesh & meshOverlap,
	OfflineMap & mapRemap
) {
	// Verify ReverseNodeArray has been calculated
	if (meshInput.edgemap.size() == 0) {
		_EXCEPTIONT("EdgeMap has not been calculated for meshInput");
	}

	// Order of triangular quadrature rule
	const int TriQuadRuleOrder = 4;

	// Triangular quadrature rule
	TriangularQuadratureRule triquadrule(TriQuadRuleOrder);

	const DataArray2D<double> & dG = triquadrule.GetG();
	const DataArray1D<double> & dW = triquadrule.GetW();

	// Get SparseMatrix represntation of the OfflineMap
	SparseMatrix<double> & smatMap = mapRemap.GetSparseMatrix();

	// kd-tree for nearest neighbor search
    kdtree * kdTarget = kd_create(3);

	// Vectors used in determining contributions from each Face
	std::vector<int> vecContributingFaceIxs;
	std::vector<double> vecContributingFaceWeights;

	// Vector of centers of the source mesh
	for (int i = 0; i < meshInput.faces.size(); i++){

		const Face & face = meshInput.faces[i];

		Node nodeCentroid = GetFaceCentroid(face, meshInput.nodes);

		kd_insert3(
			kdTarget,
			nodeCentroid.x,
			nodeCentroid.y,
			nodeCentroid.z,
			(void*)(&(meshInput.faces[i])));
	}

	// Overlap face index
	int ixOverlap = 0;

	// Loop through all source faces
	for (int ixFirst = 0; ixFirst < meshInput.faces.size(); ixFirst++) {

		// Output every 1000 overlap elements
		if (ixFirst % 1000 == 0) {
			Announce("Element %i/%i", ixFirst, meshInput.faces.size());
		}

		// This Face
		const Face & faceFirst = meshInput.faces[ixFirst];

		// Find the set of Faces that overlap faceFirst
		int ixOverlapBegin = ixOverlap;
		int ixOverlapEnd = ixOverlapBegin;
	
		for (; ixOverlapEnd < meshOverlap.faces.size(); ixOverlapEnd++) {
			if (meshOverlap.vecSourceFaceIx[ixOverlapEnd] != ixFirst) {
				break;
			}
		}

		int nOverlapFaces = ixOverlapEnd - ixOverlapBegin;

		// Loop through all overlap faces associated with this source face
		for (int j = 0; j < nOverlapFaces; j++) {
			int iTargetFace = meshOverlap.vecTargetFaceIx[ixOverlap + j];

			const Face & faceOverlap = meshOverlap.faces[ixOverlap + j];

			int nSubTriangles = faceOverlap.edges.size() - 2;

			// Jacobian at each quadrature point
			DataArray2D<double> dQuadPtWeight(nSubTriangles, dW.GetRows());
			DataArray2D<Node> dQuadPtNodes(nSubTriangles, dW.GetRows());

			// Compute quadrature area of each overlap face
			double dQuadratureArea = 0.0;

			for (int k = 0; k < nSubTriangles; k++) {
				for (int p = 0; p < dW.GetRows(); p++) {

					dQuadPtWeight(k,p) =
						CalculateSphericalTriangleJacobianBarycentric(
							meshOverlap.nodes[faceOverlap[0]],
							meshOverlap.nodes[faceOverlap[k+1]],
							meshOverlap.nodes[faceOverlap[k+2]],
							dG(p,0), dG(p,1),
							&(dQuadPtNodes(k,p)));

					dQuadPtWeight(k,p) *= dW[p];

					dQuadratureArea += dQuadPtWeight(k,p);
				}
			}

			//printf("%1.15e %1.15e %1.15e %1.15e\n", meshOverlap.vecFaceArea[ixOverlap + j], dQuadratureArea, CalculateFaceAreaKarneysMethod(faceOverlap, meshOverlap.nodes), dQuadratureArea / meshOverlap.vecFaceArea[ixOverlap + j]);

			// Loop through all sub-triangles of this overlap Face
			for (int k = 0; k < nSubTriangles; k++) {

				// Loop through all quadrature nodes on this overlap Face
				for (int p = 0; p < dW.GetRows(); p++) {

					// Get quadrature node and pointwise Jacobian
					const Node & nodeQ = dQuadPtNodes(k,p);

					// Find nearest source mesh face and add its contribution
					// to the inverse distance.
					kdres * kdresTarget =
						kd_nearest3(
							kdTarget,
							nodeQ.x,
							nodeQ.y,
							nodeQ.z);

					Face * pFace = (Face *)(kd_res_item_data(kdresTarget));

					int iNearestFace = pFace - &(meshInput.faces[0]);

					const Face & faceCurrent = meshInput.faces[iNearestFace];

					// TODO: Precompute face centroids on input mesh
					Node nodeX1 =
						GetFaceCentroid(
							meshInput.faces[iNearestFace],
							meshInput.nodes);

					Node nodeX1minusQ = nodeX1;
					nodeX1minusQ.x -= nodeQ.x;
					nodeX1minusQ.y -= nodeQ.y;
					nodeX1minusQ.z -= nodeQ.z;

					vecContributingFaceIxs.clear();
					vecContributingFaceWeights.clear();

					// TODO: Switch to using great circle distance rather than chord dist
					// TODO: Be careful about nodeCenter1 being zero
					double dInvDist1 = 1.0 / nodeX1minusQ.Magnitude();
					vecContributingFaceIxs.push_back(iNearestFace);
					vecContributingFaceWeights.push_back(dInvDist1);

					// Push neighboring face if it is on the correct side
					for (int i = 0; i < faceCurrent.edges.size(); i++) {

						const FacePair & facepair =
								meshInput.edgemap.find(faceCurrent.edges[i])->second;

						if (iNearestFace == facepair[0]){
							Node nodeX2 =
								GetFaceCentroid(
									meshInput.faces[facepair[1]],
									meshInput.nodes);

							Node nodeX1minusX2 = nodeX1;
							nodeX1minusX2.x -= nodeX2.x;
							nodeX1minusX2.y -= nodeX2.y;
							nodeX1minusX2.z -= nodeX2.z;

							if (DotProduct(nodeX1minusX2, nodeX1minusQ) > 0.0) {
								Node nodeX2minusQ = nodeX2;
								nodeX2minusQ.x -= nodeQ.x;
								nodeX2minusQ.y -= nodeQ.y;
								nodeX2minusQ.z -= nodeQ.z;

								double dInvDist2 = 1.0 / nodeX2minusQ.Magnitude();
								vecContributingFaceIxs.push_back(facepair[1]);
								vecContributingFaceWeights.push_back(dInvDist2);
							}

						} else if (iNearestFace == facepair[1]) {
							Node nodeX2 =
								GetFaceCentroid(
									meshInput.faces[facepair[0]],
									meshInput.nodes);

							Node nodeX1minusX2 = nodeX1;
							nodeX1minusX2.x -= nodeX2.x;
							nodeX1minusX2.y -= nodeX2.y;
							nodeX1minusX2.z -= nodeX2.z;

							if (DotProduct(nodeX1minusX2, nodeX1minusQ) > 0.0){
								Node nodeX2minusQ = nodeX2;
								nodeX2minusQ.x -= nodeQ.x;
								nodeX2minusQ.y -= nodeQ.y;
								nodeX2minusQ.z -= nodeQ.z;

								double dInvDist2 = 1.0 / nodeX2minusQ.Magnitude();
								vecContributingFaceIxs.push_back(facepair[0]);
								vecContributingFaceWeights.push_back(dInvDist2);
							}

						} else {
							_EXCEPTIONT("Logic error");
						}
					}

					// Total contributions
					double dInvWeightSum = 0.0;
					for (int i = 0; i < vecContributingFaceWeights.size(); i++){
						dInvWeightSum += vecContributingFaceWeights[i];
					}

					// Contribution of this quadrature point to the map
					for (int i = 0; i < vecContributingFaceIxs.size(); i++){
						smatMap(iTargetFace, vecContributingFaceIxs[i]) +=
							vecContributingFaceWeights[i]
							/ dInvWeightSum
							* dQuadPtWeight(k,p)
							* meshOverlap.vecFaceArea[ixOverlap + j]
							/ dQuadratureArea
							/ meshOutput.vecFaceArea[iTargetFace];
					}
				}
			}
		}

		// Increment the current overlap index
		ixOverlap += nOverlapFaces;
	}
}

///////////////////////////////////////////////////////////////////////////////

void ForceIntArrayConsistencyConservation(
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
		dRHS[ix] = dCoeff(i,j);
		ix++;
	}
	}

	// Consistency
	ix = 0;
	for (int i = 0; i < dCoeff.GetRows(); i++) {
		for (int j = 0; j < dCoeff.GetColumns(); j++) {
			dC(i * dCoeff.GetColumns() + j,ix) = 1.0;
		}
		dRHS[nCoeff + ix] = 1.0;
		ix++;
	}

	// Conservation
	for (int j = 0; j < dCoeff.GetColumns()-1; j++) {
		for (int i = 0; i < dCoeff.GetRows(); i++) {
			dC(i * dCoeff.GetColumns() + j,ix) = vecTargetArea[i];
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
		dCCt(i,i) = static_cast<double>(dCoeff.GetColumns());
		for (int j = 0; j < dCoeff.GetColumns()-1; j++) {
			dCCt(i,dCoeff.GetRows() + j) = vecTargetArea[i];
			dCCt(dCoeff.GetRows() + j,i) = vecTargetArea[i];
		}
	}

	for (int i = 0; i < dCoeff.GetColumns()-1; i++) {
		int ix = dCoeff.GetRows() + i;
		dCCt(ix,ix) = dP;
	}
/*
	for (int i = 0; i < nCond; i++) {
	for (int j = 0; j < nCond; j++) {
		for (int k = 0; k < nCoeff; k++) {
			dCCt(i,j) += dC(k,i) * dC(k,j);
		}
	}
	}
*/
/*
	FILE * fp = fopen("cct.dat", "w");
	for (int i = 0; i < nCond; i++) {
		for (int j = 0; j < nCond; j++) {
			fprintf(fp, "%1.15e\t", dCCt(i,j));
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
		&(dC(0,0)),
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
		&(dCCt(0,0)),
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
		&(dCCt(0,0)),
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
		&(dC(0,0)),
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
		dCoeff(i,j) = dRHS[ix];
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
			dMonoCoeff(i,j) =
				vecSourceArea[j]
				/ dTotalJacobian;
		}
		}

		// Compute scaling factor
		double dA = 0.0;
		for (int i = 0; i < dCoeff.GetRows(); i++) {
		for (int j = 0; j < dCoeff.GetColumns(); j++) {
			if (dCoeff(i,j) < 0.0) {
				double dNewA =
					- dCoeff(i,j) / fabs(dMonoCoeff(i,j) - dCoeff(i,j));

				if (dNewA > dA) {
					dA = dNewA;
				}
			}
		}
		}

		for (int i = 0; i < dCoeff.GetRows(); i++) {
		for (int j = 0; j < dCoeff.GetColumns(); j++) {
			dCoeff(i,j) = (1.0 - dA) * dCoeff(i,j) + dA * dMonoCoeff(i,j);
		}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

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
) {
	// Order of triangular quadrature rule
	const int TriQuadRuleOrder = 8;

	// Verify ReverseNodeArray has been calculated
	if (meshInput.revnodearray.size() == 0) {
		_EXCEPTIONT("ReverseNodeArray has not been calculated for meshInput");
	}
	if (meshInput.edgemap.size() == 0) {
		_EXCEPTIONT("EdgeMap has not been calculated for meshInput");
	}

	// Get SparseMatrix represntation of the OfflineMap
	SparseMatrix<double> & smatMap = mapRemap.GetSparseMatrix();

	// Fit weight exponent
	int nFitWeightsExponent = nOrder + 2;

	// Order of the finite element method
	int nP = dataGLLNodes.GetRows();

	// Mesh utilities
	MeshUtilitiesFuzzy meshutil;

	// GLL nodes
	DataArray1D<double> dG;
	DataArray1D<double> dW;

	GaussLobattoQuadrature::GetPoints(nP, 0.0, 1.0, dG, dW);

	// Triangular quadrature rule
	TriangularQuadratureRule triquadrule(TriQuadRuleOrder);

	// Number of elements needed
#ifdef RECTANGULAR_TRUNCATION
	int nCoefficients = nOrder * nOrder;
#endif
#ifdef TRIANGULAR_TRUNCATION 
	int nCoefficients = nOrder * (nOrder + 1) / 2;
#endif

	int nRequiredFaceSetSize = nCoefficients;

	// Set of found nodes
	std::set<int> setFoundNodes;

	// Loop through all faces on meshInput
	int ixOverlap = 0;

	for (int ixFirst = 0; ixFirst < meshInput.faces.size(); ixFirst++) {

		// Output every 100 elements
		if (ixFirst % 1000 == 0) {
                        Announce("Element %i/%i", ixFirst, meshInput.faces.size());
		}

		// This Face
		const Face & faceFirst = meshInput.faces[ixFirst];

		// Area of the First Face
		double dFirstArea = meshInput.vecFaceArea[ixFirst];

		// Coordinate axes
		Node nodeRef = GetFaceCentroid(faceFirst, meshInput.nodes);

		Node nodeA1 = meshInput.nodes[faceFirst[1]] - nodeRef;
		Node nodeA2 = meshInput.nodes[faceFirst[2]] - nodeRef;

		Node nodeC = CrossProduct(nodeA1, nodeA2);

		// Fit matrix
		DataArray2D<double> dFit(3,3);

		dFit(0,0) = nodeA1.x; dFit(0,1) = nodeA1.y; dFit(0,2) = nodeA1.z;
		dFit(1,0) = nodeA2.x; dFit(1,1) = nodeA2.y; dFit(1,2) = nodeA2.z;
		dFit(2,0) = nodeC.x;  dFit(2,1) = nodeC.y;  dFit(2,2) = nodeC.z;

		// Set of Faces to use in building the reconstruction and associated
		// distance metric.
		AdjacentFaceVector vecAdjFaces;

//#ifdef RECTANGULAR_TRUNCATION
//		GetAdjacentFaceVectorByNode(
//#endif
//#ifdef TRIANGULAR_TRUNCATION
		GetAdjacentFaceVectorByEdge(
//#endif
			meshInput,
			ixFirst,
			nRequiredFaceSetSize,
			vecAdjFaces);

		// Number of adjacent Faces
		int nAdjFaces = vecAdjFaces.size();

		// Blank constraint
		DataArray1D<double> dConstraint;

		// Least squares arrays
		DataArray2D<double> dFitArray;
		DataArray1D<double> dFitWeights;
		DataArray2D<double> dFitArrayPlus;

		BuildFitArray(
			meshInput,
			triquadrule,
			ixFirst,
			vecAdjFaces,
			nOrder,
			nFitWeightsExponent,
			dConstraint,
			dFitArray,
			dFitWeights
		);

		// Compute the inverse fit array
		bool fSuccess =
			InvertFitArray_Corrected(
				dConstraint,
				dFitArray,
				dFitWeights,
				dFitArrayPlus
			);

		if (!fSuccess) {
			dFitArrayPlus.Zero();
			dFitArrayPlus(0,0) = 1.0;
		}

		// Number of overlapping Faces
		int nOverlapFaces = 0;
		int ixOverlapTemp = ixOverlap;
		for (; ixOverlapTemp < meshOverlap.faces.size(); ixOverlapTemp++) {

			const Face & faceOverlap = meshOverlap.faces[ixOverlapTemp];

			if (meshOverlap.vecSourceFaceIx[ixOverlapTemp] != ixFirst) {
				break;
			}

			nOverlapFaces++;
		}

		// Loop through all Overlap Faces
		for (int i = 0; i < nOverlapFaces; i++) {

			// Quantities from the overlap Mesh
			const Face & faceOverlap = meshOverlap.faces[ixOverlap + i];

			// Quantities from the Second Mesh
			int ixSecondFace = meshOverlap.vecTargetFaceIx[ixOverlap + i];

			const NodeVector & nodesSecond = meshOutput.nodes;

			const Face & faceSecond = meshOutput.faces[ixSecondFace];

			for (int s = 0; s < nP; s++) {
			for (int t = 0; t < nP; t++) {
				
				// Determine if this Node is in faceFirst
				Node node;
				Node dDx1G;
				Node dDx2G;

				ApplyLocalMap(
					faceSecond,
					meshOutput.nodes,
					dG[s],
					dG[t],
					node,
					dDx1G,
					dDx2G);

				Face::NodeLocation loc;
				int ixLocation;

				meshutil.ContainsNode(
					faceFirst,
					meshInput.nodes,
					node,
					loc,
					ixLocation);

				if (loc == Face::NodeLocation_Exterior) {
					continue;
				}

				// Second node index
				int ixSecondNode;

				if (fContinuous) {
					ixSecondNode = dataGLLNodes(t,s,ixSecondFace) - 1;
				} else {
					ixSecondNode = ixSecondFace * nP * nP + s * nP + t;
				}

				// Avoid duplicates
				if (setFoundNodes.find(ixSecondNode) != setFoundNodes.end()) {
					continue;
				}

				setFoundNodes.insert(ixSecondNode);

				// Find the coefficients for this point
				double dX[3];
				dX[0] = node.x - nodeRef.x;
				dX[1] = node.y - nodeRef.y;
				dX[2] = node.z - nodeRef.z;
/*
				if (ixSecondNode < 10) {

					double dLon = atan2(node.y, node.x);
					if (dLon < 0.0) {
						dLon += 2.0 * M_PI;
					}
					double dLat = asin(node.z);

					dLon *= 180.0 / M_PI;
					dLat *= 180.0 / M_PI;

					printf("%i %1.15e %1.15e\n", ixSecondNode, dLon, dLat);
				}
*/
				int n = 3;
				int nrhs = 1;
				int lda = 3;
				int ipiv[3];
				int ldb = 3;
				int info;

				DataArray2D<double> dFitTemp;
				dFitTemp = dFit;
				dgesv_(
					&n, &nrhs, &(dFitTemp(0,0)), &lda, ipiv, dX, &ldb, &info);

				// Sample the reconstruction at this point
				int ixp = 0;

#ifdef RECTANGULAR_TRUNCATION
				for (int p = 0; p < nOrder; p++) {
				for (int q = 0; q < nOrder; q++) {
#endif
#ifdef TRIANGULAR_TRUNCATION 
				for (int p = 0; p < nOrder; p++) {
				for (int q = 0; q < nOrder - p; q++) {
#endif

					for (int n = 0; n < vecAdjFaces.size(); n++) {
						int ixAdjFace = vecAdjFaces[n].first;

						smatMap(ixSecondNode, ixAdjFace) +=
							  IPow(dX[0], p)
							* IPow(dX[1], q)
							* dFitArrayPlus(n,ixp);
					}

					ixp++;
				}
				}
			}
			}
		}

		// Increment the current overlap index
		ixOverlap += nOverlapFaces;

	}
}

///////////////////////////////////////////////////////////////////////////////

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
) {
	// Order of triangular quadrature rule
	const int TriQuadRuleOrder = 4;

	// Verify ReverseNodeArray has been calculated
	if (meshInput.revnodearray.size() == 0) {
		_EXCEPTIONT("ReverseNodeArray has not been calculated for meshInput");
	}

	// Triangular quadrature rule
	TriangularQuadratureRule triquadrule(TriQuadRuleOrder);

	// Fit weight exponent
	int nFitWeightsExponent = nOrder + 2;

	// Order of the finite element method
	int nP = dataGLLNodes.GetRows();

	// Gauss-Lobatto quadrature nodes and weights
	DataArray1D<double> dG;
	DataArray1D<double> dW;

	GaussLobattoQuadrature::GetPoints(nP, 0.0, 1.0, dG, dW);

	// Get SparseMatrix represntation of the OfflineMap
	SparseMatrix<double> & smatMap = mapRemap.GetSparseMatrix();

	// Number of elements needed
#ifdef RECTANGULAR_TRUNCATION
	int nCoefficients = nOrder * nOrder;
#endif
#ifdef TRIANGULAR_TRUNCATION 
	int nCoefficients = nOrder * (nOrder + 1) / 2;
#endif

//#pragma message "This should be a command-line parameter"
	int nRequiredFaceSetSize = nCoefficients;

	// Accumulated weight vector
	DataArray1D<double> dAccumW(nP+1);
	dAccumW[0] = 0.0;
	for (int i = 1; i < nP+1; i++) {
		dAccumW[i] = dAccumW[i-1] + dW[i-1];
	}
	if (fabs(dAccumW[dAccumW.GetRows()-1] - 1.0) > 1.0e-14) {
		_EXCEPTIONT("Logic error in accumulated weight");
	}

	// Create sub-element mesh and redistribution map
	Announce("Generating sub-element mesh");
	Mesh meshTargetSubElement;

	DataArray1D<double> dFiniteVolumeArea(nP * nP);
	DataArray1D<double> dQuadratureArea(nP * nP);
	std::vector< DataArray2D<double> > dRedistributionMaps;
	dRedistributionMaps.resize(meshOutput.faces.size());

	for (int ixSecond = 0; ixSecond < meshOutput.faces.size(); ixSecond++) {

		const Face & faceSecond = meshOutput.faces[ixSecond];

		const Node & nodeOutput0 = meshOutput.nodes[faceSecond[0]];
		const Node & nodeOutput1 = meshOutput.nodes[faceSecond[1]];
		const Node & nodeOutput2 = meshOutput.nodes[faceSecond[2]];
		const Node & nodeOutput3 = meshOutput.nodes[faceSecond[3]];

		for (int q = 0; q < nP; q++) {
		for (int p = 0; p < nP; p++) {

			Node node0 =
				InterpolateQuadrilateralNode(
					nodeOutput0, nodeOutput1, nodeOutput2, nodeOutput3,
					dAccumW[p], dAccumW[q]);

			Node node1 =
				InterpolateQuadrilateralNode(
					nodeOutput0, nodeOutput1, nodeOutput2, nodeOutput3,
					dAccumW[p+1], dAccumW[q]);

			Node node2 =
				InterpolateQuadrilateralNode(
					nodeOutput0, nodeOutput1, nodeOutput2, nodeOutput3,
					dAccumW[p+1], dAccumW[q+1]);

			Node node3 =
				InterpolateQuadrilateralNode(
					nodeOutput0, nodeOutput1, nodeOutput2, nodeOutput3,
					dAccumW[p], dAccumW[q+1]);

			int nNodeStart = meshTargetSubElement.nodes.size();
			meshTargetSubElement.nodes.push_back(node0);
			meshTargetSubElement.nodes.push_back(node1);
			meshTargetSubElement.nodes.push_back(node2);
			meshTargetSubElement.nodes.push_back(node3);

			Face faceNew(4);
			faceNew.SetNode(0, nNodeStart);
			faceNew.SetNode(1, nNodeStart+1);
			faceNew.SetNode(2, nNodeStart+2);
			faceNew.SetNode(3, nNodeStart+3);

			meshTargetSubElement.faces.push_back(faceNew);

			dFiniteVolumeArea[q * nP + p] =
				CalculateFaceArea(
					faceNew, meshTargetSubElement.nodes);

			dQuadratureArea[q * nP + p] =
				dataGLLJacobian(q,p,ixSecond);
		}
		}

		dRedistributionMaps[ixSecond].Allocate(nP * nP, nP * nP);
		for (int i = 0; i < nP * nP; i++) {
			dRedistributionMaps[ixSecond](i,i) = 1.0;
		}

		if (!fNoConservation) {
			ForceIntArrayConsistencyConservation(
				dFiniteVolumeArea,
				dQuadratureArea,
				dRedistributionMaps[ixSecond],
				(nMonotoneType != 0));
		}
/*
		double dSumQuadArea = 0.0;
		double dSumFVArea = 0.0;
		for (int i = 0; i < nP * nP; i++) {
			dSumQuadArea += dQuadratureArea[i];
			dSumFVArea += dFiniteVolumeArea[i];
		}
		if (fabs(dSumQuadArea - dSumFVArea) > 1.0e-14) {
			printf("%1.15e\n", dSumQuadArea - dSumFVArea);
			_EXCEPTION();
		}
*/
		for (int i = 0; i < nP * nP; i++) {
		for (int j = 0; j < nP * nP; j++) {
			dRedistributionMaps[ixSecond](i,j) *=
				dQuadratureArea[i] / dFiniteVolumeArea[j];
		}
		}
/*
		for (int i = 0; i < nP * nP; i++) {
			double dSum = 0.0;
			for (int j = 0; j < nP * nP; j++) {
				dSum += dRedistributionMaps(ixSecond,i,j) * dFiniteVolumeArea[j];
			}
			printf("%i %1.15e %1.15e\n", i, dSum, dQuadratureArea[i]);
		}
		_EXCEPTION();
*/
/*
		for (int i = 0; i < nP * nP; i++) {
		for (int j = 0; j < nP * nP; j++) {
			printf("%i %i %1.15e\n", i, j, dRedistributionMaps(ixSecond,i,j));
		}
		}
		_EXCEPTION();
*/
	}

	// Current overlap face
	int ixOverlap = 0;

	// Loop through all faces on meshInput
	for (int ixFirst = 0; ixFirst < meshInput.faces.size(); ixFirst++) {

		// Output every 100 elements
		if (ixFirst % 1000 == 0) {
                        Announce("Element %i/%i", ixFirst, meshInput.faces.size());
		}

		// This Face
		const Face & faceFirst = meshInput.faces[ixFirst];

		// Find the set of Faces that overlap faceFirst
		int ixOverlapBegin = ixOverlap;
		int ixOverlapEnd = ixOverlapBegin;
	
		for (; ixOverlapEnd < meshOverlap.faces.size(); ixOverlapEnd++) {

			if (meshOverlap.vecSourceFaceIx[ixOverlapEnd] != ixFirst) {
				break;
			}
		}

		int nOverlapFaces = ixOverlapEnd - ixOverlapBegin;

		// Create a new Mesh representing the division of target finite
		// elements associated with this finite volume.
		Mesh meshThisElement;
		meshThisElement.faces.reserve(nOverlapFaces * nP * nP);
		meshThisElement.vecTargetFaceIx.reserve(nOverlapFaces * nP * nP);

		for (int i = ixOverlapBegin; i < ixOverlapEnd; i++) {

			int iTargetFace = meshOverlap.vecTargetFaceIx[i];

			int iSubElementBegin =  iTargetFace      * nP * nP;
			int iSubElementEnd   = (iTargetFace + 1) * nP * nP;

			int iSubElement = iSubElementBegin;
			for (int p = 0; p < nP; p++) {
			for (int q = 0; q < nP; q++) {

				// Calculate overlap polygon between sub-element
				// and finite volume
				NodeVector nodevecOutput;

				GenerateOverlapFace<MeshUtilitiesFuzzy, Node>(
					meshInput,
					meshTargetSubElement,
					ixFirst,
					iSubElementBegin + p * nP + q,
					nodevecOutput);
/*
				if (nodevecOutput.size() < 3) {
					continue;
				}

				Face faceNew(nodevecOutput.size());
				for (int n = 0; n < nodevecOutput.size(); n++) {
					faceNew.SetNode(n, n);
				}
				Real dArea = CalculateFaceArea(faceNew, nodevecOutput);

				if (dArea < 1.0e-13) {
					continue;
				}
*/
/*
				if (dataGLLNodes[p][q][iTargetFace] - 1 == 3) {
					printf("%1.15e %1.15e %1.15e\n", dArea, dataGLLJacobian[p][q][iTargetFace], dataGLLNodalArea[dataGLLNodes[p][q][iTargetFace] - 1]);
				}
*/
				Face faceNew(nodevecOutput.size());
				for (int n = 0; n < nodevecOutput.size(); n++) {
					meshThisElement.nodes.push_back(nodevecOutput[n]);
					faceNew.SetNode(n, meshThisElement.nodes.size()-1);
				}

				meshThisElement.faces.push_back(faceNew);

				meshThisElement.vecTargetFaceIx.push_back(
					dataGLLNodes(p,q,iTargetFace) - 1);

				iSubElement++;
			}
			}
		}

		if (meshThisElement.faces.size() != nOverlapFaces * nP * nP) {
			_EXCEPTIONT("Logic error");
		}

		// Build integration array
		DataArray2D<double> dIntArray;

		BuildIntegrationArray(
			meshInput,
			meshThisElement,
			triquadrule,
			ixFirst,
			0,
			meshThisElement.faces.size(),
			nOrder,
			dIntArray);

		// Set of Faces to use in building the reconstruction and associated
		// distance metric.
		AdjacentFaceVector vecAdjFaces;

//#ifdef RECTANGULAR_TRUNCATION
//		GetAdjacentFaceVectorByNode(
//#endif
//#ifdef TRIANGULAR_TRUNCATION
		GetAdjacentFaceVectorByEdge(
//#endif
			meshInput,
			ixFirst,
			nRequiredFaceSetSize,
			vecAdjFaces);

		// Number of adjacent Faces
		int nAdjFaces = vecAdjFaces.size();

		// Determine the conservative constraint equation
		DataArray1D<double> dConstraint(nCoefficients);

		double dFirstArea = meshInput.vecFaceArea[ixFirst];

		for (int p = 0; p < nCoefficients; p++) {
			for (int j = 0; j < meshThisElement.faces.size(); j++) {
				dConstraint[p] += dIntArray(p,j);
			}
			dConstraint[p] /= dFirstArea;
		}

		// Least squares arrays
		DataArray2D<double> dFitArray;
		DataArray1D<double> dFitWeights;
		DataArray2D<double> dFitArrayPlus;

		BuildFitArray(
			meshInput,
			triquadrule,
			ixFirst,
			vecAdjFaces,
			nOrder,
			nFitWeightsExponent,
			dConstraint,
			dFitArray,
			dFitWeights
		);

		// Compute the pseudoinverse fit array
		bool fSuccess =
			InvertFitArray_Corrected(
				dConstraint,
				dFitArray,
				dFitWeights,
				dFitArrayPlus
			);

		// Build the composition, which maps average values in adjacent cells
		// to the integrated values of the reconstruction in overlap faces.
		DataArray2D<double> dComposedArray(nAdjFaces, meshThisElement.faces.size());
		if (fSuccess) {
			for (int i = 0; i < nAdjFaces; i++) {
			for (int j = 0; j < meshThisElement.faces.size(); j++) {
			for (int k = 0; k < nCoefficients; k++) {
				dComposedArray(i,j) += dIntArray(k,j) * dFitArrayPlus(i,k);
			}
			}
			}

		// Unable to invert fit array, drop to 1st order.  In this case
		// dFitArrayPlus(0,0) = 1 and all other entries are zero.
		} else {
			for (int j = 0; j < meshThisElement.faces.size(); j++) {
				dComposedArray(0,j) += dIntArray(0,j);
			}
		}

		// Apply redistribution operator, which redistributes mass within
		// a finite element.
		DataArray2D<double> dRedistributedArray(nAdjFaces, meshThisElement.faces.size());

		for (int i = 0; i < nAdjFaces; i++) {
		for (int j = 0; j < meshThisElement.faces.size(); j++) {
			int ixSubElement = j % (nP * nP);
			int ixElement = j / (nP * nP);

			int ixSecondFace =
				meshOverlap.vecTargetFaceIx[ixOverlap + ixElement];

			for (int k = 0; k < nP * nP; k++) {
				dRedistributedArray(i,j) +=
					dComposedArray(i,ixElement * nP * nP + k)
					* dRedistributionMaps[ixSecondFace](ixSubElement,k);
			}
		}
		}

		// Put composed array into map
		for (int i = 0; i < vecAdjFaces.size(); i++) {
		for (int j = 0; j < meshThisElement.faces.size(); j++) {
			int ixFirstFace = vecAdjFaces[i].first;
			int ixSecondNode = meshThisElement.vecTargetFaceIx[j];

			smatMap(ixSecondNode, ixFirstFace) +=
				dRedistributedArray(i,j)
				/ dataGLLNodalArea[ixSecondNode];
				// meshThisElement.vecFaceArea[j];
				// dataGLLJacobian(ixS,ixT,ixSecondElement);
		}
		}

		// Increment the current overlap index
		ixOverlap += nOverlapFaces;

		//_EXCEPTION();
	}
}

///////////////////////////////////////////////////////////////////////////////

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
) {
	// NOTE: Reducing this quadrature rule order greatly affects error norms
	// Order of triangular quadrature rule
	const int TriQuadRuleOrder = 8;

	// Verify ReverseNodeArray has been calculated
	if (meshInput.revnodearray.size() == 0) {
		_EXCEPTIONT("ReverseNodeArray has not been calculated for meshInput");
	}
	if (meshInput.edgemap.size() == 0) {
		_EXCEPTIONT("EdgeMap has not been calculated for meshInput");
	}

	// Triangular quadrature rule
	TriangularQuadratureRule triquadrule(TriQuadRuleOrder);

	const DataArray2D<double> & dG = triquadrule.GetG();
	const DataArray1D<double> & dW = triquadrule.GetW();

	// Get SparseMatrix represntation of the OfflineMap
	SparseMatrix<double> & smatMap = mapRemap.GetSparseMatrix();

	// Fit weight exponent
	int nFitWeightsExponent = nOrder + 2;

	// Order of the finite element method
	int nP = dataGLLNodes.GetRows();

	// Sample coefficients
	DataArray2D<double> dSampleCoeff(nP, nP);

	// Number of elements needed
#ifdef RECTANGULAR_TRUNCATION
	int nCoefficients = nOrder * nOrder;
#endif
#ifdef TRIANGULAR_TRUNCATION 
	int nCoefficients = nOrder * (nOrder + 1) / 2;
#endif

	int nRequiredFaceSetSize = nCoefficients;

	// Announcemnets
	Announce("Triangular quadrature rule order %i", TriQuadRuleOrder);
	Announce("Number of coefficients: %i", nCoefficients);
	Announce("Required adjacency set size: %i", nRequiredFaceSetSize);
	Announce("Fit weights exponent: %i", nFitWeightsExponent);

	// Current overlap face
	int ixOverlap = 0;

	// Build the integration array for each element on meshOverlap
	DataArray3D<double> dGlobalIntArray(
		nCoefficients,
		meshOverlap.faces.size(),
		nP * nP);
/*
	// Build the mass matrix for each element on meshOutput
	DataArray3D<double> dMassMatrix;
	dMassMatrix.Initialize(
		meshOutput.faces.size(),
		nP * nP,
		nP * nP);
*/
	// Number of overlap Faces per source Face
	DataArray1D<int> nAllOverlapFaces(meshInput.faces.size());
	DataArray1D<int> nAllTotalOverlapTriangles(meshInput.faces.size());

	for (int ixFirst = 0; ixFirst < meshInput.faces.size(); ixFirst++) {

		int ixOverlapTemp = ixOverlap;
		for (; ixOverlapTemp < meshOverlap.faces.size(); ixOverlapTemp++) {

			const Face & faceOverlap = meshOverlap.faces[ixOverlapTemp];

			if (meshOverlap.vecSourceFaceIx[ixOverlapTemp] != ixFirst) {
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

    // generic triangle used for area computation, for triangles around the center of overlap face;
    // used for overlap faces with more than 4 edges;
    // nodes array will be set for each triangle;
    // these triangles are not part of the mesh structure, they are just temporary during
    //   aforementioned decomposition.
    Face faceTri( 3 );
    NodeVector nodes( 3 );
    faceTri.SetNode( 0, 0 );
    faceTri.SetNode( 1, 1 );
    faceTri.SetNode( 2, 2 );

	for (int ixFirst = 0; ixFirst < meshInput.faces.size(); ixFirst++) {

		// Output every 100 elements
		if (ixFirst % 1000 == 0) {
                        Announce("Element %i/%i", ixFirst, meshInput.faces.size());
		}

		// This Face
		const Face & faceFirst = meshInput.faces[ixFirst];

		// Area of the First Face
		double dFirstArea = meshInput.vecFaceArea[ixFirst];

		// Coordinate axes
		Node nodeRef = GetFaceCentroid(faceFirst, meshInput.nodes);

		Node nodeA1 = meshInput.nodes[faceFirst[1]] - nodeRef;
		Node nodeA2 = meshInput.nodes[faceFirst[2]] - nodeRef;

		Node nodeC = CrossProduct(nodeA1, nodeA2);

		// Fit matrix
		DataArray2D<double> dFit(3,3);

		dFit(0,0) = nodeA1.x; dFit(0,1) = nodeA1.y; dFit(0,2) = nodeA1.z;
		dFit(1,0) = nodeA2.x; dFit(1,1) = nodeA2.y; dFit(1,2) = nodeA2.z;
		dFit(2,0) = nodeC.x;  dFit(2,1) = nodeC.y;  dFit(2,2) = nodeC.z;

		// Number of overlapping Faces and triangles
		int nOverlapFaces = nAllOverlapFaces[ixFirst];
		int nTotalOverlapTriangles = nAllTotalOverlapTriangles[ixFirst];

		// Loop through all Overlap Faces
		for (int i = 0; i < nOverlapFaces; i++) {

			// Quantities from the overlap Mesh
			const Face & faceOverlap = meshOverlap.faces[ixOverlap + i];

			const NodeVector & nodesOverlap = meshOverlap.nodes;

			// Quantities from the Second Mesh
			int ixSecond = meshOverlap.vecTargetFaceIx[ixOverlap + i];

			const NodeVector & nodesSecond = meshOutput.nodes;

			const Face & faceSecond = meshOutput.faces[ixSecond];
			int nbEdges = faceOverlap.edges.size();
            int nOverlapTriangles = 1;
            Node center; // not used if nbEdges == 3
            if (nbEdges > 3) { // decompose from center in this case
                nOverlapTriangles = nbEdges;
                for (int k = 0; k < nbEdges; k++) {
                    const Node &node = nodesOverlap[faceOverlap[k]];
                    center = center + node;
                }
                center = center / nbEdges;
                double magni = sqrt(
                        center.x * center.x + center.y * center.y
                                + center.z * center.z);
                center = center / magni; // project back on sphere of radius 1
            }

            Node node0, node1, node2;
            double dTriArea;
			// Loop over all sub-triangles of this Overlap Face
			for (int j = 0; j < nOverlapTriangles; j++) {

			    if (nbEdges == 3) // will come here only once, nOverlapTriangles == 1 in this case
                {
                    node0 = nodesOverlap[faceOverlap[0]];
                    node1 = nodesOverlap[faceOverlap[1]];
                    node2 = nodesOverlap[faceOverlap[2]];
                    dTriArea = CalculateFaceArea(faceOverlap, nodesOverlap);
                }
                else // decompose polygon in triangles around the center
                {
                    node0 = center;
                    node1 = nodesOverlap[faceOverlap[j]];
                    int j1 = (j + 1) % nbEdges;
                    node2 = nodesOverlap[faceOverlap[j1]];
                    nodes[0] = center;
                    nodes[1] = node1;
                    nodes[2] = node2;
                    dTriArea = CalculateFaceArea(faceTri, nodes);
                }

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

					dX[0] -= nodeRef.x;
					dX[1] -= nodeRef.y;
					dX[2] -= nodeRef.z;

					// Find the coefficients for this point of the polynomial
					int n = 3;
					int nrhs = 1;
					int lda = 3;
					int ipiv[3];
					int ldb = 3;
					int info;

					DataArray2D<double> dFitTemp;
					dFitTemp = dFit;
					dgesv_(
						&n, &nrhs, &(dFitTemp(0,0)), &lda, ipiv, dX, &ldb, &info);

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
/*
					// Check inverse map value
					if ((dAlpha < -1.0e-12) || (dAlpha > 1.0 + 1.0e-12) ||
						(dBeta  < -1.0e-12) || (dBeta  > 1.0 + 1.0e-12)
					) {
						_EXCEPTION2("Inverse Map out of range (%1.5e %1.5e)",
							dAlpha, dBeta);
					}
*/
					// Sample the finite element at this point
					SampleGLLFiniteElement(
						nMonotoneType,
						nP,
						dAlpha,
						dBeta,
						dSampleCoeff);

					// Sample this point in the GLL element
					int ixs = 0;
					for (int s = 0; s < nP; s++) {
					for (int t = 0; t < nP; t++) {
/*
						int ixu = 0;
						for (int u = 0; u < nP; u++) {
						for (int v = 0; v < nP; v++) {
							dMassMatrix(ixSecond,ixs,ixu) +=
								  dSampleCoeff(s,t)
								* dSampleCoeff(u,v)
								* dW[k]
								* dTriArea;

							ixu++;
						}
						}
*/
						int ixp = 0;

#ifdef RECTANGULAR_TRUNCATION
						for (int p = 0; p < nOrder; p++) {
						for (int q = 0; q < nOrder; q++) {
#endif
#ifdef TRIANGULAR_TRUNCATION 
						for (int p = 0; p < nOrder; p++) {
						for (int q = 0; q < nOrder - p; q++) {
#endif

							dGlobalIntArray(ixp,ixOverlap + i,ixs) +=
								  dSampleCoeff(s,t)
								* IPow(dX[0], p)
								* IPow(dX[1], q)
								* dW[k]
								* dTriArea
								/ dataGLLJacobian(s,t,ixSecond);

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

/*
	// Calculate inverse mass matrix over all target elements
	ixOverlap = 0;

	for (int ixSecond = 0; ixSecond < meshOutput.faces.size(); ixSecond++) {

		// Calculate inverse of the mass matrix
		int m = nP * nP;
		int n = nP * nP;
		int lda = nP * nP;
		int info;

		DataArray1D<int> iPIV;
		iPIV.Initialize(nP * nP);

		DataArray1D<double> dWork;
		dWork.Initialize(nP * nP);

		int lWork = nP * nP;

		dgetrf_(&m, &n, &(dMassMatrix(ixSecond,0,0)),
			&lda, &(iPIV[0]), &info);

		if (info != 0) {
			_EXCEPTIONT("Mass matrix triangulation error");
		}

		dgetri_(&n, &(dMassMatrix(ixSecond,0,0)),
			&lda, &(iPIV[0]), &(dWork[0]), &lWork, &info);

		if (info != 0) {
			_EXCEPTIONT("Mass matrix inversion error");
		}
	}

	// Apply inverse mass matrix
	ixOverlap = 0;

	DataArray1D<double> dTemp;
	dTemp.Initialize(nP * nP);

	for (int ixFirst = 0; ixFirst < meshInput.faces.size(); ixFirst++) {

		// Number of overlapping Faces and triangles
		int nOverlapFaces = nAllOverlapFaces[ixFirst];

		// Loop through all Overlap Faces
		for (int i = 0; i < nOverlapFaces; i++) {

			// Quantities from the Second Mesh
			int ixSecond = meshOverlap.vecTargetFaceIx[ixOverlap + i];

			for (int ixp = 0; ixp < dGlobalIntArray.GetRows(); ixp++) {

				memcpy(&(dTemp[0]),
					&(dGlobalIntArray(ixp,ixOverlap + i,0)),
					nP * nP * sizeof(double));

				memset(&(dGlobalIntArray(ixp,ixOverlap + i,0)),
					0, nP * nP * sizeof(double));

				for (int s = 0; s < nP * nP; s++) {
				for (int t = 0; t < nP * nP; t++) {
					dGlobalIntArray(ixp,ixOverlap + i,s) +=
						dMassMatrix(ixSecond,s,t) * dTemp[t];
				}
				}
			}
		}

		ixOverlap += nOverlapFaces;
	}
*/
	// Reverse map
	std::vector< std::vector<int> > vecReverseFaceIx;
	vecReverseFaceIx.resize(meshOutput.faces.size());
	for (int i = 0; i < meshOverlap.faces.size(); i++) {
		int ixSecond = meshOverlap.vecTargetFaceIx[i];

		vecReverseFaceIx[ixSecond].push_back(i);
	}
/*
	for (int ixOverlap = 0; ixOverlap < meshOverlap.faces.size(); ixOverlap++) {
		int ixSecond = meshOverlap.vecTargetFaceIx[ixOverlap];

		for (int ixp = 0; ixp < dGlobalIntArray.GetRows(); ixp++) {
		for (int ixs = 0; ixs < dGlobalIntArray.GetSubColumns(); ixs++) {
			dGlobalIntArray(ixp,ixOverlap,ixs) *=
				dataGLLJacobian(ixs/nP,ixs%nP,ixSecond)
				/ dNumericalTargetArea(ixSecond,ixs);
		}
		}
	}
*/
	// Force consistency and conservation
	for (int ixSecond = 0; ixSecond < meshOutput.faces.size(); ixSecond++) {

		if (vecReverseFaceIx[ixSecond].size() == 0) {
			continue;
		}

		DataArray2D<double> dCoeff(
			nP * nP,
			vecReverseFaceIx[ixSecond].size());

		for (int i = 0; i < vecReverseFaceIx[ixSecond].size(); i++) {
			int ixOverlap = vecReverseFaceIx[ixSecond][i];

			for (int s = 0; s < nP * nP; s++) {
				dCoeff[s][i] = dGlobalIntArray(0,ixOverlap,s);
			}
		}

		// Target areas
		DataArray1D<double> vecTargetArea(nP * nP);

		for (int s = 0; s < nP * nP; s++) {
			vecTargetArea[s] =
				dataGLLJacobian(s/nP,s%nP,ixSecond);
		}

		// Source areas
		DataArray1D<double> vecSourceArea(vecReverseFaceIx[ixSecond].size());

		for (int i = 0; i < vecReverseFaceIx[ixSecond].size(); i++) {
			int ixOverlap = vecReverseFaceIx[ixSecond][i];
			vecSourceArea[i] = meshOverlap.vecFaceArea[ixOverlap];
		}

		if (!fNoConservation) {
			ForceIntArrayConsistencyConservation(
				vecSourceArea,
				vecTargetArea,
				dCoeff,
				(nMonotoneType != 0));
		}
/*
		for (int i = 0; i < dCoeff.GetRows(); i++) {
			double dConsistency = 0.0;
			for (int j = 0; j < dCoeff.GetColumns(); j++) {
				dConsistency += dCoeff(i,j);
			}
			//printf("%1.15e\n", dConsistency);
		}
*/
		for (int i = 0; i < vecReverseFaceIx[ixSecond].size(); i++) {
			int ixOverlap = vecReverseFaceIx[ixSecond][i];

			for (int s = 0; s < nP * nP; s++) {
				//printf("%1.15e %1.15e\n", dGlobalIntArray[0][ixOverlap][s], dCoeff[s][i]);
				dGlobalIntArray(0,ixOverlap,s) = dCoeff(s,i);
			}
		}

/*
		for (int i = 0; i < dCoeff.GetRows(); i++) {
			double dConsistency = 0.0;
			for (int j = 0; j < dCoeff.GetColumns(); j++) {
				dConsistency += dCoeff(i,j);
			}
			printf("%1.15e\n", dConsistency);
		}

		for (int i = 0; i < dCoeff.GetRows(); i++) {
			int ixFirst = vecReverseFaceIx(ixSecond,i)

			for (int s = 0; s < dCoeff.GetColumns(); s++) {
				vecTargetArea[i] += dCoeff(i,s)
					* dataGLLJacobian(s/nP,s%nP,ixSecond)
					 meshInput.vecFaceArea[ixFirst];
			}
			printf("%1.15e\n", vecTargetArea[i]);
		}
*/
/*
		double dConsistency = 0.0;
		double dConservation = 0.0;

		int ixFirst = meshOverlap.vecSourceFaceIx[i];
		int ixSecond = meshOverlap.vecTargetFaceIx[i];

		for (int s = 0; s < nP * nP; s++) {
			//dConsistency += dGlobalIntArray(0,i,s)
			dConservation += dGlobalIntArray(0,i,s)
				* dataGLLJacobian(s/nP,s%nP,ixSecond)
				/ meshInput.vecFaceArea[ixFirst];

			printf("%1.15e\n", dataGLLJacobian(s/nP,s%nP,ixSecond));
		}

		//printf("Consistency: %1.15e\n", dConsistency);
		printf("Conservation: %1.15e\n", dConservation);

		_EXCEPTION();
*/
	}

/*
	// Check consistency
	DataArray2D<double> dIntSums(meshOutput.faces.size(), nP * nP);

	for (int i = 0; i < meshOverlap.faces.size(); i++) {
	for (int s = 0; s < nP * nP; s++) {
		int ixSecond = meshOverlap.vecTargetFaceIx[i];

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
	DataArray1D<double> dMassSums(meshInput.faces.size());

	for (int i = 0; i < meshOverlap.faces.size(); i++) {
	for (int s = 0; s < nP * nP; s++) {
		int ixFirst = meshOverlap.vecSourceFaceIx[i];
		int ixSecond = meshOverlap.vecTargetFaceIx[i];

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
	DataArray1D<double> dOverlapMass(meshOutput.faces.size());

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
		if (ixFirst % 1000 == 0) {
                        Announce("Element %i/%i", ixFirst, meshInput.faces.size());
		}

		// This Face
		const Face & faceFirst = meshInput.faces[ixFirst];

		// Area of the First Face
		double dFirstArea = meshInput.vecFaceArea[ixFirst];

		// Number of overlapping Faces and triangles
		int nOverlapFaces = nAllOverlapFaces[ixFirst];
		int nTotalOverlapTriangles = nAllTotalOverlapTriangles[ixFirst];
/*
		// Verify equal partition of mass in integration array
		double dTotal = 0.0;
		for (int i = 0; i < nOverlapFaces; i++) {
			int ixSecond = meshOverlap.vecTargetFaceIx[ixOverlap + i];

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
*/
		// Determine the conservative constraint equation
		DataArray1D<double> dConstraint(nCoefficients);

		for (int p = 0; p < nCoefficients; p++) {
		for (int i = 0; i < nOverlapFaces; i++) {
			int ixSecond = meshOverlap.vecTargetFaceIx[ixOverlap + i];

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

//#ifdef RECTANGULAR_TRUNCATION
//		GetAdjacentFaceVectorByNode(
//#endif
//#ifdef TRIANGULAR_TRUNCATION
		GetAdjacentFaceVectorByEdge(
//#endif
			meshInput,
			ixFirst,
			nRequiredFaceSetSize,
			vecAdjFaces);

		// Number of adjacent Faces
		int nAdjFaces = vecAdjFaces.size();

		for (int x = 0; x < nAdjFaces; x++) {
			if (vecAdjFaces[x].first == (-1)) {
				_EXCEPTION();
			}
		}

		// Build the fit operator
		DataArray2D<double> dFitArray;
		DataArray1D<double> dFitWeights;
		DataArray2D<double> dFitArrayPlus;

		BuildFitArray(
			meshInput,
			triquadrule,
			ixFirst,
			vecAdjFaces,
			nOrder,
			nFitWeightsExponent,
			dConstraint,
			dFitArray,
			dFitWeights
		);

/*
		DataArray1D<double> dRowSum;
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

		// Compute the pseudoinverse fit array
		bool fSuccess =
			InvertFitArray_Corrected(
				dConstraint,
				dFitArray,
				dFitWeights,
				dFitArrayPlus
			);

		// Build the composition, which maps average values in adjacent cells
		// to the integrated values of the reconstruction in overlap faces.
		DataArray2D<double> dComposedArray(nAdjFaces, nOverlapFaces * nP * nP);
		if (fSuccess) {
			for (int j = 0; j < nOverlapFaces; j++) {
				//int ixSecond = meshOverlap.vecTargetFaceIx[ixOverlap + j];

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

		// Unable to invert fit array, drop to 1st order.  In this case
		// dFitArrayPlus(0,0) = 1 and all other entries are zero.
		} else {
			for (int j = 0; j < nOverlapFaces; j++) {
				//int ixSecond = meshOverlap.vecTargetFaceIx[ixOverlap + j];

				for (int s = 0; s < nP * nP; s++) {
					dComposedArray[0][j * nP * nP + s] +=
						dGlobalIntArray[0][ixOverlap + j][s];
				}
			}
		}

		// Put composed array into map
		for (int i = 0; i < vecAdjFaces.size(); i++) {
		for (int j = 0; j < nOverlapFaces; j++) {
			int ixFirstFace = vecAdjFaces[i].first;
			int ixSecondFace = meshOverlap.vecTargetFaceIx[ixOverlap + j];

			for (int s = 0; s < nP; s++) {
			for (int t = 0; t < nP; t++) {

				int jx = j * nP * nP + s * nP + t;

				if (fContinuous) {
					int ixSecondNode = dataGLLNodes[s][t][ixSecondFace] - 1;

					smatMap(ixSecondNode, ixFirstFace) +=
						dComposedArray[i][jx]
						* dataGLLJacobian[s][t][ixSecondFace]
						/ dataGLLNodalArea[ixSecondNode];

				} else {
					int ixSecondNode = ixSecondFace * nP * nP + s * nP + t;

					smatMap(ixSecondNode, ixFirstFace) +=
						dComposedArray[i][jx];
				}
			}
			}
		}
		}

		// Increment the current overlap index
		ixOverlap += nOverlapFaces;
	}
}

///////////////////////////////////////////////////////////////////////////////

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
) {
	// Triangular quadrature rule
	TriangularQuadratureRule triquadrule(8);

	const DataArray2D<double> & dG = triquadrule.GetG();
	const DataArray1D<double> & dW = triquadrule.GetW();

	// Get SparseMatrix represntation of the OfflineMap
	SparseMatrix<double> & smatMap = mapRemap.GetSparseMatrix();

	// Sample coefficients
	DataArray2D<double> dSampleCoeffIn(nPin, nPin);
	DataArray2D<double> dSampleCoeffOut(nPout, nPout);

	// Build the integration array for each element on meshOverlap
	DataArray3D<double> dGlobalIntArray(
		nPin * nPin,
		meshOverlap.faces.size(),
		nPout * nPout);

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

	// Geometric area of each output node
	DataArray2D<double> dGeometricOutputArea(
		meshOutput.faces.size(), nPout * nPout);

	// Area of each overlap element in the output basis
	DataArray2D<double> dOverlapOutputArea(
		meshOverlap.faces.size(), nPout * nPout);

	// Loop through all faces on meshInput
	ixOverlap = 0;

	// generic triangle used for area computation, for triangles around the center of overlap face;
    // used for overlap faces with more than 4 edges;
    // nodes array will be set for each triangle;
    // these triangles are not part of the mesh structure, they are just temporary during
    //   aforementioned decomposition.
    Face faceTri( 3 );
    NodeVector nodes( 3 );
    faceTri.SetNode( 0, 0 );
    faceTri.SetNode( 1, 1 );
    faceTri.SetNode( 2, 2 );

	Announce("Building conservative distribution maps");
	for (int ixFirst = 0; ixFirst < meshInput.faces.size(); ixFirst++) {

		// Output every 100 elements
		if (ixFirst % 1000 == 0) {
                        Announce("Element %i/%i", ixFirst, meshInput.faces.size());
		}

		// Quantities from the First Mesh
		const Face & faceFirst = meshInput.faces[ixFirst];

		const NodeVector & nodesFirst = meshInput.nodes;

		// Number of overlapping Faces and triangles
		int nOverlapFaces = nAllOverlapFaces[ixFirst];
/*
		// Calculate total element Jacobian
		double dTotalJacobian = 0.0;
		for (int s = 0; s < nPin; s++) {
		for (int t = 0; t < nPin; t++) {
			dTotalJacobian += dataGLLJacobianIn[s][t][ixFirst];
		}
		}
*/
		// Loop through all Overlap Faces
		for (int i = 0; i < nOverlapFaces; i++) {

			// Quantities from the overlap Mesh
			const Face & faceOverlap = meshOverlap.faces[ixOverlap + i];

			const NodeVector & nodesOverlap = meshOverlap.nodes;

			// Quantities from the Second Mesh
			int ixSecond = meshOverlap.vecTargetFaceIx[ixOverlap + i];

			const NodeVector & nodesSecond = meshOutput.nodes;

			const Face & faceSecond = meshOutput.faces[ixSecond];

			int nbEdges = faceOverlap.edges.size();
            int nOverlapTriangles = 1;
            Node center; // not used if nbEdges == 3
            if (nbEdges > 3) { // decompose from center in this case
                nOverlapTriangles = nbEdges;
                for (int k = 0; k < nbEdges; k++) {
                    const Node &node = nodesOverlap[faceOverlap[k]];
                    center = center + node;
                }
                center = center / nbEdges;
                double magni = sqrt(
                        center.x * center.x + center.y * center.y
                                + center.z * center.z);
                center = center / magni; // project back on sphere of radius 1
            }

            Node node0, node1, node2;
            double dTriArea;

			// Loop over all sub-triangles of this Overlap Face
			for (int j = 0; j < nOverlapTriangles; j++) {

			    if (nbEdges == 3) // will come here only once, nOverlapTriangles == 1 in this case
                {
                    node0 = nodesOverlap[faceOverlap[0]];
                    node1 = nodesOverlap[faceOverlap[1]];
                    node2 = nodesOverlap[faceOverlap[2]];
                    dTriArea = CalculateFaceArea(faceOverlap, nodesOverlap);
                }
                else // decompose polygon in triangles around the center
                {
                    node0 = center;
                    node1 = nodesOverlap[faceOverlap[j]];
                    int j1 = (j + 1) % nbEdges;
                    node2 = nodesOverlap[faceOverlap[j1]];
                    nodes[0] = center;
                    nodes[1] = node1;
                    nodes[2] = node2;
                    dTriArea = CalculateFaceArea(faceTri, nodes);
                }

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

/*
					// Check inverse map value
					if ((dAlphaIn < 0.0) || (dAlphaIn > 1.0) ||
						(dBetaIn  < 0.0) || (dBetaIn  > 1.0)
					) {
						_EXCEPTION2("Inverse Map out of range (%1.5e %1.5e)",
							dAlphaIn, dBetaIn);
					}

					// Check inverse map value
					if ((dAlphaOut < 0.0) || (dAlphaOut > 1.0) ||
						(dBetaOut  < 0.0) || (dBetaOut  > 1.0)
					) {
						_EXCEPTION2("Inverse Map out of range (%1.5e %1.5e)",
							dAlphaOut, dBetaOut);
					}
*/
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

					// Overlap output area
					for (int s = 0; s < nPout; s++) {
					for (int t = 0; t < nPout; t++) {
						double dNodeArea =
							dSampleCoeffOut[s][t]
							* dW[k]
							* dTriArea;

						dOverlapOutputArea[ixOverlap + i][s * nPout + t] +=
							dNodeArea;

						dGeometricOutputArea[ixSecond][s * nPout + t] +=
							dNodeArea;
					}
					}

					// Compute overlap integral
					int ixp = 0;
					for (int p = 0; p < nPin; p++) {
					for (int q = 0; q < nPin; q++) {

						int ixs = 0;
						for (int s = 0; s < nPout; s++) {
						for (int t = 0; t < nPout; t++) {

							// Sample the Second finite element at this point
							dGlobalIntArray[ixp][ixOverlap + i][ixs] +=
								  dSampleCoeffOut[s][t]
								* dSampleCoeffIn[p][q]
								* dW[k]
								* dTriArea;

							ixs++;
						}
						}

						ixp++;
					}
					}
				}
			}
		}

		// Coefficients
		DataArray2D<double> dCoeff(nOverlapFaces * nPout * nPout, nPin * nPin);

		for (int i = 0; i < nOverlapFaces; i++) {

			int ixSecondFace = meshOverlap.vecTargetFaceIx[ixOverlap + i];

			int ixp = 0;
			for (int p = 0; p < nPin; p++) {
			for (int q = 0; q < nPin; q++) {

				int ixs = 0;
				for (int s = 0; s < nPout; s++) {
				for (int t = 0; t < nPout; t++) {
					if (fabs(dOverlapOutputArea[ixOverlap + i][s * nPout + t]) < ReferenceTolerance) {
						continue;
					}

					dCoeff[i * nPout * nPout + ixs][ixp] =
						dGlobalIntArray[ixp][ixOverlap + i][ixs]
						/ dOverlapOutputArea[ixOverlap + i][s * nPout + t];

					ixs++;
				}
				}

				ixp++;
			}
			}
		}

		// Source areas
		DataArray1D<double> vecSourceArea(nPin * nPin);

		for (int p = 0; p < nPin; p++) {
		for (int q = 0; q < nPin; q++) {
			vecSourceArea[p * nPin + q] =
				dataGLLJacobianIn[p][q][ixFirst];
		}
		}

		// Target areas
		DataArray1D<double> vecTargetArea(nOverlapFaces * nPout * nPout);

		for (int i = 0; i < nOverlapFaces; i++) {
			int ixSecond = meshOverlap.vecTargetFaceIx[ixOverlap + i];

			int ixs = 0;
			for (int s = 0; s < nPout; s++) {
			for (int t = 0; t < nPout; t++) {

				vecTargetArea[i * nPout * nPout + ixs] =
					dOverlapOutputArea[ixOverlap + i][nPout * s + t];

				ixs++;
			}
			}
		}

		// Force consistency and conservation
		if (!fNoConservation) {
			ForceIntArrayConsistencyConservation(
				vecSourceArea,
				vecTargetArea,
				dCoeff,
				(nMonotoneType != 0));
		}

		// Update global coefficients
		for (int i = 0; i < nOverlapFaces; i++) {

			int ixp = 0;
			for (int p = 0; p < nPin; p++) {
			for (int q = 0; q < nPin; q++) {

				int ixs = 0;
				for (int s = 0; s < nPout; s++) {
				for (int t = 0; t < nPout; t++) {

					dGlobalIntArray[ixp][ixOverlap + i][ixs] =
						dCoeff[i * nPout * nPout + ixs][ixp]
						* dOverlapOutputArea[ixOverlap + i][s * nPout + t];

					ixs++;
				}
				}

				ixp++;
			}
			}
		}
/*
		// Check column sums (conservation)
		for (int i = 0; i < nPin * nPin; i++) {
			double dColSum = 0.0;
			for (int j = 0; j < nOverlapFaces * nPout * nPout; j++) {
				dColSum += dCoeff[j][i] * vecTargetArea[j];
			}
			printf("Col %i: %1.15e\n", i, dColSum / vecSourceArea[i]);
		}

		// Check row sums (consistency)
		for (int j = 0; j < nOverlapFaces * nPout * nPout; j++) {
			double dRowSum = 0.0;
			for (int i = 0; i < nPin * nPin; i++) {
				dRowSum += dCoeff[j][i];
			}
			printf("Row %i: %1.15e\n", j, dRowSum);
		}
		_EXCEPTION();
*/

		// Increment the current overlap index
		ixOverlap += nOverlapFaces;
	}

	// Build redistribution map within target element
	Announce("Building redistribution maps on target mesh");
	DataArray1D<double> dRedistSourceArea(nPout * nPout);
	DataArray1D<double> dRedistTargetArea(nPout * nPout);
	std::vector< DataArray2D<double> > dRedistributionMaps;
	dRedistributionMaps.resize(meshOutput.faces.size());

	for (int ixSecond = 0; ixSecond < meshOutput.faces.size(); ixSecond++) {

		dRedistributionMaps[ixSecond].Allocate(
			nPout * nPout, nPout * nPout);

		for (int i = 0; i < nPout * nPout; i++) {
			dRedistributionMaps[ixSecond][i][i] = 1.0;
		}

		for (int s = 0; s < nPout * nPout; s++) {
			dRedistSourceArea[s] =
				dGeometricOutputArea[ixSecond][s];
		}

		for (int s = 0; s < nPout * nPout; s++) {
			dRedistTargetArea[s] =
				dataGLLJacobianOut[s/nPout][s%nPout][ixSecond];
		}

		if (!fNoConservation) {
			ForceIntArrayConsistencyConservation(
				dRedistSourceArea,
				dRedistTargetArea,
				dRedistributionMaps[ixSecond],
				(nMonotoneType != 0));

			for (int s = 0; s < nPout * nPout; s++) {
			for (int t = 0; t < nPout * nPout; t++) {
				dRedistributionMaps[ixSecond][s][t] *=
					dRedistTargetArea[s] / dRedistSourceArea[t];
			}
			}
		}
	}

	// Construct the total geometric area
	DataArray1D<double> dTotalGeometricArea(dataNodalAreaOut.GetRows());
	for (int ixSecond = 0; ixSecond < meshOutput.faces.size(); ixSecond++) {
		for (int s = 0; s < nPout; s++) {
		for (int t = 0; t < nPout; t++) {
			dTotalGeometricArea[dataGLLNodesOut[s][t][ixSecond] - 1]
				+= dGeometricOutputArea[ixSecond][s * nPout + t];
		}
		}
	}

	// Compose the integration operator with the output map
	ixOverlap = 0;

	Announce("Assembling map");

	// Map from source DOFs to target DOFs with redistribution applied
	DataArray2D<double> dRedistributedOp(
		nPin * nPin, nPout * nPout);

	for (int ixFirst = 0; ixFirst < meshInput.faces.size(); ixFirst++) {

		// Output every 100 elements
		if (ixFirst % 1000 == 0) {
                        Announce("Element %i/%i", ixFirst, meshInput.faces.size());
		}

		// This Face
		const Face & faceFirst = meshInput.faces[ixFirst];

		// Number of overlapping Faces and triangles
		int nOverlapFaces = nAllOverlapFaces[ixFirst];

		// Put composed array into map
		for (int j = 0; j < nOverlapFaces; j++) {
			int ixSecondFace = meshOverlap.vecTargetFaceIx[ixOverlap + j];

			dRedistributedOp.Zero();
			for (int p = 0; p < nPin * nPin; p++) {
			for (int s = 0; s < nPout * nPout; s++) {
				for (int t = 0; t < nPout * nPout; t++) {
					dRedistributedOp[p][s] +=
						dRedistributionMaps[ixSecondFace][s][t]
						* dGlobalIntArray[p][ixOverlap + j][t];
				}
			}
			}

			int ixp = 0;
			for (int p = 0; p < nPin; p++) {
			for (int q = 0; q < nPin; q++) {

				int ixFirstNode;
				if (fContinuousIn) {
					ixFirstNode = dataGLLNodesIn[p][q][ixFirst] - 1;
				} else {
					ixFirstNode = ixFirst * nPin * nPin + p * nPin + q;
				}

				int ixs = 0;
				for (int s = 0; s < nPout; s++) {
				for (int t = 0; t < nPout; t++) {

					int ixSecondNode;
					if (fContinuousOut) {
						ixSecondNode = dataGLLNodesOut[s][t][ixSecondFace] - 1;

						if (!fNoConservation) {
							smatMap(ixSecondNode, ixFirstNode) +=
								dRedistributedOp[ixp][ixs]
								/ dataNodalAreaOut[ixSecondNode];
						} else {
							smatMap(ixSecondNode, ixFirstNode) +=
								dRedistributedOp[ixp][ixs]
								/ dTotalGeometricArea[ixSecondNode];
						}

					} else {
						ixSecondNode =
							ixSecondFace * nPout * nPout + s * nPout + t;

						if (!fNoConservation) {
							smatMap(ixSecondNode, ixFirstNode) +=
								dRedistributedOp[ixp][ixs]
								/ dataGLLJacobianOut[s][t][ixSecondFace];
						} else {
							smatMap(ixSecondNode, ixFirstNode) +=
								dRedistributedOp[ixp][ixs]
								/ dGeometricOutputArea[ixSecondFace][s * nPout + t];
						}
					}

					ixs++;
				}
				}

				ixp++;
			}
			}
		}

		// Increment the current overlap index
		ixOverlap += nOverlapFaces;
	}
}

///////////////////////////////////////////////////////////////////////////////

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
) {
	// Gauss-Lobatto quadrature within Faces
	DataArray1D<double> dGL;
	DataArray1D<double> dWL;

	GaussLobattoQuadrature::GetPoints(nPout, 0.0, 1.0, dGL, dWL);

	// Utilities
	MeshUtilitiesFuzzy utils;

	// Get SparseMatrix represntation of the OfflineMap
	SparseMatrix<double> & smatMap = mapRemap.GetSparseMatrix();

	// Sample coefficients
	DataArray2D<double> dSampleCoeffIn(nPin, nPin);

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

	// Number of times this point was found
	DataArray1D<bool> fSecondNodeFound(dataNodalAreaOut.GetRows());

	// Loop through all faces on meshInput
	ixOverlap = 0;

	for (int ixFirst = 0; ixFirst < meshInput.faces.size(); ixFirst++) {

		// Output every 100 elements
		if (ixFirst % 1000 == 0) {
                        Announce("Element %i/%i", ixFirst, meshInput.faces.size());
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

			// Quantities from the Second Mesh
			int ixSecond = meshOverlap.vecTargetFaceIx[ixOverlap + i];

			const NodeVector & nodesSecond = meshOutput.nodes;

			const Face & faceSecond = meshOutput.faces[ixSecond];

			// Loop through all nodes on the second face
			int ixs = 0;
			for (int s = 0; s < nPout; s++) {
			for (int t = 0; t < nPout; t++) {

				int ixSecondNode;
				if (fContinuousOut) {
					ixSecondNode = dataGLLNodesOut[s][t][ixSecond] - 1;
				} else {
					ixSecondNode =
						ixSecond * nPout * nPout + s * nPout + t;
				}

				if (ixSecondNode >= fSecondNodeFound.GetRows()) {
					_EXCEPTIONT("Logic error");
				}

				// Check if this node has been found already
				if (fSecondNodeFound[ixSecondNode]) {
					continue;
				}

				// Check this node
				Node node;
				Node dDx1G;
				Node dDx2G;

				ApplyLocalMap(
					faceSecond,
					nodesSecond,
					dGL[t],
					dGL[s],
					node,
					dDx1G,
					dDx2G);

				// Find the components of this quadrature point in the basis
				// of the first Face.
				double dAlphaIn;
				double dBetaIn;

				ApplyInverseMap(
					faceFirst,
					nodesFirst,
					node,
					dAlphaIn,
					dBetaIn);

				// Check if this node is within the first Face
				if ((dAlphaIn < -1.0e-10) || (dAlphaIn > 1.0 + 1.0e-10) ||
					(dBetaIn  < -1.0e-10) || (dBetaIn  > 1.0 + 1.0e-10)
				) {
					continue;
				}

				// Node is within the overlap region, mark as found
				fSecondNodeFound[ixSecondNode] = true;

				// Sample the First finite element at this point
				SampleGLLFiniteElement(
					nMonotoneType,
					nPin,
					dAlphaIn,
					dBetaIn,
					dSampleCoeffIn);

				// Add to map
				int ixp = 0;
				for (int p = 0; p < nPin; p++) {
				for (int q = 0; q < nPin; q++) {

					int ixFirstNode;
					if (fContinuousIn) {
						ixFirstNode = dataGLLNodesIn[p][q][ixFirst] - 1;
					} else {
						ixFirstNode =
							ixFirst * nPin * nPin + p * nPin + q;
					}

					smatMap(ixSecondNode, ixFirstNode) +=
						dSampleCoeffIn[p][q];
				}
				}
			}
			}
		}

		// Increment the current overlap index
		ixOverlap += nOverlapFaces;
	}

	// Check for missing samples
	for (int i = 0; i < fSecondNodeFound.GetRows(); i++) {
		if (!fSecondNodeFound[i]) {
			_EXCEPTION1("Can't sample point %i", i);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

