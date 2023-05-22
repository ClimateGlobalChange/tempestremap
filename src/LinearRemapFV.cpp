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
#include "GaussLobattoQuadrature.h"
#include "TriangularQuadrature.h"
#include "MeshUtilitiesFuzzy.h"
#include "OverlapMesh.h"
#include "triangle.h"
#include "kdtree.h"
#include "DataArray3D.h"

#include "Announce.h"
#include "MathHelper.h"

#include <cstring>
#include <map>
#include <vector>

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

			// signal to not participate, because it is a ghost target
			if( ixSecondFace < 0 ) continue;  // skip and do not do anything

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

		if( nOverlapFaces == 0 ) continue;

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

			// signal to not participate, because it is a ghost target
			if( ixSecondFace < 0 ) continue;  // skip and do not do anything

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

	// Check mask size
	_ASSERT((meshInput.vecMask.size() == 0) || (meshInput.vecMask.size() == meshInput.faces.size()));

	if (meshInput.vecMask.size() != 0) {
		Announce("Source mesh contains mask information which will be used in map calculation");
	}

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

		// Check mask
		if (meshInput.vecMask.size() != 0) {
			if (meshInput.vecMask[ixFirst] == 0) {
 		               // Increment the current overlap index
                		ixOverlap += nOverlapFaces;

				continue;
			}
		}
		
		if( nOverlapFaces == 0 ) continue;

		// Loop through all overlap faces associated with this source face
		for (int j = 0; j < nOverlapFaces; j++) {
			int iTargetFace = meshOverlap.vecTargetFaceIx[ixOverlap + j];

			// signal to not participate, because it is a ghost target
			if( iTargetFace < 0 ) continue;  // skip and do not do anything

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

					// Check mask
					if (meshInput.vecMask.size() != 0) {
						if (meshInput.vecMask[iNearestFace] == 0) {
							iNearestFace = ixFirst;
						}
					}

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

							// Check mask
							if (meshInput.vecMask.size() != 0) {
								if (meshInput.vecMask[facepair[1]] == 0) {
									continue;
								}
							}

							// Add contribution
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

							// Check mask
							if (meshInput.vecMask.size() != 0) {
								if (meshInput.vecMask[facepair[0]] == 0) {
									continue;
								}
							}

							// Add contribution
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

void LinearRemapIntegratedTriangulation(
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

	// Get SparseMatrix representation of the OfflineMap
	SparseMatrix<double> & smatMap = mapRemap.GetSparseMatrix();

	// Array of global indices
	std::vector <std::vector<int>> vecGlobalIndexI(6);

	// Vector storing meshes of each panel
	std::vector<Mesh> vecMesh(6);

	// kd-tree for each panel
	std::vector<kdtree *> vecKDTreePanelI(6);

	// Arrays of panel boundaries in lat/lon coordinates (radians)
	DataArray2D<double> dPanelLat(6,2);
	DataArray2D<double> dPanelLon(6,2);

	// Width of buffer zone for each panel (degrees)
	double dBuff = 30;

	// Define equatorial panel boundaries
	for (int i = 0; i < 4; i++){

		dPanelLon(i,0) = 90*M_PI*i/180;  //Left boundary
		dPanelLon(i,1) = 90*M_PI*(i+1)/180;  //Right Boundary
		dPanelLat(i,1) = 45*M_PI/180;  //Upper Boundary
		dPanelLat(i,0) = -45*M_PI/180;  //Lower Boundary

	}

	// Define polar panels

	dPanelLon(4,0) = 0; dPanelLon(4,1) = 2*M_PI;
	dPanelLat(4,0) = 45*M_PI/180; dPanelLat(4,1) = M_PI/2;

	dPanelLon(5,0) = 0; dPanelLon(5,1) = 2*M_PI;
	dPanelLat(5,1) = -45*M_PI/180; dPanelLat(5,0) = -M_PI/2;


	// Generate the Delaunay triangulation for each of the 6 panels
	for (int i = 0; i < 6; i++){

		// Array storing the coordinates of the face centroids of panel i
		std::vector<std::vector<double>> dPanelCentroid(2);

		// Coordinates of tangent plane point of tangency for panel i
		double dLatTan;
		double dLonTan;

		if ((0 <= i) && (i <= 3)){

			dLatTan = 0;
			dLonTan = (i+1)*45*M_PI/180 + 45*i*M_PI/180;
		}
		else if (i == 4){

			dLatTan = M_PI/2;
			dLonTan = 0;

		}
		else{

			dLonTan = 0;
			dLatTan = -M_PI/2;

		}

		// Number of source faces centers for panel i; this is the size
		// of in.pointlist
		int nNodes = 0;

		for (int j = 0; j < meshInput.faces.size(); j++){

			// Get coordinates of centroid of current face

			double dLonRad0;
			double dLatRad0;

			const Face & face = meshInput.faces[j];

			Node node = GetFaceCentroid(face,meshInput.nodes);

			node = node.Normalized();
			//node.Normalized();

			XYZtoRLL_Rad(node.x,node.y,node.z,dLonRad0,dLatRad0);

			// Determine if centroid is in panel i plus buffer zone
			bool fPanelContainsCentroid = false;

			if ((0 <= i) && (i <= 3)){

				if ((dPanelLat(i,0) - dBuff*M_PI/180 <= dLatRad0) && (dLatRad0 <= dPanelLat(i,1)+dBuff*M_PI/180)){

					if ( i == 0 ){

						if ( ((2*M_PI - dBuff*M_PI/180 <= dLonRad0) && (dLonRad0 <= 2*M_PI)) ||
							 ((dPanelLon(i,0) <= dLonRad0) && (dLonRad0 <= dPanelLon(i,1) + dBuff*M_PI/180)) ) {

							fPanelContainsCentroid = true;
						}
					}

					else if ( i == 3 ){

						if ( ((0 <= dLonRad0) && (dLonRad0 <= 0 + dBuff*M_PI/180)) ||
							 ((dPanelLon(i,0) - dBuff*M_PI/180 <= dLonRad0) && (dLonRad0 <= dPanelLon(i,1) + dBuff*M_PI/180)) ){

							fPanelContainsCentroid = true;

						}

					}

					else {

						if ( (dPanelLon(i,0) - dBuff*M_PI/180 <= dLonRad0) && (dLonRad0 <= dPanelLon(i,1) + dBuff*M_PI/180)){

							fPanelContainsCentroid = true;
						}


					}

				}
			}
			else if (i == 4){

				if ((dPanelLat(i,0) - dBuff*M_PI/180 <= dLatRad0) && (dLatRad0 <= dPanelLat(i,1))){

					fPanelContainsCentroid = true;

				}

			}

			else {

				if ((dPanelLat(i,0) <= dLatRad0) && (dLatRad0 <= dPanelLat(i,1) + dBuff*M_PI/180)){
					fPanelContainsCentroid = true;

				}

			}

			if(fPanelContainsCentroid){

				vecGlobalIndexI[i].push_back(j);
				nNodes++;

				// Gnomonic projection of centroid coordinates
				double dGX;
				double dGY;
				GnomonicProjection(dLonTan,dLatTan,dLonRad0,dLatRad0,dGX,dGY);

				// Add Gnomonic coordinates to vector
				dPanelCentroid[0].push_back(dGX);
				dPanelCentroid[1].push_back(dGY);

			}
		}

		// Structures for Delaunay triangulation

		struct triangulateio in, out, vorout;

		in.numberofpoints           = nNodes;
		in.numberofpointattributes  = 0;
		in.numberofsegments         = 0;
		in.numberofholes            = 0;
		in.numberofregions          = 0;
		in.pointlist                = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));
		in.segmentlist              = (int  *) malloc(in.numberofsegments * 2 * sizeof(int));;
		in.pointattributelist       = (REAL *) NULL;
		in.pointmarkerlist          = (int  *) NULL;
		in.trianglelist             = (int  *) NULL;
		in.triangleattributelist    = (REAL *) NULL;
		in.neighborlist             = (int  *) NULL;
		in.segmentmarkerlist        = (int  *) NULL;
		in.edgelist                 = (int  *) NULL;
		in.edgemarkerlist           = (int  *) NULL;

		// initialize data structure for output triangulation
		out.pointlist               = (REAL *) NULL;
		out.pointattributelist      = (REAL *) NULL;
		out.pointmarkerlist         = (int  *) NULL;
		out.trianglelist            = (int  *) NULL;
		out.triangleattributelist   = (REAL *) NULL;
		out.neighborlist            = (int  *) NULL;
		out.segmentlist             = (int  *) NULL;
		out.segmentmarkerlist       = (int  *) NULL;
		out.edgelist                = (int  *) NULL;
		out.edgemarkerlist          = (int  *) NULL;

		// initialize data structure for output Voronoi diagram (unused)
		vorout.pointlist            = (REAL *) NULL;
		vorout.pointattributelist   = (REAL *) NULL;
		vorout.edgelist             = (int  *) NULL;
		vorout.normlist             = (REAL *) NULL;

		// Add points to in.pointlist

		for (int j = 0; j < nNodes; j++){

			in.pointlist[2*j+0] = dPanelCentroid[0][j];
			in.pointlist[2*j+1] = dPanelCentroid[1][j];

		}

		// Compute Delaunay triangulation.  Use the 'c' option so that
		// the convex hull is included

		char options[256] ="cjzenYQ";
		triangulate(options, &in, &out, &vorout);

		// Convert to mesh object by building face vector
		for (int j = 0; j < out.numberoftriangles; j++){

			Face newFace(3);

			newFace.SetNode(0, out.trianglelist[3*j+0]);
			newFace.SetNode(1, out.trianglelist[3*j+1]);
			newFace.SetNode(2, out.trianglelist[3*j+2]);
			vecMesh[i].faces.push_back(newFace);

		}

		// Build node vector.  Note that the Delaunay triangulation preserves
		// point indexing.
		for (int j = 0; j < out.numberofpoints; j++){

			Node node(out.pointlist[2*j+0],out.pointlist[2*j+1],0);
			vecMesh[i].nodes.push_back(node);

		}

		vecMesh[i].ConstructReverseNodeArray();
		vecMesh[i].RemoveCoincidentNodes();
		vecMesh[i].RemoveZeroEdges();
		vecMesh[i].ConstructEdgeMap();

		free(in.pointlist);
		free(out.pointlist);
		free(out.pointattributelist);
		free(out.pointmarkerlist);
		free(out.trianglelist);
		free(out.triangleattributelist);
		free(out.neighborlist);
		free(out.segmentlist);
		free(out.segmentmarkerlist);
		free(out.edgelist);
		free(out.edgemarkerlist);

	}

	// Construct kd-tree for each panel

	for (int i = 0; i < 6; i++){

		vecKDTreePanelI[i] = kd_create(3);

		for (int j = 0; j < vecMesh[i].faces.size(); j++){

				Face face = vecMesh[i].faces[j];

				Node nodeCenter = GetFaceCentroid(face, vecMesh[i].nodes);

				kd_insert3(
					vecKDTreePanelI[i],
					nodeCenter.x,
					nodeCenter.y,
					nodeCenter.z,
					(void*)(&(vecMesh[i].faces[j])));

		}

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
		
		if( nOverlapFaces == 0 ) continue;

		// Loop through all overlap faces associated with this source face
		for (int j = 0; j < nOverlapFaces; j++) {
			int iTargetFace = meshOverlap.vecTargetFaceIx[ixOverlap + j];

			// signal to not participate, because it is a ghost target
			if( iTargetFace < 0 ) continue;  // skip and do not do anything

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

			// Loop through all sub-triangles of this overlap Face
			for (int k = 0; k < nSubTriangles; k++) {

				// Loop through all quadrature nodes on this overlap Face
				for (int p = 0; p < dW.GetRows(); p++) {

					// Get quadrature node and pointwise Jacobian
					const Node & nodeQ = dQuadPtNodes(k,p);

					// Get lat/lon coordinates of nodeQ
					double dLatRad;
					double dLonRad;

					XYZtoRLL_Rad(nodeQ.x,nodeQ.y,nodeQ.z,dLonRad,dLatRad);

					// Determine the panel where nodeQ is located
					int iPanelIndex;

					for (int i = 0; i < 6; i++){

						if ((dPanelLat(i,0) <= dLatRad) && (dLatRad <= dPanelLat(i,1)) &&
							(dPanelLon(i,0) <= dLonRad) && (dLonRad <= dPanelLon(i,1))){

								iPanelIndex = i;
								break;

							}

					}
					// Compute Gnomonic projection onto the corresponding plane

					double dGX,dGY;
					double dLatTan;
					double dLonTan;

					if ((0 <= iPanelIndex) && (iPanelIndex <= 3)){

						dLatTan = 0;
						dLonTan = (iPanelIndex+1)*45*M_PI/180 + 45*iPanelIndex*M_PI/180;

					}
					else if (iPanelIndex == 4){

						dLatTan = M_PI/2;
						dLonTan = 0;

					}
					else{

						dLonTan = 0;
						dLatTan = -M_PI/2;

					}

					GnomonicProjection(dLonTan,dLatTan,dLonRad,dLatRad,dGX,dGY);

					// Get point closest to dGX and dGY
					kdres * kdresTarget =
						kd_nearest3(
							vecKDTreePanelI[iPanelIndex],
							dGX,
							dGY,
							0.0);

					Face * pFace = (Face *)(kd_res_item_data(kdresTarget));

					int iNearestFace = pFace - &(vecMesh[iPanelIndex].faces[0]);

					// Find triangle that contains the Gnomonic projection of nodeQ.
					// This is the face whose local index is iFaceFinal

					int iFaceFinal;

					double dA,dB;

					BarycentricCoordinates(vecMesh[iPanelIndex],iNearestFace,dGX,dGY,dA,dB);

					GetTriangleThatContainsPoint(vecMesh[iPanelIndex],iNearestFace,iFaceFinal,dGX,dGY);

					// Get global indices of the vertices of this triangle and
					// calculate corresponding centroids

					std::vector<int> vecContributingFaceI(3);
					DataArray2D<double> dataContributingCentroids(3,3);

					// Indices on the local mesh of the containing triangle

					int iLocalVertex1 = vecMesh[iPanelIndex].faces[iFaceFinal][0];
					int iLocalVertex2 = vecMesh[iPanelIndex].faces[iFaceFinal][1];
					int iLocalVertex3 = vecMesh[iPanelIndex].faces[iFaceFinal][2];

					// Indices on the source mesh of the containing triangle

					int iGlobalVertex1 = vecGlobalIndexI[iPanelIndex][iLocalVertex1];
					int iGlobalVertex2 = vecGlobalIndexI[iPanelIndex][iLocalVertex2];
					int iGlobalVertex3 = vecGlobalIndexI[iPanelIndex][iLocalVertex3];

					vecContributingFaceI[0] = iGlobalVertex1;
					vecContributingFaceI[1] = iGlobalVertex2;
					vecContributingFaceI[2] = iGlobalVertex3;

					// The centroids of the source mesh are the vertices of the containing triangle

					for (int i = 0; i < 3; i++){

						Face faceCurrent = meshInput.faces[vecContributingFaceI[i]];

						Node nodeCenter = GetFaceCentroid(faceCurrent,meshInput.nodes);

						nodeCenter = nodeCenter.Normalized();

						dataContributingCentroids(i,0) = nodeCenter.x;
						dataContributingCentroids(i,1) = nodeCenter.y;
						dataContributingCentroids(i,2) = nodeCenter.z;


					}

					// Vector of areas of subtriangles
					std::vector<double> vecSubAreas(3);

					// Area of triangle is obtained by adding up areas of the three
					// subtriangles that are formed by the sample point

					double dTriangleArea = 0;

					for (int i = 0; i < 3; i++){

						Face faceCurrent(3);

						faceCurrent.SetNode(0,0);
						faceCurrent.SetNode(1,1);
						faceCurrent.SetNode(2,2);

						NodeVector nodesCurrent;

						nodesCurrent.push_back(nodeQ);

						Node node1(dataContributingCentroids((i+1)%3,0),dataContributingCentroids((i+1)%3,1),
								   dataContributingCentroids((i+1)%3,2));

						Node node2(dataContributingCentroids((i+2)%3,0),dataContributingCentroids((i+2)%3,1),
								   dataContributingCentroids((i+2)%3,2));

						nodesCurrent.push_back(node1);
						nodesCurrent.push_back(node2);

						vecSubAreas[i] = CalculateFaceArea(faceCurrent,nodesCurrent);

						dTriangleArea += vecSubAreas[i];

					}

					// Calculate vector of weights
					std::vector<double> vecContributingFaceWeights(3);

					for (int i = 0; i < 3; i++){

						vecContributingFaceWeights[i] = vecSubAreas[i]/dTriangleArea;

					}

					// Contribution of this quadrature point to the map
					for (int i = 0; i < vecContributingFaceI.size(); i++){

						smatMap(iTargetFace, vecContributingFaceI[i]) +=
							vecContributingFaceWeights[i]
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

void LinearRemapIntegratedGeneralizedBarycentric(
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

	// Get SparseMatrix representation of the OfflineMap
	SparseMatrix<double> & smatMap = mapRemap.GetSparseMatrix();
	
	//If there are no holes on the source mesh, use the (global) dual of the source.		
	
	if (!DoesMeshHaveHoles(meshInput)) {
		
		Mesh meshInputDual = meshInput;
		
		//Construct dual mesh
		Dual(meshInputDual);
			
		//Construct edge map
		
		meshInputDual.ConstructEdgeMap();
		
		//kd-tree of dual mesh centers
	    kdtree * kdTarget = kd_create(3);
		for (int i = 0; i < meshInputDual.faces.size(); i++){

			const Face & face = meshInputDual.faces[i];

			Node nodeCentroid = GetFaceCentroid(face, meshInputDual.nodes);
		
			nodeCentroid = nodeCentroid.Normalized();

			kd_insert3(
				kdTarget,
				nodeCentroid.x,
				nodeCentroid.y,
				nodeCentroid.z,
				(void*)(&(meshInputDual.faces[i])));
			
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
	
			if( nOverlapFaces == 0 ) continue;
	
			// Loop through all overlap faces associated with this source face
			for (int j = 0; j < nOverlapFaces; j++) {
				
				int iTargetFace = meshOverlap.vecTargetFaceIx[ixOverlap + j];
	
				// signal to not participate, because it is a ghost target
				if( iTargetFace < 0 ) continue;  // skip and do not do anything
	
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
	
				// Loop through all sub-triangles of this overlap Face
				for (int k = 0; k < nSubTriangles; k++) {
	
					// Loop through all quadrature nodes on this overlap Face
					for (int p = 0; p < dW.GetRows(); p++) {
	
						// Get quadrature node and pointwise Jacobian
						Node & nodeQ = dQuadPtNodes(k,p);
	
						//Get the dual mesh face whose center is nearest to the sample point
	
						kdres * kdresTarget =
							kd_nearest3(
								kdTarget,
								nodeQ.x,
								nodeQ.y,
								nodeQ.z);
								
						
						Face * pFace = (Face *)(kd_res_item_data(kdresTarget));
		
						int iNearestFace = pFace - &(meshInputDual.faces[0]);
						
						// Find face that contains nodeQ
						// This is the face whose local index is iFaceFinal
						
						int iFaceFinal = 0;
						
						GetFaceThatContainsPoint(meshInputDual,iNearestFace,iFaceFinal,nodeQ.x,nodeQ.y,nodeQ.z);
						
						int iEdges = meshInputDual.faces[iFaceFinal].edges.size();
						
						//NodeVector nodesCurrentFace(iEdges);
						NodeVector nodesCurrentFace;
						
						for (int i = 0; i < iEdges; i++){
							
							nodesCurrentFace.push_back(meshInputDual.nodes[meshInputDual.faces[iFaceFinal][i]]);
							
						}
									
						std::vector<double> vecContributingFaceWeights;
					
						std::vector<int> vecContributingFaceI;		
							
						DataArray1D<double> dCoeffs(3);
						
						BilinearWeights(nodeQ,nodesCurrentFace,meshInputDual.faces[iFaceFinal],vecContributingFaceWeights,vecContributingFaceI);
												
						// Contribution of each point to the map
						for (int i = 0; i < vecContributingFaceI.size(); i++){
					
							if( vecContributingFaceWeights[i] < -1e-12 || vecContributingFaceWeights[i] > 1+1e-12 ){
								std::cout << "\nFound weight value = " << vecContributingFaceWeights[i] << std::endl;
								_EXCEPTIONT("Non-monotone weight");
							}
							
							smatMap(iTargetFace, vecContributingFaceI[i]) +=
								vecContributingFaceWeights[i]
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
	
	//If the source mesh has holes, use local dual faces
	
	else {
		
		//Vector of source face centers
				
		NodeVector vecSourceCenterArray;	
		
		for (int i = 0; i < meshInput.faces.size(); i++){
	
			const Face & face = meshInput.faces[i];
	
			Node nodeCentroid = GetFaceCentroid(face, meshInput.nodes);
			
			nodeCentroid = nodeCentroid.Normalized();
			
			vecSourceCenterArray.push_back(nodeCentroid);
	
				
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
			
			int iEdges = faceFirst.edges.size();
			
			//Nodes of the current face
			
			NodeVector nodesFaceFirst;
			
			for (int i = 0; i < faceFirst.edges.size(); i++){
				
				nodesFaceFirst.push_back(meshInput.nodes[faceFirst[i]]);
				
			}	
	
			// Find the set of Faces that overlap faceFirst
			int ixOverlapBegin = ixOverlap;
			
			int ixOverlapEnd = ixOverlapBegin;
	
			for (; ixOverlapEnd < meshOverlap.faces.size(); ixOverlapEnd++) {
				if (meshOverlap.vecSourceFaceIx[ixOverlapEnd] != ixFirst) {
					break;
				}
			}
			
			int nOverlapFaces = ixOverlapEnd - ixOverlapBegin;
	
			if( nOverlapFaces == 0 ) continue;
			
			std::vector<std::vector<int>> vecContributingFaceI;
			
			std::vector<std::vector<double>> vecContributingWeights;		
	
			// Loop through all overlap faces associated with this source face
			for (int j = 0; j < nOverlapFaces; j++) {
				
				int iTargetFace = meshOverlap.vecTargetFaceIx[ixOverlap + j];
	
				// signal to not participate, because it is a ghost target
				if( iTargetFace < 0 ) continue;  // skip and do not do anything
	
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
	
				// Loop through all sub-triangles of this overlap Face
				for (int k = 0; k < nSubTriangles; k++) {
	
					// Loop through all quadrature nodes on this overlap Face
					for (int p = 0; p < dW.GetRows(); p++) {
						
						vecContributingFaceI.clear();
						vecContributingWeights.clear();
						
						//Current quadrature node
						
						Node & nodeQ = dQuadPtNodes(k,p);
					
						std::vector<double> vecWeightsCurrentFace(iEdges);
						
						GeneralizedBarycentricCoordinates(nodeQ, nodesFaceFirst, vecWeightsCurrentFace);
						
						//loop through all nodes on source face
						
						for (int i = 0; i < iEdges; i++){
							
							//Global index of current node
								
							int iCurrentNode = meshInput.faces[ixFirst][i];
							
							//Number of nodes on local dual test face
							
							int iEdgesTestK = meshInput.revnodearray[iCurrentNode].size();
							
							if (iEdgesTestK < 3){
								
								std::vector<double> vecLocalWeights;
								std::vector<int> vecLocalContributingFaces;		
							
								for (auto it = meshInput.revnodearray[iCurrentNode].begin(); 
									  it != meshInput.revnodearray[iCurrentNode].end(); it++){
											
									vecLocalContributingFaces.push_back(*it);
									vecLocalWeights.push_back(vecWeightsCurrentFace[i]/iEdgesTestK);
								
								
								}
							
								vecContributingFaceI.push_back(vecLocalContributingFaces);
								vecContributingWeights.push_back(vecLocalWeights);
								
							
							}
							
							//Local dual face has at least 3 edges
							
							else{
														
								//Construct local dual face
								
								Face faceLocalDual(iEdgesTestK);
								
								NodeVector nodesLocalDual;
														
								ConstructLocalDualFace(meshInput,vecSourceCenterArray,iCurrentNode,faceLocalDual,nodesLocalDual);
												
								//Put the dual face nodes in a ccw oriented vector
						
								NodeVector nodesLocalDualReordered;
						
								for (int l = 0; l < nodesLocalDual.size(); l++){
									
									nodesLocalDualReordered.push_back(vecSourceCenterArray[faceLocalDual[l]]);
									
								}
								
								std::vector<double> vecLocalWeights(iEdgesTestK);
								std::vector<int> vecLocalContributingFaces(iEdgesTestK);
								
								//If the local dual face around the current node contains the target, then quit
								
								if (DoesFaceContainPoint(nodesLocalDualReordered, nodeQ.x, nodeQ.y, nodeQ.z)) {
									
									vecContributingWeights.clear();
									vecContributingFaceI.clear();
										
									GeneralizedBarycentricCoordinates(nodeQ, nodesLocalDualReordered, vecLocalWeights);
									
									for ( int m = 0; m < iEdgesTestK; m++){
									
										vecLocalContributingFaces[m] = faceLocalDual[m];
													
									}
									
									vecContributingFaceI.push_back(vecLocalContributingFaces);
									vecContributingWeights.push_back(vecLocalWeights);
									
									break;
										
								}
								
								else {
									
									//Use simple average or barycentric coorindates depending on whether the dual face 
									//contains the current node or not
									
									Node nodeCurrentNode = meshInput.nodes[iCurrentNode];
									
									if (DoesFaceContainPoint(nodesLocalDualReordered, nodeCurrentNode.x, nodeCurrentNode.y, 
															nodeCurrentNode.z)){
										
										//Compute generalized barycentric coordinate of current source face node wrt the local dual face
		
										GeneralizedBarycentricCoordinates(nodeCurrentNode, nodesLocalDualReordered, vecLocalWeights);
			
										for ( int m = 0; m < iEdgesTestK; m++ ){
									
											vecLocalContributingFaces[m] = faceLocalDual[m];
											vecLocalWeights[m] *= vecWeightsCurrentFace[i];
								
										}
										
										vecContributingFaceI.push_back(vecLocalContributingFaces);
										vecContributingWeights.push_back(vecLocalWeights);
										
									}
									
									
									else {
										
										//If the local dual face doesn't contain the current node, use averaging
										
										for ( int m = 0; m < iEdgesTestK; m++ ){
											
											vecLocalContributingFaces[m] = (faceLocalDual[m]);
											
											vecLocalWeights[m] = vecWeightsCurrentFace[i]*(1.0/iEdgesTestK);
											
										}
										
										vecContributingFaceI.push_back(vecLocalContributingFaces);
										vecContributingWeights.push_back(vecLocalWeights);
										
									}
									
								}
							
							
							}
										
						}
						
						//Add contributions to map
						
						for (int m = 0; m < vecContributingWeights.size(); m++){
							
							for (int l = 0; l < vecContributingWeights[m].size(); l++){
														
								int iContributingFace = vecContributingFaceI[m][l];
								
								smatMap(iTargetFace,iContributingFace) += vecContributingWeights[m][l]
											* dQuadPtWeight(k,p)
											* meshOverlap.vecFaceArea[ixOverlap + j]
											/ dQuadratureArea
											/ meshOutput.vecFaceArea[iTargetFace];
								
							}
							
						}
						
					}
					
				}
				
			}
			 
			// Increment the current overlap index
			ixOverlap += nOverlapFaces;
			
		}
		
	}

}

///////////////////////////////////////////////////////////////////////////////

void LinearRemapGeneralizedBarycentric(
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

	// Get SparseMatrix representation of the OfflineMap
	SparseMatrix<double> & smatMap = mapRemap.GetSparseMatrix();
	
	//If there are no holes on the source mesh, use the (global) dual of the source.
	
	if(!DoesMeshHaveHoles(meshInput)){
	
		//Dual mesh
		Mesh meshInputDual = meshInput;
		
		//Construct dual mesh
		Dual(meshInputDual);
		
		//Reverse node array
		meshInputDual.ConstructReverseNodeArray();	
		
		//Construct edge map
		meshInputDual.ConstructEdgeMap();
		
		//kd-tree of dual mesh centers
		kdtree * kdDual = kd_create(3);
		
		// Vector of centers of the source mesh
		for (int i = 0; i < meshInputDual.faces.size(); i++){

			const Face & face = meshInputDual.faces[i];

			Node nodeCentroid = GetFaceCentroid(face, meshInputDual.nodes);
		
			nodeCentroid = nodeCentroid.Normalized();
		
			kd_insert3(
				kdDual,
				nodeCentroid.x,
				nodeCentroid.y,
				nodeCentroid.z,
				(void*)(&(meshInputDual.faces[i])));
		
			
		}
		
		// Loop through all target faces
		for (int ixFirst = 0; ixFirst < meshOutput.faces.size(); ixFirst++) {

			// Output every 1000 overlap elements
			if (ixFirst % 1000 == 0) {
				Announce("Element %i/%i", ixFirst, meshOutput.faces.size());
			}
	
			// This Face
			const Face & faceFirst = meshOutput.faces[ixFirst];
	
			// Get node coordinates of each target face center
			Node nodeQ = GetFaceCentroid(faceFirst,meshOutput.nodes);

			nodeQ = nodeQ.Normalized();

			//Get the dual mesh face whose center is nearest to the target face

			kdres * kdresTarget =
				kd_nearest3(
					kdDual,
					nodeQ.x,
					nodeQ.y,
					nodeQ.z);

			Face * pFace = (Face *)(kd_res_item_data(kdresTarget));

			int iNearestFace = pFace - &(meshInputDual.faces[0]);

			// Find triangle that contains nodeQ.
			// This is the face whose local index is iFaceFinal

			int iFaceFinal = 0;

			GetFaceThatContainsPoint(meshInputDual,iNearestFace,iFaceFinal,nodeQ.x,nodeQ.y,nodeQ.z);
			
			int iEdges = meshInputDual.faces[iFaceFinal].edges.size();
			
			//Current dual mesh face
			
			Face faceDualFace = meshInputDual.faces[iFaceFinal];
			
			//Vector of nodes on current dual face
			
			NodeVector nodesDualFace;
			
			for (int i = 0; i < iEdges; i++){
				
				nodesDualFace.push_back(meshInputDual.nodes[faceDualFace[i]]);
				
			}
			
			//Vector of contributing weights
			
			std::vector<double> vecWeights(iEdges);
			
			//Compute generalized barycentric coordinates
			
			GeneralizedBarycentricCoordinates(nodeQ, nodesDualFace, vecWeights);	
		
			//Add weights to map
			
			for (int i = 0; i < iEdges; i++){

				int iContributingFaceI = meshInputDual.faces[iFaceFinal][i];

				smatMap(ixFirst, iContributingFaceI) = vecWeights[i];
				
			}
			
		}
		
	}
	
	//If the source mesh has holes, use local dual faces
	
	else{

		//kd-tree of dual mesh centers
		kdtree * kdTarget = kd_create(3);
    
		//Vector of source face centers
		NodeVector vecSourceCenterArray;
    
		//Vector of target face centers
		NodeVector vecTargetCenterArray;

		// Vector of centers of the source mesh
		for (int i = 0; i < meshInput.faces.size(); i++){
	
			const Face & face = meshInput.faces[i];
	
			Node nodeCentroid = GetFaceCentroid(face, meshInput.nodes);
			
			nodeCentroid = nodeCentroid.Normalized();
			
			vecSourceCenterArray.push_back(nodeCentroid);
	
				
		}
		
		for (int i = 0; i < meshOutput.faces.size(); i++){

			const Face & face = meshOutput.faces[i];

			Node nodeCentroid = GetFaceCentroid(face, meshOutput.nodes);
		
			nodeCentroid = nodeCentroid.Normalized();

			kd_insert3(
				kdTarget,
				nodeCentroid.x,
				nodeCentroid.y,
				nodeCentroid.z,
				(void*)(&(meshOutput.faces[i])));
			
			vecTargetCenterArray.push_back(nodeCentroid);
			
		}
		
		std::set<int> setAllTargets;
		
		// Loop through all source faces
		
		for (int ixFirst = 0; ixFirst < meshInput.faces.size(); ixFirst++) {
	
			// Output every 1000 source elements
			if (ixFirst % 1000 == 0) {
				
				Announce("Element %i/%i", ixFirst, meshInput.faces.size());
				
			}
			
			// This Face
			const Face & faceFirst = meshInput.faces[ixFirst];
			
			//Find the maximum distance from the current source face center to its nodes
			
			double dMaxNodeDist = 0.0;
			
			Node faceFirstCenter = vecSourceCenterArray[ixFirst];
			
			int iEdges = faceFirst.edges.size();
			
			for (int i = 0; i < iEdges; i++){
				
				//Great circle distance is just the arccos of dot product
				
				Node nodeI = meshInput.nodes[faceFirst[i]];
				
				double dDistToNodeI = acos(DotProduct(nodeI, faceFirstCenter));
				
				if ( dDistToNodeI > dMaxNodeDist ){
					
					dMaxNodeDist = dDistToNodeI;
					
				}
				
			}
			
			//Find all the target face centers that are within a given distance of the given source face center
			
			/* find points closest to the origin and within distance radius */
			//struct *presults;
										 
			struct kdres *presults = kd_nearest_range3(kdTarget,faceFirstCenter.x, faceFirstCenter.y, 
													   faceFirstCenter.z, dMaxNodeDist+.0001);
			
										   
			if( kd_res_size(presults) == 0 ){
				
				continue;
				
			}										   
				   
			//Nodes of the current face
			
			NodeVector nodesFaceFirst;
			
			for (int i = 0; i < faceFirst.edges.size(); i++){
				
				nodesFaceFirst.push_back(meshInput.nodes[faceFirst[i]]);
				
			}
						
			std::vector<int> vecTargetsInFace;
			
			while( !kd_res_end( presults ) ) {
				
				
			    /* get the data and position of the current result item */
			    Face * pFace = (Face *)kd_res_item_data( presults );
			    
			    int iNearestTarget = pFace - &(meshOutput.faces[0]);
			    			    
			    if(DoesFaceContainPoint(nodesFaceFirst, vecTargetCenterArray[iNearestTarget].x,
										vecTargetCenterArray[iNearestTarget].y,vecTargetCenterArray[iNearestTarget].z)) {
											
					
					if( setAllTargets.find(iNearestTarget) == setAllTargets.end() ){
						
						setAllTargets.insert(iNearestTarget);
						vecTargetsInFace.push_back(iNearestTarget);
						
					}
						
				}
								
			    /* go to the next entry */
			    kd_res_next( presults );
			    
			}
			
			
			kd_res_free( presults );
						
			std::vector<std::vector<int>> vecContributingFaceI;
			
			std::vector<std::vector<double>> vecContributingWeights;
												
				
			//Loop through all centers of the target faces in the given source face
			for (int i = 0; i < vecTargetsInFace.size(); i++){
				
				vecContributingFaceI.clear();
				vecContributingWeights.clear();
				
				//Current target face center
				
				Node nodeQ = vecTargetCenterArray[vecTargetsInFace[i]];
			
				std::vector<double> vecWeightsCurrentFace(iEdges);
				
				GeneralizedBarycentricCoordinates(nodeQ, nodesFaceFirst, vecWeightsCurrentFace);
				
				//loop through all nodes on source face
				
				for (int k = 0; k < iEdges; k++){
					
					//Global index of current node
						
					int iCurrentNode = meshInput.faces[ixFirst][k];
					
					//Number of nodes on local dual test face
					
					int iEdgesTestK = meshInput.revnodearray[iCurrentNode].size();
					
					if (iEdgesTestK < 3){
								
						std::vector<double> vecLocalWeights;
						std::vector<int> vecLocalContributingFaces;		
					
						for (auto it = meshInput.revnodearray[iCurrentNode].begin(); 
							  it != meshInput.revnodearray[iCurrentNode].end(); it++){
									
							vecLocalContributingFaces.push_back(*it);
							vecLocalWeights.push_back(vecWeightsCurrentFace[k]/iEdgesTestK);
									
						}
					
						vecContributingFaceI.push_back(vecLocalContributingFaces);
						vecContributingWeights.push_back(vecLocalWeights);
					
					}
					
					//Local dual face has at least three edges
					
					else{
												
						//Construct local dual face
						
						Face faceLocalDual(iEdgesTestK);
						
						NodeVector nodesLocalDual;
												
						ConstructLocalDualFace(meshInput,vecSourceCenterArray,iCurrentNode,faceLocalDual,nodesLocalDual);
										
						//Put the dual face nodes in an ccw oriented vector
				
						NodeVector nodesLocalDualReordered;
				
						for (int j = 0; j < nodesLocalDual.size(); j++){
							
							nodesLocalDualReordered.push_back(vecSourceCenterArray[faceLocalDual[j]]);
							
						}
						
						std::vector<double> vecLocalWeights(iEdgesTestK);
						std::vector<int> vecLocalContributingFaces(iEdgesTestK);
						
						
						//If the local dual face around the current node contains the target, then quit
						
						if (DoesFaceContainPoint(nodesLocalDualReordered, nodeQ.x, nodeQ.y, nodeQ.z)) {
							
								vecContributingWeights.clear();
								vecContributingFaceI.clear();
						
								GeneralizedBarycentricCoordinates(nodeQ, nodesLocalDualReordered, vecLocalWeights);
						
								for ( int m = 0; m < iEdgesTestK; m++){
							
									vecLocalContributingFaces[m] = faceLocalDual[m];
									
								}			
							
								vecContributingFaceI.push_back(vecLocalContributingFaces);
								vecContributingWeights.push_back(vecLocalWeights);
								
								break;
								
						}
						
						else {
							
							//Use simple average or barycentric coorindates depending on whether the dual face 
							//contains the current node or not
							
							Node nodeCurrentNode = meshInput.nodes[iCurrentNode];
							
							if (DoesFaceContainPoint(nodesLocalDualReordered, nodeCurrentNode.x, nodeCurrentNode.y, 
													nodeCurrentNode.z)){
								
								//Compute generalized barycentric coordinate of current source face node wrt the local dual face
							
								GeneralizedBarycentricCoordinates(nodeCurrentNode, nodesLocalDualReordered, vecLocalWeights);
								
												
								for ( int m = 0; m < iEdgesTestK; m++ ){
							
									vecLocalContributingFaces[m] = faceLocalDual[m];
									
									vecLocalWeights[m] *= vecWeightsCurrentFace[k];
						
								}
								
								vecContributingFaceI.push_back(vecLocalContributingFaces);
								vecContributingWeights.push_back(vecLocalWeights);
								
							}
							
							
							else {
								
								for ( int m = 0; m < iEdgesTestK; m++ ){
									
									vecLocalContributingFaces[m] = (faceLocalDual[m]);
									
									vecLocalWeights[m] = vecWeightsCurrentFace[k]*(1.0/iEdgesTestK);
									
								}
								
								vecContributingFaceI.push_back(vecLocalContributingFaces);
								vecContributingWeights.push_back(vecLocalWeights);
								
							}
							
						}
					
					
					}
								
				}
				
				//Add contributions to map
				
				for (int m = 0; m < vecContributingWeights.size(); m++){
					
					for (int j = 0; j < vecContributingWeights[m].size(); j++){
												
						int iContributingFace = vecContributingFaceI[m][j];
						
						smatMap(vecTargetsInFace[i],iContributingFace) += vecContributingWeights[m][j];
						
					}
					
				}
					
			}
			
		}
	}
	
}

///////////////////////////////////////////////////////////////////////////////

void LinearRemapTriangulation(
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

	// Array of global indices
	std::vector <std::vector<int>> vecGlobalIndexI(6);

	// Vector storing meshes of each panel
	std::vector<Mesh> vecMesh(6);

	// kd-tree for each panel
	std::vector<kdtree *> vecKDTreePanelI(6);

	// Arrays of panel boundaries in lat/lon coordinates (radians)
	DataArray2D<double> dPanelLat(6,2);
	DataArray2D<double> dPanelLon(6,2);

	// Width of buffer zone for each panel (degrees)
	double dBuff = 30;

	// Define equatorial panel boundaries
	for (int i = 0; i < 4; i++){

		dPanelLon(i,0) = 90*M_PI*i/180;  //Left boundary
		dPanelLon(i,1) = 90*M_PI*(i+1)/180;  //Right Boundary
		dPanelLat(i,1) = 45*M_PI/180;  //Upper Boundary
		dPanelLat(i,0) = -45*M_PI/180;  //Lower Boundary

	}

	// Define polar panels

	dPanelLon(4,0) = 0; dPanelLon(4,1) = 2*M_PI;
	dPanelLat(4,0) = 45*M_PI/180; dPanelLat(4,1) = M_PI/2;

	dPanelLon(5,0) = 0; dPanelLon(5,1) = 2*M_PI;
	dPanelLat(5,1) = -45*M_PI/180; dPanelLat(5,0) = -M_PI/2;


	// Generate the Delaunay triangulation for each of the 6 panels
	for (int i = 0; i < 6; i++){

		// Array storing the coordinates of the face centroids of panel i
		std::vector<std::vector<double>> dPanelCentroid(2);

		// Coordinates of tangent plane point of tangency for panel i
		double dLatTan;
		double dLonTan;

		if ((0 <= i) && (i <= 3)){

			dLatTan = 0;
			dLonTan = (i+1)*45*M_PI/180 + 45*i*M_PI/180;

		}
		else if (i == 4){

			dLatTan = M_PI/2;
			dLonTan = 0;

		}
		else{

			dLonTan = 0;
			dLatTan = -M_PI/2;

		}

		// Number of source faces centers for panel i; this is the size
		// of in.pointlist
		int nNodes = 0;

		for (int j = 0; j < meshInput.faces.size(); j++){

			// Get coordinates of centroid of current face

			double dLonRad0;
			double dLatRad0;

			const Face & face = meshInput.faces[j];

			//const Node & node = GetFaceCentroid(face,meshInput.nodes);
			Node node = GetFaceCentroid(face,meshInput.nodes);

			node = node.Normalized();

			XYZtoRLL_Rad(node.x,node.y,node.z,dLonRad0,dLatRad0);

			// Determine if centroid is in panel i plus buffer zone
			bool fPanelContainsCentroid = false;

			if ((0 <= i) && (i <= 3)){

				if ((dPanelLat(i,0) - dBuff*M_PI/180 <= dLatRad0) && (dLatRad0 <= dPanelLat(i,1)+dBuff*M_PI/180)){

					if ( i == 0 ){

						if ( ((2*M_PI - dBuff*M_PI/180 <= dLonRad0) && (dLonRad0 <= 2*M_PI)) ||
							 ((dPanelLon(i,0) <= dLonRad0) && (dLonRad0 <= dPanelLon(i,1) + dBuff*M_PI/180)) ) {

							fPanelContainsCentroid = true;
						}
					}

					else if ( i == 3 ){

						if ( ((0 <= dLonRad0) && (dLonRad0 <= 0 + dBuff*M_PI/180)) ||
							 ((dPanelLon(i,0) - dBuff*M_PI/180 <= dLonRad0) && (dLonRad0 <= dPanelLon(i,1) + dBuff*M_PI/180)) ){

							fPanelContainsCentroid = true;

						}

					}

					else {

						if ( (dPanelLon(i,0) - dBuff*M_PI/180 <= dLonRad0) && (dLonRad0 <= dPanelLon(i,1) + dBuff*M_PI/180)){

							fPanelContainsCentroid = true;
						}


					}

				}
			}
			else if (i == 4){

				if ((dPanelLat(i,0) - dBuff*M_PI/180 <= dLatRad0) && (dLatRad0 <= dPanelLat(i,1))){

					fPanelContainsCentroid = true;

				}

			}

			else {

				if ((dPanelLat(i,0) <= dLatRad0) && (dLatRad0 <= dPanelLat(i,1) + dBuff*M_PI/180)){
					fPanelContainsCentroid = true;

				}

			}

			if(fPanelContainsCentroid){

				vecGlobalIndexI[i].push_back(j);
				nNodes++;

				// Gnomonic projection of centroid coordinates
				double dGX;
				double dGY;
				GnomonicProjection(dLonTan,dLatTan,dLonRad0,dLatRad0,dGX,dGY);

				// Add Gnomonic coordinates to vector
				dPanelCentroid[0].push_back(dGX);
				dPanelCentroid[1].push_back(dGY);

			}
		}

		// Structures for Delaunay triangulation

		struct triangulateio in, out, vorout;

		in.numberofpoints           = nNodes;
		in.numberofpointattributes  = 0;
		in.numberofsegments         = 0;
		in.numberofholes            = 0;
		in.numberofregions          = 0;
		in.pointlist                = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));
		in.segmentlist              = (int  *) malloc(in.numberofsegments * 2 * sizeof(int));;
		in.pointattributelist       = (REAL *) NULL;
		in.pointmarkerlist          = (int  *) NULL;
		in.trianglelist             = (int  *) NULL;
		in.triangleattributelist    = (REAL *) NULL;
		in.neighborlist             = (int  *) NULL;
		in.segmentmarkerlist        = (int  *) NULL;
		in.edgelist                 = (int  *) NULL;
		in.edgemarkerlist           = (int  *) NULL;

		// initialize data structure for output triangulation
		out.pointlist               = (REAL *) NULL;
		out.pointattributelist      = (REAL *) NULL;
		out.pointmarkerlist         = (int  *) NULL;
		out.trianglelist            = (int  *) NULL;
		out.triangleattributelist   = (REAL *) NULL;
		out.neighborlist            = (int  *) NULL;
		out.segmentlist             = (int  *) NULL;
		out.segmentmarkerlist       = (int  *) NULL;
		out.edgelist                = (int  *) NULL;
		out.edgemarkerlist          = (int  *) NULL;

		// initialize data structure for output Voronoi diagram (unused)
		vorout.pointlist            = (REAL *) NULL;
		vorout.pointattributelist   = (REAL *) NULL;
		vorout.edgelist             = (int  *) NULL;
		vorout.normlist             = (REAL *) NULL;

		// Add points to in.pointlist

		for (int j = 0; j < nNodes; j++){

			in.pointlist[2*j+0] = dPanelCentroid[0][j];
			in.pointlist[2*j+1] = dPanelCentroid[1][j];

		}

		// Compute Delaunay triangulation.  Use the 'c' option so that
		// the convex hull is included

		char options[256] ="cjzenYQ";
		triangulate(options, &in, &out, &vorout);

		// Convert to mesh object by building face vector
		for (int j = 0; j < out.numberoftriangles; j++){

			Face newFace(3);

			newFace.SetNode(0, out.trianglelist[3*j+0]);
			newFace.SetNode(1, out.trianglelist[3*j+1]);
			newFace.SetNode(2, out.trianglelist[3*j+2]);
			vecMesh[i].faces.push_back(newFace);

		}

		// Build node vector.  Note that the Delaunay triangulation preserves
		// point indexing.
		for (int j = 0; j < out.numberofpoints; j++){

			Node node(out.pointlist[2*j+0],out.pointlist[2*j+1],0);
			vecMesh[i].nodes.push_back(node);

		}

		vecMesh[i].ConstructReverseNodeArray();
		vecMesh[i].RemoveCoincidentNodes();
		vecMesh[i].RemoveZeroEdges();
		vecMesh[i].ConstructEdgeMap();

		free(in.pointlist);
		free(out.pointlist);
		free(out.pointattributelist);
		free(out.pointmarkerlist);
		free(out.trianglelist);
		free(out.triangleattributelist);
		free(out.neighborlist);
		free(out.segmentlist);
		free(out.segmentmarkerlist);
		free(out.edgelist);
		free(out.edgemarkerlist);

	}

	// Construct kd-tree for each panel

	for (int i = 0; i < 6; i++){

		vecKDTreePanelI[i] = kd_create(3);

		for (int j = 0; j < vecMesh[i].faces.size(); j++){

				Face face = vecMesh[i].faces[j];

				Node nodeCenter = GetFaceCentroid(face, vecMesh[i].nodes);

				kd_insert3(
					vecKDTreePanelI[i],
					nodeCenter.x,
					nodeCenter.y,
					nodeCenter.z,
					(void*)(&(vecMesh[i].faces[j])));

		}

	}

	// Loop through all target faces
	for (int ixFirst = 0; ixFirst < meshOutput.faces.size(); ixFirst++) {

		// Output every 1000 overlap elements
		if (ixFirst % 1000 == 0) {
			Announce("Element %i/%i", ixFirst, meshOutput.faces.size());
		}

		// This Face
		const Face & faceFirst = meshOutput.faces[ixFirst];

			// Get node coordinates of each target face center
			Node nodeQ = GetFaceCentroid(faceFirst,meshOutput.nodes);

			nodeQ = nodeQ.Normalized();

			// Get lat/lon coordinates of nodeQ
			double dLatRad;
			double dLonRad;

			XYZtoRLL_Rad(nodeQ.x,nodeQ.y,nodeQ.z,dLonRad,dLatRad);

			// Determine the panel where nodeQ is located
			int iPanelIndex;

			for (int i = 0; i < 6; i++){

				if ((dPanelLat(i,0) <= dLatRad) && (dLatRad <= dPanelLat(i,1)) &&
					(dPanelLon(i,0) <= dLonRad) && (dLonRad <= dPanelLon(i,1))){

						iPanelIndex = i;
						break;

					}

			}
			// Compute Gnomonic projection onto the corresponding plane

			double dGX,dGY;
			double dLatTan;
			double dLonTan;

			if ((0 <= iPanelIndex) && (iPanelIndex <= 3)){

				dLatTan = 0;
				dLonTan = (iPanelIndex+1)*45*M_PI/180 + 45*iPanelIndex*M_PI/180;

			}
			else if (iPanelIndex == 4){

				dLatTan = M_PI/2;
				dLonTan = 0;

			}
			else{

				dLonTan = 0;
				dLatTan = -M_PI/2;

			}

			GnomonicProjection(dLonTan,dLatTan,dLonRad,dLatRad,dGX,dGY);

			// Get point closest to dGX and dGY
			kdres * kdresTarget =
				kd_nearest3(
					vecKDTreePanelI[iPanelIndex],
					dGX,
					dGY,
					0.0);

			Face * pFace = (Face *)(kd_res_item_data(kdresTarget));

			int iNearestFace = pFace - &(vecMesh[iPanelIndex].faces[0]);

			// Find triangle that contains the Gnomonic projection nodeQ.
			// This is the face whose local index is iFaceFinal

			int iFaceFinal;

			double dA,dB;

			BarycentricCoordinates(vecMesh[iPanelIndex],iNearestFace,dGX,dGY,dA,dB);

			GetTriangleThatContainsPoint(vecMesh[iPanelIndex],iNearestFace,iFaceFinal,dGX,dGY);

			// Get global indices of the vertices of this triangle and
			// calculate corresponding centroids

			std::vector<int> vecContributingFaceI(3);
			DataArray2D<double> dataContributingCentroids(3,3);

			// Indices on the local mesh of the containing triangle

			int iLocalVertex1 = vecMesh[iPanelIndex].faces[iFaceFinal][0];
			int iLocalVertex2 = vecMesh[iPanelIndex].faces[iFaceFinal][1];
			int iLocalVertex3 = vecMesh[iPanelIndex].faces[iFaceFinal][2];

			// Indices on the source mesh of the containing triangle

			int iGlobalVertex1 = vecGlobalIndexI[iPanelIndex][iLocalVertex1];
			int iGlobalVertex2 = vecGlobalIndexI[iPanelIndex][iLocalVertex2];
			int iGlobalVertex3 = vecGlobalIndexI[iPanelIndex][iLocalVertex3];

			vecContributingFaceI[0] = iGlobalVertex1;
			vecContributingFaceI[1] = iGlobalVertex2;
			vecContributingFaceI[2] = iGlobalVertex3;

			// The centroids of the source mesh are the vertices of the containing triangle

			for (int i = 0; i < 3; i++){

				Face faceCurrent = meshInput.faces[vecContributingFaceI[i]];

				Node nodeCenter = GetFaceCentroid(faceCurrent,meshInput.nodes);

				nodeCenter = nodeCenter.Normalized();

				dataContributingCentroids(i,0) = nodeCenter.x;
				dataContributingCentroids(i,1) = nodeCenter.y;
				dataContributingCentroids(i,2) = nodeCenter.z;


			}

			// Vector of areas of subtriangles
			std::vector<double> vecSubAreas(3);

			// Area of triangle is obtained by adding up areas of the three
			// subtriangles that are formed by the sample point

			double dTriangleArea = 0;

			for (int i = 0; i < 3; i++){

				Face faceCurrent(3);

				faceCurrent.SetNode(0,0);
				faceCurrent.SetNode(1,1);
				faceCurrent.SetNode(2,2);

				NodeVector nodesCurrent;

				nodesCurrent.push_back(nodeQ);

				Node node1(dataContributingCentroids((i+1)%3,0),dataContributingCentroids((i+1)%3,1),
						   dataContributingCentroids((i+1)%3,2));

				Node node2(dataContributingCentroids((i+2)%3,0),dataContributingCentroids((i+2)%3,1),
						   dataContributingCentroids((i+2)%3,2));

				nodesCurrent.push_back(node1);
				nodesCurrent.push_back(node2);

				vecSubAreas[i] = CalculateFaceArea(faceCurrent,nodesCurrent);

				dTriangleArea += vecSubAreas[i];

			}

			// Calculate vector of weights
			std::vector<double> vecContributingFaceWeights(3);

			for (int i = 0; i < 3; i++){

				vecContributingFaceWeights[i] = vecSubAreas[i]/dTriangleArea;

			}

			// Contribution of this point to the map
			for (int i = 0; i < vecContributingFaceI.size(); i++){

				smatMap(ixFirst, vecContributingFaceI[i]) = vecContributingFaceWeights[i];

			}
		}


}

///////////////////////////////////////////////////////////////////////////////

void LinearRemapIntegratedNearestNeighbor(
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

	// Get SparseMatrix representation of the OfflineMap
	SparseMatrix<double> & smatMap = mapRemap.GetSparseMatrix();
	
	//kd-tree of source mesh centers
	kdtree * kdSource = kd_create(3);
	
	for (int i = 0; i < meshInput.faces.size(); i++){

		const Face & face = meshInput.faces[i];

		Node nodeCentroid = GetFaceCentroid(face, meshInput.nodes);
	
		nodeCentroid = nodeCentroid.Normalized();

		kd_insert3(
			kdSource,
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
		
		int iEdges = faceFirst.edges.size();
		
		//Nodes of the current face
		
		NodeVector nodesFaceFirst;
		
		for (int i = 0; i < faceFirst.edges.size(); i++){
			
			nodesFaceFirst.push_back(meshInput.nodes[faceFirst[i]]);
			
		}	

		// Find the set of Faces that overlap faceFirst
		int ixOverlapBegin = ixOverlap;
		
		int ixOverlapEnd = ixOverlapBegin;

		for (; ixOverlapEnd < meshOverlap.faces.size(); ixOverlapEnd++) {
			if (meshOverlap.vecSourceFaceIx[ixOverlapEnd] != ixFirst) {
				break;
			}
		}
		
		int nOverlapFaces = ixOverlapEnd - ixOverlapBegin;

		if( nOverlapFaces == 0 ) continue;

		// Loop through all overlap faces associated with this source face
		for (int j = 0; j < nOverlapFaces; j++) {
			
			//int iTargetFace = meshOverlap.vecTargetFaceIx[ixFirst];
			int iTargetFace = meshOverlap.vecTargetFaceIx[ixOverlap + j];

			// signal to not participate, because it is a ghost target
			if( iTargetFace < 0 ) continue;  // skip and do not do anything

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

			// Loop through all sub-triangles of this overlap Face
			for (int k = 0; k < nSubTriangles; k++) {

				// Loop through all quadrature nodes on this overlap Face
				for (int p = 0; p < dW.GetRows(); p++) {

					// Get quadrature node
					Node & nodeQ = dQuadPtNodes(k,p);
					
					// Find nearest source mesh face
					kdres * kdresSource =
						kd_nearest3(
							kdSource,
							nodeQ.x,
							nodeQ.y,
							nodeQ.z);

					Face * pFace = (Face *)(kd_res_item_data(kdresSource));

					int iNearestFace = pFace - &(meshInput.faces[0]);		
			
					smatMap(iTargetFace, iNearestFace) +=
						dQuadPtWeight(k,p)
						* meshOverlap.vecFaceArea[ixOverlap + j]
						/ dQuadratureArea
						/ meshOutput.vecFaceArea[iTargetFace];	
					
				}
				
			}
			
		}
		 
		// Increment the current overlap index
		ixOverlap += nOverlapFaces;
		
	}
	
}

///////////////////////////////////////////////////////////////////////////////

void LinearRemapIntegratedBilinear(
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

	// Get SparseMatrix representation of the OfflineMap
	SparseMatrix<double> & smatMap = mapRemap.GetSparseMatrix();
	
	//If there are no holes on the source mesh, use the (global) dual of the source.	
	
	if (!DoesMeshHaveHoles(meshInput)) {
		
		Mesh meshInputDual = meshInput;
		
		//Construct dual mesh
		Dual(meshInputDual);
			
		//Construct edge map
		
		meshInputDual.ConstructEdgeMap();
		
		//kd-tree of dual mesh centers
	    kdtree * kdTarget = kd_create(3);
		for (int i = 0; i < meshInputDual.faces.size(); i++){

			const Face & face = meshInputDual.faces[i];

			Node nodeCentroid = GetFaceCentroid(face, meshInputDual.nodes);
		
			nodeCentroid = nodeCentroid.Normalized();

			kd_insert3(
				kdTarget,
				nodeCentroid.x,
				nodeCentroid.y,
				nodeCentroid.z,
				(void*)(&(meshInputDual.faces[i])));
			
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
	
			if( nOverlapFaces == 0 ){ 
				
				continue;
				
			}
	
			// Loop through all overlap faces associated with this source face
			for (int j = 0; j < nOverlapFaces; j++) {
				
				int iTargetFace = meshOverlap.vecTargetFaceIx[ixOverlap + j];
	
				// signal to not participate, because it is a ghost target
				if( iTargetFace < 0 ) continue;  // skip and do not do anything
	
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
	
				// Loop through all sub-triangles of this overlap Face
				for (int k = 0; k < nSubTriangles; k++) {
	
					// Loop through all quadrature nodes on this overlap Face
					for (int p = 0; p < dW.GetRows(); p++) {
	
						// Get quadrature node and pointwise Jacobian
						Node & nodeQ = dQuadPtNodes(k,p);
	
						//Get the dual mesh face whose center is nearest to the sample point
	
						kdres * kdresTarget =
							kd_nearest3(
								kdTarget,
								nodeQ.x,
								nodeQ.y,
								nodeQ.z);
								
						
						Face * pFace = (Face *)(kd_res_item_data(kdresTarget));
		
						int iNearestFace = pFace - &(meshInputDual.faces[0]);
						
						// Find face that contains nodeQ
						// This is the face whose local index is iFaceFinal
						
						int iFaceFinal = 0;
						
						GetFaceThatContainsPoint(meshInputDual,iNearestFace,iFaceFinal,nodeQ.x,nodeQ.y,nodeQ.z);
						
						int iEdges = meshInputDual.faces[iFaceFinal].edges.size();
						
						//NodeVector nodesCurrentFace(iEdges);
						NodeVector nodesCurrentFace;
						
						for (int i = 0; i < iEdges; i++){
							
							nodesCurrentFace.push_back(meshInputDual.nodes[meshInputDual.faces[iFaceFinal][i]]);
							
						}
									
						std::vector<double> vecContributingFaceWeights;
					
						std::vector<int> vecContributingFaceI;		
							
						DataArray1D<double> dCoeffs(3);
						
						BilinearWeights(nodeQ,nodesCurrentFace,meshInputDual.faces[iFaceFinal],vecContributingFaceWeights,vecContributingFaceI);
												
						// Contribution of each point to the map
						for (int i = 0; i < vecContributingFaceI.size(); i++){
					
							if( vecContributingFaceWeights[i] < -1e-12 || vecContributingFaceWeights[i] > 1+1e-12 ){
								std::cout << "\nFound weight value = " << vecContributingFaceWeights[i] << std::endl;
								_EXCEPTIONT("Non-monotone weight");
							}
							
							smatMap(iTargetFace, vecContributingFaceI[i]) +=
								vecContributingFaceWeights[i]
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
	
	//If the source mesh has holes, use local dual faces
	
	else {
		
		//Vector of source face centers
		
		NodeVector vecSourceCenterArray;	
		
		for (int i = 0; i < meshInput.faces.size(); i++){
	
			const Face & face = meshInput.faces[i];
	
			Node nodeCentroid = GetFaceCentroid(face, meshInput.nodes);
			
			nodeCentroid = nodeCentroid.Normalized();
			
			vecSourceCenterArray.push_back(nodeCentroid);
	
				
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
			
			int iEdges = faceFirst.edges.size();
			
			//Nodes of the current face
			
			NodeVector nodesFaceFirst;
			
			for (int i = 0; i < faceFirst.edges.size(); i++){
				
				nodesFaceFirst.push_back(meshInput.nodes[faceFirst[i]]);
				
			}	
	
			// Find the set of Faces that overlap faceFirst
			int ixOverlapBegin = ixOverlap;
			
			int ixOverlapEnd = ixOverlapBegin;
	
			for (; ixOverlapEnd < meshOverlap.faces.size(); ixOverlapEnd++) {
				if (meshOverlap.vecSourceFaceIx[ixOverlapEnd] != ixFirst) {
					break;
				}
			}
			
			int nOverlapFaces = ixOverlapEnd - ixOverlapBegin;
	
			if( nOverlapFaces == 0 ) continue;
	
			// Loop through all overlap faces associated with this source face
			for (int j = 0; j < nOverlapFaces; j++) {
				
				int iTargetFace = meshOverlap.vecTargetFaceIx[ixOverlap + j];
	
				// signal to not participate, because it is a ghost target
				if( iTargetFace < 0 ) continue;  // skip and do not do anything
	
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
	
				// Loop through all sub-triangles of this overlap Face
				for (int k = 0; k < nSubTriangles; k++) {
	
					// Loop through all quadrature nodes on this overlap Face
					for (int p = 0; p < dW.GetRows(); p++) {
	
						//These are the weights that are input into the sparse matrix map
				
						std::vector<std::vector<int>> vecContributingFaceI;
			
						std::vector<std::vector<double>> vecContributingFaceWeights;
					
						//These are the temporary weights associated with each source face node.
						//They are used if the current target face isn't in any of the local dual faces
				
						std::vector<std::vector<int>> vecLocalFacesTemp;
				
						std::vector<std::vector<double>> vecLocalWeightsTemp;
	
						// Get quadrature node and pointwise Jacobian
						Node & nodeQ = dQuadPtNodes(k,p);
		
						//Bilinear weights of current target center wrt to nodes of current face
			
						std::vector<double> vecWeightsCurrentFace(iEdges);
										
						//Number of nodes we've visited
						int iCountK = 0;
						
						//loop through all nodes on source face
						
						for (int i = 0; i < iEdges; i++){
							
							//Global index of current node
								
							int iCurrentNode = meshInput.faces[ixFirst][i];
							
							//Number of nodes on local dual test face
							
							int iEdgesTestK = meshInput.revnodearray[iCurrentNode].size();
							
							if (iEdgesTestK < 3){
								
								std::vector<double> vecLocalWeights;
								std::vector<int> vecLocalContributingFaces;		
							
								for (auto it = meshInput.revnodearray[iCurrentNode].begin(); 
									  it != meshInput.revnodearray[iCurrentNode].end(); it++){
											
									vecLocalContributingFaces.push_back(*it);
									vecLocalWeights.push_back(1.0/iEdgesTestK);
								
								
								}
							
								vecLocalWeightsTemp.push_back(vecLocalWeights);
								vecLocalFacesTemp.push_back(vecLocalContributingFaces);
							
								iCountK += 1;
							
							}
							
							//Local dual face has at least 3 edges, i.e. current source node shared by more than 3 faces
							
							else{
				
								Face faceLocalDual(iEdgesTestK);
								
								NodeVector nodesLocalDual;
														
								ConstructLocalDualFace(meshInput,vecSourceCenterArray,iCurrentNode,faceLocalDual,nodesLocalDual);
												
								//Put the dual face nodes in an ccw oriented vector
						
								NodeVector nodesLocalDualReordered;
						
								for (int l = 0; l < nodesLocalDual.size(); l++){
									
									nodesLocalDualReordered.push_back(vecSourceCenterArray[faceLocalDual[l]]);
									
								}
								
								std::vector<double> vecLocalWeights;
								std::vector<int> vecLocalContributingFaces;
								
								//If the local dual face around the current node contains the target, then quit
								
								if (DoesFaceContainPoint(nodesLocalDualReordered, nodeQ.x, nodeQ.y, nodeQ.z)) {
									
										vecContributingFaceWeights.clear();
										vecContributingFaceI.clear();
										
										//Bilinear weights
										
										BilinearWeights(nodeQ,nodesLocalDualReordered,faceLocalDual,vecLocalWeights,vecLocalContributingFaces);
										
										vecContributingFaceI.push_back(vecLocalContributingFaces);
										vecContributingFaceWeights.push_back(vecLocalWeights);								
										
										break;
												
								}
								
								else {
									
									//Add to count
									iCountK += 1;
									
									//Use simple average or bilinear depending on whether the dual face 
									//contains the current node or not
									
									Node nodeCurrentNode = meshInput.nodes[iCurrentNode];
									
									if (DoesFaceContainPoint(nodesLocalDualReordered, nodeCurrentNode.x, nodeCurrentNode.y, 
															nodeCurrentNode.z)){
										
										//Compute bilinear weights of current source face node wrt the local dual face
									
										BilinearWeights(nodeCurrentNode,nodesLocalDualReordered,faceLocalDual,vecLocalWeights,
														vecLocalContributingFaces);
										
										
										vecLocalFacesTemp.push_back(vecLocalContributingFaces);
										vecLocalWeightsTemp.push_back(vecLocalWeights);
										
									}
									
									else {
										
										//Local dual face does not contain current node
										
										vecLocalWeights.resize(iEdgesTestK);
										vecLocalContributingFaces.resize(iEdgesTestK);
										
										for ( int m = 0; m < iEdgesTestK; m++ ){
											
											vecLocalContributingFaces[m] = faceLocalDual[m];
											
											vecLocalWeights[m] = (1.0/iEdgesTestK);
											
										}
										
										vecLocalFacesTemp.push_back(vecLocalContributingFaces);
										vecLocalWeightsTemp.push_back(vecLocalWeights);	
										
										
									}
									
									
								}
											
							}
							
						}
						
						//If we've visited very node of the current source face, we need to do another interpolation
						//from these nodes to the current target face center
									
						if( iCountK == iEdges ){
							
							//New vectors to store bilinear weights of current source face wrt current quadrature node
							
							std::vector<int> vecLocalFaces;
							std::vector<double> vecLocalWeights;
							
							//New face with local ordering
							Face faceCurrentLocal(iEdges);
							
							for ( int r = 0; r < iEdges; r++ ) {
								
								faceCurrentLocal.SetNode(r,r);
								
							}
							
							//compute bilinear weights
							
							BilinearWeights(nodeQ, nodesFaceFirst, faceCurrentLocal, vecLocalWeights, vecLocalFaces);
							
							for (int r = 0; r < vecLocalFaces.size(); r++) {
									
									for (int q = 0; q < vecLocalWeightsTemp[vecLocalFaces[r]].size(); q++ ){
										
										vecLocalWeightsTemp[vecLocalFaces[r]][q] *= vecLocalWeights[r];
										
									}
									
									//Push back to weight and face vectors
									vecContributingFaceI.push_back(vecLocalFacesTemp[vecLocalFaces[r]]);
									vecContributingFaceWeights.push_back(vecLocalWeightsTemp[vecLocalFaces[r]]);
								
							}
						
						}
						
						for (int m = 0; m < vecContributingFaceWeights.size(); m++){
		
							for (int l = 0; l < vecContributingFaceWeights[m].size(); l++){
								
								
								int iContributingFace = vecContributingFaceI[m][l];
			
								if( vecContributingFaceWeights[m][l] < -1e-12 || 
									vecContributingFaceWeights[m][l] > 1+1e-12 ){
										
										std::cout << "\nFound weight value = " << 
										vecContributingFaceWeights[m][l] << std::endl;
										
										_EXCEPTIONT("Non-monotone weight");
									
								}
								
								smatMap(iTargetFace, iContributingFace) +=
											
									vecContributingFaceWeights[m][l]
									* dQuadPtWeight(k,p)
									* meshOverlap.vecFaceArea[ixOverlap + j]
									/ dQuadratureArea
									/ meshOutput.vecFaceArea[iTargetFace];
			
							}
		
						}
						
					}
					
				}
				
			}
			 
			// Increment the current overlap index
			ixOverlap += nOverlapFaces;
			
		}
		
	}

}

///////////////////////////////////////////////////////////////////////////////

void LinearRemapBilinear(	
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

	// Get SparseMatrix representation of the OfflineMap
	SparseMatrix<double> & smatMap = mapRemap.GetSparseMatrix();
	
	//If there are no holes on the source mesh, use the (global) dual of the source.
	
	if (!DoesMeshHaveHoles(meshInput)) {
		
		Mesh meshInputDual = meshInput;
		
		//Construct dual mesh
		Dual(meshInputDual);
			
		//Construct edge map
		
		meshInputDual.ConstructEdgeMap();
		
		//kd-tree of dual mesh centers
	  kdtree * kdTarget = kd_create(3);
		for (int i = 0; i < meshInputDual.faces.size(); i++){

			const Face & face = meshInputDual.faces[i];

			Node nodeCentroid = GetFaceCentroid(face, meshInputDual.nodes);
		
			nodeCentroid = nodeCentroid.Normalized();

			kd_insert3(
				kdTarget,
				nodeCentroid.x,
				nodeCentroid.y,
				nodeCentroid.z,
				(void*)(&(meshInputDual.faces[i])));
			
		}
		
		// Loop through all target faces
		for (int ixFirst = 0; ixFirst < meshOutput.faces.size(); ixFirst++) {
	
			// Output every 1000 overlap elements
			if (ixFirst % 1000 == 0) {
				Announce("Element %i/%i", ixFirst, meshOutput.faces.size());
			}
	
			// This Face
			const Face & faceFirst = meshOutput.faces[ixFirst];
			
			// Get node coordinates of each target face center
			Node nodeQ = GetFaceCentroid(faceFirst,meshOutput.nodes);
			
			nodeQ = nodeQ.Normalized();
					
			//Get the dual mesh face whose center is nearest to the target face
			
			kdres * kdresTarget =
				kd_nearest3(
					kdTarget,
					nodeQ.x,
					nodeQ.y,
					nodeQ.z);
					
			Face * pFace = (Face *)(kd_res_item_data(kdresTarget));
	
			int iNearestFace = pFace - &(meshInputDual.faces[0]);
			
			// Find face that contains nodeQ
			// This is the face whose local index is iFaceFinal
			
			int iFaceFinal = 0;
			
			GetFaceThatContainsPoint(meshInputDual,iNearestFace,iFaceFinal,nodeQ.x,nodeQ.y,nodeQ.z);
			
			int iEdges = meshInputDual.faces[iFaceFinal].edges.size();
			
			//NodeVector nodesCurrentFace(iEdges);
			NodeVector nodesCurrentFace;
			
			for (int i = 0; i < iEdges; i++){
				
				nodesCurrentFace.push_back(meshInputDual.nodes[meshInputDual.faces[iFaceFinal][i]]);
				
			}
						
			std::vector<double> vecContributingFaceWeights;
		
			std::vector<int> vecContributingFaceI;		
				
			DataArray1D<double> dCoeffs(3);
			
			BilinearWeights(nodeQ,nodesCurrentFace,meshInputDual.faces[iFaceFinal],vecContributingFaceWeights,vecContributingFaceI);
			
			// Contribution of each point to the map
			for (int i = 0; i < vecContributingFaceI.size(); i++){
				
			  if( vecContributingFaceWeights[i] < -1e-12 || vecContributingFaceWeights[i] > 1+1e-12 ){
              std::cout << "\nFound weight value = " << vecContributingFaceWeights[i] << std::endl;
						 _EXCEPTIONT("Non-monotone weight");
        }
	
				int iContributingFaceI = vecContributingFaceI[i];
				
				smatMap(ixFirst, iContributingFaceI) = vecContributingFaceWeights[i];
				
			}					
											
		}
		
	}
	
	//If the source mesh has holes, use local dual faces
	
	else {

		//kd-tree of dual mesh centers
		kdtree * kdTarget = kd_create(3);
    
		//Vector of source face centers
		NodeVector vecSourceCenterArray;
    
		//Vector of target face centers
		NodeVector vecTargetCenterArray;

		// Vector of centers of the source mesh
		for (int i = 0; i < meshInput.faces.size(); i++){
	
			const Face & face = meshInput.faces[i];
	
			Node nodeCentroid = GetFaceCentroid(face, meshInput.nodes);
			
			nodeCentroid = nodeCentroid.Normalized();
			
			vecSourceCenterArray.push_back(nodeCentroid);
	
				
		}
		
		for (int i = 0; i < meshOutput.faces.size(); i++){

			const Face & face = meshOutput.faces[i];

			Node nodeCentroid = GetFaceCentroid(face, meshOutput.nodes);
		
			nodeCentroid = nodeCentroid.Normalized();

			kd_insert3(
				kdTarget,
				nodeCentroid.x,
				nodeCentroid.y,
				nodeCentroid.z,
				(void*)(&(meshOutput.faces[i])));
			
			vecTargetCenterArray.push_back(nodeCentroid);
			
		}		
		
		std::set<int> setAllTargets;
		
		// Loop through all source faces
		for (int ixFirst = 0; ixFirst < meshInput.faces.size(); ixFirst++) {
	
			// Output every 1000 elements
			if (ixFirst % 1000 == 0) {
				Announce("Element %i/%i", ixFirst, meshInput.faces.size());
			}
			
			// This Face
			const Face & faceFirst = meshInput.faces[ixFirst];
			
			//Find the maximum distance from the current source face center to its nodes
			
			double dMaxNodeDist = 0.0;
			
			Node faceFirstCenter = vecSourceCenterArray[ixFirst];
			
			int iEdges = faceFirst.edges.size();
			
			for (int i = 0; i < iEdges; i++){
				
				//Great circle distance is just the arccos of dot product
				
				Node nodeI = meshInput.nodes[faceFirst[i]];
				
				double dDistToNodeI = acos(DotProduct(nodeI, faceFirstCenter));
				
				if ( dDistToNodeI > dMaxNodeDist ){
					
					dMaxNodeDist = dDistToNodeI;
					
				}
				
			}
			
			//Find all the target face centers that are within a given distance of the given source face center
			
			/* find points closest to the origin and within distance radius */
			//struct *presults;
										 
			struct kdres *presults = kd_nearest_range3(kdTarget,faceFirstCenter.x, faceFirstCenter.y, 
													   faceFirstCenter.z, dMaxNodeDist+.0001);
			
										   
			if( kd_res_size(presults) == 0 ){
				
				continue;
				
			}										   
					   
			//Nodes of the current face
			
			NodeVector nodesFaceFirst;
			
			for (int i = 0; i < faceFirst.edges.size(); i++){
				
				nodesFaceFirst.push_back(meshInput.nodes[faceFirst[i]]);
				
			}
			
			//set of target faces that have already been used.  This is necessry to avoid double
			//counting of target face centers that lie on a source face edge/node
			
			std::vector<int> vecTargetsInFace;
			
			while( !kd_res_end( presults ) ) {
				
				
			    /* get the data and position of the current result item */
			    Face * pFace = (Face *)kd_res_item_data( presults );
			    
			    int iNearestTarget = pFace - &(meshOutput.faces[0]);
			    			    
			    if(DoesFaceContainPoint(nodesFaceFirst, vecTargetCenterArray[iNearestTarget].x,
										vecTargetCenterArray[iNearestTarget].y,vecTargetCenterArray[iNearestTarget].z)) {
										
					if( setAllTargets.find(iNearestTarget) == setAllTargets.end() ){
						
						setAllTargets.insert(iNearestTarget);
						vecTargetsInFace.push_back(iNearestTarget);
						
					}
						
				}
								
			    /* go to the next entry */
			    kd_res_next( presults );
			    
			}
			
			kd_res_free( presults );

			//Loop through all centers of the target faces in the given source face
			for (int i = 0; i < vecTargetsInFace.size(); i++){
				
				//These are the weights that are input into the spares matrix map
				
				std::vector<std::vector<int>> vecContributingFaceI;
			
				std::vector<std::vector<double>> vecContributingFaceWeights;
				
				//These are the temporary weights associated with each source face node.
				//They are used if the current target face isn't in any of the local dual faces
				
				std::vector<std::vector<int>> vecLocalFacesTemp;
				
				std::vector<std::vector<double>> vecLocalWeightsTemp;
				
				//Current target face center
				
				Node nodeQ = vecTargetCenterArray[vecTargetsInFace[i]];
				
				//Bilinear weights of current target center wrt to nodes of current face
			
				std::vector<double> vecWeightsCurrentFace(iEdges);
								
				//Number of nodes we've visited
				int iCountK = 0;
				
				//loop through all nodes on source face
				
				for (int k = 0; k < iEdges; k++){
					
					//Global index of current node
						
					int iCurrentNode = meshInput.faces[ixFirst][k];
					
					//Number of nodes on local dual test face
					
					int iEdgesTestK = meshInput.revnodearray[iCurrentNode].size();
					
					if (iEdgesTestK < 3){
						
						std::vector<double> vecLocalWeights;
						std::vector<int> vecLocalContributingFaces;		
						
						for (auto it = meshInput.revnodearray[iCurrentNode].begin(); 
								  it != meshInput.revnodearray[iCurrentNode].end(); it++){
									  
							vecLocalContributingFaces.push_back(*it);
							vecLocalWeights.push_back(1.0/iEdgesTestK);
							
						
						}
						
						vecLocalWeightsTemp.push_back(vecLocalWeights);
						vecLocalFacesTemp.push_back(vecLocalContributingFaces);
						
						iCountK += 1;
						
					}
					
					//Local dual face has at least three edges, i.e. current source node shared by more than 3 faces
					
					else{
		
						Face faceLocalDual(iEdgesTestK);
						
						NodeVector nodesLocalDual;
												
						ConstructLocalDualFace(meshInput,vecSourceCenterArray,iCurrentNode,faceLocalDual,nodesLocalDual);
										
						//Put the dual face nodes in an ccw oriented vector
				
						NodeVector nodesLocalDualReordered;
				
						for (int j = 0; j < nodesLocalDual.size(); j++){
							
							nodesLocalDualReordered.push_back(vecSourceCenterArray[faceLocalDual[j]]);
							
						}
						
						std::vector<double> vecLocalWeights;
						std::vector<int> vecLocalContributingFaces;
						
						//If the local dual face around the current node contains the target, then quit
						
						if (DoesFaceContainPoint(nodesLocalDualReordered, nodeQ.x, nodeQ.y, nodeQ.z)) {
							
								vecContributingFaceWeights.clear();
								vecContributingFaceI.clear();
								
								//Bilinear weights
								
								BilinearWeights(nodeQ,nodesLocalDualReordered,faceLocalDual,vecLocalWeights,vecLocalContributingFaces);
								
								vecContributingFaceI.push_back(vecLocalContributingFaces);
								vecContributingFaceWeights.push_back(vecLocalWeights);								
								
								break;
										
						}
						
						else {
							
							//Add to count
							iCountK += 1;
							
							//Use simple average or barycentric coorindates depending on whether the dual face 
							//contains the current node or not
							
							Node nodeCurrentNode = meshInput.nodes[iCurrentNode];
							
							if (DoesFaceContainPoint(nodesLocalDualReordered, nodeCurrentNode.x, nodeCurrentNode.y, 
													nodeCurrentNode.z)){
								
								//Compute bilinear weights of current source face node wrt the local dual face
							
								BilinearWeights(nodeCurrentNode,nodesLocalDualReordered,faceLocalDual,vecLocalWeights,
												vecLocalContributingFaces);
								
								
								vecLocalFacesTemp.push_back(vecLocalContributingFaces);
								vecLocalWeightsTemp.push_back(vecLocalWeights);
								
							}
							
							else {
								
								vecLocalWeights.resize(iEdgesTestK);
								vecLocalContributingFaces.resize(iEdgesTestK);
								
								for ( int m = 0; m < iEdgesTestK; m++ ){
									
									vecLocalContributingFaces[m] = (faceLocalDual[m]);
									
									vecLocalWeights[m] = (1.0/iEdgesTestK);
									
								}
								
								vecLocalFacesTemp.push_back(vecLocalContributingFaces);
								vecLocalWeightsTemp.push_back(vecLocalWeights);	
								
								
							}					
							
						}
									
					}
					
				}
				
				//If we've visited very node of the current source face, we need to do another interpolation
				//from these nodes to the current target face center
				
				if( iCountK == iEdges ){
					
					//New vectors to store bilinear weights of current source face wrt current target face center
					
					std::vector<int> vecLocalFaces;
					std::vector<double> vecLocalWeights;
					
					//New face with local ordering
					Face faceCurrentLocal(iEdges);
					
					for ( int r = 0; r < iEdges; r++ ) {
						
						faceCurrentLocal.SetNode(r,r);
						
					}
					
					//compute bilinear weights
					
					BilinearWeights(nodeQ, nodesFaceFirst, faceCurrentLocal, vecLocalWeights, vecLocalFaces);
					
					for (int r = 0; r < vecLocalFaces.size(); r++) {
							
							for (int q = 0; q < vecLocalWeightsTemp[vecLocalFaces[r]].size(); q++ ){
								
								vecLocalWeightsTemp[vecLocalFaces[r]][q] *= vecLocalWeights[r];
								
							}
							
							//Push back to weight and face vectors
							vecContributingFaceI.push_back(vecLocalFacesTemp[vecLocalFaces[r]]);
							vecContributingFaceWeights.push_back(vecLocalWeightsTemp[vecLocalFaces[r]]);
						
					}
				
				}
			
				for (int m = 0; m < vecContributingFaceWeights.size(); m++){
					
					for (int j = 0; j < vecContributingFaceWeights[m].size(); j++){
												
						int iContributingFace = vecContributingFaceI[m][j];
						
						if( vecContributingFaceWeights[m][j] < -1e-12 || vecContributingFaceWeights[m][j] > 1+1e-12 ){
              std::cout << "\nFound weight value = " << vecContributingFaceWeights[m][j] << std::endl;
						 _EXCEPTIONT("Non-monotone weight");
								
						}
						
						smatMap(vecTargetsInFace[i],iContributingFace) += vecContributingFaceWeights[m][j];
						
					}
					
				}
			
			}
		
		}

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

		if( nOverlapFaces == 0 ) continue;

		// Loop through all Overlap Faces
		for (int i = 0; i < nOverlapFaces; i++) {

			// Quantities from the overlap Mesh
			const Face & faceOverlap = meshOverlap.faces[ixOverlap + i];

			// Quantities from the Second Mesh
			int ixSecondFace = meshOverlap.vecTargetFaceIx[ixOverlap + i];

			// signal to not participate, because it is a ghost target
			if( ixSecondFace < 0 ) continue;  // skip and do not do anything

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

			// signal to not participate, because it is a ghost target
			if( iTargetFace < 0 ) continue;  // skip and do not do anything

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

			// signal to not participate, because it is a ghost target
			if( ixSecondFace < 0 ) continue;  // skip and do not do anything

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

		if( nOverlapFaces == 0 ) continue;

		// Loop through all Overlap Faces
		for (int i = 0; i < nOverlapFaces; i++) {

			// Quantities from the overlap Mesh
			const Face & faceOverlap = meshOverlap.faces[ixOverlap + i];

			const NodeVector & nodesOverlap = meshOverlap.nodes;

			// Quantities from the Second Mesh
			int ixSecond = meshOverlap.vecTargetFaceIx[ixOverlap + i];

			// signal to not participate, because it is a ghost target
			if( ixSecond < 0 ) continue;  // skip and do not do anything

			const NodeVector & nodesSecond = meshOutput.nodes;

			const Face & faceSecond = meshOutput.faces[ixSecond];
			int nbEdges = faceOverlap.edges.size();
			int nOverlapTriangles = 1;
			Node nodeCenter; // not used if nbEdges == 3
			if (nbEdges > 3) { // decompose from nodeCenter in this case
				nOverlapTriangles = nbEdges;
				for (int k = 0; k < nbEdges; k++) {
					const Node &node = nodesOverlap[faceOverlap[k]];
					nodeCenter = nodeCenter + node;
				}
				nodeCenter = nodeCenter / nbEdges;
				double dMagni = sqrt(
						nodeCenter.x * nodeCenter.x + nodeCenter.y * nodeCenter.y
								+ nodeCenter.z * nodeCenter.z);
				nodeCenter = nodeCenter / dMagni; // project back on sphere of radius 1
			}

			Node node0, node1, node2;
			double dTriArea;
			// Loop over all sub-triangles of this Overlap Face
			for (int j = 0; j < nOverlapTriangles; j++) {

				if (nbEdges == 3) { // will come here only once, nOverlapTriangles == 1 in this case
					node0 = nodesOverlap[faceOverlap[0]];
					node1 = nodesOverlap[faceOverlap[1]];
					node2 = nodesOverlap[faceOverlap[2]];
				}
				else { // decompose polygon in triangles around the nodeCenter
					node0 = nodeCenter;
					node1 = nodesOverlap[faceOverlap[j]];
					int j1 = (j + 1) % nbEdges;
					node2 = nodesOverlap[faceOverlap[j1]];
				}
				dTriArea = CalculateTriangleAreaQuadratureMethod(node0, node1,
						node2);
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

			// signal to not participate, because it is a ghost target
			if( ixSecond < 0 ) continue;  // skip and do not do anything

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

		// signal to not participate, because it is a ghost target
		if( ixSecond < 0 ) continue;  // skip and do not do anything

		vecReverseFaceIx[ixSecond].push_back(i);
	}
/*
	for (int ixOverlap = 0; ixOverlap < meshOverlap.faces.size(); ixOverlap++) {
		int ixSecond = meshOverlap.vecTargetFaceIx[ixOverlap];

    // signal to not participate, because it is a ghost target
    if( ixSecond < 0 ) continue;  // skip and do not do anything

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

    // signal to not participate, because it is a ghost target
    if( ixSecond < 0 ) continue;  // skip and do not do anything

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

    // signal to not participate, because it is a ghost target
    if( ixSecond < 0 ) continue;  // skip and do not do anything

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

    // signal to not participate, because it is a ghost target
    if( ixSecond < 0 ) continue;  // skip and do not do anything

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

			// signal to not participate, because it is a ghost target
			if( ixSecond < 0 ) continue;  // skip and do not do anything

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

			// signal to not participate, because it is a ghost target
			if( ixSecond < 0 ) continue;  // skip and do not do anything

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

			// signal to not participate, because it is a ghost target
			if( ixSecondFace < 0 ) continue;  // skip and do not do anything

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

		if( nOverlapFaces == 0 ) continue;
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

			// signal to not participate, because it is a ghost target
			if( ixSecond < 0 ) continue;  // skip and do not do anything

			const NodeVector & nodesSecond = meshOutput.nodes;
			const Face & faceSecond = meshOutput.faces[ixSecond];
			int nbEdges = faceOverlap.edges.size();
			int nOverlapTriangles = 1;
			Node nodeCenter; // not used if nbEdges == 3
			if (nbEdges > 3) { // decompose from nodeCenter in this case
				nOverlapTriangles = nbEdges;
				for (int k = 0; k < nbEdges; k++) {
					const Node &node = nodesOverlap[faceOverlap[k]];
					nodeCenter = nodeCenter + node;
				}
				nodeCenter = nodeCenter / nbEdges;
				double dMagni = sqrt(
						nodeCenter.x * nodeCenter.x + nodeCenter.y * nodeCenter.y
								+ nodeCenter.z * nodeCenter.z);
				nodeCenter = nodeCenter / dMagni; // project back on sphere of radius 1
			}

			Node node0, node1, node2;
			double dTriArea;

			// Loop over all sub-triangles of this Overlap Face
			for (int j = 0; j < nOverlapTriangles; j++) {

				if (nbEdges == 3) { // will come here only once, nOverlapTriangles == 1 in this case
					node0 = nodesOverlap[faceOverlap[0]];
					node1 = nodesOverlap[faceOverlap[1]];
					node2 = nodesOverlap[faceOverlap[2]];
				}
				else { // decompose polygon in triangles around the nodeCenter
					node0 = nodeCenter;
					node1 = nodesOverlap[faceOverlap[j]];
					int j1 = (j + 1) % nbEdges;
					node2 = nodesOverlap[faceOverlap[j1]];
				}
				dTriArea = CalculateTriangleAreaQuadratureMethod(node0, node1,
						node2);

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

			// signal to not participate, because it is a ghost target
			if( ixSecondFace < 0 ) continue;  // skip and do not do anything

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

			// signal to not participate, because it is a ghost target
			if( ixSecond < 0 ) continue;  // skip and do not do anything

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

		if( nOverlapFaces == 0 ) continue;

		// Put composed array into map
		for (int j = 0; j < nOverlapFaces; j++) {
			int ixSecondFace = meshOverlap.vecTargetFaceIx[ixOverlap + j];

			// signal to not participate, because it is a ghost target
			if( ixSecondFace < 0 ) continue;  // skip and do not do anything

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

		if( nOverlapFaces == 0 ) continue;

		// Loop through all Overlap Faces
		for (int i = 0; i < nOverlapFaces; i++) {

			// Quantities from the overlap Mesh
			const Face & faceOverlap = meshOverlap.faces[ixOverlap + i];

			const NodeVector & nodesOverlap = meshOverlap.nodes;

			// Quantities from the Second Mesh
			int ixSecond = meshOverlap.vecTargetFaceIx[ixOverlap + i];

			// signal to not participate, because it is a ghost target
			if( ixSecond < 0 ) continue;  // skip and do not do anything

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

