///////////////////////////////////////////////////////////////////////////////
///
///	\file    FiniteVolumeTools.cpp
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

#include "Defines.h"
#include "FiniteVolumeTools.h"
#include "MathHelper.h"
#include "Announce.h"
#include "Exception.h"
#include "CoordTransforms.h"

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

void GetAdjacentFaceVectorByEdge(
	const Mesh & mesh,
	int iFaceInitial,
	int nRequiredFaceSetSize,
	AdjacentFaceVector & vecFaces
) {
	// Ensure the ReverseNodeArray has been constructed
	if (mesh.edgemap.size() == 0) {
		_EXCEPTIONT("EdgeMap is required");
	}

	// Insert the initial Face
	vecFaces.push_back(FaceDistancePair(iFaceInitial,1));

	// Set of all Faces
	std::set<int> setAllFaces;

	// Set of current Faces
	std::set<int> setCurrentFaces;

	// Insert initial Face
	setAllFaces.insert(iFaceInitial);
	setCurrentFaces.insert(iFaceInitial);

	// Generate the set of faces
	int iDistance = 1;
	for (;;) {

		if (vecFaces.size() >= nRequiredFaceSetSize) {
			break;
		}

		// Increment distance metric
		iDistance++;

		// Set of Faces to examine next
		std::set<int> setNextFaces;

		// Loop through all Faces adjacent to Faces in setCurrentFaces
		std::set<int>::const_iterator iterCurrentFace = setCurrentFaces.begin();
		for (; iterCurrentFace != setCurrentFaces.end(); iterCurrentFace++) {

			const Face & faceCurrent = mesh.faces[*iterCurrentFace];
			for (int i = 0; i < faceCurrent.edges.size(); i++) {
				const FacePair & facepair =
					mesh.edgemap.find(faceCurrent.edges[i])->second;

				// New face index
				int iNewFace;
				if (facepair[0] == *iterCurrentFace) {
					iNewFace = facepair[1];
				} else if (facepair[1] == *iterCurrentFace) {
					iNewFace = facepair[0];
				} else {
					_EXCEPTIONT("Logic error");
				}

				if (iNewFace == InvalidFace) {
					continue;
				}

				// If this is a new Face, insert into vector
				if (setAllFaces.find(iNewFace) == setAllFaces.end()) {
					vecFaces.push_back(FaceDistancePair(iNewFace, iDistance));
					setAllFaces.insert(iNewFace);
					setNextFaces.insert(iNewFace);
				}
			}
		}

		setCurrentFaces = setNextFaces;
	}

	if (vecFaces.size() < nRequiredFaceSetSize) {
		_EXCEPTION1("Unable to find enough adjacent faces to Face %i to meet required size", iFaceInitial);
	}
}

///////////////////////////////////////////////////////////////////////////////

void GetTriangleThatContainsPoint(
	const Mesh & mesh,
	int iFaceInitial,
	int & iFaceFinal,
	double dX,
	double dY
) {
	// Ensure the ReverseNodeArray has been constructed
	if (mesh.edgemap.size() == 0) {
		_EXCEPTIONT("EdgeMap is required");
	}
	
	// First check if point is in initial face
	if (DoesTriangleContainPoint(mesh,iFaceInitial,dX,dY)){
		
		iFaceFinal = iFaceInitial;
		return;
						
	}
	
	// Set of all Faces
	std::set<int> setAllFaces;

	// Set of current Faces
	std::set<int> setCurrentFaces;
	
	// Insert initial Face
	setAllFaces.insert(iFaceInitial);
	setCurrentFaces.insert(iFaceInitial);
	
	while (setAllFaces.size() < mesh.faces.size()) {

		// Set of Faces to examine next
		std::set<int> setNextFaces;

		// Loop through all Faces adjacent to Faces in setCurrentFaces
		std::set<int>::const_iterator iterCurrentFace = setCurrentFaces.begin();
		for (; iterCurrentFace != setCurrentFaces.end(); iterCurrentFace++) {

			const Face & faceCurrent = mesh.faces[*iterCurrentFace];
			for (int i = 0; i < faceCurrent.edges.size(); i++) {
				const FacePair & facepair =
					mesh.edgemap.find(faceCurrent.edges[i])->second;
				
				// New face index
				int iNewFace;
				if (facepair[0] == *iterCurrentFace) {
					iNewFace = facepair[1];
				} else if (facepair[1] == *iterCurrentFace) {
					iNewFace = facepair[0];
				} else {
					_EXCEPTIONT("Logic error");
				}

				if (iNewFace == InvalidFace) {
					continue;
				}
				
				// If this is a new Face, check whether it contains the point
				if (setAllFaces.find(iNewFace) == setAllFaces.end()) {
					if(DoesTriangleContainPoint(mesh,iNewFace,dX,dY)){
						
						iFaceFinal = iNewFace;
						return;
						
					}
					else{
					
						setAllFaces.insert(iNewFace);
						setNextFaces.insert(iNewFace);
						
					}
				}
			}
		}

		setCurrentFaces = setNextFaces;
	}
	_EXCEPTIONT("Unable to find a triangle that contains the point");
}

///////////////////////////////////////////////////////////////////////////////

void GetFaceThatContainsPoint(
	const Mesh & mesh,
	int iFaceInitial,
	int & iFaceFinal,
	double dX,
	double dY,
	double dZ

){
	
	// Ensure the ReverseNodeArray has been constructed
	if (mesh.edgemap.size() == 0) {
		_EXCEPTIONT("EdgeMap is required");
	}
	
	// First check if point is in initial face
	if (DoesFaceContainPoint(mesh,iFaceInitial,dX,dY,dZ)){
		
		iFaceFinal = iFaceInitial;
		return;
						
	}
	
	// Set of all Faces
	std::set<int> setAllFaces;

	// Set of current Faces
	std::set<int> setCurrentFaces;
	
	// Insert initial Face
	setAllFaces.insert(iFaceInitial);
	setCurrentFaces.insert(iFaceInitial);
	
	while (setAllFaces.size() < mesh.faces.size()) {

		// Set of Faces to examine next
		std::set<int> setNextFaces;

		// Loop through all Faces adjacent to Faces in setCurrentFaces
		std::set<int>::const_iterator iterCurrentFace = setCurrentFaces.begin();
		for (; iterCurrentFace != setCurrentFaces.end(); iterCurrentFace++) {

			const Face & faceCurrent = mesh.faces[*iterCurrentFace];
			for (int i = 0; i < faceCurrent.edges.size(); i++) {
				const FacePair & facepair =
					mesh.edgemap.find(faceCurrent.edges[i])->second;
				
				// New face index
				int iNewFace;
				if (facepair[0] == *iterCurrentFace) {
					iNewFace = facepair[1];
				} else if (facepair[1] == *iterCurrentFace) {
					iNewFace = facepair[0];
				} else {
					_EXCEPTIONT("Logic error");
				}

				if (iNewFace == InvalidFace) {
					continue;
				}
				
				// If this is a new Face, check whether it contains the point
				if (setAllFaces.find(iNewFace) == setAllFaces.end()) {
					if(DoesFaceContainPoint(mesh,iNewFace,dX,dY,dZ)){
						
						iFaceFinal = iNewFace;
						return;
						
					}
					else{
					
						setAllFaces.insert(iNewFace);
						setNextFaces.insert(iNewFace);
						
					}
				}
			}
		}

		setCurrentFaces = setNextFaces;
	}
	_EXCEPTIONT("Unable to find a face that contains the point");
}

///////////////////////////////////////////////////////////////////////////////

bool DoesFaceContainPoint(
	const Mesh & mesh,
	int iFace,
	double dX,
	double dY,
	double dZ

){
	
	//Convert sample point to lat/lon coordinates
	
	double dLonRad0 = 0.0;
	double dLatRad0 = 0.0;
	
	XYZtoRLL_Rad(dX,dY,dZ,dLonRad0,dLatRad0);
	
	int iEdges = mesh.faces[iFace].edges.size();
	
	NodeVector nodesPlane;
	
	double dCenterX = 0;
	double dCenterY = 0;
	
	// Project face nodes onto plane tangent to the sphere at the sample point whose 
	// coordinates are dX, dY, dZ
	for (int i = 0; i < iEdges; i++){
		
		Node nodeCurrent = mesh.nodes[mesh.faces[iFace][i]];
		
		double dLonRad = 0.0;
		double dLatRad = 0.0;
		
		//Convert to lat/lon coordinates
		
		XYZtoRLL_Rad(nodeCurrent.x,nodeCurrent.y,nodeCurrent.z,dLonRad,dLatRad);
		
		//Project on tangent plane at sample point
		
		double dGX = 0.0;
		double dGY = 0.0;
		
		GnomonicProjection(dLonRad0,dLatRad0,dLonRad,dLatRad,dGX,dGY);
		
		Node nodeI(dGX,dGY,0.0);
		
		nodesPlane.push_back(nodeI);
		
		//add contributions to planar face centroid
		dCenterX += dGX;
		dCenterY += dGY;
		
	}
	
	//Compute coordinates of planar face centroid and turn it into a Node
	
	dCenterX /= iEdges;
	dCenterY /= iEdges;
		
	Node nodeCenter(dCenterX,dCenterY,0);
	
	//loop over all edges
	for (int i = 0; i < iEdges; i++){
		
		Node nodeI = nodesPlane[i];
		
		Node nodeIPlusOne = nodesPlane[(i+1)%iEdges];
		
		Node nodeEdge(nodeIPlusOne.x - nodeI.x,
					  nodeIPlusOne.y - nodeI.y,
					  0.0);
		
		Node nodeEdgeNormal(-nodeEdge.y, nodeEdge.x, 0.0);
		
		
		Node nodeCenterMinusI(nodeCenter.x - nodeI.x,
							  nodeCenter.y - nodeI.y,
							  0.0);
		
		//This is the projected sampled point (whose coordinates are zero in the gnomonic plane)
		//minus nodeI
						  
		Node nodeQ(-nodeI.x,-nodeI.y,0.0);
					
		//If the signs of the dot products are different, than the sample point is outside the polygon		  
		
		if((DotProduct(nodeCenterMinusI,nodeEdgeNormal) > 0 && DotProduct(nodeQ,nodeEdgeNormal) < 0) ||
			(DotProduct(nodeCenterMinusI,nodeEdgeNormal) < 0 && DotProduct(nodeQ,nodeEdgeNormal) > 0)){
				
				return false;
				
		}
		
	}
	
	return true;
	
}

///////////////////////////////////////////////////////////////////////////////

bool DoesFaceContainPoint(
	const NodeVector & nodesP,
	double dX,
	double dY,
	double dZ

){
	
	//Convert sample point to lat/lon coordinates
	
	double dLonRad0 = 0.0;
	double dLatRad0 = 0.0;
	
	XYZtoRLL_Rad(dX,dY,dZ,dLonRad0,dLatRad0);
	
	int iEdges = nodesP.size();
	
	NodeVector nodesPlane;
	
	double dCenterX = 0;
	double dCenterY = 0;
	
	// Project face nodes onto plane tangent to the sphere at the sample point whose 
	// coordinates are dX, dY, dZ
	for (int i = 0; i < iEdges; i++){
		
		Node nodeCurrent = nodesP[i];
		
		double dLonRad = 0.0;
		double dLatRad = 0.0;
		
		//Convert to lat/lon coordinates
		
		XYZtoRLL_Rad(nodeCurrent.x,nodeCurrent.y,nodeCurrent.z,dLonRad,dLatRad);
		
		//Project on tangent plane at sample point
		
		double dGX = 0.0;
		double dGY = 0.0;
		
		GnomonicProjection(dLonRad0,dLatRad0,dLonRad,dLatRad,dGX,dGY);
		
		Node nodeI(dGX,dGY,0.0);
		
		nodesPlane.push_back(nodeI);
		
		//add contributions to planar face centroid
		dCenterX += dGX;
		dCenterY += dGY;
		
	}
	
	//Compute coordinates of planar face centroid and turn it into a Node
	
	dCenterX /= iEdges;
	dCenterY /= iEdges;
		
	Node nodeCenter(dCenterX,dCenterY,0);
	
	//loop over all edges
	for (int i = 0; i < iEdges; i++){
		
		Node nodeI = nodesPlane[i];
		
		Node nodeIPlusOne = nodesPlane[(i+1)%iEdges];
		
		Node nodeEdge(nodeIPlusOne.x - nodeI.x,
					  nodeIPlusOne.y - nodeI.y,
					  0.0);
		
		Node nodeEdgeNormal(-nodeEdge.y, nodeEdge.x, 0.0);
		
		
		Node nodeCenterMinusI(nodeCenter.x - nodeI.x,
							  nodeCenter.y - nodeI.y,
							  0.0);
		
		//This is the projected sampled point (whose coordinates are zero in the gnomonic plane)
		//minus nodeI
						  
		Node nodeQ(-nodeI.x,-nodeI.y,0.0);
					
		//If the signs of the dot products are different, than the sample point is outside the polygon		  
		
		if((DotProduct(nodeCenterMinusI,nodeEdgeNormal) > 0.0 && DotProduct(nodeQ,nodeEdgeNormal) < 0.0) ||
			(DotProduct(nodeCenterMinusI,nodeEdgeNormal) < 0.0 && DotProduct(nodeQ,nodeEdgeNormal) > 0.0)){
				
				return false;
				
		}
		
	}
	
	return true;
	
}

///////////////////////////////////////////////////////////////////////////////

void BarycentricCoordinates(
	const Mesh & mesh,
	int iFace,
	double dX,
	double dY,
	double & dA,
	double & dB
){
	
	Face face = mesh.faces[iFace];
	
	// The input face has to be a triangle
	if(face.edges.size() != 3){
		_EXCEPTIONT("The input face must be a triangle");
	}
	
	Node nodeV1 = mesh.nodes[face[0]];
	Node nodeV2 = mesh.nodes[face[1]];
	Node nodeV3 = mesh.nodes[face[2]];
	
	double dX1 = nodeV1.x;
	double dY1 = nodeV1.y;
	
	double dX2 = nodeV2.x;
	double dY2 = nodeV2.y;
	
	double dX3 = nodeV3.x;
	double dY3 = nodeV3.y;
	
	double dDenom = (dY2 - dY3)*(dX1 - dX3) + (dX3 - dX2)*(dY1 - dY3);
	
	if( abs(dDenom) > 1e-12 ){
		
		_EXCEPTIONT("Points are close to colinear");
		
	}
	
	dA = ((dY2 - dY3)*(dX - dX3) + (dX3 - dX2)*(dY - dY3))/dDenom;
		 
	dB = ((dY3 - dY1)*(dX - dX3) + (dX1 - dX3)*(dY - dY3))/dDenom;
		
}

///////////////////////////////////////////////////////////////////////////////

bool DoesTriangleContainPoint(
	const Mesh & mesh,
	int iFace,
	double dX,
	double dY
){
	
	double dA;
	double dB;
	
	BarycentricCoordinates(mesh,iFace,dX,dY,dA,dB);
		
	if ((0.0 <= dA) && (0.0 <= dB) && (dA + dB <= 1+ReferenceTolerance)){
		return true;
	}
	else{
		return false;
	}
	
}

///////////////////////////////////////////////////////////////////////////////

void GetAdjacentFaceVectorByNode(
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

void MatVectorMult(const DataArray2D<double> & dMat,
				   DataArray1D<double> & dRHS,
				   DataArray1D<double> & dOutput)
{

		dOutput(0) = (dMat(0,0))*dRHS(0) + (dMat(0,1))*dRHS(1) + (dMat(0,2))*dRHS(2);
		dOutput(1) = (dMat(1,0))*dRHS(0) + (dMat(1,1))*dRHS(1) + (dMat(1,2))*dRHS(2);
		dOutput(2) = (dMat(2,0))*dRHS(0) + (dMat(2,1))*dRHS(1) + (dMat(2,2))*dRHS(2);
	
}

///////////////////////////////////////////////////////////////////////////////

void TriangleLineIntersection(
	Node & nodeQ,
	NodeVector & nodesP,
	DataArray1D<double> & dCoeffs,
	double & dCond
) {
	
	_ASSERT(dCoeffs.GetRows() == 3);
	
	//Setup up columns of 3x3 matrix
	DataArray2D<double> dInterpMat(3,3);
	
	dInterpMat(0,0) = nodeQ.x;
	dInterpMat(1,0) = nodeQ.y;
	dInterpMat(2,0) = nodeQ.z;
	
	dInterpMat(0,1) = nodesP[1].x - nodesP[0].x;						  
	dInterpMat(1,1) = nodesP[1].y - nodesP[0].y;
	dInterpMat(2,1) = nodesP[1].z - nodesP[0].z;
	dInterpMat(0,2) = nodesP[2].x - nodesP[0].x;						  
	dInterpMat(1,2) = nodesP[2].y - nodesP[0].y;
	dInterpMat(2,2) = nodesP[2].z - nodesP[0].z;	
					  			
	int m = 3;
	int n = 3;
	int lda = 3;
	int info;
	
	DataArray1D<int> iPIV;
	iPIV.Allocate(3);
	
	DataArray1D<double> dWork;
	dWork.Allocate(3);
	
	int lWork = 3;
	
	//Column sums for A and A inverse
	DataArray1D<double> dColSumA(3);
	
	for (int j = 0; j < 3; j++) {
		
		for (int k = 0; k < 3; k++) {
			
			dColSumA[j] += fabs(dInterpMat(j,k));
			
		}
		
	}
	
	DataArray1D<double> dColSumAInv(3);
	
	dgetrf_(&m, &n, &(dInterpMat(0,0)),
			&lda, &(iPIV[0]), &info);

		
	dgetri_(&n, &(dInterpMat(0,0)),
			&lda, &(iPIV[0]), &(dWork[0]), &lWork, &info);
			

	//A inverse column sums
	for (int j = 0; j < 3; j++) {
	
		for (int k = 0; k < 3; k++) {
		
			dColSumAInv[j] += fabs(dInterpMat(j,k));
		
		}
	
	}
	
	
	//max column sums of A and A inverse
	double dMaxColSumA = dColSumA[0];
	
	double dMaxColSumAInv = dColSumAInv[0];
	
	for (int k = 1; k < 3; k++) {
		
		if (dColSumA[k] > dMaxColSumA) {
			dMaxColSumA = dColSumA[k];
		}
		if (dColSumAInv[k] > dMaxColSumAInv) {
			dMaxColSumAInv = dColSumAInv[k];
		}
	}		
	
	dCond = dMaxColSumAInv * dMaxColSumA;
		
	if (info < 0) {
		_EXCEPTION1("dgetrf_ reports matrix had an illegal value (%i)", info);
	}
	if (info > 0) {
		Announce("WARNING: Singular matrix detected in fit (likely colinear elements)");
	}
	
	//Set right hand side of linear system
	
	DataArray1D<double> dRHS(3);
	
	dRHS(0) = nodeQ.x - nodesP[0].x;
	dRHS(1) = nodeQ.y - nodesP[0].y;
	dRHS(2) = nodeQ.z - nodesP[0].z;
	
	//Solve linear system with matrix multiply
	
	MatVectorMult(dInterpMat, dRHS, dCoeffs);
	
	
}

///////////////////////////////////////////////////////////////////////////////

void NewtonQuadrilateral(
	Node & nodeQ,
	NodeVector & nodesP,
	DataArray1D<double> & dCoeffs,
	bool & fConverged
) {
	
	_ASSERT(dCoeffs.GetRows() == 3);
	
	int iMaxIterations = 100;
	
	//This algorithm is essentially the one that is used in ESMF
	
	// Four corners of quadrilateral
	Node nodeQ0 = nodesP[0];
	Node nodeQ1 = nodesP[1];
	Node nodeQ2 = nodesP[2];
	Node nodeQ3 = nodesP[3];
	
	// Jacobian
	DataArray2D<double> dJacobian(3,3);
	
	DataArray1D<double> A(3), B(3), C(3), D(3), E(3), F(3);
	
	A[0] = nodeQ0.x - nodeQ1.x + nodeQ2.x - nodeQ3.x;
	A[1] = nodeQ0.y - nodeQ1.y + nodeQ2.y - nodeQ3.y;
	A[2] = nodeQ0.z - nodeQ1.z + nodeQ2.z - nodeQ3.z;
	
	B[0] = nodeQ1.x - nodeQ0.x;
	B[1] = nodeQ1.y - nodeQ0.y;
	B[2] = nodeQ1.z - nodeQ0.z;

	C[0] = nodeQ3.x-nodeQ0.x;
	C[1] = nodeQ3.y-nodeQ0.y;
	C[2] = nodeQ3.z-nodeQ0.z;
	
	D[0] = nodeQ.x;
	D[1] = nodeQ.y;
	D[2] = nodeQ.z;
	
	E[0] = nodeQ0.x - nodeQ.x;
	E[1] = nodeQ0.y - nodeQ.y;
	E[2] = nodeQ0.z - nodeQ.z;
	
	for (int i = 0; i < iMaxIterations; i++){
		
		double dA = dCoeffs[0];
		double dB = dCoeffs[1];
		double dC = dCoeffs[2];
		
		// Calculate Value of function at X
		F[0] = dA*dB*A[0] + dA*B[0] + dB*C[0] + dC*D[0] + E[0];
		F[1] = dA*dB*A[1] + dA*B[1] + dB*C[1] + dC*D[1] + E[1];
		F[2] = dA*dB*A[2] + dA*B[2] + dB*C[2] + dC*D[2] + E[2];
		
	    // If we're close enough to 0.0 then exit
		if (F[0]*F[0]+F[1]*F[1]+F[2]*F[2] < 1.0E-15) {
			
			fConverged = true;
			break;
			
		}

		// Construct Jacobian
		dJacobian(0,0) = A[0]*dB + B[0];
		dJacobian(1,0) = A[1]*dB + B[1];
		dJacobian(2,0) = A[2]*dB + B[2];
		
		dJacobian(0,1) = A[0]*dA + C[0];
		dJacobian(1,1) = A[1]*dA + C[1];
		dJacobian(2,1) = A[2]*dA + C[2];

		dJacobian(0,2) = D[0];
		dJacobian(1,2) = D[1];
		dJacobian(2,2) = D[2];
		
		//Invert jacobian
		int m = 3;
		int n = 3;
		int lda = 3;
		int info;
		
		DataArray1D<int> iPIV;
		iPIV.Allocate(3);
		
		DataArray1D<double> dWork;
		dWork.Allocate(3);
		
		int lWork = 3;
		
		//LU decomposition
		dgetrf_(&m, &n, &(dJacobian(0,0)),
				&lda, &(iPIV[0]), &info);
	
		//Invert LU decomposition
		dgetri_(&n, &(dJacobian(0,0)),
				&lda, &(iPIV[0]), &(dWork[0]), &lWork, &info);
			
		if (info != 0) {
			_EXCEPTIONT("Mass matrix inversion error");
		}
		
		DataArray1D<double> dDeltaX(3);
		
		//Multiply solution vector by inverse jacobian
		MatVectorMult(dJacobian, F, dDeltaX);
		
		//Update solution
		dCoeffs[0] = dCoeffs[0] - dDeltaX[0];
		dCoeffs[1] = dCoeffs[1] - dDeltaX[1];
		dCoeffs[2] = dCoeffs[2] - dDeltaX[2];
		
		
	}
	
	
}

///////////////////////////////////////////////////////////////////////////////

void BuildIntegrationArray(
	const Mesh & meshInput,
	const Mesh & meshOverlap,
	const TriangularQuadratureRule & triquadrule,
	int ixFirstFace,
	int ixOverlapBegin,
	int ixOverlapEnd,
	int nOrder,
	DataArray2D<double> & dIntArray
) {
	// Number of coefficients needed at this order
#ifdef RECTANGULAR_TRUNCATION
	int nCoefficients = nOrder * nOrder;
#endif
#ifdef TRIANGULAR_TRUNCATION
	int nCoefficients = nOrder * (nOrder + 1) / 2;
#endif

	// Triangular quadrature rule
	const DataArray2D<double> & dG = triquadrule.GetG();
	const DataArray1D<double> & dW = triquadrule.GetW();

	// This Face
	const Face & faceFirst = meshInput.faces[ixFirstFace];

	// Coordinate axes
	Node nodeRef = GetFaceCentroid(faceFirst, meshInput.nodes);

#if defined(USE_STEREOGRAPHIC_FITS)
	Node nodeA1, nodeA2;
	Node nodeC = nodeRef;
	GetTangentBasis(nodeRef, nodeA1, nodeA2);
#else
	Node nodeA1 = meshInput.nodes[faceFirst[1]] - nodeRef;
	Node nodeA2 = meshInput.nodes[faceFirst[2]] - nodeRef;
	Node nodeC = CrossProduct(nodeA1, nodeA2);
#endif

	// Fit matrix
	DataArray2D<double> dFit(3,3);

	dFit(0,0) = nodeA1.x; dFit(0,1) = nodeA1.y; dFit(0,2) = nodeA1.z;
	dFit(1,0) = nodeA2.x; dFit(1,1) = nodeA2.y; dFit(1,2) = nodeA2.z;
	dFit(2,0) = nodeC.x;  dFit(2,1) = nodeC.y;  dFit(2,2) = nodeC.z;

	DataArray2D<double> dFitTemp;

	// Number of overlapping faces and triangles
	int nOverlapFaces = ixOverlapEnd - ixOverlapBegin;

	// Build integration array
	dIntArray.Allocate(nCoefficients, nOverlapFaces);

	// Loop through all overlap Faces
	for (int i = 0; i < nOverlapFaces; i++) {
		const Face &faceOverlap = meshOverlap.faces[ixOverlapBegin + i];

		const NodeVector &nodesOverlap = meshOverlap.nodes;

		int nbEdges = faceOverlap.edges.size();
		int nOverlapTriangles = 1;
		Node nodeCenter; // not used if nbEdges == 3
		if (nbEdges > 3) { // decompose from center in this case
			nOverlapTriangles = nbEdges;
			for (int k = 0; k < nbEdges; k++) {
				const Node & currNode = nodesOverlap[faceOverlap[k]];
				nodeCenter = nodeCenter + currNode;
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

			if (nbEdges == 3) {
				node0 = nodesOverlap[faceOverlap[0]];
				node1 = nodesOverlap[faceOverlap[1]];
				node2 = nodesOverlap[faceOverlap[2]];
			}
			else { // decompose polygon in triangles around the center
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

				dFitTemp = dFit;
				dgesv_(
					&n, &nrhs, &(dFitTemp(0,0)), &lda, ipiv, dX, &ldb, &info);

				// Sample this point
				int ixp = 0;

#ifdef RECTANGULAR_TRUNCATION
				for (int p = 0; p < nOrder; p++) {
				for (int q = 0; q < nOrder; q++) {
#endif
#ifdef TRIANGULAR_TRUNCATION 
				for (int p = 0; p < nOrder; p++) {
				for (int q = 0; q < nOrder - p; q++) {
#endif
					dIntArray(ixp,i) +=
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
}

///////////////////////////////////////////////////////////////////////////////

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
) {

	// Reference to active Face
	const Face & faceFirst = mesh.faces[ixFirst];

	// Number of coefficients needed at this order
#ifdef RECTANGULAR_TRUNCATION
	int nCoefficients = nOrder * nOrder;
#endif
#ifdef TRIANGULAR_TRUNCATION 
	int nCoefficients = nOrder * (nOrder + 1) / 2;
#endif

	// Number of adjacent Faces
	int nAdjFaces = vecAdjFaces.size();

	// Initialize arrays,
	dFitArray.Allocate(nCoefficients, nAdjFaces);
	dFitWeights.Allocate(nAdjFaces);

	// Triangular quadrature rule
	const DataArray2D<double> & dG = triquadrule.GetG();
	const DataArray1D<double> & dW = triquadrule.GetW();

	// Coordinate axes
	Node nodeRef = GetFaceCentroid(faceFirst, mesh.nodes);

#if defined(USE_STEREOGRAPHIC_FITS)
	Node nodeA1, nodeA2;
	Node nodeC = nodeRef;
	GetTangentBasis(nodeRef, nodeA1, nodeA2);
#else
	Node nodeA1 = mesh.nodes[faceFirst[1]] - nodeRef;
	Node nodeA2 = mesh.nodes[faceFirst[2]] - nodeRef;
	Node nodeC = CrossProduct(nodeA1, nodeA2);
#endif

	// Fit matrix
	DataArray2D<double> dFit(3,3);

	dFit(0,0) = nodeA1.x; dFit(0,1) = nodeA1.y; dFit(0,2) = nodeA1.z;
	dFit(1,0) = nodeA2.x; dFit(1,1) = nodeA2.y; dFit(1,2) = nodeA2.z;
	dFit(2,0) = nodeC.x;  dFit(2,1) = nodeC.y;  dFit(2,2) = nodeC.z;

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

//#pragma message "Pre-triangularize this matrix for efficiency"
				DataArray2D<double> dFitTemp;
				dFitTemp = dFit;
				dgesv_(
					&n, &nrhs, &(dFitTemp(0,0)), &lda, ipiv, dX, &ldb, &info);

				if (info != 0) {
					_EXCEPTIONT("Solve failure in dgesv");
				}

				// Loop through all coefficients
				int ixp = 0;

#ifdef RECTANGULAR_TRUNCATION
				for (int p = 0; p < nOrder; p++) {
				for (int q = 0; q < nOrder; q++) {
#endif
#ifdef TRIANGULAR_TRUNCATION 
				for (int p = 0; p < nOrder; p++) {
				for (int q = 0; q < nOrder - p; q++) {
#endif
					dFitArray(ixp,iAdjFace) +=
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

		// Integrate locally using the constraint
		if ((dConstraint.GetRows() != 0) && (iAdjFace == 0)) {
			for (int p = 0; p < nCoefficients; p++) {
				dFitArray(p,0) = dConstraint[p];
			}
		}

		// Reweight the fit array
		dFitWeights[iAdjFace] =
			pow(static_cast<double>(fdp.second),
			    - static_cast<double>(nFitWeightsExponent));

		for (int j = 0; j < dFitArray.GetRows(); j++) {
			dFitArray(j,iAdjFace) *= dFitWeights[iAdjFace];
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

bool InvertFitArray_Corrected(
	const DataArray1D<double> & dConstraint,
	DataArray2D<double> & dFitArray,
	DataArray1D<double> & dFitWeights,
	DataArray2D<double> & dFitArrayPlus
) {
	// Dimensions of the fit operator
	int nCoefficients = dFitArray.GetRows();
	int nAdjFaces = dFitArray.GetColumns();

	// Check parameters
	if (dConstraint.GetRows() != 0) {
	if (dConstraint.GetRows() != nCoefficients) {
		_EXCEPTIONT("Dimension mismatch between dConstraint and dFitWeights");
	}
	}

	if (dFitWeights.GetRows() != nAdjFaces) {
		_EXCEPTIONT("Dimension mismatch between dFitArray and dFitWeights");
	}

	// Allocate inverse
	dFitArrayPlus.Allocate(nAdjFaces, nCoefficients);

/*
	// Invert the fit array using Moore-Penrose pseudoinverse
	int m = nCoefficients;
	int n = nAdjFaces;
	int lda = nCoefficients;
	int ldb = nAdjFaces;
	int nrhs = nCoefficients;
	double rcond = 0.0;
	int rank = 0;
	int lwork = -1;
	int info = 0;

	DataArray1D<double> dS(nAdjFaces);

	DataArray1D<double> dWork(1);

	DataArray2D<double> dA(nAdjFaces, nCoefficients);
	for (int i = 0; i < nAdjFaces; i++) {
	for (int j = 0; j < nCoefficients; j++) {
		dA(i,j) = dFitArray(j,i);
	}
	}

	DataArray2D<double> dB(nCoefficients, nAdjFaces);
	for (int i = 0; i < nCoefficients; i++) {
		dB(i,i) = 1.0;
	}

	dgelss_(
		&m,
		&n,
		&nrhs,
		&(dA(0,0)),
		&lda,
		&(dB(0,0)),
		&ldb,
		&(dS[0]),
		&rcond,
		&rank,
		&(dWork[0]),
		&lwork,
		&info);

	if (info != 0) {
		_EXCEPTION();
	}

	lwork = static_cast<int>(dWork[0]);
	dWork.Initialize(lwork);

	dgelss_(
		&m,
		&n,
		&nrhs,
		&(dA(0,0)),
		&lda,
		&(dB(0,0)),
		&ldb,
		&(dS[0]),
		&rcond,
		&rank,
		&(dWork[0]),
		&lwork,
		&info);

	if (info != 0) {
		_EXCEPTION();
	}

	for (int i = 0; i < nAdjFaces; i++) {
	for (int j = 0; j < nCoefficients; j++) {
		dFitArrayPlus(i,j) = dB(j,i);
	}
	}
*/

	// Used for calculation of the 1-norm condition number using max column sums
	DataArray1D<double> dColSumA(nCoefficients);
	DataArray1D<double> dColSumAT(nCoefficients);

	// Compute Moore-Penrose pseudoinverse via inv(A^T*A)*A^T
	DataArray2D<double> dFit2(nCoefficients, nCoefficients);

	for (int j = 0; j < nCoefficients; j++) {
	for (int k = 0; k < nCoefficients; k++) {
	for (int l = 0; l < nAdjFaces; l++) {
		dFit2(j,k) += dFitArray(j,l) * dFitArray(k,l);
	}
	}
	}

	// Column sums of dFit2
	for (int j = 0; j < nCoefficients; j++) {
	for (int k = 0; k < nCoefficients; k++) {
		dColSumA[j] += fabs(dFit2(j,k));
	}
	}

	// Calculate pseudoinverse of FitArray
	int m = nCoefficients;
	int n = nCoefficients;
	int lda = nCoefficients;
	int info;

	DataArray1D<int> iPIV(nCoefficients);

	DataArray1D<double> dWork(nCoefficients);

	int lWork = nCoefficients;

	// Compute LU factorization
	dgetrf_(&m, &n, &(dFit2(0,0)), &lda, &(iPIV[0]), &info);
	if (info < 0) {
		_EXCEPTION1("dgetrf_ reports matrix had an illegal value (%i)", info);
	}
	if (info > 0) {
		Announce("WARNING: Singular matrix detected in fit (likely colinear elements)");
		return false;
	}

	// Matrix inverse using LU factorization
	dgetri_(&n, &(dFit2(0,0)), &lda, &(iPIV[0]), &(dWork[0]), &lWork, &info);
	if (info < 0) {
		_EXCEPTION1("dgetri_ reports matrix had an illegal value (%i)", info);
	}
	if (info > 0) {
		Announce("WARNING: Singular matrix detected in fit (likely colinear elements)");
		return false;
	}

	// Column sums of dFit2
	for (int j = 0; j < nCoefficients; j++) {
	for (int k = 0; k < nCoefficients; k++) {
		dColSumAT[j] += fabs(dFit2(j,k));
	}
	}

	// Warning if condition number estimate is too high
	double dMaxColSumA = dColSumA[0];
	double dMaxColSumAT = dColSumAT[0];
	for (int k = 1; k < nCoefficients; k++) {
		if (dColSumA[k] > dMaxColSumA) {
			dMaxColSumA = dColSumA[k];
		}
		if (dColSumAT[k] > dMaxColSumAT) {
			dMaxColSumAT = dColSumAT[k];
		}
	}
	if (dMaxColSumA * dMaxColSumAT > FVConditionNumberThreshold) {
		Announce("WARNING: Poor conditioning in fit matrix (%1.15e); dropping to 1st order",
			dMaxColSumA * dMaxColSumAT);
		return false;
	}

	// Calculate pseudoinverse
	for (int j = 0; j < nAdjFaces; j++) {
	for (int k = 0; k < nCoefficients; k++) {
	for (int l = 0; l < nCoefficients; l++) {
		dFitArrayPlus(j,k) += dFitArray(l,j) * dFit2(l,k);
	}
	}
	}

	for (int j = 0; j < nAdjFaces; j++) {
	for (int k = 0; k < nCoefficients; k++) {
		dFitArrayPlus(j,k) *= dFitWeights[j];
	}
	}

	// Apply correction to constant mode
	if (dConstraint.GetRows() != 0) {
		for (int i = 0; i < nAdjFaces; i++) {
			double dSum = 0.0;
			for (int p = 1; p < nCoefficients; p++) {
				dSum += dConstraint[p] * dFitArrayPlus(i,p);
			}
			if (i == 0) {
				dFitArrayPlus(i,0) = 1.0 - dSum;
			} else {
				dFitArrayPlus(i,0) = - dSum;
			}
		}
	}

/*
	// Apply correction to highest degree mode
	for (int i = 0; i < nAdjFaces; i++) {
		double dSum = 0.0;
		for (int p = 0; p < nCoefficients-1; p++) {
			dSum += dConstraint[p] * dFitArrayPlus(i,p);
		}
		if (i == 0) {
			dFitArrayPlus(i,nCoefficients-1) =
				(1.0 - dSum) / dConstraint[nCoefficients-1];
		} else {
			dFitArrayPlus(i,nCoefficients-1) =
				- dSum / dConstraint[nCoefficients-1];
		}
	}
*/

	return true;

}

///////////////////////////////////////////////////////////////////////////////

void InvertFitArray_LeastSquares(
	const DataArray1D<double> & dConstraint,
	DataArray2D<double> & dFitArray,
	DataArray1D<double> & dFitWeights,
	DataArray2D<double> & dFitArrayPlus
) {
	// Dimensions of the fit operator
	int nCoefficients = dFitArray.GetRows();
	int nAdjFaces = dFitArray.GetColumns();

	// Check parameters
	if (dConstraint.GetRows() != nCoefficients) {
		_EXCEPTIONT("Dimension mismatch between dConstraint and dFitWeights");
	}
	if (dFitWeights.GetRows() != nAdjFaces) {
		_EXCEPTIONT("Dimension mismatch between dFitArray and dFitWeights");
	}

	// Allocate inverse
	dFitArrayPlus.Allocate(nAdjFaces, nCoefficients);

	// Special case: First order
	if (nCoefficients == 1) {
		dFitArrayPlus(0,0) = 1.0;
		return;
	}

	// Compute QR factorization of the constraint
	DataArray2D<double> dQ(nCoefficients, nCoefficients);

	double dR;

	{
		int m = nCoefficients;
		int n = 1;
		int lda = m;

		double tau;

		DataArray1D<double> dWork(nCoefficients);
		int lwork = nCoefficients;

		int info;

		memcpy(&(dQ(0,0)), &(dConstraint[0]), nCoefficients * sizeof(double));

		dgeqrf_(&m, &n, &(dQ(0,0)), &lda, &tau, &(dWork[0]), &lwork, &info);
		if (info != 0) {
			_EXCEPTION1("Error in dgeqrf: %i", info);
		}

		dR = dQ(0,0);

		int k = 1;
		n = nCoefficients;
		dorgqr_(&m, &n, &k, &(dQ(0,0)), &lda, &tau, &(dWork[0]), &lwork, &info);
		if (info != 0) {
			_EXCEPTION1("Error in dorgqr: %i", info);
		}
	}

	// Calculate G = F * Q 
	DataArray2D<double> dGG(nCoefficients, nAdjFaces);

	for (int i = 0; i < nCoefficients; i++) {
	for (int j = 0; j < nAdjFaces; j++) {
		for (int k = 0; k < nCoefficients; k++) {
			dGG(i,j) += dFitArray(k,j) * dQ(i,k);
		}
	}
	}

	// Calculate Moore-Penrose pseudoinverse of G(:,2:p)
	DataArray2D<double> dGxPlus(nAdjFaces, nCoefficients-1);

	{
		// Gx2 = G(:, 2:p)^T * G(:,2:p)
		DataArray2D<double> dGx2(nCoefficients-1, nCoefficients-1);

		for (int i = 0; i < nCoefficients-1; i++) {
		for (int j = 0; j < nCoefficients-1; j++) {
			for (int k = 0; k < nAdjFaces; k++) {
				dGx2(i,j) += dGG(i+1,k) * dGG(j+1,k);
			}
		}
		}

		int m = nCoefficients-1;
		int n = nCoefficients-1;
		int lda = nCoefficients-1;
		int info;

		DataArray1D<int> iPIV(nCoefficients-1);

		DataArray1D<double> dWork(nCoefficients-1);

		int lWork = nCoefficients-1;

		dgetrf_(&m, &n, &(dGx2(0,0)), &lda, &(iPIV[0]), &info);

		dgetri_(&n, &(dGx2(0,0)), &lda, &(iPIV[0]), &(dWork[0]), &lWork, &info);

		// Calculate pseudoinverse
		for (int i = 0; i < nAdjFaces; i++) {
		for (int j = 0; j < nCoefficients-1; j++) {
			for (int k = 0; k < nCoefficients-1; k++) {
				dGxPlus(i,j) += dGG(k+1,i) * dGx2(k,j);
			}
		}
		}
	}

	// Z = G+ * (I - G(:,1) * R^{-1} * e0T)
	DataArray2D<double> dZ(nAdjFaces, nCoefficients-1);

	{
		DataArray2D<double> dSubZ(nAdjFaces, nAdjFaces);

		for (int i = 0; i < nAdjFaces; i++) {
			dSubZ(i,i) = 1.0;
		}

		for (int i = 0; i < nAdjFaces; i++) {
			dSubZ(0,i) -= dGG(0,i) / dR;
		}

		for (int i = 0; i < nAdjFaces; i++) {
		for (int j = 0; j < nCoefficients-1; j++) {
			for (int k = 0; k < nAdjFaces; k++) {
				dZ(i,j) += dGxPlus(k,j) * dSubZ(i,k);
			}
		}
		}
	}

	// Fhat = Q(:,1) * R^{-1} * e0T + Q(:,2:p) * Z
	for (int i = 0; i < nCoefficients; i++) {
		dFitArrayPlus(0,i) += dQ(0,i) / dR;
	}
	for (int i = 0; i < nAdjFaces; i++) {
	for (int j = 0; j < nCoefficients; j++) {
		for (int k = 0; k < nCoefficients-1; k++) {
			dFitArrayPlus(i,j) += dQ(k+1,j) * dZ(i,k);
		}
	}
	}

	for (int i = 0; i < nAdjFaces; i++) {
	for (int j = 0; j < nCoefficients; j++) {
		dFitArrayPlus(i,j) *= dFitWeights[i];
	}
	}
}

///////////////////////////////////////////////////////////////////////////////
