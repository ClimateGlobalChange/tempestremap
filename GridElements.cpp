///////////////////////////////////////////////////////////////////////////////
///
///	\file    GridElements.cpp
///	\author  Paul Ullrich
///	\version March 7, 2014
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
#include "GridElements.h"

#include "Announce.h"

#include <cmath>
#include <netcdfcpp.h>

///////////////////////////////////////////////////////////////////////////////
/// Face
///////////////////////////////////////////////////////////////////////////////

int Face::GetEdgeIndex(
	const Edge & edge
) const {
	for (int i = 0; i < edges.size(); i++) {
		if (edges[i] == edge) {
			return i;
		}
	}
	_EXCEPTIONT("Edge not found on Face");
}

///////////////////////////////////////////////////////////////////////////////

void Face::RemoveZeroEdges() {

	// Loop through all edges of this face
	for (int i = 0; i < edges.size(); i++) {

		// Remove zero edges
		if (edges[i][0] == edges[i][1]) {
			edges.erase(edges.begin()+i);
			i--;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
/// Mesh
///////////////////////////////////////////////////////////////////////////////

void Mesh::Clear() {
	nodes.clear();
	faces.clear();
	edgemap.clear();
	revnodearray.clear();
}

///////////////////////////////////////////////////////////////////////////////

void Mesh::ConstructEdgeMap() {

	// Construct the edge map
	edgemap.clear();
	for (int i = 0; i < faces.size(); i++) {
		for (int k = 0; k < 4; k++) {
			if (faces[i][k] == faces[i][(k+1)%4]) {
				continue;
			}

			Edge edge(faces[i][k], faces[i][(k+1)%4]);
			FacePair facepair;

			EdgeMapIterator iter =
				edgemap.insert(EdgeMapPair(edge, facepair)).first;

			iter->second.AddFace(i);
		}
	}

	Announce("Mesh size: Edges [%i]", edgemap.size());
}

///////////////////////////////////////////////////////////////////////////////

void Mesh::ConstructReverseNodeArray() {

	// Initialize the object
	revnodearray.resize(nodes.size());
	for (int i = 0; i < revnodearray.size(); i++) {
		revnodearray[i].clear();
	}

	// Build set for each node
	for (int i = 0; i < faces.size(); i++) {
		for (int k = 0; k < faces[i].edges.size(); k++) {
			int ixNode = faces[i].edges[k][0];
			revnodearray[ixNode].insert(i);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

Real Mesh::CalculateFaceArea(
	int iFace
) const {
	const Real Tolerance = ReferenceTolerance;

	const Face & face = faces[iFace];

	double dAccumulatedExcess = 0.0;

	int nEdges = static_cast<int>(face.edges.size());

	double dFaceArea = 0.0;

	for (int j = 0; j < nEdges; j++) {

		// Get Nodes bounding this Node
		int jPrev = (j + nEdges - 1) % nEdges;
		int jNext = (j + 1) % nEdges;

		const Node & node0 = nodes[face[jPrev]];
		const Node & node1 = nodes[face[j]];
		const Node & node2 = nodes[face[jNext]];

		// Calculate area using Karney's method
		// http://osgeo-org.1560.x6.nabble.com/Area-of-a-spherical-polygon-td3841625.html
		double dLon1 = atan2(node1.y, node1.x);
		double dLat1 = asin(node1.z);

		double dLon2 = atan2(node2.y, node2.x);
		double dLat2 = asin(node2.z);

		if ((dLon1 < -0.5 * M_PI) && (dLon2 > 0.5 * M_PI)) {
			dLon1 += 2.0 * M_PI;
		}
		if ((dLon2 < -0.5 * M_PI) && (dLon1 > 0.5 * M_PI)) {
			dLon2 += 2.0 * M_PI;
		}

		double dLamLat1 = 2.0 * atanh(tan(dLat1 / 2.0));
		double dLamLat2 = 2.0 * atanh(tan(dLat2 / 2.0));

		double dS = tan(0.5 * (dLon2 - dLon1))
			* tanh(0.5 * (dLamLat1 + dLamLat2));

		dFaceArea -= 2.0 * atan(dS);
/*
		// Calculate angle between Nodes
		Real dDotN0N1 = DotProduct(node0, node1);
		Real dDotN0N2 = DotProduct(node0, node2);
		Real dDotN2N1 = DotProduct(node2, node1);

		Real dDotN0N0 = DotProduct(node0, node0);
		Real dDotN1N1 = DotProduct(node1, node1);
		Real dDotN2N2 = DotProduct(node2, node2);

		Real dDenom0 = dDotN0N0 * dDotN1N1 - dDotN0N1 * dDotN0N1;
		Real dDenom2 = dDotN2N2 * dDotN1N1 - dDotN2N1 * dDotN2N1;

		if ((dDenom0 < Tolerance) || (dDenom2 < Tolerance)) {
			continue;
		}

		// Angle between nodes
		double dDenom = sqrt(dDenom0 * dDenom2);

		double dRatio =
			(dDotN1N1 * dDotN0N2 - dDotN0N1 * dDotN2N1) / dDenom;

		if (dRatio > 1.0) {
			dRatio = 1.0;
		} else if (dRatio < -1.0) {
			dRatio = -1.0;
		}

		dAccumulatedExcess += acos(dRatio);
*/
		//printf("[%1.15e %1.15e %1.15e],\n", node1.x, node1.y, node1.z);
		//printf("%1.15e %1.15e\n",  dDotN1N1, dDotN2N2);

	}

	//double dPlanarSum = static_cast<double>(nEdges-2) * M_PI;

	if (dFaceArea < 0.0) {
		dFaceArea += 2.0 * M_PI;
		//printf("%1.15e %1.15e\n", dFaceArea, dAccumulatedExcess - dPlanarSum);
		//_EXCEPTIONT("Negative area element detected");
	}
/*
	if (dPlanarSum > dAccumulatedExcess) {
		printf("%1.15e %1.15e\n", dPlanarSum, dAccumulatedExcess);
		_EXCEPTIONT("Negative area element detected");
	}
*/
	return dFaceArea;
}

///////////////////////////////////////////////////////////////////////////////

Real Mesh::CalculateFaceAreas() {

	// Calculate the area of each Face
	vecFaceArea.Initialize(faces.size());
	for (int i = 0; i < faces.size(); i++) {
		vecFaceArea[i] = CalculateFaceArea(i);
	}

	// Calculate accumulated area carefully
	static const int Jump = 10;
	std::vector<double> vecFaceAreaBak;
	vecFaceAreaBak.resize(vecFaceArea.GetRows());
	memcpy(&(vecFaceAreaBak[0]), &(vecFaceArea[0]),
		vecFaceArea.GetRows() * sizeof(double));

	for (;;) {
		if (vecFaceAreaBak.size() == 1) {
			break;
		}
		for (int i = 0; i <= (vecFaceAreaBak.size()-1) / Jump; i++) {
			int ixRef = Jump * i;
			vecFaceAreaBak[i] = vecFaceAreaBak[ixRef];
			for (int j = 1; j < Jump; j++) {
				if (ixRef + j >= vecFaceAreaBak.size()) {
					break;
				}
				vecFaceAreaBak[i] += vecFaceAreaBak[ixRef + j];
			}
		}
		vecFaceAreaBak.resize((vecFaceAreaBak.size()-1) / Jump + 1);
	}

	return vecFaceAreaBak[0];
}

///////////////////////////////////////////////////////////////////////////////

void Mesh::Write(const std::string & strFile) const {
	const int ParamFour = 4;
	const int ParamLenString = 33;

	// Determine the maximum number of nodes per element
	int nNodesPerElement = 0;
	for (int i = 0; i < faces.size(); i++) {
		if (faces[i].edges.size() > nNodesPerElement) {
			nNodesPerElement = faces[i].edges.size();
		}
	}
	Announce("Max nodes per element: %i", nNodesPerElement);

	// Output to a NetCDF Exodus file
	NcFile ncOut(strFile.c_str(), NcFile::Replace);

	// Random Exodus dimensions
	NcDim * dimLenString = ncOut.add_dim("len_string", ParamLenString);
	NcDim * dimLenLine = ncOut.add_dim("len_line", 81);
	NcDim * dimFour = ncOut.add_dim("four", ParamFour);
	NcDim * dimTime = ncOut.add_dim("time_step");
	NcDim * dimDimension = ncOut.add_dim("num_dim", 3);

	// Number of nodes
	int nNodeCount = nodes.size();
	NcDim * dimNodes = ncOut.add_dim("num_nodes", nNodeCount);

	// Number of elements
	int nElementCount = faces.size();
	NcDim * dimElements = ncOut.add_dim("num_elem", nElementCount);

	// Other dimensions
	NcDim * dimNumElementBlocks = ncOut.add_dim("num_el_blk", 1);
	NcDim * dimNumQARec = ncOut.add_dim("num_qa_rec", 1);
	NcDim * dimElementBlock1 = ncOut.add_dim("num_el_in_blk1", nElementCount);
	NcDim * dimNodesPerElement =
		ncOut.add_dim("num_nod_per_el1", nNodesPerElement);
	NcDim * dimAttBlock1 = ncOut.add_dim("num_att_in_blk1", 1);

	// Global attributes
	ncOut.add_att("api_version", 4.98f);
	ncOut.add_att("version", 4.98f);
	ncOut.add_att("floating_point_word_size", 8);
	ncOut.add_att("file_size", 0);

	char szTitle[128];
	sprintf(szTitle, "tempest(%s) 01/01/2013: 00:00:00", strFile.c_str());
	ncOut.add_att("title", szTitle);

	// Time_whole (unused)
	ncOut.add_var("time_whole", ncDouble, dimTime);

	// QA records
	char szQARecord[ParamFour][ParamLenString] = {
		"Tempest", "13.0", "01/01/2013", "00:00:00"};

	NcVar * varQARecords =
		ncOut.add_var("qa_records", ncChar, dimNumQARec, dimFour, dimLenString);
	varQARecords->set_cur(0, 0, 0);
	varQARecords->put(&(szQARecord[0][0]), 1, 4, ParamLenString);

	// Coordinate names
	char szCoordNames[3][ParamLenString] = {"x", "y", "z"};
	
	NcVar * varCoordNames =
		ncOut.add_var("coor_names", ncChar, dimDimension, dimLenString);
	varCoordNames->set_cur(0, 0, 0);
	varCoordNames->put(&(szCoordNames[0][0]), 3, ParamLenString);

	// Element block names
	NcVar * varElementBlockNames =
		ncOut.add_var("eb_names", ncChar, dimNumElementBlocks, dimLenString);

	// Element map
	int * nElementMap = new int[nElementCount];
	for (int i = 0; i < nElementCount; i++) {
		nElementMap[i] = i+1;
	}

	NcVar * varElementMap =
		ncOut.add_var("elem_map", ncInt, dimElements);
	varElementMap->put(nElementMap, nElementCount);

	delete[] nElementMap;

	// Element block status
	int nOne = 1;

	NcVar * varElementBlockStatus =
		ncOut.add_var("eb_status", ncInt, dimNumElementBlocks);
	varElementBlockStatus->put(&nOne, 1);

	NcVar * varElementProperty =
		ncOut.add_var("eb_prop1", ncInt, dimNumElementBlocks);
	varElementProperty->put(&nOne, 1);
	varElementProperty->add_att("name", "ID");

	// Attributes
	double * dAttrib1 = new double[nElementCount];
	for (int i = 0; i < nElementCount; i++) {
		dAttrib1[i] = 1.0;
	}

	NcVar * varAttrib1 =
		ncOut.add_var("attrib1", ncDouble, dimElementBlock1, dimAttBlock1);
	varAttrib1->put(dAttrib1, nElementCount, 1);
	delete[] dAttrib1;

	// Face nodes (1-indexed)
	NcVar * varFaces =
		ncOut.add_var("connect1", ncInt, dimElementBlock1, dimNodesPerElement);

	varFaces->add_att("elem_type", "SHELL4");

	int * nConnect = new int[nNodesPerElement];
	for (int i = 0; i < nElementCount; i++) {
		int nEdges = faces[i].edges.size();
		int k = 0;
		for (; k < nEdges; k++) {
			nConnect[k] = faces[i][k] + 1;
		}
		for (; k < nNodesPerElement; k++) {
			nConnect[k] = nConnect[nEdges-1];
		}

		varFaces->set_cur(i, 0);
		varFaces->put(nConnect, 1, nNodesPerElement);
	}

	delete[] nConnect;

	// Node list
	NcVar * varNodes =
		ncOut.add_var("coord", ncDouble, dimDimension, dimNodes);

	double * dCoord = new double[nNodeCount];
	for (int i = 0; i < nNodeCount; i++) {
		dCoord[i] = static_cast<double>(nodes[i].x);
	}
	varNodes->set_cur(0, 0);
	varNodes->put(dCoord, 1, nNodeCount);
	for (int i = 0; i < nNodeCount; i++) {
		dCoord[i] = static_cast<double>(nodes[i].y);
	}
	varNodes->set_cur(1, 0);
	varNodes->put(dCoord, 1, nNodeCount);
	for (int i = 0; i < nNodeCount; i++) {
		dCoord[i] = static_cast<double>(nodes[i].z);
	}
	varNodes->set_cur(2, 0);
	varNodes->put(dCoord, 1, nNodeCount);
	delete[] dCoord;

	// Edge types
	NcVar * varEdgeTypes =
		ncOut.add_var("edge_type", ncInt,
			dimElementBlock1, dimNodesPerElement);

	int * nEdgeType = new int[nNodesPerElement];
	for (int i = 0; i < nElementCount; i++) {
		for (int k = 0; k < nNodesPerElement; k++) {
			nEdgeType[k] = static_cast<int>(faces[i].edges[k].type);
		}
		varEdgeTypes->set_cur(i, 0);
		varEdgeTypes->put(nEdgeType, 1, nNodesPerElement);
	}
	delete[] nEdgeType;

	// Source elements from mesh 1
	if (vecFirstFaceIx.size() != 0) {
		if (vecFirstFaceIx.size() != nElementCount) {
			_EXCEPTIONT("Incorrect size of vecFirstFaceIx");
		}

		NcVar * varFirstMeshSourceFace =
			ncOut.add_var("face_source_1", ncInt, dimElementBlock1);

		varFirstMeshSourceFace->set_cur((long)0);
		varFirstMeshSourceFace->put(&(vecFirstFaceIx[0]), nElementCount);
	}

	// Source elements from mesh 2
	if (vecSecondFaceIx.size() != 0) {
		if (vecSecondFaceIx.size() != nElementCount) {
			_EXCEPTIONT("Incorrect size of vecSecondFaceIx");
		}

		NcVar * varSecondMeshSourceFace =
			ncOut.add_var("face_source_2", ncInt, dimElementBlock1);

		varSecondMeshSourceFace->set_cur((long)0);
		varSecondMeshSourceFace->put(&(vecSecondFaceIx[0]), nElementCount);
	}
}

///////////////////////////////////////////////////////////////////////////////

void Mesh::Read(const std::string & strFile) {

	// Input from a NetCDF Exodus file
	NcFile ncFile(strFile.c_str(), NcFile::ReadOnly);

	// Determine number of nodes per element
	NcDim * dimNodesPerElement = ncFile.get_dim("num_nod_per_el1");
	int nNodesPerElement = dimNodesPerElement->size();

	// Number of nodes
	NcDim * dimNodes = ncFile.get_dim("num_nodes");
	int nNodeCount = dimNodes->size();

	// Number of elements
	NcDim * dimElements = ncFile.get_dim("num_elem");
	int nElementCount = dimElements->size();

	// Output size
	Announce("Mesh size: Nodes [%i] Elements [%i]", nNodeCount, nElementCount);

	// Load in node array
	nodes.resize(nNodeCount);

	NcVar * varNodes = ncFile.get_var("coord");

	double dCoord[3];
	for (int i = 0; i < nNodeCount; i++) {
		varNodes->set_cur(0, i);
		varNodes->get(dCoord, 3, 1);
		nodes[i].x = static_cast<Real>(dCoord[0]);
		nodes[i].y = static_cast<Real>(dCoord[1]);
		nodes[i].z = static_cast<Real>(dCoord[2]);

#ifdef USE_EXACT_ARITHMETIC
		nodes[i].fx.Set(nodes[i].x);
		nodes[i].fy.Set(nodes[i].y);
		nodes[i].fz.Set(nodes[i].z);
/*
		printf("%1.15e : ", nodes[i].x); nodes[i].fx.Print(); printf("\n");
		printf("%1.15e : ", nodes[i].y); nodes[i].fy.Print(); printf("\n");
		printf("%1.15e : ", nodes[i].z); nodes[i].fz.Print(); printf("\n");
*/
#endif
	}

	// Load in face array
	faces.resize(nElementCount, Face(nNodesPerElement));

	NcVar * varFaces = ncFile.get_var("connect1");

	int * nNodes = new int[nNodesPerElement];
	for (int i = 0; i < nElementCount; i++) {
		varFaces->set_cur(i, 0);
		varFaces->get(nNodes, 1, nNodesPerElement);

		for (int j = 0; j < nNodesPerElement; j++) {
			faces[i].SetNode(j, nNodes[j]-1);
		}
	}
	delete[] nNodes;

	// Check for variables
	bool fHasEdgeType = false;
	bool fHasFirstMeshSourceFace = false;
	bool fHasSecondMeshSourceFace = false;

	for (int v = 0; v < ncFile.num_vars(); v++) {
		if (strcmp(ncFile.get_var(v)->name(), "edge_type") == 0) {
			fHasEdgeType = true;
		}
		if (strcmp(ncFile.get_var(v)->name(), "face_source_1") == 0) {
			fHasFirstMeshSourceFace = true;
		}
		if (strcmp(ncFile.get_var(v)->name(), "face_source_2") == 0) {
			fHasSecondMeshSourceFace = true;
		}
	}


	// Load in edge type array
	if (fHasEdgeType) {
		NcVar * varEdgeTypes = ncFile.get_var("edge_type");

		int * nEdgeTypes = new int[nNodesPerElement];
		for (int i = 0; i < nElementCount; i++) {
			varEdgeTypes->set_cur(i, 0);
			varEdgeTypes->get(nEdgeTypes, 1, nNodesPerElement);
			faces[i].edges[0].type = static_cast<Edge::Type>(nEdgeTypes[0]);
			faces[i].edges[1].type = static_cast<Edge::Type>(nEdgeTypes[1]);
			faces[i].edges[2].type = static_cast<Edge::Type>(nEdgeTypes[2]);
			faces[i].edges[3].type = static_cast<Edge::Type>(nEdgeTypes[3]);
		}
		delete[] nEdgeTypes;
	}

	// Load in first mesh source face ix
	if (fHasFirstMeshSourceFace) {
		NcVar * varFaceSource1 = ncFile.get_var("face_source_1");
		vecFirstFaceIx.resize(nElementCount);
		varFaceSource1->set_cur((long)0);
		varFaceSource1->get(&(vecFirstFaceIx[0]), nElementCount);
	}

	// Load in second mesh source face ix
	if (fHasSecondMeshSourceFace) {
		NcVar * varFaceSource2 = ncFile.get_var("face_source_2");
		vecSecondFaceIx.resize(nElementCount);
		varFaceSource2->set_cur((long)0);
		varFaceSource2->get(&(vecSecondFaceIx[0]), nElementCount);
	}
}

///////////////////////////////////////////////////////////////////////////////

void Mesh::RemoveZeroEdges() {

	// Remove zero edges from all Faces
	for (int i = 0; i < faces.size(); i++) {
		faces[i].RemoveZeroEdges();
	}
}

///////////////////////////////////////////////////////////////////////////////

void Mesh::Validate() const {

	// Validate that edges are oriented counter-clockwise
	for (int i = 0; i < faces.size(); i++) {
		const Face & face = faces[i];

		const int nEdges = face.edges.size();

		for (int j = 0; j < nEdges; j++) {

			// Check for zero edges
			for(;;) {
				if (face.edges[j][0] == face.edges[j][1]) {
					j++;
				} else {
					break;
				}
				if (j == nEdges) {
					break;
				}
			}

			if (j == nEdges) {
				break;
			}

			// Find the next non-zero edge
			int jNext = (j + 1) % nEdges;

			for(;;) {
				if (face.edges[jNext][0] == face.edges[jNext][1]) {
					jNext++;
				} else {
					break;
				}
				if (jNext == nEdges) {
					jNext = 0;
				}
				if (jNext == ((j + 1) % nEdges)) {
					_EXCEPTIONT("Mesh validation failed: "
						"No edge information on Face");
				}
			}

			// Get edges
			const Edge & edge0 = face.edges[j];
			const Edge & edge1 = face.edges[(j + 1) % nEdges];

			if (edge0[1] != edge1[0]) {
				_EXCEPTIONT("Mesh validation failed: Edge cyclicity error");
			}

			const Node & node0 = nodes[edge0[0]];
			const Node & node1 = nodes[edge0[1]];
			const Node & node2 = nodes[edge1[1]];

			// Vectors along edges
			Node nodeD1 = node0 - node1;
			Node nodeD2 = node2 - node1;
/*
			node0.fx.Print(); printf("\n");
			node1.fx.Print(); printf("\n");
			nodeD1.fx.Print(); printf("\n");
			_EXCEPTION();
*/
			// Compute cross-product
			Node nodeCross(CrossProductIX(nodeD1, nodeD2));

			// Dot cross product with radial vector
			Real dDot = DotProduct(node1, nodeCross);

#ifdef USE_EXACT_ARITHMETIC
			FixedPoint dDotX = DotProductX(node1, nodeCross);

			printf("%1.15e : ", nodeCross.x); nodeCross.fx.Print(); printf("\n");

			if (fabs(nodeCross.x - nodeCross.fx.ToReal()) > ReferenceTolerance) {
				printf("X0: %1.15e : ", node0.x); node0.fx.Print(); printf("\n");
				printf("Y0: %1.15e : ", node0.y); node0.fy.Print(); printf("\n");
				printf("Z0: %1.15e : ", node0.z); node0.fz.Print(); printf("\n");
				printf("X1: %1.15e : ", node1.x); node1.fx.Print(); printf("\n");
				printf("Y1: %1.15e : ", node1.y); node1.fy.Print(); printf("\n");
				printf("Z1: %1.15e : ", node1.z); node1.fz.Print(); printf("\n");
				printf("X2: %1.15e : ", node2.x); node2.fx.Print(); printf("\n");
				printf("Y2: %1.15e : ", node2.y); node2.fy.Print(); printf("\n");
				printf("Z2: %1.15e : ", node2.z); node2.fz.Print(); printf("\n");

				printf("X1: %1.15e : ", nodeD1.x); nodeD1.fx.Print(); printf("\n");
				printf("Y1: %1.15e : ", nodeD1.y); nodeD1.fy.Print(); printf("\n");
				printf("Z1: %1.15e : ", nodeD1.z); nodeD1.fz.Print(); printf("\n");
				printf("X2: %1.15e : ", nodeD2.x); nodeD2.fx.Print(); printf("\n");
				printf("Y2: %1.15e : ", nodeD2.y); nodeD2.fy.Print(); printf("\n");
				printf("Z2: %1.15e : ", nodeD2.z); nodeD2.fz.Print(); printf("\n");
				_EXCEPTIONT("FixedPoint mismatch (X)");
			}
			if (fabs(nodeCross.y - nodeCross.fy.ToReal()) > ReferenceTolerance) {
				_EXCEPTIONT("FixedPoint mismatch (Y)");
			}
			if (fabs(nodeCross.z - nodeCross.fz.ToReal()) > ReferenceTolerance) {
				_EXCEPTIONT("FixedPoint mismatch (Z)");
			}

#endif

			if (dDot > 0.0) {
				printf("\nError detected (orientation):\n");
				printf("  Face %i, Edge %i, Orientation %1.5e\n",
					i, j, dDot);

				printf("  (x,y,z):\n");
				printf("    n0: %1.5e %1.5e %1.5e\n", node0.x, node0.y, node0.z);
				printf("    n1: %1.5e %1.5e %1.5e\n", node1.x, node1.y, node1.z);
				printf("    n2: %1.5e %1.5e %1.5e\n", node2.x, node2.y, node2.z);

				Real dR0 = sqrt(
					node0.x * node0.x + node0.y * node0.y + node0.z * node0.z);
				Real dLat0 = asin(node0.z / dR0);
				Real dLon0 = atan2(node0.y, node0.x);

				Real dR1 = sqrt(
					node1.x * node1.x + node1.y * node1.y + node1.z * node1.z);
				Real dLat1 = asin(node1.z / dR1);
				Real dLon1 = atan2(node1.y, node1.x);

				Real dR2 = sqrt(
					node2.x * node2.x + node2.y * node2.y + node2.z * node2.z);
				Real dLat2 = asin(node2.z / dR2);
				Real dLon2 = atan2(node2.y, node2.x);

				printf("  (lambda, phi):\n");
				printf("    n0: %1.5e %1.5e\n", dLon0, dLat0);
				printf("    n1: %1.5e %1.5e\n", dLon1, dLat1);
				printf("    n2: %1.5e %1.5e\n", dLon2, dLat2);

				printf("  X-Product:\n");
				printf("    %1.5e %1.5e %1.5e\n",
					nodeCross.x, nodeCross.y, nodeCross.z);

				_EXCEPTIONT(
					"Mesh validation failed: Clockwise element detected");
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
// General purpose functions
///////////////////////////////////////////////////////////////////////////////

bool IsPositivelyOrientedEdge(
	const Node & nodeBegin,
	const Node & nodeEnd
) {
	const Real Tolerance = ReferenceTolerance;

	if ((fabs(nodeBegin.x - nodeEnd.x) < Tolerance) &&
		(fabs(nodeBegin.y - nodeEnd.y) < Tolerance) &&
		(fabs(nodeBegin.z - nodeEnd.z) < Tolerance)
	) {
		_EXCEPTIONT("Latitude line of zero length");
	}

	// Both nodes in positive y half-plane
	if ((nodeBegin.y >= 0.0) && (nodeEnd.y >= 0.0)) {
		if (nodeEnd.x < nodeBegin.x) {
			return true;
		} else {
			return false;
		}

	// Both nodes in negative y half-plane
	} else if ((nodeBegin.y <= 0.0) && (nodeEnd.y <= 0.0)) {
		if (nodeEnd.x > nodeBegin.x) {
			return true;
		} else {
			return false;
		}

	// Both nodes in positive x half-plane
	} else if ((nodeBegin.x >= 0.0) && (nodeEnd.x >= 0.0)) {
		if (nodeEnd.y > nodeBegin.y) {
			return true;
		} else {
			return false;
		}

	// Both nodes in negative x half-plane
	} else if ((nodeBegin.x <= 0.0) && (nodeEnd.x <= 0.0)) {
		if (nodeEnd.y < nodeBegin.y) {
			return true;
		} else {
			return false;
		}

	// Arc length too large
	} else {
		_EXCEPTIONT("Arc length too large to determine orientation.");
	}
}

///////////////////////////////////////////////////////////////////////////////

void GetLocalDirection(
	const Node & nodeBegin,
	const Node & nodeEnd,
	const Node & nodeRef,
	const Edge::Type edgetype,
	Node & nodeDir
) {

	// Direction along a great circle arc
	if (edgetype == Edge::Type_GreatCircleArc) {

		// Cartesian direction
		nodeDir = nodeEnd - nodeBegin;

		// Project onto surface of the sphere
		Real dDotDirBegin = DotProduct(nodeDir, nodeRef);
		Real dNormNodeBegin = DotProduct(nodeRef, nodeRef);

		nodeDir.x -= dDotDirBegin / dNormNodeBegin * nodeRef.x;
		nodeDir.y -= dDotDirBegin / dNormNodeBegin * nodeRef.y;
		nodeDir.z -= dDotDirBegin / dNormNodeBegin * nodeRef.z;

	// Direction along a line of constant latitude
	} else if (edgetype == Edge::Type_ConstantLatitude) {
		nodeDir.z = 0.0;

		if (IsPositivelyOrientedEdge(nodeBegin, nodeEnd)) {
			nodeDir.x = - nodeBegin.y;
			nodeDir.y = + nodeBegin.x;
		} else {
			nodeDir.x = + nodeBegin.y;
			nodeDir.y = - nodeBegin.x;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void GetLocalDirection(
	const Node & nodeBegin,
	const Node & nodeEnd,
	const Edge::Type edgetype,
	Node & nodeDir
) {

	// Direction along a great circle arc
	if (edgetype == Edge::Type_GreatCircleArc) {

		// Cartesian direction
		nodeDir = nodeEnd - nodeBegin;

		// Project onto surface of the sphere
		Real dDotDirBegin   = DotProduct(nodeDir, nodeBegin);
		Real dNormNodeBegin = DotProduct(nodeBegin, nodeBegin);

		nodeDir.x -= dDotDirBegin / dNormNodeBegin * nodeBegin.x;
		nodeDir.y -= dDotDirBegin / dNormNodeBegin * nodeBegin.y;
		nodeDir.z -= dDotDirBegin / dNormNodeBegin * nodeBegin.z;

	// Direction along a line of constant latitude
	} else if (edgetype == Edge::Type_ConstantLatitude) {
		nodeDir.z = 0.0;

		if (IsPositivelyOrientedEdge(nodeBegin, nodeEnd)) {
			nodeDir.x = - nodeBegin.y;
			nodeDir.y = + nodeBegin.x;
		} else {
			nodeDir.x = + nodeBegin.y;
			nodeDir.y = - nodeBegin.x;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void EqualizeCoincidentNodes(
	const Mesh & meshFirst,
	Mesh & meshSecond
) {
	int nCoincidentNodes = 0;

	// Sort nodes
	std::map<Node, int> setSortedFirstNodes;
	for (int i = 0; i < meshFirst.nodes.size(); i++) {
		setSortedFirstNodes.insert(
			std::pair<Node, int>(meshFirst.nodes[i], i));
	}

	// For each node in meshSecond determine if a corresponding node
	// exists in meshFirst.
	for (int i = 0; i < meshSecond.nodes.size(); i++) {
		std::map<Node, int>::const_iterator iter =
			setSortedFirstNodes.find(meshSecond.nodes[i]);

		if (iter != setSortedFirstNodes.end()) {
			meshSecond.nodes[i] = iter->first;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void EqualizeCoincidentNodes(
	Mesh & mesh
) {
	int nCoincidentNodes = 0;

	// Sort nodes
	std::map<Node, int> mapSortedNodes;
	for (int i = 0; i < mesh.nodes.size(); i++) {
		std::map<Node, int>::const_iterator iter =
			mapSortedNodes.find(mesh.nodes[i]);

		if (iter != mapSortedNodes.end()) {
			nCoincidentNodes++;
			mesh.nodes[i] = iter->first;
		} else {
			mapSortedNodes.insert(
				std::pair<Node, int>(mesh.nodes[i], i));
		}
	}

	printf("Coincident nodes: %i\n", nCoincidentNodes);
}

///////////////////////////////////////////////////////////////////////////////

int BuildCoincidentNodeVector(
	const Mesh & meshFirst,
	const Mesh & meshSecond,
	std::vector<int> & vecSecondToFirstCoincident
) {
	int nCoincidentNodes = 0;

	// Sort nodes
	std::map<Node, int> setSortedFirstNodes;
	for (int i = 0; i < meshFirst.nodes.size(); i++) {
		setSortedFirstNodes.insert(
			std::pair<Node, int>(meshFirst.nodes[i], i));
	}

	// Resize array
	vecSecondToFirstCoincident.resize(meshSecond.nodes.size(), InvalidNode);

	// For each node in meshSecond determine if a corresponding node
	// exists in meshFirst.
	for (int i = 0; i < meshSecond.nodes.size(); i++) {
		std::map<Node, int>::const_iterator iter =
			setSortedFirstNodes.find(meshSecond.nodes[i]);

		if (iter != setSortedFirstNodes.end()) {
			vecSecondToFirstCoincident[i] = iter->second;
			nCoincidentNodes++;
		}
	}

	return nCoincidentNodes;
}

///////////////////////////////////////////////////////////////////////////////

