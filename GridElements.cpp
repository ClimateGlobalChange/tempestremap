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

void Face::ContainsNode(
	const NodeVector & nodevec,
	const Node & node,
	Face::NodeLocation & loc,
	int & ixLocation
) const {
	static const Real Tolerance = ReferenceTolerance;

	// Set of edges which "contain" this node
	std::set<int> setContainedEdgeIx;

	// Loop through all edges of this face
	for (int i = 0; i < edges.size(); i++) {

		// Paired edges means edge is invalid
		if (edges[i][0] == edges[i][1]) {
			_EXCEPTIONT("Zero Edge detected");
		}

		// Check which side of the face this edge is on
		const Node & na = nodevec[edges[i][0]];
		const Node & nb = nodevec[edges[i][1]];

		if (edges[i].type == Edge::Type_GreatCircleArc) {
			Real dDotNorm = DotProduct(CrossProduct(na, nb), node);

			//printf("Norm1: %i %1.5e\n", i, dDotNorm);

			if (dDotNorm < - Tolerance) {
				loc = NodeLocation_Exterior;
				ixLocation = 0;
				return;
			}
			if (dDotNorm < Tolerance) {
				setContainedEdgeIx.insert(i);
			}

		} else if (edges[i].type == Edge::Type_ConstantLatitude) {
			Real dAlignment = (na.x * nb.y - nb.x * na.y);
			Real dDotNorm = dAlignment / fabs(dAlignment) * (node.z - na.z);

			//printf("Norm2: %i %1.5e %1.5e %1.5e %1.5e\n", i, dAlignment, node.z, na.z, dDotNorm);

			if (dDotNorm < - Tolerance) {
				loc = NodeLocation_Exterior;
				ixLocation = 0;
				return;
			}
			if (dDotNorm < Tolerance) {
				setContainedEdgeIx.insert(i);
			}

		} else {
			_EXCEPTIONT("Invalid EdgeType");
		}
	}

	// Check if the node is contained on an edge
	if (setContainedEdgeIx.size() == 1) {
		loc = NodeLocation_Edge;
		ixLocation = *(setContainedEdgeIx.begin());
		return;
	}

	// Node is coincident with a corner of this face
	if (setContainedEdgeIx.size() == 2) {

		std::set<int>::iterator iter;

		iter = setContainedEdgeIx.begin();
		int ix0 = *(iter);

		iter++;
		int ix1 = *(iter);

		if ((ix0 == 0) && (ix1 != 1)) {
			ixLocation = 0;
		} else {
			ixLocation = ix1;
		}

		loc = NodeLocation_Corner;
		return;
	}

	// Node occurs in more than two edges; error.
	if (setContainedEdgeIx.size() > 2) {
		_EXCEPTIONT("Logic error: Node occurs in more than two edges");
	}

	// Default; node occurs in the interior of the face
	loc = NodeLocation_Interior;
	ixLocation = 0;
	return;
}

///////////////////////////////////////////////////////////////////////////////

#ifdef USE_EXACT_ARITHMETIC

void Face::ContainsNodeX(
	const NodeVector & nodevec,
	const Node & node,
	Face::NodeLocation & loc,
	int & ixLocation
) const {

	// Set of edges which "contain" this node
	std::set<int> setContainedEdgeIx;

	// Loop through all Edges of this face
	for (int i = 0; i < edges.size(); i++) {

		// Paired Edges means Edge is invalid
		if (edges[i][0] == edges[i][1]) {
			_EXCEPTIONT("Zero Edge detected");
		}

		// Check which side of the Face this Edge is on
		const Node & na = nodevec[edges[i][0]];
		const Node & nb = nodevec[edges[i][1]];

		if (edges[i].type == Edge::Type_GreatCircleArc) {
			FixedPoint fpDotNorm = DotProductX(CrossProductX(na, nb), node);

			if (fpDotNorm.IsNegative()) {
				loc = NodeLocation_Exterior;
				ixLocation = 0;
				return;
			}
			if (fpDotNorm.IsZero()) {
				setContainedEdgeIx.insert(i);
			}

		} else if (edges[i].type == Edge::Type_ConstantLatitude) {
			_EXCEPTIONT("Unimplemented");

		} else {
			_EXCEPTIONT("Invalid EdgeType");
		}
	}

	// Check if the node is contained on an edge
	if (setContainedEdgeIx.size() == 1) {
		loc = NodeLocation_Edge;
		ixLocation = *(setContainedEdgeIx.begin());
		return;
	}

	// Node is coincident with a corner of this face
	if (setContainedEdgeIx.size() == 2) {

		std::set<int>::iterator iter;

		iter = setContainedEdgeIx.begin();
		int ix0 = *(iter);

		iter++;
		int ix1 = *(iter);

		if ((ix0 == 0) && (ix1 != 1)) {
			ixLocation = 0;
		} else {
			ixLocation = ix1;
		}

		loc = NodeLocation_Corner;
		return;
	}

	// Node occurs in more than two edges; error.
	if (setContainedEdgeIx.size() > 2) {
		_EXCEPTIONT("Logic error: Node occurs in more than two edges");
	}

	// Default; node occurs in the interior of the face
	loc = NodeLocation_Interior;
	ixLocation = 0;
	return;
}

#endif

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
		printf("%1.15Le : ", nodes[i].x); nodes[i].fx.Print(); printf("\n");
		printf("%1.15Le : ", nodes[i].y); nodes[i].fy.Print(); printf("\n");
		printf("%1.15Le : ", nodes[i].z); nodes[i].fz.Print(); printf("\n");
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

	// Load in edge type array
	bool fHasEdgeType = false;

	for (int v = 0; v < ncFile.num_vars(); v++) {
		if (strcmp(ncFile.get_var(v)->name(), "edge_type") == 0) {
			fHasEdgeType = true;
			break;
		}
	}

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
			Node vecCross(CrossProductIX(nodeD1, nodeD2));

			// Dot cross product with radial vector
			Real dDot = DotProduct(node1, vecCross);

#ifdef USE_EXACT_ARITHMETIC
			FixedPoint dDotX = DotProductX(node1, vecCross);
#endif
			printf("%1.15Le : ", vecCross.x); vecCross.fx.Print(); printf("\n");

			if (dDot > 0.0) {
				printf("\nError detected (orientation):\n");
				printf("  Face %i, Edge %i, Orientation %1.5Le\n",
					i, j, dDot);

				printf("  (x,y,z):\n");
				printf("    n0: %1.5Le %1.5Le %1.5Le\n", node0.x, node0.y, node0.z);
				printf("    n1: %1.5Le %1.5Le %1.5Le\n", node1.x, node1.y, node1.z);
				printf("    n2: %1.5Le %1.5Le %1.5Le\n", node2.x, node2.y, node2.z);

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
				printf("    n0: %1.5Le %1.5Le\n", dLon0, dLat0);
				printf("    n1: %1.5Le %1.5Le\n", dLon1, dLat1);
				printf("    n2: %1.5Le %1.5Le\n", dLon2, dLat2);

				printf("  X-Product:\n");
				printf("    %1.5Le %1.5Le %1.5Le\n",
					vecCross.x, vecCross.y, vecCross.z);

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

void NudgeAlongEdge(
	const Node & nodeBegin,
	const Node & nodeEnd,
	const Edge::Type type,
	Node & nodeNudged,
	Real Nudge
) {
	//static const Real Nudge = 1.0e-6;

	Node nodeDelta(
		nodeEnd.x - nodeBegin.x,
		nodeEnd.y - nodeBegin.y,
		nodeEnd.z - nodeBegin.z);

	Real dModNudge = Nudge / nodeDelta.Magnitude();

	if (fabs(dModNudge) < 1.0e-12) {
		_EXCEPTIONT("Coincident Begin and End nodes");
	}

	// Edge is a great circle arc
	if (type == Edge::Type_GreatCircleArc) {

		nodeNudged.x = nodeBegin.x * (1.0 - dModNudge) + dModNudge * nodeEnd.x;
		nodeNudged.y = nodeBegin.y * (1.0 - dModNudge) + dModNudge * nodeEnd.y;
		nodeNudged.z = nodeBegin.z * (1.0 - dModNudge) + dModNudge * nodeEnd.z;

		Real dAbsNodeNudge = nodeNudged.Magnitude();

		nodeNudged.x /= dAbsNodeNudge;
		nodeNudged.y /= dAbsNodeNudge;
		nodeNudged.z /= dAbsNodeNudge;

	// Edge is a line of constant latitude
	} else if (type == Edge::Type_ConstantLatitude) {
		nodeNudged.x = nodeBegin.x * (1.0 - dModNudge) + dModNudge * nodeEnd.x;
		nodeNudged.y = nodeBegin.y * (1.0 - dModNudge) + dModNudge * nodeEnd.y;
		nodeNudged.z = nodeBegin.z;

		Real dAbsNodeNudge =
			sqrt(nodeNudged.x * nodeNudged.x + nodeNudged.y * nodeNudged.y);

		Real dRadius = sqrt(1.0 - nodeNudged.z * nodeNudged.z);

		nodeNudged.x *= dRadius / dAbsNodeNudge;
		nodeNudged.y *= dRadius / dAbsNodeNudge;

	} else {
		_EXCEPTIONT("Invalid edge");
	}
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

