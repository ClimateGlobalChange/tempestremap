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
	static const double Tolerance = 1.0e-12;

	// Set of edges which "contain" this node
	std::set<int> setContainedEdgeIx;

	// Loop through all edges of this face
	for (int i = 0; i < edges.size(); i++) {

		// Paired edges means edge is invalid
		if (edges[i][0] == edges[i][1]) {
			continue;
		}

		// Check which side of the face this edge is on
		const Node & na = nodevec[edges[i][0]];
		const Node & nb = nodevec[edges[i][1]];

		if (edges[i].type == Edge::Type_GreatCircleArc) {
			double dDotNorm =
				  (na.y * nb.z - nb.y * na.z) * node.x
				+ (nb.x * na.z - na.x * nb.z) * node.y
				+ (na.x * nb.y - nb.x * na.y) * node.z;

			if (dDotNorm > Tolerance) {
				loc = NodeLocation_Exterior;
				ixLocation = 0;
				return;
			}
			if (dDotNorm > - Tolerance) {
				setContainedEdgeIx.insert(i);
			}

		} else if (edges[i].type == Edge::Type_ConstantLatitude) {
			double dDotNorm =
				  (na.x * nb.y - nb.x * na.y) * (node.z - na.z);

			if (dDotNorm > Tolerance) {
				loc = NodeLocation_Exterior;
				ixLocation = 0;
				return;
			}
			if (dDotNorm > - Tolerance) {
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

		loc = NodeLocation_Corner;
		ixLocation = ix1;
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
		dCoord[i] = nodes[i].x;
	}
	varNodes->set_cur(0, 0);
	varNodes->put(dCoord, 1, nNodeCount);
	for (int i = 0; i < nNodeCount; i++) {
		dCoord[i] = nodes[i].y;
	}
	varNodes->set_cur(1, 0);
	varNodes->put(dCoord, 1, nNodeCount);
	for (int i = 0; i < nNodeCount; i++) {
		dCoord[i] = nodes[i].z;
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
		nodes[i].x = dCoord[0];
		nodes[i].y = dCoord[1];
		nodes[i].z = dCoord[2];
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

///////////////////////////////////////////////////////////////////////////////
// General purpose functions
///////////////////////////////////////////////////////////////////////////////

bool CalculateEdgeIntersections(
	const Node & nodeFirstBegin,
	const Node & nodeFirstEnd,
	const Edge::Type typeFirst,
	const Node & nodeSecondBegin,
	const Node & nodeSecondEnd,
	const Edge::Type typeSecond,
	std::vector<Node> & nodeIntersections,
	bool fIncludeFirstBeginNode
) {
	static const double Tolerance = 1.0e-12;

	// Make a locally modifyable version of the Nodes
	Node node11;
	Node node12;
	Node node21;
	Node node22;

	// Second edge is a line of constant latitude; first is a great circle arc
	if ((typeFirst  == Edge::Type_ConstantLatitude) &&
		(typeSecond == Edge::Type_GreatCircleArc)
	) {
		node11 = nodeSecondBegin;
		node12 = nodeSecondEnd;
		node21 = nodeFirstBegin;
		node22 = nodeFirstEnd;

	} else {
		node11 = nodeFirstBegin;
		node12 = nodeFirstEnd;
		node21 = nodeSecondBegin;
		node22 = nodeSecondEnd;
	}

	// Check for coincident nodes
	if (node11 == node12) {
		_EXCEPTIONT("Coincident nodes used to define edge");
	}
	if (node21 == node22) {
		_EXCEPTIONT("Coincident nodes used to define edge");
	}

	// Clear the intersection vector
	nodeIntersections.clear();

	// Both edges are great circle arcs
	if ((typeFirst  == Edge::Type_GreatCircleArc) &&
		(typeSecond == Edge::Type_GreatCircleArc)
	) {
		// n11 dot n12
		double dN11oN12 =
			+ node11.x * node12.x
			+ node11.y * node12.y
			+ node11.z * node12.z;

		// Cross product of second vectors
		Node nodeN21xN22(
			+ node21.y * node22.z - node21.z * node22.y,
			- node21.x * node22.z + node21.z * node22.x,
			+ node21.x * node22.y - node21.y * node22.x);

		// Other Cross products
		Node nodeN12xN21(
			+ node12.y * node21.z - node12.z * node21.y,
			- node12.x * node21.z + node12.z * node21.x,
			+ node12.x * node21.y - node12.y * node21.x);

		Node nodeN12xN22(
			+ node12.y * node22.z - node12.z * node22.y,
			- node12.x * node22.z + node12.z * node22.x,
			+ node12.x * node22.y - node12.y * node22.x);

		// n12 dot (n21 cross n22)
		double dN12oN21xN22 =
			+ node12.x * nodeN21xN22.x
			+ node12.y * nodeN21xN22.y
			+ node12.z * nodeN21xN22.z;

		// n11 dot (n21 cross n22)
		double dN11oN21xN22 =
			+ node11.x * nodeN21xN22.x
			+ node11.y * nodeN21xN22.y
			+ node11.z * nodeN21xN22.z;

		// Check if all four vectors lay in the same plane (coincident lines)
		if ((fabs(dN11oN21xN22) < Tolerance) &&
		    (fabs(dN12oN21xN22) < Tolerance)
		) {
			// If coincident, check for intersections
			Node nodeN11xN12(
				+ node11.y * node12.z - node11.z * node12.y,
				- node11.x * node12.z + node11.z * node12.x,
				+ node11.x * node12.y - node11.y * node12.x);

			// Use the largest element of the cross product
			double dA0;
			double dA1;
			double dB0;
			double dB1;

			if ((fabs(nodeN11xN12.x) > fabs(nodeN11xN12.y)) &&
				(fabs(nodeN11xN12.x) > fabs(nodeN11xN12.z))
			) {
				dA0 = - node12.y * node21.z + node12.z * node21.y;
				dA1 = - node12.y * node22.z + node12.z * node22.y;

				dB0 = + node11.y * node21.z - node11.z * node21.y;
				dB1 = + node11.y * node22.z - node11.z * node22.y;

				dA0 /= - nodeN11xN12.x;
				dA1 /= - nodeN11xN12.x;
				dB0 /= - nodeN11xN12.x;
				dB1 /= - nodeN11xN12.x;

			} else if (
				(fabs(nodeN11xN12.y) > fabs(nodeN11xN12.x)) &&
				(fabs(nodeN11xN12.y) > fabs(nodeN11xN12.z))
			) {
				dA0 = - node12.x * node21.z + node12.z * node21.x;
				dA1 = - node12.x * node22.z + node12.z * node22.x;

				dB0 = + node11.x * node21.z - node11.z * node21.x;
				dB1 = + node11.x * node22.z - node11.z * node22.x;

				dA0 /= - nodeN11xN12.y;
				dA1 /= - nodeN11xN12.y;
				dB0 /= - nodeN11xN12.y;
				dB1 /= - nodeN11xN12.y;

			} else {
				dA0 = - node12.x * node21.y + node12.y * node21.x;
				dA1 = - node12.x * node22.y + node12.y * node22.x;

				dB0 = + node11.x * node21.y - node11.y * node21.x;
				dB1 = + node11.x * node22.y - node11.y * node22.x;

				dA0 /= nodeN11xN12.z;
				dA1 /= nodeN11xN12.z;
				dB0 /= nodeN11xN12.z;
				dB1 /= nodeN11xN12.z;
			}

			// Insert in order from FirstBegin to FirstEnd
			if (dA0 < dA1) {
				if ((dA0 > -Tolerance) && (dB0 > -Tolerance)) {
					nodeIntersections.push_back(node21);
				}
				if ((dA1 > -Tolerance) && (dB1 > -Tolerance)) {
					nodeIntersections.push_back(node22);
				}
			} else {
				if ((dA1 > -Tolerance) && (dB1 > -Tolerance)) {
					nodeIntersections.push_back(node22);
				}
				if ((dA0 > -Tolerance) && (dB0 > -Tolerance)) {
					nodeIntersections.push_back(node21);
				}
			}

			// Remove intersections that coincide with the first begin node
			if (!fIncludeFirstBeginNode) {
				for (int i = 0; i < nodeIntersections.size(); i++) {
					if (nodeIntersections[i] == nodeFirstBegin) {
						nodeIntersections.erase(nodeIntersections.begin()+i);
					}
				}
			}

			return true;
		}

		// Solution coefficients
		double dA0;
		double dB0;
		double dC0;
		double dD0;

		double dNumerC =
			+ node11.x * nodeN12xN22.x
			+ node11.y * nodeN12xN22.y
			+ node11.z * nodeN12xN22.z;

		double dNumerD =
			+ node11.x * nodeN12xN21.x
			+ node11.y * nodeN12xN21.y
			+ node11.z * nodeN12xN21.z;

		// node12 is farthest from the plane defining node21 and node22
		if (fabs(dN11oN21xN22) < fabs(dN12oN21xN22)) {

			double dNumerB =
				+ node11.x * nodeN21xN22.x
				+ node11.y * nodeN21xN22.y
				+ node11.z * nodeN21xN22.z;

			double dMB = - dNumerB / dN12oN21xN22;
			double dMC = - dNumerC / dN12oN21xN22;
			double dMD = + dNumerD / dN12oN21xN22;

			double dDenom = dMB * dMB + 2.0 * dMB * dN11oN12 + 1.0;

			if (fabs(dDenom) < Tolerance) {
				_EXCEPTIONT("Unknown case");
			}

			dA0 = 1.0 / sqrt(dDenom);
			dB0 = dA0 * dMB;
			dC0 = dA0 * dMC;
			dD0 = dA0 * dMD;

		// node11 is farthest from the plane defining node21 and node22
		} else {

			double dNumerA =
				+ node12.x * nodeN21xN22.x
				+ node12.y * nodeN21xN22.y
				+ node12.z * nodeN21xN22.z;

			double dMA = - dNumerA / dN11oN21xN22;
			double dMC = + dNumerC / dN11oN21xN22;
			double dMD = - dNumerD / dN11oN21xN22;

			double dDenom = dMA * dMA + 2.0 * dMA * dN11oN12 + 1.0;

			if (fabs(dDenom) < Tolerance) {
				_EXCEPTIONT("Unknown case");
			}

			dB0 = 1.0 / sqrt(dDenom);
			dA0 = dB0 * dMA;
			dC0 = dB0 * dMC;
			dD0 = dB0 * dMD;
		}

		// Check if first solution lies within interval
		if ((dA0 > -Tolerance) &&
			(dB0 > -Tolerance) &&
			(dC0 > -Tolerance) &&
			(dD0 > -Tolerance)
		) {
			nodeIntersections.push_back(Node(
				(node11.x * dA0 + node12.x * dB0),
				(node11.y * dA0 + node12.y * dB0),
				(node11.z * dA0 + node12.z * dB0)));

		// Check if second solution lies within interval
		} else if (
			(dA0 < Tolerance) &&
			(dB0 < Tolerance) &&
			(dC0 < Tolerance) &&
			(dD0 < Tolerance)
		) {
			nodeIntersections.push_back(Node(
				- (node11.x * dA0 + node12.x * dB0),
				- (node11.y * dA0 + node12.y * dB0),
				- (node11.z * dA0 + node12.z * dB0)));
		}

	// First edge is a line of constant latitude; second is a great circle arc
	} else if (
		(typeFirst  == Edge::Type_GreatCircleArc) &&
		(typeSecond == Edge::Type_ConstantLatitude)
	) {

		// Cross product of basis vectors for great circle plane
		double dCrossX = + node11.y * node12.z - node11.z * node12.y;
		double dCrossY = - node11.x * node12.z + node11.z * node12.x;
		double dCrossZ = + node11.x * node12.y - node11.y * node12.x;

		// node12.z is larger than node11.z
		if (fabs(node11.z) < fabs(node12.z)) {

			// Check for equatorial great circle arc
			if (fabs(node12.z) < Tolerance) {
				if (fabs(node21.z) < Tolerance) {
					_EXCEPTIONT("Not implemented");
					return true;
				} else {
					return false;
				}
			}

			// Quadratic coefficients, used to solve for A
			double dDTermA = (dCrossY * dCrossY + dCrossX * dCrossX)
				/ (node12.z * node12.z);

			double dDTermB = + 2.0 * node21.z / (node12.z * node12.z) * (
				- node12.x * dCrossY + node12.y * dCrossX);

			double dDTermC = node21.z * node21.z / (node12.z * node12.z) - 1.0;

			double dDisc = dDTermB * dDTermB - 4.0 * dDTermA * dDTermC;

			double dCross2 = node21.x * node22.y - node21.y * node22.x;

			// Only one solution
			if (fabs(dDisc) < Tolerance) {

				// Components of intersection in node1 basis
				double dA = - dDTermB / (2.0 * dDTermA);

				double dB = (-dA * node11.z + node21.z) / node12.z;

				// Components of intersection in (1,0)-(0,1) basis
				double dC = (-dA * dCrossY + node12.x * node21.z) / node12.z;

				double dD = ( dA * dCrossX + node12.y * node21.z) / node12.z;

				// Components of intersection in node2 basis
				double dE = ( dC * node22.y - dD * node22.x) / dCross2;

				double dF = (-dC * node21.y + dD * node21.x) / dCross2;

				if ((dA > -Tolerance) &&
					(dB > -Tolerance) &&
					(dE > -Tolerance) &&
					(dF > -Tolerance)
				) {
					nodeIntersections.resize(1);
					nodeIntersections[0].x = dC;
					nodeIntersections[0].y = dD;
					nodeIntersections[0].z = node21.z;
				}

			} else {
				double dSqrtDisc =
					sqrt(dDTermB * dDTermB - 4.0 * dDTermA * dDTermC);

				// Components of intersection in node1 basis
				double dA0 = (- dDTermB + dSqrtDisc) / (2.0 * dDTermA);
				double dA1 = (- dDTermB - dSqrtDisc) / (2.0 * dDTermA);

				double dB0 = (-dA0 * node11.z + node21.z) / node12.z;
				double dB1 = (-dA1 * node11.z + node21.z) / node12.z;

				// Components of intersection in (1,0,0)-(0,1,0) basis
				double dC0 = (-dA0 * dCrossY + node12.x * node21.z) / node12.z;
				double dC1 = (-dA1 * dCrossY + node12.x * node21.z) / node12.z;

				double dD0 = ( dA0 * dCrossX + node12.y * node21.z) / node12.z;
				double dD1 = ( dA1 * dCrossX + node12.y * node21.z) / node12.z;

				// Components of intersection in node2 basis
				double dE0 = ( dC0 * node22.y - dD0 * node22.x) / dCross2;
				double dE1 = ( dC1 * node22.y - dD1 * node22.x) / dCross2;

				double dF0 = (-dC0 * node21.y + dD0 * node21.x) / dCross2;
				double dF1 = (-dC1 * node21.y + dD1 * node21.x) / dCross2;

				if ((dA0 > -Tolerance) &&
					(dB0 > -Tolerance) &&
					(dE0 > -Tolerance) &&
					(dF0 > -Tolerance)
				) {
					nodeIntersections.push_back(
						Node(dC0, dD0, node21.z));
				}
				if ((dA1 > -Tolerance) &&
					(dB1 > -Tolerance) &&
					(dE1 > -Tolerance) &&
					(dF1 > -Tolerance)
				) {
					nodeIntersections.push_back(
						Node(dC1, dD1, node21.z));
				}
			}

		// node11.z is larger than node12.z
		} else {

			// Quadratic coefficients, used to solve for B
			double dDTermA = (dCrossY * dCrossY + dCrossX * dCrossX)
				/ (node11.z * node11.z);

			double dDTermB = - 2.0 * node21.z / (node11.z * node11.z) * (
				- node11.x * dCrossY + node11.y * dCrossX);

			double dDTermC = node21.z * node21.z / (node11.z * node11.z) - 1.0;

			double dDisc = dDTermB * dDTermB - 4.0 * dDTermA * dDTermC;

			double dCross2 = node21.x * node22.y - node21.y * node22.x;

			// Only one solution
			if (fabs(dDisc) < Tolerance) {

				// Components of intersection in node1 basis
				double dB = - dDTermB / (2.0 * dDTermA);

				double dA = (-dB * node12.z + node21.z) / node11.z;

				// Components of intersection in (1,0,0)-(0,1,0) basis
				double dC = (dB * dCrossY + node11.x * node21.z) / node11.z;

				double dD = (-dB * dCrossX + node11.y * node21.z) / node11.z;

				// Components of intersection in node2 basis
				double dE = ( dC * node22.y - dD * node22.x) / dCross2;

				double dF = (-dC * node21.y + dD * node21.x) / dCross2;

				if ((dA > -Tolerance) &&
					(dB > -Tolerance) &&
					(dE > -Tolerance) &&
					(dF > -Tolerance)
				) {
					nodeIntersections.resize(1);
					nodeIntersections[0].x = dC;
					nodeIntersections[0].y = dD;
					nodeIntersections[0].z = node21.z;
				}

			// Two solutions
			} else {
				double dSqrtDisc = sqrt(dDisc);

				// Components of intersection in node1 basis
				double dB0 = (- dDTermB + dSqrtDisc) / (2.0 * dDTermA);
				double dB1 = (- dDTermB - dSqrtDisc) / (2.0 * dDTermA);

				double dA0 = (-dB0 * node12.z + node21.z) / node11.z;
				double dA1 = (-dB1 * node12.z + node21.z) / node11.z;

				// Components of intersection in (1,0,0)-(0,1,0) basis
				double dC0 = (dB0 * dCrossY + node11.x * node21.z) / node11.z;
				double dC1 = (dB1 * dCrossY + node11.x * node21.z) / node11.z;

				double dD0 = (-dB0 * dCrossX + node11.y * node21.z) / node11.z;
				double dD1 = (-dB1 * dCrossX + node11.y * node21.z) / node11.z;

				// Components of intersection in node2 basis
				double dE0 = ( dC0 * node22.y - dD0 * node22.x) / dCross2;
				double dE1 = ( dC1 * node22.y - dD1 * node22.x) / dCross2;

				double dF0 = (-dC0 * node21.y + dD0 * node21.x) / dCross2;
				double dF1 = (-dC1 * node21.y + dD1 * node21.x) / dCross2;

				if ((dA0 > -Tolerance) &&
					(dB0 > -Tolerance) &&
					(dE0 > -Tolerance) &&
					(dF0 > -Tolerance)
				) {
					nodeIntersections.push_back(
						Node(dC0, dD0, node21.z));
				}
				if ((dA1 > -Tolerance) &&
					(dB1 > -Tolerance) &&
					(dE1 > -Tolerance) &&
					(dF1 > -Tolerance)
				) {
					nodeIntersections.push_back(
						Node(dC1, dD1, node21.z));
				}
			}
		}

	// Both edges are lines of constant latitude
	} else if (
		(typeFirst  == Edge::Type_ConstantLatitude) &&
		(typeSecond == Edge::Type_ConstantLatitude)
	) {
		if (fabs(node11.z - node21.z) < Tolerance) {
			_EXCEPTIONT("Not implemented");
			return true;
		} else {
			return false;
		}

	// Unknown
	} else {
		_EXCEPTIONT("Invalid Edge::Type");
	}

	// If the begin node is not to be included, erase from intersections
	if (!fIncludeFirstBeginNode) {
		for (int i = 0; i < nodeIntersections.size(); i++) {
			if (nodeIntersections[i] == nodeFirstBegin) {
				nodeIntersections.erase(nodeIntersections.begin()+i);
			}
		}
	}

	return false;
}

///////////////////////////////////////////////////////////////////////////////

void NudgeAlongEdge(
	const Node & nodeBegin,
	const Node & nodeEnd,
	const Edge::Type type,
	Node & nodeNudged
) {
	static const double Nudge = 1.0e-5;

	Node nodeDelta(
		nodeEnd.x - nodeBegin.x,
		nodeEnd.y - nodeBegin.y,
		nodeEnd.z - nodeBegin.z);

	double dModNudge = Nudge / nodeDelta.Magnitude();

	if (fabs(dModNudge) < 1.0e-12) {
		_EXCEPTIONT("Coincident Begin and End nodes");
	}

	// Edge is a great circle arc
	if (type == Edge::Type_GreatCircleArc) {

		nodeNudged.x = nodeBegin.x * (1.0 - dModNudge) + dModNudge * nodeEnd.x;
		nodeNudged.y = nodeBegin.y * (1.0 - dModNudge) + dModNudge * nodeEnd.y;
		nodeNudged.z = nodeBegin.z * (1.0 - dModNudge) + dModNudge * nodeEnd.z;

		double dAbsNodeNudge = nodeNudged.Magnitude();

		nodeNudged.x /= dAbsNodeNudge;
		nodeNudged.y /= dAbsNodeNudge;
		nodeNudged.z /= dAbsNodeNudge;

	// Edge is a line of constant latitude
	} else if (type == Edge::Type_ConstantLatitude) {
		nodeNudged.x = nodeBegin.x * (1.0 - dModNudge) + dModNudge * nodeEnd.x;
		nodeNudged.y = nodeBegin.y * (1.0 - dModNudge) + dModNudge * nodeEnd.y;
		nodeNudged.z = nodeBegin.z;

		double dAbsNodeNudge =
			sqrt(nodeNudged.x * nodeNudged.x + nodeNudged.y * nodeNudged.y);

		double dRadius = sqrt(1.0 - nodeNudged.z * nodeNudged.z);

		nodeNudged.x *= dRadius / dAbsNodeNudge;
		nodeNudged.y *= dRadius / dAbsNodeNudge;

	} else {
		_EXCEPTIONT("Invalid edge");
	}
}

///////////////////////////////////////////////////////////////////////////////

