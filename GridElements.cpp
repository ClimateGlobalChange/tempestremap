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
	const int NodesPerElement = 4;

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
	NcDim * dimNodesPerElement = ncOut.add_dim("num_nod_per_el1", 4);
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

	int * nConnect = new int[NodesPerElement];
	for (int i = 0; i < nElementCount; i++) {
		for (int k = 0; k < NodesPerElement; k++) {
			nConnect[k] = faces[i][k] + 1;
		}
		varFaces->set_cur(i, 0);
		varFaces->put(nConnect, 1, NodesPerElement);
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

	int * nEdgeType = new int[NodesPerElement];
	for (int i = 0; i < nElementCount; i++) {
		for (int k = 0; k < NodesPerElement; k++) {
			nEdgeType[k] = static_cast<int>(faces[i].GetEdge(k).type);
		}
		varEdgeTypes->set_cur(i, 0);
		varEdgeTypes->put(nEdgeType, 1, NodesPerElement);
	}
	delete[] nEdgeType;
}

///////////////////////////////////////////////////////////////////////////////

void Mesh::Read(const std::string & strFile) {

	const int NodesPerElement = 4;

	// Input from a NetCDF Exodus file
	NcFile ncFile(strFile.c_str(), NcFile::ReadOnly);

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
	faces.resize(nElementCount, Face(4));

	NcVar * varFaces = ncFile.get_var("connect1");

	int nNodes[NodesPerElement];
	for (int i = 0; i < nElementCount; i++) {
		varFaces->set_cur(i, 0);
		varFaces->get(nNodes, 1, NodesPerElement);
		faces[i].SetNode(0, nNodes[0]-1);
		faces[i].SetNode(1, nNodes[1]-1);
		faces[i].SetNode(2, nNodes[2]-1);
		faces[i].SetNode(3, nNodes[3]-1);
	}

	// Load in edge type array
	NcVar * varEdgeTypes = ncFile.get_var("edge_type");

	int nEdgeTypes[NodesPerElement];
	for (int i = 0; i < nElementCount; i++) {
		varEdgeTypes->set_cur(i, 0);
		varEdgeTypes->get(nEdgeTypes, 1, NodesPerElement);
		faces[i].GetEdge(0).type = static_cast<Edge::Type>(nEdgeTypes[0]);
		faces[i].GetEdge(1).type = static_cast<Edge::Type>(nEdgeTypes[1]);
		faces[i].GetEdge(2).type = static_cast<Edge::Type>(nEdgeTypes[2]);
		faces[i].GetEdge(3).type = static_cast<Edge::Type>(nEdgeTypes[3]);
	}
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
	std::vector<Node> & nodeIntersections
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

	// Clear the intersection vector
	nodeIntersections.clear();

	// Both edges are great circle arcs
	if ((typeFirst  == Edge::Type_GreatCircleArc) &&
		(typeSecond == Edge::Type_GreatCircleArc)
	) {
		double dDenom =
			+ node11.x * (node12.y * node21.z - node12.z * node21.y)
			- node11.y * (node12.x * node21.z - node12.z * node21.x)
			+ node11.z * (node12.x * node21.y - node12.y * node21.x);

		double dNumerA =
			+ node12.x * (node21.y * node22.z - node21.z * node22.y)
			- node12.y * (node21.x * node22.z - node21.z * node22.x)
			+ node12.z * (node21.x * node22.y - node21.y * node22.x);

		double dNumerB =
			- node11.x * (node21.y * node22.z - node21.z * node22.y)
			+ node11.y * (node21.x * node22.z - node21.z * node22.x)
			- node11.z * (node21.x * node22.y - node21.y * node22.x);

		double dNumerC =
			- node11.x * (node12.y * node22.z - node12.z * node22.y)
			+ node11.y * (node12.x * node22.z - node12.z * node22.x)
			- node11.z * (node12.x * node22.y - node12.y * node22.x);

		if (dDenom == 0.0) {
			_EXCEPTIONT("Coincident nodes");
		}

		double dMA = dNumerA / dDenom;
		double dMB = dNumerB / dDenom;
		double dMC = dNumerC / dDenom;

		double dDotNode =
			node11.x * node12.x + node11.y * node12.y + node11.z * node12.z;

		double dDenomD =
			+ dMA * dMA + 2.0 * dMA * dMB * dDotNode + dMB * dMB;

		if (fabs(dDenomD) < Tolerance) {
			return false;
		}

		double dD0 = + 1.0 / sqrt(dDenomD);
		double dD1 = - dD0;

		if ((dMA > -Tolerance) &&
			(dMB > -Tolerance) &&
			(dMC > -Tolerance) &&
			(dD0 > -Tolerance)
		) {
			nodeIntersections.resize(1);

			nodeIntersections[0].x =
				(node11.x * dMA + node12.x * dMB) * dD0;
			nodeIntersections[0].y =
				(node11.y * dMA + node12.y * dMB) * dD0;
			nodeIntersections[0].z =
				(node11.z * dMA + node12.z * dMB) * dD0;
		}
/*
		if ((dMA < 0.0) && (dMB < 0.0)) {
			nodeIntersections[1].x =
				(node11.x * dMA + node12.x * dMB) * dD1;
			nodeIntersections[1].y =
				(node11.y * dMA + node12.y * dMB) * dD1;
			nodeIntersections[1].z =
				(node11.z * dMA + node12.z * dMB) * dD1;
		}
*/
		return true;

	// First edge is a line of constant latitude; second is a great circle arc
	} else if (
		(typeFirst  == Edge::Type_GreatCircleArc) &&
		(typeSecond == Edge::Type_ConstantLatitude)
	) {

		// Cross product of basis vectors for great circle plane
		double dCrossX = + node11.y * node12.z - node11.z * node12.y;
		double dCrossY = - node11.x * node12.z + node11.z * node12.x;
		double dCrossZ = + node11.x * node12.y - node11.y * node12.x;

		// Case when node11.z is on the equator
		if (fabs(node11.z) < Tolerance) {

			// Check for equatorial great circle arc
			if (fabs(node12.z) < Tolerance) {
				if (fabs(node21.z) < Tolerance) {
					return false;
				} else {
					return true;
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

		// Case when node11.z is not on the equator
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

		return true;

	// Both edges are lines of constant latitude
	} else if (
		(typeFirst  == Edge::Type_ConstantLatitude) &&
		(typeSecond == Edge::Type_ConstantLatitude)
	) {
		if (fabs(node11.z - node21.z) < Tolerance) {
			return false;
		} else {
			return true;
		}

	// Unknown
	} else {
		_EXCEPTIONT("Invalid Edge::Type");
	}

	_EXCEPTION();
}

///////////////////////////////////////////////////////////////////////////////
