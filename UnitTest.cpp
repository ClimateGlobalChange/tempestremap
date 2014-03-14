///////////////////////////////////////////////////////////////////////////////
///
///	\file    UnitTest.cpp
///	\author  Paul Ullrich
///	\version March 11, 2014
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

#include "Announce.h"
#include "CommandLine.h"
#include "Exception.h"
#include "GridElements.h"
#include "OverlapMesh.h"

#include <cmath>

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Error tolerance.
///	</summary>
static const double Tolerance = 1.0e-12;

///////////////////////////////////////////////////////////////////////////////

void TestEdgeIntersections01() {

	double dXX = 1.0/sqrt(2.0);
	double dYY = 1.0/sqrt(3.0);

	Node node11(0.0, -1.0, 0.0);
	Node node12(1.0, 0.0, 0.0);
	Node node21(dYY, -dYY, dYY);
	Node node22(dYY, -dYY, -dYY);

	std::vector<Node> nodeIntersections;

	printf("(01) Expected solution:\n");
	printf("(a) 1 intersection found\n");
	printf("Int0: 7.07107e-01 -7.07107e-01 0.00000e+00\n");
	printf("(b) 1 intersection found\n");
	printf("Int0: 7.07107e-01 -7.07107e-01 0.00000e+00\n");
	printf("(c) 0 intersection found\n");
	printf("(d) Coincident lines\n");

	printf("Computed Solution:\n");

	bool fFound;
	double dError;

	// Test 01a
	printf("(a) %i intersection found\n",
		static_cast<int>(nodeIntersections.size()));

	fFound = CalculateEdgeIntersections(
		node11, node12, Edge::Type_GreatCircleArc,
		node21, node22, Edge::Type_GreatCircleArc,
		nodeIntersections);

	if (nodeIntersections.size() != 1) {
		_EXCEPTIONT("TEST FAILED");
	}

	printf("Int0: %1.5e %1.5e %1.5e\n",
		nodeIntersections[0].x,
		nodeIntersections[0].y,
		nodeIntersections[0].z);

	dError =
		+ fabs(nodeIntersections[0].x - dXX)
		+ fabs(nodeIntersections[0].y + dXX)
		+ fabs(nodeIntersections[0].z);

	if (dError > Tolerance) {
		_EXCEPTIONT("TEST FAILED");
	}

	// Test 01b
	printf("(b) %i intersection found\n",
		static_cast<int>(nodeIntersections.size()));

	fFound = CalculateEdgeIntersections(
		node21, node22, Edge::Type_GreatCircleArc,
		node11, node12, Edge::Type_GreatCircleArc,
		nodeIntersections);

	if (nodeIntersections.size() != 1) {
		_EXCEPTIONT("TEST FAILED");
	}

	printf("Int0: %1.5e %1.5e %1.5e\n",
		nodeIntersections[0].x,
		nodeIntersections[0].y,
		nodeIntersections[0].z);

	dError =
		+ fabs(nodeIntersections[0].x - dXX)
		+ fabs(nodeIntersections[0].y + dXX)
		+ fabs(nodeIntersections[0].z);

	if (dError > Tolerance) {
		_EXCEPTIONT("TEST FAILED");
	}

	// Test 01c
	node12.x = cos(M_PI/4.0 + 0.01);
	node12.y = -sin(M_PI/4.0 + 0.01);

	fFound = CalculateEdgeIntersections(
		node11, node12, Edge::Type_GreatCircleArc,
		node21, node22, Edge::Type_GreatCircleArc,
		nodeIntersections);

	printf("(c) %i intersection found\n",
		static_cast<int>(nodeIntersections.size()));

	if (nodeIntersections.size() != 0) {
		_EXCEPTIONT("TEST FAILED");
	}

	// Test 01d
	node11 = Node(0.0, 0.0, 1.0);
	node12 = Node(1.0, 0.0, 0.0);
	node21 = Node(0.0, 0.0, -1.0);
	node22 = Node(1.0, 0.0, 0.0);

	fFound = CalculateEdgeIntersections(
		node11, node12, Edge::Type_GreatCircleArc,
		node21, node22, Edge::Type_GreatCircleArc,
		nodeIntersections);

	if (!fFound) {
		printf("(d) Coincident lines\n");
	} else {
		printf("(d) Non-coincident lines\n");
		_EXCEPTIONT("TEST FAILED\n");
	}

	printf("TEST PASSED\n");
}

///////////////////////////////////////////////////////////////////////////////

void TestEdgeIntersections02() {
	double dYY = 1.0/sqrt(3.0);

	printf("(02) Expected solution:\n");
	printf("2 intersections found\n");
	printf("Int0: -7.07107e-01 0.00000e+00 7.07107e-01\n");
	printf("Int1: 7.07107e-01 0.00000e+00 7.07107e-01\n");

	Node node11(-1.0, 0.0, 0.0);
	Node node12(0.0, 0.0, 1.0);
	Node node21(-dYY, dYY, dYY);
	Node node22(-dYY, -dYY, dYY);

	std::vector<Node> nodeIntersections;

	bool fFound = CalculateEdgeIntersections(
		node11, node12, Edge::Type_GreatCircleArc,
		node21, node22, Edge::Type_ConstantLatitude,
		nodeIntersections);

	printf("Computed Solution:\n");
	printf("%i intersections found\n",
		static_cast<int>(nodeIntersections.size()));

	if (nodeIntersections.size() != 1) {
		_EXCEPTIONT("TEST FAILED");
	}

	printf("Int: %1.5e %1.5e %1.5e\n",
		nodeIntersections[0].x,
		nodeIntersections[0].y,
		nodeIntersections[0].z);

	double dError =
		+ fabs(nodeIntersections[0].x + sqrt(2.0/3.0))
		+ fabs(nodeIntersections[0].y - 0.0)
		+ fabs(nodeIntersections[0].z - dYY);

	if (dError > Tolerance) {
		_EXCEPTIONT("TEST FAILED");
	}

	printf("TEST PASSED\n");
}

///////////////////////////////////////////////////////////////////////////////

void TestEdgeIntersections03() {
	double dXX = 1.0/sqrt(2.0);
	double dYY = 1.0/sqrt(3.0);

	Node node11(-dYY, dYY, dYY);
	Node node12(+dYY, dYY, dYY);
	Node node21(0.0, dXX, dXX);
	Node node22(dXX, 0.0, dXX);

	std::vector<Node> nodeIntersections;

	bool fFound = CalculateEdgeIntersections(
		node11, node12, Edge::Type_GreatCircleArc,
		node21, node22, Edge::Type_ConstantLatitude,
		nodeIntersections);

	printf("%i intersections found\n",
		static_cast<int>(nodeIntersections.size()));

	if (nodeIntersections.size() != 1) {
		_EXCEPTIONT("TEST FAILED");
	}

	printf("Computed Solution:\n");
	printf("Int: %1.5e %1.5e %1.5e\n",
		nodeIntersections[0].x,
		nodeIntersections[0].y,
		nodeIntersections[0].z);
}

///////////////////////////////////////////////////////////////////////////////

void TestEdgeIntersections04() {
	double dXX = 1.0/sqrt(2.0);
	double dYY = 1.0/sqrt(3.0);
	double dAA = sqrt(2.0/5.0);
	double dBB = sqrt(3.0/10.0);

	Node node11(-dYY, dYY, dYY);
	Node node12(+dYY, dYY, dYY);
	Node node21(-dBB, dBB, dAA);
	Node node22(dBB, dBB, dAA);

	std::vector<Node> nodeIntersections;

	bool fFound = CalculateEdgeIntersections(
		node11, node12, Edge::Type_GreatCircleArc,
		node21, node22, Edge::Type_ConstantLatitude,
		nodeIntersections);

	printf("%i intersections found\n",
		static_cast<int>(nodeIntersections.size()));

	if (nodeIntersections.size() != 2) {
		_EXCEPTIONT("TEST FAILED");
	}

	printf("Computed Solution:\n");
	printf("Int0: %1.5e %1.5e %1.5e\n",
		nodeIntersections[0].x,
		nodeIntersections[0].y,
		nodeIntersections[0].z);
	printf("Int1: %1.5e %1.5e %1.5e\n",
		nodeIntersections[1].x,
		nodeIntersections[1].y,
		nodeIntersections[1].z);
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

try {

	// Test edge intersections
	AnnounceBanner();
	TestEdgeIntersections01();

	AnnounceBanner();
	TestEdgeIntersections02();

	AnnounceBanner();
	TestEdgeIntersections03();

	AnnounceBanner();
	TestEdgeIntersections04();

	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}
}

///////////////////////////////////////////////////////////////////////////////


