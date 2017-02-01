///////////////////////////////////////////////////////////////////////////////
///
///	\file    GenerateOverlapMesh.cpp
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

#include "Announce.h"
#include "CommandLine.h"
#include "STLStringHelper.h"
#include "Exception.h"
#include "GridElements.h"
#include "OverlapMesh.h"

#include "netcdfcpp.h"

#include <cmath>

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	NcError error(NcError::silent_nonfatal);

try {

	// Input mesh A
	std::string strMeshA;

	// Input mesh B
	std::string strMeshB;

	// Output mesh file
	std::string strOverlapMesh;

	// Overlap grid generation method
	std::string strMethod;

	// No validation of the meshes
	bool fNoValidate;

	// Concave elements may be present in mesh A
	bool fHasConcaveFacesA;

	// Concave elements may be present in mesh B
	bool fHasConcaveFacesB;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strMeshA, "a", "");
		CommandLineString(strMeshB, "b", "");
		CommandLineString(strOverlapMesh, "out", "overlap.g");
		CommandLineStringD(strMethod, "method", "fuzzy", "(fuzzy|exact|mixed)");
		CommandLineBool(fNoValidate, "novalidate");
		CommandLineBool(fHasConcaveFacesA, "concavea");
		CommandLineBool(fHasConcaveFacesB, "concaveb");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Method string
	OverlapMeshMethod method;
	STLStringHelper::ToLower(strMethod);
	if (strMethod == "fuzzy") {
		method = OverlapMeshMethod_Fuzzy;
	} else if (strMethod == "exact") {
		method = OverlapMeshMethod_Exact;
	} else if (strMethod == "mixed") {
		method = OverlapMeshMethod_Mixed;
	} else {
		_EXCEPTIONT("Invalid \"method\" value");
	}

	// Load input mesh
	AnnounceStartBlock("Loading mesh A");
	Mesh meshA(strMeshA);
	meshA.RemoveZeroEdges();
	AnnounceEndBlock(NULL);

	// Convexify Mesh
	if (fHasConcaveFacesA) {
		Mesh meshTemp = meshA;
		ConvexifyMesh(meshTemp, meshA);
	}

	// Validate mesh
	if (!fNoValidate) {
		AnnounceStartBlock("Validate mesh A");
		meshA.Validate();
		AnnounceEndBlock(NULL);
	}

	// Load output mesh
	AnnounceStartBlock("Loading mesh B");
	Mesh meshB(strMeshB);
	meshB.RemoveZeroEdges();
	AnnounceEndBlock(NULL);

	// Convexify Mesh
	if (fHasConcaveFacesB) {
		Mesh meshTemp = meshB;
		ConvexifyMesh(meshTemp, meshB);
	}

	// Validate mesh
	if (!fNoValidate) {
		AnnounceStartBlock("Validate mesh B");
		meshB.Validate();
		AnnounceEndBlock(NULL);
	}

	// Construct the edge map on both meshes
	AnnounceStartBlock("Constructing edge map on mesh A");
	meshA.ConstructEdgeMap();
	AnnounceEndBlock(NULL);

	AnnounceStartBlock("Constructing edge map on mesh B");
	meshB.ConstructEdgeMap();
	AnnounceEndBlock(NULL);

	// 
	Mesh meshOverlap;
/*
	GenerateOverlapMeshFromFace(
		meshA,
		meshB,
		0,
		meshOverlap,
		method);
*/
	AnnounceStartBlock("Construct overlap mesh\n");
	GenerateOverlapMesh_v2(meshA, meshB, meshOverlap, method);
	AnnounceEndBlock(NULL);
/*
	// Construct the reverse node array on both meshes
	AnnounceStartBlock("Constructing reverse node array on input mesh");
	meshA.ConstructReverseNodeArray();
	AnnounceEndBlock(NULL);

	AnnounceStartBlock("Constructing reverse node array on output mesh");
	meshB.ConstructReverseNodeArray();
	AnnounceEndBlock(NULL);

	// Equalize nearly coincident nodes on these Meshes
	AnnounceStartBlock("Equalize coicident Nodes");
	EqualizeCoincidentNodes(meshA, meshB);
	AnnounceEndBlock(NULL);

	// Construct the overlap mesh
	Mesh meshOverlap;

	AnnounceStartBlock("Construct overlap mesh");
	GenerateOverlapMesh(meshA, meshB, meshOverlap, method);
	AnnounceEndBlock(NULL);
*/
	// Write the overlap mesh
	AnnounceStartBlock("Writing overlap mesh");
	meshOverlap.Write(strOverlapMesh.c_str());
	AnnounceEndBlock(NULL);

	AnnounceBanner();

	return (0);

} catch(Exception & e) {
	Announce(e.ToString().c_str());
	return (-1);

} catch(...) {
	return (-2);
}
}

///////////////////////////////////////////////////////////////////////////////

