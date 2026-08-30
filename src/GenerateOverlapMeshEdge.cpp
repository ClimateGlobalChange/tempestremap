///////////////////////////////////////////////////////////////////////////////
///
///	\file    GenerateOverlapMeshEdge.cpp
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
#include "NetCDFUtilities.h"

#include <cmath>

///////////////////////////////////////////////////////////////////////////////

extern "C"
int GenerateOverlapMeshEdge(
	std::string strMeshA,
	std::string strMeshB,
	Mesh& meshOverlap,
	std::string strOverlapMesh,
	std::string strOutputFormat,
	std::string strMethod,
	bool fNoValidate
) {

	NcError error(NcError::silent_nonfatal);

try {
	// Check command line parameters (data type arguments)
	STLStringHelper::ToLower(strOutputFormat);

	NcFile::FileFormat eOutputFormat =
		GetNcFileFormatFromString(strOutputFormat);
	if (eOutputFormat == NcFile::BadFormat) {
		_EXCEPTION1("Invalid \"out_format\" value (%s), "
			"expected [Classic|Offset64Bits|Netcdf4|Netcdf4Classic]",
			strOutputFormat.c_str());
	}

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

	AnnounceStartBlock("Construct overlap mesh");
	GenerateOverlapMeshEdge(meshA, meshB, meshOverlap, method);
	AnnounceEndBlock(NULL);

	// Write the overlap mesh
	AnnounceStartBlock("Writing overlap mesh");
	meshOverlap.Write(strOverlapMesh.c_str(), eOutputFormat);
	AnnounceEndBlock(NULL);

} catch(Exception & e) {
	Announce(e.ToString().c_str());
	return (0);

} catch(...) {
	return (0);
}
	return (0);
}

///////////////////////////////////////////////////////////////////////////////
