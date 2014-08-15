///////////////////////////////////////////////////////////////////////////////
///
///	\file    gecore2.cpp
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
#include "Exception.h"
#include "GridElements.h"
#include "OverlapMesh.h"
#include "DataMatrix3D.h"
#include "FiniteElementTools.h"
#include "SparseMatrix.h"

#include "OfflineMap.h"
#include "LinearRemapSE0.h"

#include "netcdfcpp.h"
#include <cmath>

///////////////////////////////////////////////////////////////////////////////

void ParseVariableList(
	const std::string & strVariables,
	std::vector< std::string > & vecVariableStrings
) {
	int iVarBegin = 0;
	int iVarCurrent = 0;

	// Parse variable name
	for (;;) {
		if ((iVarCurrent >= strVariables.length()) ||
			(strVariables[iVarCurrent] == ',') ||
			(strVariables[iVarCurrent] == ' ')
		) {
			if (iVarCurrent == iVarBegin) {
				if (iVarCurrent >= strVariables.length()) {
					break;
				}

				continue;
			}

			vecVariableStrings.push_back(
				strVariables.substr(iVarBegin, iVarCurrent - iVarBegin));

			iVarBegin = iVarCurrent + 1;
		}

		iVarCurrent++;
	}
}

///////////////////////////////////////////////////////////////////////////////

void LoadMetaDataFile(
	const std::string & strMetaFile,
	DataMatrix3D<int> & dataGLLNodes,
	DataMatrix3D<double> & dataGLLJacobian
) {
	NcFile ncMeta(strMetaFile.c_str(), NcFile::ReadOnly);

	NcDim * dimNp = ncMeta.get_dim("np");
	NcDim * dimNelem = ncMeta.get_dim("nelem");

	NcVar * varGLLNodes = ncMeta.get_var("GLLnodes");
	NcVar * varGLLJacobian = ncMeta.get_var("J");

	int nP = dimNp->size();
	int nElem = dimNelem->size();

	dataGLLNodes.Initialize(nP, nP, nElem);
	dataGLLJacobian.Initialize(nP, nP, nElem);

	varGLLNodes->get(&(dataGLLNodes[0][0][0]), nP, nP, nElem);
	varGLLJacobian->get(&(dataGLLJacobian[0][0][0]), nP, nP, nElem);
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

try {

	// Method
	std::string strMethod;

	// Input mesh file
	std::string strInputMesh;

	// Metadata file
	std::string strMetaFile;

	// Order of polynomial in each element
	int nP;

	// Use bubble on interior of spectral element nodes
	bool fBubble;

	// Output mesh file
	std::string strOutputMesh;

	// Overlap mesh file
	std::string strOverlapMesh;

	// Variable list
	std::string strVariables;

	// Output weights file
	std::string strOutputWeights;

	// Input data file
	std::string strInputData;

	// Output data file
	std::string strOutputData;

	// Parse the command line
	BeginCommandLine()
		//CommandLineStringD(strMethod, "method", "", "[se]");
		CommandLineString(strInputMesh, "in_mesh", "");
		CommandLineString(strOutputMesh, "out_mesh", "");
		CommandLineString(strMetaFile, "in_meta", "");
		CommandLineInt(nP, "np", 4);
		CommandLineBool(fBubble, "bubble");
		CommandLineString(strOverlapMesh, "ov_mesh", "");
		CommandLineString(strVariables, "var", "");
		CommandLineString(strOutputWeights, "out_weights", "");
		CommandLineString(strInputData, "in_data", "");
		CommandLineString(strOutputData, "out_data", "");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Check command line parameters
	if ((strInputData != "") && (strOutputData == "")) {
		_EXCEPTIONT("in_data specified without out_data");
	}
	if ((strInputData == "") && (strOutputData != "")) {
		_EXCEPTIONT("out_data specified without in_data");
	}

	// Parse variable list
	std::vector< std::string > vecVariableStrings;
	ParseVariableList(strVariables, vecVariableStrings);

	if (vecVariableStrings.size() == 0) {
		_EXCEPTIONT("No variables specified");
	}

	// Load input mesh
	AnnounceStartBlock("Loading input mesh");
	Mesh meshInput(strInputMesh);
	meshInput.RemoveZeroEdges();
	AnnounceEndBlock(NULL);

	// Calculate Face areas
	AnnounceStartBlock("Calculating input mesh Face areas");
	double dTotalAreaInput = meshInput.CalculateFaceAreas();
	Announce("Input Mesh Geometric Area: %1.15e", dTotalAreaInput);
	AnnounceEndBlock(NULL);

	// Load metadata file
	DataMatrix3D<int> dataGLLNodes;
	DataMatrix3D<double> dataGLLJacobian;

	if (strMetaFile != "") {
		AnnounceStartBlock("Loading meta data file");
		LoadMetaDataFile(strMetaFile, dataGLLNodes, dataGLLJacobian);
		AnnounceEndBlock(NULL);

	} else {
		AnnounceStartBlock("Generating mesh meta data");
		double dNumericalArea =
			GenerateMetaData(
				meshInput,
				nP,
				fBubble,
				dataGLLNodes,
				dataGLLJacobian);

		Announce("Input Mesh Numerical Area: %1.15e", dNumericalArea);
		AnnounceEndBlock(NULL);

		if (fabs(dNumericalArea - dTotalAreaInput) > 1.0e-12) {
			Announce("WARNING: Significant mismatch between numerical area "
				"and geometric area\n\t(correct with --bubble)");
		}
	}

	if (dataGLLNodes.GetSubColumns() != meshInput.faces.size()) {
		_EXCEPTIONT("Number of element does not match between metadata and "
			"input mesh");
	}

	// Generate the unique Jacobian
	DataVector<double> vecInputAreas;
	GenerateUniqueJacobian(
		dataGLLNodes,
		dataGLLJacobian,
		vecInputAreas);

	// Load output mesh
	AnnounceStartBlock("Loading output mesh");
	Mesh meshOutput(strOutputMesh);
	meshOutput.RemoveZeroEdges();
	AnnounceEndBlock(NULL);

	// Calculate Face areas
	AnnounceStartBlock("Calculating output mesh Face areas");
	Real dTotalAreaOutput = meshOutput.CalculateFaceAreas();
	Announce("Output Mesh Geometric Area: %1.15e", dTotalAreaOutput);
	AnnounceEndBlock(NULL);

	// Load overlap mesh
	AnnounceStartBlock("Loading overlap mesh");
	Mesh meshOverlap(strOverlapMesh);
	meshOverlap.RemoveZeroEdges();
	AnnounceEndBlock(NULL);

	// Calculate Face areas
	AnnounceStartBlock("Calculating overlap mesh Face areas");
	Real dTotalAreaOverlap = meshOverlap.CalculateFaceAreas();
	Announce("Overlap Mesh Area: %1.15e", dTotalAreaOverlap);
	AnnounceEndBlock(NULL);

	// Generate weights file
	AnnounceStartBlock("Calculating offline map");
	OfflineMap mapRemap;
	mapRemap.InitializeOutputDimensionsFromFile(strOutputMesh);

	LinearRemapSE0(
		meshInput,
		meshOutput,
		meshOverlap,
		dataGLLNodes,
		dataGLLJacobian,
		mapRemap
	);

	// Determine first-order and conservative properties of map
	mapRemap.IsFirstOrder(1.0e-8);
	mapRemap.IsConservative(vecInputAreas, meshOutput.vecFaceArea, 1.0e-8);

	AnnounceEndBlock(NULL);

	if (strInputData != "") {
		AnnounceStartBlock("Applying offline map to data");
		mapRemap.Apply(
			vecInputAreas,
			meshOutput.vecFaceArea,
			strInputData,
			strOutputData,
			vecVariableStrings);
		AnnounceEndBlock(NULL);
	}
/*
		DataVector<double> dInput;
		dInput.Initialize(meshInput.faces.size());

		for (int i = 0; i < dInput.GetRows(); i++) {
			dInput[i] = 1.0;
		}

		DataVector<double> dOutput;
		dOutput.Initialize(meshOutput.faces.size());

		mapRemap.Apply(dInput, dOutput);

		for (int i = 0; i < 10; i++) {
			printf("%1.10e\n", dOutput[i]);
		}
*/

	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}
}

///////////////////////////////////////////////////////////////////////////////

