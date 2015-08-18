///////////////////////////////////////////////////////////////////////////////
///
///	\file    GenerateOfflineMap.cpp
///	\author  Paul Ullrich
///	\version June 29, 2015
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
#include "STLStringHelper.h"

#include "OverlapMesh.h"
#include "OfflineMap.h"
#include "LinearRemapSE0.h"
#include "LinearRemapFV.h"

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
	const std::string & strInputMeta,
	DataMatrix3D<int> & dataGLLNodes,
	DataMatrix3D<double> & dataGLLJacobian
) {
	NcFile ncMeta(strInputMeta.c_str(), NcFile::ReadOnly);

	NcDim * dimNp = ncMeta.get_dim("np");
	if (dimNp == NULL) {
		_EXCEPTIONT("Dimension \"np\" missing from metadata file");
	}

	NcDim * dimNelem = ncMeta.get_dim("nelem");
	if (dimNelem == NULL) {
		_EXCEPTIONT("Dimension \"nelem\" missing from metadata file");
	}

	NcVar * varGLLNodes = ncMeta.get_var("GLLnodes");
	if (dimNelem == NULL) {
		_EXCEPTIONT("Variable \"GLLnodes\" missing from metadata file");
	}

	NcVar * varGLLJacobian = ncMeta.get_var("J");
	if (dimNelem == NULL) {
		_EXCEPTIONT("Variable \"J\" missing from metadata file");
	}

	int nP = dimNp->size();
	int nElem = dimNelem->size();

	DataMatrix3D<int> dataGLLNodes_tmp;
	DataMatrix3D<double> dataGLLJacobian_tmp;
 
	dataGLLNodes.Initialize(nP, nP, nElem);
	dataGLLJacobian.Initialize(nP, nP, nElem);
	dataGLLNodes_tmp.Initialize(nP, nP, nElem);
	dataGLLJacobian_tmp.Initialize(nP, nP, nElem);
 
	varGLLNodes->get(&(dataGLLNodes_tmp[0][0][0]), nP, nP, nElem);
 	varGLLJacobian->get(&(dataGLLJacobian_tmp[0][0][0]), nP, nP, nElem);

	for (int i = 0; i < nP; i++) {
		for (int j = 0; j < nP; j++) {
			for (int k = 0; k < nElem; k++) {
				dataGLLNodes[i][j][k] = dataGLLNodes_tmp[j][i][k];
				dataGLLJacobian[i][j][k] = dataGLLJacobian_tmp[j][i][k];
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	NcError error(NcError::silent_nonfatal);

try {

	// Input / Output types
	enum DiscretizationType {
		DiscretizationType_FV,
		DiscretizationType_CGLL,
		DiscretizationType_DGLL
	};

	// Input mesh file
	std::string strInputMesh;

	// Overlap mesh file
	std::string strOverlapMesh;

	// Input metadata file
	std::string strInputMeta;

	// Output metadata file
	std::string strOutputMeta;

	// Input data type
	std::string strInputType;

	// Output data type
	std::string strOutputType;

	// Order of polynomial in each element
	int nPin;

	// Order of polynomial in each output element
	int nPout;

	// Use bubble on interior of spectral element nodes
	bool fBubble;

	// Enforce monotonicity
	bool fMonotoneType1;

	// Enforce monotonicity
	bool fMonotoneType2;

	// Enforce monotonicity
	bool fMonotoneType3;

	// Volumetric remapping
	bool fVolumetric;

	// No conservation
	bool fNoConservation;

	// Turn off checking for conservation / consistency
	bool fNoCheck;

	// Output mesh file
	std::string strOutputMesh;

	// Variable list
	std::string strVariables;

	// Output map file
	std::string strOutputMap;

	// Input data file
	std::string strInputData;

	// Output data file
	std::string strOutputData;

	// Name of the ncol variable
	std::string strNColName;

	// Output as double
	bool fOutputDouble;

	// List of variables to preserve
	std::string strPreserveVariables;

	// Preserve all non-remapped variables
	bool fPreserveAll;

	// Fill value override
	double dFillValueOverride;

	// Parse the command line
	BeginCommandLine()
		//CommandLineStringD(strMethod, "method", "", "[se]");
		CommandLineString(strInputMesh, "in_mesh", "");
		CommandLineString(strOutputMesh, "out_mesh", "");
		CommandLineString(strOverlapMesh, "ov_mesh", "");
		CommandLineString(strInputMeta, "in_meta", "");
		CommandLineString(strOutputMeta, "out_meta", "");
		CommandLineStringD(strInputType, "in_type", "fv", "[fv|cgll|dgll]");
		CommandLineStringD(strOutputType, "out_type", "fv", "[fv|cgll|dgll]");
		CommandLineInt(nPin, "in_np", 4);
		CommandLineInt(nPout, "out_np", 4);
		CommandLineBool(fBubble, "bubble");
		CommandLineBool(fMonotoneType1, "mono");
		CommandLineBool(fMonotoneType2, "mono2");
		CommandLineBool(fMonotoneType3, "mono3");
		CommandLineBool(fVolumetric, "volumetric");
		CommandLineBool(fNoConservation, "noconserve");
		CommandLineBool(fNoCheck, "nocheck");
		CommandLineString(strVariables, "var", "");
		CommandLineString(strOutputMap, "out_map", "");
		CommandLineString(strInputData, "in_data", "");
		CommandLineString(strOutputData, "out_data", "");
		CommandLineString(strNColName, "ncol_name", "ncol");
		CommandLineBool(fOutputDouble, "out_double");
		CommandLineString(strPreserveVariables, "preserve", "");
		CommandLineBool(fPreserveAll, "preserveall");
		CommandLineDouble(dFillValueOverride, "fillvalue", 0.0);

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Check command line parameters (mesh arguments)
	if (strInputMesh == "") {
		_EXCEPTIONT("No input mesh (--in_mesh) specified");
	}
	if (strOutputMesh == "") {
		_EXCEPTIONT("No output mesh (--out_mesh) specified");
	}

	// Overlap mesh
	if (strOverlapMesh == "") {
		_EXCEPTIONT("No overlap mesh specified");
	}

	// Check command line parameters (data arguments)
	if ((strInputData != "") && (strOutputData == "")) {
		_EXCEPTIONT("--in_data specified without --out_data");
	}
	if ((strInputData == "") && (strOutputData != "")) {
		_EXCEPTIONT("--out_data specified without --in_data");
	}

	// Check metadata parameters
	if ((strInputMeta != "") && (strInputType == "fv")) {
		_EXCEPTIONT("--in_meta cannot be used with --in_type fv");
	}
	if ((strOutputMeta != "") && (strOutputType == "fv")) {
		_EXCEPTIONT("--out_meta cannot be used with --out_type fv");
	}

	// Check command line parameters (data type arguments)
	STLStringHelper::ToLower(strInputType);
	STLStringHelper::ToLower(strOutputType);

	DiscretizationType eInputType;
	DiscretizationType eOutputType;

	if (strInputType == "fv") {
		eInputType = DiscretizationType_FV;
	} else if (strInputType == "cgll") {
		eInputType = DiscretizationType_CGLL;
	} else if (strInputType == "dgll") {
		eInputType = DiscretizationType_DGLL;
	} else {
		_EXCEPTION1("Invalid \"in_type\" value (%s), expected [fv|cgll|dgll]",
			strInputType.c_str());
	}

	if (strOutputType == "fv") {
		eOutputType = DiscretizationType_FV;
	} else if (strOutputType == "cgll") {
		eOutputType = DiscretizationType_CGLL;
	} else if (strOutputType == "dgll") {
		eOutputType = DiscretizationType_DGLL;
	} else {
		_EXCEPTION1("Invalid \"out_type\" value (%s), expected [fv|cgll|dgll]",
			strOutputType.c_str());
	}

	// Monotonicity flags
	int nMonotoneType = 0;
	if (fMonotoneType1) {
		nMonotoneType = 1;
	}
	if (fMonotoneType2) {
		if (nMonotoneType != 0) {
			_EXCEPTIONT("Only one of --mono, --mono2 and --mono3 may be set");
		}
		nMonotoneType = 2;
	}
	if (fMonotoneType3) {
		if (nMonotoneType != 0) {
			_EXCEPTIONT("Only one of --mono, --mono2 and --mono3 may be set");
		}
		nMonotoneType = 3;
	}
/*
	// Volumetric
	if (fVolumetric && (nMonotoneType != 0)) {
		_EXCEPTIONT("--volumetric cannot be used in conjunction with --mono#");
	}
*/
	// Create Offline Map
	OfflineMap mapRemap;

	// Initialize dimension information from file
	AnnounceStartBlock("Initializing dimensions of map");
	Announce("Input mesh");
	mapRemap.InitializeSourceDimensionsFromFile(strInputMesh);
	Announce("Output mesh");
	mapRemap.InitializeTargetDimensionsFromFile(strOutputMesh);
	AnnounceEndBlock(NULL);

	// Parse variable list
	std::vector< std::string > vecVariableStrings;
	ParseVariableList(strVariables, vecVariableStrings);

	// Parse preserve variable list
	std::vector< std::string > vecPreserveVariableStrings;
	ParseVariableList(strPreserveVariables, vecPreserveVariableStrings);

	if (fPreserveAll && (vecPreserveVariableStrings.size() != 0)) {
		_EXCEPTIONT("--preserveall and --preserve cannot both be specified");
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

	// Input mesh areas
	if (eInputType == DiscretizationType_FV) {
		mapRemap.SetSourceAreas(meshInput.vecFaceArea);
	}

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

	// Output mesh areas
	if (eOutputType == DiscretizationType_FV) {
		mapRemap.SetTargetAreas(meshOutput.vecFaceArea);
	}

	// Load overlap mesh
	AnnounceStartBlock("Loading overlap mesh");
	Mesh meshOverlap(strOverlapMesh);
	meshOverlap.RemoveZeroEdges();

	// Verify that overlap mesh is in the correct order
	int ixSourceFaceMax = (-1);
	int ixTargetFaceMax = (-1);

	if (meshOverlap.vecSourceFaceIx.size() !=
		meshOverlap.vecTargetFaceIx.size()
	) {
		_EXCEPTIONT("Invalid overlap mesh:\n"
			"    Possible mesh file corruption?");
	}

	for (int i = 0; i < meshOverlap.vecSourceFaceIx.size(); i++) {
		if (meshOverlap.vecSourceFaceIx[i] + 1 > ixSourceFaceMax) {
			ixSourceFaceMax = meshOverlap.vecSourceFaceIx[i] + 1;
		}
		if (meshOverlap.vecTargetFaceIx[i] + 1 > ixTargetFaceMax) {
			ixTargetFaceMax = meshOverlap.vecTargetFaceIx[i] + 1;
		}
	}

	// Check for forward correspondence in overlap mesh
	if (ixSourceFaceMax == meshInput.faces.size() //&&
		//(ixTargetFaceMax == meshOutput.faces.size())
	) {
		Announce("Overlap mesh forward correspondence found");

	// Check for reverse correspondence in overlap mesh
	} else if (
		ixSourceFaceMax == meshOutput.faces.size() //&&
		//(ixTargetFaceMax == meshInput.faces.size())
	) {
		Announce("Overlap mesh reverse correspondence found (reversing)");

		// Reorder overlap mesh
		meshOverlap.ExchangeFirstAndSecondMesh();

	// No correspondence found
	} else {
		_EXCEPTION2("Invalid overlap mesh:\n"
			"    No correspondence found with input and output meshes (%i,%i)",
			ixSourceFaceMax, ixTargetFaceMax);
	}

	AnnounceEndBlock(NULL);

	// Calculate Face areas
	AnnounceStartBlock("Calculating overlap mesh Face areas");
	Real dTotalAreaOverlap = meshOverlap.CalculateFaceAreas();
	Announce("Overlap Mesh Area: %1.15e", dTotalAreaOverlap);
	AnnounceEndBlock(NULL);

	// Partial cover
	if (fabs(dTotalAreaOverlap - dTotalAreaInput) > 1.0e-10) {
		if (!fNoCheck) {
			Announce("WARNING: Significant mismatch between overlap mesh area "
				"and input mesh area.\n  Automatically enabling --nocheck");
			fNoCheck = true;
		}
	}

/*
	// Recalculate input mesh area from overlap mesh
	if (fabs(dTotalAreaOverlap - dTotalAreaInput) > 1.0e-10) {
		AnnounceStartBlock("Overlap mesh only covers a sub-area of the sphere");
		Announce("Recalculating source mesh areas");
		dTotalAreaInput = meshInput.CalculateFaceAreasFromOverlap(meshOverlap);
		Announce("New Input Mesh Geometric Area: %1.15e", dTotalAreaInput);
		AnnounceEndBlock(NULL);
	}
*/
	// Finite volume input / Finite volume output
	if ((eInputType  == DiscretizationType_FV) &&
		(eOutputType == DiscretizationType_FV)
	) {

		// Generate reverse node array and edge map
		meshInput.ConstructReverseNodeArray();
		meshInput.ConstructEdgeMap();

		// Initialize coordinates for map
		mapRemap.InitializeSourceCoordinatesFromMeshFV(meshInput);
		mapRemap.InitializeTargetCoordinatesFromMeshFV(meshOutput);

		// Construct OfflineMap
		AnnounceStartBlock("Calculating offline map");
		LinearRemapFVtoFV(
			meshInput, meshOutput, meshOverlap, nPin, mapRemap);

	// Finite volume input / Finite element output
	} else if (eInputType == DiscretizationType_FV) {
		DataMatrix3D<int> dataGLLNodes;
		DataMatrix3D<double> dataGLLJacobian;

		if (strOutputMeta != "") {
			AnnounceStartBlock("Loading meta data file");
			LoadMetaDataFile(strOutputMeta, dataGLLNodes, dataGLLJacobian);
			AnnounceEndBlock(NULL);

		} else {
			AnnounceStartBlock("Generating output mesh meta data");
			double dNumericalArea =
				GenerateMetaData(
					meshOutput,
					nPout,
					fBubble,
					dataGLLNodes,
					dataGLLJacobian);

			Announce("Output Mesh Numerical Area: %1.15e", dNumericalArea);
			AnnounceEndBlock(NULL);
		}

		// Initialize coordinates for map
		mapRemap.InitializeSourceCoordinatesFromMeshFV(meshInput);
		mapRemap.InitializeTargetCoordinatesFromMeshFE(
			meshOutput, nPout, dataGLLNodes);

		// Generate the continuous Jacobian
		bool fContinuous = (eOutputType == DiscretizationType_CGLL);

		if (eOutputType == DiscretizationType_CGLL) {
			GenerateUniqueJacobian(
				dataGLLNodes,
				dataGLLJacobian,
				mapRemap.GetTargetAreas());

		} else {
			GenerateDiscontinuousJacobian(
				dataGLLJacobian,
				mapRemap.GetTargetAreas());
		}

		// Generate reverse node array and edge map
		meshInput.ConstructReverseNodeArray();
		meshInput.ConstructEdgeMap();

		// Generate remap weights
		AnnounceStartBlock("Calculating offline map");

		if (fVolumetric) {
			LinearRemapFVtoGLL_Volumetric(
				meshInput,
				meshOutput,
				meshOverlap,
				dataGLLNodes,
				dataGLLJacobian,
				mapRemap.GetTargetAreas(),
				nPin,
				mapRemap,
				nMonotoneType,
				fContinuous,
				fNoConservation);

		} else {
			LinearRemapFVtoGLL(
				meshInput,
				meshOutput,
				meshOverlap,
				dataGLLNodes,
				dataGLLJacobian,
				mapRemap.GetTargetAreas(),
				nPin,
				mapRemap,
				nMonotoneType,
				fContinuous,
				fNoConservation);
		}

	// Finite element input / Finite volume output
	} else if (
		(eInputType != DiscretizationType_FV) &&
		(eOutputType == DiscretizationType_FV)
	) {
		DataMatrix3D<int> dataGLLNodes;
		DataMatrix3D<double> dataGLLJacobian;

		if (strInputMeta != "") {
			AnnounceStartBlock("Loading meta data file");
			LoadMetaDataFile(strInputMeta, dataGLLNodes, dataGLLJacobian);
			AnnounceEndBlock(NULL);

		} else {
			AnnounceStartBlock("Generating input mesh meta data");
			double dNumericalArea =
				GenerateMetaData(
					meshInput,
					nPin,
					fBubble,
					dataGLLNodes,
					dataGLLJacobian);

			Announce("Input Mesh Numerical Area: %1.15e", dNumericalArea);
			AnnounceEndBlock(NULL);

			if (fabs(dNumericalArea - dTotalAreaInput) > 1.0e-12) {
				Announce("WARNING: Significant mismatch between input mesh "
					"numerical area and geometric area");
			}
		}

		if (dataGLLNodes.GetSubColumns() != meshInput.faces.size()) {
			_EXCEPTIONT("Number of element does not match between metadata and "
				"input mesh");
		}

		// Initialize coordinates for map
		mapRemap.InitializeSourceCoordinatesFromMeshFE(
			meshInput, nPin, dataGLLNodes);
		mapRemap.InitializeTargetCoordinatesFromMeshFV(meshOutput);

		// Generate the continuous Jacobian for input mesh
		bool fContinuousIn = (eInputType == DiscretizationType_CGLL);

		if (eInputType == DiscretizationType_CGLL) {
			GenerateUniqueJacobian(
				dataGLLNodes,
				dataGLLJacobian,
				mapRemap.GetSourceAreas());

		} else {
			GenerateDiscontinuousJacobian(
				dataGLLJacobian,
				mapRemap.GetSourceAreas());
		}

		// Generate offline map
		AnnounceStartBlock("Calculating offline map");

		if (fVolumetric) {
			_EXCEPTIONT("Unimplemented: Volumetric currently unavailable for"
				"GLL input mesh");
		}

		LinearRemapSE4(
			meshInput,
			meshOutput,
			meshOverlap,
			dataGLLNodes,
			dataGLLJacobian,
			nMonotoneType,
			fContinuousIn,
			fNoConservation,
			mapRemap
		);

	// Finite element input / Finite element output
	} else if (
		(eInputType  != DiscretizationType_FV) &&
		(eOutputType != DiscretizationType_FV)
	) {
		DataMatrix3D<int> dataGLLNodesIn;
		DataMatrix3D<double> dataGLLJacobianIn;

		DataMatrix3D<int> dataGLLNodesOut;
		DataMatrix3D<double> dataGLLJacobianOut;

		// Input metadata
		if (strInputMeta != "") {
			AnnounceStartBlock("Loading input meta data file");
			LoadMetaDataFile(
				strInputMeta, dataGLLNodesIn, dataGLLJacobianIn);
			AnnounceEndBlock(NULL);

		} else {
			AnnounceStartBlock("Generating input mesh meta data");
			double dNumericalAreaIn =
				GenerateMetaData(
					meshInput,
					nPin,
					fBubble,
					dataGLLNodesIn,
					dataGLLJacobianIn);

			Announce("Input Mesh Numerical Area: %1.15e", dNumericalAreaIn);
			AnnounceEndBlock(NULL);

			if (fabs(dNumericalAreaIn - dTotalAreaInput) > 1.0e-12) {
				Announce("WARNING: Significant mismatch between input mesh "
					"numerical area and geometric area");
			}
		}

		// Output metadata
		if (strOutputMeta != "") {
			AnnounceStartBlock("Loading output meta data file");
			LoadMetaDataFile(
				strOutputMeta, dataGLLNodesOut, dataGLLJacobianOut);
			AnnounceEndBlock(NULL);

		} else {
			AnnounceStartBlock("Generating output mesh meta data");
			double dNumericalAreaOut =
				GenerateMetaData(
					meshOutput,
					nPout,
					fBubble,
					dataGLLNodesOut,
					dataGLLJacobianOut);

			Announce("Output Mesh Numerical Area: %1.15e", dNumericalAreaOut);
			AnnounceEndBlock(NULL);

			if (fabs(dNumericalAreaOut - dTotalAreaOutput) > 1.0e-12) {
				Announce("WARNING: Significant mismatch between output mesh "
					"numerical area and geometric area");
			}
		}

		// Initialize coordinates for map
		mapRemap.InitializeSourceCoordinatesFromMeshFE(
			meshInput, nPin, dataGLLNodesIn);
		mapRemap.InitializeTargetCoordinatesFromMeshFE(
			meshOutput, nPout, dataGLLNodesOut);

		// Generate the continuous Jacobian for input mesh
		bool fContinuousIn = (eInputType == DiscretizationType_CGLL);

		if (eInputType == DiscretizationType_CGLL) {
			GenerateUniqueJacobian(
				dataGLLNodesIn,
				dataGLLJacobianIn,
				mapRemap.GetSourceAreas());

		} else {
			GenerateDiscontinuousJacobian(
				dataGLLJacobianIn,
				mapRemap.GetSourceAreas());
		}

		// Generate the continuous Jacobian for output mesh
		bool fContinuousOut = (eOutputType == DiscretizationType_CGLL);

		if (eOutputType == DiscretizationType_CGLL) {
			GenerateUniqueJacobian(
				dataGLLNodesOut,
				dataGLLJacobianOut,
				mapRemap.GetTargetAreas());

		} else {
			GenerateDiscontinuousJacobian(
				dataGLLJacobianOut,
				mapRemap.GetTargetAreas());
		}

		// Generate offline map
		AnnounceStartBlock("Calculating offline map");

		LinearRemapGLLtoGLL2(
			meshInput,
			meshOutput,
			meshOverlap,
			dataGLLNodesIn,
			dataGLLJacobianIn,
			dataGLLNodesOut,
			dataGLLJacobianOut,
			mapRemap.GetTargetAreas(),
			nPin,
			nPout,
			nMonotoneType,
			fContinuousIn,
			fContinuousOut,
			fNoConservation,
			mapRemap
		);

	} else {
		_EXCEPTIONT("Not implemented");
	}

//#pragma warning "NOTE: VERIFICATION DISABLED"

	// Verify consistency, conservation and monotonicity
	if (!fNoCheck) {
		AnnounceStartBlock("Verifying map");
		mapRemap.IsConsistent(1.0e-8);
		mapRemap.IsConservative(1.0e-8);

		if (nMonotoneType != 0) {
			mapRemap.IsMonotone(1.0e-12);
		}
		AnnounceEndBlock(NULL);
	}

	AnnounceEndBlock(NULL);

	// Initialize element dimensions from input/output Mesh
	AnnounceStartBlock("Writing output");

	// Output the Offline Map
	if (strOutputMap != "") {
		AnnounceStartBlock("Writing offline map");
		mapRemap.Write(strOutputMap);
		AnnounceEndBlock(NULL);
	}

	// Apply Offline Map to data
	if (strInputData != "") {
		AnnounceStartBlock("Applying offline map to data");

		mapRemap.SetFillValueOverride(static_cast<float>(dFillValueOverride));
		mapRemap.Apply(
			strInputData,
			strOutputData,
			vecVariableStrings,
			strNColName,
			fOutputDouble,
			false);
		AnnounceEndBlock(NULL);
	}
	AnnounceEndBlock(NULL);

	// Copy variables from input file to output file
	if ((strInputData != "") && (strOutputData != "")) {
		if (fPreserveAll) {
			AnnounceStartBlock("Preserving variables");
			mapRemap.PreserveAllVariables(strInputData, strOutputData);
			AnnounceEndBlock(NULL);

		} else if (vecPreserveVariableStrings.size() != 0) {
			AnnounceStartBlock("Preserving variables");
			mapRemap.PreserveVariables(
				strInputData,
				strOutputData,
				vecPreserveVariableStrings);
			AnnounceEndBlock(NULL);
		}
	}

	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}
}

///////////////////////////////////////////////////////////////////////////////

