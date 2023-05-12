///////////////////////////////////////////////////////////////////////////////
///
///	\file	GenerateOfflineMap.cpp
///	\author  Paul Ullrich
///	\version November 26, 2017
///
///	<remarks>
///		Copyright 2000-2017 Paul Ullrich
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
#include "TempestRemapAPI.h"
#include "GridElements.h"
#include "OverlapMesh.h"
#include "DataArray3D.h"
#include "FiniteElementTools.h"
#include "SparseMatrix.h"
#include "STLStringHelper.h"
#include "NetCDFUtilities.h"
#include "triangle.h"
#include "FiniteVolumeTools.h"

#include "OverlapMesh.h"
#include "OfflineMap.h"
#include "LinearRemapSE0.h"
#include "LinearRemapFV.h"

#include "netcdfcpp.h"
#include <cmath>

///////////////////////////////////////////////////////////////////////////////

std::string g_strVersion = "GenerateOfflineMap 2.6 : 2022-10-11";

///////////////////////////////////////////////////////////////////////////////

static void ParseVariableList(
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
	const std::string & strSourceMeta,
	DataArray3D<int> & dataGLLNodes,
	DataArray3D<double> & dataGLLJacobian
) {
	NcFile ncMeta(strSourceMeta.c_str(), NcFile::ReadOnly);

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

	DataArray3D<int> dataGLLNodes_tmp;
	DataArray3D<double> dataGLLJacobian_tmp;
 
	dataGLLNodes.Allocate(nP, nP, nElem);
	dataGLLJacobian.Allocate(nP, nP, nElem);
	dataGLLNodes_tmp.Allocate(nP, nP, nElem);
	dataGLLJacobian_tmp.Allocate(nP, nP, nElem);
 
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

extern "C"
int GenerateOfflineMapWithMeshes (
	Mesh & meshSource,
	Mesh & meshTarget,
	Mesh & meshOverlap,
	std::string strSourceType,
	std::string strTargetType,
	const GenerateOfflineMapAlgorithmOptions & optsAlg,
	OfflineMap & mapRemap
) {
	NcError error(NcError::silent_nonfatal);

try {

	// Input / Output types
	enum DiscretizationType {
		DiscretizationType_FV,
		DiscretizationType_CGLL,
		DiscretizationType_DGLL
	};

	// Check metadata parameters
	if ((optsAlg.strSourceMeta != "") && (strSourceType == "fv")) {
		_EXCEPTIONT("--in_meta cannot be used with --in_type fv");
	}
	if ((optsAlg.strTargetMeta != "") && (strTargetType == "fv")) {
		_EXCEPTIONT("--out_meta cannot be used with --out_type fv");
	}

	// Check data type arguments
	std::string strNetCDFFormat = optsAlg.strOutputFormat;
	STLStringHelper::ToLower(strNetCDFFormat);

	NcFile::FileFormat eOutputFormat =
		GetNcFileFormatFromString(strNetCDFFormat);
	if (eOutputFormat == NcFile::BadFormat) {
		_EXCEPTION1("Invalid \"out_format\" value (%s), "
			"expected [Classic|Offset64Bits|Netcdf4|Netcdf4Classic]",
			optsAlg.strOutputFormat.c_str());
	}
   
	STLStringHelper::ToLower(strSourceType);
	STLStringHelper::ToLower(strTargetType);

	DiscretizationType eSourceType;
	DiscretizationType eTargetType;

	if (strSourceType == "fv") {
		eSourceType = DiscretizationType_FV;
	} else if (strSourceType == "cgll") {
		eSourceType = DiscretizationType_CGLL;
	} else if (strSourceType == "dgll") {
		eSourceType = DiscretizationType_DGLL;
	} else {
		_EXCEPTION1("Invalid \"in_type\" value (%s), expected [fv|cgll|dgll]",
			strSourceType.c_str());
	}

	if (strTargetType == "fv") {
		eTargetType = DiscretizationType_FV;
	} else if (strTargetType == "cgll") {
		eTargetType = DiscretizationType_CGLL;
	} else if (strTargetType == "dgll") {
		eTargetType = DiscretizationType_DGLL;
	} else {
		_EXCEPTION1("Invalid \"out_type\" value (%s), expected [fv|cgll|dgll]",
			strTargetType.c_str());
	}

	// Make an index of method arguments
	std::set<std::string> setMethodStrings;
	{
		int iLast = 0;
		for (int i = 0; i <= optsAlg.strMethod.length(); i++) {
			if ((i == optsAlg.strMethod.length()) || (optsAlg.strMethod[i] == ';')) {
				std::string strMethodString =
					optsAlg.strMethod.substr(iLast, i-iLast);
				STLStringHelper::RemoveWhitespaceInPlace(strMethodString);
				if (strMethodString.length() > 0) {
					setMethodStrings.insert(strMethodString);
				}
				iLast = i+1;
			}
		}
	}

	// Method flags
	std::string strMapAlgorithm("");
	int nMonotoneType = (optsAlg.fMonotone)?(1):(0);
	bool fNoConservation = optsAlg.fNoConservation;
	bool fNormalize = false;

	for (auto it : setMethodStrings) {

		// Piecewise constant monotonicity
		if (it == "mono2") {
			if (nMonotoneType != 0) {
				_EXCEPTIONT("Multiple monotonicity specifications found (--mono) or (--method \"mono#\")");
			}
			if ((eSourceType == DiscretizationType_FV) && (eTargetType == DiscretizationType_FV)) {
				_EXCEPTIONT("--method \"mono2\" is only used when remapping to/from CGLL or DGLL grids");
			}
			nMonotoneType = 2;

		// Piecewise linear monotonicity
		} else if (it == "mono3") {
			if (nMonotoneType != 0) {
				_EXCEPTIONT("Multiple monotonicity specifications found (--mono) or (--method \"mono#\")");
			}
			if ((eSourceType == DiscretizationType_FV) && (eTargetType == DiscretizationType_FV)) {
				_EXCEPTIONT("--method \"mono3\" is only used when remapping to/from CGLL or DGLL grids");
			}
			nMonotoneType = 3;

		// Volumetric remapping from FV to GLL
		} else if (it == "volumetric") {
			if ((eSourceType != DiscretizationType_FV) || (eTargetType == DiscretizationType_FV)) {
				_EXCEPTIONT("--method \"volumetric\" may only be used for FV->CGLL or FV->DGLL remapping");
			}
			strMapAlgorithm = "volumetric";

		// Inverse distance mapping
		} else if (it == "invdist") {
			if ((eSourceType != DiscretizationType_FV) || (eTargetType != DiscretizationType_FV)) {
				_EXCEPTIONT("--method \"invdist\" may only be used for FV->FV remapping");
			}
			strMapAlgorithm = "invdist";

		// Delaunay triangulation mapping
		} else if (it == "delaunay") {
			if ((eSourceType != DiscretizationType_FV) || (eTargetType != DiscretizationType_FV)) {
				_EXCEPTIONT("--method \"delaunay\" may only be used for FV->FV remapping");
			}
			strMapAlgorithm = "delaunay";

		// Bilinear
		} else if (it == "bilin") {
			if ((eSourceType != DiscretizationType_FV) || (eTargetType != DiscretizationType_FV)) {
				_EXCEPTIONT("--method \"bilin\" may only be used for FV->FV remapping");
			}
			strMapAlgorithm = "fvbilin";

		// Integrated bilinear (same as mono3 when source grid is CGLL/DGLL)
		} else if (it == "intbilin") {
			if (eTargetType != DiscretizationType_FV) {
				_EXCEPTIONT("--method \"intbilin\" may only be used when mapping to FV.");
			}
			if (eSourceType == DiscretizationType_FV) {
				strMapAlgorithm = "fvintbilin";
			} else {
				strMapAlgorithm = "mono3";
			}

		// Normalize rows
		} else if (it == "normalize") {
			fNormalize = true;

		// Integrated bilinear with generalized Barycentric coordinates
		} else if (it == "intbilingb") {
			if ((eSourceType != DiscretizationType_FV) || (eTargetType != DiscretizationType_FV)) {
				_EXCEPTIONT("--method \"intbilingb\" may only be used for FV->FV remapping");
			}
			strMapAlgorithm = "fvintbilingb";

		} else {
			_EXCEPTION1("Invalid --method argument \"%s\"", it.c_str());
		}
	}

	// Initialize dimension information from file
	if (!mapRemap.AreDimensionsInitialized()) {
		AnnounceStartBlock("Initializing dimensions of map");
		std::vector<std::string> srcDimNames, tgtDimNames;
		std::vector<int> srcDimSizes, tgtDimSizes;
		Announce("Input mesh");
		srcDimNames.push_back("num_elem");
		srcDimSizes.push_back(meshSource.faces.size());
		mapRemap.InitializeSourceDimensions(srcDimNames, srcDimSizes);
		Announce("Output mesh");
		tgtDimNames.push_back("num_elem");
		tgtDimSizes.push_back(meshTarget.faces.size());
		mapRemap.InitializeTargetDimensions(tgtDimNames, tgtDimSizes);
		AnnounceEndBlock(NULL);
	}

	// Calculate Face areas
	AnnounceStartBlock("Calculating input mesh Face areas");
	double dTotalAreaInput = meshSource.CalculateFaceAreas(optsAlg.fSourceConcave);
	Announce("Input Mesh Geometric Area: %1.15e (%1.15e)", dTotalAreaInput, dTotalAreaInput / (4.0 * M_PI));
	AnnounceEndBlock(NULL);

	// Calculate Face areas
	AnnounceStartBlock("Calculating output mesh Face areas");
	Real dTotalAreaOutput = meshTarget.CalculateFaceAreas(optsAlg.fTargetConcave);
	Announce("Output Mesh Geometric Area: %1.15e (%1.15e)", dTotalAreaOutput, dTotalAreaOutput / (4.0 * M_PI));
	AnnounceEndBlock(NULL);

	// Verify that overlap mesh is in the correct order
	int ixSourceFaceMax = (-1);
	int ixTargetFaceMax = (-1);

	if (meshOverlap.vecSourceFaceIx.size() !=
		meshOverlap.vecTargetFaceIx.size()
	) {
		_EXCEPTIONT("Invalid overlap mesh:\n"
			"	Possible mesh file corruption?");
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
	if (ixSourceFaceMax == meshSource.faces.size() //&&
		//(ixTargetFaceMax == meshTarget.faces.size())
	) {
		Announce("Overlap mesh forward correspondence found");

	// Check for reverse correspondence in overlap mesh
	} else if (
		ixSourceFaceMax == meshTarget.faces.size() //&&
		//(ixTargetFaceMax == meshSource.faces.size())
	) {
		Announce("Overlap mesh reverse correspondence found (reversing)");

		// Reorder overlap mesh
		meshOverlap.ExchangeFirstAndSecondMesh();

	// No correspondence found
	} else {
		_EXCEPTION2("Invalid overlap mesh:\n"
			"	No correspondence found with input and output meshes (%i,%i)",
			ixSourceFaceMax, ixTargetFaceMax);
	}

	AnnounceEndBlock(NULL);

	// Calculate Face areas
	AnnounceStartBlock("Calculating overlap mesh Face areas");
	Real dTotalAreaOverlap = meshOverlap.CalculateFaceAreas(false);
	Announce("Overlap Mesh Area: %1.15e (%1.15e)", dTotalAreaOverlap, dTotalAreaOverlap / (4.0 * M_PI));
	AnnounceEndBlock(NULL);

	// Correct areas to match the areas calculated in the overlap mesh
	if (!optsAlg.fNoCorrectAreas) {
		AnnounceStartBlock("Correcting source/target areas to overlap mesh areas");
		DataArray1D<double> dSourceArea(meshSource.faces.size());
		DataArray1D<double> dTargetArea(meshTarget.faces.size());

		_ASSERT(meshOverlap.vecSourceFaceIx.size() == meshOverlap.faces.size());
		_ASSERT(meshOverlap.vecTargetFaceIx.size() == meshOverlap.faces.size());
		_ASSERT(meshOverlap.vecFaceArea.GetRows() == meshOverlap.faces.size());

		_ASSERT(meshSource.vecFaceArea.GetRows() == meshSource.faces.size());
		_ASSERT(meshTarget.vecFaceArea.GetRows() == meshTarget.faces.size());

		for (int i = 0; i < meshOverlap.faces.size(); i++) {
			dSourceArea[ meshOverlap.vecSourceFaceIx[i] ] += meshOverlap.vecFaceArea[i];
			dTargetArea[ meshOverlap.vecTargetFaceIx[i] ] += meshOverlap.vecFaceArea[i];
		}

		for (int i = 0; i < meshSource.faces.size(); i++) {
			if (fabs(dSourceArea[i] - meshSource.vecFaceArea[i]) < 1.0e-10) {
				meshSource.vecFaceArea[i] = dSourceArea[i];
			}
		}
		for (int i = 0; i < meshTarget.faces.size(); i++) {
			if (fabs(dTargetArea[i] - meshTarget.vecFaceArea[i]) < 1.0e-10) {
				meshTarget.vecFaceArea[i] = dTargetArea[i];
			}
		}
		AnnounceEndBlock(NULL);
	}

	// Set source mesh areas in map
	if (eSourceType == DiscretizationType_FV) {
		mapRemap.SetSourceAreas(meshSource.vecFaceArea);
		if (meshSource.vecMask.size() != 0) {
			mapRemap.SetSourceMask(meshSource.vecMask);
		}
	}

	// Set target mesh areas in map
	if (eTargetType == DiscretizationType_FV) {
		mapRemap.SetTargetAreas(meshTarget.vecFaceArea);
		if (meshTarget.vecMask.size() != 0) {
			mapRemap.SetTargetMask(meshTarget.vecMask);
		}
	}

	// Checks
	bool fCheckConsistency = !optsAlg.fNoCheck;
	bool fCheckConservation = !optsAlg.fNoCheck;
	bool fCheckMonotonicity = (!optsAlg.fNoCheck) && (nMonotoneType != 0);

	// Partial cover on input
	if (fabs(dTotalAreaOverlap - dTotalAreaInput) > 1.0e-10) {
		if (fCheckConsistency) {
			Announce("WARNING: Significant mismatch between overlap mesh area "
				"and input mesh area.\n  Disabling checks for consistency.");
			fCheckConsistency = false;
		}
	}

/*
	// Recalculate input mesh area from overlap mesh
	if (fabs(dTotalAreaOverlap - dTotalAreaInput) > 1.0e-10) {
		AnnounceStartBlock("Overlap mesh only covers a sub-area of the sphere");
		Announce("Recalculating source mesh areas");
		dTotalAreaInput = meshSource.CalculateFaceAreasFromOverlap(meshOverlap);
		Announce("New Input Mesh Geometric Area: %1.15e", dTotalAreaInput);
		AnnounceEndBlock(NULL);
	}
*/
	// Finite volume input / Finite volume output
	if ((eSourceType  == DiscretizationType_FV) &&
		(eTargetType == DiscretizationType_FV)
	) {

		// Generate reverse node array and edge map
		meshSource.ConstructReverseNodeArray();
		meshSource.ConstructEdgeMap();

		// Initialize coordinates for map
		mapRemap.InitializeSourceCoordinatesFromMeshFV(meshSource);
		mapRemap.InitializeTargetCoordinatesFromMeshFV(meshTarget);

		// Construct OfflineMap
		if (strMapAlgorithm == "invdist") {
			AnnounceStartBlock("Calculating offline map (invdist)");
			LinearRemapFVtoFVInvDist(
				meshSource,
				meshTarget,
				meshOverlap,
				mapRemap);

		} else if (strMapAlgorithm == "delaunay") {
			AnnounceStartBlock("Calculating offline map (delaunay)");
			LinearRemapTriangulation(
				meshSource,
				meshTarget,
				meshOverlap,
				mapRemap);

		} else if (strMapAlgorithm == "fvintbilin") {
			AnnounceStartBlock("Calculating offline map (intbilin)");
			LinearRemapIntegratedBilinear(
				meshSource,
				meshTarget,
				meshOverlap,
				mapRemap);

		} else if (strMapAlgorithm == "fvintbilingb") {
			AnnounceStartBlock("Calculating offline map (intbilingb)");
			LinearRemapIntegratedGeneralizedBarycentric(
				meshSource,
				meshTarget,
				meshOverlap,
				mapRemap);

		} else if (strMapAlgorithm == "fvbilin") {
			AnnounceStartBlock("Calculating offline map (bilin)");
			LinearRemapBilinear(
				meshSource,
				meshTarget,
				meshOverlap,
				mapRemap);

		} else {
			AnnounceStartBlock("Calculating offline map (default)");
			LinearRemapFVtoFV(
				meshSource,
				meshTarget,
				meshOverlap,
				(optsAlg.fMonotone)?(1):(optsAlg.nPin),
				mapRemap);
				
				//To run the non-conservative monotone remapping schemes, uncomment the following and replace with
				//LinearRemapGeneralizedBarycentric, LinearRemapTriangulation, LinearRemapBilinear, 
				//LinearRemapIntegratedGeneralizedBarycentric, LinearRemapIntegratedTriangulation, 
				//or LinearRemapIntegratedBilinear
				
				//LinearRemapIntegratedGeneralizedBarycentric(
					//meshSource,
					//meshTarget,
					//meshOverlap,
					//mapRemap);
				
		}

	// Finite volume input / Finite element output
	} else if (eSourceType == DiscretizationType_FV) {
		DataArray3D<int> dataGLLNodes;
		DataArray3D<double> dataGLLJacobian;

		if (optsAlg.strTargetMeta != "") {
			AnnounceStartBlock("Loading meta data file");
			LoadMetaDataFile(optsAlg.strTargetMeta, dataGLLNodes, dataGLLJacobian);
			AnnounceEndBlock(NULL);

		} else {
			AnnounceStartBlock("Generating output mesh meta data");
			double dNumericalArea =
				GenerateMetaData(
					meshTarget,
					optsAlg.nPout,
					optsAlg.fNoBubble,
					dataGLLNodes,
					dataGLLJacobian,
					(eTargetType == DiscretizationType_CGLL));

			Announce("Output Mesh Numerical Area: %1.15e", dNumericalArea);
			AnnounceEndBlock(NULL);
		}

		// Initialize coordinates for map
		mapRemap.InitializeSourceCoordinatesFromMeshFV(meshSource);
		mapRemap.InitializeTargetCoordinatesFromMeshFE(
			meshTarget, optsAlg.nPout, dataGLLNodes);

		// Generate the continuous Jacobian
		bool fContinuous = (eTargetType == DiscretizationType_CGLL);

		if (eTargetType == DiscretizationType_CGLL) {
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
		meshSource.ConstructReverseNodeArray();
		meshSource.ConstructEdgeMap();

		// Generate remap weights
		if (strMapAlgorithm == "volumetric") {
			AnnounceStartBlock("Calculating offline map (volumetric)");
			LinearRemapFVtoGLL_Volumetric(
				meshSource,
				meshTarget,
				meshOverlap,
				dataGLLNodes,
				dataGLLJacobian,
				mapRemap.GetTargetAreas(),
				(optsAlg.fMonotone)?(1):(optsAlg.nPin),
				mapRemap,
				nMonotoneType,
				fContinuous,
				fNoConservation);

		} else {
			AnnounceStartBlock("Calculating offline map");
			LinearRemapFVtoGLL(
				meshSource,
				meshTarget,
				meshOverlap,
				dataGLLNodes,
				dataGLLJacobian,
				mapRemap.GetTargetAreas(),
				(optsAlg.fMonotone)?(1):(optsAlg.nPin),
				mapRemap,
				nMonotoneType,
				fContinuous,
				fNoConservation);
		}

	// Finite element input / Finite volume output
	} else if (
		(eSourceType != DiscretizationType_FV) &&
		(eTargetType == DiscretizationType_FV)
	) {
		DataArray3D<int> dataGLLNodes;
		DataArray3D<double> dataGLLJacobian;

		if (optsAlg.strSourceMeta != "") {
			AnnounceStartBlock("Loading meta data file");
			LoadMetaDataFile(optsAlg.strSourceMeta, dataGLLNodes, dataGLLJacobian);
			AnnounceEndBlock(NULL);

		} else {
			AnnounceStartBlock("Generating input mesh meta data");
			double dNumericalArea =
				GenerateMetaData(
					meshSource,
					optsAlg.nPin,
					optsAlg.fNoBubble,
					dataGLLNodes,
					dataGLLJacobian,
					(eSourceType == DiscretizationType_CGLL));

			Announce("Input Mesh Numerical Area: %1.15e", dNumericalArea);
			AnnounceEndBlock(NULL);

			if (fabs(dNumericalArea - dTotalAreaInput) > 1.0e-12) {
				Announce("WARNING: Significant mismatch between input mesh "
					"numerical area and geometric area");
			}
		}

		if (dataGLLNodes.GetSubColumns() != meshSource.faces.size()) {
			_EXCEPTIONT("Number of element does not match between metadata and "
				"input mesh");
		}

		// Initialize coordinates for map
		mapRemap.InitializeSourceCoordinatesFromMeshFE(
			meshSource, optsAlg.nPin, dataGLLNodes);
		mapRemap.InitializeTargetCoordinatesFromMeshFV(meshTarget);

		// Generate the continuous Jacobian for input mesh
		bool fContinuousIn = (eSourceType == DiscretizationType_CGLL);

		if (eSourceType == DiscretizationType_CGLL) {
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

		LinearRemapSE4(
			meshSource,
			meshTarget,
			meshOverlap,
			dataGLLNodes,
			dataGLLJacobian,
			nMonotoneType,
			fContinuousIn,
			fNoConservation,
			optsAlg.fSparseConstraints,
			mapRemap
		);

	// Finite element input / Finite element output
	} else if (
		(eSourceType  != DiscretizationType_FV) &&
		(eTargetType != DiscretizationType_FV)
	) {
		DataArray3D<int> dataGLLNodesIn;
		DataArray3D<double> dataGLLJacobianIn;

		DataArray3D<int> dataGLLNodesOut;
		DataArray3D<double> dataGLLJacobianOut;

		// Input metadata
		if (optsAlg.strSourceMeta != "") {
			AnnounceStartBlock("Loading input meta data file");
			LoadMetaDataFile(
				optsAlg.strSourceMeta, dataGLLNodesIn, dataGLLJacobianIn);
			AnnounceEndBlock(NULL);

		} else {
			AnnounceStartBlock("Generating input mesh meta data");
			double dNumericalAreaIn =
				GenerateMetaData(
					meshSource,
					optsAlg.nPin,
					optsAlg.fNoBubble,
					dataGLLNodesIn,
					dataGLLJacobianIn,
					(eSourceType == DiscretizationType_CGLL));

			Announce("Input Mesh Numerical Area: %1.15e", dNumericalAreaIn);
			AnnounceEndBlock(NULL);

			if (fabs(dNumericalAreaIn - dTotalAreaInput) > 1.0e-12) {
				Announce("WARNING: Significant mismatch between input mesh "
					"numerical area and geometric area");
			}
		}

		// Output metadata
		if (optsAlg.strTargetMeta != "") {
			AnnounceStartBlock("Loading output meta data file");
			LoadMetaDataFile(
				optsAlg.strTargetMeta, dataGLLNodesOut, dataGLLJacobianOut);
			AnnounceEndBlock(NULL);

		} else {
			AnnounceStartBlock("Generating output mesh meta data");
			double dNumericalAreaOut =
				GenerateMetaData(
					meshTarget,
					optsAlg.nPout,
					optsAlg.fNoBubble,
					dataGLLNodesOut,
					dataGLLJacobianOut,
					(eTargetType == DiscretizationType_CGLL));

			Announce("Output Mesh Numerical Area: %1.15e", dNumericalAreaOut);
			AnnounceEndBlock(NULL);

			if (fabs(dNumericalAreaOut - dTotalAreaOutput) > 1.0e-12) {
				Announce("WARNING: Significant mismatch between output mesh "
					"numerical area and geometric area");
			}
		}

		// Initialize coordinates for map
		mapRemap.InitializeSourceCoordinatesFromMeshFE(
			meshSource, optsAlg.nPin, dataGLLNodesIn);
		mapRemap.InitializeTargetCoordinatesFromMeshFE(
			meshTarget, optsAlg.nPout, dataGLLNodesOut);

		// Generate the continuous Jacobian for input mesh
		bool fContinuousIn = (eSourceType == DiscretizationType_CGLL);

		if (eSourceType == DiscretizationType_CGLL) {
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
		bool fContinuousOut = (eTargetType == DiscretizationType_CGLL);

		if (eTargetType == DiscretizationType_CGLL) {
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
			meshSource,
			meshTarget,
			meshOverlap,
			dataGLLNodesIn,
			dataGLLJacobianIn,
			dataGLLNodesOut,
			dataGLLJacobianOut,
			mapRemap.GetTargetAreas(),
			optsAlg.nPin,
			optsAlg.nPout,
			nMonotoneType,
			fContinuousIn,
			fContinuousOut,
			fNoConservation,
			mapRemap
		);

	} else {
		_EXCEPTIONT("Not implemented");
	}

	if (fNormalize) {
		AnnounceStartBlock("Normalizing weights");
		mapRemap.EnforceConsistency();
		AnnounceEndBlock("Done");
	}

	Announce("Map generation complete");
	AnnounceEndBlock(NULL);

	// Verify consistency, conservation and monotonicity
	if (!optsAlg.fNoCheck) {
		mapRemap.CheckMap(
			fCheckConsistency,
			fCheckConservation,
			fCheckMonotonicity,
			1.0e-8,			   // Normal tolerance
			1.0e-12,			  // Strict tolerance
			dTotalAreaOverlap);
	}

	// Output the Offline Map
	if (optsAlg.strOutputMapFile != "") {
		AnnounceStartBlock("Writing offline map");

		typedef std::map<std::string, std::string> AttributeMap;
		typedef AttributeMap::value_type AttributePair;

		AttributeMap mapAttributes;

		mapAttributes.insert(AttributePair("domain_a", meshSource.strFileName));
		mapAttributes.insert(AttributePair("domain_b", meshTarget.strFileName));
		mapAttributes.insert(AttributePair("grid_file_src", meshSource.strFileName));
		mapAttributes.insert(AttributePair("grid_file_dst", meshTarget.strFileName));
		mapAttributes.insert(AttributePair("grid_file_ovr", meshOverlap.strFileName));
		mapAttributes.insert(AttributePair("concave_src", (optsAlg.fSourceConcave)?("true"):("false")));
		mapAttributes.insert(AttributePair("concave_dst", (optsAlg.fTargetConcave)?("true"):("false")));
		if (optsAlg.strSourceMeta != "") {
			mapAttributes.insert(AttributePair("meta_src", optsAlg.strSourceMeta));
		}
		if (optsAlg.strTargetMeta != "") {
			mapAttributes.insert(AttributePair("meta_dst", optsAlg.strTargetMeta));
		}
		mapAttributes.insert(AttributePair("type_src", strSourceType));
		mapAttributes.insert(AttributePair("type_dst", strTargetType));
		mapAttributes.insert(AttributePair("np_src", std::to_string((long long)optsAlg.nPin)));
		mapAttributes.insert(AttributePair("np_dst", std::to_string((long long)optsAlg.nPout)));
		mapAttributes.insert(AttributePair("mono", (optsAlg.fMonotone)?("true"):("false")));
		mapAttributes.insert(AttributePair("nobubble", (optsAlg.fNoBubble)?("true"):("false")));
		mapAttributes.insert(AttributePair("nocorrectareas", (optsAlg.fNoCorrectAreas)?("true"):("false")));
		mapAttributes.insert(AttributePair("noconserve", (fNoConservation)?("true"):("false")));
		mapAttributes.insert(AttributePair("sparse_constraints", (optsAlg.fSparseConstraints)?("true"):("false")));
		mapAttributes.insert(AttributePair("method", optsAlg.strMethod));
		mapAttributes.insert(AttributePair("version", g_strVersion));

		mapRemap.Write(optsAlg.strOutputMapFile, mapAttributes, eOutputFormat);
		AnnounceEndBlock("Done");
		AnnounceBanner();
	}

	return (0);

} catch(Exception & e) {
	Announce(e.ToString().c_str());
	return (0);

} catch(...) {
	return (0);
}
}

///////////////////////////////////////////////////////////////////////////////

extern "C" 
int GenerateOfflineMap (
	std::string strSourceMesh,
	std::string strTargetMesh,
	std::string strOverlapMesh,
	std::string strSourceType,
	std::string strTargetType,
	const GenerateOfflineMapAlgorithmOptions & optsAlg,
	OfflineMap & mapRemap
) {
	NcError error(NcError::silent_nonfatal);

try {

	// Input / Output types
	enum DiscretizationType {
		DiscretizationType_FV,
		DiscretizationType_CGLL,
		DiscretizationType_DGLL
	};

	// Check command line parameters (mesh arguments)
	if (strSourceMesh == "") {
		_EXCEPTIONT("No input mesh (--in_mesh) specified");
	}
	if (strTargetMesh == "") {
		_EXCEPTIONT("No output mesh (--out_mesh) specified");
	}

	// Overlap mesh
	if (strOverlapMesh == "") {
		_EXCEPTIONT("No overlap mesh specified");
	}

	// Initialize dimension information from file
	AnnounceStartBlock("Initializing dimensions of map");
	Announce("Input mesh");
	mapRemap.InitializeSourceDimensionsFromFile(strSourceMesh);
	Announce("Output mesh");
	mapRemap.InitializeTargetDimensionsFromFile(strTargetMesh);
	AnnounceEndBlock(NULL);

	// Load input mesh
	AnnounceStartBlock("Loading input mesh");
	Mesh meshSource(strSourceMesh);
	meshSource.RemoveZeroEdges();
	AnnounceEndBlock(NULL);

	// Load output mesh
	AnnounceStartBlock("Loading output mesh");
	Mesh meshTarget(strTargetMesh);
	meshTarget.RemoveZeroEdges();
	AnnounceEndBlock(NULL);

	// Load overlap mesh
	AnnounceStartBlock("Loading overlap mesh");
	Mesh meshOverlap(strOverlapMesh);
	meshOverlap.RemoveZeroEdges();

	int err =
		GenerateOfflineMapWithMeshes(
			meshSource,
			meshTarget,
			meshOverlap,
			strSourceType,
			strTargetType,
			optsAlg,
			mapRemap);

	return err;

} catch(Exception & e) {
	Announce(e.ToString().c_str());
	return (0);

} catch(...) {
	return (0);
}
}

///////////////////////////////////////////////////////////////////////////////

extern "C"
int GenerateOfflineMapAndApply (
	std::string strSourceMesh,
	std::string strTargetMesh,
	std::string strOverlapMesh,
	std::string strSourceType,
	std::string strTargetType,
	const GenerateOfflineMapAlgorithmOptions & optsAlg,
	const ApplyOfflineMapOptions & optsApply,
	OfflineMap & mapRemap
) {
	NcError error(NcError::silent_nonfatal);

try {

	// Parse variable list
	std::vector< std::string > vecVariableStrings;
	ParseVariableList(optsApply.strVariables, vecVariableStrings);

	// Parse preserve variable list
	std::vector< std::string > vecPreserveVariableStrings;
	ParseVariableList(optsApply.strPreserveVariables, vecPreserveVariableStrings);

	if (optsApply.fPreserveAll && (vecPreserveVariableStrings.size() != 0)) {
		_EXCEPTIONT("--preserveall and --preserve cannot both be specified");
	}

	// Check data arguments
	if ((optsApply.strInputData != "") && (optsApply.strOutputData == "")) {
		_EXCEPTIONT("--in_data specified without --out_data");
	}
	if ((optsApply.strInputData == "") && (optsApply.strOutputData != "")) {
		_EXCEPTIONT("--out_data specified without --in_data");
	}

	// Generate OfflineMap
	OfflineMap mapRemap;
	int err =
		GenerateOfflineMap(
			strSourceMesh,
			strTargetMesh,
			strOverlapMesh,
			strSourceType,
			strTargetType,
			optsAlg,
			mapRemap);

	if (err != 0) return err;

	// Apply OfflineMap to data
	if (optsApply.strInputData != "") {
		AnnounceStartBlock("Applying offline map to data");

		mapRemap.SetFillValueOverrideDbl(optsApply.dFillValueOverride);
		mapRemap.SetFillValueOverride(static_cast<float>(optsApply.dFillValueOverride));
		mapRemap.Apply(
			optsApply.strInputData,
			optsApply.strOutputData,
			vecVariableStrings,
			optsApply.strNColName,
			optsApply.fOutputDouble,
			false);

		AnnounceEndBlock("Done");
	}
	AnnounceEndBlock(NULL);

	// Copy variables from input file to output file
	if ((optsApply.strInputData != "") && (optsApply.strOutputData != "")) {
		if (optsApply.fPreserveAll) {
			AnnounceStartBlock("Preserving variables");
			mapRemap.PreserveAllVariables(
				optsApply.strInputData,
				optsApply.strOutputData);
			AnnounceEndBlock("Done");

		} else if (vecPreserveVariableStrings.size() != 0) {
			AnnounceStartBlock("Preserving variables");
			mapRemap.PreserveVariables(
				optsApply.strInputData,
				optsApply.strOutputData,
				vecPreserveVariableStrings);
			AnnounceEndBlock("Done");
		}
	}

	return (0);

} catch(Exception & e) {
	Announce(e.ToString().c_str());
	return (-1);

} catch(...) {
	return (-2);
}

}

///////////////////////////////////////////////////////////////////////////////

#ifdef TEMPEST_DRIVER_MODE

int main(int argc, char** argv) {

	// Input mesh file
	std::string strSourceMesh;

	// Output mesh file
	std::string strTargetMesh;

	// Overlap mesh file
	std::string strOverlapMesh;

	// Input data type
	std::string strSourceType;

	// Output data type
	std::string strTargetType;

	// Algorithm options
	GenerateOfflineMapAlgorithmOptions optsAlg;

	// Apply options
	ApplyOfflineMapOptions optsApply;

	// NetCDF output format
	std::string strOutputFormat;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strSourceMesh, "in_mesh", "");
		CommandLineString(strTargetMesh, "out_mesh", "");
		CommandLineString(strOverlapMesh, "ov_mesh", "");
		CommandLineStringD(strSourceType, "in_type", "fv", "[fv|cgll|dgll]");
		CommandLineStringD(strTargetType, "out_type", "fv", "[fv|cgll|dgll]");

		// Optional algorithm arguments
		CommandLineString(optsAlg.strOutputMapFile, "out_map", "");
		CommandLineString(optsAlg.strSourceMeta, "in_meta", "");
		CommandLineString(optsAlg.strTargetMeta, "out_meta", "");
		CommandLineBool(optsAlg.fSourceConcave, "in_concave");
		CommandLineBool(optsAlg.fTargetConcave, "out_concave");
		CommandLineInt(optsAlg.nPin, "in_np", 4);
		CommandLineInt(optsAlg.nPout, "out_np", 4);
		CommandLineString(optsAlg.strMethod, "method", "");
		CommandLineBool(optsAlg.fMonotone, "mono");
		CommandLineBool(optsAlg.fNoBubble, "nobubble");
		CommandLineBool(optsAlg.fNoCorrectAreas, "nocorrectareas");
		CommandLineBool(optsAlg.fNoConservation, "noconserve");
		CommandLineBool(optsAlg.fNoCheck, "nocheck");
		CommandLineBool(optsAlg.fSparseConstraints, "sparse_constraints");

		// Absorbed into --method
		//CommandLineBool(fVolumetric, "volumetric");
		//CommandLineBool(fMonotoneType2, "mono2");
		//CommandLineBool(fMonotoneType3, "mono3");

		// Optional apply arguments
		CommandLineString(optsApply.strInputData, "in_data", "");
		CommandLineString(optsApply.strOutputData, "out_data", "");
		//CommandLineString(optsApply.strInputDataList, "in_data_list", "");
		//CommandLineString(optsApply.strOutputDataList, "out_data_list", "");
		CommandLineString(optsApply.strVariables, "var", "");
		CommandLineString(optsApply.strNColName, "ncol_name", "ncol");
		CommandLineBool(optsApply.fOutputDouble, "out_double");
		CommandLineString(optsApply.strPreserveVariables, "preserve", "");
		CommandLineBool(optsApply.fPreserveAll, "preserveall");
		CommandLineDouble(optsApply.dFillValueOverride, "fillvalue", 0.0);
		//CommandLineString(optsApply.strLogDir, "logdir", "");

		// Optional output format
		CommandLineStringD(strOutputFormat, "out_format","Netcdf4","[Classic|Offset64Bits|Netcdf4|Netcdf4Classic]");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Store NetCDF output format in both option lists
	optsAlg.strOutputFormat = strOutputFormat;
	optsApply.strOutputFormat = strOutputFormat;

	// Call the actual mesh generator
	OfflineMap mapRemap;
	int err =
		GenerateOfflineMapAndApply(
			strSourceMesh,
			strTargetMesh,
			strOverlapMesh,
			strSourceType,
			strTargetType,
			optsAlg,
			optsApply,
			mapRemap);

	if (err) exit(err);

	return 0;
}

#endif

///////////////////////////////////////////////////////////////////////////////
