///////////////////////////////////////////////////////////////////////////////
///
///	\file    GenerateOfflineMap.cpp
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
#include "GridElements.h"
#include "OverlapMesh.h"
#include "DataArray3D.h"
#include "FiniteElementTools.h"
#include "SparseMatrix.h"
#include "STLStringHelper.h"
#include "NetCDFUtilities.h"

#include "OverlapMesh.h"
#include "OfflineMap.h"
#include "LinearRemapSE0.h"
#include "LinearRemapFV.h"

#include "netcdfcpp.h"
#include <cmath>

///////////////////////////////////////////////////////////////////////////////

std::string g_strVersion = "GenerateOfflineMap 2.0 : 2017-11-26";

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
int GenerateOfflineMapWithMeshes(
	OfflineMap& mapRemap,
	Mesh& meshSource,
	Mesh& meshTarget,
	Mesh& meshOverlap,
	std::string strSourceMeta,
	std::string strTargetMeta,
	std::string strSourceType, std::string strTargetType,
	int nPin, int nPout,
	bool fNoBubble, bool fCorrectAreas, int fMonotoneTypeID,
	bool fVolumetric, bool fNoConservation, bool fNoCheck,
	std::string strVariables, std::string strOutputMap,
	std::string strInputData, std::string strOutputData,
	std::string strNColName, bool fOutputDouble,
	std::string strOutputFormat,
	std::string strPreserveVariables, bool fPreserveAll, double dFillValueOverride,
	bool fSourceConcave, bool fTargetConcave, double lb, double ub, bool fCAAS, bool fCAASLocal
) {
	NcError error(NcError::silent_nonfatal);

try {

    // Input / Output types
    enum DiscretizationType {
        DiscretizationType_FV,
        DiscretizationType_CGLL,
        DiscretizationType_DGLL
    };

    // Check command line parameters (data arguments)
    if ((strInputData != "") && (strOutputData == "")) {
        _EXCEPTIONT("--in_data specified without --out_data");
    }
    if ((strInputData == "") && (strOutputData != "")) {
        _EXCEPTIONT("--out_data specified without --in_data");
    }

    // Check metadata parameters
    if ((strSourceMeta != "") && (strSourceType == "fv")) {
        _EXCEPTIONT("--in_meta cannot be used with --in_type fv");
    }
    if ((strTargetMeta != "") && (strTargetType == "fv")) {
        _EXCEPTIONT("--out_meta cannot be used with --out_type fv");
    }

    // Check command line parameters (data type arguments)
    STLStringHelper::ToLower(strOutputFormat);

	NcFile::FileFormat eOutputFormat =
		GetNcFileFormatFromString(strOutputFormat);
	if (eOutputFormat == NcFile::BadFormat) {
		_EXCEPTION1("Invalid \"out_format\" value (%s), "
			"expected [Classic|Offset64Bits|Netcdf4|Netcdf4Classic]",
			strOutputFormat.c_str());
	}
   
    STLStringHelper::ToLower(strSourceType);
    STLStringHelper::ToLower(strTargetType);

    DiscretizationType eInputType;
    DiscretizationType eOutputType;

    if (strSourceType == "fv") {
        eInputType = DiscretizationType_FV;
    } else if (strSourceType == "cgll") {
        eInputType = DiscretizationType_CGLL;
    } else if (strSourceType == "dgll") {
        eInputType = DiscretizationType_DGLL;
    } else {
        _EXCEPTION1("Invalid \"in_type\" value (%s), expected [fv|cgll|dgll]",
            strSourceType.c_str());
    }

    if (strTargetType == "fv") {
        eOutputType = DiscretizationType_FV;
    } else if (strTargetType == "cgll") {
        eOutputType = DiscretizationType_CGLL;
    } else if (strTargetType == "dgll") {
        eOutputType = DiscretizationType_DGLL;
    } else {
        _EXCEPTION1("Invalid \"out_type\" value (%s), expected [fv|cgll|dgll]",
            strTargetType.c_str());
    }

    // Monotonicity flags
    int nMonotoneType = fMonotoneTypeID;

/*
    // Volumetric
    if (fVolumetric && (nMonotoneType != 0)) {
        _EXCEPTIONT("--volumetric cannot be used in conjunction with --mono#");
    }
*/

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

    // Parse variable list
    std::vector< std::string > vecVariableStrings;
    ParseVariableList(strVariables, vecVariableStrings);

    // Parse preserve variable list
    std::vector< std::string > vecPreserveVariableStrings;
    ParseVariableList(strPreserveVariables, vecPreserveVariableStrings);

    if (fPreserveAll && (vecPreserveVariableStrings.size() != 0)) {
        _EXCEPTIONT("--preserveall and --preserve cannot both be specified");
    }

    // Calculate Face areas
    AnnounceStartBlock("Calculating input mesh Face areas");
    double dTotalAreaInput = meshSource.CalculateFaceAreas(fSourceConcave);
    Announce("Input Mesh Geometric Area: %1.15e (%1.15e)", dTotalAreaInput, dTotalAreaInput / (4.0 * M_PI));
    AnnounceEndBlock(NULL);

    // Calculate Face areas
    AnnounceStartBlock("Calculating output mesh Face areas");
    Real dTotalAreaOutput = meshTarget.CalculateFaceAreas(fTargetConcave);
    Announce("Output Mesh Geometric Area: %1.15e (%1.15e)", dTotalAreaOutput, dTotalAreaOutput / (4.0 * M_PI));
    AnnounceEndBlock(NULL);

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
            "    No correspondence found with input and output meshes (%i,%i)",
            ixSourceFaceMax, ixTargetFaceMax);
    }

    AnnounceEndBlock(NULL);

    // Calculate Face areas
    AnnounceStartBlock("Calculating overlap mesh Face areas");
    Real dTotalAreaOverlap = meshOverlap.CalculateFaceAreas(false);
    Announce("Overlap Mesh Area: %1.15e (%1.15e)", dTotalAreaOverlap, dTotalAreaOverlap / (4.0 * M_PI));
    AnnounceEndBlock(NULL);

	// Correct areas to match the areas calculated in the overlap mesh
	if (fCorrectAreas) {
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
    if (eInputType == DiscretizationType_FV) {
        mapRemap.SetSourceAreas(meshSource.vecFaceArea);
		if (meshSource.vecMask.IsAttached()) {
			mapRemap.SetSourceMask(meshSource.vecMask);
		}
    }

    // Set target mesh areas in map
    if (eOutputType == DiscretizationType_FV) {
        mapRemap.SetTargetAreas(meshTarget.vecFaceArea);
		if (meshTarget.vecMask.IsAttached()) {
			mapRemap.SetTargetMask(meshTarget.vecMask);
		}
    }

	// Checks
	bool fCheckConsistency = !fNoCheck;
	bool fCheckConservation = !fNoCheck;
	bool fCheckMonotonicity = (!fNoCheck) && (nMonotoneType != 0);

    // Partial cover on input
    if (fabs(dTotalAreaOverlap - dTotalAreaInput) > 1.0e-10) {
        if (fCheckConsistency) {
            Announce("WARNING: Significant mismatch between overlap mesh area "
                "and input mesh area.\n  Disabling checks for consistency.");
            fCheckConsistency = false;
        }
    }

        DataArray3D<int> dataGLLNodes;
		DataArray3D<int> dataGLLNodesIn;
		DataArray3D<int> dataGLLNodesOut;

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
    if ((eInputType  == DiscretizationType_FV) &&
        (eOutputType == DiscretizationType_FV)
    ) {

        // Generate reverse node array and edge map
        meshSource.ConstructReverseNodeArray();
        meshSource.ConstructEdgeMap();

        // Initialize coordinates for map
        mapRemap.InitializeSourceCoordinatesFromMeshFV(meshSource);
        mapRemap.InitializeTargetCoordinatesFromMeshFV(meshTarget);

        // Construct OfflineMap
        AnnounceStartBlock("Calculating offline map");
        LinearRemapFVtoFV(meshSource, meshTarget, meshOverlap, nPin, mapRemap);

    // Finite volume input / Finite element output
    } else if (eInputType == DiscretizationType_FV) {
        //DataArray3D<int> dataGLLNodes;
        DataArray3D<double> dataGLLJacobian;

        if (strTargetMeta != "") {
            AnnounceStartBlock("Loading meta data file");
            LoadMetaDataFile(strTargetMeta, dataGLLNodes, dataGLLJacobian);
            AnnounceEndBlock(NULL);

        } else {
            AnnounceStartBlock("Generating output mesh meta data");
            double dNumericalArea =
                GenerateMetaData(
                    meshTarget,
                    nPout,
                    fNoBubble,
                    dataGLLNodes,
                    dataGLLJacobian);

            Announce("Output Mesh Numerical Area: %1.15e", dNumericalArea);
            AnnounceEndBlock(NULL);
            
        }

        // Initialize coordinates for map
        mapRemap.InitializeSourceCoordinatesFromMeshFV(meshSource);
        mapRemap.InitializeTargetCoordinatesFromMeshFE(
            meshTarget, nPout, dataGLLNodes);

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
        meshSource.ConstructReverseNodeArray();
        meshSource.ConstructEdgeMap();

        // Generate remap weights
        AnnounceStartBlock("Calculating offline map");

        if (fVolumetric) {
            LinearRemapFVtoGLL_Volumetric(
                meshSource,
                meshTarget,
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
                meshSource,
                meshTarget,
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
        //DataArray3D<int> dataGLLNodes;
        DataArray3D<double> dataGLLJacobian;

        if (strSourceMeta != "") {
            AnnounceStartBlock("Loading meta data file");
            LoadMetaDataFile(strSourceMeta, dataGLLNodes, dataGLLJacobian);
            AnnounceEndBlock(NULL);

        } else {
            AnnounceStartBlock("Generating input mesh meta data");
            double dNumericalArea =
                GenerateMetaData(
                    meshSource,
                    nPin,
                    fNoBubble,
                    dataGLLNodes,
                    dataGLLJacobian);

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
            meshSource, nPin, dataGLLNodes);
        mapRemap.InitializeTargetCoordinatesFromMeshFV(meshTarget);

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
            meshSource,
            meshTarget,
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
        //DataArray3D<int> dataGLLNodesIn;
        DataArray3D<double> dataGLLJacobianIn;

        //DataArray3D<int> dataGLLNodesOut;
        DataArray3D<double> dataGLLJacobianOut;

        // Input metadata
        if (strSourceMeta != "") {
            AnnounceStartBlock("Loading input meta data file");
            LoadMetaDataFile(
                strSourceMeta, dataGLLNodesIn, dataGLLJacobianIn);
            AnnounceEndBlock(NULL);

        } else {
            AnnounceStartBlock("Generating input mesh meta data");
            double dNumericalAreaIn =
                GenerateMetaData(
                    meshSource,
                    nPin,
                    fNoBubble,
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
        if (strTargetMeta != "") {
            AnnounceStartBlock("Loading output meta data file");
            LoadMetaDataFile(
                strTargetMeta, dataGLLNodesOut, dataGLLJacobianOut);
            AnnounceEndBlock(NULL);

        } else {
            AnnounceStartBlock("Generating output mesh meta data");
            double dNumericalAreaOut =
                GenerateMetaData(
                    meshTarget,
                    nPout,
                    fNoBubble,
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
            meshSource, nPin, dataGLLNodesIn);
        mapRemap.InitializeTargetCoordinatesFromMeshFE(
            meshTarget, nPout, dataGLLNodesOut);

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
            meshSource,
            meshTarget,
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

	Announce("Map generation complete");
    AnnounceEndBlock(NULL);

    // Verify consistency, conservation and monotonicity
	if (!fNoCheck) {
		mapRemap.CheckMap(
			fCheckConsistency,
			fCheckConservation,
			fCheckMonotonicity,
			1.0e-8,               // Normal tolerance
			1.0e-12,              // Strict tolerance
			dTotalAreaOverlap);
	}

    // Initialize element dimensions from input/output Mesh
    AnnounceStartBlock("Writing output");

    // Output the Offline Map
    if (strOutputMap != "") {
        AnnounceStartBlock("Writing offline map");

		typedef std::map<std::string, std::string> AttributeMap;
		typedef AttributeMap::value_type AttributePair;

		AttributeMap mapAttributes;

		mapAttributes.insert(AttributePair("domain_a", meshSource.strFileName));
		mapAttributes.insert(AttributePair("domain_b", meshTarget.strFileName));
		mapAttributes.insert(AttributePair("grid_file_src", meshSource.strFileName));
		mapAttributes.insert(AttributePair("grid_file_dst", meshTarget.strFileName));
		mapAttributes.insert(AttributePair("grid_file_ovr", meshOverlap.strFileName));
		if (strSourceMeta != "") {
			mapAttributes.insert(AttributePair("meta_src", strSourceMeta));
		}
		if (strTargetMeta != "") {
			mapAttributes.insert(AttributePair("meta_dst", strTargetMeta));
		}
		mapAttributes.insert(AttributePair("type_src", strSourceType));
		mapAttributes.insert(AttributePair("type_dst", strTargetType));
		mapAttributes.insert(AttributePair("np_src", std::to_string((long long)nPin)));
		mapAttributes.insert(AttributePair("np_dst", std::to_string((long long)nPout)));
		mapAttributes.insert(AttributePair("bubble", (!fNoBubble)?("true"):("false")));
		mapAttributes.insert(AttributePair("mono_type", std::to_string((long long)fMonotoneTypeID)));
		if (fVolumetric) {
			mapAttributes.insert(AttributePair("volumetric", "true"));
		}
		if (fNoConservation) {
			mapAttributes.insert(AttributePair("no_conserve", "true"));
		}
		mapAttributes.insert(AttributePair("concave_src", (fSourceConcave)?("true"):("false")));
		mapAttributes.insert(AttributePair("concave_dst", (fTargetConcave)?("true"):("false")));
		mapAttributes.insert(AttributePair("version", g_strVersion));

        mapRemap.Write(strOutputMap, mapAttributes, eOutputFormat);
        AnnounceEndBlock("Done");
    }

    // Apply Offline Map to data
    if (strInputData != "") {
        AnnounceStartBlock("Applying offline map to data");

        mapRemap.SetFillValueOverride(static_cast<float>(dFillValueOverride));
        
        //DAVE ADDS NEXT TWO CONDITIONALS
        if((eInputType  == DiscretizationType_FV) &&
        (eOutputType != DiscretizationType_FV)){
			dataGLLNodesOut = dataGLLNodes;
		}
		
		else if((eInputType  != DiscretizationType_FV) &&
        (eOutputType == DiscretizationType_FV)){
			dataGLLNodesIn = dataGLLNodes;
		}
        
        mapRemap.Apply(
            strInputData,
            strOutputData,
            vecVariableStrings,
            strNColName,
            meshOverlap,
            meshSource,
            nPin,
            dataGLLNodesIn,
            dataGLLNodesOut,
            lb,
            ub,
            fOutputDouble,
            false,
            fCAAS,
            fCAASLocal);
        AnnounceEndBlock("Done");
    }
    AnnounceEndBlock(NULL);

    // Copy variables from input file to output file
    if ((strInputData != "") && (strOutputData != "")) {
        if (fPreserveAll) {
            AnnounceStartBlock("Preserving variables");
            mapRemap.PreserveAllVariables(strInputData, strOutputData);
            AnnounceEndBlock("Done");

        } else if (vecPreserveVariableStrings.size() != 0) {
            AnnounceStartBlock("Preserving variables");
            mapRemap.PreserveVariables(
                strInputData,
                strOutputData,
                vecPreserveVariableStrings);
            AnnounceEndBlock("Done");
        }
    }

	AnnounceBanner();

    return (0);

} catch(Exception & e) {
    Announce(e.ToString().c_str());
    return (0);

} catch(...) {
    return (0);
}
}


extern "C" 
int GenerateOfflineMap(
	OfflineMap& mapRemap,
	std::string strInputMesh, std::string strOutputMesh,
	std::string strOverlapMesh,
	std::string strSourceMeta, std::string strTargetMeta,
	std::string strSourceType, std::string strTargetType,
	int nPin, int nPout,
	bool fNoBubble,
	bool fCorrectAreas,
	int fMonotoneTypeID,
	bool fVolumetric,
	bool fNoConservation,
	bool fNoCheck,
	std::string strVariables,
	std::string strOutputMap,
	std::string strInputData, std::string strOutputData,
	std::string strNColName,
	bool fOutputDouble,
	std::string strOutputFormat,
	std::string strPreserveVariables,
	bool fPreserveAll,
	double dFillValueOverride,
	bool fSourceConcave, bool fTargetConcave, double lb, double ub, bool fCAAS, bool fCAASLocal
	 )
{
	NcError error(NcError::silent_nonfatal);

try {

	// Input / Output types
	enum DiscretizationType {
		DiscretizationType_FV,
		DiscretizationType_CGLL,
		DiscretizationType_DGLL
	};

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

	// Initialize dimension information from file
	AnnounceStartBlock("Initializing dimensions of map");
	Announce("Input mesh");
	mapRemap.InitializeSourceDimensionsFromFile(strInputMesh);
	Announce("Output mesh");
	mapRemap.InitializeTargetDimensionsFromFile(strOutputMesh);
	AnnounceEndBlock(NULL);

	// Load input mesh
	AnnounceStartBlock("Loading input mesh");
	Mesh meshSource(strInputMesh);
	meshSource.RemoveZeroEdges();
	AnnounceEndBlock(NULL);

	// Load output mesh
	AnnounceStartBlock("Loading output mesh");
	Mesh meshTarget(strOutputMesh);
	meshTarget.RemoveZeroEdges();
	AnnounceEndBlock(NULL);

	// Load overlap mesh
	AnnounceStartBlock("Loading overlap mesh");
	Mesh meshOverlap(strOverlapMesh);
	meshOverlap.RemoveZeroEdges();

    int err = GenerateOfflineMapWithMeshes(mapRemap, meshSource, meshTarget, meshOverlap,
                                            strSourceMeta, strTargetMeta,
                                            strSourceType, strTargetType,
                                            nPin, nPout,
                                            fNoBubble, fCorrectAreas, fMonotoneTypeID,
                                            fVolumetric, fNoConservation, fNoCheck,
                                            strVariables, strOutputMap,
                                            strInputData, strOutputData,
                                            strNColName, fOutputDouble, strOutputFormat,
                                            strPreserveVariables, fPreserveAll, dFillValueOverride,
                                            fSourceConcave, fTargetConcave, lb, ub, fCAAS, fCAASLocal );

    return err;

} catch(Exception & e) {
	Announce(e.ToString().c_str());
	return (0);

} catch(...) {
	return (0);
}
}

///////////////////////////////////////////////////////////////////////////////

#ifdef TEMPEST_DRIVER_MODE

int main(int argc, char** argv) {

	// Input mesh file
	std::string strInputMesh;

	// Overlap mesh file
	std::string strOverlapMesh;

	// Input metadata file
	std::string strSourceMeta;

	// Output metadata file
	std::string strTargetMeta;

	// Input data type
	std::string strSourceType;

	// Output data type
	std::string strTargetType;

	// Output mesh file
	std::string strOutputMesh;

	// Order of polynomial in each element
	int nPin;

	// Order of polynomial in each output element
	int nPout;

	// Use bubble on interior of spectral element nodes
	bool fNoBubble;

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

	// Output data Netcdf format
	std::string strOutputFormat;

	// List of variables to preserve
	std::string strPreserveVariables;

	// Preserve all non-remapped variables
	bool fPreserveAll;

	// Fill value override
	double dFillValueOverride;

	// Input mesh contains concave elements
	bool fSourceConcave;

	// Output mesh contains concave elements
	bool fTargetConcave;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputMesh, "in_mesh", "");
		CommandLineString(strOutputMesh, "out_mesh", "");
		CommandLineString(strOverlapMesh, "ov_mesh", "");
		CommandLineString(strSourceMeta, "in_meta", "");
		CommandLineString(strTargetMeta, "out_meta", "");
		CommandLineStringD(strSourceType, "in_type", "fv", "[fv|cgll|dgll]");
		CommandLineStringD(strTargetType, "out_type", "fv", "[fv|cgll|dgll]");

		// Optional arguments
		CommandLineInt(nPin, "in_np", 4);
		CommandLineInt(nPout, "out_np", 4);
		CommandLineBool(fNoBubble, "no_bubble");
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
                CommandLineString(strOutputFormat, "out_format","Classic","[Classic|Offset64Bits|Netcdf4|Netcdf4Classic]");
		CommandLineString(strPreserveVariables, "preserve", "");
		CommandLineBool(fPreserveAll, "preserveall");
		CommandLineDouble(dFillValueOverride, "fillvalue", 0.0);
		CommandLineBool(fSourceConcave, "in_concave");
		CommandLineBool(fTargetConcave, "out_concave");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	int fMonotoneTypeID=0;
	if (fMonotoneType1) fMonotoneTypeID=1;
	if (fMonotoneType2) fMonotoneTypeID=2;
	if (fMonotoneType3) fMonotoneTypeID=3;

	// Call the actual mesh generator
    OfflineMap mapRemap;
	int err = GenerateOfflineMap(  mapRemap, strInputMesh, strOutputMesh, strOverlapMesh,
                                    strSourceMeta, strTargetMeta,
                                    strSourceType, strTargetType,
                                    nPin, nPout,
                                    fNoBubble, fMonotoneTypeID,
                                    fVolumetric, fNoConservation, fNoCheck,
                                    strVariables, strOutputMap, strInputData, strOutputData,
                                    strNColName, fOutputDouble, strOutputFormat, strPreserveVariables, fPreserveAll, dFillValueOverride,
                                    fSourceConcave, fTargetConcave );

	if (err) exit(err);

	return 0;
}

#endif

///////////////////////////////////////////////////////////////////////////////
