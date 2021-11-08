///////////////////////////////////////////////////////////////////////////////
///
///	\file    GenerateOfflineMapExe.cpp
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
#include "OfflineMap.h"

#include "TempestRemapAPI.h"


///////////////////////////////////////////////////////////////////////////////

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

///////////////////////////////////////////////////////////////////////////////
