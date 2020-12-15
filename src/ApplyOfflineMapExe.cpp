///////////////////////////////////////////////////////////////////////////////
///
///	\file    ApplyOfflineMapExe.cpp
///	\author  Paul Ullrich
///	\version February 14, 2020
///
///	<remarks>
///		Copyright 2020 Paul Ullrich
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
#include "OfflineMap.h"

#include "TempestRemapAPI.h"

#if defined(TEMPEST_MPIOMP)
#include <mpi.h>
#endif

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

#if defined(TEMPEST_MPIOMP)
	// Initialize MPI
	MPI_Init(&argc, &argv);

	// Enable output only on rank zero
	AnnounceOnlyOutputOnRankZero();
#endif

	// Input data file
	std::string strInputData;

	// Input data file list
	std::string strInputDataList;

	// Input map file
	std::string strInputMap;

	// List of variables
	std::string strVariables;

	// Input data file (second instance)
	//std::string strInputData2;

	// Input map file (second instance)
	//std::string strInputMap2;

	// List of variables (second instance)
	//std::string strVariables2;

	// Output data file
	std::string strOutputData;

	// Output data file list
	std::string strOutputDataList;

	// Name of the ncol variable
	std::string strNColName;
	
	std::string strInputMesh;
	
	std::string strOutputMesh;
	
	std::string strOverlapMesh;
	
	int nPin;
	
	int nPout;

	// Output as double
	bool fOutputDouble;

	// List of variables to preserve
	std::string strPreserveVariables;

	// Preserve all non-remapped variables
	bool fPreserveAll;

	// Fill value override
	double dFillValueOverride;

	// Log directory
	std::string strLogDir;
	
	double lb;
	
	double ub;
	
	bool fCAAS;
	
	bool fCAASLocal;
	
	bool fContainsConcaveFaces;
	
	bool fGLL;
	
	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputData, "in_data", "");
		CommandLineString(strInputDataList, "in_data_list", "");
		CommandLineString(strInputMap, "map", "");
		CommandLineString(strVariables, "var", "");
		//CommandLineString(strInputData2, "in_data2", "");
		//CommandLineString(strInputMap2, "map2", "");
		//CommandLineString(strVariables2, "var2", "");
		CommandLineString(strOutputData, "out_data", "");
		CommandLineString(strOutputDataList, "out_data_list", "");
		CommandLineString(strNColName, "ncol_name", "ncol");
		CommandLineString(strInputMesh,"mesh_in","");
		CommandLineString(strOutputMesh,"mesh_out","");
		CommandLineString(strOverlapMesh,"mesh_ov","");
		CommandLineInt(nPin, "in_np", 0);
		CommandLineInt(nPout, "out_np", 0);
		CommandLineBool(fOutputDouble, "out_double");
		CommandLineString(strPreserveVariables, "preserve", "");
		CommandLineBool(fPreserveAll, "preserveall");
		CommandLineDouble(dFillValueOverride, "fillvalue", 0.0);
		CommandLineString(strLogDir, "logdir", "");
		CommandLineDouble(lb, "lb", 0.0);
		CommandLineDouble(ub, "ub", 1.0);
		CommandLineBool(fCAAS,"f_CAAS");
		CommandLineBool(fCAASLocal,"CAAS_local");
		CommandLineBool(fContainsConcaveFaces,"concave");
		CommandLineBool(fGLL,"gll");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Calculate metadata
	int err = ApplyOfflineMap(
		strInputData,
		strInputDataList,
		strInputMap,
		strVariables,
		strOutputData,
		strOutputDataList,
		strNColName,
		strInputMesh,
		strOutputMesh,
		strOverlapMesh,
		nPin,
		nPout,
		fOutputDouble,
		strPreserveVariables,
		fPreserveAll,
		dFillValueOverride,
		strLogDir,
		lb,
		ub,
		fCAAS,
		fCAASLocal,
		fContainsConcaveFaces,
		fGLL );

	// Done
	AnnounceBanner();

#if defined(TEMPEST_MPIOMP)
	// Deinitialize MPI
	MPI_Finalize();
#endif

	if (err) exit(err);

	return 0;
}

///////////////////////////////////////////////////////////////////////////////
