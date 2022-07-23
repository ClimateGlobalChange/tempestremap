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

	// Input map file
	std::string strInputMap;

	// Options
	ApplyOfflineMapOptions optsApply;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputMap, "map", "");
		CommandLineString(optsApply.strInputData, "in_data", "");
		CommandLineString(optsApply.strOutputData, "out_data", "");
		CommandLineString(optsApply.strInputDataList, "in_data_list", "");
		CommandLineString(optsApply.strOutputDataList, "out_data_list", "");
		CommandLineString(optsApply.strVariables, "var", "");
		//CommandLineString(strInputData2, "in_data2", "");
		//CommandLineString(strInputMap2, "map2", "");
		//CommandLineString(strVariables2, "var2", "");
		CommandLineString(optsApply.strNColName, "ncol_name", "ncol");
		CommandLineString(optsApply.strEnforceBounds, "bounds", "");
		CommandLineBool(optsApply.fOutputDouble, "out_double");
		CommandLineString(optsApply.strPreserveVariables, "preserve", "");
		CommandLineBool(optsApply.fPreserveAll, "preserveall");
		CommandLineDouble(optsApply.dFillValueOverride, "fillvalue", 0.0);
		CommandLineString(optsApply.strLogDir, "logdir", "");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Calculate metadata
	int err = ApplyOfflineMap(
		strInputMap,
		optsApply );

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
