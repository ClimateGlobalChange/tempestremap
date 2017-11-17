///////////////////////////////////////////////////////////////////////////////
///
///	\file    ApplyOfflineMapExe.cpp
///	\author  Paul Ullrich
///	\version September 15, 2014
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
#include "OfflineMap.h"

#include "TempestRemapAPI.h"

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	// Input data file
	std::string strInputData;

	// Input map file
	std::string strInputMap;

	// List of variables
	std::string strVariables;

	// Input data file (second instance)
	std::string strInputData2;

	// Input map file (second instance)
	std::string strInputMap2;

	// List of variables (second instance)
	std::string strVariables2;

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
		CommandLineString(strInputData, "in_data", "");
		CommandLineString(strInputMap, "map", "");
		CommandLineString(strVariables, "var", "");
		CommandLineString(strInputData2, "in_data2", "");
		CommandLineString(strInputMap2, "map2", "");
		CommandLineString(strVariables2, "var2", "");
		CommandLineString(strOutputData, "out_data", "");
		CommandLineString(strNColName, "ncol_name", "ncol");
		CommandLineBool(fOutputDouble, "out_double");
		CommandLineString(strPreserveVariables, "preserve", "");
		CommandLineBool(fPreserveAll, "preserveall");
		CommandLineDouble(dFillValueOverride, "fillvalue", 0.0);

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Calculate metadata
	int err = ApplyOfflineMap ( strInputData, strInputMap, strVariables, strInputData2, 
								strInputMap2, strVariables2, strOutputData, strNColName, 
								fOutputDouble, strPreserveVariables, fPreserveAll, dFillValueOverride );
	if (err) exit(err);

	// Done
	AnnounceBanner();

	return 0;
}

///////////////////////////////////////////////////////////////////////////////
