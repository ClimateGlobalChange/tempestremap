///////////////////////////////////////////////////////////////////////////////
///
///	\file    ApplyOfflineMap.cpp
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

int main(int argc, char** argv) {

try {

	// Input map file
	std::string strInputMap;

	// List of variables
	std::string strVariables;

	// Input data file
	std::string strInputData;

	// Output data file
	std::string strOutputData;

	// Name of the ncol variable
	std::string strNColName;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputMap, "map", "");
		CommandLineString(strVariables, "var", "");
		CommandLineString(strInputData, "in_data", "");
		CommandLineString(strOutputData, "out_data", "");
		CommandLineString(strNColName, "ncol_name", "ncol");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Check parameters
	if (strInputMap == "") {
		_EXCEPTIONT("No map specified");
	}
	if (strInputData == "") {
		_EXCEPTIONT("No input data specified");
	}
	if (strOutputData == "") {
		_EXCEPTIONT("No output data specified");
	}

	// Parse variable list
	std::vector< std::string > vecVariableStrings;
	ParseVariableList(strVariables, vecVariableStrings);

	if ((strInputData != "") && (vecVariableStrings.size() == 0)) {
		_EXCEPTIONT("No variables specified");
	}

	// OfflineMap
	OfflineMap mapRemap;
	mapRemap.Read(strInputMap);

	// Apply OfflineMap to data
	DataVector<double> vecDummyAreas;

	AnnounceStartBlock("Applying offline map to data");
	mapRemap.Apply(
		vecDummyAreas,
		vecDummyAreas,
		strInputData,
		strOutputData,
		vecVariableStrings,
		strNColName,
		false);
	AnnounceEndBlock(NULL);

	AnnounceBanner();

} catch(Exception & e) {
	Announce(e.ToString().c_str());
}
}

///////////////////////////////////////////////////////////////////////////////


