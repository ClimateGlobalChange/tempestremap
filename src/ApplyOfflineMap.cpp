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
#include "netcdfcpp.h"

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

	NcError error(NcError::silent_nonfatal);

try {

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

	// Parse preserve variable list
	std::vector< std::string > vecPreserveVariableStrings;
	ParseVariableList(strPreserveVariables, vecPreserveVariableStrings);

	if (fPreserveAll && (vecPreserveVariableStrings.size() != 0)) {
		_EXCEPTIONT("--preserveall and --preserve cannot both be specified");
	}

	// Second input data file
	std::vector< std::string > vecVariableStrings2;
	if (strInputData2 != "") {
		ParseVariableList(strVariables2, vecVariableStrings2);
		if (vecVariableStrings2.size() == 0) {
			_EXCEPTIONT("No variables specified for --in_data2");
		}
		if (strInputMap2 == "") {
			_EXCEPTIONT("No map specified for --in_data2");
		}
	}
	if ((strInputMap2 != "") && (strInputData2 == "")) {
		_EXCEPTIONT("No input data specified for --map2");
	}

	// Apply OfflineMap to data
	if (strInputMap2 == "") {
		AnnounceStartBlock("Applying offline map to data");
	} else {
		AnnounceStartBlock("Applying first offline map to data");
	}

	// OfflineMap
	OfflineMap mapRemap;
	mapRemap.Read(strInputMap);
	mapRemap.SetFillValueOverride(static_cast<float>(dFillValueOverride));

	mapRemap.Apply(
		strInputData,
		strOutputData,
		vecVariableStrings,
		strNColName,
		fOutputDouble,
		false);
	AnnounceEndBlock(NULL);

	if (strInputMap2 != "") {
		AnnounceStartBlock("Applying second offline map to data");

		// OfflineMap
		OfflineMap mapRemap2;
		mapRemap2.Read(strInputMap2);

		// Verify consistency of maps
		SparseMatrix<double> & smatRemap  = mapRemap .GetSparseMatrix();
		SparseMatrix<double> & smatRemap2 = mapRemap2.GetSparseMatrix();
		if ((smatRemap.GetRows() != smatRemap2.GetRows()) ||
			(smatRemap.GetColumns() != smatRemap2.GetColumns())
		) {
			_EXCEPTIONT("Mismatch in dimensions of input maps "
				"--map and --map2");
		}

		mapRemap2.Apply(
			strInputData2,
			strOutputData,
			vecVariableStrings2,
			strNColName,
			false,
			true);

		AnnounceEndBlock(NULL);
	}

	// Copy variables from input file to output file
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

	AnnounceBanner();

	return (0);

} catch(Exception & e) {
	Announce(e.ToString().c_str());
	return (-1);

} catch(...) {
	return (-2);
}
}

///////////////////////////////////////////////////////////////////////////////


