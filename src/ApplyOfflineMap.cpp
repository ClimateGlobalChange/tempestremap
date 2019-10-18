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
#include "Exception.h"
#include "OfflineMap.h"
#include "netcdfcpp.h"

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

extern "C" 
int ApplyOfflineMap(
	std::string strInputData,
	std::string strInputMap,
	std::string strVariables,
	std::string strInputData2, 
	std::string strInputMap2,
	std::string strVariables2,
	std::string strOutputData,
	std::string strNColName, 
	bool fOutputDouble,
	std::string strPreserveVariables,
	bool fPreserveAll,
	double dFillValueOverride
) {

	NcError error(NcError::silent_nonfatal);

try {

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
	mapRemap.SetFillValueOverrideDbl(dFillValueOverride);
	mapRemap.SetFillValueOverride(static_cast<float>(dFillValueOverride));

	mapRemap.Apply(
		strInputData,
		strOutputData,
		vecVariableStrings,
		strNColName,
		fOutputDouble,
		false);
	AnnounceEndBlock(NULL);

	if (strInputMap2.size()) {
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
			fOutputDouble,
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


} catch(Exception & e) {
	Announce(e.ToString().c_str());
	return (-1);

} catch(...) {
	return (-2);
}
	return 0;
}

///////////////////////////////////////////////////////////////////////////////
