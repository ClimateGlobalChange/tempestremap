///////////////////////////////////////////////////////////////////////////////
///
///	\file    AnalyzeMap.cpp
///	\author  Paul Ullrich
///	\version November 14, 2019
///
///	<remarks>
///		Copyright 2019 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#include "CommandLine.h"
#include "OfflineMap.h"
#include "Announce.h"

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	NcError error(NcError::silent_nonfatal);

try {
	// Input filename
	std::string strInputMap;

	// Don't check individual entries
	bool fNoCheck;

	// Monotonic map
	bool fMono;

	// Normal tolerance
	double dNormalTolerance;

	// Strict tolerance
	double dStrictTolerance;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputMap,  "map",  "");
		CommandLineBool(fMono, "mono");
		CommandLineBool(fNoCheck, "nocheck");
		CommandLineDouble(dNormalTolerance, "tol", 1.0e-8);
		CommandLineDouble(dStrictTolerance, "stricttol", 1.0e-12);
		
		ParseCommandLine(argc, argv);
	EndCommandLine(argv)
	AnnounceBanner();

	// Check parameters
	if (strInputMap == "") {
		_EXCEPTIONT("No map specified");
	}

	// OfflineMap
	OfflineMap mapRemap;
	AnnounceStartBlock("Reading map");
	mapRemap.Read(strInputMap);
	AnnounceEndBlock("Done");
	mapRemap.CheckMap(
		!fNoCheck,
		!fNoCheck,
		!fNoCheck && fMono,
		dNormalTolerance,
		dStrictTolerance);

} catch(Exception & e) {
	Announce(e.ToString().c_str());
	return (-1);

} catch(...) {
	return (-2);
}
	return 0;
}

///////////////////////////////////////////////////////////////////////////////
