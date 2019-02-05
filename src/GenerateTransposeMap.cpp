///////////////////////////////////////////////////////////////////////////////
///
///	\file    GenerateTransposeMap.cpp
///	\author  Paul Ullrich
///	\version November 16, 2018
///
///	<remarks>
///		Copyright 2000-2018 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#include "CommandLine.h"
#include "Announce.h"
#include "OfflineMap.h"

#include "netcdfcpp.h"

#include <cmath>
#include <iostream>

///////////////////////////////////////////////////////////////////////////////

typedef std::map<std::string, std::string> AttributeMap;
typedef AttributeMap::value_type AttributePair;

///////////////////////////////////////////////////////////////////////////////

void SwapAttributeNames(
	AttributeMap & mapAttributes,
	const std::string & strFirst,
	const std::string & strSecond
) {
	const int nExt = strFirst.length();
	if (nExt != strSecond.length()) {
		_EXCEPTIONT("Attribute extensions must have identical length");
	}

	// Swap src and dst attributes
	AttributeMap::iterator iterAtt = mapAttributes.begin();
	for (; iterAtt != mapAttributes.end(); iterAtt++) {
		const std::string & strName = iterAtt->first;
		if (strName.length() > nExt) {
			if (strName.substr(strName.length()-nExt) == strFirst) {
				const std::string strDstName =
					strName.substr(0, strName.length()-nExt) + strSecond;

				AttributeMap::iterator iterAttDst =
					mapAttributes.find(strDstName);

				if (iterAttDst == mapAttributes.end()) {
					mapAttributes.insert(
						AttributePair(
							strDstName,
							iterAtt->second));

				} else {
					const std::string strValue = iterAttDst->second;
					iterAttDst->second = iterAtt->second;
					iterAtt->second = strValue;
				}

			} else if (strName.substr(strName.length()-nExt) == strSecond) {
				const std::string strSrcName =
					strName.substr(0, strName.length()-nExt) + strFirst;

				AttributeMap::iterator iterAttSrc =
					mapAttributes.find(strSrcName);

				if (iterAttSrc == mapAttributes.end()) {
					mapAttributes.insert(
						AttributePair(
							strSrcName,
							iterAtt->second));
				}
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	NcError error(NcError::silent_nonfatal);

try {
	// Map file for input
	std::string strInputMapFile;

	// Overlap mesh file
	//std::string strOverlapMesh;

	// Map file for output
	std::string strOutputMapFile;

	// Do not verify the mesh
	bool fNoCheck;

	// Check monotonicity
	bool fCheckMonotone;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputMapFile, "in", "");
		//CommandLineString(strOverlapMesh, "ov_mesh", "");
		CommandLineString(strOutputMapFile, "out", "");
		CommandLineBool(fNoCheck, "nocheck");
		CommandLineBool(fCheckMonotone, "checkmono");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	// Check arguments
	if (strInputMapFile == "") {
		_EXCEPTIONT("Input map file (--in) must be specified");
	}
	//if (strOverlapMesh == "") {
	//	_EXCEPTIONT("Overlap mesh file (--in) must be specified");
	//}
	if (strOutputMapFile == "") {
		_EXCEPTIONT("Output map file (--out) must be specified");
	}

	// Atribute map
	AttributeMap mapAttributes;

	// Load map from file
	AnnounceStartBlock("Loading input map");
	OfflineMap mapIn;
        NcFile::FileFormat eFileFormat;
	mapIn.Read(strInputMapFile, eFileFormat, &mapAttributes);
	AnnounceEndBlock("Done");

	// Generate transpose map
	AnnounceStartBlock("Generating transpose map");
	OfflineMap mapOut;
	mapOut.SetTranspose(mapIn);
	AnnounceEndBlock("Done");

	// Verify map
	if (!fNoCheck) {
		AnnounceStartBlock("Verifying map");
		mapOut.IsConsistent(1.0e-8);
		mapOut.IsConservative(1.0e-8);

		if (fCheckMonotone) {
			mapOut.IsMonotone(1.0e-12);
		}
		AnnounceEndBlock("Done");
	}

	// Swap attribute names
	SwapAttributeNames(mapAttributes, "_src", "_dst");
	SwapAttributeNames(mapAttributes, "_a", "_b");

	// Find version name
	AttributeMap::iterator iterVersion = mapAttributes.find("version");
	if (iterVersion == mapAttributes.end()) {
		mapAttributes.insert(
			AttributePair("version", "GenerateTransposeMap 2.0 : 2018-11-16"));
	} else {
		iterVersion->second =
			"GenerateTransposeMap 2.0 : 2018-11-16 :: " + iterVersion->second;
	}

	// Write map to file
	AnnounceStartBlock("Writing transpose map");
	mapOut.Write(strOutputMapFile, eFileFormat, mapAttributes);
	AnnounceEndBlock("Done");

	return (0);

} catch(Exception & e) {
	Announce(e.ToString().c_str());
	return (-1);

} catch(...) {
	return (-2);
}
}


///////////////////////////////////////////////////////////////////////////////

