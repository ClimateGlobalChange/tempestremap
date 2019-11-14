///////////////////////////////////////////////////////////////////////////////
///
///   \file    ConvertSCRIPtoExodus.cpp
///   \author  Paul Ullrich
///   \version November 14, 2019
///
///   <remarks>
///      Copyright 2019 Paul Ullrich
///
///      This file is distributed as part of the Tempest source code package.
///      Permission is granted to use, copy, modify and distribute this
///      source code and its documentation under the terms of the GNU General
///      Public License.  This software is provided "as is" without express
///      or implied warranty.
///   </remarks>

#include "CommandLine.h"
#include "GridElements.h"
#include "FiniteElementTools.h"
#include "STLStringHelper.h"
#include "NetCDFUtilities.h"
#include "Exception.h"
#include "Announce.h"

#include "netcdfcpp.h"

#include <cmath>
#include <iostream>

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	NcError error(NcError::silent_nonfatal);

try {
	// Input filename
	std::string strInputFile;

	// Output scrip filename
	std::string strOutputFile;

	// Output format
	std::string strOutputFormat;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile,  "in",  "");
		CommandLineString(strOutputFile, "out", "");
		CommandLineString(strOutputFormat, "out_format", "netcdf4");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	// Check file names
	if (strInputFile  == "") {
		_EXCEPTIONT("No input file specified");  
	}
	if (strOutputFile == "") {
		_EXCEPTIONT("No output file specified"); 
	}

	//---------------------------------------------------------------------------
	//---------------------------------------------------------------------------

	// Load input mesh
	std::cout << std::endl;
	std::cout << "..Loading input mesh" << std::endl;

	Mesh meshIn(strInputFile);
	meshIn.RemoveZeroEdges(); 		// Do we need this?

	// Generate new mesh
	// std::cout << "..Generating SCRIP format data" << std::endl;
	// Mesh meshOut;

	// Write the mesh	
	std::cout << "..Writing mesh" << std::endl;

	// Output format
	STLStringHelper::ToLower(strOutputFormat);

	NcFile::FileFormat eOutputFormat =
		GetNcFileFormatFromString(strOutputFormat);
	if (eOutputFormat == NcFile::BadFormat) {
		_EXCEPTION1("Invalid \"out_format\" value (%s), "
			"expected [Classic|Offset64Bits|Netcdf4|Netcdf4Classic]",
			strOutputFormat.c_str());
	}

	meshIn.Write(strOutputFile, eOutputFormat);

	std::cout << "..Done writing" << std::endl;


	// Announce
	std::cout << "..Mesh converter exited successfully" << std::endl;
	std::cout << "=========================================================";
	std::cout << std::endl;

	return (0);

} catch(Exception & e) {
	Announce(e.ToString().c_str());
	return (-1);

} catch(...) {
	return (-2);
}
}
