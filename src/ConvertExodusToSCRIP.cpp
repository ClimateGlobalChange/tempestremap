///////////////////////////////////////////////////////////////////////////////
///
///   \file    ConvertExodusToSCRIP.cpp
///   \author  Paul Ullrich
///   \version September 18, 2015
///
///   <remarks>
///      Copyright 2000-2014 Paul Ullrich
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

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile,  "in",  "");
		CommandLineString(strOutputFile, "out", "");
		
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
	NcFile::FileFormat eOutputFormat = NcFile::Offset64Bits;
	meshIn.WriteScrip(strOutputFile, eOutputFormat);

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