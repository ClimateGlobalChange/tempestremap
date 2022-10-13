///////////////////////////////////////////////////////////////////////////////
///
///	\file    ApplyOfflineMap.cpp
///	\author  Paul Ullrich
///	\version Feburary 14, 2020
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
#include "Exception.h"
#include "TempestRemapAPI.h"
#include "OfflineMap.h"
#include "netcdfcpp.h"
#include "GridElements.h"
#include "FiniteElementTools.h"

#include <fstream>

#if defined(TEMPEST_MPIOMP)
#include <mpi.h>
#endif

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

///	<summary>
///		Parse the list of input files.
///	</summary>
void ParseInputFiles(
	const std::string & strInputFile,
	std::vector<NcFile *> & vecFiles
) {
	int iLast = 0;
	for (int i = 0; i <= strInputFile.length(); i++) {
		if ((i == strInputFile.length()) ||
		    (strInputFile[i] == ';')
		) {
			std::string strFile =
				strInputFile.substr(iLast, i - iLast);

			NcFile * pNewFile = new NcFile(strFile.c_str());

			if (!pNewFile->is_valid()) {
				_EXCEPTION1("Cannot open input file \"%s\"",
					strFile.c_str());
			}

			vecFiles.push_back(pNewFile);
			iLast = i+1;
		}
	}

	if (vecFiles.size() == 0) {
		_EXCEPTION1("No input files found in \"%s\"",
			strInputFile.c_str());
	}
}

///////////////////////////////////////////////////////////////////////////////

extern "C" 
int ApplyOfflineMap(
	std::string strInputMap,
	const ApplyOfflineMapOptions & optsApply
) {

	NcError error(NcError::silent_nonfatal);

try {

	// Check parameters
	if (strInputMap == "") {
		_EXCEPTIONT("No map specified");
	}
	if ((optsApply.strInputData == "") && (optsApply.strInputDataList == "")) {
		_EXCEPTIONT("No input data (--in_data) or (--in_data_list) specified");
	}
	if ((optsApply.strInputData != "") && (optsApply.strInputDataList != "")) {
		_EXCEPTIONT("Only one of --in_data or --in_data_list may be specified");
	}
	if ((optsApply.strOutputData == "") && (optsApply.strOutputDataList == "")) {
		_EXCEPTIONT("No output data (--out_data) or (--out_data_list)");
	}
	if ((optsApply.strOutputData != "") && (optsApply.strOutputDataList != "")) {
		_EXCEPTIONT("Only one of --out_data or --out_data_list may be specified");
	}
	if ((optsApply.strInputData != "") && (optsApply.strOutputData == "")) {
		_EXCEPTIONT("If --in_data is specified then --out_data must also be specified");
	}
	if ((optsApply.strInputDataList != "") && (optsApply.strOutputDataList == "")) {
		_EXCEPTIONT("If --in_data_list is specified then --out_data_list must also be specified");
	}

	// Load input file list
	std::vector<std::string> vecInputDataFiles;

	if (optsApply.strInputData.length() != 0) {
		vecInputDataFiles.push_back(optsApply.strInputData);

	} else {
		std::ifstream ifInputFileList(optsApply.strInputDataList.c_str());
		if (!ifInputFileList.is_open()) {
			_EXCEPTION1("Unable to open file \"%s\"",
				optsApply.strInputDataList.c_str());
		}
		std::string strFileLine;
		while (std::getline(ifInputFileList, strFileLine)) {
			if (strFileLine.length() == 0) {
				continue;
			}
			if (strFileLine[0] == '#') {
				continue;
			}
			vecInputDataFiles.push_back(strFileLine);
		}
	}

	// Load output file list
	std::vector<std::string> vecOutputDataFiles;

	if (optsApply.strOutputData.length() != 0) {
		vecOutputDataFiles.push_back(optsApply.strOutputData);

	} else {
		std::ifstream ifOutputFileList(optsApply.strOutputDataList.c_str());
		if (!ifOutputFileList.is_open()) {
			_EXCEPTION1("Unable to open file \"%s\"",
				optsApply.strOutputDataList.c_str());
		}
		std::string strFileLine;
		while (std::getline(ifOutputFileList, strFileLine)) {
			if (strFileLine.length() == 0) {
				continue;
			}
			if (strFileLine[0] == '#') {
				continue;
			}
			vecOutputDataFiles.push_back(strFileLine);
		}
	}

	// Check length
	if (vecInputDataFiles.size() != vecOutputDataFiles.size()) {
		_EXCEPTIONT("Mistmatch in --in_data_list and --out_data_list file length");
	}

	// Parse variable list
	std::vector< std::string > vecVariableStrings;
	ParseVariableList(optsApply.strVariables, vecVariableStrings);

	// Parse preserve variable list
	std::vector< std::string > vecPreserveVariableStrings;
	ParseVariableList(optsApply.strPreserveVariables, vecPreserveVariableStrings);

	if (optsApply.fPreserveAll && (vecPreserveVariableStrings.size() != 0)) {
		_EXCEPTIONT("--preserveall and --preserve cannot both be specified");
	}

#if defined(TEMPEST_MPIOMP)
	// Spread files across nodes
	int nMPIRank;
	MPI_Comm_rank(MPI_COMM_WORLD, &nMPIRank);

	int nMPISize;
	MPI_Comm_size(MPI_COMM_WORLD, &nMPISize);

	Announce("Executing ApplyOfflineMap with %i threads over %i files",
		nMPISize, vecInputDataFiles.size());
#endif

#if defined(TEMPEST_MPIOMP)
	// Set up logging
	if (optsApply.strLogDir == "") {
		Announce("Reporting only enabled on thread 0 (if reporting desired, use --logdir)");
	} else {
		Announce("Logs will be written to directory \"%s\"", optsApply.strLogDir.c_str());
	}

	// Open log file
	if (optsApply.strLogDir != "") {
		std::string strLogFile = optsApply.strLogDir;
		if (strLogFile[strLogFile.length()-1] != '/') {
			strLogFile += "/";
		}
		char szTemp[20];
		sprintf(szTemp, "log%06i.txt", nMPIRank);
		strLogFile += szTemp;

		FILE * fp = fopen(strLogFile.c_str(), "w");
		if (fp == NULL) {
			_EXCEPTION1("Unable to open output log file \"%s\"",
				strLogFile.c_str());
		}
		AnnounceSetOutputBuffer(fp);
		AnnounceOutputOnAllRanks();
	}
#else
	if (optsApply.strLogDir != "") {
		std::string strLogFile = optsApply.strLogDir;
		if (strLogFile[strLogFile.length()-1] != '/') {
			strLogFile += "/";
		}
		strLogFile += "log000000.txt";

		Announce("Log will be written to file \"%s\"", strLogFile.c_str());

		FILE * fp = fopen(strLogFile.c_str(), "w");
		if (fp == NULL) {
			_EXCEPTION1("Unable to open output log file \"%s\"",
				strLogFile.c_str());
		}
		AnnounceSetOutputBuffer(fp);
	}
#endif
/*
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
*/

	DataArray3D<int> dataGLLNodesIn;
	DataArray3D<double> dataGLLJacobianIn;
	
	DataArray3D<int> dataGLLNodesOut;
	DataArray3D<double> dataGLLJacobianOut;
	
	DataArray3D<int> * pdataGLLNodesIn = NULL;
			
	DataArray3D<int> * pdataGLLNodesOut = NULL;
							
	Mesh meshOverlap;
	Mesh meshSource;
	Mesh meshTarget;
	
	Mesh * pmeshSource = NULL;
	Mesh * pmeshOverlap = NULL;
		
		
	if(optsApply.strInputMesh != ""){
		
		//If the source mesh is finite volume, we need the source mesh for local p bounds preservation
		
		meshSource.Read(optsApply.strInputMesh);
		meshSource.ConstructReverseNodeArray();
		meshSource.ConstructEdgeMap();
		pmeshSource = &meshSource;
		
		
		//If the source mesh is finite element, we need dataGLLNodes for local bounds preservation
		if(optsApply.fgll){
			
			double dTotalAreaInput = meshSource.CalculateFaceAreas(optsApply.fContainsConcaveFaces);
			double dNumericalAreaIn = GenerateMetaData(meshSource,optsApply.nPin,false,dataGLLNodesIn,dataGLLJacobianIn);
			
			pdataGLLNodesIn = &dataGLLNodesIn;
			
		}
		
	}
	
	if(optsApply.strOverlapMesh != ""){
		
		meshOverlap.Read(optsApply.strOverlapMesh);
		meshOverlap.RemoveZeroEdges();
		pmeshOverlap = &meshOverlap;
		
	}
	
	//If the target mesh is finite element, get dataGLLNodesOut
	if(optsApply.strOutputMesh != ""){
		
		meshTarget.Read(optsApply.strOutputMesh);
		meshTarget.ConstructReverseNodeArray();
		meshTarget.ConstructEdgeMap();
		
		double dTotalAreaOutput = meshTarget.CalculateFaceAreas(optsApply.fContainsConcaveFaces);
		double dNumericalAreaOut = GenerateMetaData(meshTarget,optsApply.nPout,false,dataGLLNodesOut,dataGLLJacobianOut);
		
		pdataGLLNodesOut = &dataGLLNodesOut;
					
	}
		
	// OfflineMap
	AnnounceStartBlock("Loading offline map");

	OfflineMap mapRemap;
	mapRemap.Read(strInputMap);
	mapRemap.SetFillValueOverrideDbl(optsApply.dFillValueOverride);
	mapRemap.SetFillValueOverride(static_cast<float>(optsApply.dFillValueOverride));
	mapRemap.SetEnforcementBounds(optsApply.strEnforceBounds,pmeshSource,pmeshOverlap,pdataGLLNodesIn,
								  pdataGLLNodesOut,optsApply.nPin);
	//mapRemap.SetEnforcementBounds(optsApply.strEnforceBounds);

	AnnounceEndBlock("Done");

	for (int f = 0; f < vecInputDataFiles.size(); f++) {

#if defined(TEMPEST_MPIOMP)
		if (f % nMPISize != nMPIRank) {
			continue;
		}
#endif
		AnnounceStartBlock("Processing \"%s\"", vecInputDataFiles[f].c_str());

		// Apply the map
		mapRemap.Apply(
			vecInputDataFiles[f],
			vecOutputDataFiles[f],
			vecVariableStrings,
			optsApply.strNColName,
			optsApply.fOutputDouble,
			false);
	
		// Copy variables from input file to output file
		if (optsApply.fPreserveAll) {
			AnnounceStartBlock("Preserving variables");
			mapRemap.PreserveAllVariables(
				vecInputDataFiles[f],
				vecOutputDataFiles[f]);
			AnnounceEndBlock(NULL);

		} else if (vecPreserveVariableStrings.size() != 0) {
			AnnounceStartBlock("Preserving variables");
			mapRemap.PreserveVariables(
				vecInputDataFiles[f],
				vecOutputDataFiles[f],
				vecPreserveVariableStrings);
			AnnounceEndBlock(NULL);
		}

		AnnounceEndBlock("Done");
	}
	AnnounceEndBlock(NULL);
/*
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
*/

	AnnounceOnlyOutputOnRankZero();
	AnnounceSetOutputBuffer(stdout);

} catch(Exception & e) {
	Announce(e.ToString().c_str());
	return (-1);

} catch(...) {
	return (-2);
}

	return 0;
}

///////////////////////////////////////////////////////////////////////////////
