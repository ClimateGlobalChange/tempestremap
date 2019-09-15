///////////////////////////////////////////////////////////////////////////////
///
///	\file    RemapIterationsExe.cpp
///	\author  Vijay Mahadevan
///	\version September 4, 2019
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

///
/// Example usage:
///
///  ./RemapIterations --iterations 10 --src_mesh source.g --tgt_mesh target.g --ov_mesh overlap.g --fwdmap fwdRemapWeights.nc --revmap revRemapWeights.nc \
///                    --var "Temperature" --out_data outputdata.nc
///

#include "Announce.h"
#include "CommandLine.h"
#include "Exception.h"
#include "OfflineMap.h"
#include "netcdfcpp.h"
#include "NetCDFUtilities.h"

#include "TempestRemapAPI.h"
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

int main(int argc, char** argv) {

	// Source mesh file
	std::string strSrcMesh;

  // Target mesh file
  std::string strTgtMesh;

  // Overlap mesh file
  std::string strOvMesh;

	// Input map file
	std::string strForwardMap;

  // Reverse map file
  std::string strReverseMap;

  // Total number of remap iterations
  int iNRemapIterations;

  // Total number of remap iterations
  int iOutputFrequency;

	// List of variables
	std::string strVariables;

	// Output data file
	std::string strOutputData;

	// Name of the ncol variable
	std::string strNColName;

  // Fill value override
  double dFillValueOverride=0.0;

  int nPin, nPout;

	// Parse the command line
	BeginCommandLine()
    CommandLineString(strSrcMesh, "src_mesh", "");
    CommandLineString(strTgtMesh, "tgt_mesh", "");
    CommandLineString(strOvMesh, "ov_mesh", "");
    CommandLineInt(nPin, "np_in", 1);
    CommandLineInt(nPout, "np_out", 1);
    CommandLineString(strForwardMap, "fwdmap", "");
    CommandLineString(strReverseMap, "revmap", "");

		CommandLineInt(iNRemapIterations, "iterations", 1);
    CommandLineInt(iOutputFrequency, "frequency", 1);
		CommandLineString(strVariables, "var", "");
		
		CommandLineString(strOutputData, "out_data", "");
		CommandLineString(strNColName, "ncol_name", "ncol");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	AnnounceBanner();

	NcError error(NcError::silent_nonfatal);

  try {

    // Check parameters
    if (strSrcMesh == "") {
      _EXCEPTIONT("No source mesh specified");
    }
    if (strTgtMesh == "") {
      _EXCEPTIONT("No target mesh specified");
    }
    if (strOvMesh == "") strOvMesh = "overlapMesh.nc";
    if (strForwardMap == "") strForwardMap = "fwdRemapWeights.nc";
    if (strReverseMap == "") strReverseMap = "revRemapWeights.nc";
    if (strOutputData == "") {
      _EXCEPTIONT("No output data specified");
    }

    // Parse variable list
    std::vector< std::string > vecVariableStrings;
    ParseVariableList(strVariables, vecVariableStrings);
    const unsigned nvars = vecVariableStrings.size();

    // Apply OfflineMap to data
    if (strReverseMap == "") {
      AnnounceStartBlock("Applying offline map to data");
    } else {
      AnnounceStartBlock("Applying first offline map to data");
    }

    // Load source and target grids

    // Generate overlap meshes
    Mesh meshOverlap;
    int ierr = GenerateOverlapMesh ( strSrcMesh, strTgtMesh, 
                                     meshOverlap, strOvMesh,
                                     "Netcdf4", "exact", false );

    OfflineMap mapFwdRemap, mapRevRemap;
    // Compute the Forward map
    ierr = GenerateOfflineMap ( mapFwdRemap, strSrcMesh, strTgtMesh, strOvMesh, "", "", 
                                "fv", "fv", nPin, nPout,
                                /* fBubble */ false, 
                                /* fMonotoneTypeID */ 2,
                                /* fVolumetric */ false,
                                /* fNoConservation */ false, 
                                /* fNoCheck */ false,
                                "", strForwardMap);

    // Compute the Reverse map
    ierr = GenerateOfflineMap ( mapRevRemap, strTgtMesh, strSrcMesh, strOvMesh, "", "", 
                                "fv", "fv", nPout, nPin,
                                /* fBubble */ false, 
                                /* fMonotoneTypeID */ 2,
                                /* fVolumetric */ false,
                                /* fNoConservation */ false, 
                                /* fNoCheck */ false,
                                "", strReverseMap);

    // OfflineMap
    std::vector<DataArray1D<double> > source_solutions, target_solutions;
    // mapFwdRemap.Read(strForwardMap);
    mapFwdRemap.SetFillValueOverride(static_cast<float>(dFillValueOverride));
    mapFwdRemap.RetrieveFieldData("source", strSrcMesh, vecVariableStrings, strNColName, source_solutions);

    // mapRevRemap.Read(strReverseMap);
    mapRevRemap.SetFillValueOverride(static_cast<float>(dFillValueOverride));
    mapFwdRemap.RetrieveFieldData("target", strTgtMesh, vecVariableStrings, strNColName, target_solutions);

    // Verify consistency of maps
    SparseMatrix<double> & smatRemap  = mapFwdRemap .GetSparseMatrix();
    SparseMatrix<double> & smatRemap2 = mapRevRemap.GetSparseMatrix();
    std::cout << "FwdRemap: " << smatRemap.GetRows() << " " << smatRemap.GetColumns() << std::endl;
    std::cout << "RevRemap: " << smatRemap2.GetRows() << " " << smatRemap2.GetColumns() << std::endl;
    if ((smatRemap.GetRows() != smatRemap2.GetColumns()) ||
      (smatRemap.GetColumns() != smatRemap2.GetRows())
    ) {
      _EXCEPTIONT("Mismatch in dimensions of input maps "
        "--map and --revmap");
    }

    // Open source data file
    NcFile ncOutput(strOutputData.c_str(), NcFile::Replace, NULL, 0, NcFile::Offset64Bits);
    if (!ncOutput.is_valid()) {
      _EXCEPTION1("Cannot open output data file \"%s\"",
        strOutputData.c_str());
    }
    std::cout << "\nWriting output solution vectors to " << strOutputData << "\n";

    // Attributes
    ncOutput.add_att("Title", "TempestRemap Repetitive Remapper");

    // Map dimensions
    int nA = (int)(smatRemap.GetRows());
    int nB = (int)(smatRemap.GetColumns());

    NcDim * dimSrcGridRank = ncOutput.add_dim("n_src_dofs", nB);
    NcDim * dimDstGridRank = ncOutput.add_dim("n_dst_dofs", nA);
    NcDim * dimOutputFrequency = ncOutput.add_dim("output_frequency", iOutputFrequency);
    NcDim * dimLevels      = ncOutput.add_dim("time", iNRemapIterations/iOutputFrequency+1);
    std::vector<NcVar *> varFields(nvars*2);

    std::stringstream sstr;
    for (unsigned ivar=0; ivar < nvars; ++ivar)
    {
      sstr.str("");
      sstr << vecVariableStrings[ivar] << "_remap_src";
      std::cout << ivar << ") Source solution: " << sstr.str() << " " << source_solutions[ivar].GetRows() << std::endl;
      assert(source_solutions[ivar].GetRows() == nB);
      varFields[ivar]       = ncOutput.add_var(sstr.str().c_str(), ncDouble, dimSrcGridRank, dimLevels);
      sstr.str("");
      sstr << vecVariableStrings[ivar] << "_remap_tgt";
      std::cout << ivar << ") Target solution: " << sstr.str() << " " << target_solutions[ivar].GetRows() << std::endl;
      assert(target_solutions[ivar].GetRows() == nA);
      varFields[nvars+ivar] = ncOutput.add_var(sstr.str().c_str(), ncDouble, dimDstGridRank, dimLevels);
    }

    AnnounceEndBlock(NULL);

    // Write out the original source data to file
    {
      for (unsigned ivar=0; ivar < nvars; ++ivar)
      {
        varFields[ivar]->set_cur(0, 0);
        varFields[ivar]->put(source_solutions[ivar], nB, 1);
      }
    }

    for (int iter=0; iter < iNRemapIterations; ++iter)
    {
      std::cout << "Iteration " << iter << " ...\n";
      
      std::cout << "\tApplying forward offline map to data\n";
      for (unsigned ivar=0; ivar < nvars; ++ivar)
      {
        // Apply the offline map to the data
        smatRemap.Apply(source_solutions[ivar], target_solutions[ivar]);
        if (iter % iOutputFrequency == 0) { // Put the data in NetCDF file as needed
          varFields[nvars+ivar]->set_cur(0, iter);
          varFields[nvars+ivar]->put(target_solutions[ivar], nA, 1);
        }
      }

      std::cout << "\tApplying reverse offline map to forward data\n";
      for (unsigned ivar=0; ivar < nvars; ++ivar)
      {
        // Apply the offline map to the data
        smatRemap2.Apply(target_solutions[ivar], source_solutions[ivar]);
        if (iter % iOutputFrequency == 0) { // Put the data in NetCDF file as needed
          varFields[ivar]->set_cur(0, iter+1);
          varFields[ivar]->put(source_solutions[ivar], nB, 1);
        }
      }

      if ( (iter/iOutputFrequency) % 10) // After every 5 iterations, flush to disk
        ncOutput.sync();

    }

    ncOutput.close();

  } catch(Exception & e) {
    Announce(e.ToString().c_str());
    return (-1);

  } catch(...) {
    return (-2);
  }
	return 0;

	// Done
	AnnounceBanner();

	return 0;
}

///////////////////////////////////////////////////////////////////////////////
