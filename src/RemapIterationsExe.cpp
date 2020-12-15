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

// #define ENABLE_EIGEN3

#include "Announce.h"
#include "CommandLine.h"
#include "Exception.h"
#include "OfflineMap.h"
#include "netcdfcpp.h"
#include "NetCDFUtilities.h"
#include "FiniteVolumeTools.h"

#ifdef ENABLE_EIGEN3
#include <Eigen/Dense>
#include <Eigen/Sparse>
#endif

#include "TempestRemapAPI.h"

#ifdef ENABLE_EIGEN3

typedef Eigen::Map< Eigen::Matrix< double, 1, Eigen::Dynamic > > WeightRowVector;
typedef Eigen::Map< Eigen::Matrix< double, Eigen::Dynamic, 1 > > WeightColVector;
typedef Eigen::SparseMatrix< double, Eigen::RowMajor > WeightMatrix;

void copy_tempest_sparsemat_to_eigen3(SparseMatrix<double>& trmat, WeightMatrix& emat)
{
    /* Should the columns be the global size of the matrix ? */
    emat.resize( trmat.GetRows(), trmat.GetColumns() );

    DataArray1D< int > lrows;
    DataArray1D< int > lcols;
    DataArray1D< double > lvals;
    trmat.GetEntries( lrows, lcols, lvals );
    unsigned locvals = lvals.GetRows();

    emat.reserve( locvals );
    for( unsigned iv = 0; iv < locvals; iv++ )
    {
        // std::cout << "Row = " << row_ldofmap[lrows[iv]] << ", Col = " << col_ldofmap[lcols[iv]]
        // << ", DATA = " << lvals[iv] << std::endl; std::cout << "Row = " << lrows[iv] << ", Col =
        // " << lcols[iv] << ", DATA = " << lvals[iv] << std::endl;
        emat.insert( lrows[iv], lcols[iv] ) = lvals[iv];
    }

    emat.makeCompressed();

    return;
}

#endif

// For CAAS
double ApplyCAASLimiting( OfflineMap& mapOperator, Mesh& meshInput, Mesh& meshOverlap, const int nPin,
                        DataArray1D< double >& dataInDouble, DataArray1D< double >& dataOutDouble, bool useCAAS,
                        bool useCAASLocal );


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


void print_vec(DataArray1D<double>& vec)
{
  printf("Vec: %f %f %f %f %f\n", vec[0], vec[1], vec[2], vec[3], vec[4]);
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

  // Specify if we are restarting starting from provided index
  int iRestartOffset;

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

  bool skipMapGen = false;

  bool useCAAS = false;

  int nPin, nPout;

	// Parse the command line
	BeginCommandLine()
    CommandLineString(strSrcMesh, "src_mesh", "");
    CommandLineString(strTgtMesh, "tgt_mesh", "");
    CommandLineString(strOvMesh, "ov_mesh", "");
    CommandLineBool(skipMapGen, "skip");
    CommandLineInt(nPin, "np_in", 1);
    CommandLineInt(nPout, "np_out", 1);
    CommandLineString(strForwardMap, "fwdmap", "");
    CommandLineString(strReverseMap, "revmap", "");
		CommandLineInt(iRestartOffset, "restart", 0);
    CommandLineBool(useCAAS, "CAAS");

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

    // Optionally load the input and output meshes
    Mesh meshInput, meshOutput;

    // Generate overlap meshes
    Mesh meshOverlap;
    if (!skipMapGen)
    {
      int ierr = GenerateOverlapMesh ( strSrcMesh, strTgtMesh, 
                                       meshOverlap, strOvMesh,
                                       "Netcdf4", "exact", false );
    }

    // Load source and target grids
    if (useCAAS)
    {
      meshInput.Read ( strSrcMesh );
      meshOutput.Read ( strTgtMesh );
      if (skipMapGen) {
        meshOverlap.Read ( strOvMesh );
      }
    }

    OfflineMap mapFwdRemap, mapRevRemap;
    if (skipMapGen)
    {
      mapFwdRemap.Read(strForwardMap);
      mapRevRemap.Read(strReverseMap);
    }
    else
    {
      // Compute the Forward map
      int ierr = GenerateOfflineMap ( mapFwdRemap, strSrcMesh, strTgtMesh, strOvMesh, "", "", 
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
    }

    // Verify consistency of maps
    SparseMatrix<double> & smatRemap  = mapFwdRemap.GetSparseMatrix();
    SparseMatrix<double> & smatRemapReverse = mapRevRemap.GetSparseMatrix();
    std::cout << "FwdRemap: " << smatRemap.GetRows() << " " << smatRemap.GetColumns() << std::endl;
    std::cout << "RevRemap: " << smatRemapReverse.GetRows() << " " << smatRemapReverse.GetColumns() << std::endl;
    if ((smatRemap.GetRows() != smatRemapReverse.GetColumns()) ||
      (smatRemap.GetColumns() != smatRemapReverse.GetRows())
    ) {
      _EXCEPTIONT("Mismatch in dimensions of input maps "
        "--map and --revmap");
    }

#ifdef ENABLE_EIGEN3
    WeightMatrix eigFwdMap, eigRevMap;
    copy_tempest_sparsemat_to_eigen3(smatRemap, eigFwdMap);
    copy_tempest_sparsemat_to_eigen3(smatRemapReverse, eigRevMap);
#endif

    int nA = (int)(smatRemap.GetRows());
    int nB = (int)(smatRemap.GetColumns());

    // OfflineMap
    std::vector<DataArray1D<double> > source_solutions, target_solutions;
    // mapFwdRemap.Read(strForwardMap);
    mapFwdRemap.SetFillValueOverride(static_cast<double>(dFillValueOverride));
    mapRevRemap.SetFillValueOverride(static_cast<double>(dFillValueOverride));
    
    NcFile::FileMode fmode = NcFile::Replace;
    NcFile::FileFormat fformat = NcFile::/*Offset64Bits*//*Netcdf4*/Netcdf4Classic;
    // populate starting vectors for source and target
    if (iRestartOffset == 0) {
       // Read the solution data from source grid file
       mapFwdRemap.RetrieveFieldData("source", strSrcMesh, vecVariableStrings, strNColName, source_solutions);

       // Read the solution data from target grid file
       mapFwdRemap.RetrieveFieldData("target", strTgtMesh, vecVariableStrings, strNColName, target_solutions);
    }
    else { // performing a restart
      fmode = NcFile::Write;
      
      // Allocate memory
      source_solutions.resize(nvars);
      target_solutions.resize(nvars);
      for (unsigned ivar=0; ivar < nvars; ++ivar) {
         source_solutions[ivar].Allocate(nB);
         target_solutions[ivar].Allocate(nA);
      }
    }

    // Open source data file
    NcFile ncOutput(strOutputData.c_str(), fmode, NULL, 0, fformat);
    if (!ncOutput.is_valid()) {
      _EXCEPTION1("Cannot open output data file \"%s\"", strOutputData.c_str());
    }
    std::cout << "\nWriting output solution vectors to " << strOutputData << "\n";

    // Attributes
    ncOutput.add_att("Title", "TempestRemap Repetitive Remapper");

    // Map dimensions
    NcDim * dimSrcGridRank = (iRestartOffset == 0 ? ncOutput.add_dim("n_src_dofs", nB) : ncOutput.get_dim("n_src_dofs"));
    NcDim * dimDstGridRank = (iRestartOffset == 0 ? ncOutput.add_dim("n_dst_dofs", nA) : ncOutput.get_dim("n_dst_dofs"));
    NcDim * dimOutputFrequency = (iRestartOffset == 0 ? ncOutput.add_dim("output_frequency", iOutputFrequency) : ncOutput.get_dim("output_frequency"));
    //NcDim * dimLevels      = (iRestartOffset == 0 ? ncOutput.add_dim("time", iNRemapIterations/iOutputFrequency+1) : ncOutput.get_dim("time"));
    NcDim * dimLevels      = (iRestartOffset == 0 ? ncOutput.add_dim("time", NC_UNLIMITED) : ncOutput.get_dim("time"));
    std::vector<NcVar *> varFields(nvars*2);

    std::stringstream sstr;
    for (unsigned ivar=0; ivar < nvars; ++ivar)
    {
      sstr.str("");
      sstr << vecVariableStrings[ivar] << "_remap_src";
      std::cout << ivar << ") Source solution: " << sstr.str() << " " << source_solutions[ivar].GetRows() << std::endl;
      if (iRestartOffset == 0) {
         varFields[ivar]       = ncOutput.add_var(sstr.str().c_str(), ncDouble, dimLevels, dimSrcGridRank);
         print_vec(source_solutions[ivar]);
         assert(source_solutions[ivar].GetRows() == nB);

      }
      else {
         varFields[ivar]       = ncOutput.get_var(sstr.str().c_str());
      }
      sstr.str("");
      sstr << vecVariableStrings[ivar] << "_remap_tgt";
      std::cout << ivar << ") Target solution: " << sstr.str() << " " << target_solutions[ivar].GetRows() << std::endl;
      if (iRestartOffset == 0) {
         varFields[nvars+ivar] = ncOutput.add_var(sstr.str().c_str(), ncDouble, dimLevels, dimDstGridRank);
         print_vec(target_solutions[ivar]);
         assert(target_solutions[ivar].GetRows() == nA);
      }
      else {
         varFields[nvars+ivar] = ncOutput.get_var(sstr.str().c_str());
      }
    }

    AnnounceEndBlock(NULL);

    // Write out the original source data to file
    if (iRestartOffset == 0)
    {
      for (unsigned ivar=0; ivar < nvars; ++ivar)
      {
        /* Write initial source sampled solution data */
        varFields[ivar]->set_cur(0, 0);
        varFields[ivar]->put(source_solutions[ivar], 1, nB);
        
        /* Write initial target sampled solution data */
        varFields[nvars+ivar]->set_cur(0, 0);
        varFields[nvars+ivar]->put(target_solutions[ivar], 1, nA);
      }
    }
    else {
      for (unsigned ivar=0; ivar < nvars; ++ivar)
      {
        /* Write initial source sampled solution data */
        varFields[ivar]->set_cur(iRestartOffset, 0);
        varFields[ivar]->get(source_solutions[ivar], 1, nB);
        print_vec(source_solutions[ivar]);
        
        /* Write initial target sampled solution data */
        varFields[nvars+ivar]->set_cur(iRestartOffset, 0);
        varFields[nvars+ivar]->get(target_solutions[ivar], 1, nA);
        print_vec(target_solutions[ivar]);
      }
    }

    for (int iter=iRestartOffset+1; iter <= iNRemapIterations; ++iter)
    {
      std::cout << "Iteration " << iter << " ...\n";
      
      // std::cout << "\tApplying forward offline map to data\n";
      for (unsigned ivar=0; ivar < nvars; ++ivar)
      {
        // Apply the offline map to the data
#ifdef ENABLE_EIGEN3
        WeightColVector cVarSVec(source_solutions[ivar], nB);
        // WeightColVector cVarTVec(target_solutions[ivar], nA);
        Eigen::VectorXd cVarTVec = eigFwdMap * cVarSVec;
        //for(unsigned il=0; il < nA; ++il) target_solutions[ivar][il] = cVarTVec[il];
        const double* tdata = cVarTVec.data();
#else

        smatRemap.Apply(source_solutions[ivar], target_solutions[ivar]);

        if (useCAAS) {
          // double ApplyCAASLimiting( OfflineMap& mapOperator, Mesh& meshInput, Mesh& meshOverlap, const int nPin,
          //                DataArray1D< double >& dataInDouble, DataArray1D< double >& dataOutDouble, bool useCAAS,
          //                bool useCAASLocal );
          double glb_error_fwd = ApplyCAASLimiting ( mapFwdRemap, meshInput, meshOverlap, nPin, 
                                                     source_solutions[ivar], target_solutions[ivar],
                                                     false, true );
        }


        const double* tdata = (double*)(target_solutions[ivar]);
#endif
        if (iter == 1 || iter % iOutputFrequency == 0) { // Put the data in NetCDF file as needed
          if (ivar == 0) std::cout << iter << " -- Writing target solution at " << (iter)/iOutputFrequency+1 << "\n";
          if (iter == 1 || iOutputFrequency == 1) varFields[nvars+ivar]->set_cur(iter, 0);
          else varFields[nvars+ivar]->set_cur((iter)/iOutputFrequency+1, 0);
          varFields[nvars+ivar]->put(tdata, 1, nA);
        }
      }

      // std::cout << "\tApplying reverse offline map to forward data\n";
      for (unsigned ivar=0; ivar < nvars; ++ivar)
      {
        // Apply the offline map to the data
#ifdef ENABLE_EIGEN3
        // WeightColVector cVarSVec(source_solutions[ivar], nB);
        WeightColVector cVarTVec(target_solutions[ivar], nA);
        Eigen::VectorXd cVarSVec = eigRevMap * cVarTVec;
        //for(unsigned il=0; il < nB; ++il) source_solutions[ivar][il] = cVarSVec[il];
        const double* sdata = cVarSVec.data();
#else
        // Apply the reverse map
        smatRemapReverse.Apply(target_solutions[ivar], source_solutions[ivar]);

        if (useCAAS) {
          // double ApplyCAASLimiting( OfflineMap& mapOperator, Mesh& meshInput, Mesh& meshOverlap, const int nPin,
          //                DataArray1D< double >& dataInDouble, DataArray1D< double >& dataOutDouble, bool useCAAS,
          //                bool useCAASLocal );
          double glb_error_fwd = ApplyCAASLimiting ( mapRevRemap, meshOutput, meshOverlap, nPin, 
                                                     target_solutions[ivar], source_solutions[ivar],
                                                     false, true );
        }

        const double* sdata = (double*)(source_solutions[ivar]);
#endif
        if (iter == 1 || iter % iOutputFrequency == 0) { // Put the data in NetCDF file as needed
          if (ivar == 0) std::cout << iter << " -- Writing source solution at " << (iter)/iOutputFrequency+1 << "\n";
          if (iter == 1 || iOutputFrequency == 1) varFields[ivar]->set_cur(iter, 0);
          else varFields[ivar]->set_cur((iter)/iOutputFrequency+1, 0);
          varFields[ivar]->put(sdata, 1, nB);
        }
      }

      if ( (iter/iOutputFrequency) % 5) // After every 5 iterations, flush to disk
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


double ApplyCAASLimiting( OfflineMap& mapOperator, Mesh& meshInput, Mesh& meshOverlap, const int nPin,
                        DataArray1D< double >& dataInDouble, DataArray1D< double >& dataOutDouble, bool useCAAS,
                        bool useCAASLocal )
{

    const int nSourceCount                      = dataInDouble.GetRows();
    const int nTargetCount                      = dataOutDouble.GetRows();
    const DataArray1D< double >& m_dSourceAreas = mapOperator.GetSourceAreas();
    const DataArray1D< double >& m_dTargetAreas = mapOperator.GetTargetAreas();

    // Announce input mass
    double dSourceMass = 0.0;
    double dSourceMin  = dataInDouble[0];
    double dSourceMax  = dataInDouble[0];
    for( int i = 0; i < nSourceCount; i++ )
    {
        dSourceMass += dataInDouble[i] * m_dSourceAreas[i];
        if( dataInDouble[i] < dSourceMin ) { dSourceMin = dataInDouble[i]; }
        if( dataInDouble[i] > dSourceMax ) { dSourceMax = dataInDouble[i]; }
    }

    // Apply the offline map to the data

    if( useCAASLocal || useCAAS )
    {
        DataArray1D< double > l = dataOutDouble;
        DataArray1D< double > u = dataOutDouble;
        DataArray1D< double > x( nTargetCount );
        double b = dSourceMass;

        for( int i = 0; i < l.GetRows(); i++ )
        {
            b -= dataOutDouble[i] * m_dTargetAreas[i];
        }

        if( useCAASLocal )
        {
            int GLLSizeIn  = 0; // FV
            int GLLSizeOut = 0; // FV
            // int pOut       = dataGLLNodesOut.GetSize( 0 );
            // int qOut       = dataGLLNodesOut.GetSize( 1 );
            // int pIn        = dataGLLNodesIn.GetSize( 0 );
            // int qIn        = dataGLLNodesIn.GetSize( 1 );
            int pIn          = nPin;
            double f_maxI  = 0.0;
            double f_minI  = 0.0;
            int nTargetFaces = nTargetCount;

            std::vector< std::vector< int > > SourceOvTarget( nTargetFaces );

            for( int i = 0; i < meshOverlap.faces.size(); i++ )
            {
                int ixT = meshOverlap.vecTargetFaceIx[i];
                int ixS = meshOverlap.vecSourceFaceIx[i];
                SourceOvTarget[ixT].push_back( ixS );
            }

            std::vector< double > local_UB( nTargetCount );
            std::vector< double > local_LB( nTargetCount );

            for( int i = 0; i < nTargetCount; i++ )
            {
                AdjacentFaceVector vecAdjFaces;

                GetAdjacentFaceVectorByEdge( meshInput, SourceOvTarget[i][0], ( nPin + 1 ) * ( nPin + 1 ),
                                             vecAdjFaces );
                f_maxI = dataInDouble[vecAdjFaces[0].first];
                f_minI = dataInDouble[vecAdjFaces[0].first];
                for( int j = 0; j < vecAdjFaces.size(); j++ )
                {
                    int k  = vecAdjFaces[j].first;
                    f_maxI = fmax( f_maxI, dataInDouble[k] );
                    f_minI = fmin( f_minI, dataInDouble[k] );
                }

                for( int j = 0; j < SourceOvTarget[i].size(); j++ )
                {

                    int k  = SourceOvTarget[i][j];
                    f_maxI = fmax( f_maxI, dataInDouble[k] );
                    f_minI = fmin( f_minI, dataInDouble[k] );
                }

                // f_minI=fmax(f_minI,0.0);

                local_UB[i] = f_maxI;
                local_LB[i] = f_minI;
            }

            double mt = 0.0;
            for( int i = 0; i < nTargetCount; i++ )
            {
                mt += m_dTargetAreas[i] * ( local_LB[i] - dataOutDouble[i] );
            }

            for( int i = 0; i < l.GetRows(); i++ )
            {
                l[i] = local_LB[i] - l[i];
                u[i] = local_UB[i] - u[i];
            }

            // Adjust mass of lower bound if greater than b
            double mL = 0.0;

            for( int i = 0; i < nTargetCount; i++ )
            {
                mL += m_dTargetAreas[i] * l[i];
            }
            if( mL > b )
            {
                for( int i = 0; i < nTargetCount; i++ )
                {
                    mL   = mL - m_dTargetAreas[i] * l[i] + m_dTargetAreas[i] * ( dSourceMin - dataOutDouble[i] );
                    l[i] = dSourceMin - dataOutDouble[i];
                    if( mL < b ) { break; }
                }
            }

            // Adjust mass of upper bound if less than b
            double mU = 0.0;

            if( mU < b )
            {
                for( int i = 0; i < nTargetCount; i++ )
                {
                    mU   = mU - m_dTargetAreas[i] * u[i] + m_dTargetAreas[i] * ( dSourceMax - dataOutDouble[i] );
                    u[i] = dSourceMax - dataOutDouble[i];
                    if( mU > b ) { break; }
                }
            }
        } // if( useCAASLocal )
        else if( useCAAS )
        {
            for( int i = 0; i < l.GetRows(); i++ )
            {
                l[i] = lb - l[i];
                u[i] = ub - u[i];
            }
        }

        // Invoke CAAS application on the offline map
        mapOperator.CAAS( x, l, u, b );

        // Add correction
        for( int i = 0; i < l.GetRows(); i++ )
        {
            dataOutDouble[i] += x[i];
        }

    }

    // Announce output mass
    double dTargetMass = 0.0;
    double dTargetMin  = dataOutDouble[0];
    double dTargetMax  = dataOutDouble[0];
    for( int i = 0; i < nTargetCount; i++ )
    {
        dTargetMass += dataOutDouble[i] * m_dTargetAreas[i];
        if( dataOutDouble[i] < dTargetMin ) { dTargetMin = dataOutDouble[i]; }
        if( dataOutDouble[i] > dTargetMax ) { dTargetMax = dataOutDouble[i]; }
    }

    return ( dTargetMass - dSourceMass );
}

///////////////////////////////////////////////////////////////////////////////////////////

