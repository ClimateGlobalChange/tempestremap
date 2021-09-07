///////////////////////////////////////////////////////////////////////////////
///
/// \file	GenerateOverlapMeshLint.cpp
/// \author  Paul Ullrich
/// \version August 19, 2021
///
/// <remarks>
///	 Copyright 2021 Paul Ullrich
///
///	 This file is distributed as part of the Tempest source code package.
///	 Permission is granted to use, copy, modify and distribute this
///	 source code and its documentation under the terms of the GNU General
///	 Public License.  This software is provided "as is" without express
///	 or implied warranty.
/// </remarks>

#include "Announce.h"
#include "CommandLine.h"
#include "STLStringHelper.h"
#include "Exception.h"
#include "GridElements.h"
#include "OverlapMesh.h"

#include "netcdfcpp.h"
#include "NetCDFUtilities.h"

#include <cmath>

#if defined(TEMPEST_MPIOMP)
#include <mpi.h>
#endif

///////////////////////////////////////////////////////////////////////////////

extern "C"
int GenerateOverlapWithMeshesLint (
	Mesh & meshA,
	Mesh & meshB,
	Mesh & meshOverlap,
	std::string strOverlapMesh,
	std::string strOutputFormat,
	std::string strMethod,
	const bool fNoValidate,
	const bool fHasConcaveFacesA,
	const bool fHasConcaveFacesB,
	const bool fParallel,
	std::string strTempDir,
	const bool fVerbose
) {

	NcError error ( NcError::silent_nonfatal );

	try
	{
		// Check command line parameters (data type arguments)
		STLStringHelper::ToLower(strOutputFormat);

		NcFile::FileFormat eOutputFormat =
			GetNcFileFormatFromString(strOutputFormat);
		if (eOutputFormat == NcFile::BadFormat) {
			_EXCEPTION1("Invalid \"out_format\" value (%s), "
				"expected [Classic|Offset64Bits|Netcdf4|Netcdf4Classic]",
				strOutputFormat.c_str());
		}

		// Method string
		OverlapMeshMethod method;
		STLStringHelper::ToLower ( strMethod );

		if (strMethod == "fuzzy") {
			method = OverlapMeshMethod_Fuzzy;

		} else if (strMethod == "exact") {
			method = OverlapMeshMethod_Exact;

		} else if (strMethod == "mixed") {
			method = OverlapMeshMethod_Mixed;

		} else {
			_EXCEPTIONT("Invalid \"method\" value, expected \"fuzzy\", \"exact\" or \"mixed\"");
		}

		// Convexify Mesh A
		if (fHasConcaveFacesA) {
			Mesh meshTemp = meshA;
			ConvexifyMesh(meshTemp, meshA, fVerbose);
		}

		// Validate Mesh A
		if (!fNoValidate) {
			AnnounceStartBlock("Validate mesh A");
			meshA.Validate();
			AnnounceEndBlock(NULL);
		}

		// Convexify Mesh B
		if (fHasConcaveFacesB) {
			Mesh meshTemp = meshB;
			ConvexifyMesh(meshTemp, meshB, fVerbose);
		}

		// Validate Mesh B
		if (!fNoValidate) {
			AnnounceStartBlock("Validate mesh B");
			meshB.Validate();
			AnnounceEndBlock(NULL);
		}

		// Generate meshes
		meshOverlap.type = Mesh::MeshType_Overlap;

		AnnounceStartBlock ( "Construct overlap mesh" );
		GenerateOverlapMeshLint (
			meshA, meshB,
			meshOverlap,
			method,
			fParallel,
			strTempDir,
			fVerbose );
		AnnounceEndBlock ( NULL );

		/*
			GenerateOverlapMeshFromFace(
				meshA,
				meshB,
				0,
				meshOverlap,
				method);

			// Construct the reverse node array on both meshes
			AnnounceStartBlock("Constructing reverse node array on input mesh");
			meshA.ConstructReverseNodeArray();
			AnnounceEndBlock(NULL);

			AnnounceStartBlock("Constructing reverse node array on output mesh");
			meshB.ConstructReverseNodeArray();
			AnnounceEndBlock(NULL);

			// Equalize nearly coincident nodes on these Meshes
			AnnounceStartBlock("Equalize coicident Nodes");
			EqualizeCoincidentNodes(meshA, meshB);
			AnnounceEndBlock(NULL);

			// Construct the overlap mesh
			Mesh meshOverlap;

			AnnounceStartBlock("Construct overlap mesh");
			GenerateOverlapMesh_v1(meshA, meshB, meshOverlap, method);
			AnnounceEndBlock(NULL);
		*/

		// Write the overlap mesh
		if ( strOverlapMesh.size() )
		{
			AnnounceStartBlock("Writing overlap mesh");
			meshOverlap.Write(strOverlapMesh.c_str(), eOutputFormat);
			AnnounceEndBlock(NULL);
		}

		return 0;

	}
	catch ( Exception& e )
	{
		Announce ( e.ToString().c_str() );
		return ( 0 );

	}
	catch ( ... )
	{
		return ( 0 );
	}
}

///////////////////////////////////////////////////////////////////////////////

extern "C"
int GenerateOverlapMeshLint(
	std::string strMeshA,
	std::string strMeshB,
	Mesh & meshOverlap,
	std::string strOverlapMesh,
	std::string strOutputFormat,
	std::string strMethod,
	const bool fNoValidate,
	const bool fHasConcaveFacesA,
	const bool fHasConcaveFacesB,
	const bool fParallel,
	std::string strTempDir,
	const bool fVerbose
) {

	NcError error ( NcError::silent_nonfatal );

	try
	{
		// Set output format
		STLStringHelper::ToLower(strOutputFormat);

		NcFile::FileFormat eOutputFormat =
			GetNcFileFormatFromString(strOutputFormat);
		if (eOutputFormat == NcFile::BadFormat) {
			_EXCEPTION1("Invalid \"out_format\" value (%s), "
				"expected [Classic|Offset64Bits|Netcdf4|Netcdf4Classic]",
				strOutputFormat.c_str());
		}

		// Prepare for parallel computation
		int nMPIRank = 0;
		int nMPISize = 1;

#if defined(TEMPEST_MPIOMP)
		if (fParallel) {
			MPI_Comm_rank(MPI_COMM_WORLD, &nMPIRank);
			MPI_Comm_size(MPI_COMM_WORLD, &nMPISize);
		}
#endif

		// Input meshes
		Mesh meshA;
		Mesh meshB;

		// Run in serial
		if (nMPISize == 1) {
			AnnounceStartBlock("Loading mesh A");
			meshA.Read(strMeshA);
			meshA.RemoveZeroEdges();
			AnnounceEndBlock(NULL);

			AnnounceStartBlock("Loading mesh B");
			meshB.Read(strMeshB);
			meshB.RemoveZeroEdges();
			AnnounceEndBlock(NULL);

			int err =
				GenerateOverlapWithMeshesLint (
					meshA, meshB,
					meshOverlap,
					strOverlapMesh,
					strOutputFormat,
					strMethod,
					fNoValidate,
					fHasConcaveFacesA,
					fHasConcaveFacesB,
					fParallel, strTempDir,
					fVerbose);

			return err;

		// Run in parallel
		} else {

#if !defined(TEMPEST_MPIOMP)
			_EXCEPTIONT("MPI not enabled: Cannot run in parallel");
#else
			// Number of pieces (used in parallel operation)
			int nPiecesA = sqrt(nMPISize);
			int nPiecesB = nMPISize / nPiecesA;
			int nPiecesTotal = nPiecesA * nPiecesB;

			// Split mesh in serial (on MPI Rank 0)
			if (nMPIRank == 0) {
				AnnounceStartBlock("Loading mesh A");
				meshA.Read(strMeshA);
				meshA.RemoveZeroEdges();
				AnnounceEndBlock(NULL);

				Announce("Splitting mesh A into %i pieces", nPiecesA);
				meshA.SplitAndWrite(strTempDir + std::string("/meshA"), nPiecesA);

				AnnounceStartBlock("Loading mesh B");
				meshB.Read(strMeshB);
				meshB.RemoveZeroEdges();
				AnnounceEndBlock(NULL);

				AnnounceStartBlock("Splitting mesh B into %i pieces", nPiecesB);
				meshB.SplitAndWrite(strTempDir + std::string("/meshB"), nPiecesB);
				AnnounceEndBlock(NULL);
			}

			// Synchronize tasks
			if (fParallel) {
				MPI_Barrier(MPI_COMM_WORLD);
			}

			// Computer overlap meshes in parallel
			AnnounceStartBlock("Computing overlap meshes (in parallel)");

			int iPieceAB = 0;
			for (int iPieceA = 0; iPieceA < nPiecesA; iPieceA++) {
			for (int iPieceB = 0; iPieceB < nPiecesB; iPieceB++) {

				if (iPieceAB % nMPISize == nMPIRank) {

					char szIndexA[10];
					snprintf(szIndexA, 10, "%06i", iPieceA);
					char szIndexB[10];
					snprintf(szIndexB, 10, "%06i", iPieceB);

					std::string strMeshAFile = strTempDir + std::string("/meshA") + szIndexA + std::string(".g");
					std::string strMeshBFile = strTempDir + std::string("/meshB") + szIndexB + std::string(".g");
					std::string strMeshOutFile = strTempDir + std::string("/meshOvA") + szIndexA + std::string("B") + szIndexB + std::string(".g");
					std::string strMeshLogFile = strTempDir + std::string("/logA") + szIndexA + std::string("B") + szIndexB + std::string(".txt");

					FILE * fpLog = fopen(strMeshLogFile.c_str(), "w");
					if (fpLog == NULL) {
						_EXCEPTION1("Unable to open log file \"%s\" for writing", strMeshLogFile.c_str());
					}
					AnnounceSetOutputBuffer(fpLog);
					AnnounceOutputOnAllRanks();

					Announce("Processor %i computing overlap on mesh pieces A%i and B%i\n", nMPIRank, iPieceA, iPieceB);

					Mesh meshA(strMeshAFile);
					Mesh meshB(strMeshBFile);
					Mesh meshOverlapSub;

					int err =
						GenerateOverlapWithMeshesLint (
							meshA, meshB,
							meshOverlapSub,
							strMeshOutFile,
							strOutputFormat,
							strMethod,
							fNoValidate,
							fHasConcaveFacesA,
							fHasConcaveFacesB,
							fParallel, strTempDir,
							fVerbose);

					AnnounceOnlyOutputOnRankZero();
					AnnounceSetOutputBuffer(stdout);

					if (err) {
						MPI_Abort(MPI_COMM_WORLD, err);
						return err;
					}
				}

				iPieceAB++;
			}
			}

			AnnounceEndBlock("Done");

			// Synchronize tasks
			if (fParallel) {
				MPI_Barrier(MPI_COMM_WORLD);
			}

			// Merge overlap meshes in serial (on MPI Rank 0)
			if (nMPIRank == 0) {
				AnnounceStartBlock("Merging meshes");

				Mesh meshCombined;
				meshCombined.Read(strTempDir + std::string("/meshOvA000000B000000.g"));
				meshCombined.BeginAppend();

				int iPieceAB = 0;
				for (int iPieceA = 0; iPieceA < nPiecesA; iPieceA++) {
				for (int iPieceB = 0; iPieceB < nPiecesB; iPieceB++) {
					if ((iPieceA == 0) && (iPieceB == 0)) {
						continue;
					}

					char szIndexA[10];
					snprintf(szIndexA, 10, "%06i", iPieceA);
					char szIndexB[10];
					snprintf(szIndexB, 10, "%06i", iPieceB);
					std::string strMeshOvFile = strTempDir + std::string("/meshOvA") + szIndexA + std::string("B") + szIndexB + std::string(".g");

					Mesh meshOther;
					meshOther.Read(strMeshOvFile);
					meshCombined.Append(meshOther);
				}
				}

				meshCombined.EndAppend();

				AnnounceEndBlock("Done");

				AnnounceStartBlock("Writing final overlap mesh");
				meshCombined.Write(strOverlapMesh);
				AnnounceEndBlock("Done");
			}

			// Clean-up files
			if (nMPIRank == 0) {
				AnnounceStartBlock("Cleanup temporary files");

				for (int iPieceA = 0; iPieceA < nPiecesA; iPieceA++) {
					char szIndexA[10];
					snprintf(szIndexA, 10, "%06i", iPieceA);
					std::string strMeshAFile = strTempDir + std::string("/meshA") + szIndexA + std::string(".g");
					Announce("\"%s\"", strMeshAFile.c_str());
					remove(strMeshAFile.c_str());
				}
				for (int iPieceB = 0; iPieceB < nPiecesB; iPieceB++) {
					char szIndexB[10];
					snprintf(szIndexB, 10, "%06i", iPieceB);
					std::string strMeshBFile = strTempDir + std::string("/meshB") + szIndexB + std::string(".g");
					Announce("\"%s\"", strMeshBFile.c_str());
					remove(strMeshBFile.c_str());
				}

				for (int iPieceA = 0; iPieceA < nPiecesA; iPieceA++) {
				for (int iPieceB = 0; iPieceB < nPiecesB; iPieceB++) {

					char szIndexA[10];
					snprintf(szIndexA, 10, "%06i", iPieceA);
					char szIndexB[10];
					snprintf(szIndexB, 10, "%06i", iPieceB);
					std::string strMeshOvFile = strTempDir + std::string("/meshOvA") + szIndexA + std::string("B") + szIndexB + std::string(".g");
					std::string strMeshLogFile = strTempDir + std::string("/logA") + szIndexA + std::string("B") + szIndexB + std::string(".txt");

					Announce("\"%s\"", strMeshOvFile.c_str());
					remove(strMeshOvFile.c_str());
					//Announce("\"%s\"", strMeshLogFile.c_str());
					//remove(strMeshLogFile.c_str());
				}
				}

				AnnounceEndBlock("Done");
			}
#endif
		}

		return 0;
	}
	catch ( Exception& e )
	{
		Announce ( e.ToString().c_str() );
#if defined(TEMPEST_MPIOMP)
		int flag;
		MPI_Initialized(&flag);
		if (flag) {
			MPI_Abort(MPI_COMM_WORLD, -1);
		}
#endif
		return ( 0 );

	}
	catch ( ... )
	{
		return ( 0 );
	}
}

///////////////////////////////////////////////////////////////////////////////
