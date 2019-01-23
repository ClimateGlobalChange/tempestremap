///////////////////////////////////////////////////////////////////////////////
///
///	\file    CalculateDiffNorms.cpp
///	\author  Paul Ullrich
///	\version August 31st, 2014
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

#include "CommandLine.h"
#include "GridElements.h"
#include "FiniteElementTools.h"
#include "GaussLobattoQuadrature.h"
#include "Exception.h"
#include "Announce.h"
#include "DataArray1D.h"
#include "DataArray3D.h"

#include "netcdfcpp.h"

#include <cmath>
#include <iostream>

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	NcError error(NcError::silent_nonfatal);

try {
	// First data file
	std::string strFileA;

	// Second data file
	std::string strFileB;

	// Variable to compare
	std::string strVariableName;

	// Mesh file to use
	std::string strMeshFile;

	// Input data is on GLL nodes
	bool fGLL;

	// Degree of polynomial
	int nP;

	// Use bubble function
	bool fBubble;

	// Output filename
	std::string strOutputFile;

	// Contains concave faces
	bool fContainsConcaveFaces;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strFileA, "a", "");
		CommandLineString(strFileB, "b", "");
		CommandLineString(strVariableName, "var", "Psi");
		CommandLineBool(fGLL, "gll");
		CommandLineInt(nP, "np", 4);
		CommandLineBool(fBubble, "bubble");
		CommandLineString(strMeshFile, "mesh", "");
		CommandLineString(strOutputFile, "outfile", "");
		CommandLineBool(fContainsConcaveFaces, "concave");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	// Input mesh
	AnnounceStartBlock("Loading Mesh");
	Mesh mesh(strMeshFile);

	// Check for rectilinear Mesh
	NcFile ncMesh(strMeshFile.c_str(), NcFile::ReadOnly);

	// Check for rectilinear attribute
	bool fRectilinear = false;
	for (int a = 0; a < ncMesh.num_atts(); a++) {
		if (strcmp(ncMesh.get_att(a)->name(), "rectilinear") == 0) {
			fRectilinear = true;
			break;
		}
	}

	if (fGLL && fRectilinear) {
		_EXCEPTIONT("--gll cannot be used with rectilinear grids");
	}

	std::vector<long> vecOutputDimSizes;
	std::vector<std::string> vecOutputDimNames;

	// Non-rectilinear output
	if (!fRectilinear) {
		vecOutputDimSizes.resize(1);
		vecOutputDimSizes[0] = mesh.faces.size();

		vecOutputDimNames.resize(1);
		vecOutputDimNames[0] = "ncol";

		Announce("Non-rectilinear mesh detected");

	// Obtain rectilinear attributes
	} else {
		int nDim0Size = ncMesh.get_att("rectilinear_dim0_size")->as_int(0);
		int nDim1Size = ncMesh.get_att("rectilinear_dim1_size")->as_int(0);

		std::string strDim0Name =
			ncMesh.get_att("rectilinear_dim0_name")->as_string(0);
		std::string strDim1Name =
			ncMesh.get_att("rectilinear_dim1_name")->as_string(0);

		vecOutputDimSizes.resize(2);
		vecOutputDimSizes[0] = nDim0Size;
		vecOutputDimSizes[1] = nDim1Size;

		vecOutputDimNames.resize(2);
		vecOutputDimNames[0] = strDim0Name;
		vecOutputDimNames[1] = strDim1Name;

		Announce("Rectilinear mesh detected");
	}

	AnnounceEndBlock("Done");

	// Get Mesh information
	int nTotalDataSize = 1;

	DataArray1D<double> dataUniqueJacobian;
	
	// Finite volumes
	if (!fGLL) {
		// Total data size
		for (int d = 0; d < vecOutputDimSizes.size(); d++) {
			nTotalDataSize *= vecOutputDimSizes[d];
		}

	// Finite elements
	} else {
		// Calculate pointwise areas
		DataArray3D<int> dataGLLnodes;
		DataArray3D<double> dataGLLJacobian;

		GenerateMetaData(mesh, nP, fBubble, dataGLLnodes, dataGLLJacobian);

		GenerateUniqueJacobian(
			dataGLLnodes, dataGLLJacobian, dataUniqueJacobian);

		for (int s = 0; s < nP; s++) {
		for (int t = 0; t < nP; t++) {
		for (int i = 0; i < dataGLLnodes.GetSubColumns(); i++) {
			if (dataGLLnodes[s][t][i] > vecOutputDimSizes[0]) {
				vecOutputDimSizes[0] = dataGLLnodes[s][t][i];
			}
		}
		}
		}

		nTotalDataSize = vecOutputDimSizes[0];
	}

	// Load data from file A
	AnnounceStartBlock("Loading data from file A");

	DataArray1D<double> dDataA(nTotalDataSize);

	NcFile ncFileA(strFileA.c_str(), NcFile::ReadOnly);
	if (!ncFileA.is_valid()) {
		_EXCEPTION1("Unable to open file \"%s\"", strFileA.c_str());
	}

	NcVar * varA = ncFileA.get_var(strVariableName.c_str());
	if (varA == NULL) {
		_EXCEPTION2("File \"%s\" does not contain variable \"%s\"",
			strFileA.c_str(), strVariableName.c_str());
	}

	varA->get(&(dDataA[0]), &(vecOutputDimSizes[0]));

	int nVarASize = varA->get_dim(0)->size();

	ncFileA.close();

	AnnounceEndBlock("Done");

	// Load data from file B
	AnnounceStartBlock("Loading data from file B");

	DataArray1D<double> dDataB(nTotalDataSize);

	NcFile ncFileB(strFileB.c_str(), NcFile::ReadOnly);
	if (!ncFileB.is_valid()) {
		_EXCEPTION1("Unable to open file \"%s\"", strFileB.c_str());
	}

	NcVar * varB = ncFileB.get_var(strVariableName.c_str());
	if (varB == NULL) {
		_EXCEPTION2("File \"%s\" does not contain variable \"%s\"",
			strFileB.c_str(), strVariableName.c_str());
	}

	varB->get(&(dDataB[0]), &(vecOutputDimSizes[0]));

	int nVarBSize = varB->get_dim(0)->size();

	ncFileB.close();

	AnnounceEndBlock("Done");

	// Check sizes
	if (nVarASize != nVarBSize) {
		_EXCEPTION2("Variable size mismatch [%i,%i]", nVarASize, nVarBSize);
	}

	// Calculate the difference and store in dDataA
	double dMaxA = dDataA[0];
	double dMinA = dDataA[0];

	double dMaxB = dDataB[0];
	double dMinB = dDataB[0];

	for (int i = 0; i < nTotalDataSize; i++) {
		if (dDataA[i] > dMaxA) {
			dMaxA = dDataA[i];
		}
		if (dDataA[i] < dMinA) {
			dMinA = dDataA[i];
		}

		if (dMinA == 0.0) {
			_EXCEPTIONT("Zero minimum field value");
		}

		if (dDataB[i] > dMaxB) {
			dMaxB = dDataB[i];
		}
		if (dDataB[i] < dMinB) {
			dMinB = dDataB[i];
		}

		dDataA[i] = fabs(dDataA[i] - dDataB[i]);
	}

	// Min / Max Norm
	double dNormLmin;
	double dNormLmax;

	if (dMaxB == dMinB) {
		dNormLmin = dMinB - dMinA;
		dNormLmax = dMaxA - dMaxB;
	} else {
		dNormLmin = (dMinB - dMinA) / (dMaxB - dMinB);
		dNormLmax = (dMaxA - dMaxB) / (dMaxB - dMinB);
	}

	// Norms and Sums
	double dNormL1 = 0.0;
	double dNormL2 = 0.0;
	double dNormLi = 0.0;

	double dSumL1 = 0.0;
	double dSumL2 = 0.0;
	double dSumLi = 0.0;

	// Use face areas
	if (!fGLL) {
		mesh.CalculateFaceAreas(fContainsConcaveFaces);

		for (int i = 0; i < nTotalDataSize; i++) {
			dNormL1 += dDataA[i] * mesh.vecFaceArea[i];
			dNormL2 += dDataA[i] * dDataA[i] * mesh.vecFaceArea[i];

			if (dDataA[i] > dNormLi) {
				dNormLi = dDataA[i];
			}

			dSumL1 += dDataB[i] * mesh.vecFaceArea[i];
			dSumL2 += dDataB[i] * dDataB[i] * mesh.vecFaceArea[i];

			if (dDataB[i] > dSumLi) {
				dSumLi = dDataB[i];
			}
		}

		dNormL1 = dNormL1 / dSumL1;
		dNormL2 = sqrt(dNormL2 / dSumL2);
		dNormLi = dNormLi / dSumLi;

	} else {

		// Calculate differences
		for (int i = 0; i < nTotalDataSize; i++) {
			dNormL1 += dDataA[i] * dataUniqueJacobian[i];
			dNormL2 += dDataA[i] * dDataA[i] * dataUniqueJacobian[i];

			if (dDataA[i] > dNormLi) {
				dNormLi = dDataA[i];
			}

			dSumL1 += dDataB[i] * dataUniqueJacobian[i];
			dSumL2 += dDataB[i] * dDataB[i] * dataUniqueJacobian[i];

			if (dDataB[i] > dSumLi) {
				dSumLi = dDataB[i];
			}
		}

		dNormL1 = dNormL1 / dSumL1;
		dNormL2 = sqrt(dNormL2 / dSumL2);
		dNormLi = dNormLi / dSumLi;
	}

	// Announce results
	AnnounceStartBlock("Results:");
	Announce("L1:   %1.15e | %1.15e", dNormL1, dSumL1);
	Announce("L2:   %1.15e | %1.15e", dNormL2, dSumL2);
	Announce("Li:   %1.15e | %1.15e", dNormLi, dSumLi);
	Announce("Lmin: %1.15e | %1.5e %1.5e", dNormLmin, dMinA, dMinB);
	Announce("Lmax: %1.15e | %1.5e %1.5e", dNormLmax, dMaxA, dMaxB);
	AnnounceEndBlock(NULL);

	// Print results to file
	if (strOutputFile != "") {
		FILE * fp = fopen(strOutputFile.c_str(), "a");
		fprintf(fp, "%1.15e %1.15e %1.15e %1.15e %1.15e\n",
			dNormL1, dNormL2, dNormLi, dNormLmin, dNormLmax);
		fclose(fp);
	}

	return (0);

} catch(Exception & e) {
	Announce(e.ToString().c_str());
	return (-1);

} catch(...) {
	return (-2);
}
}


///////////////////////////////////////////////////////////////////////////////

