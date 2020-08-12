///////////////////////////////////////////////////////////////////////////////
///
///	\file    GenerateStereographicMesh.cpp
///	\author  Paul Ullrich
///	\version July 20, 2020
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

#include "GridElements.h"
#include "Exception.h"
#include "Announce.h"
#include "STLStringHelper.h"
#include "NetCDFUtilities.h"

#include "netcdfcpp.h"

#include <cmath>
#include <iostream>

///////////////////////////////////////////////////////////////////////////////
// 
// Input Parameters:
// Output filename:  std::string strOutputFile;
// 
// Output Parameters: Mesh*
// 
extern "C" 
int GenerateStereographicMesh(
	Mesh & mesh,
	double dLonDegP,
	double dLatDegP,
	double dLonDeg0,
	double dLatDeg0,
	double dLonDeg1,
	double dLatDeg1,
	int nXElements,
	int nYElements,
	bool fCentroids,
	std::string strOutputFile,
	std::string strOutputFormat
) {

	NcError error(NcError::silent_nonfatal);

try {

	// Announce
	std::cout << "=========================================================";
	std::cout << std::endl;

	// Check arguments
	if (fCentroids) {
		if (nXElements < 2) {
			_EXCEPTIONT("At least two X grid elements expected");
		}
		if (nYElements < 2) {
			_EXCEPTIONT("At least two Y grid elements expected");
		}

	} else {
		if (nXElements < 1) {
			_EXCEPTIONT("At least one X grid element expected");
		}
		if (nYElements < 1) {
			_EXCEPTIONT("At least one Y grid element expected");
		}
	}

    // Check command line parameters (data type arguments)
    STLStringHelper::ToLower(strOutputFormat);

	NcFile::FileFormat eOutputFormat =
		GetNcFileFormatFromString(strOutputFormat);
	if (eOutputFormat == NcFile::BadFormat) {
		_EXCEPTION1("Invalid \"out_format\" value (%s), "
			"expected [Classic|Offset64Bits|Netcdf4|Netcdf4Classic]",
			strOutputFormat.c_str());
	}

	// Announce
	std::cout << "..Generating polar stereographic mesh" << std::endl;

	NodeVector & nodes = mesh.nodes;
	FaceVector & faces = mesh.faces;
    mesh.type = Mesh::MeshType_Transect;

	// Convert degrees to radians
	double dLonRadP = dLonDegP * M_PI / 180.0;
	double dLatRadP = dLatDegP * M_PI / 180.0;

	double dLonRad0 = dLonDeg0 * M_PI / 180.0;
	double dLatRad0 = dLatDeg0 * M_PI / 180.0;

	double dLonRad1 = dLonDeg1 * M_PI / 180.0;
	double dLatRad1 = dLatDeg1 * M_PI / 180.0;

	// Calculate corner points in stereographic projection
	double dA0;
	double dB0;

	StereographicProjection(
		dLonRadP,
		dLatRadP,
		dLonRad0,
		dLatRad0,
		dA0,
		dB0);

	double dA1;
	double dB1;

	StereographicProjection(
		dLonRadP,
		dLatRadP,
		dLonRad1,
		dLatRad1,
		dA1,
		dB1);

	// Calculate other corner points
	double dLonRad2;
	double dLatRad2;

	double dLonRad3;
	double dLatRad3;

	StereographicProjectionInv(
		dLonRadP,
		dLatRadP,
		dA1,
		dB0,
		dLonRad2,
		dLatRad2);

	StereographicProjectionInv(
		dLonRadP,
		dLatRadP,
		dA0,
		dB1,
		dLonRad3,
		dLatRad3);

	printf("Corners:\n");
	printf("  [%2.14f, %2.14f]\n", RadToDeg(dLonRad0), RadToDeg(dLatRad0));
	printf("  [%2.14f, %2.14f]\n", RadToDeg(dLonRad2), RadToDeg(dLatRad2));
	printf("  [%2.14f, %2.14f]\n", RadToDeg(dLonRad1), RadToDeg(dLatRad1));
	printf("  [%2.14f, %2.14f]\n", RadToDeg(dLonRad3), RadToDeg(dLatRad3));

	// Calculate grid spacing in stereographic projection
	double dDeltaA;
	double dDeltaB;

	if (fCentroids) {
		dDeltaA = (dA1 - dA0) / static_cast<double>(nXElements-1);
		dDeltaB = (dB1 - dB0) / static_cast<double>(nYElements-1);

		dA0 -= 0.5 * dDeltaA;
		dB0 -= 0.5 * dDeltaB;

	} else {
		dDeltaA = (dA1 - dA0) / static_cast<double>(nXElements);
		dDeltaB = (dB1 - dB0) / static_cast<double>(nYElements);
	}

	// Insert vertices into mesh
	for (int j = 0; j <= nYElements; j++) {
		double dB = dB0 + dDeltaB * static_cast<double>(j);

		for (int i = 0; i <= nXElements; i++) {
			double dA = dA0 + dDeltaA * i;

			double dLonRad;
			double dLatRad;
			StereographicProjectionInv(
				dLonRadP,
				dLatRadP,
				dA,
				dB,
				dLonRad,
				dLatRad);

			double dX = cos(dLonRad) * cos(dLatRad);
			double dY = sin(dLonRad) * cos(dLatRad);
			double dZ = sin(dLatRad);

			_ASSERT(fabs(dX * dX + dY * dY + dZ * dZ - 1.0) < ReferenceTolerance);

			nodes.push_back(Node(dX,dY,dZ));
		}
	}

	std::cout << "..Inserting faces" << std::endl;

	// Insert faces
	for (int j = 0; j < nYElements; j++) {
	for (int i = 0; i < nXElements; i++) {
		Face face(4);
		face.SetNode(0,  i      +  j      * (nXElements+1));
		face.SetNode(1, (i + 1) +  j      * (nXElements+1));
		face.SetNode(2, (i + 1) + (j + 1) * (nXElements+1));
		face.SetNode(3,  i      + (j + 1) * (nXElements+1));
		faces.push_back(face);
	}
	}

	// Output the mesh
	if (strOutputFile.size()) {

		// Announce
		std::cout << "..Writing mesh to file [" << strOutputFile.c_str() << "] ";
		std::cout << std::endl;

		mesh.Write(strOutputFile, eOutputFormat);

		NcFile ncOutput(strOutputFile.c_str(), NcFile::Write);
		ncOutput.add_att("rectilinear", "true");
		ncOutput.add_att("rectilinear_dim0_size", nYElements);
		ncOutput.add_att("rectilinear_dim1_size", nXElements);
		ncOutput.add_att("rectilinear_dim0_name", "Y");
		ncOutput.add_att("rectilinear_dim1_name", "X");
	}

	// Announce
	std::cout << "..Mesh generator exited successfully" << std::endl;
	std::cout << "=========================================================";
	std::cout << std::endl;

} catch(Exception & e) {
	Announce(e.ToString().c_str());
	return (0);

} catch(...) {
	return (0);
}
    return 0;
}

///////////////////////////////////////////////////////////////////////////////
