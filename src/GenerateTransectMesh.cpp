///////////////////////////////////////////////////////////////////////////////
///
///	\file    GenerateTransectMesh.cpp
///	\author  Paul Ullrich
///	\version February 13, 2020
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
// Number of elements in mesh: int nResolution;
// Alternate arrangement: bool fAlt;
// Output filename:  std::string strOutputFile;
// 
// Output Parameters: Mesh*
// 
extern "C" 
int GenerateTransectMesh(
	Mesh & mesh,
	double dLonDeg0,
	double dLatDeg0,
	double dLonDeg1,
	double dLatDeg1,
	double dPerpDtheta,
	int nParaElements,
	int nPerpElements,
	std::string strOutputFile,
	std::string strOutputFormat
) {

	NcError error(NcError::silent_nonfatal);

try {

	// Check arguments
	if (nParaElements < 1) {
		_EXCEPTIONT("At least one grid element expected along transect");
	}
	if (nPerpElements < 1) {
		_EXCEPTIONT("At least one grid element expected perpendicular to transect");
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
	std::cout << "=========================================================";
	std::cout << std::endl;
	std::cout << "..Generating transect mesh" << std::endl;

	NodeVector & nodes = mesh.nodes;
	FaceVector & faces = mesh.faces;
    mesh.type = Mesh::MeshType_Transect;

	// Calculate angle between initial and final points
	double dLonRad0 = dLonDeg0 * M_PI / 180.0;
	double dLatRad0 = dLatDeg0 * M_PI / 180.0;

	double dLonRad1 = dLonDeg1 * M_PI / 180.0;
	double dLatRad1 = dLatDeg1 * M_PI / 180.0;

	double dX0 = cos(dLonRad0) * cos(dLatRad0);
	double dY0 = sin(dLonRad0) * cos(dLatRad0);
	double dZ0 = sin(dLatRad0);

	double dX1 = cos(dLonRad1) * cos(dLatRad1);
	double dY1 = sin(dLonRad1) * cos(dLatRad1);
	double dZ1 = sin(dLatRad1);

	double dDot0 = dX0 * dX0 + dY0 * dY0 + dZ0 * dZ0;
	double dDot1 = dX1 * dX1 + dY1 * dY1 + dZ1 * dZ1;

	if (fabs(dDot0 - 1.0) > 1.0e-12) {
		_EXCEPTIONT("Logic error");
	}
	if (fabs(dDot1 - 1.0) > 1.0e-12) {
		_EXCEPTIONT("Logic error");
	}

	// Dot product
	double dDot = dX0 * dX1 + dY0 * dY1 + dZ0 * dZ1;
	if (fabs(dDot) >= 1.0) {
		_EXCEPTIONT("Transect uses coincident endpoints");
	}

	double dParaDtheta = acos(dDot) / static_cast<double>(nParaElements);

	if (dPerpDtheta <= 0.0) {
		dPerpDtheta = dParaDtheta;
	}

	printf("..Parallel resolution: %1.5f degrees (GCD)\n", dParaDtheta * 180.0 / M_PI);
	printf("..Perpendicular resolution: %1.5f degrees (GCD)\n", dPerpDtheta * 180.0 / M_PI);

	// Vector pointing from p0 to p1
	double dXt = dX1 - dX0;
	double dYt = dY1 - dY0;
	double dZt = dZ1 - dZ0;

	double dMagt = sqrt(dXt * dXt + dYt * dYt + dZt * dZt);

	// Perpendicular vector
	double dXp = dY1 * dZ0 - dZ1 * dY0;
	double dYp = dZ1 * dX0 - dX1 * dZ0;
	double dZp = dX1 * dY0 - dY1 * dX0;

	double dMagp = sqrt(dXp * dXp + dYp * dYp + dZp * dZp);

	std::cout << "..Inserting vertices" << std::endl;

	// Perpendicular angle start
	double dPerpTheta0 = - 0.5 * static_cast<double>(nPerpElements) * dPerpDtheta;

	// Insert vertices of transect
	for (int j = 0; j <= nParaElements; j++) {
		//double dS = (cos(dParaDtheta * static_cast<double>(j)) - 1.0) / (dDot - 1.0);
		double dDeltaTheta = dParaDtheta * static_cast<double>(j);
		double dS = sin(dDeltaTheta) / (cos(dDeltaTheta - 0.5 * acos(dDot)) * dMagt);

		double dXn = dX0 + dXt * dS;
		double dYn = dY0 + dYt * dS;
		double dZn = dZ0 + dZt * dS;

		double dMag = sqrt(dXn * dXn + dYn * dYn + dZn * dZn);

		dXn /= dMag;
		dYn /= dMag;
		dZn /= dMag;

		for (int i = 0; i <= nPerpElements; i++) {
			double dTanTheta =
				tan(dPerpTheta0 + static_cast<double>(i) * dPerpDtheta);

			double dX = dXn + dXp * dTanTheta / dMagp;
			double dY = dYn + dYp * dTanTheta / dMagp;
			double dZ = dZn + dZp * dTanTheta / dMagp;
			dMag = sqrt(dX * dX + dY * dY + dZ * dZ);

			dX /= dMag;
			dY /= dMag;
			dZ /= dMag;

			nodes.push_back(Node(dX,dY,dZ));
		}
	}

	std::cout << "..Inserting faces" << std::endl;


	// Insert faces
	for (int j = 0; j < nParaElements; j++) {
	for (int i = 0; i < nPerpElements; i++) {
		Face face(4);
		face.SetNode(0,  i      +  j      * (nPerpElements+1));
		face.SetNode(1, (i + 1) +  j      * (nPerpElements+1));
		face.SetNode(2, (i + 1) + (j + 1) * (nPerpElements+1));
		face.SetNode(3,  i      + (j + 1) * (nPerpElements+1));
		faces.push_back(face);
	}
	}

	// Output the mesh
	if (strOutputFile.size()) {

		// Announce
		std::cout << "..Writing mesh to file [" << strOutputFile.c_str() << "] ";
		std::cout << std::endl;

		mesh.Write(strOutputFile, eOutputFormat);
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
