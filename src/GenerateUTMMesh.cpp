///////////////////////////////////////////////////////////////////////////////
///
///	\file    GenerateUTMMesh.cpp
///	\author  Paul Ullrich
///	\version July 2, 2018
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

#include "GridElements.h"
#include "Exception.h"
#include "Announce.h"
#include "DataMatrix.h"

#include <cmath>
#include <iostream>
#include <complex>
#include <cfloat>

#include "netcdfcpp.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Convert a coordinate from UTM to RLL.
///	</summary>
///	<license>
///	Source code adapted from UTM2LL, available from:
///	https://www.mathworks.com/matlabcentral/fileexchange/45699-ll2utm-and-utm2ll
///	
///	UTM2LL license agreement as follows:
///	
///	Copyright (c) 2001-2015, Franois Beauducel, covered by BSD License.
///	All rights reserved.
///
///	Redistribution and use in source and binary forms, with or without 
///	modification, are permitted provided that the following conditions are 
///	met:
///
///	   * Redistributions of source code must retain the above copyright 
///	     notice, this list of conditions and the following disclaimer.
///	   * Redistributions in binary form must reproduce the above copyright 
///	     notice, this list of conditions and the following disclaimer in 
///	     the documentation and/or other materials provided with the distribution
///	                           
///	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
///	AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
///	IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
///	ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
///	LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
///	CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
///	SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
///	INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
///	CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
///	ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
///	POSSIBILITY OF SUCH DAMAGE.
///	</license>

///	<summary>
///		Calculate UTM to RLL projection coefficients
///	</summary>
///	<param name="dE">
///		First ellipsoid eccentricity.
///	</param>
///	<param name="iM">
///		iM = 0 for transverse mercator
///		iM = 1 for transverse mercator reverse coefficients
///		iM = 2 for merdian arc
///	</param>
///	<param name="dC">
///		Pointer to an array of length 5.
///	</param>
void ConvertUTMtoRLL_Coeff(
	double dE,
	int iM,
	double * dC
) {
	DataMatrix<double> dC0(5,9);

	dC0[0][0] = -175./16384.;
	dC0[0][2] = -5./256.;
	dC0[0][4] = -3./64.;
	dC0[0][6] = -1./4.;
	dC0[0][8] = 1.;

	// Transverse mercator
	if (iM == 0) {
		dC0[1][0] = -105./4096.;
		dC0[1][2] = -45./1024.;
		dC0[1][4] = -3./32.;
		dC0[1][6] = -3./8.;

		dC0[2][0] = 525./16384.;
		dC0[2][2] = 45./1024.;
		dC0[2][4] = 15./256.;
		
		dC0[3][0] = -175./12288.;
		dC0[3][2] = -35./3072.;
		
		dC0[4][0] = 315./131072.;

	// Transverse mercator reverse coefficients
	} else if (iM == 1) {
		dC0[1][0] = 1./61440.;
		dC0[1][2] = 7./2048.;
		dC0[1][4] = 1./48.;
		dC0[1][6] = 1./8.;

		dC0[2][0] = 559./368640.;
		dC0[2][2] = 3./1280.;
		dC0[2][4] = 1./768.;

		dC0[3][0] = 283./430080.;
		dC0[3][2] = 17./30720.;

		dC0[4][0] = 4397./41287680.;

	// Median arc
	} else if (iM == 2) {
		dC0[1][0] = -901./184320.;
		dC0[1][2] = -9./1024.;
		dC0[1][4] = -1./96.;
		dC0[1][6] = 1./8.;

		dC0[2][0] = -311./737280.;
		dC0[2][2] = 17./5120.;
		dC0[2][4] = 13./768.;

		dC0[3][0] = 899./430080.;
		dC0[3][2] = 61./15360.;

		dC0[4][0] = 49561./41287680.;

	} else {
		_EXCEPTION();
	}

	// Evaluate the polynomial
	for (int i = 0; i < 5; i++) {
		dC[i] = 0.0;
		double dEx = 1.0;
		for (int j = 0; j < 9; j++) {
			dC[i] += dC0[i][8-j] * dEx;
			dEx = dEx * dE;
		}
	}
}

///	<summary>
///		Convert a coordinate from UTM to RLL.
///	</summary>
void ConvertUTMtoRLL(
	int nZone,
	double dX,
	double dY,
	double & dLon,
	double & dLat
) {
	// Conversion from rad to degree
	static const double dD0 = 180.0 / M_PI;

	// Semi-major axis (in meters) from WGS84 standard
	static const double dA1 = 6378137.0;

	// Flattening of the ellipsoid from WGS84 standard
	static const double dF1 = 298.257223563;

	// Maximum iteration for latitude computation
	static const int nMaxIter = 100;

	// Minimum residue for latitude computation
	static const double dEps = 1.0e-11;

	// UTM scale factor
	static const double dK0 = 0.9996;

	// UTM false East (m)
	static const double dX0 = 500000.0;

	// UTM false North (m)
	const double dY0 = 1.0e7 * static_cast<double>(nZone < 0);

	// UTM origin latitude (rad)
	static const double dP0 = 0.0;

	// UTM origin longitude (rad)
	const double dL0 = (6.0 * fabs(static_cast<double>(nZone)) - 183.0) / dD0;

	// Ellipsoid eccentricity
	const double dE1x = dA1 * (1.0 - 1.0/dF1);
	const double dE1 = sqrt((dA1 * dA1 - dE1x * dE1x)/(dA1 * dA1));

	const double dN = dK0 * dA1;

	// Computing parameters for Mercator Transverse projection
	double dC[5];

	ConvertUTMtoRLL_Coeff(dE1, 0, dC);

	const double dYS =
		dY0 - dN * (
			  dC[0] * dP0
			+ dC[1] * sin(2.0 * dP0)
			+ dC[2] * sin(4.0 * dP0)
			+ dC[3] * sin(6.0 * dP0)
			+ dC[4] * sin(8.0 * dP0));

	ConvertUTMtoRLL_Coeff(dE1, 1, dC);

	std::complex<double> dZT((dY - dYS)/dN/dC[0], (dX - dX0)/dN/dC[0]);

	std::complex<double> dZ =
		dZT
		- dC[1] * sin(2.0 * dZT)
		- dC[2] * sin(4.0 * dZT)
		- dC[3] * sin(6.0 * dZT)
		- dC[4] * sin(8.0 * dZT);

	double dL = dZ.real();
	double dLS = dZ.imag();

	double dl = dL0 + atan(sinh(dLS) / cos(dL));
	double dp = asin(sin(dL) / cosh(dLS));

	dL = log(tan(M_PI/4.0 + dp/2.0));

	dp = 2.0 * atan(exp(dL)) - M_PI/2.0;
	
	int i = 0;
	double dp0 = DBL_MAX;
	while (((dp0 == DBL_MAX) || (fabs(dp - dp0) > dEps)) && (i < nMaxIter)) {
		dp0 = dp;
		double dES = dE1 * sin(dp0);
		dp = pow(2.0 * atan((1.0 + dES) / (1.0 - dES)), dE1/2.0) * exp(dL) - M_PI/2.0;
		i++;
	}
	if (i == nMaxIter) {
		_EXCEPTIONT("Convergence failure");
	}

	dLat = dp;
	dLon = dl;

/*
	% constants
D0 = 180/pi;	% conversion rad to deg
maxiter = 100;	% m
eps = 1e-11;	% m

K0 = 0.9996;					% 
X0 = 500000;					% 
Y0 = 1e7*(f < 0);				% UTM false North (m)
P0 = 0;						% 
L0 = (6*abs(f) - 183)/D0;			% 
E1 = sqrt((A1^2 - (A1*(1 - 1/F1))^2)/A1^2);	% ellpsoid excentricity
N = K0*A1;

% computing parameters for Mercator Transverse projection
C = coef(E1,0);
YS = Y0 - N*(C(1)*P0 + C(2)*sin(2*P0) + C(3)*sin(4*P0) + C(4)*sin(6*P0) + C(5)*sin(8*P0));

C = coef(E1,1);
zt = complex((y - YS)/N/C(1),(x - X0)/N/C(1));
z = zt - C(2)*sin(2*zt) - C(3)*sin(4*zt) - C(4)*sin(6*zt) - C(5)*sin(8*zt);
L = real(z);
LS = imag(z);

l = L0 + atan(sinh(LS)./cos(L));
p = asin(sin(L)./cosh(LS));

L = log(tan(pi/4 + p/2));

% calculates latitude from the isometric latitude
p = 2*atan(exp(L)) - pi/2;
p0 = NaN;
n = 0;
while any(isnan(p0) | abs(p - p0) > eps) && n < maxiter
	p0 = p;
	es = E1*sin(p0);
	p = 2*atan(((1 + es)./(1 - es)).^(E1/2).*exp(L)) - pi/2;
	n = n + 1;
end

if nargout < 2
	lat = D0*[p(:),l(:)];
else
	lat = p*D0;
	lon = l*D0;
end
*/
}

///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// 
// Input Parameters:
// UTM Zone:  int nZone;
// Number of columns in mesh: int nCols;
// Number of rows in mesh: int nRows;
// XLL Corner of the mesh: double dXLLCorner;
// YLL Corner of the mesh: double dYLLCorner;
// Cell size of the mesh: double dCellSize;
// Output filename: std::string strOutputFile;
// 
// Output Parameters: Mesh*
// 
extern "C" 
int GenerateUTMMesh(
	Mesh & mesh,
	int nZone,
	int nCols,
	int nRows,
	double dXLLCorner,
	double dYLLCorner,
	double dCellSize,
	std::string strOutputFile,
	bool fVerbose
) {

	NcError error(NcError::silent_nonfatal);

try {

	// Clear the mesh
	mesh.Clear();
    mesh.type = Mesh::MeshType_UTM;

	NodeVector & nodes = mesh.nodes;
	FaceVector & faces = mesh.faces;

	// Loop through all nodes
	for (int j = 0; j < nRows+1; j++) {
	for (int i = 0; i < nCols+1; i++) {
		double dXLL = dXLLCorner + static_cast<double>(i) * dCellSize;
		double dYLL = dYLLCorner + static_cast<double>(j) * dCellSize;

		double dLon;
		double dLat;

		ConvertUTMtoRLL(nZone, dXLL, dYLL, dLon, dLat);

		//printf("%i %1.5e %1.5e\n", nodes.size(), dLon, dLat);

		double dX = cos(dLat) * cos(dLon);
		double dY = cos(dLat) * sin(dLon);
		double dZ = sin(dLat);

		nodes.push_back(Node(dX, dY, dZ));
	}
	}

	// Generate faces
	for (int j = 0; j < nRows; j++) {

		int iThisLatNodeIx =  j    * (nCols + 1);
		int iNextLatNodeIx = (j+1) * (nCols + 1);

		for (int i = 0; i < nCols; i++) {
			Face face(4);
			face.SetNode(0, iThisLatNodeIx + i);
			face.SetNode(1, iThisLatNodeIx + (i + 1) % (nCols + 1));
			face.SetNode(2, iNextLatNodeIx + (i + 1) % (nCols + 1));
			face.SetNode(3, iNextLatNodeIx + i);

			faces.push_back(face);
		}
	}

/*
	// Reorder the faces
	if (fFlipLatLon) {
		FaceVector faceTemp = mesh.faces;
		mesh.faces.clear();
		for (int i = 0; i < nLongitudes; i++) {
		for (int j = 0; j < nLatitudes; j++) {
			mesh.faces.push_back(faceTemp[j * nLongitudes + i]);
		}
		}
	}
*/
	// Output the mesh
	if (strOutputFile.size()) {

		// Announce
		std::cout << "..Writing mesh to file [" << strOutputFile.c_str() << "] ";
		std::cout << std::endl;

		mesh.Write(strOutputFile);

		NcFile ncOutput(strOutputFile.c_str(), NcFile::Write);
		ncOutput.add_att("rectilinear", "true");

		ncOutput.add_att("rectilinear_dim0_size", nRows);
		ncOutput.add_att("rectilinear_dim1_size", nCols);
		ncOutput.add_att("rectilinear_dim0_name", "rows");
		ncOutput.add_att("rectilinear_dim1_name", "cols");

		ncOutput.close();
	}

	// Announce
	std::cout << "..Mesh generator exited successfully" << std::endl;
	std::cout << "=========================================================";
	std::cout << std::endl;

  return 0;

} catch(Exception & e) {
	Announce(e.ToString().c_str());
	return (0);

} catch(...) {
	return (0);
}
}

///////////////////////////////////////////////////////////////////////////////
