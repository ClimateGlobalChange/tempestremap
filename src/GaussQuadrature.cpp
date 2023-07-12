///////////////////////////////////////////////////////////////////////////////
///
///	\file    GaussQuadrature.cpp
///	\author  Paul Ullrich
///	\version January 1, 2015
///
///	<remarks>
///		Copyright 2000-2010 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#include "GaussQuadrature.h"
#include "LegendrePolynomial.h"
#include "Exception.h"

///////////////////////////////////////////////////////////////////////////////

void GaussQuadrature::GetPoints(
	int nCount,
	DataArray1D<double> & dG,
	DataArray1D<double> & dW
) {
	// Check for valid range
	if (nCount < 1) {
		_EXCEPTION1("Invalid count (%i): Minimum count 1", nCount);
	}

	// Initialize the arrays
	dG.Allocate(nCount);
	dW.Allocate(nCount);

	// Degree 1
	if (nCount == 1) {
		dG[0] =  0.0;
		dW[0] = +2.0;

	// Degree 2
	} else if (nCount == 2) {
		dG[0] = -0.5773502691896257;
		dG[1] = +0.5773502691896257;

		dW[0] = +1.0;
		dW[1] = +1.0;

	// Degree 3
	} else if (nCount == 3) {
		dG[0] = -0.7745966692414834;
		dG[1] =  0.0;
		dG[2] = +0.7745966692414834;

		dW[0] = +0.5555555555555556;
		dW[1] = +0.8888888888888888;
		dW[2] = +0.5555555555555556;

	// Degree 4
	} else if (nCount == 4) {
		dG[0] = -0.8611363115940526;
		dG[1] = -0.3399810435848563;
		dG[2] = +0.3399810435848563;
		dG[3] = +0.8611363115940526;

		dW[0] = 0.3478548451374538;
		dW[1] = 0.6521451548625461;
		dW[2] = 0.6521451548625461;
		dW[3] = 0.3478548451374538;

	// Degree 5
	} else if (nCount == 5) {
		dG[0] = -0.9061798459386640;
		dG[1] = -0.5384693101056831;
		dG[2] =  0.0;
		dG[3] = +0.5384693101056831;
		dG[4] = +0.9061798459386640;

		dW[0] = 0.2369268850561891;
		dW[1] = 0.4786286704993665;
		dW[2] = 0.5688888888888889;
		dW[3] = 0.4786286704993665;
		dW[4] = 0.2369268850561891;

	// Degree 6
	} else if (nCount == 6) {
		dG[0] = -0.9324695142031521;
		dG[1] = -0.6612093864662645;
		dG[2] = -0.2386191860831969;
		dG[3] = +0.2386191860831969;
		dG[4] = +0.6612093864662645;
		dG[5] = +0.9324695142031521;

		dW[0] = 0.1713244923791704;
		dW[1] = 0.3607615730481386;
		dW[2] = 0.4679139345726910;
		dW[3] = 0.4679139345726910;
		dW[4] = 0.3607615730481386;
		dW[5] = 0.1713244923791704;

	// Degree 7
	} else if (nCount == 7) {
		dG[0] = -0.9491079123427585;
		dG[1] = -0.7415311855993945;
		dG[2] = -0.4058451513773972;
		dG[3] =  0.0;
		dG[4] = +0.4058451513773972;
		dG[5] = +0.7415311855993945;
		dG[6] = +0.9491079123427585;

		dW[0] = 0.1294849661688697;
		dW[1] = 0.2797053914892766;
		dW[2] = 0.3818300505051189;
		dW[3] = 0.4179591836734694;
		dW[4] = 0.3818300505051189;
		dW[5] = 0.2797053914892766;
		dW[6] = 0.1294849661688697;

	// Degree 8
	} else if (nCount == 8) {
		dG[0] = -0.9602898564975363;
		dG[1] = -0.7966664774136267;
		dG[2] = -0.5255324099163290;
		dG[3] = -0.1834346424956498;
		dG[4] = +0.1834346424956498;
		dG[5] = +0.5255324099163290;
		dG[6] = +0.7966664774136267;
		dG[7] = +0.9602898564975363;

		dW[0] = 0.1012285362903763;
		dW[1] = 0.2223810344533745;
		dW[2] = 0.3137066458778873;
		dW[3] = 0.3626837833783620;
		dW[4] = 0.3626837833783620;
		dW[5] = 0.3137066458778873;
		dW[6] = 0.2223810344533745;
		dW[7] = 0.1012285362903763;

	// Degree 9
	} else if (nCount == 9) {
		dG[0] = -1.0;
		dG[1] = -0.899757995411460;
		dG[2] = -0.677186279510738;
		dG[3] = -0.363117463826178;
		dG[4] =  0.0;
		dG[5] = +0.363117463826178;
		dG[6] = +0.677186279510738;
		dG[7] = +0.899757995411460;
		dG[8] = +1.0;

		dW[0] = 0.0277777777777778;
		dW[1] = 0.1654953615608055;
		dW[2] = 0.2745387125001617;
		dW[3] = 0.3464285109730463;
		dW[4] = 0.3715192743764172;
		dW[5] = 0.3464285109730463;
		dW[6] = 0.2745387125001617;
		dW[7] = 0.1654953615608055;
		dW[8] = 0.0277777777777778;

	// Degree 10
	} else if (nCount == 10) {
		dG[0] = -0.9739065285171717;
		dG[1] = -0.8650633666889845;
		dG[2] = -0.6794095682990244;
		dG[3] = -0.4333953941292472;
		dG[4] = -0.1488743389816312;
		dG[5] = +0.1488743389816312;
		dG[6] = +0.4333953941292472;
		dG[7] = +0.6794095682990244;
		dG[8] = +0.8650633666889845;
		dG[9] = +0.9739065285171717;

		dW[0] = 0.0666713443086881;
		dW[1] = 0.1494513491505806;
		dW[2] = 0.2190863625159820;
		dW[3] = 0.2692667193099963;
		dW[4] = 0.2955242247147529;
		dW[5] = 0.2955242247147529;
		dW[6] = 0.2692667193099963;
		dW[7] = 0.2190863625159820;
		dW[8] = 0.1494513491505806;
		dW[9] = 0.0666713443086881;

	// Higher degrees
	} else {
		LegendrePolynomial::AllRoots(nCount, dG);

		for (int k = 0; k < nCount; k++) {
			double dDeriv =
				LegendrePolynomial::EvaluateDerivative(nCount, dG[k]);

			dW[k] = 2.0 / ((1.0 - dG[k] * dG[k]) * dDeriv * dDeriv);
		}

	}
}

///////////////////////////////////////////////////////////////////////////////

void GaussQuadrature::GetPoints(
	int nCount,
	double dXi0,
	double dXi1,
	DataArray1D<double> & dG,
	DataArray1D<double> & dW
) {
	// Get quadrature points in the [-1, 1] reference element
	GetPoints(nCount, dG, dW);

	// Scale quadrature points
	for (int i = 0; i < nCount; i++) {
		dG[i] = dXi0 + 0.5 * (dXi1 - dXi0) * (dG[i] + 1.0);
		dW[i] = 0.5 * (dXi1 - dXi0) * dW[i];
	}
}

///////////////////////////////////////////////////////////////////////////////

