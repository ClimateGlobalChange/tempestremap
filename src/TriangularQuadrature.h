///////////////////////////////////////////////////////////////////////////////
///
///	\file    TriangularQuadrature.h
///	\author  Paul Ullrich
///	\version September 1, 2014
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

#ifndef _TRIANGULARQUADRATURE_H_
#define _TRIANGULARQUADRATURE_H_

#include "DataArray1D.h"
#include "DataArray2D.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A wrapper for various triangular quadrature rules.
///	</summary>
///	<note>
///		Triangular quadrature to use for integration
///		Dunavant, D.A. "High Degree Efficient Symmetrical Gaussian Quadrature
///		Rules for the Triangle."  J. Numer. Meth. Eng., 21, pp 1129-1148.
///	</note>
class TriangularQuadratureRule {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	TriangularQuadratureRule(
		int nOrder
	);

	///	<summary>
	///		Get the number of points in the rule.
	///	</summary>
	int GetPoints() const {
		return m_nPoints;
	}

	///	<summary>
	///		Get the DataArray2D of coordinates.
	///	</summary>
	const DataArray2D<double> & GetG() const {
		return m_dG;
	}

	///	<summary>
	///		Get the DataArray1D of weights.
	///	</summary>
	const DataArray1D<double> & GetW() const {
		return m_dW;
	}

protected:
	///	<summary>
	///		Number of points associated with triangular quadrature rule.
	///	</summary>
	int m_nPoints;

	///	<summary>
	///		Coordinates of triangular quadrature rule.
	///	</summary>
	DataArray2D<double> m_dG;

	///	<summary>
	///		Weights of triangular quadrature rule.
	///	</summary>
	DataArray1D<double> m_dW;
};

///////////////////////////////////////////////////////////////////////////////

#endif

