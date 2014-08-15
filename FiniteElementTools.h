///////////////////////////////////////////////////////////////////////////////
///
///	\file    FiniteElementTools.h
///	\author  Paul Ullrich
///	\version August 14, 2014
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

#include "Defines.h"
#include "DataVector.h"
#include "DataMatrix3D.h"

///////////////////////////////////////////////////////////////////////////////

class Mesh;

///////////////////////////////////////////////////////////////////////////////

double GenerateMetaData(
	const Mesh & mesh,
	int nP,
	bool fBubble,
	DataMatrix3D<int> & dataGLLnodes,
	DataMatrix3D<double> & dataGLLJacobian
);

///////////////////////////////////////////////////////////////////////////////

void GenerateUniqueJacobian(
	const DataMatrix3D<int> & dataGLLnodes,
	const DataMatrix3D<double> & dataGLLJacobian,
	DataVector<double> & dataUniqueJacobian
);

///////////////////////////////////////////////////////////////////////////////

