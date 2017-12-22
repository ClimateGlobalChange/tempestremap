///////////////////////////////////////////////////////////////////////////////
///
///	\file    OverlapMesh.h
///	\author  Paul Ullrich
///	\version March 7, 2014
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

#ifndef _DEFINES_H_
#define _DEFINES_H_

///////////////////////////////////////////////////////////////////////////////

#define USE_EXACT_ARITHMETIC

// defines required by Triangle package
#define ANSI_DECLARATORS
#define VOID int

///////////////////////////////////////////////////////////////////////////////

typedef double Real;
typedef Real   REAL;  // for use with triangles.h

static const Real HighTolerance      = 1.0e-10;
static const Real ReferenceTolerance = 1.0e-12;

///////////////////////////////////////////////////////////////////////////////

#endif

