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

typedef double Real;
typedef Real   REAL;  // for use with triangles.h

///////////////////////////////////////////////////////////////////////////////
//
// Number of warning messages for individual cells that should be printed
// by AnalyzeMap.
//
static const int OfflineMapWarningMessageCount = 10;

///////////////////////////////////////////////////////////////////////////////
//
// Defines the threshold for the 1-norm condition number, above which
// the finite volume reconstruction drops to first order accuracy.
//
static const Real FVConditionNumberThreshold = 1.0e7;

///////////////////////////////////////////////////////////////////////////////
//
// Defines what type of finite volume reconstruction to use when generating
// maps whose source grid is of type finite volume.
//
#define RECTANGULAR_TRUNCATION
//#define TRIANGULAR_TRUNCATION

///////////////////////////////////////////////////////////////////////////////
//
// Defines for floating point tolerance.
//
static const Real HighTolerance      = 1.0e-10;
static const Real ReferenceTolerance = 1.0e-12;

///////////////////////////////////////////////////////////////////////////////
//
// These defines determine the behavior of GenerateOverlapMesh.
//
// If OVERLAPMESH_RETAIN_REPEATED_NODES is specified this function will make
// no effort to remove repeated nodes during the overlap mesh generation
// calculation.
//
// If OVERLAPMESH_USE_UNSORTED_MAP is specified node removal will use an
// std::unsorted_map() which may nonetheless produce some coincident nodes
// if some very unlikely conditions are met.
//
// If OVERLAPMESH_USE_NODE_MULTIMAP is specified node removal will use the
// node_multimap_3d which is guaranteed to produce no coincident nodes (but
// is the slowest).
//
//#define OVERLAPMESH_RETAIN_REPEATED_NODES
#define OVERLAPMESH_USE_UNSORTED_MAP
//#define OVERLAPMESH_USE_NODE_MULTIMAP

///////////////////////////////////////////////////////////////////////////////
//
// This define specifies the bin width for the std::unsorted_map() and
// node_multimap_3d.
//
#define OVERLAPMESH_BIN_WIDTH 1.0e-1

///////////////////////////////////////////////////////////////////////////////
//
// This define specifies that exact arithmetic should be used in the overlap
// mesh calculation.
//
#define USE_EXACT_ARITHMETIC

///////////////////////////////////////////////////////////////////////////////
//
// Defines required by Triangle package
//
#define ANSI_DECLARATORS
#define VOID int

///////////////////////////////////////////////////////////////////////////////
//
// Use stereographic fit and integration arrays when building a finite volume
// reconstruction.
//
//#define USE_STEREOGRAPHIC_FITS

///////////////////////////////////////////////////////////////////////////////
//
// Maximum number of target faces per source face when imposing
// consistency and conservation.
//
#define FORCECC_MAX_TARGET_FACES (-1)

///////////////////////////////////////////////////////////////////////////////
//
// Maximum number of faces to search before giving up in GenerateOverlapMesh.
//
static const int OverlapFaceSearchMaximumFaces = (-1);

///////////////////////////////////////////////////////////////////////////////

#endif

