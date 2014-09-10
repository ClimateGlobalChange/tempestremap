///////////////////////////////////////////////////////////////////////////////
///
///	\file    OfflineMap.h
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

#ifndef _OFFLINEMAP_H_
#define _OFFLINEMAP_H_

#include "SparseMatrix.h"
#include <string>
#include <vector>

class Mesh;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		An offline map between two Meshes.
///	</summary>
class OfflineMap {

public:
	///	<summary>
	///		Initialize the array of output dimensions from a file.
	///	</summary>
	void InitializeOutputDimensionsFromFile(
		const std::string & strOutputMesh
	);

public:
	///	<summary>
	///		Apply the offline map to a data file.
	///	</summary>
	void Apply(
		const DataVector<double> & vecAreaInput,
		const DataVector<double> & vecAreaOutput,
		const std::string & strInputDataFile,
		const std::string & strOutputDataFile,
		const std::vector<std::string> & vecVariables,
		const std::string & strNColName
	);

	///	<summary>
	///		Read the OfflineMap from a NetCDF file.
	///	</summary>
	void Read(
		const std::string & strInput
	);

	///	<summary>
	///		Write the OfflineMap to a NetCDF file.
	///	</summary>
	void Write(
		const std::string & strOutput
	);

public:
	///	<summary>
	///		Determine if the map is first-order accurate.
	///	</summary>
	bool IsConsistent(
		double dTolerance
	);

	///	<summary>
	///		Determine if the map is conservative.
	///	</summary>
	bool IsConservative(
		const DataVector<double> & vecInputAreas,
		const DataVector<double> & vecOutputAreas,
		double dTolerance
	);

	///	<summary>
	///		Determine if the map is monotone.
	///	</summary>
	bool IsMonotone(
		double dTolerance
	);

public:
	///	<summary>
	///		Get the SparseMatrix representation of the OfflineMap.
	///	</summary>
	const SparseMatrix<double> & GetSparseMatrix() const {
		return m_mapRemap;
	}

	///	<summary>
	///		Get the SparseMatrix representation of the OfflineMap.
	///	</summary>
	SparseMatrix<double> & GetSparseMatrix() {
		return m_mapRemap;
	}

protected:
	///	<summary>
	///		The SparseMatrix representing this operator.
	///	</summary>
	SparseMatrix<double> m_mapRemap;

	///	<summary>
	///		Vector of dimension sizes for Output.
	///	</summary>
	std::vector<int> m_vecOutputDimSizes;

	///	<summary>
	///		Vector of dimension names for Output.
	///	</summary>
	std::vector<std::string> m_vecOutputDimNames;
};

///////////////////////////////////////////////////////////////////////////////

#endif

