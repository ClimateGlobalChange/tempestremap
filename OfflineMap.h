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
#include "DataVector.h"
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
	///		Initialize the array of input dimensions from a file.
	///	</summary>
	void InitializeInputDimensionsFromFile(
		const std::string & strInputMesh
	);

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
		const std::string & strInputDataFile,
		const std::string & strOutputDataFile,
		const std::vector<std::string> & vecVariables,
		const std::string & strNColName,
		bool fOutputDouble = false,
		bool fAppend = false
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
	///		Get the vector of areas associated with degrees of freedom on
	///		the input Mesh.
	///	</summary>
	DataVector<double> & GetInputAreas() {
		return m_dInputAreas;
	}

	///	<summary>
	///		Get the vector of areas associated with degrees of freedom on
	///		the input Mesh.
	///	</summary>
	const DataVector<double> & GetInputAreas() const {
		return m_dInputAreas;
	}

	///	<summary>
	///		Get the vector of areas associated with degrees of freedom on
	///		the output Mesh.
	///	</summary>
	DataVector<double> & GetOutputAreas() {
		return m_dOutputAreas;
	}

	///	<summary>
	///		Get the vector of areas associated with degrees of freedom on
	///		the output Mesh.
	///	</summary>
	const DataVector<double> & GetOutputAreas() const {
		return m_dOutputAreas;
	}

	///	<summary>
	///		Set the vector of areas associated with degrees of freedom on
	///		the input Mesh.
	///	</summary>
	void SetInputAreas(const DataVector<double> & dInputAreas) {
		m_dInputAreas = dInputAreas;
	}

	///	<summary>
	///		Set the vector of areas associated with degrees of freedom on
	///		the output Mesh.
	///	</summary>
	void SetOutputAreas(const DataVector<double> & dOutputAreas) {
		m_dOutputAreas = dOutputAreas;
	}

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

public:
	///	<summary>
	///		Set the fill value override.
	///	</summary>
	void SetFillValueOverride(float flFillValueOverride) {
		m_flFillValueOverride = flFillValueOverride;
	}

protected:
	///	<summary>
	///		The SparseMatrix representing this operator.
	///	</summary>
	SparseMatrix<double> m_mapRemap;

	///	<summary>
	///		Vector of areas associated with input degrees of freedom.
	///	</summary>
	DataVector<double> m_dInputAreas;

	///	<summary>
	///		Vector of areas associated with output degrees of freedom.
	///	</summary>
	DataVector<double> m_dOutputAreas;

	///	<summary>
	///		Vector of dimension sizes for Input.
	///	</summary>
	std::vector<int> m_vecInputDimSizes;

	///	<summary>
	///		Vector of dimension names for Input.
	///	</summary>
	std::vector<std::string> m_vecInputDimNames;

	///	<summary>
	///		Vector of dimension sizes for Output.
	///	</summary>
	std::vector<int> m_vecOutputDimSizes;

	///	<summary>
	///		Vector of dimension names for Output.
	///	</summary>
	std::vector<std::string> m_vecOutputDimNames;

	///	<summary>
	///		The fill value override.
	///	</summary>
	float m_flFillValueOverride;
};

///////////////////////////////////////////////////////////////////////////////

#endif

