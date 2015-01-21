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
#include "DataMatrix.h"
#include "DataMatrix3D.h"
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
	void InitializeSourceDimensionsFromFile(
		const std::string & strSourceMesh
	);

	///	<summary>
	///		Initialize the array of output dimensions from a file.
	///	</summary>
	void InitializeTargetDimensionsFromFile(
		const std::string & strTargetMesh
	);

private:
	///	<summary>
	///		Initialize the coordinate arrays for a finite-volume mesh.
	///	</summary>
	void InitializeCoordinatesFromMeshFV(
		const Mesh & mesh,
		DataVector<double> & dCenterLon,
		DataVector<double> & dCenterLat,
		DataMatrix<double> & dVertexLon,
		DataMatrix<double> & dVertexLat
	);

	///	<summary>
	///		Initialize the coordinate arrays for a finite-element mesh.
	///	</summary>
	void InitializeCoordinatesFromMeshFE(
		const Mesh & mesh,
		int nP,
		const DataMatrix3D<int> & dataGLLnodes,
		DataVector<double> & dCenterLon,
		DataVector<double> & dCenterLat,
		DataMatrix<double> & dVertexLon,
		DataMatrix<double> & dVertexLat
	);

public:
	///	<summary>
	///		Initialize the source coordinate arrays for a finite-volume mesh.
	///	</summary>
	void InitializeSourceCoordinatesFromMeshFV(
		const Mesh & meshSource
	);

	///	<summary>
	///		Initialize the target coordinate arrays for a finite-volume mesh.
	///	</summary>
	void InitializeTargetCoordinatesFromMeshFV(
		const Mesh & meshTarget
	);

	///	<summary>
	///		Initialize the source coordinate arrays for a finite-element mesh.
	///	</summary>
	void InitializeSourceCoordinatesFromMeshFE(
		const Mesh & meshSource,
		int nP,
		const DataMatrix3D<int> & dataGLLnodesSource
	);

	///	<summary>
	///		Initialize the target coordinate arrays for a finite-element mesh.
	///	</summary>
	void InitializeTargetCoordinatesFromMeshFE(
		const Mesh & meshTarget,
		int nP,
		const DataMatrix3D<int> & dataGLLnodesTarget
	);

public:
	///	<summary>
	///		Apply the offline map to a data file.
	///	</summary>
	void Apply(
		const std::string & strSourceDataFile,
		const std::string & strTargetDataFile,
		const std::vector<std::string> & vecVariables,
		const std::string & strNColName,
		bool fTargetDouble = false,
		bool fAppend = false
	);

	///	<summary>
	///		Read the OfflineMap from a NetCDF file.
	///	</summary>
	void Read(
		const std::string & strSource
	);

	///	<summary>
	///		Write the OfflineMap to a NetCDF file.
	///	</summary>
	void Write(
		const std::string & strTarget
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
	DataVector<double> & GetSourceAreas() {
		return m_dSourceAreas;
	}

	///	<summary>
	///		Get the vector of areas associated with degrees of freedom on
	///		the input Mesh.
	///	</summary>
	const DataVector<double> & GetSourceAreas() const {
		return m_dSourceAreas;
	}

	///	<summary>
	///		Get the vector of areas associated with degrees of freedom on
	///		the output Mesh.
	///	</summary>
	DataVector<double> & GetTargetAreas() {
		return m_dTargetAreas;
	}

	///	<summary>
	///		Get the vector of areas associated with degrees of freedom on
	///		the output Mesh.
	///	</summary>
	const DataVector<double> & GetTargetAreas() const {
		return m_dTargetAreas;
	}

	///	<summary>
	///		Set the vector of areas associated with degrees of freedom on
	///		the input Mesh.
	///	</summary>
	void SetSourceAreas(const DataVector<double> & dSourceAreas) {
		m_dSourceAreas = dSourceAreas;
	}

	///	<summary>
	///		Set the vector of areas associated with degrees of freedom on
	///		the output Mesh.
	///	</summary>
	void SetTargetAreas(const DataVector<double> & dTargetAreas) {
		m_dTargetAreas = dTargetAreas;
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
	DataVector<double> m_dSourceAreas;

	///	<summary>
	///		Vector of areas associated with output degrees of freedom.
	///	</summary>
	DataVector<double> m_dTargetAreas;

protected:
	///	<summary>
	///		Vector of cell center longitudes on source grid.
	///	</summary>
	DataVector<double> m_dSourceCenterLon;

	///	<summary>
	///		Vector of cell center latitudes on source grid.
	///	</summary>
	DataVector<double> m_dSourceCenterLat;

	///	<summary>
	///		Vector of cell center longitudes on target grid.
	///	</summary>
	DataVector<double> m_dTargetCenterLon;

	///	<summary>
	///		Vector of cell center latitudes on target grid.
	///	</summary>
	DataVector<double> m_dTargetCenterLat;

	///	<summary>
	///		Matrix of vertex longitudes on source grid.
	///	</summary>
	DataMatrix<double> m_dSourceVertexLon;

	///	<summary>
	///		Matrix of vertex latitudes on source grid.
	///	</summary>
	DataMatrix<double> m_dSourceVertexLat;

	///	<summary>
	///		Matrix of vertex longitudes on source grid.
	///	</summary>
	DataMatrix<double> m_dTargetVertexLon;

	///	<summary>
	///		Matrix of vertex latitudes on source grid.
	///	</summary>
	DataMatrix<double> m_dTargetVertexLat;

protected:
	///	<summary>
	///		Vector of dimension sizes for source.
	///	</summary>
	std::vector<int> m_vecSourceDimSizes;

	///	<summary>
	///		Vector of dimension names for source.
	///	</summary>
	std::vector<std::string> m_vecSourceDimNames;

	///	<summary>
	///		Vector of dimension sizes for target.
	///	</summary>
	std::vector<int> m_vecTargetDimSizes;

	///	<summary>
	///		Vector of dimension names for target.
	///	</summary>
	std::vector<std::string> m_vecTargetDimNames;

	///	<summary>
	///		The fill value override.
	///	</summary>
	float m_flFillValueOverride;
};

///////////////////////////////////////////////////////////////////////////////

#endif

