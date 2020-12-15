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
#include "DataArray1D.h"
#include "DataArray2D.h"
#include "DataArray3D.h"
#include "netcdfcpp.h"
#include <string>
#include <vector>
#include <cfloat>

class Mesh;

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		An offline map between two Meshes.
///	</summary>
class OfflineMap {

public:
	///	<summary>
	///		Constructor.
	///	</summary>
	OfflineMap() :
		m_flFillValueOverride(FLT_MAX),
		m_dFillValueOverride(DBL_MAX)
	{ }

	///	<summary>
	///		An empty virtual destructor.
	///	</summary>
	virtual ~OfflineMap()
	{ }

public:
	///	<summary>
	///		Initialize the array of dimensions from a file.
	///	</summary>
	static void InitializeDimensionsFromMeshFile(
		const std::string & strMeshFile,
		std::vector<std::string> & vecDimNames,
		std::vector<int> & vecDimSizes,
		DataArray1D<double> & dCenterLon,
		DataArray1D<double> & dCenterLat,
		DataArray2D<double> & dVertexLon,
		DataArray2D<double> & dVertexLat
	);

	///	<summary>
	///		Initialize the array of input dimensions from a file.
	///	</summary>
	void InitializeSourceDimensionsFromFile(
		const std::string & strSourceMesh
	);

    ///	<summary>
    ///		Initialize the array of input dimensions from a mesh.
    ///	</summary>
    void InitializeSourceDimensions(
        const std::vector<std::string>& p_srcDimNames,
        const std::vector<int>& p_srcDimSizes
    );

	///	<summary>
	///		Initialize the array of output dimensions from a file.
	///	</summary>
	void InitializeTargetDimensionsFromFile(
		const std::string & strTargetMesh
	);

    ///	<summary>
    ///		Initialize the array of output dimensions from a mesh.
    ///	</summary>
    void InitializeTargetDimensions(
        const std::vector<std::string>& p_tgtDimNames,
        const std::vector<int>& p_tgtDimSizes
    );

protected:

	///	Clip and assured sum function.
	void CAAS(
		DataArray1D<double> & x,
		DataArray1D<double> & l,
		DataArray1D<double> & u,
		double & b
		);

	//double CAAS(
		//double & l,
		//double & u,
		//double & b,
		//int N
	//);

	///	<summary>
	///		Initialize the coordinate arrays for a finite-volume mesh.
	///	</summary>
	void InitializeCoordinatesFromMeshFV(
		const Mesh & mesh,
		DataArray1D<double> & dCenterLon,
		DataArray1D<double> & dCenterLat,
		DataArray2D<double> & dVertexLon,
		DataArray2D<double> & dVertexLat,
		bool fLatLon,
		int nNodesPerFace = 0
	);

	///	<summary>
	///		Initialize the coordinate arrays for a finite-element mesh.
	///	</summary>
	void InitializeCoordinatesFromMeshFE(
		const Mesh & mesh,
		int nP,
		const DataArray3D<int> & dataGLLnodes,
		DataArray1D<double> & dCenterLon,
		DataArray1D<double> & dCenterLat,
		DataArray2D<double> & dVertexLon,
		DataArray2D<double> & dVertexLat
	);

	///	<summary>
	///		Initialize the rectilinear coordinate vectors.
	///	</summary>
	void InitializeRectilinearCoordinateVector(
		int nLon,
		int nLat,
		const DataArray2D<double> & dVertexLon,
		const DataArray2D<double> & dVertexLat,
		bool fLonFirst,
		DataArray1D<double> & dCenterLon,
		DataArray1D<double> & dCenterLat,
		DataArray1D<double> & dVectorCenterLon,
		DataArray1D<double> & dVectorCenterLat,
		DataArray2D<double> & dVectorBoundsLon,
		DataArray2D<double> & dVectorBoundsLat
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
		const DataArray3D<int> & dataGLLnodesSource
	);

	///	<summary>
	///		Initialize the target coordinate arrays for a finite-element mesh.
	///	</summary>
	void InitializeTargetCoordinatesFromMeshFE(
		const Mesh & meshTarget,
		int nP,
		const DataArray3D<int> & dataGLLnodesTarget
	);

public:

	///	<summary>
	///		Get the writable reference to the source element center longitude vector.
	///	</summary>
	DataArray1D<double>& GetSourceCenterLon()
	{
		return m_dSourceCenterLon;
	}

	///	<summary>
	///		Get the const reference to the source element center longitude vector.
	///	</summary>
	const DataArray1D<double>& GetSourceCenterLon() const
	{
		return m_dSourceCenterLon;
	}

	///	<summary>
	///		Get the writable reference to the source element center latitude vector.
	///	</summary>
	DataArray1D<double>& GetSourceCenterLat()
	{
		return m_dSourceCenterLat;
	}

	///	<summary>
	///		Get the const reference to the source element center latitude vector.
	///	</summary>
	const DataArray1D<double>& GetSourceCenterLat() const
	{
		return m_dSourceCenterLat;
	}

	///	<summary>
	///		Get the writable reference to the target element center longitude vector.
	///	</summary>
	DataArray1D<double>& GetTargetCenterLon()
	{
		return m_dTargetCenterLon;
	}

	///	<summary>
	///		Get the const reference to the target element center longitude vector.
	///	</summary>
	const DataArray1D<double>& GetTargetCenterLon() const
	{
		return m_dTargetCenterLon;
	}

	///	<summary>
	///		Get the writable reference to the target element center latitude vector.
	///	</summary>
	DataArray1D<double>& GetTargetCenterLat()
	{
		return m_dTargetCenterLat;
	}

	///	<summary>
	///		Get the const reference to the target element center latitude vector.
	///	</summary>
	const DataArray1D<double>& GetTargetCenterLat() const
	{
		return m_dTargetCenterLat;
	}

	///	<summary>
	///		Get the writable reference to the source vertex longitude matrix.
	///	</summary>
	DataArray2D<double>& GetSourceVertexLon()
	{
		return m_dSourceVertexLon;
	}

	///	<summary>
	///		Get the const reference to the source vertex longitude matrix.
	///	</summary>
	const DataArray2D<double>& GetSourceVertexLon() const
	{
		return m_dSourceVertexLon;
	}

	///	<summary>
	///		Get the writable reference to the target vertex longitude matrix.
	///	</summary>
	DataArray2D<double>& GetTargetVertexLon()
	{
		return m_dTargetVertexLon;
	}

	///	<summary>
	///		Get the const reference to the target vertex longitude matrix.
	///	</summary>
	const DataArray2D<double>& GetTargetVertexLon() const
	{
		return m_dTargetVertexLon;
	}

	///	<summary>
	///		Get the writable reference to the source vertex latitude matrix.
	///	</summary>
	DataArray2D<double>& GetSourceVertexLat()
	{
		return m_dSourceVertexLat;
	}

	///	<summary>
	///		Get the const reference to the source vertex latitude matrix.
	///	</summary>
	const DataArray2D<double>& GetSourceVertexLat() const
	{
		return m_dSourceVertexLat;
	}

	///	<summary>
	///		Get the writable reference to the target vertex latitude natrix.
	///	</summary>
	DataArray2D<double>& GetTargetVertexLat()
	{
		return m_dTargetVertexLat;
	}

	///	<summary>
	///		Get the const reference to the target vertex latitude matrix.
	///	</summary>
	const DataArray2D<double>& GetTargetVertexLat() const
	{
		return m_dTargetVertexLat;
	}

public:
	///	<summary>
	///		Copy a list of variables from a source file to target file.
	///	</summary>
	void PreserveVariables(
		const std::string & strSourceDataFile,
		const std::string & strTargetDataFile,
		const std::vector<std::string> & vecPreserveVariables
	);

	///	<summary>
	///		Copy all non-remapped, non-dimensional variables from source file
	///		to target file.
	///	</summary>
	void PreserveAllVariables(
		const std::string & strSourceDataFile,
		const std::string & strTargetDataFile
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
		Mesh & meshOverlap,
		Mesh & meshInput,
		int & nPin,
		DataArray3D<int> & dataGLLNodesIn,
		DataArray3D<int> & dataGLLNodesOut,
		double lb,
		double ub,
		bool fTargetDouble = false,
		bool fAppend = false,
		bool fCAAS = false,
		bool fCAASLocal = false
	);

	///	<summary>
	///		Allocate the solution data correctly and load field
	///     values from a source file.
	///	</summary>
	void RetrieveFieldData(
		const std::string context, /* "source" or "target" */
		const std::string & strSourceDataFile,
		const std::vector<std::string> & vecVariables,
		const std::string & strNColName,
		std::vector<DataArray1D<double> > & vecSolutions
	);

	///	<summary>
	///		Read the OfflineMap from a NetCDF file.
	///	</summary>
	virtual void Read(
		const std::string & strSource,
		std::map<std::string, std::string> * pmapAttributes = NULL,
		NcFile::FileFormat * peFileFormat = NULL
	);

	///	<summary>
	///		Write the OfflineMap to a NetCDF file, with attribute map.
	///	</summary>
	virtual void Write(
		const std::string & strTarget,
		const std::map<std::string, std::string> & mapAttributes,
		NcFile::FileFormat eFileFormat = NcFile::Classic
	);

	///	<summary>
	///		Initialize a map that is the transverse of the given map.
	///	</summary>
	void SetTranspose(
		const OfflineMap & mapIn
	);

private:
	///	<summary>
	///		Determine if the map is first-order accurate.
	///	</summary>
	virtual int IsConsistent(
		double dTolerance,
		const DataArray1D<int> & dataRows,
		const DataArray1D<int> & dataCols,
		const DataArray1D<double> & dataEntries,
		DataArray1D<double> * pdRowSums = NULL
	);

	///	<summary>
	///		Determine if the map is conservative.
	///	</summary>
	virtual int IsConservative(
		double dTolerance,
		const DataArray1D<int> & dataRows,
		const DataArray1D<int> & dataCols,
		const DataArray1D<double> & dataEntries,
		DataArray1D<double> * pdColSums = NULL
	);

	///	<summary>
	///		Determine if the map is monotone.
	///	</summary>
	virtual int IsMonotone(
		double dTolerance,
		const DataArray1D<int> & dataRows,
		const DataArray1D<int> & dataCols,
		const DataArray1D<double> & dataEntries
	);

public:
	///	<summary>
	///		Determine if the map is first-order accurate.
	///	</summary>
	virtual int IsConsistent(
		double dTolerance
	);

	///	<summary>
	///		Determine if the map is conservative.
	///	</summary>
	virtual int IsConservative(
		double dTolerance
	);

	///	<summary>
	///		Determine if the map is monotone.
	///	</summary>
	virtual int IsMonotone(
		double dTolerance
	);

	///	<summary>
	///		Determine if the map is sane.
	///	</summary>
	virtual bool CheckMap(
		bool fCheckConsistency,
		bool fCheckConservation,
		bool fCheckMonotonicity,
		double dNormalTolerance,
		double dStrictTolerance,
		double dTotalOverlapArea = 0.0
	);

public:
	///	<summary>
	///		Get the vector of areas associated with the source mesh.
	///	</summary>
	DataArray1D<double> & GetSourceAreas() {
		return m_dSourceAreas;
	}

	///	<summary>
	///		Get the vector of areas associated with the source mesh.
	///	</summary>
	const DataArray1D<double> & GetSourceAreas() const {
		return m_dSourceAreas;
	}

	///	<summary>
	///		Get the vector of areas associated with the target mesh.
	///	</summary>
	DataArray1D<double> & GetTargetAreas() {
		return m_dTargetAreas;
	}

	///	<summary>
	///		Get the vector of areas associated with the target mesh.
	///	</summary>
	const DataArray1D<double> & GetTargetAreas() const {
		return m_dTargetAreas;
	}

	///	<summary>
	///		Set the vector of areas associated with the source mesh.
	///	</summary>
	void SetSourceAreas(const DataArray1D<double> & dSourceAreas) {
		m_dSourceAreas = dSourceAreas;
	}

	///	<summary>
	///		Set the vector of areas associated with the target mesh.
	///	</summary>
	void SetTargetAreas(const DataArray1D<double> & dTargetAreas) {
		m_dTargetAreas = dTargetAreas;
	}

public:
	///	<summary>
	///		Get the mask vector associated with the source grid.
	///	</summary>
	const DataArray1D<int> & GetSourceMask() const {
		return m_iSourceMask;
	}

	///	<summary>
	///		Get the mask vector associated with the target grid.
	///	</summary>
	const DataArray1D<int> & GetTargetMask() const {
		return m_iTargetMask;
	}

	///	<summary>
	///		Set the mask vector associated with the source grid.
	///	</summary>
	void SetSourceMask(const DataArray1D<int> & iSourceMask) {
		m_iSourceMask = iSourceMask;
	}

	///	<summary>
	///		Set the mask vector associated with the target grid.
	///	</summary>
	void SetTargetMask(const DataArray1D<int> & iTargetMask) {
		m_iTargetMask = iTargetMask;
	}

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

	/// <summary>
	///		Determine if dimension information has been initialized.
	///	</summary>
	bool AreDimensionsInitialized() const {
		if (m_vecSourceDimSizes.size() == 0) {
			return false;
		}
		if (m_vecTargetDimSizes.size() == 0) {
			return false;
		}
		if (m_vecSourceDimNames.size() != m_vecSourceDimSizes.size()) {
			_EXCEPTIONT("Invalid dimension initialization");
		}
		if (m_vecTargetDimNames.size() != m_vecTargetDimSizes.size()) {
			_EXCEPTIONT("Invalid dimension initialization");
		}
		return true;
	}

public:
	///	<summary>
	///		Set the fill value override (float).
	///	</summary>
	void SetFillValueOverride(float flFillValueOverride) {
		m_flFillValueOverride = flFillValueOverride;
	}

	///	<summary>
	///		Set the fill value override (double).
	///	</summary>
	void SetFillValueOverrideDbl(double dFillValueOverride) {
		m_dFillValueOverride = dFillValueOverride;
	}

protected:
	///	<summary>
	///		The SparseMatrix representing this operator.
	///	</summary>
	SparseMatrix<double> m_mapRemap;

	///	<summary>
	///		Vector of areas associated with source mesh.
	///	</summary>
	DataArray1D<double> m_dSourceAreas;

	///	<summary>
	///		Vector of areas associated with target mesh.
	///	</summary>
	DataArray1D<double> m_dTargetAreas;

	///	<summary>
	///		Vector containing grid mask associated with source mesh.
	///	</summary>
	DataArray1D<int> m_iSourceMask;

	///	<summary>
	///		Vector containing grid mask associated with target mesh.
	///	</summary>
	DataArray1D<int> m_iTargetMask;

protected:
	///	<summary>
	///		Vector of cell center longitudes on source grid.
	///	</summary>
	DataArray1D<double> m_dSourceCenterLon;

	///	<summary>
	///		Vector of cell center latitudes on source grid.
	///	</summary>
	DataArray1D<double> m_dSourceCenterLat;

	///	<summary>
	///		Vector of cell center longitudes on target grid.
	///	</summary>
	DataArray1D<double> m_dTargetCenterLon;

	///	<summary>
	///		Vector of cell center latitudes on target grid.
	///	</summary>
	DataArray1D<double> m_dTargetCenterLat;

	///	<summary>
	///		Matrix of vertex longitudes on source grid.
	///	</summary>
	DataArray2D<double> m_dSourceVertexLon;

	///	<summary>
	///		Matrix of vertex latitudes on source grid.
	///	</summary>
	DataArray2D<double> m_dSourceVertexLat;

	///	<summary>
	///		Vector containing cell center longitude along "lon" dimension.
	///	</sumamry>
	DataArray1D<double> m_dVectorSourceCenterLon;

	///	<summary>
	///		Vector containing cell center latitude along "lat" dimension.
	///	</sumamry>
	DataArray1D<double> m_dVectorSourceCenterLat;

	///	<summary>
	///		Vector containing bounds for longitude along "lon" dimension.
	///	</summary>
	DataArray2D<double> m_dVectorSourceBoundsLon;

	///	<summary>
	///		Vector containing bounds for latitude along "lat" dimension.
	///	</summary>
	DataArray2D<double> m_dVectorSourceBoundsLat;

	///	<summary>
	///		Matrix of vertex longitudes on source grid.
	///	</summary>
	DataArray2D<double> m_dTargetVertexLon;

	///	<summary>
	///		Matrix of vertex latitudes on source grid.
	///	</summary>
	DataArray2D<double> m_dTargetVertexLat;
	
	///	<summary>
	///		Vector containing cell center longitude along "lon" dimension.
	///	</sumamry>
	DataArray1D<double> m_dVectorTargetCenterLon;

	///	<summary>
	///		Vector containing cell center latitude along "lat" dimension.
	///	</sumamry>
	DataArray1D<double> m_dVectorTargetCenterLat;

	///	<summary>
	///		Vector containing bounds for longitude along "lon" dimension.
	///	</summary>
	DataArray2D<double> m_dVectorTargetBoundsLon;

	///	<summary>
	///		Vector containing bounds for latitude along "lat" dimension.
	///	</summary>
	DataArray2D<double> m_dVectorTargetBoundsLat;

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
	///		The fill value override (float).
	///	</summary>
	float m_flFillValueOverride;

	///	<summary>
	///		The fill value override (double).
	///	</summary>
	double m_dFillValueOverride;
};

///////////////////////////////////////////////////////////////////////////////

#endif

