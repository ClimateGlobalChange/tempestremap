#ifndef TEMPESTREMAP_API_H
#define TEMPESTREMAP_API_H

#include "DataArray3D.h"
#include "GridElements.h"
#include "OfflineMap.h"
#include "netcdfcpp.h"
#include <string>

extern "C" {

	// Generate a Cubed-Sphere mesh
	int GenerateCSMesh (
		Mesh & meshOut,
		int nResolution,
		std::string strOutputFile,
		std::string strOutputFormat );

	// Generate a transect mesh
	int GenerateTransectMesh(
		Mesh & mesh,
		double dLonDeg0,
		double dLatDeg0,
		double dLonDeg1,
		double dLatDeg1,
		double dPerpDtheta,
		int nParaElements,
		int nPerpElements,
		std::string strOutputFile,
		std::string strOutputFormat );

	// Generate a polar stereographic mesh
	int GenerateStereographicMesh(
		Mesh & mesh,
		double dLonDegP,
		double dLatDegP,
		double dLonDeg0,
		double dLatDeg0,
		double dLonDeg1,
		double dLatDeg1,
		int nXElements,
		int nYElements,
		bool fCentroids,
		std::string strOutputFile,
		std::string strOutputFormat );

	// Generate a Latitude-Longitude mesh
	int GenerateRLLMesh (
		Mesh & meshOut,
		int nLongitudes,
		int nLatitudes,
		double dLonBegin,
		double dLonEnd,
		double dLatBegin,
		double dLatEnd,
		bool fGlobalCap,
		bool fFlipLatLon,
		bool fForceGlobal,
		std::string strInputFile,
		std::string strInputFileLonName,
		std::string strInputFileLatName,
		std::string strOutputFile,
		std::string strOutputFormat,
		bool fVerbose );

	// Generate a Latitude-Longitude mesh
	int GenerateUTMMesh (
		Mesh & meshOut,
		int nZone,
		int nCols,
		int nRows,
		double dXLLCorner,
		double dYLLCorner,
		double dCellSize,
		std::string strOutputFile,
		bool fVerbose );

	// Generate a Icosahedral-Sphere mesh
	int GenerateICOMesh (
		Mesh & meshOut,
		int nResolution,
		bool fDual,
		std::string strOutputFile,
		std::string strOutputFormat );

	// Generate Lambert-Conic mesh
	int GenerateLambertConfConicMesh (
		Mesh & meshOut,
		int nNCol,
		int nNRow,
		double dLon0,
		double dLat0,
		double dLat1,
		double dLat2,
		double dXLL,
		double dYLL,
		double dDX,
		std::string strOutputFile );

	// Compute the overlap mesh given a source and target mesh file names
	int GenerateOverlapMesh (
		std::string strMeshA,
		std::string strMeshB,
		Mesh & meshOverlap,
		std::string strOverlapMesh,
		std::string strOutputFormat,
		std::string strMethod,
		bool fNoValidate,
		bool fHasConcaveFacesA = false,
		bool fHasConcaveFacesB = false,
		bool fAllowNoOverlap = false,
		bool fVerbose = true );

	// Compute the overlap mesh given a source and target mesh objects
	// An overload method which takes as arguments the source and target
	// meshes that are pre-loaded into memory
	int GenerateOverlapWithMeshes (
		Mesh & meshA,
		Mesh & meshB,
		Mesh& meshOverlap,
		std::string strOverlapMesh,
		std::string strOutputFormat,
		std::string strMethod,
		bool fHasConcaveFacesA = false,
		bool fHasConcaveFacesB = false,
		bool fAllowNoOverlap = false,
		bool fVerbose = true );

	// Old version of the implementation to compute the overlap mesh
	// given a source and target mesh file names
	int GenerateOverlapMesh_v1 (
		std::string strMeshA,
		std::string strMeshB,
		Mesh& meshOverlap,
		std::string strOverlapMesh,
		std::string strMethod,
		const bool fNoValidate = true );

	// Generate the Gauss-Lobatto-Legendre metadata for the given Mesh
	int GenerateGLLMetaData (
		std::string strMesh,
		Mesh & meshOut,
		int nP,
		bool fNoBubble,
		std::string strOutput,
		DataArray3D<int> & dataGLLnodes,
		DataArray3D<double> & dataGLLJacobian );

	// Generate the OfflineMap between input and output meshes (read from file)
	int GenerateOfflineMap (
		OfflineMap & mapOut,
		std::string strInputMesh,
		std::string strOutputMesh,
		std::string strOverlapMesh,
		std::string strInputMeta,
		std::string strOutputMeta,
		std::string strInputType,
		std::string strOutputType,
		int nPin = 4,
		int nPout = 4,
		bool fNoBubble = false,
		bool fCorrectAreas = false,
		int fMonotoneTypeID = 0,
		bool fVolumetric = false,
		bool fNoConservation = false,
		bool fNoCheck = false,
		std::string strVariables = "",
		std::string strOutputMap = "",
		std::string strInputData = "",
		std::string strOutputData = "",
		std::string strNColName = "",
		bool fOutputDouble = false,
		std::string strOutputFormat = "Classic",
		std::string strPreserveVariables = "",
		bool fPreserveAll = false,
		double dFillValueOverride = 0.0,
		bool fInputConcave = false,
		bool fOutputConcave = false,
		double lb = 0.0,
		double ub = 1.0,
		bool fCAAS = false);

	// Generate the OfflineMap between input and output meshes
	int GenerateOfflineMapWithMeshes (
		OfflineMap & mapRemap,
		Mesh & meshInput,
		Mesh & meshOutput,
		Mesh & meshOverlap,
		std::string strInputMeta,
		std::string strOutputMeta,
		std::string strInputType,
		std::string strOutputType,
		int nPin = 4,
		int nPout = 4,
		bool fBubble = false,
		bool fCorrectAreas = false,
		int fMonotoneTypeID = 0,
		bool fVolumetric = false,
		bool fNoConservation = false,
		bool fNoCheck = false,
		std::string strVariables = "",
		std::string strOutputMap = "",
		std::string strInputData = "",
		std::string strOutputData = "",
		std::string strNColName = "",
		bool fOutputDouble = false,
		std::string strOutputFormat = "Netcdf4",
		std::string strPreserveVariables = "",
		bool fPreserveAll = false,
		double dFillValueOverride = 0.0,
		bool fInputConcave = false,
		bool fOutputConcave = false );

	// Apply an offline map to a datafile
	int ApplyOfflineMap(
		std::string strInputData,
		std::string strInputDataList,
		std::string strInputMap,
		std::string strVariables,
		std::string strOutputData,
		std::string strOutputDataList,
		std::string strNColName, 
		bool fOutputDouble,
		std::string strPreserveVariables,
		bool fPreserveAll,
		double dFillValueOverride,
		std::string strLogDir,
		double lb,
		double ub,
		bool fCAAS );

	// Generate the connectivity data for faces of the given Mesh
	int GenerateConnectivityData(
		const Mesh & meshIn,
		std::vector< std::set<int> > & vecConnectivity);

}

#endif // TEMPESTREMAP_API_H
