#ifndef TEMPESTREMAP_API_H
#define TEMPESTREMAP_API_H

#include "DataArray3D.h"
#include "GridElements.h"
#include "OfflineMap.h"
#include "netcdfcpp.h"
#include <string>

extern "C" {

	///	<summary>
	///		Generate a cubed-sphere mesh.
	///	</summary>
	int GenerateCSMesh (
		Mesh & meshOut,
		int nResolution,
		std::string strOutputFile,
		std::string strOutputFormat );

	///	<summary>
	///		Generate a transect mesh.
	///	</summary>
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

	///	<summary>
	///		Generate a polar stereographic mesh.
	///	</summary>
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

	///	<summary>
	///		Generate a Latitude-Longitude mesh.
	///	</summary>
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

	///	<summary>
	///		Generate a rectilinear mesh from a file.
	///	</summary>
	int GenerateRectilinearMeshFromFile(
		Mesh & mesh, 
		std::string strInputFile,
		std::string strInputFileLonName,
		std::string strInputFileLatName,
		std::string strOutputFile, 
		std::string strOutputFormat,
		bool fVerbose
	);

	///	<summary>
	///		Restructure 2D data into 1D or vice versa.
	///	</summary>
	int RestructureData(
		std::string strInputFile,
		std::string strVariable,
		std::string strFillValue,
		std::string strRefFile,
		std::string strRefFileLonName,
		std::string strRefFileLatName,
		std::string strOutputFile, 
		std::string strOutputFormat,
		bool fVerbose
	);

	///	<summary>
	///		Generate a Universal Transverse Mercator mesh.
	///	</summary>
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

	///	<summary>
	///		Generate a Icosahedral-Sphere mesh.
	///	</summary>
	int GenerateICOMesh (
		Mesh & meshOut,
		int nResolution,
		bool fDual,
		std::string strOutputFile,
		std::string strOutputFormat );

	///	<summary>
	///		Generate a Lambert-Conic mesh.
	///	</summary>
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

	///	<summary>
	///		Compute the overlap mesh given a source and target mesh file names.
	///	</summary>
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

	///	<summary>
	///		Compute the overlap mesh given two mesh objects.
	///		This is an overloaded method which takes as arguments the source and target
	///		meshes that are pre-loaded into memory.
	///	</summary>
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

	///	<summary>
	///		A structure containing optional arguments for GenerateOfflineMap.
	///	</summary>
	struct GenerateOfflineMapAlgorithmOptions {

	public:
		///	<summary>
		///		Constructor.
		///	</summary>
		GenerateOfflineMapAlgorithmOptions() :
			strOutputMapFile(""),
			strOutputFormat("Netcdf4"),
			strSourceMeta(""),
			strTargetMeta(""),
			fSourceConcave(false),
			fTargetConcave(false),
			nPin(4),
			nPout(4),
			strMethod(""),
			fMonotone(false),
			fNoBubble(false),
			fNoCorrectAreas(false),
			fNoConservation(false),
			fNoCheck(false),
			fSparseConstraints(false)
		{ }

	public:
		///	<summary>
		///		A filename for the output map after its generation.
		///	</summary>
		std::string strOutputMapFile;

		///	<summary>
		///		NetCDF format to use for output.
		///	</summary>
		std::string strOutputFormat;

		///	<summary>
		///		A filename containing source mesh metadata.
		///	</summary>
		std::string strSourceMeta;

		///	<summary>
		///		A filename containing target mesh metadata.
		///	</summary>
		std::string strTargetMeta;

		///	<summary>
		///		Source mesh contains concave Faces.
		///	</summary>
		bool fSourceConcave;

		///	<summary>
		///		Target mesh contains concave Faces.
		///	</summary>
		bool fTargetConcave;
		
		///	<summary>
		///		Input polynomial order.
		///	</summary>
		int nPin;

		///	<summary>
		///		Output polynomial order (only used for finite element output).
		///	</summary>
		int nPout;

		///	<summary>
		///		Method arguments.
		///	</summary>
		std::string strMethod;

		///	<summary>
		///		Generate a monotone map.
		///	</summary>
		bool fMonotone;

		///	<summary>
		///		Do not use a bubble correction for finite elements.
		///	</summary>
		bool fNoBubble;

		///	<summary>
		///		Do not correct the Face areas on the input and output meshes to match
		///		the Face areas on the overlap mesh.
		///	</summary>
		bool fNoCorrectAreas;

		///	<summary>
		///		Do not correct conservation errors.
		///	</summary>
		bool fNoConservation;

		///	<summary>
		///		Do not check the final map.
		///	</summary>
		bool fNoCheck;

		///	<summary>
		///		Use sparse constraints.
		///	</summary>
		bool fSparseConstraints;
	};

	///	<summary>
	///		Generate the OfflineMap between input and output meshes.
	///	</summary>
	int GenerateOfflineMapWithMeshes (
		Mesh & meshSource,
		Mesh & meshTarget,
		Mesh & meshOverlap,
		std::string strSourceType,
		std::string strTargetType,
		const GenerateOfflineMapAlgorithmOptions & optsAlg,
		OfflineMap & mapRemap );

	///	<summary>
	///		Generate the OfflineMap between input and output meshes.
	///	</summary>
	int GenerateOfflineMap (
		std::string strSourceMesh,
		std::string strTargetMesh,
		std::string strOverlapMesh,
		std::string strSourceType,
		std::string strTargetType,
		const GenerateOfflineMapAlgorithmOptions & optsAlg,
		OfflineMap & mapRemap );
/*
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
		std::string strOutputFormat = "Netcdf4",
		std::string strPreserveVariables = "",
		bool fPreserveAll = false,
		double dFillValueOverride = 0.0,
		bool fSourceConcave = false,
		bool fTargetConcave = false,
		bool fSparseConstraints = false);
*/

	///	<summary>
	///		A structure containing optional arguments for outputs from GenerateOfflineMap.
	///	</summary>
	struct ApplyOfflineMapOptions {

	public:
		///	<summary>
		///		Constructor.
		///	</summary>
		ApplyOfflineMapOptions() :
			strInputData(""),
			strOutputData(""),
			strInputDataList(""),
			strOutputDataList(""),
			strVariables(""),
			strNColName("ncol"),
			fOutputDouble(false),
			strOutputFormat("Netcdf4"),
			strPreserveVariables(""),
			fPreserveAll(false),
			dFillValueOverride(0.0),
			strLogDir("")
		{ }

	public:
		///	<summary>
		///		The input data file.
		///	</summary>
		std::string strInputData;

		///	<summary>
		///		The output data file.
		///	</summary>
		std::string strOutputData;

		///	<summary>
		///		A text file containing a list of input data files.
		///	</summary>
		std::string strInputDataList;

		///	<summary>
		///		A text file containing a list of output data files.
		///	</summary>
		std::string strOutputDataList;

		///	<summary>
		///		A list of variables to operate on.
		///	</summary>
		std::string strVariables;

		///	<summary>
		///		The name of the unstructured dimension in the data.
		///	</summary>
		std::string strNColName;

		///	<summary>
		///		A string describing how bounds enforcement should be performed.
		///	</summary>
		std::string strEnforceBounds;

		///	<summary>
		///		Output data using double precision.
		///	</summary>
		bool fOutputDouble;

		///	<summary>
		///		NetCDF format to use for output.
		///	</summary>
		std::string strOutputFormat;

		///	<summary>
		///		List of variables to preserve.
		///	</summary>
		std::string strPreserveVariables;

		///	<summary>
		///		Preserve all output variables.
		///	</summary>
		bool fPreserveAll;

		///	<summary>
		///		Fill value to use for output data.
		///	</summary>
		double dFillValueOverride;

		///	<summary>
		///		A directory for writing log files.
		///	</summary>
		std::string strLogDir;
	};

	///	<summary>
	///		Generate the OfflineMap between input and output meshes.
	///	</summary>
	int GenerateOfflineMapAndApply (
		std::string strSourceMesh,
		std::string strTargetMesh,
		std::string strOverlapMesh,
		std::string strSourceType,
		std::string strTargetType,
		const GenerateOfflineMapAlgorithmOptions & optsAlg,
		const ApplyOfflineMapOptions & optsApply,
		OfflineMap & mapRemap );
/*
	///	<summary>
	///		Generate the OfflineMap between input and output meshes.
	///	</summary>
	int GenerateOfflineMap (
		OfflineMap & mapRemap,
		std::string strSourceMesh,
		std::string strTargetMesh,
		std::string strOverlapMesh,
		std::string strSourceMeta,
		std::string strTargetMeta,
		std::string strSourceType,
		std::string strTargetType,
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
		std::string strOutputFormat = "Netcdf4",
		std::string strPreserveVariables = "",
		bool fPreserveAll = false,
		double dFillValueOverride = 0.0,
		bool fSourceConcave = false,
		bool fTargetConcave = false,
		bool fSparseConstraints = false);

	///	<summary>
	///		Generate the OfflineMap between input and output meshes.
	///	</summary>
	int GenerateOfflineMapWithMeshes (
		OfflineMap & mapRemap,
		Mesh & meshSource,
		Mesh & meshTarget,
		Mesh & meshOverlap,
		std::string strSourceMeta,
		std::string strTargetMeta,
		std::string strSourceType,
		std::string strTargetType,
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
		bool fSourceConcave = false,
		bool fTargetConcave = false,
		bool fSparseConstraints = false);

	///	<summary>
	///		Apply an OfflineMap to a datafile.
	///	</summary>
	int ApplyOfflineMap(
		std::string strInputData,
		std::string strInputDataList,
		std::string strInputMap,
		std::string strVariables,
		std::string strOutputData,
		std::string strOutputDataList,
		std::string strNColName, 
		std::string strEnforceBounds,
		bool fOutputDouble,
		std::string strPreserveVariables,
		bool fPreserveAll,
		double dFillValueOverride,
		std::string strLogDir );
*/

	///	<summary>
	///		Apply an OfflineMap to a datafile.
	///	</summary>
	int ApplyOfflineMap(
		std::string strInputMap,
		const ApplyOfflineMapOptions & optsApply);

	///	<summary>
	///		Generate the connectivity data for a given input file.
	///	</summary>
	int GenerateConnectivityData(
		const Mesh & meshIn,
		std::vector< std::set<int> > & vecConnectivity );

}

#endif // TEMPESTREMAP_API_H

