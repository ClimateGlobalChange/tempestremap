#include "TempestConfig.h"
#include "DataMatrix3D.h"
#include "GridElements.h"
#include "OfflineMap.h"
#include <string>

extern "C" {

	// Generate a Cubed-Sphere mesh
	Mesh* GenerateCSMesh(int nResolution, bool fAlt, std::string strOutputFile);

	// Generate a Latitude-Longitude mesh
	Mesh* GenerateRLLMesh(  int nLongitudes, int nLatitudes, 
							double dLonBegin, double dLonEnd, 
							double dLatBegin, double dLatEnd, 
							bool fFlipLatLon, std::string strOutputFile);

	// Generate a Icosahedral-Sphere mesh
	Mesh* GenerateICOMesh(int nResolution, bool fDual, std::string strOutputFile);

	// Generate Lambert-Conic mesh
	Mesh* GenerateLambertConfConicMesh(  int nNCol, int nNRow, 
												double dLon0, double dLat0, 
												double dLat1, double dLat2, 
												double dXLL, double dYLL, double dDX, 
												std::string strOutputFile);

	// Compute the overlap mesh given a source and target mesh file names
	Mesh* GenerateOverlapMesh(std::string strMeshA, std::string strMeshB, std::string strOverlapMesh, std::string strMethod, bool fNoValidate);

	// Compute the overlap mesh given a source and target mesh objects
	// An overload method which takes as arguments the source and target meshes that are pre-loaded into memory
    Mesh* GenerateOverlapWithMeshes(Mesh& meshA, Mesh& meshB, std::string strOverlapMesh, std::string strMethod, bool fNoValidate);

	// New version of the implementation to compute the overlap mesh given a source and target mesh file names
	Mesh* GenerateOverlapMesh_v1(std::string strMeshA, std::string strMeshB, std::string strOverlapMesh, std::string strMethod, bool fNoValidate);

	Mesh* GenerateGLLMetaData(std::string strMesh, int nP, std::string strOutput, DataMatrix3D<int>& dataGLLnodes, DataMatrix3D<double>& dataGLLJacobian);

	OfflineMap* GenerateOfflineMap( std::string strInputMesh, std::string strOutputMesh, std::string strOverlapMesh, 
									std::string strInputMeta, std::string strOutputMeta, 
									std::string strInputType, std::string strOutputType,
									int nPin=4, int nPout=4, 
									bool fBubble=false, int fMonotoneTypeID=0, 
									bool fVolumetric=false, bool fNoConservation=false, bool fNoCheck=false,
									std::string strVariables="", std::string strOutputMap="", 
									std::string strInputData="", std::string strOutputData="",
									std::string strNColName="", bool fOutputDouble=false, 
									std::string strPreserveVariables="", bool fPreserveAll=false, double dFillValueOverride=0.0 );

	int ApplyOfflineMap(std::string strInputData, std::string strInputMap, std::string strVariables, std::string strInputData2, 
						std::string strInputMap2, std::string strVariables2, std::string strOutputData, std::string strNColName, 
						bool fOutputDouble, std::string strPreserveVariables, bool fPreserveAll, double dFillValueOverride);

}
