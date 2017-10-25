#ifndef TEMPESTREMAP_API_H
#define TEMPESTREMAP_API_H

#include "TempestConfig.h"
#include "DataMatrix3D.h"
#include "GridElements.h"
#include "OfflineMap.h"
#include <string>

extern "C" {

    // Generate a Cubed-Sphere mesh
    int GenerateCSMesh ( Mesh& meshOut,
                         int nResolution,
                         bool fAlt,
                         std::string strOutputFile );

    // Generate a Latitude-Longitude mesh
    int GenerateRLLMesh (  Mesh& meshOut,
                           int nLongitudes, int nLatitudes,
                           double dLonBegin, double dLonEnd,
                           double dLatBegin, double dLatEnd,
                           bool fFlipLatLon, bool fForceGlobal,
                           std::string strInputFile, std::string strOutputFile,
                           bool fVerbose );

    // Generate a Icosahedral-Sphere mesh
    int GenerateICOMesh ( Mesh& meshOut,
                          int nResolution,
                          bool fDual,
                          std::string strOutputFile );

    // Generate Lambert-Conic mesh
    int GenerateLambertConfConicMesh ( Mesh& meshOut,
                                       int nNCol, int nNRow,
                                       double dLon0, double dLat0,
                                       double dLat1, double dLat2,
                                       double dXLL, double dYLL, double dDX,
                                       std::string strOutputFile );

    // Compute the overlap mesh given a source and target mesh file names
    int GenerateOverlapMesh ( std::string strMeshA, std::string strMeshB,
                              Mesh& meshOverlap, std::string strOverlapMesh,
                              std::string strMethod, bool fNoValidate,
                              bool fHasConcaveFacesA = false, bool fHasConcaveFacesB = false,
                              bool verbose = true );

    // Compute the overlap mesh given a source and target mesh objects
    // An overload method which takes as arguments the source and target meshes that are pre-loaded into memory
    int GenerateOverlapWithMeshes ( Mesh& meshA, Mesh& meshB,
                                    Mesh& meshOverlap, std::string strOverlapMesh,
                                    std::string strMethod,
                                    bool fHasConcaveFacesA = false, bool fHasConcaveFacesB = false,
                                    bool verbose = true );

    // New version of the implementation to compute the overlap mesh given a source and target mesh file names
    int GenerateOverlapMesh_v1 ( std::string strMeshA, std::string strMeshB,
                                 Mesh& meshOverlap, std::string strOverlapMesh,
                                 std::string strMethod,
                                 const bool fNoValidate = true );

    int GenerateGLLMetaData ( std::string strMesh, Mesh& meshOut,
                              int nP, bool fBubble, std::string strOutput,
                              DataMatrix3D<int>& dataGLLnodes,
                              DataMatrix3D<double>& dataGLLJacobian );

    int GenerateOfflineMap ( OfflineMap& offlineMapOut, std::string strInputMesh, std::string strOutputMesh, std::string strOverlapMesh,
                             std::string strInputMeta, std::string strOutputMeta,
                             std::string strInputType, std::string strOutputType,
                             int nPin = 4, int nPout = 4,
                             bool fBubble = false, int fMonotoneTypeID = 0,
                             bool fVolumetric = false, bool fNoConservation = false, bool fNoCheck = false,
                             std::string strVariables = "", std::string strOutputMap = "",
                             std::string strInputData = "", std::string strOutputData = "",
                             std::string strNColName = "", bool fOutputDouble = false,
                             std::string strPreserveVariables = "", bool fPreserveAll = false, double dFillValueOverride = 0.0,
                             bool fInputConcave = false, bool fOutputConcave = false );

    int GenerateOfflineMapWithMeshes ( OfflineMap& mapRemap,
                                       Mesh& meshInput, Mesh& meshOutput, Mesh& meshOverlap,
                                       std::string strInputMeta, std::string strOutputMeta,
                                       std::string strInputType, std::string strOutputType,
                                       int nPin = 4, int nPout = 4,
                                       bool fBubble = false, int fMonotoneTypeID = 0,
                                       bool fVolumetric = false, bool fNoConservation = false, bool fNoCheck = false,
                                       std::string strVariables = "", std::string strOutputMap = "",
                                       std::string strInputData = "", std::string strOutputData = "",
                                       std::string strNColName = "", bool fOutputDouble = false,
                                       std::string strPreserveVariables = "", bool fPreserveAll = false, double dFillValueOverride = 0.0,
                                       bool fInputConcave = false, bool fOutputConcave = false );

    int ApplyOfflineMap ( std::string strInputData, std::string strInputMap, std::string strVariables, std::string strInputData2,
                          std::string strInputMap2, std::string strVariables2, std::string strOutputData, std::string strNColName,
                          bool fOutputDouble, std::string strPreserveVariables, bool fPreserveAll, double dFillValueOverride );

    int GenerateConnectivityData ( Mesh& meshIn, std::vector< std::set<int> >& vecConnectivity );

}

#endif // TEMPESTREMAP_API_H
