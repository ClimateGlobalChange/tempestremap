# Don't require GNU-standard files (Changelog, README, etc.)
AUTOMAKE_OPTIONS = foreign subdir-objects

# Don't loose the autoconf include path
ACLOCAL_AMFLAGS = -I config

##############################################################################
lib_LTLIBRARIES = libTempestRemap.la

libTempestRemap_la_includedir = $(includedir)

# The list of header files which are to be installed 
libTempestRemap_la_include_HEADERS = \
	src/Announce.h \
	src/Subscript.h \
	src/CoordTransforms.h \
	src/DataArray1D.h \
	src/FixedPoint.h \
	src/GridElements.h \
	src/LinearRemapSE0.h \
	src/MeshUtilities.h \
	src/OverlapMesh.h \
	src/TriangularQuadrature.h \
	src/CommandLine.h \
	src/Defines.h \
	src/GaussLobattoQuadrature.h \
	src/kdtree.h \
	src/order32.h \
	src/MathHelper.h \
	src/NetCDFUtilities.h \
	src/ncvalues.h \
	src/netcdf.hh \
	src/netcdfcpp.h \
	src/PolynomialInterp.h \
	src/DataArray3D.h \
	src/Exception.h \
	src/GaussQuadrature.h \
	src/LegendrePolynomial.h \
	src/MeshUtilitiesExact.h \
	src/OfflineMap.h \
	src/SparseMatrix.h \
	src/DataArray2D.h \
	src/FiniteElementTools.h \
	src/FiniteVolumeTools.h \
	src/GridElementsExact.h \
	src/LinearRemapFV.h \
	src/MeshUtilitiesFuzzy.h \
	src/OverlapFace.h \
	src/STLStringHelper.h \
	src/TempestRemapAPI.h \
	src/TempestConfig.h \
	src/triangle.h

# The list of source files, and any header files that do not need to be installed
libTempestRemap_la_SOURCES = \
	src/Announce.cpp \
	src/PolynomialInterp.cpp \
	src/GridElements.cpp \
	src/MeshUtilities.cpp \
	src/MeshUtilitiesFuzzy.cpp \
	src/MeshUtilitiesExact.cpp \
	src/GenerateCSMesh.cpp \
	src/GenerateTransectMesh.cpp \
	src/GenerateStereographicMesh.cpp \
	src/GenerateICOMesh.cpp \
	src/GenerateRLLMesh.cpp \
	src/GenerateRectilinearMeshFromFile.cpp \
	src/GenerateUTMMesh.cpp \
	src/GenerateGLLMetaData.cpp \
	src/GenerateLambertConfConicMesh.cpp \
	src/GenerateOverlapMesh.cpp \
	src/GenerateOverlapMesh_v1.cpp \
	src/GaussQuadrature.cpp \
	src/GaussLobattoQuadrature.cpp \
	src/LegendrePolynomial.cpp \
	src/FiniteElementTools.cpp \
	src/FiniteVolumeTools.cpp \
	src/NetCDFUtilities.cpp \
	src/RestructureData.cpp \
	src/ncvalues.cpp \
	src/netcdf.cpp \
	src/OverlapMesh.cpp \
	src/OfflineMap.cpp \
	src/LinearRemapSE0.cpp \
	src/LinearRemapFV.cpp \
	src/TriangularQuadrature.cpp \
	src/ApplyOfflineMap.cpp \
	src/GenerateOfflineMap.cpp \
	src/GenerateConnectivityData.cpp \
	src/kdtree.cpp \
	src/triangle.cpp \
	src/node_multimap_3d.h

libTempestRemap_la_LIBADD = ${NETCDF_LIBS} ${LAPACK_LIBS} ${BLAS_LIBS} ${LIBS}

# Load system-specific defaults
AM_CPPFLAGS = -I$(srcdir)/src -I$(builddir)/src ${NETCDF_CPPFLAGS}
AM_LDFLAGS = ${LDFLAGS} ${NETCDF_LDFLAGS}
LDADD = libTempestRemap.la ${NETCDF_LIBS} ${LAPACK_LIBS} ${BLAS_LIBS} ${LIBS}

# Mesh generation drivers
GenerateTestData_SOURCES = src/GenerateTestData.cpp
GenerateRLLMesh_SOURCES = src/GenerateRLLMeshExe.cpp
GenerateRectilinearMeshFromFile_SOURCES = src/GenerateRectilinearMeshFromFileExe.cpp
GenerateUTMMesh_SOURCES = src/GenerateUTMMeshExe.cpp
GenerateCSMesh_SOURCES = src/GenerateCSMeshExe.cpp
GenerateTransectMesh_SOURCES = src/GenerateTransectMeshExe.cpp
GenerateStereographicMesh_SOURCES = src/GenerateStereographicMeshExe.cpp
GenerateICOMesh_SOURCES = src/GenerateICOMeshExe.cpp
GenerateLambertConfConicMesh_SOURCES = src/GenerateLambertConfConicMeshExe.cpp
GenerateVolumetricMesh_SOURCES = src/GenerateVolumetricMesh.cpp

# Overlap mesh generators
GenerateOverlapMesh_SOURCES = src/GenerateOverlapMeshExe.cpp
GenerateOverlapMesh_v1_SOURCES = src/GenerateOverlapMeshExe_v1.cpp

# Compute and apply the offline mapping weights
ApplyOfflineMap_SOURCES = src/ApplyOfflineMapExe.cpp
GenerateOfflineMap_SOURCES = src/GenerateOfflineMapExe.cpp 
# GenerateOfflineMap_v1_SOURCES = src/GenerateOfflineMap_v1.cpp

GenerateGLLMetaData_SOURCES = src/GenerateGLLMetaDataExe.cpp
GenerateTransposeMap_SOURCES = src/GenerateTransposeMap.cpp
CoarsenRectilinearData_SOURCES = src/CoarsenRectilinearData.cpp
CalculateDiffNorms_SOURCES = src/CalculateDiffNorms.cpp
AnalyzeMap_SOURCES = src/AnalyzeMap.cpp
VerticalInterpolate_SOURCES = src/VerticalInterpolate.cpp
RestructureData_SOURCES = src/RestructureDataExe.cpp

MeshToTxt_SOURCES = src/MeshToTxt.cpp
ShpToMesh_SOURCES = src/ShpToMesh.cpp
ConvertMeshToUGRID_SOURCES = src/ConvertMeshToUGRID.cpp
ConvertMeshToSCRIP_SOURCES = src/ConvertMeshToSCRIP.cpp
ConvertMeshToExodus_SOURCES = src/ConvertMeshToExodus.cpp

bin_PROGRAMS = GenerateTestData \
				GenerateCSMesh GenerateTransectMesh GenerateStereographicMesh GenerateRLLMesh \
				GenerateUTMMesh GenerateICOMesh GenerateRectilinearMeshFromFile \
				GenerateVolumetricMesh GenerateLambertConfConicMesh \
				GenerateOverlapMesh GenerateOverlapMesh_v1 \
				ApplyOfflineMap GenerateOfflineMap \
				CalculateDiffNorms GenerateGLLMetaData \
				GenerateTransposeMap CoarsenRectilinearData \
				MeshToTxt ShpToMesh ConvertMeshToUGRID ConvertMeshToSCRIP ConvertMeshToExodus \
				AnalyzeMap VerticalInterpolate RestructureData


# Utility target: build but don't run tests
build-check:
	$(MAKE) 'TESTS_ENVIRONMENT=: ' check

doc_DATA = README.md doc/GreatCircleArcIntersections.mw doc/GreatCircleLatitudeIntersections.mw

nobase_dist_bin_SCRIPTS = \
	test/areaint.m \
	test/run.sh \
	test/run_fvtogll_cs.sh \
	test/run_fvtogll_ico_diffnorms.sh \
	test/run_glltofv_rll_diffnorms.sh \
	test/testdata_diffnorms_mono.sh \
	test/generate_testdata.sh \
	test/run_fvtofv_ico.sh \
	test/run_fvtogll_cs_diffnorms.sh \
	test/run_fvtogll_rll.sh \
	test/run_glltofv_rll_m.sh \
	test/optimtest.m \
	test/run_fvtofv_ico_diffnorms.sh \
	test/run_fvtogll_cs_m.sh \
	test/run_fvtogll_rll_diffnorms.sh \
	test/run_glltofv_rll_m_diffnorms.sh \
	test/plot_meshes.m \
	test/run_fvtofv_rll.sh \
	test/run_fvtogll_cs_m_diffnorms.sh \
	test/run_glltofv_ico.sh \
	test/run_glltogll_cs.sh \
	test/plot_nodes.m \
	test/run_fvtofv_rll_diffnorms.sh \
	test/run_fvtogll_ico.sh \
	test/run_glltofv_rll.sh \
	test/run_glltogll_cs_diffnorms.sh \
	test/too_many_targets/too_many_targets.sh

GMAKE_FILES = mk/config.make  mk/defs.make  mk/framework.make  mk/rules.make  mk/system.make \
							mk/system/agri.make  mk/system/cheyenne.make  mk/system/cori.make  mk/system/example.make  \
							mk/system/macosx.make

EXTRA_DIST = $(doc_DATA) $(GMAKE_FILES) LICENSE Makefile.gmake src/Makefile.gmake \
						 src/netcdf-cxx-4.2.COPYRIGHT src/netcdf-cxx-4.2.README src/netcdf-cxx-4.2.VERSION

CLEANFILES = configs.sed
