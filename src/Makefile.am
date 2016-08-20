##############################################################################
lib_LTLIBRARIES = libTempestRemap.la

# Library files to include
LDFILES= -lnetcdf -lnetcdf_c++ -framework accelerate

libTempestRemap_la_includedir = $(includedir)

# The list of source files, and any header files that do not need to be installed
libTempestRemap_la_SOURCES = \
	Announce.cpp \
	PolynomialInterp.cpp \
	GridElements.cpp \
	MeshUtilities.cpp \
	MeshUtilitiesFuzzy.cpp \
	MeshUtilitiesExact.cpp \
	GaussQuadrature.cpp \
	GaussLobattoQuadrature.cpp \
	LegendrePolynomial.cpp \
	FiniteElementTools.cpp \
	NetCDFUtilities.cpp \
	OfflineMap.cpp \
	TriangularQuadrature.cpp

# The list of header files which are to be installed 
nobase_libTempestRemap_la_include_HEADERS = \
	Announce.h \
	DataVector.h \
	FixedPoint.h \
	GridElements.h \
	LinearRemapSE0.h \
	MeshUtilities.h \
	OverlapMesh.h \
	TriangularQuadrature.h \
	CommandLine.h \
	Defines.h \
	GaussLobattoQuadrature.h \
	kdtree.h \
	MathHelper.h \
	NetCDFUtilities.h \
	PolynomialInterp.h \
	DataMatrix3D.h \
	Exception.h \
	GaussQuadrature.h \
	LegendrePolynomial.h \
	MeshUtilitiesExact.h \
	OfflineMap.h \
	SparseMatrix.h \
	DataMatrix.h \
	FiniteElementTools.h \
	GridElementsExact.h \
	LinearRemapFV.h \
	MeshUtilitiesFuzzy.h \
	OverlapFace.h \
	STLStringHelper.h

# Load system-specific defaults
#CFLAGS+= -I$(NETCDF_INCLUDEDIR)
#LDFLAGS+= -L$(NETCDF_LIBDIR)
AM_CPPFLAGS = -I$(srcdir)/src -I$(builddir) -I$(builddir)/src ${NETCDF_CPPFLAGS}
AM_LDFLAGS = ${LDFLAGS} ${NETCDF_LDFLAGS}
LDADD = libTempestRemap.la ${NETCDF_LIBS} ${LAPACK_LIBS} ${BLAS_LIBS} ${LIBS}

GenerateRLLMesh_SOURCES = GenerateRLLMesh.cpp
GenerateCSMesh_SOURCES = GenerateCSMesh.cpp
GenerateICOMesh_SOURCES = GenerateICOMesh.cpp
GenerateLambertConfConicMesh_SOURCES = GenerateLambertConfConicMesh.cpp
GenerateOverlapMesh_SOURCES = GenerateOverlapMesh.cpp OverlapMesh.cpp kdtree.c
GenerateOverlapMesh_v1_SOURCES = GenerateOverlapMesh_v1.cpp OverlapMesh.cpp kdtree.c
GenerateGLLMetaData_SOURCES = GenerateGLLMetaData.cpp
MeshToTxt_SOURCES = MeshToTxt.cpp
GenerateTestData_SOURCES = GenerateTestData.cpp
CalculateDiffNorms_SOURCES = CalculateDiffNorms.cpp
ApplyOfflineMap_SOURCES = ApplyOfflineMap.cpp
GenerateOfflineMap_SOURCES = GenerateOfflineMap.cpp LinearRemapSE0.cpp LinearRemapFV.cpp OverlapMesh.cpp kdtree.c
GenerateOfflineMap_v1_SOURCES = GenerateOfflineMap_v1.cpp LinearRemapSE0.cpp LinearRemapFV.cpp

#include Make.defs

bin_PROGRAMS = GenerateRLLMesh$(EXEEXT) GenerateCSMesh$(EXEEXT) GenerateICOMesh$(EXEEXT) \
							 GenerateLambertConfConicMesh$(EXEEXT) GenerateOverlapMesh$(EXEEXT) \
							 GenerateOverlapMesh_v1$(EXEEXT) GenerateGLLMetaData$(EXEEXT) MeshToTxt$(EXEEXT) \
							 GenerateTestData$(EXEEXT) CalculateDiffNorms$(EXEEXT) \
							 ApplyOfflineMap$(EXEEXT) GenerateOfflineMap$(EXEEXT) GenerateOfflineMap_v1$(EXEEXT)

# DO NOT DELETE

