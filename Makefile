##############################################################################
# Compiler and flags
CC= g++
F90= gfortran
CFLAGS= -O3
FORTFLAGS= -std=legacy

USBLAS= True

# NETCDF library directories
NETCDF_INCLUDEDIR=/usr/local/include
NETCDF_LIBDIR=/usr/local/lib

##############################################################################
# DO NOT MODIFY BELOW THIS LINE
##############################################################################

# Library files to include
LDFILES= -lgfortran -lnetcdf -lnetcdf_c++ -framework accelerate

# Local files
FILES= Announce.cpp \
       PolynomialInterp.cpp \
       GridElements.cpp \
       OverlapMesh.cpp \
       MeshUtilities.cpp \
       MeshUtilitiesFuzzy.cpp \
       MeshUtilitiesExact.cpp \
	   GaussLobattoQuadrature.cpp \
       FiniteElementTools.cpp \
	   NetCDFUtilities.cpp \
       OfflineMap.cpp

#FORTRAN_FILES= dqed.f90

GENERATERLLMESH_FILES= GenerateRLLMesh.cpp $(FILES)

GENERATECSMESH_FILES= GenerateCSMesh.cpp $(FILES)

GENERATEOVERLAPMESH_FILES= GenerateOverlapMesh.cpp $(FILES)

GENERATEGLLMETADATA_FILES= GenerateGLLMetaData.cpp $(FILES)

MESHTOTXT_FILES= MeshToTxt.cpp $(FILES)

GENERATETESTDATA_FILES= GenerateTestData.cpp $(FILES)

GECORE2_FILES= gecore2.cpp LinearRemapSE0.cpp $(FILES)

# Load system-specific defaults
CFLAGS+= -I$(NETCDF_INCLUDEDIR)
LDFLAGS+= -L$(NETCDF_LIBDIR)

include Make.defs

##
## Build instructions
##
all: GenerateRLLMesh GenerateCSMesh GenerateOverlapMesh GenerateGLLMetaData MeshToTxt GenerateTestData gecore2

GenerateRLLMesh: $(GENERATERLLMESH_FILES:%.cpp=$(BUILDDIR)/%.o)
	$(CC) $(LDFLAGS) -o $@ $(GENERATERLLMESH_FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

GenerateCSMesh: $(GENERATECSMESH_FILES:%.cpp=$(BUILDDIR)/%.o)
	$(CC) $(LDFLAGS) -o $@ $(GENERATECSMESH_FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

GenerateOverlapMesh: $(GENERATEOVERLAPMESH_FILES:%.cpp=$(BUILDDIR)/%.o)
	$(CC) $(LDFLAGS) -o $@ $(GENERATEOVERLAPMESH_FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

GenerateGLLMetaData: $(GENERATEGLLMETADATA_FILES:%.cpp=$(BUILDDIR)/%.o)
	$(CC) $(LDFLAGS) -o $@ $(GENERATEGLLMETADATA_FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

MeshToTxt: $(MESHTOTXT_FILES:%.cpp=$(BUILDDIR)/%.o)
	$(CC) $(LDFLAGS) -o $@ $(MESHTOTXT_FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

GenerateTestData: $(GENERATETESTDATA_FILES:%.cpp=$(BUILDDIR)/%.o)
	$(CC) $(LDFLAGS) -o $@ $(GENERATETESTDATA_FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

gecore2: $(GECORE2_FILES:%.cpp=$(BUILDDIR)/%.o) $(FORTRAN_FILES:%.f90=$(BUILDDIR)/%.o)
	$(CC) $(LDFLAGS) -o $@ $(GECORE2_FILES:%.cpp=$(BUILDDIR)/%.o) $(FORTRAN_FILES:%.f90=$(BUILDDIR)/%.o) $(LDFILES)


##
## Clean
##
clean:
	rm -f GenerateRLLMesh GenerateCSMesh GenerateOverlapMesh GenerateGLLMetaData MeshToTxt GenerateTestData gecore2 *.o
	rm -rf $(DEPDIR)
	rm -rf $(BUILDDIR)

##
## Include dependencies
##
include $(FILES:%.cpp=$(DEPDIR)/%.d)

include $(FORTRAN_FILES:%.f90=$(DEPDIR)/%.d)

# DO NOT DELETE

