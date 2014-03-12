##############################################################################
# Compiler and flags
CC= g++
CFLAGS= -O3

# NETCDF library directories
NETCDF_INCLUDEDIR=/usr/local/include
NETCDF_LIBDIR=/usr/local/lib

##############################################################################
# DO NOT MODIFY BELOW THIS LINE
##############################################################################

# Library files to include
LDFILES= -lnetcdf -lnetcdf_c++

# Local files
FILES= Announce.cpp \
       GridElements.cpp \
       OverlapMesh.cpp

GENERATERLLMESH_FILES= GenerateRLLMesh.cpp $(FILES)

GENERATECSMESH_FILES= GenerateCSMesh.cpp $(FILES)

MESHTOTXT_FILES= MeshToTxt.cpp $(FILES)

UNITTEST_FILES= UnitTest.cpp $(FILES)

GECORE2_FILES= gecore2.cpp $(FILES)

# Load system-specific defaults
CFLAGS+= -I$(NETCDF_INCLUDEDIR)
LDFLAGS+= -L$(NETCDF_LIBDIR)

include Make.defs

##
## Build instructions
##
all: GenerateRLLMesh GenerateCSMesh MeshToTxt UnitTest gecore2

GenerateRLLMesh: $(GENERATERLLMESH_FILES:%.cpp=$(BUILDDIR)/%.o)
	$(CC) $(LDFLAGS) -o $@ $(GENERATERLLMESH_FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

GenerateCSMesh: $(GENERATECSMESH_FILES:%.cpp=$(BUILDDIR)/%.o)
	$(CC) $(LDFLAGS) -o $@ $(GENERATECSMESH_FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

MeshToTxt: $(MESHTOTXT_FILES:%.cpp=$(BUILDDIR)/%.o)
	$(CC) $(LDFLAGS) -o $@ $(MESHTOTXT_FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

UnitTest: $(UNITTEST_FILES:%.cpp=$(BUILDDIR)/%.o)
	$(CC) $(LDFLAGS) -o $@ $(UNITTEST_FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)

gecore2: $(GECORE2_FILES:%.cpp=$(BUILDDIR)/%.o)
	$(CC) $(LDFLAGS) -o $@ $(GECORE2_FILES:%.cpp=$(BUILDDIR)/%.o) $(LDFILES)


##
## Clean
##
clean:
	rm -f GenerateRLLMesh GenerateCSMesh MeshToTxt UnitTest gecore2 *.o
	rm -rf $(DEPDIR)
	rm -rf $(BUILDDIR)

##
## Include dependencies
##
include $(FILES:%.cpp=$(DEPDIR)/%.d)

# DO NOT DELETE

