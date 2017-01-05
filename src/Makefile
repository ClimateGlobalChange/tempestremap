# Copyright (c) 2016      Bryce Adelstein-Lelbach aka wash
# Copyright (c) 2000-2016 Paul Ullrich 
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying 
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

# Base directory.
TEMPESTREMAPDIR=..

# Load Makefile framework. 
include $(TEMPESTREMAPDIR)/mk/framework.make

UTIL_FILES= Announce.cpp \
            FiniteElementTools.cpp \
			GaussLobattoQuadrature.cpp \
			GaussQuadrature.cpp \
			GridElements.cpp \
			LegendrePolynomial.cpp \
			LinearRemapFV.cpp \
			LinearRemapSE0.cpp \
			MeshUtilities.cpp \
			MeshUtilitiesExact.cpp \
			MeshUtilitiesFuzzy.cpp \
            NetCDFUtilities.cpp \
            OfflineMap.cpp \
			OverlapMesh.cpp \
			PolynomialInterp.cpp \
			TriangularQuadrature.cpp \
			kdtree.cpp

EXEC_FILES= ApplyOfflineMap.cpp \
            CalculateDiffNorms.cpp \
			GenerateConnectivityFile.cpp \
            GenerateCSMesh.cpp \
            GenerateGLLMetaData.cpp \
            GenerateICOMesh.cpp \
            GenerateLambertConfConicMesh.cpp \
            GenerateOfflineMap.cpp \
            GenerateOverlapMesh.cpp \
            GenerateOverlapMesh_v1.cpp \
            GenerateRLLMesh.cpp \
            GenerateTestData.cpp \
            GenerateVolumetricMesh.cpp \
			MeshToTxt.cpp \
			ShpToMesh.cpp

EXEC_TARGETS= $(EXEC_FILES:%.cpp=%)

FILES= $(UTIL_FILES) $(EXEC_FILES)

.PHONY: all clean

# Build rules. 
all: $(EXEC_TARGETS)

$(EXEC_TARGETS): %: $(UTIL_FILES:%.cpp=$(BUILDDIR)/%.o) $(BUILDDIR)/%.o $(TEMPESTREMAPBASELIB)
	$(CXX) $(LDFLAGS) -o $@ $(UTIL_FILES:%.cpp=$(BUILDDIR)/%.o) $(BUILDDIR)/$*.o $(TEMPESTREMAPBASELIB) $(LIBRARIES)
	mv $@ $(TEMPESTREMAPDIR)/bin

$(TEMPESTREMAPBASELIB): %:
	cd $(TEMPESTREMAPDIR)/src/base; $(MAKE)

# Clean rules.
clean:
	rm -rf $(DEPDIR)
	rm -rf $(BUILDDIR)

# Include dependencies.
-include $(FILES:%.cpp=$(DEPDIR)/%.d)

# DO NOT DELETE
