///////////////////////////////////////////////////////////////////////////////
///
///	\file    GenerateConnectivityFile.cpp
///	\author  Paul Ullrich
///	\version January 4, 2017
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

#include "CommandLine.h"
#include "GridElements.h"
#include "FiniteElementTools.h"
#include "GaussLobattoQuadrature.h"
#include "Exception.h"
#include "Announce.h"

#include "netcdfcpp.h"

#include <cmath>
#include <iostream>

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	NcError error(NcError::silent_nonfatal);

try {
	// Input filename
	std::string strInputFile;

	// Output mesh filename
	std::string strOutputFile;

	// Polynomial degree per element
	int nP = 2;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in", "");
		CommandLineString(strOutputFile, "out", "");
		//CommandLineInt(nP, "np", 2);
		//CommandLineBool(fCGLL, "cgll");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	// Check file names
	if (strInputFile == "") {
		std::cout << "ERROR: No input file specified" << std::endl;
		return (-1);
	}
	if (strOutputFile == "") {
		std::cout << "ERROR: No output file specified" << std::endl;
		return (-1);
	}
	if (nP < 1) {
		std::cout << "ERROR: --np must be >= 2" << std::endl;
		return (-1);
	}

	AnnounceBanner();

	// Load input mesh
	AnnounceStartBlock("Loading input mesh");

	Mesh meshIn(strInputFile);
	meshIn.RemoveZeroEdges();

	AnnounceEndBlock("Done");

	// Construct edge map
	AnnounceStartBlock("Constructing edge map");

	meshIn.ConstructEdgeMap();

	AnnounceEndBlock("Done");

	// Number of elements
	int nElements = meshIn.faces.size();

	// Build connectivity vector using edge map
	AnnounceStartBlock("Constructing connectivity");

	std::vector< std::set<int> > vecConnectivity;
	vecConnectivity.resize(nElements);

	EdgeMapConstIterator iter = meshIn.edgemap.begin();
	for (; iter != meshIn.edgemap.end(); iter++) {

		if ((iter->second[0] != InvalidFace) &&
		    (iter->second[1] != InvalidFace)
		) {
			if ((iter->second[0] < 0) || (iter->second[0] >= nElements)) {
				_EXCEPTION1("Face index (%i) out of range",
					iter->second[0]);
			}
			if ((iter->second[1] < 0) || (iter->second[1] >= nElements)) {
				_EXCEPTION1("Face index (%i) out of range",
					iter->second[1]);
			}

			vecConnectivity[iter->second[0]].insert(iter->second[1]+1);
			vecConnectivity[iter->second[1]].insert(iter->second[0]+1);
		}
	}

	AnnounceEndBlock("Done");	

	// Open output file
	AnnounceStartBlock("Writing connectivity file");

	NcFile ncmesh(strInputFile.c_str(), NcFile::ReadOnly);

	NcVar * varLat = ncmesh.get_var("grid_center_lat");
	NcVar * varLon = ncmesh.get_var("grid_center_lon");

	// Check if center latitudes and longitudes are already available
	DataVector<double> dAllLats;
	DataVector<double> dAllLons;

	bool fConvertLatToDegrees = true;
	bool fConvertLonToDegrees = true;

	if ((varLat == NULL) || (varLon == NULL)) {
		Announce("grid_center_lat not found, recalculating face centers");
	} else {
		Announce("grid_center_lat found in file, loading values");

		if (varLat->get_dim(0)->size() != vecConnectivity.size()) {
			_EXCEPTIONT("grid_center_lat dimension mismatch");
		}
		if (varLon->get_dim(0)->size() != vecConnectivity.size()) {
			_EXCEPTIONT("grid_center_lon dimension mismatch");
		}

		dAllLats.Initialize(vecConnectivity.size());
		varLat->set_cur((long)0);
		varLat->get(dAllLats, vecConnectivity.size());

		NcAtt * attLatUnits = varLat->get_att("units");
		std::string strLatUnits = attLatUnits->as_string(0);
		if (strLatUnits == "degrees") {
			fConvertLatToDegrees = false;
		}

		dAllLons.Initialize(vecConnectivity.size());
		varLon->set_cur((long)0);
		varLon->get(dAllLons, vecConnectivity.size());

		NcAtt * attLonUnits = varLon->get_att("units");
		std::string strLonUnits = attLonUnits->as_string(0);
		if (strLonUnits == "degrees") {
			fConvertLonToDegrees = false;
		}

	}

	// Write connectiivty file
	FILE * fp = fopen(strOutputFile.c_str(), "w");
	fprintf(fp, "%lu\n", vecConnectivity.size());
	for (size_t f = 0; f < vecConnectivity.size(); f++) {
		double dLon;
		double dLat;

		if ((varLat == NULL) || (varLon == NULL)) {
			Node nodeCentroid;
			for (int i = 0; i < meshIn.faces[f].edges.size(); i++) {
				nodeCentroid.x += meshIn.nodes[meshIn.faces[f][i]].x;
				nodeCentroid.y += meshIn.nodes[meshIn.faces[f][i]].y;
				nodeCentroid.z += meshIn.nodes[meshIn.faces[f][i]].z;
			}
			double dMagnitude = nodeCentroid.Magnitude();

			nodeCentroid.x /= dMagnitude;
			nodeCentroid.y /= dMagnitude;
			nodeCentroid.z /= dMagnitude;

			dLon = atan2(nodeCentroid.y, nodeCentroid.x);
			dLat = asin(nodeCentroid.z);

			if (dLon < 0.0) {
				dLon += 2.0 * M_PI;
			}

		} else {
			dLon = dAllLons[f];
			dLat = dAllLats[f];
		}

		if (fConvertLonToDegrees) {
			dLon *= 180.0 / M_PI;
		}
		if (fConvertLatToDegrees) {
			dLat *= 180.0 / M_PI;
		}

		fprintf(fp, "%1.14f,", dLon);
		fprintf(fp, "%1.14f,", dLat);
		fprintf(fp, "%lu", vecConnectivity[f].size());

		std::set<int>::const_iterator iter = vecConnectivity[f].begin();
		for (; iter != vecConnectivity[f].end(); iter++) {
			fprintf(fp, ",%i", *iter);
		}
		if (f != vecConnectivity.size()-1) {
			fprintf(fp,"\n");
		}
	}
	fclose(fp);

	AnnounceEndBlock("Done");

	// Announce
	AnnounceBanner();

	return (0);

} catch(Exception & e) {
	Announce(e.ToString().c_str());
	return (-1);

} catch(...) {
	return (-2);
}
}

///////////////////////////////////////////////////////////////////////////////


