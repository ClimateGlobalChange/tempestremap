///////////////////////////////////////////////////////////////////////////////
///
///	\file    ShpToMesh.cpp
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
#include "Exception.h"
#include "Announce.h"
#include "order32.h"
#include "MeshUtilitiesFuzzy.h"

#include "netcdfcpp.h"

#include <cmath>
#include <iostream>

///////////////////////////////////////////////////////////////////////////////

static const int32_t SHPFileCodeRef = 0x0000270a;

static const int32_t SHPVersionRef = 1000;

static const int32_t SHPPolygonType = 5;

struct SHPHeader {
	int32_t iFileCode;
	int32_t iUnused[5];
	int32_t iFileLength;
	int32_t iVersion;
	int32_t iShapeType;
};

struct SHPBounds {
	double dXmin;
	double dYmin;
	double dXmax;
	double dYmax;

	double dZmin;
	double dZmax;

	double dMmin;
	double dMmax;
};

struct SHPRecordHeader {
	int32_t iNumber;
	int32_t nLength;
};

struct SHPPolygonHeader {
	double dXmin;
	double dYmin;
	double dXmax;
	double dYmax;
	int32_t nNumParts;
	int32_t nNumPoints;
};

///////////////////////////////////////////////////////////////////////////////

int32_t SwapEndianInt32(const int32_t num) {

	int32_t res;
	const char * pnum = (const char *)(&num);
	char * pres = (char *)(&res);

	pres[0] = pnum[3];
	pres[1] = pnum[2];
	pres[2] = pnum[1];
	pres[3] = pnum[0];

	return res;
}

///////////////////////////////////////////////////////////////////////////////

double SwapEndianDouble(const double num) {

	double res;
  	const char * pnum = (const char *)(&num);
	char * pres = (char *)(&res);

	pres[0] = pnum[7];
	pres[1] = pnum[6];
	pres[2] = pnum[5];
	pres[3] = pnum[4];
	pres[4] = pnum[3];
	pres[5] = pnum[2];
	pres[6] = pnum[1];
	pres[7] = pnum[0];

	return res;
}

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Convert a concave Face into a Convex face.
///	</summary>
///	<returns>
///		true if the Face is convex and has been removed from the mesh.faces
///		vector.
///	</returns>
bool ConvexifyFace(
	Mesh & mesh,
	int iFace
) {
	if ((iFace < 0) || (iFace > mesh.faces.size())) {
		_EXCEPTIONT("Face index out of range");
	}

	Face & face = mesh.faces[iFace];

	const int nEdges = face.edges.size();

	MeshUtilitiesFuzzy meshutils;

	bool fHasReflexNodes = false;

	// Search for reflex nodes on Face
	for (int i = 0; i < nEdges; i++) {
		
		int ixLast = (i + nEdges - 1) % nEdges;
		int ixCurr = i;
		int ixNext = (i + 1) % nEdges;

		const Node & nodeLast = mesh.nodes[face[ixLast]];
		const Node & nodeCurr = mesh.nodes[face[ixCurr]];
		const Node & nodeNext = mesh.nodes[face[ixNext]];

		int iSide = meshutils.FindNodeEdgeSide(
			nodeLast,
			nodeCurr,
			Edge::Type_GreatCircleArc,
			nodeNext);

		if (iSide != (-1)) {
			continue;
		}

		Announce("Reflex node found: %i", face[ixCurr]);

		// Reflex node found; divide mesh at this node
		int ixDividingNode = (-1);
		double dMinDist = (-1.0);
		for (int j = 0; j < nEdges; j++) {
			if ((j == ixLast) || (j == ixCurr) || (j == ixNext)) {
				continue;
			}

			const Node & nodeCandidate = mesh.nodes[face[j]];

			// Check that this Node is in the range of the reflex node
			const int iSide0 = meshutils.FindNodeEdgeSide(
				nodeLast,
				nodeCurr,
				Edge::Type_GreatCircleArc,
				nodeCandidate);

			if (iSide0 == (-1)) {
				continue;
			}

			const int iSide1 = meshutils.FindNodeEdgeSide(
				nodeCurr,
				nodeNext,
				Edge::Type_GreatCircleArc,
				nodeCandidate);

			if (iSide1 == (-1)) {
				continue;
			}

			// Check that this Node is of minimum distance
			Node nodeDelta = nodeCandidate - nodeCurr;

			double dDist = nodeDelta.Magnitude();

			//printf("%i ", face[j]);

			if ((dMinDist < 0.0) || (dDist < dMinDist)) {
				ixDividingNode = j;
				dMinDist = dDist;
			}
		}

		// No dividing node found -- add a Steiner vertex
		if (ixDividingNode == (-1)) {

			Announce("No dividing node found -- adding a Steiner vertex");

			// Find a node that bisects the two angles
			Node nodeBisect;
			nodeBisect.x = 3.0 * nodeCurr.x - nodeLast.x - nodeNext.x;
			nodeBisect.y = 3.0 * nodeCurr.y - nodeLast.y - nodeNext.y;
			nodeBisect.z = 3.0 * nodeCurr.z - nodeLast.z - nodeNext.z;

			double dBisectMag = nodeBisect.Magnitude();
			nodeBisect.x /= dBisectMag;
			nodeBisect.y /= dBisectMag;
			nodeBisect.z /= dBisectMag;

			int ixIntersectEdge;
			Node nodeClosestIntersect;
			dMinDist = (-1.0);

			for (int j = 0; j < nEdges; j++) {
				if ((j == ixCurr) || (j == ixLast)) {
					continue;
				}

				std::vector<Node> vecIntersections;

				bool fCoincident = meshutils.CalculateEdgeIntersectionsSemiClip(
					mesh.nodes[face[j]],
					mesh.nodes[face[(j+1)%nEdges]],
					Edge::Type_GreatCircleArc,
					nodeCurr,
					nodeBisect,
					Edge::Type_GreatCircleArc,
					vecIntersections);

				if (fCoincident) {
					_EXCEPTIONT("Coincident lines detected");
				}

				if (vecIntersections.size() > 1) {
					_EXCEPTIONT("Logic error");

				} else if (vecIntersections.size() == 1) {
					Node nodeDelta = vecIntersections[0] - nodeCurr;
					double dDist = nodeDelta.Magnitude();

					if ((dMinDist == -1.0) || (dDist < dMinDist)) {
						dMinDist = dDist;
						ixIntersectEdge = j;
						nodeClosestIntersect = vecIntersections[0];
					}
				}
			}

			Edge edgeIntersect = face.edges[ixIntersectEdge];

			if (dMinDist == -1.0) {
				_EXCEPTIONT("Logic error: No intersecting lines found");
			}

			int ixNewIntersectNode = mesh.nodes.size();
			mesh.nodes.push_back(nodeClosestIntersect);

			int iFaceSize1 = 1;
			for (int k = (i+1)%nEdges; face[k] != edgeIntersect[0]; k = (k+1)%nEdges) {
				iFaceSize1++;
			}

			Face faceNew1(iFaceSize1 + 2);
			Face faceNew2(nEdges - iFaceSize1 + 1);

			for (int k = 0; k < iFaceSize1 + 1; k++) {
				faceNew1.SetNode(k, face[(i+k)%nEdges]);
				//printf("%i\n", (i+k)%nEdges);
			}
			faceNew1.SetNode(iFaceSize1 + 1, ixNewIntersectNode);
			for (int k = 0; k < nEdges - iFaceSize1; k++) {
				faceNew2.SetNode(k, face[(ixIntersectEdge+k+1)%nEdges]);
				//printf("%i\n", (ixIntersectEdge+k+1)%nEdges);
			}
			faceNew2.SetNode(nEdges - iFaceSize1, ixNewIntersectNode);

			mesh.faces.push_back(faceNew1);
			mesh.faces.push_back(faceNew2);
			mesh.faces.erase(mesh.faces.begin() + iFace);

			int nFaces = mesh.faces.size();
			ConvexifyFace(mesh, nFaces-1);
			ConvexifyFace(mesh, nFaces-2);

		// Divide the mesh
		} else {
			Announce("Dividing node found %i", face[ixDividingNode]);

			int iFaceSize1 = 0;
			for (int k = (i+1)%nEdges; k != ixDividingNode; k = (k+1)%nEdges) {
				iFaceSize1++;
			}

			Face faceNew1(iFaceSize1 + 2);
			Face faceNew2(nEdges - iFaceSize1);

			for (int k = 0; k < iFaceSize1 + 2; k++) {
				faceNew1.SetNode(k, face[(i+k)%nEdges]);
				//printf("%i\n", (i+k)%nEdges);
			}
			for (int k = 0; k < nEdges - iFaceSize1; k++) {
				faceNew2.SetNode(k, face[(ixDividingNode+k)%nEdges]);
				//printf("%i\n", (ixDividingNode+k)%nEdges);
			}

			mesh.faces.push_back(faceNew1);
			mesh.faces.push_back(faceNew2);
			mesh.faces.erase(mesh.faces.begin() + iFace);

			int nFaces = mesh.faces.size();
			ConvexifyFace(mesh, nFaces-1);
			ConvexifyFace(mesh, nFaces-2);
		}

		fHasReflexNodes = true;
		break;
	}

	return fHasReflexNodes;
}

///////////////////////////////////////////////////////////////////////////////

void ConvexifyMesh(
	Mesh & mesh
) {
	char szBuffer[256];

	// Loop through all Faces in the Mesh
	int nFaces = mesh.faces.size();
	for (int f = 0; f < nFaces; f++) {
		sprintf(szBuffer, "Face %i", f);
		AnnounceStartBlock(szBuffer);

		// Adjust current Face index
		bool fConcaveFaceRemoved = ConvexifyFace(mesh, f);
		if (fConcaveFaceRemoved) {
			f--;
			nFaces--;
		}

		AnnounceEndBlock("Done");
	}
}

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {

	NcError error(NcError::silent_nonfatal);

try {
	// Input filename
	std::string strInputFile;

	// Output mesh filename
	std::string strOutputFile;

	// Convexify the mesh
	bool fConvexify;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strInputFile, "in", "");
		CommandLineString(strOutputFile, "out", "");
		CommandLineBool(fConvexify, "convexify");

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

	AnnounceBanner();

	// Load shapefile
	AnnounceStartBlock("Loading shapefile");

	std::ifstream shpfile(
		strInputFile.c_str(), std::ios::in | std::ios::binary);

	SHPHeader shphead;
	shpfile.read((char*)(&shphead), sizeof(SHPHeader));

	if (O32_HOST_ORDER == O32_LITTLE_ENDIAN) {
		shphead.iFileCode = SwapEndianInt32(shphead.iFileCode);
		shphead.iFileLength = SwapEndianInt32(shphead.iFileLength);
	} else if (O32_HOST_ORDER == O32_BIG_ENDIAN) {
		shphead.iVersion = SwapEndianInt32(shphead.iVersion);
		shphead.iShapeType = SwapEndianInt32(shphead.iShapeType);
	} else {
		_EXCEPTIONT("Invalid system Endian");
	}

	if (shphead.iFileCode != SHPFileCodeRef) {
		_EXCEPTIONT("Input file does not appear to be a ESRI Shapefile: "
			"File code mismatch");
	}
	if (shphead.iVersion != SHPVersionRef) {
		_EXCEPTIONT("Input file error: Version mismatch");
	}
	if (shphead.iShapeType != SHPPolygonType) {
		_EXCEPTIONT("Input file error: Polygon type expected");
	}

	SHPBounds shpbounds;
	shpfile.read((char*)(&shpbounds), sizeof(SHPBounds));

	if (O32_HOST_ORDER == O32_BIG_ENDIAN) {
		shpbounds.dXmin = SwapEndianDouble(shpbounds.dXmin);
		shpbounds.dYmin = SwapEndianDouble(shpbounds.dYmin);
		shpbounds.dXmax = SwapEndianDouble(shpbounds.dXmax);
		shpbounds.dYmax = SwapEndianDouble(shpbounds.dXmax);

		shpbounds.dZmin = SwapEndianDouble(shpbounds.dZmin);
		shpbounds.dZmax = SwapEndianDouble(shpbounds.dZmax);

		shpbounds.dMmin = SwapEndianDouble(shpbounds.dMmin);
		shpbounds.dMmax = SwapEndianDouble(shpbounds.dMmax);
	}

	// Current position (in 16-bit words)
	int32_t iCurrentPosition = 50;

	int32_t iPolygonIx = 1;

	// Exodus mesh
	Mesh mesh;

	// Load records
	while (iCurrentPosition < shphead.iFileLength) {

		// Read the record header
		SHPRecordHeader shprechead;
		shpfile.read((char*)(&shprechead), sizeof(SHPRecordHeader));
		if (shpfile.eof()) {
			break;
		}

		if (O32_HOST_ORDER == O32_LITTLE_ENDIAN) {
			shprechead.iNumber = SwapEndianInt32(shprechead.iNumber);
			shprechead.nLength = SwapEndianInt32(shprechead.nLength);
		}

		char szBuffer[128];
		sprintf(szBuffer, "Polygon %i", shprechead.iNumber);
		iPolygonIx++;
		AnnounceStartBlock(szBuffer);

		iCurrentPosition += shprechead.nLength;

		// Read the shape type
		int32_t iShapeType;
		shpfile.read((char*)(&iShapeType), sizeof(int32_t));
		if (shpfile.eof()) {
			break;
		}

		if (O32_HOST_ORDER == O32_BIG_ENDIAN) {
			iShapeType = SwapEndianInt32(iShapeType);
		}
		if (iShapeType != SHPPolygonType) {
			_EXCEPTIONT("Input file error: Record Polygon type expected");
		}

		// Read the polygon header
		SHPPolygonHeader shppolyhead;
		shpfile.read((char*)(&shppolyhead), sizeof(SHPPolygonHeader));
		if (shpfile.eof()) {
			break;
		}

		if (O32_HOST_ORDER == O32_BIG_ENDIAN) {
			shppolyhead.dXmin = SwapEndianDouble(shppolyhead.dXmin);
			shppolyhead.dYmin = SwapEndianDouble(shppolyhead.dYmin);
			shppolyhead.dXmax = SwapEndianDouble(shppolyhead.dXmax);
			shppolyhead.dYmax = SwapEndianDouble(shppolyhead.dYmax);
			shppolyhead.nNumParts = SwapEndianInt32(shppolyhead.nNumParts);
			shppolyhead.nNumPoints = SwapEndianInt32(shppolyhead.nNumPoints);
		}

		// Sanity check
		if (shppolyhead.nNumParts > 0x1000000) {
			_EXCEPTION1("Polygon NumParts exceeds sanity bound (%i)",
				shppolyhead.nNumParts);
		}
		if (shppolyhead.nNumPoints > 0x1000000) {
			_EXCEPTION1("Polygon NumPoints exceeds sanity bound (%i)",
				shppolyhead.nNumPoints);
		}
		Announce("containing %i part(s) with %i points",
			shppolyhead.nNumParts,
			shppolyhead.nNumPoints);
		Announce("Xmin: %3.5f", shppolyhead.dXmin);
		Announce("Ymin: %3.5f", shppolyhead.dYmin);
		Announce("Xmax: %3.5f", shppolyhead.dXmax);
		Announce("Ymax: %3.5f", shppolyhead.dYmax);

		if (shppolyhead.nNumParts != 1) {
			_EXCEPTIONT("Only polygons with 1 part currently supported"
				" in Exodus format");
		}

		DataVector<int32_t> iParts(shppolyhead.nNumParts);
		shpfile.read((char*)&(iParts[0]),
			shppolyhead.nNumParts * sizeof(int32_t));
		if (shpfile.eof()) {
			break;
		}

		DataVector<double> dPoints(shppolyhead.nNumPoints * 2);
		shpfile.read((char*)&(dPoints[0]),
			shppolyhead.nNumPoints * 2 * sizeof(double));
		if (shpfile.eof()) {
			break;
		}

		if (O32_HOST_ORDER == O32_BIG_ENDIAN) {
			for (int i = 0; i < shppolyhead.nNumParts; i++) {
				iParts[i] = SwapEndianInt32(iParts[i]);
			}
			for (int i = 0; i < shppolyhead.nNumPoints * 2; i++) {
				dPoints[i] = SwapEndianDouble(dPoints[i]);
			}
		}

		// Convert to Exodus mesh.  Note that shapefile polygons are specified
		// in clockwise order, whereas Exodus files request polygons to be
		// specified in counter-clockwise order.  Hence we need to reorient
		// the alignment of Faces.
		int nFaces = mesh.faces.size();
		int nNodes = mesh.nodes.size();
		mesh.faces.resize(nFaces+1);
		mesh.nodes.resize(nNodes + shppolyhead.nNumPoints);

		mesh.faces[nFaces] = Face(shppolyhead.nNumPoints);
		for (int i = 0; i < shppolyhead.nNumPoints; i++) {
			double dLonRad = dPoints[2*i] / 180.0 * M_PI;
			double dLatRad = dPoints[2*i+1] / 180.0 * M_PI;

			mesh.nodes[nNodes+i].x = cos(dLatRad) * cos(dLonRad);
			mesh.nodes[nNodes+i].y = cos(dLatRad) * sin(dLonRad);
			mesh.nodes[nNodes+i].z = sin(dLatRad);

			mesh.faces[nFaces].SetNode(
				shppolyhead.nNumPoints - i - 1, nNodes + i);
		}
/*
		Face face5(5);
		face5.SetNode(0, 37);
		face5.SetNode(1, 32);
		face5.SetNode(2, 28);
		face5.SetNode(3, 27);
		face5.SetNode(4, 26);

		mesh.faces[0] = face5;
*/
		AnnounceEndBlock("Done");
	}

	AnnounceEndBlock("Done");

	// Convexify the mesh
	if (fConvexify) {
		AnnounceStartBlock("Convexify mesh");
		ConvexifyMesh(mesh);
		AnnounceEndBlock("Done");
	}

	// Write to file
	AnnounceStartBlock("Write Exodus mesh");

	mesh.Write(strOutputFile);

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


