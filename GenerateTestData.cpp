///////////////////////////////////////////////////////////////////////////////
///
///	\file    GenerateTestData.cpp
///	\author  Paul Ullrich
///	\version August 31st, 2014
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
#include "DataMatrix3D.h"

#include "netcdfcpp.h"

#include <cmath>
#include <iostream>

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A class wrapper for a test function.
///	</summary>
class TestFunction {

public:
	///	<summary>
	///		Virtual destructor.
	///	</summary>
	virtual ~TestFunction()
	{ }

	///	<summary>
	///		Evaluate the test function.
	///	</summary>
	virtual double operator()(
		double dLon,
		double dLat
	) = 0;
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A relatively smooth low-order harmonic.
///	</summary>
class TestFunctionY2b2 : public TestFunction {

public:
	///	<summary>
	///		Evaluate the test function.
	///	<summary>
	virtual double operator()(
		double dLon,
		double dLat
	) {
		return (2.0 + cos(dLat) * cos(dLat) * cos(2.0 * dLon));
	}
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A high frequency spherical harmonic.
///	</summary>
class TestFunctionY16b32 : public TestFunction {

public:
	///	<summary>
	///		Evaluate the test function.
	///	</summary>
	virtual double operator()(
		double dLon,
		double dLat
	) {
		return (2.0 + pow(sin(2.0 * dLat), 16.0) * cos(16.0 * dLon));
	}
};

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Stationary vortex fields.
///	</summary>
class TestFunctionVortex : public TestFunction {

public:
	///	<summary>
	///		Find the rotated longitude and latitude of a point on a sphere
	///		with pole at (dLonC, dLatC).
	///	</summary>
	void RotatedSphereCoord(
		double dLonC,
		double dLatC,
		double & dLonT,
		double & dLatT
	) {
		double dSinC = sin(dLatC);
		double dCosC = cos(dLatC);
		double dCosT = cos(dLatT);
		double dSinT = sin(dLatT);
		
		double dTrm  = dCosT * cos(dLonT - dLonC);
		double dX = dSinC * dTrm - dCosC * dSinT;
		double dY = dCosT * sin(dLonT - dLonC);
		double dZ = dSinC * dSinT + dCosC * dTrm;

		dLonT = atan2(dY, dX);
		if (dLonT < 0.0) {
			dLonT += 2.0 * M_PI;
		}
		dLatT = asin(dZ);
	}

	///	<summary>
	///		Evaluate the test function.
	///	</summary>
	virtual double operator()(
		double dLon,
		double dLat
	) {
		const double dLon0 = 0.0;
		const double dLat0 = 0.6;
		const double dR0 = 3.0;
		const double dD = 5.0;
		const double dT = 6.0;

		RotatedSphereCoord(dLon0, dLat0, dLon, dLat);

		double dRho = dR0 * cos(dLat);
		double dVt = 3.0 * sqrt(3.0) / 2.0
			/ cosh(dRho) / cosh(dRho) * tanh(dRho);

		double dOmega;
		if (dRho == 0.0) {
			dOmega = 0.0;
		} else {
			dOmega = dVt / dRho;
		}

		return (1.0 - tanh(dRho / dD * sin(dLon - dOmega * dT)));
	}
};

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {
try {
	// Mesh file to use
	std::string strMeshFile;

	// Test data to use
	int iTestData;

	// Output on GLL grid
	bool fGLL;

	// Degree of polynomial
	int nP;

	// Output variable name
	std::string strVariableName;

	// Output filename
	std::string strTestData;

	// Parse the command line
	BeginCommandLine()
		CommandLineString(strMeshFile, "mesh", "");
		CommandLineInt(iTestData, "test", 1);
		CommandLineBool(fGLL, "gll");
		CommandLineInt(nP, "np", 4);
		CommandLineString(strVariableName, "var", "Psi");
		CommandLineString(strTestData, "out", "testdata.nc");

		ParseCommandLine(argc, argv);
	EndCommandLine(argv)

	// Triangular quadrature nodes (4th order accuracy)
	const int TriQuadraturePoints = 6;

	const double TriQuadratureG[6][3] = {
		{0.108103018168070, 0.445948490915965, 0.445948490915965},
		{0.445948490915965, 0.108103018168070, 0.445948490915965},
		{0.445948490915965, 0.445948490915965, 0.108103018168070},
		{0.816847572980458, 0.091576213509771, 0.091576213509771},
		{0.091576213509771, 0.816847572980458, 0.091576213509771},
		{0.091576213509771, 0.091576213509771, 0.816847572980458}};

	const double TriQuadratureW[6] =
		{0.223381589678011, 0.223381589678011, 0.223381589678011,
		 0.109951743655322, 0.109951743655322, 0.109951743655322};

	// Test data
	TestFunction * pTest;
	if (iTestData == 1) {
		pTest = new TestFunctionY2b2;
	} else if (iTestData == 2) {
		pTest = new TestFunctionY16b32;
	} else if (iTestData == 3) {
		pTest = new TestFunctionVortex;
	} else {
		_EXCEPTIONT("Test index out of range; expected [1,2,3]");
	}

	// Announce
	Announce("=========================================================");

	// Input mesh
	AnnounceStartBlock("Loading Mesh");
	Mesh mesh(strMeshFile);

	// Check for rectilinear Mesh
	NcFile ncMesh(strMeshFile.c_str(), NcFile::ReadOnly);

	// Check for rectilinear attribute
	bool fRectilinear = false;
	for (int a = 0; a < ncMesh.num_atts(); a++) {
		if (strcmp(ncMesh.get_att(a)->name(), "rectilinear") == 0) {
			fRectilinear = true;
			break;
		}
	}

	if (fGLL && fRectilinear) {
		_EXCEPTIONT("--gll cannot be used with rectilinear grids");
	}

	std::vector<long> vecOutputDimSizes;
	std::vector<std::string> vecOutputDimNames;

	// Non-rectilinear output
	if (!fRectilinear) {
		vecOutputDimSizes.resize(1);
		vecOutputDimSizes[0] = mesh.faces.size();

		vecOutputDimNames.resize(1);
		vecOutputDimNames[0] = "ncol";

		Announce("Non-rectilinear mesh detected");

	// Obtain rectilinear attributes
	} else {
		int nDim0Size = ncMesh.get_att("rectilinear_dim0_size")->as_int(0);
		int nDim1Size = ncMesh.get_att("rectilinear_dim1_size")->as_int(0);

		std::string strDim0Name =
			ncMesh.get_att("rectilinear_dim0_name")->as_string(0);
		std::string strDim1Name =
			ncMesh.get_att("rectilinear_dim1_name")->as_string(0);

		vecOutputDimSizes.resize(2);
		vecOutputDimSizes[0] = nDim0Size;
		vecOutputDimSizes[1] = nDim1Size;

		vecOutputDimNames.resize(2);
		vecOutputDimNames[0] = strDim0Name;
		vecOutputDimNames[1] = strDim1Name;

		Announce("Rectilinear mesh detected");
	}

	AnnounceEndBlock("Done");

	// Generate test data
	AnnounceStartBlock("Generating test data");

	// Output data
	DataVector<double> dVar;

	// Sample as element averages
	if (!fGLL) {

		// Calculate element areas
		mesh.CalculateFaceAreas();

		// Resize the array
		dVar.Initialize(mesh.faces.size());

		// Loop through all Faces
		for (int i = 0; i < mesh.faces.size(); i++) {

			const Face & face = mesh.faces[i];

			// Loop through all sub-triangles
			for (int j = 0; j < face.edges.size()-2; j++) {

				const Node & node0 = mesh.nodes[face[0]];
				const Node & node1 = mesh.nodes[face[j+1]];
				const Node & node2 = mesh.nodes[face[j+2]];

				// Triangle area
				Face faceTri(3);
				faceTri.SetNode(0, face[0]);
				faceTri.SetNode(1, face[j+1]);
				faceTri.SetNode(2, face[j+2]);

				double dTriangleArea = CalculateFaceArea(faceTri, mesh.nodes);

				// Calculate the element average
				double dTotalSample = 0.0;

				// Loop through all quadrature points
				for (int k = 0; k < TriQuadraturePoints; k++) {
					Node node(
						  TriQuadratureG[k][0] * node0.x
						+ TriQuadratureG[k][1] * node1.x
						+ TriQuadratureG[k][2] * node2.x,
						  TriQuadratureG[k][0] * node0.y
						+ TriQuadratureG[k][1] * node1.y
						+ TriQuadratureG[k][2] * node2.y,
						  TriQuadratureG[k][0] * node0.z
						+ TriQuadratureG[k][1] * node1.z
						+ TriQuadratureG[k][2] * node2.z);

					double dMagnitude = node.Magnitude();
					node.x /= dMagnitude;
					node.y /= dMagnitude;
					node.z /= dMagnitude;

					double dLon = atan2(node.y, node.x);
					if (dLon < 0.0) {
						dLon += 2.0 * M_PI;
					}
					double dLat = asin(node.z);

					double dSample = (*pTest)(dLon, dLat);

					dTotalSample += dSample * TriQuadratureW[k] * dTriangleArea;
				}

				dVar[i] += dTotalSample / mesh.vecFaceArea[i];
			}
		}

	// Sample at GLL nodes
	} else {

		// Generate grid metadata
		DataMatrix3D<int> dataGLLNodes;
		DataMatrix3D<double> dataGLLJacobian;

		GenerateMetaData(mesh, nP, false, dataGLLNodes, dataGLLJacobian);

		// Number of elements
		int nElements = mesh.faces.size();

		// Verify all elements are quadrilaterals
		for (int k = 0; k < nElements; k++) {
			const Face & face = mesh.faces[k];

			if (face.edges.size() != 4) {
				_EXCEPTIONT("Non-quadrilateral face detected; "
					"incompatible with --gll");
			}
		}

		// Number of unique nodes
		int iMaxNode = 0;
		for (int i = 0; i < nP; i++) {
		for (int j = 0; j < nP; j++) {
		for (int k = 0; k < nElements; k++) {
			if (dataGLLNodes[i][j][k] > iMaxNode) {
				iMaxNode = dataGLLNodes[i][j][k];
			}
		}
		}
		}

		// Resize output array
		vecOutputDimSizes[0] = iMaxNode;

		// Get Gauss-Lobatto quadrature nodes
		DataVector<double> dG;
		DataVector<double> dW;

		GaussLobattoQuadrature::GetPoints(nP, 0.0, 1.0, dG, dW);

		// Allocate data
		dVar.Initialize(iMaxNode+1);

		// Sample data
		for (int k = 0; k < nElements; k++) {

			const Face & face = mesh.faces[k];

			for (int i = 0; i < nP; i++) {
			for (int j = 0; j < nP; j++) {

				// Apply local map
				Node node;
				Node dDx1G;
				Node dDx2G;

				ApplyLocalMap(
					face,
					mesh.nodes,
					dG[i],
					dG[j],
					node,
					dDx1G,
					dDx2G);

				// Sample data at this point
				double dLon = atan2(node.y, node.x);
				if (dLon < 0.0) {
					dLon += 2.0 * M_PI;
				}
				double dLat = asin(node.z);

				double dSample = (*pTest)(dLon, dLat);

				dVar[dataGLLNodes[j][i][k]-1] = dSample;
			}
			}
		}
	}

	AnnounceEndBlock("Done");

	// Output file
	AnnounceStartBlock("Writing results");

	NcFile ncOut(strTestData.c_str(), NcFile::Replace);

	// Add dimensions
	std::vector<NcDim *> vecDimOut;
	for (int d = 0; d < vecOutputDimSizes.size(); d++) {
		vecDimOut.push_back(
			ncOut.add_dim(
				vecOutputDimNames[d].c_str(),
				vecOutputDimSizes[d]));
	}

	// Add variable
	NcVar * varOut =
		ncOut.add_var(
			strVariableName.c_str(),
			ncDouble,
			static_cast<int>(vecOutputDimSizes.size()),
			(const NcDim**)&(vecDimOut[0]));

	// Output data
	varOut->put(&(dVar[0]), &(vecOutputDimSizes[0]));

	AnnounceEndBlock("Done");

	// Delete the test
	delete pTest;

} catch(Exception & e) {
	std::cout << e.ToString() << std::endl;
}
}


///////////////////////////////////////////////////////////////////////////////

