///////////////////////////////////////////////////////////////////////////////
///
///	\file    RestructureData.cpp
///	\author  Paul Ullrich
///	\version March 31, 2021
///
///	<remarks>
///		Copyright 2021 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#include "GridElements.h"
#include "DataArray2D.h"
#include "Exception.h"
#include "Announce.h"
#include "STLStringHelper.h"
#include "NetCDFUtilities.h"
#include "PolynomialInterp.h"

#include <cmath>
#include <cfloat>
#include <iostream>

#include "netcdfcpp.h"

///////////////////////////////////////////////////////////////////////////////

extern "C" 
int RestructureData(
	std::string strInputFile,
	std::string strVariable,
	std::string strFillValue,
	std::string strRefFile,
	std::string strRefFileLonName,
	std::string strRefFileLatName,
	std::string strOutputFile, 
	std::string strOutputFormat,
	bool fVerbose
) {

	NcError error(NcError::silent_nonfatal);

try {

    // Output format
    STLStringHelper::ToLower(strOutputFormat);

	NcFile::FileFormat eOutputFormat =
		GetNcFileFormatFromString(strOutputFormat);
	if (eOutputFormat == NcFile::BadFormat) {
		_EXCEPTION1("Invalid \"out_format\" value (%s), "
			"expected [Classic|Offset64Bits|Netcdf4|Netcdf4Classic]",
			strOutputFormat.c_str());
	}

	// No variable specified
	bool fVariable = (strVariable != "");

	// If blank reference file, use input file
	if (strRefFile == "") {
		strRefFile = strInputFile;
	}

	AnnounceStartBlock("Loading Data");

	// Load reference file
	NcFile ncreffile(strRefFile.c_str(), NcFile::ReadOnly);
	if (!ncreffile.is_valid()) {
		_EXCEPTION1("Unable to load input file \"%s\"", strInputFile.c_str());
	}

	NcVar * varLon = ncreffile.get_var(strRefFileLonName.c_str());
	if (varLon == NULL) {
		_EXCEPTION1("Unable to find variable \"%s\" in input file", strRefFileLonName.c_str());
	}

	NcVar * varLat = ncreffile.get_var(strRefFileLatName.c_str());
	if (varLat == NULL) {
		_EXCEPTION1("Unable to find variable \"%s\" in input file", strRefFileLatName.c_str());
	}

	// Load input file
	NcFile ncinfile(strInputFile.c_str(), NcFile::ReadOnly);
	if (!ncinfile.is_valid()) {
		_EXCEPTION1("Unable to load input file \"%s\"", strInputFile.c_str());
	}

	// Load variable in input file
	NcVar * varIn = NULL;
	if (fVariable) {
		varIn = ncinfile.get_var(strVariable.c_str());
		if (varIn == NULL) {
			_EXCEPTION1("Unable to find variable \"%s\" in input file", strVariable.c_str());
		}
	}

	// Output file
	NcFile ncoutfile(strOutputFile.c_str(), NcFile::Replace, NULL, 0, eOutputFormat);
	if (!ncoutfile.is_valid()) {
		_EXCEPTION1("Unable to open output file \"%s\" for writing",
			strOutputFile.c_str());
	}

	// Load spatial dimension information
	long lSpatialSize = 0;
	long lSpatialDims = varLon->num_dims();
	NcDim * dimSpatial0 = NULL;
	NcDim * dimSpatial1 = NULL;

	if (lSpatialDims == 1) {
		dimSpatial0 = varLon->get_dim(lSpatialDims-1);
		lSpatialSize = dimSpatial0->size();
	} else if (lSpatialDims == 2) {
		dimSpatial0 = varLon->get_dim(lSpatialDims-2);
		dimSpatial1 = varLon->get_dim(lSpatialDims-1);
		lSpatialSize = dimSpatial0->size() * dimSpatial1->size();
	} else {
		_EXCEPTION1("Invalid number of dimensions in longitude variable \"%s\" -- expected 1 or 2 dimensions",
			varLon->name());
	}

	// Verify latitude array is same size
	if (lSpatialDims != varLat->num_dims()) {
		_EXCEPTION2("Longitude variable \"%s\" and latitude variable \"%s\" have mismatched size",
			varLon->name(), varLat->name());
	}
	for (long d = 0; d < lSpatialDims; d++) {
		if (varLon->get_dim(d)->size() != varLat->get_dim(d)->size()) {
			_EXCEPTION3("Longitude variable \"%s\" and latitude variable \"%s\" have mismatched dimension %li size",
				varLon->name(), varLat->name(), d);
		}
	}

	// Load longitude and latitude data
	DataArray1D<double> dLon(lSpatialSize);
	DataArray1D<double> dLat(lSpatialSize);

	if (lSpatialDims == 1) {
		varLon->get(&(dLon[0]), lSpatialSize);
		varLat->get(&(dLat[0]), lSpatialSize);
	} else {
		varLon->get(&(dLon[0]), dimSpatial0->size(), dimSpatial1->size());
		varLat->get(&(dLat[0]), dimSpatial0->size(), dimSpatial1->size());
	}

	// Load _FillValue (if exists)
	double dFillValueLon = DBL_MAX;
	NcAtt * attFillValueLon = varLon->get_att("_FillValue");
	if (attFillValueLon != NULL) {
		dFillValueLon = attFillValueLon->as_double(0);
	}

	double dFillValueLat = DBL_MAX;
	NcAtt * attFillValueLat = varLat->get_att("_FillValue");
	if (attFillValueLat != NULL) {
		dFillValueLat = attFillValueLat->as_double(0);
	}

	double dFillValueVar = DBL_MAX;
	if (fVariable) {
		NcAtt * attFillValueVar = varIn->get_att("_FillValue");
		if (attFillValueVar != NULL) {
			if (strFillValue != "") {
				_EXCEPTIONT("--fillvalue can only be used when variable _FillValue attribute is not already defined");
			}
			dFillValueVar = attFillValueVar->as_double(0);
		} else if (strFillValue != "") {
			dFillValueVar = std::stod(strFillValue);
		}
	}

	// Variable
	long lVariableDims = lSpatialDims;
	NcDim * dimVar0 = NULL;
	NcDim * dimVar1 = NULL;

	if (fVariable) {
		lVariableDims = varIn->num_dims();
		if (lVariableDims == 1) {
			dimVar0 = varIn->get_dim(lVariableDims-1);
		} else if (lVariableDims == 2) {
			dimVar0 = varIn->get_dim(lVariableDims-2);
			dimVar1 = varIn->get_dim(lVariableDims-1);
		} else {
			_EXCEPTIONT("Auxiliary dimensions in variable: Not implemented yet!");
		}
	}

	AnnounceEndBlock("Done");

	// 2D latitude and longitude arrays
	if (lSpatialDims == 2) {

		// 1D variable array -- convert to 2D
		if (lVariableDims == 1) {

			AnnounceStartBlock("Converting 1D data to 2D");

			DataArray1D<float> dData(lSpatialSize);
			if (fVariable) {
				varIn->get(&(dData[0]), dimVar0->size());
			}

			CopyNcVar(ncreffile, ncoutfile, varLon->name());
			CopyNcVar(ncreffile, ncoutfile, varLat->name());

			NcVar * varLonOut = ncoutfile.get_var(varLon->name());
			if (varLonOut == NULL) {
				_EXCEPTION2("Unable to create variable \"%s\" in file \"%s\"",
					strRefFileLonName.c_str(), strOutputFile.c_str());
			}
			NcVar * varLatOut = ncoutfile.get_var(varLat->name());
			if (varLatOut == NULL) {
				_EXCEPTION2("Unable to create variable \"%s\" in file \"%s\"",
					strRefFileLatName.c_str(), strOutputFile.c_str());
			}

			_ASSERT(varLonOut->num_dims() == 2);
			_ASSERT(varLatOut->num_dims() == 2);

			NcDim * dimX = varLonOut->get_dim(0);
			NcDim * dimY = varLonOut->get_dim(1);

			_ASSERT(dimX != NULL);
			_ASSERT(dimY != NULL);

			if (fVariable) {
				NcVar * varOut = ncoutfile.add_var(strVariable.c_str(), ncFloat, dimX, dimY);
				if (varOut == NULL) {
					_EXCEPTION2("Unable to create variable \"%s\" in file \"%s\"",
						strVariable.c_str(), strOutputFile.c_str());
				}

				CopyNcVarAttributes(varIn, varOut);

				if (strFillValue != "") {
					NcAtt * attFillValue = varOut->get_att("_FillValue");
					if (attFillValue == NULL) {
						varOut->add_att("_FillValue", static_cast<float>(dFillValueVar));
					}
				}

				DataArray1D<float> dDataOut(dimX->size() * dimY->size());

				long ixDOF = 0;
				for (long i = 0; i < lSpatialSize; i++) {
					if (dLon[i] != dFillValueLon) {
						dDataOut[i] = dData[ixDOF];
						ixDOF++;
					} else{
						dDataOut[i] = static_cast<float>(dFillValueVar);
					}
				}

				Announce("%li / %li degrees of freedom found", ixDOF, dimX->size() * dimY->size());

				varOut->put(&(dDataOut[0]), dimX->size(), dimY->size());
			}

			AnnounceEndBlock("Done");

		// 2D variable array -- convert to 1D
		} else if (lVariableDims == 2) {

			AnnounceStartBlock("Converting 2D data to 1D");

			DataArray1D<float> dData(lSpatialSize);
			if (fVariable) {
				varIn->get(&(dData[0]), dimVar0->size(), dimVar1->size());
			}

			// Count number of non-fillvalue elements
			long lCount = 0;
			for (long i = 0; i < lSpatialSize; i++) {
				if (dLon[i] != dFillValueLon) {
					lCount++;
				}
			}

			Announce("%li degrees of freedom found", lCount);

			NcDim * dimNCol = ncoutfile.add_dim("ncol", lCount);
			if (dimNCol == NULL) {
				_EXCEPTION1("Unable to create dimension \"ncol\" in file \"%s\"", strOutputFile.c_str());
			}
			NcVar * varLonOut = ncoutfile.add_var(strRefFileLonName.c_str(), ncDouble, dimNCol);
			if (varLonOut == NULL) {
				_EXCEPTION2("Unable to create variable \"%s\" in file \"%s\"",
					strRefFileLonName.c_str(), strOutputFile.c_str());
			}
			NcVar * varLatOut = ncoutfile.add_var(strRefFileLatName.c_str(), ncDouble, dimNCol);
			if (varLatOut == NULL) {
				_EXCEPTION2("Unable to create variable \"%s\" in file \"%s\"",
					strRefFileLatName.c_str(), strOutputFile.c_str());
			}
			CopyNcVarAttributes(varLon, varLonOut);
			CopyNcVarAttributes(varLat, varLatOut);

			NcVar * varOut = NULL;
			if (strVariable != "") {
				varOut = ncoutfile.add_var(strVariable.c_str(), ncFloat, dimNCol);
				if (varOut == NULL) {
					_EXCEPTION2("Unable to create variable \"%s\" in file \"%s\"",
						strVariable.c_str(), strOutputFile.c_str());
				}
				CopyNcVarAttributes(varIn, varOut);
			}

			if (lCount == lSpatialSize) {
				varLonOut->put(&(dLon[0]), lCount);
				varLatOut->put(&(dLat[0]), lCount);
				if (fVariable) {
					varOut->put(&(dData[0]), lCount);
				}

			} else {
				DataArray1D<double> dLonVectorized(lCount);
				DataArray1D<double> dLatVectorized(lCount);
				DataArray1D<float> dDataVectorized(lCount);

				long ix = 0;
				for (long i = 0; i < lSpatialSize; i++) {
					if (dLon[i] != dFillValueLon) {
						dLonVectorized[ix] = dLon[i];
						dLatVectorized[ix] = dLat[i];
						if (fVariable) {
							dDataVectorized[ix] = dData[i];
						}
					}
				}
				varLonOut->put(&(dLonVectorized[0]), lCount);
				varLatOut->put(&(dLatVectorized[0]), lCount);

				if (fVariable) {
					varOut->put(&(dDataVectorized[0]), lCount);
				}
			}

			AnnounceEndBlock("Done");
		}

	// Not implemented
	} else {
		_EXCEPTIONT("Not implemented");
	}

	AnnounceBanner();

	return 0;

} catch(Exception & e) {
	Announce(e.ToString().c_str());
	return (0);

} catch(...) {
	return (0);
}
}

///////////////////////////////////////////////////////////////////////////////

