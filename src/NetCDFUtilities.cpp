///////////////////////////////////////////////////////////////////////////////
///
///	\file    NetCDFUtilities.cpp
///	\author  Paul Ullrich
///	\version August 14, 2014
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

#include "NetCDFUtilities.h"
#include "Exception.h"
#include "DataArray1D.h"
#include "netcdfcpp.h"

#include <vector>

////////////////////////////////////////////////////////////////////////////////

 NcFile::FileFormat GetNcFileFormatFromString(
	const std::string & strFormat
) {
	NcFile::FileFormat eOutputFormat = NcFile::BadFormat;
    if (strFormat == "classic") {
        eOutputFormat = NcFile::Classic;
    } else if (strFormat == "offset64bits") {
        eOutputFormat = NcFile::Offset64Bits;
    } else if (strFormat == "netcdf4") {
        eOutputFormat = NcFile::Netcdf4;
    } else if (strFormat == "netcdf4classic") {
        eOutputFormat = NcFile::Netcdf4Classic;
    }
	return eOutputFormat;
}

////////////////////////////////////////////////////////////////////////////////

void CopyNcFileAttributes(
	NcFile * fileIn,
	NcFile * fileOut
) {
	for (int a = 0; a < fileIn->num_atts(); a++) {
		NcAtt * att = fileIn->get_att(a);
		long num_vals = att->num_vals();

		NcValues * pValues = att->values();

		if (att->type() == ncByte) {
			fileOut->add_att(att->name(), num_vals,
				(const ncbyte*)(pValues->base()));

		} else if (att->type() == ncChar) {
			fileOut->add_att(att->name(), num_vals,
				(const char*)(pValues->base()));

		} else if (att->type() == ncShort) {
			fileOut->add_att(att->name(), num_vals,
				(const short*)(pValues->base()));

		} else if (att->type() == ncInt) {
			fileOut->add_att(att->name(), num_vals,
				(const int*)(pValues->base()));

		} else if (att->type() == ncFloat) {
			fileOut->add_att(att->name(), num_vals,
				(const float*)(pValues->base()));

		} else if (att->type() == ncDouble) {
			fileOut->add_att(att->name(), num_vals,
				(const double*)(pValues->base()));

		} else if (att->type() == ncInt64) {
			fileOut->add_att(att->name(), num_vals,
				(const ncint64*)(pValues->base()));

		} else if (att->type() == ncUInt64) {
			fileOut->add_att(att->name(), num_vals,
				(const ncuint64*)(pValues->base()));

		} else {
			_EXCEPTIONT("Invalid attribute type");
		}

		delete pValues;
	}
}

////////////////////////////////////////////////////////////////////////////////

void CopyNcVarAttributes(
	NcVar * varIn,
	NcVar * varOut
) {
	bool fSuccess;

	for (int a = 0; a < varIn->num_atts(); a++) {
		NcAtt * att = varIn->get_att(a);
		long num_vals = att->num_vals();

		NcValues * pValues = att->values();
		if (pValues == NULL) {
			_EXCEPTION2("Invalid attribute type \"%s::%s\"",
				varIn->name(), att->name());
		}

		// Change _FillValue type to match variable type, if needed
		if (strcmp(att->name(), "_FillValue") == 0) {

			if ((att->type() == ncFloat) && (varOut->type() == ncDouble)) {
				double dFillValue = static_cast<double>(*((const float *)pValues->base()));
				fSuccess = varOut->add_att("_FillValue", dFillValue);
				if (!fSuccess) {
					_EXCEPTION1("Error creating attribute \"_FillValue\" for variable \"%s\" in output file", varOut->name());
				}

				delete pValues;
				continue;

			} else if ((att->type() == ncDouble) && (varOut->type() == ncFloat)) {
				float dFillValue = static_cast<float>(*((const double *)pValues->base()));
				fSuccess = varOut->add_att("_FillValue", dFillValue);
				if (!fSuccess) {
					_EXCEPTION1("Error creating attribute \"_FillValue\" for variable \"%s\" in output file", varOut->name());
				}

				delete pValues;
				continue;
			}
		}

		if (att->type() == ncByte) {
			fSuccess = varOut->add_att(att->name(), num_vals,
				(const ncbyte*)(pValues->base()));

		} else if (att->type() == ncChar) {
			fSuccess = varOut->add_att(att->name(), num_vals,
				(const char*)(pValues->base()));

		} else if (att->type() == ncShort) {
			fSuccess = varOut->add_att(att->name(), num_vals,
				(const short*)(pValues->base()));

		} else if (att->type() == ncInt) {
			fSuccess = varOut->add_att(att->name(), num_vals,
				(const int*)(pValues->base()));

		} else if (att->type() == ncFloat) {
			fSuccess = varOut->add_att(att->name(), num_vals,
				(const float*)(pValues->base()));

		} else if (att->type() == ncDouble) {
			fSuccess = varOut->add_att(att->name(), num_vals,
				(const double*)(pValues->base()));

		} else if (att->type() == ncInt64) {
			fSuccess = varOut->add_att(att->name(), num_vals,
				(const ncint64*)(pValues->base()));

		} else if (att->type() == ncUInt64) {
			fSuccess = varOut->add_att(att->name(), num_vals,
				(const ncuint64*)(pValues->base()));

		} else {
			_EXCEPTIONT("Invalid attribute type");
		}

		if (!fSuccess) {
			_EXCEPTION2("Error creating attribute \"%s\" for variable \"%s\" in output file",
				att->name(), varOut->name());
		}

		delete pValues;
	}
}

////////////////////////////////////////////////////////////////////////////////

void CopyNcVar(
	NcFile & ncIn,
	NcFile & ncOut,
	const std::string & strVarName,
	bool fCopyAttributes
) {
	if (!ncIn.is_valid()) {
		_EXCEPTIONT("Invalid input file specified");
	}
	if (!ncOut.is_valid()) {
		_EXCEPTIONT("Invalid output file specified");
	}
	NcVar * var = ncIn.get_var(strVarName.c_str());
	if (var == NULL) {
		_EXCEPTION1("NetCDF file does not contain variable \"%s\"",
			strVarName.c_str());
	}

	NcVar * varOut;

	std::vector<NcDim *> dimOut;
	dimOut.resize(var->num_dims());

	std::vector<long> counts;
	counts.resize(var->num_dims());

	long nDataSize = 1;

	for (int d = 0; d < var->num_dims(); d++) {
		NcDim * dimA = var->get_dim(d);

		dimOut[d] = ncOut.get_dim(dimA->name());

		if (dimOut[d] == NULL) {
			if (dimA->is_unlimited()) {
				dimOut[d] = ncOut.add_dim(dimA->name());
			} else {
				dimOut[d] = ncOut.add_dim(dimA->name(), dimA->size());
			}

			if (dimOut[d] == NULL) {
				_EXCEPTION2("Failed to add dimension \"%s\" (%i) to file",
					dimA->name(), dimA->size());
			}
		}
		if (dimOut[d]->size() != dimA->size()) {
			if (dimA->is_unlimited() && !dimOut[d]->is_unlimited()) {
				_EXCEPTION2("Mismatch between input file dimension \"%s\" and "
					"output file dimension (UNLIMITED / %i)",
					dimA->name(), dimOut[d]->size());
			} else if (!dimA->is_unlimited() && dimOut[d]->is_unlimited()) {
				_EXCEPTION2("Mismatch between input file dimension \"%s\" and "
					"output file dimension (%i / UNLIMITED)",
					dimA->name(), dimA->size());
			} else if (!dimA->is_unlimited() && !dimOut[d]->is_unlimited()) {
				_EXCEPTION3("Mismatch between input file dimension \"%s\" and "
					"output file dimension (%i / %i)",
					dimA->name(), dimA->size(), dimOut[d]->size());
			}
		}

		counts[d] = dimOut[d]->size();
		nDataSize *= counts[d];
	}

	// ncByte / ncChar type
	if ((var->type() == ncByte) || (var->type() == ncChar)) {
		DataArray1D<char> data(nDataSize);

		varOut =
			ncOut.add_var(
				var->name(), var->type(),
				dimOut.size(), (const NcDim**)&(dimOut[0]));

		if (varOut == NULL) {
			_EXCEPTION1("Cannot create variable \"%s\"", var->name());
		}

		if (fCopyAttributes) {
			CopyNcVarAttributes(var, varOut);
		}

		var->get(&(data[0]), &(counts[0]));
		varOut->put(&(data[0]), &(counts[0]));
	}

	// ncShort type
	if (var->type() == ncShort) {
		DataArray1D<short> data(nDataSize);

		varOut =
			ncOut.add_var(
				var->name(), var->type(),
				dimOut.size(), (const NcDim**)&(dimOut[0]));

		if (varOut == NULL) {
			_EXCEPTION1("Cannot create variable \"%s\"", var->name());
		}

		if (fCopyAttributes) {
			CopyNcVarAttributes(var, varOut);
		}

		var->get(&(data[0]), &(counts[0]));
		varOut->put(&(data[0]), &(counts[0]));
	}

	// ncInt type
	if (var->type() == ncInt) {
		DataArray1D<int> data(nDataSize);

		varOut =
			ncOut.add_var(
				var->name(), var->type(),
				dimOut.size(), (const NcDim**)&(dimOut[0]));

		if (varOut == NULL) {
			_EXCEPTION1("Cannot create variable \"%s\"", var->name());
		}

		if (fCopyAttributes) {
			CopyNcVarAttributes(var, varOut);
		}

		var->get(&(data[0]), &(counts[0]));
		varOut->put(&(data[0]), &(counts[0]));
	}

	// ncFloat type
	if (var->type() == ncFloat) {
		DataArray1D<float> data(nDataSize);

		varOut =
			ncOut.add_var(
				var->name(), var->type(),
				dimOut.size(), (const NcDim**)&(dimOut[0]));

		if (varOut == NULL) {
			_EXCEPTION1("Cannot create variable \"%s\"", var->name());
		}

		if (fCopyAttributes) {
			CopyNcVarAttributes(var, varOut);
		}

		var->get(&(data[0]), &(counts[0]));
		varOut->put(&(data[0]), &(counts[0]));
	}


	// ncDouble type
	if (var->type() == ncDouble) {
		DataArray1D<double> data(nDataSize);

		varOut =
			ncOut.add_var(
				var->name(), var->type(),
				dimOut.size(), (const NcDim**)&(dimOut[0]));

		if (varOut == NULL) {
			_EXCEPTION1("Cannot create variable \"%s\"", var->name());
		}

		if (fCopyAttributes) {
			CopyNcVarAttributes(var, varOut);
		}

		var->get(&(data[0]), &(counts[0]));
		varOut->put(&(data[0]), &(counts[0]));
	}

	// ncInt64 type
	if (var->type() == ncInt64) {
		DataArray1D<ncint64> data(nDataSize);

		varOut =
			ncOut.add_var(
				var->name(), var->type(),
				dimOut.size(), (const NcDim**)&(dimOut[0]));

		if (varOut == NULL) {
			_EXCEPTION1("Cannot create variable \"%s\"", var->name());
		}

		if (fCopyAttributes) {
			CopyNcVarAttributes(var, varOut);
		}

		var->get(&(data[0]), &(counts[0]));
		varOut->put(&(data[0]), &(counts[0]));
	}

	// ncUInt64 type
	if (var->type() == ncUInt64) {
		DataArray1D<ncuint64> data(nDataSize);

		varOut =
			ncOut.add_var(
				var->name(), var->type(),
				dimOut.size(), (const NcDim**)&(dimOut[0]));

		if (varOut == NULL) {
			_EXCEPTION1("Cannot create variable \"%s\"", var->name());
		}

		if (fCopyAttributes) {
			CopyNcVarAttributes(var, varOut);
		}

		var->get(&(data[0]), &(counts[0]));
		varOut->put(&(data[0]), &(counts[0]));
	}

	// Check output variable exists
	if (varOut == NULL) {
		_EXCEPTION1("Unable to create output variable \"%s\"",
			var->name());
	}
}

////////////////////////////////////////////////////////////////////////////////

void CopyNcVarIfExists(
	NcFile & ncIn,
	NcFile & ncOut,
	const std::string & strVarName,
	bool fCopyAttributes
) {
	// Turn off fatal errors in NetCDF
	NcError error(NcError::silent_nonfatal);

	NcVar * var = ncIn.get_var(strVarName.c_str());
	if (var != NULL) {
		CopyNcVar(ncIn, ncOut, strVarName, fCopyAttributes);
	}
}

////////////////////////////////////////////////////////////////////////////////


