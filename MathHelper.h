///////////////////////////////////////////////////////////////////////////////
///
///	\file    MathHelper.h
///	\author  Paul Ullrich
///	\version July 26, 2010
///
///	<summary>
///		This header file provides access to various mathematical functions
///		that are not normally exposed in the standard C++ libraries.
///	</summary>
///	<remarks>
///		Copyright 2000-2010 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#ifndef _MATHHELPER_H_
#define _MATHHELPER_H_

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		Calculate the maximum of two values.
///	</summary>
template <typename _T>
_T Max(_T x1, _T x2) {
	return (x1>x2)?(x1):(x2);
}

///	<summary>
///		Calculate the minimum of two values.
///	</summary>
template <typename _T>
_T Min(_T x1, _T x2) {
	return (x1<x2)?(x1):(x2);
}

///	<summary>
///		Calculate the sign of a value.
///	</summary>
template <typename _T>
_T Sign(_T x1) {
	return (x1 < static_cast<_T>(0))?(static_cast<_T>(-1)):(static_cast<_T>(1));
}

///	<summary>
///		Clamp a value to be within a given range.
///	</summary>
template <typename _T>
_T Clamp(_T y, _T x1, _T x2) {
	return (y>x2)?(x2):((y<x1)?(x1):(y));
}

///	<summary>
///		Calculate the integer square root.
///	</summary>
///	<remarks>
///		Source: Crenshaw, Jack.  Integer square roots.
///		http://www.embedded.com/98/9802fe2.htm
///	</remarks>
inline unsigned int ISqrt(unsigned int a) {
	unsigned int irem = 0;
	unsigned int iroot = 0;
	for (int i = 0; i < 16; i++) {
		iroot <<= 1;
		irem = ((irem << 2) + (a >> 30));
		a <<= 2;
		iroot++;
		if (iroot <= irem) {
			irem -= iroot;
			iroot++;
		} else {
			iroot--;
		}
	}
	return (static_cast<unsigned int>(iroot >> 1));
}

///	<summary>
///		Calculate the integer power of integers function.
///	</summary>
inline int IntPow(int d, unsigned int p) {
	if (p == 0) {
		return 1;
	}

	unsigned int q;

	int iPow = d;
	for (q = 1; q < p; q++) {
		iPow *= d;
	}
	return iPow;
}

///	<summary>
///		Calculate the integer power function.
///	</summary>
inline double IPow(double d, unsigned int p) {
	unsigned int q;

	double dPow = 1.0;

	for (q = 0; q < p; q++) {
		dPow *= d;
	}
	return dPow;
}

///	<summary>
///		Calculate the integer factorial function.
///	</summary>
inline unsigned int IFact(unsigned int p) {
	unsigned int q;
	unsigned int iFact = 1;

	for (q = 2; q <= p; q++) {
		iFact *= q;
	}

	return iFact;
}

///////////////////////////////////////////////////////////////////////////////

#endif
