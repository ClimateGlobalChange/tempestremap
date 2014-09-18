///////////////////////////////////////////////////////////////////////////////
///
///	\file    FixedPoint.h
///	\author  Paul Ullrich
///	\version August 6, 2014
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

#ifndef _FIXEDPOINT_H_
#define _FIXEDPOINT_H_

///////////////////////////////////////////////////////////////////////////////

#include <cstdint>
#include <vector>
#include <iostream>
#include <cstring>

#include "Defines.h"
#include "Exception.h"

///////////////////////////////////////////////////////////////////////////////

///	<summary>
///		A fixed point number, used for exact arithmetic.
///	</summary>
class FixedPoint {

public:
	///	<summary>
	///		Number of Digits to use.
	///	</summary>
	static const int Digits = 8;

	///	<summary>
	///		Maximum value in each Digit.
	///	</summary>
	static const int64_t MaximumDigit = 10000000000000000;

	///	<summary>
	///		Half maximum value in each Digit (log scale).
	///	</summary>
	static const int64_t HalfMaximumDigit = 100000000;

public:
	///	<summary>
	///		Default constructor.
	///	</summary>
	FixedPoint() :
		m_iSign(0),
		m_iDecimal(0)
	{
		Zero();
	}

	///	<summary>
	///		Constructor from a double.
	///	</summary>
	FixedPoint(double d) {
		Set(d);
	}

public:
	///	<summary>
	///		Normalize the FixedPoint number.
	///	</summary>
	inline void Normalize() {
		uint64_t nCarryover;
		for (int i = 0; i < Digits; i++) {
			nCarryover = m_vecDigits[i] / MaximumDigit;
			m_vecDigits[i] = m_vecDigits[i] % MaximumDigit;
			m_vecDigits[i+1] += nCarryover;
		}
		if (nCarryover != 0) {
			_EXCEPTIONT("FixedPoint overflow");
		}

		// Check for zero
		bool fIsZero = true;
		for (int i = 0; i < Digits; i++) {
			if (m_vecDigits[i] != 0) {
				fIsZero = false;
				break;
			}
		}
		if (fIsZero) {
			m_iSign = 0;
			m_iDecimal = 1;
		}
	}

public:
	///	<summary>
	///		Set equal to zero.
	///	</summary>
	inline void Zero() {
		m_iSign = 0;
		m_iDecimal = 1;
		memset(m_vecDigits, 0, sizeof(int64_t) * Digits);
	}

	///	<summary>
	///		Set the FixedPoint number.
	///	</summary>
	inline void Set(double d) {

		Zero();

		if (d > 0.0) {
			m_vecDigits[0] = static_cast<int64_t>( d * 1.0e16);
			m_iSign = +1;
		} else if (d == 0.0) {
			m_iSign = 0;
		} else {
			m_vecDigits[0] = static_cast<int64_t>(-d * 1.0e16);
			m_iSign = -1;
		}

		if (m_vecDigits[0] == MaximumDigit) {
			m_vecDigits[0] = MaximumDigit-1;
		} else if (m_vecDigits[0] > MaximumDigit) {
			_EXCEPTIONT("FixedPoint cannot be set by value larger than 1");
		}

		m_iDecimal = 1;

		Normalize();
	}

private:
	///	<summary>
	///		Member sum/difference operator
	///	</summary>
	inline void SumDifference(
		const FixedPoint & fp,
		bool fDifference
	) {
		// Check for zero
		if (fp.m_iSign == 0) {
			return;
		}

		// Shift decimal
		int iShift;
		if (m_iDecimal < fp.m_iDecimal) {
			int i;
			for (i = Digits-1; i >= fp.m_iDecimal - m_iDecimal; i--) {
				m_vecDigits[i] = m_vecDigits[i - fp.m_iDecimal + m_iDecimal];
			}
			for (; i >= 0; i--) {
				m_vecDigits[i] = 0;
			}
			m_iDecimal = fp.m_iDecimal;
			iShift = 0;

		} else {
			iShift = m_iDecimal - fp.m_iDecimal;
		}

		// Sign of the second term
		int iFPSign = fp.m_iSign;
		if (fDifference) {
			iFPSign = - iFPSign;
		}

		// Signs are the same, add values
		if (((m_iSign >= 0) && (iFPSign >= 0)) ||
			((m_iSign <= 0) && (iFPSign <= 0))
		) {
			for (int i = iShift; i < Digits; i++) {
				m_vecDigits[i] += fp.m_vecDigits[i - iShift];

				if (m_vecDigits[i] >= MaximumDigit) {
					if (i == Digits-1) {
						_EXCEPTIONT("Integer overflow detected");
					} else {
						m_vecDigits[i] -= MaximumDigit;
						m_vecDigits[i+1] += 1;
					}
				}
			}

			if (m_iSign == 0) {
				m_iSign = iFPSign;
			}

		// Signs are different, difference values
		} else {

			// Determine the first non-zero digit
			int iEnd;
			for (iEnd = Digits-1; iEnd >= iShift; iEnd--) {
				if ((m_vecDigits[iEnd] != 0) ||
					(fp.m_vecDigits[iEnd - iShift] != 0)
				) {
					break;
				}
			}

			// Both values are zero
			if (iEnd == iShift - 1) {
				_EXCEPTIONT("Logic error");
			}

			if (iShift != 0) {
				printf("Shift: %i\n", iShift);
			}

			// Determine which value is larger in magnitude
			bool fIAmLarger = true;
			for (int i = iEnd; i >= iShift; i--) {
				if (m_vecDigits[i] > fp.m_vecDigits[i - iShift]) {
					break;
				} else if (m_vecDigits[i] < fp.m_vecDigits[i - iShift]) {
					fIAmLarger = false;
					break;
				}
			}

			// This value is larger in magnitude
			if (fIAmLarger) {
				for (int i = iShift; i <= iEnd; i++) {
					m_vecDigits[i] -= fp.m_vecDigits[i - iShift];

					if (m_vecDigits[i] < 0) {
						m_vecDigits[i] += MaximumDigit;
						m_vecDigits[i+1] -= 1;
					}
				}

			// fp is larger in magnitude
			} else {
				m_iSign = - m_iSign;
				for (int i = iShift; i <= iEnd; i++) {
					m_vecDigits[i] = fp.m_vecDigits[i - iShift] - m_vecDigits[i];

					if (m_vecDigits[i] < 0) {
						m_vecDigits[i] += MaximumDigit;
						m_vecDigits[i+1] += 1;
					}
				}
			}
		}

		Normalize();
	}

	///	<summary>
	///		Product operator.
	///	</summary>
	inline void Product(
		const FixedPoint & fp
	) {
		FixedPoint fpTemp(*this);

		Zero();

		for (int i = 0; i < Digits/2; i++) {
			int64_t iS0 = fpTemp.m_vecDigits[i] % HalfMaximumDigit;
			int64_t iB0 = (fpTemp.m_vecDigits[i] - iS0) / HalfMaximumDigit;

			for (int j = 0; j < Digits/2; j++) {
				int64_t iS1 = fp.m_vecDigits[j] % HalfMaximumDigit;
				int64_t iB1 = (fp.m_vecDigits[j] - iS1) / HalfMaximumDigit;

				int64_t iSS = iS0 * iS1;
				int64_t iSB = iS0 * iB1;
				int64_t iBS = iB0 * iS1;
				int64_t iBB = iB0 * iB1;

				m_vecDigits[i+j  ] += iSS;
				m_vecDigits[i+j  ] +=
					((iSB + iBS) % HalfMaximumDigit) * HalfMaximumDigit;
				m_vecDigits[i+j+1] += (iSB + iBS) / HalfMaximumDigit;
				m_vecDigits[i+j+1] += iBB;
			}
		}

		m_iDecimal = fpTemp.m_iDecimal + fp.m_iDecimal;

		if (m_iDecimal > Digits) {
			Print(); printf("\n");
			fp.Print(); printf("\n");
			_EXCEPTIONT("FixedPoint overflow");
		}

		m_iSign = fpTemp.m_iSign * fp.m_iSign;

		Normalize();
	}

public:
	///	<summary>
	///		Assignment operator.
	///	</summary>
	inline const FixedPoint & operator=(const FixedPoint & fp) {
		m_iSign = fp.m_iSign;
		m_iDecimal = fp.m_iDecimal;
		memcpy(m_vecDigits, fp.m_vecDigits, Digits * sizeof(uint64_t));

		return (*this);
	}

	///	<summary>
	///		Inline sum operator.
	///	</summary>
	inline FixedPoint & operator+=(const FixedPoint & fp) {
		SumDifference(fp, false);
		return (*this);
	}

	///	<summary>
	///		Sum operator.
	///	</summary>
	inline FixedPoint operator+(const FixedPoint & fp) const {
		FixedPoint fpOut(*this);
		fpOut.SumDifference(fp, false);
		return fpOut;
	}

	///	<summary>
	///		Inline difference operator.
	///	</summary>
	inline FixedPoint & operator-=(const FixedPoint & fp) {
		SumDifference(fp, true);
		return (*this);
	}

	///	<summary>
	///		Difference operator.
	///	</summary>
	inline FixedPoint operator-(const FixedPoint & fp) const {
		FixedPoint fpOut(*this);
		fpOut.SumDifference(fp, true);
		return fpOut;
	}

	///	<summary>
	///		Inline product operator.
	///	</summary>
	inline FixedPoint & operator*=(const FixedPoint & fp) {
		Product(fp);
		return (*this);
	}

	///	<summary>
	///		Product operator.
	///	</summary>
	inline FixedPoint operator*(const FixedPoint & fp) const {
		FixedPoint fpOut(*this);
		fpOut.Product(fp);
		return fpOut;
	}

	///	<summary>
	///		Negation operator.
	///	</summary>
	inline const FixedPoint & Negate() {
		m_iSign = - m_iSign;
		return (*this);
	}

	///	<summary>
	///		Returns true if this number is positive.
	///	</summary>
	bool IsPositive() const {
		return (m_iSign > 0);
	}

	///	<summary>
	///		Returns true if this number is negative.
	///	</summary>
	bool IsNegative() const {
		return (m_iSign < 0);
	}

	///	<summary>
	///		Returns true if this number is nonpositive.
	///	</summary>
	bool IsNonPositive() const {
		return (m_iSign <= 0);
	}

	///	<summary>
	///		Returns true if this number is nonnegative.
	///	</summary>
	bool IsNonNegative() const {
		return (m_iSign >= 0);
	}

	///	<summary>
	///		Returns true if this number is zero.
	///	</summary>
	bool IsZero() const {
		return (m_iSign == 0);
	}

public:
	///	<summary>
	///		Convert this number to a double.
	///	</summary>
	Real ToReal() const {
		if (m_iSign == 0) {
			return 0.0;
		}
		if (m_iDecimal == 0) {
			_EXCEPTIONT("Invalid value of m_iDecimal");
		}

		const Real dInv = 1.0 / static_cast<double>(MaximumDigit);

		Real dCurrentInv = dInv;
		Real dOut = 0.0;
		for (int i = 1; i <= m_iDecimal; i++) {
			dOut += static_cast<Real>(m_vecDigits[m_iDecimal-i]) * dCurrentInv;
			dCurrentInv *= dInv;
		}

		if (m_iSign == -1) {
			dOut *= -1.0;
		}

		return dOut;
	}

	///	<summary>
	///		Print this number
	///	</summary>
	void Print() const {
		//printf("(%i) ", m_iDecimal);
		if (m_iSign < 0) {
			printf("-");
		}
		int i = Digits-1;
		for (; i > 0; i--) {
			if (m_vecDigits[i] != 0) {
				break;
			}
			if (i+1 == m_iDecimal) {
				break;
			}
		}
		for (; i >= 0; i--) {
			if (i+1 == m_iDecimal) {
				printf(".");
			}
			printf("%016llu", m_vecDigits[i]);
		}
	}

protected:
	///	<summary>
	///		Sign of this expression (-1, 0 or +1).
	///	</summary>
	int m_iSign;

	///	<summary>
	///		Decimal point location.
	///	</summary>
	int m_iDecimal;

	///	<summary>
	///		Array of 15-digit integers representing the floating point number.
	///	</summary>
	int64_t m_vecDigits[Digits];
};

///////////////////////////////////////////////////////////////////////////////

#endif

