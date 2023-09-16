// This file is part of Alpertron Calculators.
//
// Copyright 2015-2021 Dario Alejandro Alpern
//
// Alpertron Calculators is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Alpertron Calculators is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Alpertron Calculators.  If not, see <http://www.gnu.org/licenses/>.

#ifndef _BIGNBR_H
#define _BIGNBR_H

#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>

#include "enums.h"

// Tom multiplied these constants by 1000
// (It should really be MAX_LEN^2, right?)
#define MAX_LEN_MULT  25000000  // 200000 digits
// well, not this one. it doesn't compile with 2500000
#define MAX_LEN       25000  // 200000 digits
#define BITS_PER_GROUP         31
#define BITS_PER_GROUP_MINUS_1 30
#define HALF_INT_RANGE        0x40000000
#define HALF_INT_RANGE_U      0x40000000U
#define FOURTH_INT_RANGE      0x20000000
#define FOURTH_INT_RANGE_U    0x20000000U
#define MAX_INT_NBR           0x7FFFFFFF
#define MAX_INT_NBR_U         0x7FFFFFFFU
#define LIMB_RANGE            0x80000000U
#define MAX_VALUE_LIMB        0x7FFFFFFFU
#define SMALL_NUMBER_BOUND    32768

struct mylimb {
  int x;
};

typedef struct mylimb limb;

enum eSign {
  SIGN_POSITIVE = 0,
  SIGN_NEGATIVE,
};

enum eNbrCached {
  NBR_NOT_CACHED = 0,
  NBR_READY_TO_BE_CACHED,
  NBR_CACHED,
};

typedef struct BigInteger {
  int nbrLimbs;
  enum eSign sign;
  limb limbs[MAX_LEN];
} BigInteger;

extern limb MontgomeryMultR2[MAX_LEN];
extern limb MontgomeryMultR1[MAX_LEN];

void multiply(const limb* factor1, const limb* factor2, limb* result,
              int len, int* pResultLen);

void multiplyWithBothLen(const limb* factor1, const limb* factor2, limb* result,
                         int len1, int len2, int* pResultLen);

void squareRoot(const limb *argument, limb *sqRoot, int len, int *pLenSqRoot);
void AddBigInt(const limb *pAddend1, const limb *pAddend2, limb *pSum, int nbrLimbs);
bool BigIntIsZero(const BigInteger *value);
bool BigIntIsOne(const BigInteger* value);
bool BigIntEqual(const BigInteger *value1, const BigInteger *value2);
void BigIntChSign(BigInteger *value);
void BigIntAdd(const BigInteger *pAddend1, const BigInteger *pAddend2, BigInteger *pSum);
void BigIntSubt(const BigInteger* pMinuend, const BigInteger* pSubtrahend, BigInteger* pDifference);
void BigIntNegate(const BigInteger *pSrc, BigInteger *pDest);
enum eExprErr BigIntDivide(const BigInteger *pDividend, const BigInteger *pDivisor, BigInteger *pQuotient);
enum eExprErr BigIntMultiply(const BigInteger *pFactor1, const BigInteger *pFactor2, BigInteger *pProduct);
enum eExprErr BigIntRemainder(const BigInteger* pDividend,
  const BigInteger* pDivisor, BigInteger* pRemainder);
enum eExprErr BigIntPower(const BigInteger *pBase, const BigInteger *pExponent, BigInteger *pPower);
enum eExprErr BigIntPowerIntExp(const BigInteger *pBase, int exponent, BigInteger *pPower);

void floordiv(const BigInteger *num, const BigInteger *den, BigInteger *result);
enum eExprErr BigIntMultiplyPower2(BigInteger* pArg, int powerOf2);
void BigIntMultiplyBy2(BigInteger *nbr);
void BigIntDivideBy2(BigInteger *nbr);
void BigIntGcd(const BigInteger *pArg1, const BigInteger *pArg2, BigInteger *pResult);
void multint(BigInteger *pResult, const BigInteger *pMult, int factor);
void multadd(BigInteger *pResult, int iMult, const BigInteger *pMult, int addend);
void BigIntPowerOf2(BigInteger *pResult, int expon);
void subtractdivide(BigInteger *pBigInt, int subt, int divisor);
void addbigint(BigInteger *pResult, int addend);
void CopyBigInt(BigInteger *pDest, const BigInteger *pSrc);
int getNbrLimbs(const limb *bigNbr);
void BigIntDivide2(BigInteger *pArg);
void BigIntAnd(const BigInteger *firstArg, const BigInteger *secondArg, BigInteger *result);
void IntArray2BigInteger(
    int number_length, const int *ptrValues, /*@out@*/BigInteger *bigint);
void BigInteger2IntArray(int number_length,
                         /*@out@*/int *ptrValues, const BigInteger *bigint);
void UncompressLimbsBigInteger(int number_length,
                               const limb *ptrValues, /*@out@*/BigInteger *bigint);
void CompressLimbsBigInteger(int number_length,
                             /*@out@*/limb *ptrValues, const BigInteger *bigint);
void NbrToLimbs(int nbr, /*@out@*/limb *limbs, int len);

void intToBigInteger(BigInteger *bigint, int value);
double getMantissa(const limb *ptrLimb, int nbrLimbs);
void DivideBigNbrByMaxPowerOf2(int *pShRight, limb *number, int *pNbrLimbs);

void ChSignBigNbr(limb *nbr, int length);
void AddBigNbr(const limb *pNbr1, const limb *pNbr2, limb *pSum, int nbrLen);
void SubtractBigNbr(const limb *pNbr1, const limb *pNbr2, limb *pDiff, int nbrLen);
void DivBigNbrByInt(const limb *pDividend, int divisor, limb *pQuotient, int nbrLen);
int RemDivBigNbrByInt(const limb *pDividend, int divisor, int nbrLen);
void MultBigNbr(const limb *pFactor1, const limb *pFactor2, limb *pProd, int nbrLen);
void IntToBigNbr(int value, limb *bigNbr, int nbrLength);
int BigNbrToBigInt(const BigInteger *pBigNbr, limb *pBigInt);
void BigIntToBigNbr(BigInteger *pBigNbr, const limb *pBigInt, int nbrLenBigInt);
void AdjustBigIntModN(limb *Nbr, const limb *Mod, int nbrLen);
void IntToBigNbr(int value, limb *bigNbr, int nbrLength);
int BigIntJacobiSymbol(const BigInteger *upper, const BigInteger *lower);

static inline int UintToInt(unsigned int value) {
  return (int)value;
}
static inline int Uint64ToInt(uint64_t value) {
  return (int)value;
}
#endif
