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

#include <string.h>
#include <math.h>

#include "bignbr.h"

void ChSignBigNbr(limb *nbr, int length) {
  int carry = 0;
  const limb *ptrEndNbr = nbr + length;
  for (limb *ptrNbr = nbr; ptrNbr < ptrEndNbr; ptrNbr++)
  {
    carry -= ptrNbr->x;
    ptrNbr->x = carry & MAX_INT_NBR;
    carry >>= BITS_PER_GROUP;
  }
}

void AddBigNbr(const limb *pNbr1, const limb*pNbr2, limb*pSum, int nbrLen) {
  unsigned int carry = 0U;
  const limb*ptrNbr1 = pNbr1;
  const limb*ptrNbr2 = pNbr2;
  const limb*ptrEndSum = pSum + nbrLen;
  for (limb*ptrSum = pSum; ptrSum < ptrEndSum; ptrSum++)
  {
    unsigned int tmp;
    carry = (carry >> BITS_PER_GROUP) + (unsigned int)ptrNbr1->x +
      (unsigned int)ptrNbr2->x;
    tmp = carry & MAX_INT_NBR_U;
    ptrSum->x = (int)tmp;
    ptrNbr1++;
    ptrNbr2++;
  }
}

void SubtractBigNbr(const limb *pNbr1, const limb*pNbr2, limb*pDiff, int nbrLen)
{
  unsigned int borrow = 0U;
  const limb*ptrNbr1 = pNbr1;
  const limb*ptrNbr2 = pNbr2;
  const limb*ptrEndDiff = pDiff + nbrLen;
  for (limb*ptrDiff = pDiff; ptrDiff < ptrEndDiff; ptrDiff++)
  {
    unsigned int tmp;
    borrow = (unsigned int)ptrNbr1->x - (unsigned int)ptrNbr2->x -
      (borrow >> BITS_PER_GROUP);
    tmp = borrow & MAX_INT_NBR_U;
    ptrDiff->x = (int)tmp;
    ptrNbr1++;
    ptrNbr2++;
  }
}

void DivBigNbrByInt(const limb *pDividend, int divisor, limb *pQuotient, int nbrLen)
{
  const limb* ptrDividend = pDividend;
  limb* ptrQuotient = pQuotient;
  unsigned int remainder = 0U;
  double dDivisor = (double)divisor;
  double dLimb = 2147483648.0;
  int nbrLenMinus1 = nbrLen - 1;
  ptrDividend += nbrLenMinus1;
  ptrQuotient += nbrLenMinus1;
  for (int ctr = nbrLenMinus1; ctr >= 0; ctr--)
  {
    unsigned int dividend = (remainder << BITS_PER_GROUP) +
      (unsigned int)ptrDividend->x;
    double dDividend = ((double)remainder * dLimb) + (double)ptrDividend->x;
    // quotient has correct value or 1 more.
    unsigned int quotient = (unsigned int)((dDividend / dDivisor) + 0.5);
    remainder = dividend - (quotient * (unsigned int)divisor);
    if (remainder >= (unsigned int)divisor)
    {     // remainder not in range 0 <= remainder < divisor. Adjust.
      quotient--;
      remainder += (unsigned int)divisor;
    }
    ptrQuotient->x = (int)quotient;
    ptrQuotient--;
    ptrDividend--;
  }
}

void MultBigNbr(const limb *pFactor1, const limb *pFactor2, limb *pProd, int nbrLen)
{
  limb* ptrProd = pProd;
  double dRangeLimb = (double)(1U << BITS_PER_GROUP);
  double dInvRangeLimb = 1.0 / dRangeLimb;
  int low = 0;
  int factor1;
  int factor2;
  double dAccumulator = 0.0;
  for (int i = 0; i < nbrLen; i++)
  {
    for (int j = 0; j <= i; j++)
    {
      factor1 = (pFactor1 + j)->x;
      factor2 = (pFactor2 + i - j)->x;
      low += factor1*factor2;
      dAccumulator += (double)factor1 * (double)factor2;
    }
    low &= MAX_INT_NBR;    // Trim extra bits.
    ptrProd->x = low;
    ptrProd++;
    // Subtract or add 0x20000000 so the multiplication by dVal is not nearly an integer.
    // In that case, there would be an error of +/- 1.
    if (low < HALF_INT_RANGE)
    {
      dAccumulator = floor((dAccumulator + (double)FOURTH_INT_RANGE)*dInvRangeLimb);
    }
    else
    {
      dAccumulator = floor((dAccumulator - (double)FOURTH_INT_RANGE)*dInvRangeLimb);
    }
    low = (int)(dAccumulator - floor(dAccumulator * dInvRangeLimb) * dRangeLimb);
  }
  ptrProd->x = low;
  (ptrProd+1)->x = (int)floor(dAccumulator/dRangeLimb);
}

// On input: pFactor1 and pFactor2: pointers to factors.
//           pProd: pointer to product (length = 2*nbrLen)
//           nbrLen: number of limbs of factors.
void MultBigNbrComplete(const limb *pFactor1, const limb *pFactor2, limb *pProd, int nbrLen)
{
  limb* ptrProd = pProd;
  double dRangeLimb = (double)(1U << BITS_PER_GROUP);
  double dInvRangeLimb = 1.0 / dRangeLimb;
  int low = 0;
  int i;
  int j;
  int factor1;
  int factor2;
  double dAccumulator = 0.0;
  for (i = 0; i < nbrLen; i++)
  {
    for (j = 0; j <= i; j++)
    {
      factor1 = (pFactor1 + j)->x;
      factor2 = (pFactor2 + i - j)->x;
      low += factor1 * factor2;
      dAccumulator += (double)factor1 * (double)factor2;
    }
    low &= MAX_INT_NBR;    // Trim extra bits.
    ptrProd->x = low;
    ptrProd++;
    // Subtract or add 0x20000000 so the multiplication by dVal is not nearly an integer.
    // In that case, there would be an error of +/- 1.
    if (low < HALF_INT_RANGE)
    {
      dAccumulator = floor((dAccumulator + (double)FOURTH_INT_RANGE)*dInvRangeLimb);
    }
    else
    {
      dAccumulator = floor((dAccumulator - (double)FOURTH_INT_RANGE)*dInvRangeLimb);
    }
    low = (int)(dAccumulator - floor(dAccumulator * dInvRangeLimb) * dRangeLimb);
  }
  for (; i < (2*nbrLen); i++)
  {
    for (j = i-nbrLen+1; j < nbrLen; j++)
    {
      factor1 = (pFactor1 + j)->x;
      factor2 = (pFactor2 + i - j)->x;
      low += factor1 * factor2;
      dAccumulator += (double)factor1 * (double)factor2;
    }
    low &= MAX_INT_NBR;    // Trim extra bits.
    ptrProd->x = low;
    ptrProd++;
    // Subtract or add 0x20000000 so the multiplication by dVal is not nearly an integer.
    // In that case, there would be an error of +/- 1.
    if (low < HALF_INT_RANGE)
    {
      dAccumulator = floor((dAccumulator + (double)FOURTH_INT_RANGE)*dInvRangeLimb);
    }
    else
    {
      dAccumulator = floor((dAccumulator - (double)FOURTH_INT_RANGE)*dInvRangeLimb);
    }
    low = (int)(dAccumulator - floor(dAccumulator * dInvRangeLimb) * dRangeLimb);
  }
  ptrProd->x = low;
  (ptrProd + 1)->x = (int)floor(dAccumulator / dRangeLimb);
}
