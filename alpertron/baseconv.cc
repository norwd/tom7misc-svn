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
#include <stdint.h>
#include "bignbr.h"
#include "globals.h"
#include "karatsuba.h"

#define DIGITS_PER_LIMB 9
#define MAX_LIMB_CONVERSION 1000000000

static void add(const limb *addend1, const limb *addend2, limb *sum, int length);

  // Convert number to little-endian.
void Dec2Bin(const char *decimal, limb *binary, int digits, int *bitGroups) {
  limb power10000[MAX_LEN*2];
  limb temp[MAX_LEN];

  // First step: separate in groups of DIGITS_PER_LIMB digits.
  const char *ptrSrc;
  limb *ptrDest;
  limb *ptrBinary;
  int digit;
  int lenBytes;
  int multiplier;
  int nbrGroups = 1;
  while ((nbrGroups * DIGITS_PER_LIMB) < digits)
  {
    nbrGroups *= 2;
  }
  lenBytes = nbrGroups * (int)sizeof(limb);
  (void)memset(binary, 0, lenBytes);
  (void)memset(power10000, 0, lenBytes);
  power10000[0].x = MAX_LIMB_CONVERSION;
  ptrDest = binary;
  for (ptrSrc = decimal + digits - 1; ptrSrc >= (decimal + DIGITS_PER_LIMB-1); ptrSrc -= DIGITS_PER_LIMB)
  {
    int limbContents = 0;
    for (digit = DIGITS_PER_LIMB-1; digit >= 0; digit--)
    {
      limbContents = (limbContents * 10) + *(ptrSrc - digit) - '0';
    }
    ptrDest->x = limbContents;
    ptrDest++;
  }
  digit = 0;
  multiplier = 1;
  for (; ptrSrc >= decimal; ptrSrc--)
  {
    digit += multiplier * (*ptrSrc - '0');
    multiplier *= 10;
  }
  ptrDest->x = digit;
  for (int outerGroup = 1; outerGroup < nbrGroups; outerGroup += outerGroup)
  {
    for (int innerGroup = 0; innerGroup < nbrGroups; innerGroup += 2*outerGroup)
    {
      ptrBinary = binary + innerGroup;
      multiply(power10000, ptrBinary + outerGroup, temp, outerGroup, NULL);
      lenBytes = outerGroup * (int)sizeof(limb);
      (void)memset(ptrBinary + outerGroup, 0, lenBytes);
      add(temp, ptrBinary, ptrBinary, 2*outerGroup);
    }
    if ((outerGroup * 2) < nbrGroups)
    {    // Square power10000.
      multiply(power10000, power10000, temp, outerGroup, NULL);
      lenBytes = (outerGroup * 2) * (int)sizeof(limb);
      (void)memcpy(power10000, temp, lenBytes);
    }
  }
  // Determine first non-significant group.
  *bitGroups = 1;    // Initialize number of groups in advance.
  for (int groupNbr = nbrGroups-1; groupNbr > 0; groupNbr--)
  {
    if ((binary + groupNbr)->x != 0)
    {            // First non-significant group found.
      *bitGroups = groupNbr + 1;
      break;
    }
  }
}

void int2dec(char **pOutput, int nbr) {
  char *ptrOutput = *pOutput;
  bool significantZero = false;
  unsigned int div = 1000000000U;
  unsigned int value = (unsigned int)nbr;
  while (div > 0U)
  {
    int digit;

    digit = value/div;
    if ((digit > 0) || significantZero)
    {
      significantZero = true;
      *ptrOutput = (char)(digit + '0');
      ptrOutput++;
    }
    value %= div;
    div /= 10;
  }
  if (!significantZero)
  {
    *ptrOutput = '0';
    ptrOutput++;
  }
  *pOutput = ptrOutput;
}

static void Bin2DecLoop(char** ppDest, bool *pSignificantZero,
  int valueGrp, int grpLen, int *pGroupCtr, int *pDigits, bool last)
{
  int digit[DIGITS_PER_LIMB];
  int digits = *pDigits;
  int count;
  int value = valueGrp;
  char* ptrDest = *ppDest;
  int groupCtr = *pGroupCtr;
  bool significantZero = *pSignificantZero;
  for (count = 0; count < DIGITS_PER_LIMB; count++)
  {
    digit[count] = value % 10;
    value /= 10;
  }
  for (count = DIGITS_PER_LIMB - 1; count >= 0; count--)
  {
    if ((digit[count] != 0) || significantZero)
    {
      digits++;
      *ptrDest = (char)(digit[count] + '0');
      ptrDest++;
      if ((groupCtr == 1) && (!last || (count != 0)))
      {    // Do not insert space at the end of number.
        *ptrDest = ' ';
        ptrDest++;
      }
      significantZero = true;
    }
    groupCtr--;
    if (groupCtr == 0)
    {
      groupCtr = grpLen;
    }
  }
  *ppDest = ptrDest;
  *pGroupCtr = groupCtr;
  *pSignificantZero = significantZero;
  *pDigits = digits;
}

// Convert little-endian number to a string with space every groupLen digits.
  // In order to perform a faster conversion, use groups of DIGITS_PER_LIMB digits.
void Bin2Dec(char **ppDecimal, const limb *binary, int nbrLimbs, int groupLength)
{
  limb power10000[MAX_LEN*2];

  int grpLen = groupLength;
  int len;
  int index;
  int index2;
  const limb *ptrSrc = binary + nbrLimbs - 1;
  char *ptrDest;
  bool significantZero = false;
  int groupCtr;
  int digits=0;
  bool showDigitsText = true;
  unsigned int firstMult;
  unsigned int secondMult;
  firstMult = (unsigned int)BITS_PER_GROUP / 2U;
  firstMult = 1U << firstMult;
  secondMult = LIMB_RANGE / firstMult;

  if (grpLen <= 0)
  {
    grpLen = -grpLen;
    showDigitsText = false;
  }
  power10000[0].x = ptrSrc->x % MAX_LIMB_CONVERSION;
  power10000[1].x = ptrSrc->x / MAX_LIMB_CONVERSION;
  len = ((power10000[1].x == 0)? 1 : 2); // Initialize array length.
  for (index = nbrLimbs - 2; index >= 0; index--)
  {
    double dCarry;
    double dQuotient;
    limb *ptrPower;

    // Multiply by firstMult and then by secondMult, so there is never
    // more than 53 bits in the product.

    ptrPower = power10000;
    dQuotient = 0;
    for (index2 = 0; index2 < len; index2++)
    {
      dCarry = dQuotient + ((double)ptrPower->x * (double)firstMult);
      dQuotient = floor(dCarry / (double)MAX_LIMB_CONVERSION);
      ptrPower->x = (int)(dCarry - (dQuotient * (double)MAX_LIMB_CONVERSION));
      ptrPower++;
    }
    if (dQuotient != 0.0)
    {
      ptrPower->x = (int)dQuotient;
      len++;
    }
    ptrPower = power10000;
    ptrSrc--;
    dQuotient = ptrSrc->x;
    for (index2 = 0; index2 < len; index2++)
    {
      dCarry = dQuotient + ((double)ptrPower->x * (double)secondMult);
      dQuotient = floor(dCarry / (double)MAX_LIMB_CONVERSION);
      ptrPower->x = (int)(dCarry - (dQuotient * (double)MAX_LIMB_CONVERSION));
      ptrPower++;
    }
    if (dQuotient != 0.0)
    {
      ptrPower->x = (int)dQuotient;
      ptrPower++;
      len++;
    }
  }
  // At this moment the array power10000 has the representation
  // of the number in base 10000 in little-endian. Convert to
  // ASCII separating every groupLength characters.
  ptrDest = *ppDecimal;
  ptrSrc = &power10000[len-1];
  groupCtr = len * DIGITS_PER_LIMB;
  if (grpLen != 0)
  {
    groupCtr %= grpLen;
    if (groupCtr == 0)
    {
      groupCtr = grpLen;
    }
  }
  for (index = len; index > 0; index--)
  {
    int value = ptrSrc->x;
    ptrSrc--;
    Bin2DecLoop(&ptrDest, &significantZero, value, grpLen,
      &groupCtr, &digits, index == 1);
  }
  if (!significantZero)
  {     // Number is zero.
    *ptrDest = '0';
    ptrDest++;
  }
  else
  {
    if ((digits > 30) && showDigitsText)
    {
      *ptrDest = ' ';
      ptrDest++;
      *ptrDest = '(';
      ptrDest++;
      int2dec(&ptrDest, digits);
      copyStr(&ptrDest, (lang ? " dígitos)" : " digits)"));
    }
  }
  *ptrDest = '\0';             // Add terminator.
  *ppDecimal = ptrDest;
}

static void add(const limb *addend1, const limb *addend2, limb *sum, int length)
{
  const limb* ptrAddend1 = addend1;
  const limb* ptrAddend2 = addend2;
  limb* ptrSum = sum;
  unsigned int carry;
  carry = 0;
  for (int i = 0; i < length; i++)
  {
    carry += (unsigned int)ptrAddend1->x + (unsigned int)ptrAddend2->x;
    ptrAddend1++;
    ptrAddend2++;
    if (carry >= LIMB_RANGE)
    {
      ptrSum->x = (int)carry - (int)LIMB_RANGE;
      ptrSum++;
      carry = 1;
    }
    else
    {
      ptrSum->x = (int)carry;
      ptrSum++;
      carry = 0;
    }
  }
  return;
}

void BigInteger2Dec(char **ppDecimal, const BigInteger *pBigInt, int groupLength)
{
  char* ptrDecimal = *ppDecimal;
  if (pBigInt->sign == SIGN_NEGATIVE)
  {
    *ptrDecimal = '-';
    ptrDecimal++;
  }
  Bin2Dec(&ptrDecimal, pBigInt->limbs, pBigInt->nbrLimbs, groupLength);
  *ppDecimal = ptrDecimal;
}
