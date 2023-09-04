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

#ifndef _FACTOR_H
#define _FACTOR_H

#include "showtime.h"

#define MAX_FACTORS 5000

#define TYP_AURIF    100000000
#define TYP_TABLE    150000000
#define TYP_SIQS     200000000
#define TYP_LEHMAN   250000000
#define TYP_RABIN    300000000
#define TYP_DIVISION 350000000
#define TYP_EC       400000000

// I think the first entry is the residual (maybe composite?). The
// multiplicity is the number of distinct factors after that.
struct sFactors
{
  // I think that this is a length followed by limbs. -tom7
  // (See SortFactors)
  int *ptrFactor;
  int multiplicity;
  int upperBound;
  int type;
};

void factor(const BigInteger *toFactor, const int *number, int *factors,
            struct sFactors *pstFactors);

extern struct sFactors astFactorsMod[MAX_FACTORS];
extern int factorsMod[20000];

#endif
