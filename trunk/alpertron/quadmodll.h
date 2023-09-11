// This file is part of Alpertron Calculators.
//
// Copyright 2017-2021 Dario Alejandro Alpern
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

#ifndef _QUADMODLL_H
#define _QUADMODLL_H

#include "bignbr.h"
#include <functional>

using ShowSolutionsModPrimeFn = std::function<void(int factorIndex, int expon,
                                                   const BigInteger* pIncrement)>;
using SolutionFn = std::function<void(BigInteger *value)>;
using ShowNoSolsModPrimeFn = std::function<void(int expon)>;

struct QuadModLLResult {
  BigInteger Solution1[400];
  BigInteger Solution2[400];
  BigInteger Increment[400];

  BigInteger prime;
};

void SolveEquation(
    const SolutionFn &solutionCback,
    const ShowSolutionsModPrimeFn &showSolutionsModPrime,
    const ShowNoSolsModPrimeFn &showNoSolsModPrime,
    BigInteger *pValA, const BigInteger* pValB,
    const BigInteger* pValC, BigInteger* pValN,
    BigInteger *GcdAll, BigInteger *pValNn,
    QuadModLLResult *result);

#endif