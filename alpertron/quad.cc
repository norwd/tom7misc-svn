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

#include "quad.h"

#include <memory>
#include <vector>
#include <optional>
#include <utility>

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "factor.h"
#include "quadmodll.h"
#include "modmult.h"
#include "bigconv.h"

#include "base/stringprintf.h"
#include "base/logging.h"
#include "bignum/big.h"
#include "bignum/big-overloads.h"

using namespace std;

static constexpr bool VERBOSE = false;

namespace {

enum class LinSolType {
  // Solutions of the form (Xlin * t + Xind, Ylin * t + Yind).
  SOLUTION_FOUND,
  NO_SOLUTIONS,
  // (x, y) for any x,y.
  INFINITE_SOLUTIONS,
};

// Result of the linear equation solver; might not actually be
// a solution.
struct LinSol {
  LinSolType type = LinSolType::NO_SOLUTIONS;

  LinSol(LinSolType type) : type(type) {}

  // Only meaningful when SOLUTION_FOUND.
  BigInt Xlin, Xind;
  BigInt Ylin, Yind;

  void SwapXY() {
    std::swap(Xlin, Ylin);
    std::swap(Xind, Yind);
  }
};

enum class SolutionNumber {
  FIRST,
  SECOND,
};

enum class QmodCallbackType {
  PARABOLIC = 0,
  ELLIPTIC,
  HYPERBOLIC,
};

[[maybe_unused]]
static std::string CallbackString(QmodCallbackType t) {
  switch (t) {
  case QmodCallbackType::PARABOLIC: return "PARABOLIC";
  case QmodCallbackType::ELLIPTIC: return "ELLIPTIC";
  case QmodCallbackType::HYPERBOLIC: return "HYPERBOLIC";
  default: return "???";
  }
}

// This is just a display heuristic. Tests whether the
// BigInteger representation would be more than two limbs.
static bool IsBig(const BigInt &bg) {
  return BigIntNumLimbs(bg) > 2;
  // size_t s = 1 + (mpz_sizeinbase(bg.GetRep(), 2) / BITS_PER_GROUP);
  // fprintf(stderr, "%s has size %d in %d groups\n",
  // bg.ToString().c_str(), (int)s, BITS_PER_GROUP);
  // return s > 2;
}

static LinSol LinearEq(BigInt coeffX, BigInt coeffY, BigInt coeffInd) {
  if (VERBOSE) {
    printf("LinearEq %s %s %s\n",
           coeffX.ToString().c_str(),
           coeffY.ToString().c_str(),
           coeffInd.ToString().c_str());
  }
  // A linear equation. X + Y + IND = 0.

  if (coeffX == 0) {
    if (coeffY == 0) {
      if (coeffInd != 0) {
        return LinSol(LinSolType::NO_SOLUTIONS);
      } else {
        return LinSol(LinSolType::INFINITE_SOLUTIONS);
      }
    }

    if (!BigInt::DivisibleBy(coeffInd, coeffY)) {
      return LinSol(LinSolType::NO_SOLUTIONS);
    } else {
      LinSol sol(LinSolType::SOLUTION_FOUND);
      sol.Xind = BigInt(0);
      sol.Xlin = BigInt(1);
      sol.Yind = -BigInt::DivExact(coeffInd, coeffY);
      sol.Ylin = BigInt(0);
      return sol;
    }
  }

  if (coeffY == 0) {

    const auto [qq, rr] = BigInt::QuotRem(coeffInd, coeffX);

    if (rr != 0) {
      return LinSol(LinSolType::NO_SOLUTIONS);
    } else {
      LinSol sol(LinSolType::SOLUTION_FOUND);
      sol.Yind = BigInt(0);
      sol.Ylin = BigInt(1);
      sol.Xind = -qq;
      sol.Xlin = BigInt(0);
      return sol;
    }
  }

  const BigInt gcdxy = BigInt::GCD(coeffX, coeffY);

  if (gcdxy != 1) {
    // GCD is not 1.
    // To solve it, we first find the greatest common divisor of the
    // linear coefficients, that is: gcd(coeffX, coeffY) = gcdxy.

    if (!BigInt::DivisibleBy(coeffInd, gcdxy)) {
      // The independent coefficient is not a multiple of
      // the gcd, so there are no solutions.
      return LinSol(LinSolType::NO_SOLUTIONS);
    }

    // Divide all coefficients by the gcd.
    if (gcdxy != 0) {
      coeffX = BigInt::DivExact(coeffX, gcdxy);
      coeffY = BigInt::DivExact(coeffY, gcdxy);
      coeffInd = BigInt::DivExact(coeffInd, gcdxy);
    }
  }

  // Now the generalized Euclidean algorithm.
  // (BigInt may have an implementation of this?)

  BigInt U1(1);
  BigInt U2(0);
  BigInt U3 = coeffX;
  BigInt V1(0);
  BigInt V2(1);
  BigInt V3 = coeffY;

  BigInt q;

  while (V3 != 0) {

    if (VERBOSE) {
      printf("Euclid Step: %s %s %s %s %s %s\n",
             U1.ToString().c_str(),
             U2.ToString().c_str(),
             U3.ToString().c_str(),
             V1.ToString().c_str(),
             V2.ToString().c_str(),
             V3.ToString().c_str());
    }

    // q <- floor(U3 / V3).
    q = BigInt::DivFloor(U3, V3);

    {
      // T <- U1 - q * V1
      BigInt T = U1 - q * V1;
      U1 = std::move(V1);
      V1 = std::move(T);
    }

    {
      BigInt T = U2 - q * V2;
      U2 = std::move(V2);
      V2 = std::move(T);
    }

    {
      BigInt T = U3 - q * V3;
      U3 = std::move(V3);
      V3 = std::move(T);
    }
  }

  CHECK(U3 != 0);
  // Compute q as -coeffInd / U3
  q = -coeffInd / U3;

  // Compute Xind as -U1 * coeffInd / U3
  BigInt xind = U1 * q;

  BigInt xlin = coeffY;

  // Compute Yind as -U2 * coeffInd / U3
  BigInt yind = U2 * q;

  BigInt ylin = -coeffX;

  if (VERBOSE) {
    printf("Step: %s %s %s %s %s %s | %s %s %s %s\n",
           U1.ToString().c_str(),
           U2.ToString().c_str(),
           U3.ToString().c_str(),
           V1.ToString().c_str(),
           V2.ToString().c_str(),
           V3.ToString().c_str(),

           coeffX.ToString().c_str(),
           coeffY.ToString().c_str(),
           xind.ToString().c_str(),
           yind.ToString().c_str());
  }

  // Substitute variables so the independent coefficients can be minimized.
  // Reuse variables U1, U2, U3, V1, V2, V3.

  // U1 <- coeffX^2 + coeffY^2
  U1 = coeffX * coeffX + coeffY * coeffY;

  // U2 <- (coeffX^2 + coeffY^2)/2
  U2 = U1 >> 1;

  U2 += coeffX * yind;
  U2 -= coeffY * xind;

  // U1 <- delta to add to t'
  U1 = BigInt::DivFloor(U2, U1);

  if (VERBOSE) {
    printf("After subst: %s %s %s %s %s %s\n",
           U1.ToString().c_str(),
           U2.ToString().c_str(),
           U3.ToString().c_str(),
           V1.ToString().c_str(),
           V2.ToString().c_str(),
           V3.ToString().c_str());
  }

  // Xind <- Xind + coeffY * delta
  q = U1 * coeffY;
  xind += q;

  // Yind <- Yind - coeffX * delta
  q = U1 * coeffX;
  yind -= q;

  if (xlin < 0 && ylin < 0) {
    // If both coefficients are negative, make them positive.
    // printf("negate_coeff coverage\n");
    xlin = -xlin;
    ylin = -ylin;
  }

  LinSol sol(LinSolType::SOLUTION_FOUND);

  sol.Xlin = std::move(xlin);
  sol.Xind = std::move(xind);
  sol.Ylin = std::move(ylin);
  sol.Yind = std::move(yind);

  return sol;
}

// Output:
// nullopt: There are no solutions because gcd(P, Q, R) > 1
// some(P, Q, R) with gcd(P, Q, R) = 1.
std::optional<std::tuple<BigInt, BigInt, BigInt>>
PerformTransformation(
    const BigInt &A, const BigInt &B, const BigInt &C, const BigInt &K,
    const BigInt &Value) {
  // writes: P, Q, R, H, I

  const BigInt VA = A * Value;

  // Compute P as (at^2+bt+c)/K
  const BigInt P = ((VA + B) * Value + C) / K;

  // Compute Q <- -(2at + b).
  const BigInt Q = -((VA << 1) + B);

  // Compute R <- aK
  const BigInt R = A * K;

  // Compute gcd of P, Q and R.

  // Note: Used to write H and I as temporaries, but I think they're dead.
  const BigInt I = BigInt::GCD(BigInt::GCD(P, Q), R);
  if (I == 1) {
    // Gcd equals 1.
    return {std::make_tuple(P, Q, R)};
  }

  // No solutions because gcd(P, Q, R) > 1.
  return std::nullopt;
}

// Accumulates the minimum solution in xdst, ydst. Would be cleaner
// if we passed around some kind of explicit state.
void AccumulateXY(const BigInt &X, const BigInt &Y,
                  std::optional<std::pair<BigInt, BigInt>> *sol) {

  CHECK(sol != nullptr);

  if (!sol->has_value()) {
    sol->emplace(X, Y);
  } else {
    const auto &[BX, BY] = sol->value();
    // Take the smallest by the sum of absolute values.

    if (BigInt::Abs(X) + BigInt::Abs(Y) <=
        BigInt::Abs(BX) + BigInt::Abs(BY)) {
      sol->emplace(X, Y);
    }
  }
}

// Returns true if solution found.
bool AccumulatePoint(
    const BigInt &X, const BigInt &Y,
    const BigInt &Alpha, const BigInt &Beta,
    const BigInt &Div,
    std::optional<std::pair<BigInt, BigInt>> *sol) {

  CHECK(sol != nullptr);

  // (I think this should actually be impossible because Div comes from
  // the GCD of the coefficients.)
  CHECK(Div != 0) << "Might be shenanigans with divisibility by zero";

  // Check first that (X+alpha) and (Y+beta) are multiple of D.
  BigInt tmp1 = X + Alpha;
  if (BigInt::DivisibleBy(tmp1, Div)) {
    BigInt tmp2 = Y + Beta;
    if (BigInt::DivisibleBy(tmp2, Div)) {

      tmp1 = BigInt::DivExact(tmp1, Div);
      tmp2 = BigInt::DivExact(tmp2, Div);

      // HYPERBOLIC.
      AccumulateXY(tmp1, tmp2, sol);
      return true;
    }
  }

  return false;
}


// First check: |u| < g.
// Second check: |u+g| > |v|
// Third check: |u-g| < |v|
// On input G = floor(g), g > 0.
// g is not an integer number.
static bool CheckStartOfContinuedFractionPeriod(const BigInt &U,
                                                const BigInt &V,
                                                const BigInt &G) {
  if (G >= BigInt::Abs(U)) {
    // First check |u| < g passed.
    // Set Tmp1 to |v|
    BigInt Tmp1 = BigInt::Abs(V);
    // Compute Tmp2 as u + floor(g) which equals floor(u+g)
    BigInt Tmp2 = U + G;

    if (Tmp2 < 0) {
      // Round to number nearer to zero.
      Tmp2 += 1;
      // addbigint(&Tmp2, 1);
    }

    Tmp2 = BigInt::Abs(Tmp2);

    // Compute Tmp2 as floor(|u+g|)
    // Compute bigTmp as floor(|u+g|) - |v|
    if (Tmp2 >= Tmp1) {
      // Second check |u+g| > |v| passed.
      // Compute Tmp2 as u - floor(g)
      Tmp2 = U - G;

      if (Tmp2 <= 0) {
        // Round down number to integer.
        Tmp2 -= 1;
      }

      Tmp2 = BigInt::Abs(Tmp2);

      // Compute Tmp2 as floor(|u-g|)
      // Compute Tmp2 as |v| - floor(|u-g|)
      if (Tmp1 >= Tmp2) {
        // Third check |u-g| < |v| passed.
        // Save U and V to check period end.
        return true;
      }
    }
  }
  return false;
}

// Returns Temp0, Temp1
static std::pair<BigInt, BigInt>
UnimodularSubstitution(const BigInt &M,
                       const BigInt &Z,
                       const BigInt &O) {
  BigInt Temp0, Temp1;
  if (M < 0) {
    // Perform the substitution: x = X + Y, y = (|m|-1)X + |m|Y
    Temp0 = (Z + O);
    Temp1 = Temp0 * -M - Z;
  } else if (M == 0) {
    Temp0 = Z;
    Temp1 = O;
  } else {
    // Perform the substitution: x = mX + (m-1)Y, y = X + Y
    Temp1 = Z + O;
    Temp0 = Temp1 * M - O;
  }
  return std::make_pair(Temp0, Temp1);
}


struct Quad {
  // Solutions accumulated here.
  Solutions solutions;
  // True if we had a solution on a hyperbolic curve and so we should
  // print recursive solutions as well.
  bool hyperbolic_recursive_solution = false;
  // This is essentially debugging output now.
  std::string *output = nullptr;

  void ShowText(const std::string &text) {
    if (output != nullptr)
      *output += text;
  }

  inline void ShowChar(char c) {
    if (output != nullptr)
      output->push_back(c);
  }

  void showMinus() {
    if (output != nullptr)
      *output += " - ";
  }

  void ShowBigInt(const BigInt &value) {
    if (output != nullptr) {
      *output += value.ToString();
    }
  }

  void showInt(int value) {
    if (output != nullptr) {
      StringAppendF(output, "%d", value);
    }
  }

  void SolutionHeader() {
    ShowText("\nSolution:");
  }

  void ShowLin(const BigInt &coeffX, const BigInt &coeffY,
               const BigInt &coeffInd,
               const char *x, const char *y) {
    LinSolType t;
    t = Show(coeffX, x, LinSolType::SOLUTION_FOUND);
    t = Show(coeffY, y, t);
    Show1(coeffInd, t);
  }

  void ShowLinInd(const BigInt &lin, const BigInt &ind) {
    if (ind == 0 && lin == 0) {
      ShowText("0");
    }
    if (ind != 0) {
      ShowBigInt(ind);
    }
    ShowChar(' ');

    if (lin < 0) {
      showMinus();
    } else if (lin != 0 && ind != 0) {
      ShowChar('+');
    } else {
      // Nothing to do.
    }
    ShowChar(' ');
    if (lin != 0) {
      if (BigInt::Abs(lin) != 1) {
        // abs(lin) is not 1
        // CopyBigInt(&Aux0, lin);
        // Do not show negative sign twice.
        // Aux0.sign = SIGN_POSITIVE;
        ShowBigInt(BigInt::Abs(lin));
      }
      ShowText(" t");
    }
  }


  void PrintLinear(const LinSol &sol) {
    switch (sol.type) {
    default:
    case LinSolType::NO_SOLUTIONS:
      return;

    case LinSolType::INFINITE_SOLUTIONS:
      SolutionHeader();
      solutions.any_integers = true;
      ShowText("\nx, y: any integer");
      return;

    case LinSolType::SOLUTION_FOUND:
      // Port note: This used to actually have the effect of swapping
      // xind/yind xlin/ylin.
      SolutionHeader();
      ShowText("\nx = ");
      ShowLinInd(sol.Xlin, sol.Xind);
      ShowText("\ny = ");
      ShowLinInd(sol.Ylin, sol.Yind);
      LinearSolution s;
      s.MX = sol.Xlin;
      s.BX = sol.Xind;
      s.MY = sol.Ylin;
      s.BY = sol.Yind;
      solutions.linear.emplace_back(std::move(s));
      return;
    }
  }

  void ShowSolutionXY(const BigInt &x, const BigInt &y) {
    ShowText("\nx = ");
    ShowBigInt(x);
    ShowText("\ny = ");
    ShowBigInt(y);
    solutions.points.emplace_back(PointSolution{.X = x, .Y = y});
  }

  void PrintQuad(const BigInt &T2, const BigInt &T,
                 const BigInt &Ind) {
    const char *var1 = "k";
    if (BigInt::Abs(T2) == 1) {
      // abs(coeffT2) = 1
      if (T2 == 1) {
        // coeffT2 = 1
        ShowChar(' ');
      } else {
        // coeffT2 = -1
        showMinus();
      }
      ShowText(var1);
      ShowText("^2");
    } else if (T2 != 0) {
      // coeffT2 is not zero.
      ShowBigInt(T2);
      ShowChar(' ');
      ShowText(var1);
      ShowText("^2");
    } else {
      // Nothing to do.
    }

    if (T < 0) {
      ShowText(" - ");
    } else if (T != 0 && T2 != 0) {
      ShowText(" + ");
    } else {
      // Nothing to do.
    }

    if (BigInt::Abs(T) == 1) {
      // abs(coeffT) = 1
      ShowText(var1);
      ShowText(" * ");
    } else if (T != 0) {
      // Port note: original called showlimbs if negative, which I
      // think is just a way of printing the absolute value without
      // any copying.
      ShowBigInt(BigInt::Abs(T));
      ShowText(" ");
      ShowText(var1);
    } else {
      // Nothing to do.
    }

    if (Ind != 0) {
      if (T != 0 || T2 != 0) {
        if (Ind < 0) {
          ShowText(" - ");
        } else {
          ShowText(" + ");
        }
      } else if (Ind < 0) {
        ShowText(" -");
      } else {
        // Nothing to do.
      }

      // Same trick for abs value.
      ShowBigInt(BigInt::Abs(Ind));
    }
  }


  // XXX why does this take/return "linear solution type" ?
  LinSolType Show(const BigInt &num, const string &str,
                  LinSolType t) {
    LinSolType tOut = t;
    if (num != 0) {
      // num is not zero.
      if (t == LinSolType::NO_SOLUTIONS && num >= 0) {
        ShowText(" +");
      }

      if (num < 0) {
        ShowText(" -");
      }

      if (num != 1 && num != -1) {
        // num is not 1 or -1.
        ShowChar(' ');
        ShowBigInt(BigInt::Abs(num));
        ShowText(" * ");
      } else {
        ShowText(" ");
      }

      if (output != nullptr)
        *output += str;

      if (t == LinSolType::SOLUTION_FOUND) {
        tOut = LinSolType::NO_SOLUTIONS;
      }
    }
    return tOut;
  }

  void Show1(const BigInt &num, LinSolType t) {
    const LinSolType u = Show(num, "", t);
    ShowChar(' ');
    // Port note: This used to test u & 1 as a "trick" for detecting
    // NO_SOLUTIONS?
    if (u != LinSolType::NO_SOLUTIONS || num == 1 || num == -1) {
      // Show absolute value of num.
      ShowBigInt(BigInt::Abs(num));
    }
  }

  // Print the original quadratic equation.
  void ShowEq(
      const BigInt &coeffA, const BigInt &coeffB,
      const BigInt &coeffC, const BigInt &coeffD,
      const BigInt &coeffE, const BigInt &coeffF,
      const char *x, const char *y) {

    LinSolType t;
    string vxx = StringPrintf("%s^2", x);
    t = Show(coeffA, vxx, LinSolType::SOLUTION_FOUND);

    string vxy = StringPrintf("%s * %s", x, y);
    t = Show(coeffB, vxy, t);

    string vyy = StringPrintf("%s^2", y);
    t = Show(coeffC, vyy, t);

    t = Show(coeffD, x, t);

    t = Show(coeffE, y, t);

    if (coeffF < 0) {
      ShowText(" - ");
      ShowBigInt(BigInt::Abs(coeffF));
    } else {
      ShowText(" + ");
      ShowBigInt(coeffF);
    }
  }

  void ShowRecSol(char variable,
                  const BigInt &cx,
                  const BigInt &cy,
                  const BigInt &ci) {
    ShowChar(variable);
    ShowText("_n+1 = ");
    LinSolType t = Show(cx, "x_n",
                        LinSolType::SOLUTION_FOUND);
    t = Show(cy, "y_n", t);
    Show1(ci, t);
  }

  void ShowResult(const char *text, const BigInt &value) {
    ShowText(text);
    ShowText(" = ");
    ShowBigInt(value);
    ShowText("\n");
  }


  // Compute coefficients of x: V3 * w^2 + V2 * w + V1
  // Returns V3, V2, V1
  std::tuple<BigInt, BigInt, BigInt> ComputeXDiscrZero(
      const BigInt &A, const BigInt &B,
      const BigInt &C, const BigInt &D,
      const BigInt &E, const BigInt &Z,
      const BigInt &J, const BigInt &K,
      const BigInt &U2) {
    // Let m = 2be - 4cd
    // U3 <- m
    BigInt U3 = (B * E - ((C * D) << 1)) << 1;
    // Compute V1 <- (x'j - k)/2a + mx'^2
    BigInt V1 = (((U2 * J) - K) / A) >> 1;
    V1 += U3 * U2 * U2;
    // Compute V2 <- (j/2a + 2mx')z
    BigInt V2 = (J / A) >> 1;
    V2 += ((U3 * U2) << 1);
    V2 *= Z;
    // Compute V3 as m*z^2
    BigInt V3 = U3 * Z * Z;
    return std::make_tuple(V3, V2, V1);
  }

  // Compute coefficients of y: V3 * w^2 + V2 * w + V1
  // Returns V3, V2, V1
  std::tuple<BigInt, BigInt, BigInt> ComputeYDiscrZero(
      const BigInt &U, const BigInt &U2,
      const BigInt &S, const BigInt &R,
      const BigInt &Z) {

    // Compute V1 <- r + sx' + ux'^2
    BigInt V1 = (U * U2 + S) * U2 + R;

    // Compute V2 <- (s + 2ux')z
    BigInt V2 = (((U * U2) << 1) + S) * Z;

    // Compute V3 <- uz^2
    BigInt V3 = U * Z * Z;

    return std::make_tuple(V3, V2, V1);
  }


  // Only the parabolic mode will swap x and y.
  void CallbackQuadModParabolic(
      bool swap_xy,
      const BigInt &A, const BigInt &B, const BigInt &C,
      const BigInt &D, const BigInt &E,
      const BigInt &U, const BigInt &V,
      const BigInt &Value) {
    // The argument of this function is T. t = T - d + uk (k arbitrary).
    // Compute R <- (T^2 - v)/u

    BigInt R = ((Value * Value) - V) / U;
    // Compute S as 2*T
    BigInt S = Value << 1;

    // Find k from the congruence
    //  jk = K (mod z) where j = u-bs, K = d+br-T, z = 2a
    // Compute j <- u-bs
    BigInt J = U - B * S;
    // Compute K <- d+br-T
    BigInt K = (D + B * R) - Value;
    // Compute z <- 2a
    BigInt Z = A << 1;
    // If K is not multiple of gcd(j, z) there is no solution.
    BigInt BigTmp = BigInt::GCD(J, Z);
    CHECK(BigTmp != 0);
    if (!BigInt::DivisibleBy(K, BigTmp)) {
      return;
    }

    // Compute g = gcd(j, K, z), then recalculate j <- j/g, K <- K/g, z <- z/g
    BigInt U1 = BigInt::GCD(BigTmp, K);
    CHECK(U1 != 0);

    // Divide by J, K, Z by their GCD.
    // U2 <- j
    BigInt U2 = BigInt::DivExact(J, U1);
    // U3 <- K
    BigInt U3 = BigInt::DivExact(K, U1);
    // Use positive sign for modulus.
    Z = BigInt::Abs(BigInt::DivExact(Z, U1));

    if (Z != 0) U2 %= Z;
    // PERF: Can just use Mod?
    if (U2 < 0) U2 += Z;

    if (Z != 0) U3 %= Z;
    if (U3 < 0) U3 += Z;

    if (U2 == 0) {
      // M and N equal zero.
      // In this case 0*k = 0 (mod z) means any k is valid.
      Z = BigInt(1);
    } else {
      // U2 <- x'
      printf("GeneralModularDivision(%s,%s,%s) coverage\n",
              U2.ToString().c_str(), U3.ToString().c_str(),
              Z.ToString().c_str());
      solutions.interesting_coverage = true;
      // XXX This might be wrong for some inputs? See test.
      U2 = GeneralModularDivision(U2, U3, Z);
    }


    QuadraticSolution qsol;

    {
      auto coeff_x = ComputeXDiscrZero(A, B, C, D, E, Z, J, K, U2);
      auto coeff_y = ComputeYDiscrZero(U, U2, S, R, Z);

      if (swap_xy) {
        std::tie(qsol.VX, qsol.MX, qsol.BX) = std::move(coeff_y);
        std::tie(qsol.VY, qsol.MY, qsol.BY) = std::move(coeff_x);
      } else {
        std::tie(qsol.VX, qsol.MX, qsol.BX) = std::move(coeff_x);
        std::tie(qsol.VY, qsol.MY, qsol.BY) = std::move(coeff_y);
      }
    }

    SolutionHeader();

    // Result box:
    ShowText("\nx = ");
    PrintQuad(qsol.VX, qsol.MX, qsol.BX);

    ShowText("\ny = ");
    PrintQuad(qsol.VY, qsol.MY, qsol.BY);

    solutions.quadratic.emplace_back(std::move(qsol));
  }

  // Obtain next convergent of continued fraction of U/V
  // Previous convergents U1/V1, U2/V2, U3/V3.
  std::tuple<BigInt, BigInt, BigInt,
             BigInt, BigInt, BigInt> GetNextConvergent(
                 BigInt U, BigInt U1, BigInt U2,
                 BigInt V, BigInt V1, BigInt V2) {
    BigInt Tmp = BigInt::DivFloor(U, V);

    // Compute new value of U and V.
    BigInt Tmp2 = U - Tmp * V;
    U = std::move(V);
    V = Tmp2;

    // Compute new convergents: h_n = a_n*h_{n-1} + h_{n-2}
    // and also k_n = k_n*k_{n-1} + k_{n-2}
    BigInt Tmp3 = Tmp * U1 + U2;

    BigInt U3 = std::move(U2);
    U2 = std::move(U1);
    U1 = Tmp3;

    Tmp *= V1;
    Tmp += V2;
    BigInt V3 = std::move(V2);
    V2 = std::move(V1);
    V1 = Tmp;
    return std::make_tuple(U, U1, U2,
                           V, V1, V2);
  }

  void ShowAllRecSols(
      BigInt P, BigInt Q,
      BigInt R, BigInt S,
      BigInt K, BigInt L,
      const BigInt &Alpha, const BigInt &Beta) {

    RecursiveSolution rsol1;
    rsol1.P = P;
    rsol1.Q = Q;
    rsol1.R = R;
    rsol1.S = S;
    rsol1.K = K;
    rsol1.L = L;

    if (IsBig(P) || IsBig(Q)) {
      if (Alpha == 0 && Beta == 0) {
        ShowText("x_n+1 = P * x_n + Q * y_n\n"
                 "y_n+1 = R * x_n + S * y_n\n");
      } else {
        ShowText("x_n+1 = P * x_n + Q * y_n + K\n"
                 "y_n+1 = R * x_n + S * y_n + L\n");
      }
      ShowText("where:\n");
      ShowResult("P", P);
      ShowResult("Q", Q);
      if (Alpha != 0 || Beta != 0) {
        ShowResult("K", K);
      }
      ShowResult("R", R);
      ShowResult("S", S);
      if (Alpha != 0 || Beta != 0) {
        ShowResult("L", L);
      }
    } else {
      ShowText("\n");
      ShowRecSol('x', P, Q, K);
      ShowText("\n");
      ShowRecSol('y', R, S, L);
    }

    // Compute the related solution. Port note: The original
    // code made these modifications in place.

    // Compute x_{n-1} from x_n and y_n
    // Compute new value of K and L as:
    //   Knew <- L*Q - K*S and Lnew <- K*R - L*P
    BigInt Tmp1 = L * Q - K * S;
    L = K * R - L * P;
    K = std::move(Tmp1);

    // Compute new values of P, Q, R and S as:
    // Pnew <- S, Qnew <- -Q, Rnew <- -R, Snew <- P
    Q = -std::move(Q);
    R = -std::move(R);

    BigInt Tmp = P;
    P = S;
    S = std::move(Tmp);

    ShowText("\n... and also:\n");
    if (IsBig(P) || IsBig(Q)) {
      ShowResult("P", P);
      ShowResult("Q", Q);
      if (Alpha != 0 || Beta != 0) {
        ShowResult("K", K);
      }
      ShowResult("R", R);
      ShowResult("S", S);
      if (Alpha != 0 || Beta != 0) {
        ShowResult("L", L);
      }
    } else {
      ShowRecSol('x', P, Q, K);
      ShowText("\n");
      ShowRecSol('y', R, S, L);
    }

    RecursiveSolution rsol2;
    rsol2.P = std::move(P);
    rsol2.Q = std::move(Q);
    rsol2.R = std::move(R);
    rsol2.S = std::move(S);
    rsol2.K = std::move(K);
    rsol2.L = std::move(L);

    solutions.recursive.emplace_back(std::move(rsol1),
                                     std::move(rsol2));
  }

  bool SolutionFoundFromContFraction(bool isBeven,
                                     int V,
                                     const BigInt &Alpha,
                                     const BigInt &Beta,
                                     const BigInt &A,
                                     const BigInt &B,
                                     const BigInt &C,
                                     const BigInt &Discr,
                                     const BigInt &U1,
                                     const BigInt &V1) {
    BigInt P, S;
    if (isBeven) {
      // P <- r - (b/2)s
      // S <- r + (b/2)s
      BigInt tmp = (B >> 1) * V1;
      P = U1 - tmp;
      S = U1 + tmp;
    } else {
      // P <- r - bs
      // S <- r + bs
      BigInt tmp = B * V1;
      P = U1 - tmp;
      S = U1 + tmp;
      if (V == 4) {
        // P <- (r - bs)/2
        // S <- (r + bs)/2
        P >>= 1;
        S >>= 1;
      }
    }

    // Q <- -cs
    BigInt Q = -(C * V1);
    // R <- as
    BigInt R = A * V1;

    if (!isBeven && V == 1) {
      // Q <- -2cs
      Q <<= 1;
      // R <- 2as
      R <<= 1;
    }

    BigInt K = (Alpha * P) + (Beta * Q);
    BigInt L = (Alpha * R) + (Beta * S);

    if (VERBOSE) {
      printf("contf: %s %s %s %s | %s %s %s %s | %s %s | %s %s\n",
             A.ToString().c_str(),
             B.ToString().c_str(),
             C.ToString().c_str(),
             Discr.ToString().c_str(),
             P.ToString().c_str(),
             Q.ToString().c_str(),
             R.ToString().c_str(),
             S.ToString().c_str(),
             Alpha.ToString().c_str(),
             Beta.ToString().c_str(),
             K.ToString().c_str(),
             L.ToString().c_str());
    }

    CHECK(Discr != 0) << "Original code may have had shenanigans "
      "with dividing by zero";

    // Check whether alpha - K and beta - L are multiple of discriminant.
    const BigInt AlphaMinusK = Alpha - K;
    const BigInt BetaMinusL = Beta - L;
    if (BigInt::DivisibleBy(AlphaMinusK, Discr) &&
        BigInt::DivisibleBy(BetaMinusL, Discr)) {
      // Solution found.
      K = BigInt::DivExact(AlphaMinusK, Discr);
      L = BigInt::DivExact(BetaMinusL, Discr);
      ShowAllRecSols(P, Q, R, S,
                     K, L, Alpha, Beta);
      return true;
    }

    // Check whether alpha + K and beta + L are multiple of discriminant.
    const BigInt AlphaPlusK = Alpha + K;
    const BigInt BetaPlusL = Beta + L;
    if (BigInt::DivisibleBy(AlphaPlusK, Discr) &&
        BigInt::DivisibleBy(BetaPlusL, Discr)) {
      // Solution found.
      K = BigInt::DivExact(AlphaPlusK, Discr);
      L = BigInt::DivExact(BetaPlusL, Discr);

      ShowAllRecSols(-P, -Q, -R, -S,
                     K, L, Alpha, Beta);
      return true;
    }
    return false;
  }

  void ShowXYOne(const BigInt &X, const BigInt &Y) {
    SolutionHeader();
    ShowSolutionXY(X, Y);
  }


  // Use continued fraction of sqrt(B^2-4AC)
  // If the discriminant is 5, the method does not work: use 3, 1 and 7, 3.
  // If the convergent is r/s we get:
  // x(n+1) = Px(n) + Qy(n) + K
  // y(n+1) = Rx(n) + Sy(n) + L
  // where if b is odd:
  //        P = (r - bs)/2, Q = -cs, R = as, S = (r + bs)/2,
  // if b is even:
  //        P = r - (b/2)s, Q = -cs, R = as, S = r + (b/2)s,
  // in any case:
  //        K = (alpha*(1-P) - beta*Q) / D, L = (-alpha*R + beta*(1-S)) / D.
  void GetRecursiveSolution(
      BigInt A, BigInt B, BigInt C,
      const BigInt &ABack, const BigInt &BBack, const BigInt &CBack,
      const BigInt &Alpha, const BigInt &Beta,
      const BigInt &GcdHomog, BigInt Discr) {

    BigInt H = Discr;

    const bool isBeven = B.IsEven();
    if (isBeven) {
      H >>= 2;
    }

    // Obtain original discriminant.
    Discr *= GcdHomog;
    Discr *= GcdHomog;

    std::optional<int64_t> gcdo = GcdHomog.ToInt();
    CHECK(gcdo.has_value()) << "Original code seems to assume this, "
      "accessing the first limb directly.";
    const int64_t gcd_homog = gcdo.value();


    if (Discr == 5) {
      // Discriminant is 5.
      // Do not use continued fraction because it does not work.
      ShowText("\nRecursive solutions (discr 5):\n");

      // 3,1 is first solution to U1^2 - 5*V1^2 = 4
      if (SolutionFoundFromContFraction(isBeven, 4,
                                        Alpha, Beta,
                                        A, B, C,
                                        Discr,
                                        BigInt(3),
                                        BigInt(1))) {
        return;
      }

      // 9,4 is first solution to U1^2 - 5*V1^2 = 1
      (void)SolutionFoundFromContFraction(isBeven, 1,
                                          Alpha, Beta,
                                          A, B, C,
                                          Discr,
                                          BigInt(9),
                                          BigInt(4));
      return;
    }

    // g <- sqrt(discr).
    BigInt G = BigInt::Sqrt(H);
    // Port note: Was explicit SIGN_POSITIVE here in original, but I
    // think that was just because it was manipulating the limbs
    // directly? Sqrt is always non-negative...
    CHECK(G >= 0);

    int periodLength = 1;

    BigInt U(0);
    BigInt V(1);

    BigInt UU2(0);
    BigInt UU1(1);

    BigInt VV2(1);
    BigInt VV1(0);

    BigInt UBak, VBak;

    if (gcd_homog != 1) {
      periodLength = -1;
      do {
        // fprintf(stderr, "rec_gcdhomognotone coverage H=%s\n",
        // H.ToString().c_str());
        BigInt BigTmp = U + G;
        if (V < 0) {
          // If denominator is negative, round square root upwards.
          BigTmp += 1;
        }

        // Tmp1 = Term of continued fraction.
        BigInt Tmp1 = BigInt::DivFloor(BigTmp, V);

        // U <- a*V - U
        U = Tmp1 * V - U;

        // Port note: In alpertron, this uses L, which is used as
        // the discriminant (or discriminant/4) in ContFrac and may
        // have the same value here. But I think this was a bug; H
        // is discr (or discr/4) in this one. Passing around L
        // rather than setting a global results it in being
        // uninitialized here. Unfortunately, not good coverage of
        // this loop.

        // V <- (D - U^2)/V
        V = (H - U * U) / V;

        if (periodLength < 0) {
          UBak = U;
          VBak = V;
        }
        periodLength++;
      } while (periodLength == 1 || U != UBak || V != VBak);
      // Reset values of U and V.
      U = BigInt{0};
      V = BigInt{1};
    }

    if (periodLength > 1) {
      // quad.exe 6301 1575 2 7199 -1 -114995928
      printf("nonzeroperiod coverage (%d)\n", periodLength);
      solutions.interesting_coverage = true;
    }

    ShowText("\nRecursive solutions:\n");

    int periodNbr = 0;
    enum eSign sign = SIGN_POSITIVE;
    for (;;) {
      BigInt BigTmp = U + G;
      if (V < 0) {
        // If denominator is negative, round square root upwards.
        BigTmp += 1;
      }
      // Tmp1 = Term of continued fraction.
      BigInt Tmp1 = BigInt::DivFloor(BigTmp, V);

      // U3 <- U2, U2 <- U1, U1 <- a*U2 + U3
      BigInt UU3 = UU2;
      UU2 = UU1;
      UU1 = Tmp1 * UU2 + UU3;

      // V3 <- V2, V2 <- V1, V1 <- a*V2 + V3
      BigInt VV3 = VV2;
      VV2 = VV1;
      VV1 = Tmp1 * VV2 + VV3;

      U = Tmp1 * V - U;
      V = (H - U * U) / V;

      if (sign == SIGN_POSITIVE) {
        sign = SIGN_NEGATIVE;
      } else {
        sign = SIGN_POSITIVE;
      }

      if (VERBOSE)
      printf("FS: %c %s %s %s %d\n",
             isBeven ? 'e' : 'o',
             V.ToString().c_str(),
             Alpha.ToString().c_str(),
             Beta.ToString().c_str(),
             periodNbr);

      // V must have the correct sign.
      if ((sign == SIGN_NEGATIVE) ? V >= 0 : V < 0) {
        continue;
      }

      // Expecting denominator to be 1 (B even or odd)
      // or 4 (B odd) with correct sign.
      if (BigInt::Abs(V) != 1 &&
          (isBeven || BigInt::Abs(V) != 4)) {
        continue;
      }

      periodNbr++;
      if ((periodNbr * periodLength) % gcd_homog != 0) {
        continue;
      }


      // Found solution from continued fraction.
      if (SolutionFoundFromContFraction(isBeven,
                                        BigInt::Abs(V).ToInt().value(),
                                        Alpha, Beta,
                                        ABack, BBack, CBack,
                                        Discr,
                                        UU1, VV1)) {
        return;
      }
    }
  }


  // On input: H: value of u, I: value of v.
  // Output: ((tu - nv)*E, u*E) and ((-tu + nv)*E, -u*E)
  // If m is greater than zero, perform the substitution:
  //    x = mX + (m-1)Y, y = X + Y
  // If m is less than zero, perform the substitution:
  //    x = X + Y, y = (|m|-1)X + |m|Y
  // Do not substitute if m equals zero.
  // Returns true if solution found.
  bool NonSquareDiscrSolutionOne(
      const BigInt &M, const BigInt &E, const BigInt &K,
      const BigInt &Alpha, const BigInt &Beta, const BigInt &Div,
      const BigInt &H, const BigInt &I,
      const BigInt &Value) {

    // Port note: This used to modify the value of K based on the
    // callback type, but now we do that at the call site. (Also there
    // was something suspicious in here where it flipped the sign and
    // then set it negative.)

    // X = (tu - Kv)*E
    const BigInt Z = (Value * H - K * I) * E;
    // Y = u*E
    const BigInt O = H * E;

    // (we get here with both values for two_solutions)

    bool sol_found = false;
    // Undo unimodular substitution
    {
      const auto &[Temp0, Temp1] =
        UnimodularSubstitution(M, Z, O);
      sol_found = ShowPointOne(Temp0, Temp1, Alpha, Beta, Div);
    }

    // Z: (-tu - Kv)*E
    // O: -u*E

    // Undo unimodular substitution
    {
      const auto &[Temp0, Temp1] =
        UnimodularSubstitution(M, -Z, -O);
      sol_found = ShowPointOne(Temp0, Temp1, Alpha, Beta, Div) ||
        sol_found;
    }
    return sol_found;
  }

  // Same, when there are two solutions.
  // Returns true if a solution found (solFound = true in original code).
  // XXX this should return them, rather than modifying xminus/yminus
  bool NonSquareDiscrSolutionTwo(
      const BigInt &M, const BigInt &E, const BigInt &K,
      const BigInt &Alpha, const BigInt &Beta, const BigInt &Div,
      const BigInt &H, const BigInt &I,
      const BigInt &Value,
      std::optional<std::pair<BigInt, BigInt>> *sol_plus,
      std::optional<std::pair<BigInt, BigInt>> *sol_minus) {

    CHECK(sol_plus != nullptr);
    CHECK(sol_minus != nullptr);

    // Port note: This used to modify the value of K based on the callback
    // type, but now we do that at the call site. (Also there was something
    // suspicious in here where it flipped the sign and then set it negative.)

    // X = (tu - Kv)*E
    const BigInt Z = (Value * H - K * I) * E;
    // Y = u*E
    const BigInt O = H * E;


    bool sol_found = false;
    // Undo unimodular substitution
    {
      const auto &[Temp0, Temp1] =
        UnimodularSubstitution(M, Z, O);
      sol_found = AccumulatePoint(Temp0, Temp1,
                                  Alpha, Beta, Div,
                                  sol_plus);
    }

    // Z: (-tu - Kv)*E
    // O: -u*E

    // Undo unimodular substitution
    {
      const auto &[Temp0, Temp1] =
        UnimodularSubstitution(M, -Z, -O);
      sol_found = AccumulatePoint(Temp0, Temp1,
                                  Alpha, Beta, Div,
                                  sol_minus) || sol_found;
    }

    return sol_found;
  }


  // Returns true if we found a solution (and so should show
  // the recursive solution).
  bool ShowPointOne(const BigInt &X, const BigInt &Y,
                    const BigInt &Alpha, const BigInt &Beta,
                    const BigInt &Div) {

    // Check first that (X+alpha) and (Y+beta) are multiple of D.
    BigInt tmp1 = X + Alpha;
    BigInt tmp2 = Y + Beta;

    // (I think this should actually be impossible because Div comes from
    // the GCD of the coefficients.)
    CHECK(Div != 0) << "Might be shenanigans with divisibility by zero";

    if (BigInt::DivisibleBy(tmp1, Div) &&
        BigInt::DivisibleBy(tmp2, Div)) {

      if (Div != 0) {
        tmp1 = BigInt::DivExact(tmp1, Div);
        tmp2 = BigInt::DivExact(tmp2, Div);
      }

      // Not HYPERBOLIC.
      // Result box:
      ShowXYOne(tmp1, tmp2);
      return true;
    }
    return false;
  }


  void CheckSolutionSquareDiscr(
      const BigInt &CurrentFactor,
      const BigInt &H, const BigInt &I, const BigInt &L,
      const BigInt &M, const BigInt &Z,
      const BigInt &Alpha, const BigInt &Beta, const BigInt &Div) {

    CHECK(CurrentFactor != 0);
    BigInt N = Z / CurrentFactor;

    // (IL - HM)X = NI - cM
    // (IL - HM)Y = cL - NH

    // O = Denominator.
    BigInt O = I * L - H * M;

    // P = Numerator of X.
    BigInt P = N * I - CurrentFactor * M;

    if (VERBOSE) {
      printf("CheckSolutionSquareDiscr %s %s %s %s %s %s\n",
             P.ToString().c_str(),
             O.ToString().c_str(),
             CurrentFactor.ToString().c_str(),
             L.ToString().c_str(),
             N.ToString().c_str(),
             H.ToString().c_str());
    }

    CHECK(O != 0) << "Might have been shenanigans with O = 0?";
    if (BigInt::DivisibleBy(P, O)) {
      // X found.
      BigInt U1 = BigInt::DivExact(P, O);
      // ValP = Numerator of Y.
      P = CurrentFactor * L - N * H;

      CHECK(O != 0);
      if (BigInt::DivisibleBy(P, O)) {
        // Y found.
        BigInt U2 = BigInt::DivExact(P, O);
        // Show results.

        ShowPointOne(U1, U2, Alpha, Beta, Div);
        return;
      }
    }

    // The system of two equations does not have integer solutions.
    // No solution found.
  }

  // TODO: Try to make this dispatch (callbackQuadModType) static.
  template<QmodCallbackType QMOD_CALLBACK>
  void SolutionX(bool origin_translated, bool swap_xy,
                 BigInt Value, const BigInt &Modulus,
                 const BigInt &A, const BigInt &B, const BigInt &C,
                 const BigInt &D, const BigInt &E,
                 const BigInt &M, const BigInt &K,
                 const BigInt &U, const BigInt &V,
                 const BigInt &Alpha, const BigInt &Beta,
                 const BigInt &Div, const BigInt &Discr) {

    if (VERBOSE) {
      printf("SolutionX(%s, %s)\n",
             Value.ToString().c_str(),
             Modulus.ToString().c_str());
    }

    // If 2*value is greater than modulus, subtract modulus.
    if ((Value << 1) > Modulus) {
      Value -= Modulus;
    }

    if (VERBOSE) {
      printf("  with %s %s %s %s %s | %s %s | %s %s\n",
             A.ToString().c_str(),
             B.ToString().c_str(),
             C.ToString().c_str(),
             D.ToString().c_str(),
             E.ToString().c_str(),
             M.ToString().c_str(),
             K.ToString().c_str(),
             U.ToString().c_str(),
             V.ToString().c_str());
    }

    switch (QMOD_CALLBACK) {
    case QmodCallbackType::PARABOLIC:
      // Real uses of swap_xy
      CallbackQuadModParabolic(swap_xy,
                               A, B, C, D, E,
                               U, V, Value);
      break;

    case QmodCallbackType::ELLIPTIC:
      CallbackQuadModElliptic(A, B, C, E, M, K,
                              Alpha, Beta, Div, Discr,
                              Value);
      break;

    case QmodCallbackType::HYPERBOLIC:
      CallbackQuadModHyperbolic(origin_translated,
                                A, B, C, K, E, M, Alpha, Beta, Div,
                                Discr, Value);
      break;

    default:
      break;
    }
  }

  // Solve congruence an^2 + bn + c = 0 (mod n) where n is different from zero.
  template<QmodCallbackType QMOD_CALLBACK>
  void SolveQuadModEquation(
      bool origin_translated,
      bool swap_xy,
      const BigInt &coeffQuadr,
      const BigInt &coeffLinear,
      const BigInt &coeffIndep,
      BigInt Modulus,
      const BigInt &A, const BigInt &B, const BigInt &C, const BigInt &D, const BigInt &E,
      const BigInt &M, const BigInt &K, const BigInt &U, const BigInt &V,
      const BigInt &Alpha, const BigInt &Beta, const BigInt &Div, const BigInt &Discr) {

    if (VERBOSE) {
      printf("[SQME] %s %s %s %s\n",
             coeffQuadr.ToString().c_str(),
             coeffLinear.ToString().c_str(),
             coeffIndep.ToString().c_str(),
             Modulus.ToString().c_str());
    }

    CHECK(Modulus > 0);

    BigInt coeff_quadr = BigInt::CMod(coeffQuadr, Modulus);
    if (coeff_quadr < 0) coeff_quadr += Modulus;

    BigInt coeff_linear = BigInt::CMod(coeffLinear, Modulus);
    if (coeff_linear < 0) coeff_linear += Modulus;

    BigInt coeff_indep = BigInt::CMod(coeffIndep, Modulus);
    if (coeff_indep < 0) coeff_indep += Modulus;

    BigInt GcdAll = BigInt::GCD(coeff_indep,
                                BigInt::GCD(coeff_quadr, coeff_linear));

    // For a GCD of zero here, original code would cause and ignore
    // a division by zero, then read 0 from the temporary.

    if (GcdAll != 0 && !BigInt::DivisibleBy(coeff_indep, GcdAll)) {
      // C must be multiple of gcd(A, B).
      // Otherwise go out because there are no solutions.
      return;
    }

    GcdAll = BigInt::GCD(Modulus, GcdAll);

    // Divide all coefficients by gcd(A, B).
    if (GcdAll != 0) {
      coeff_quadr = BigInt::DivExact(coeff_quadr, GcdAll);
      coeff_linear = BigInt::DivExact(coeff_linear, GcdAll);
      coeff_indep = BigInt::DivExact(coeff_indep, GcdAll);
      Modulus = BigInt::DivExact(Modulus, GcdAll);
    }

    // PERF just substitute this through
    const BigInt &ValNn = Modulus;

    if (ValNn == 1) {
      // All values from 0 to GcdAll - 1 are solutions.
      if (GcdAll > 5) {
        printf("allvalues coverage\n");
        solutions.interesting_coverage = true;

        // XXX should we call SolutionX here?
        // Seems like we would want to do so until we find
        // a solution, at least?

        ShowText("\nAll values of x between 0 and ");

        // XXX Suspicious that this modifies GcdAll in place (I
        // think just to display it?) but uses it again below.
        GcdAll -= 1;
        ShowBigInt(GcdAll);
        ShowText(" are solutions.");
      } else {
        // must succeed; is < 5 and non-negative

        const int n = GcdAll.ToInt().value();
        for (int ctr = 0; ctr < n; ctr++) {
          SolutionX<QMOD_CALLBACK>(origin_translated,
                                   swap_xy,
                                   BigInt(ctr), Modulus,
                                   A, B, C, D, E,
                                   M, K,
                                   U, V,
                                   Alpha, Beta, Div, Discr);
        }
      }
      return;
    }

    if (BigInt::DivisibleBy(coeff_quadr, Modulus)) {
      // Linear equation.
      printf("linear-eq coverage\n");
      solutions.interesting_coverage = true;

      if (BigInt::GCD(coeff_linear, Modulus) != 1) {
        // B and N are not coprime. Go out.
        return;
      }

      // Calculate z <- -C / B (mod N)

      // Is it worth it to convert to montgomery form for one division??
      const std::unique_ptr<MontgomeryParams> params =
        GetMontgomeryParams(Modulus);

      BigInt z =
        BigIntModularDivision(*params, coeff_indep, coeff_linear, Modulus);

      if (z != 0) {
        // not covered by cov.sh :(
        printf("new coverage z != 0\n");
        solutions.interesting_coverage = true;

        // XXX is this a typo for ValNn in the original?
        // ValN is only set in CheckSolutionSquareDiscr.
        // Since it was static, it would usually be zero here.
        // z = ValN - z;
        z = 0 - z;
      }
      BigInt Temp0 = ValNn * GcdAll;

      for (;;) {
        // also not covered :(
        printf("new coverage: loop zz");
        solutions.interesting_coverage = true;
        SolutionX<QMOD_CALLBACK>(origin_translated,
                                 swap_xy,
                                 z, Modulus,
                                 A, B, C, D, E,
                                 M, K,
                                 U, V,
                                 Alpha, Beta, Div, Discr);
        z += Modulus;
        if (z < Temp0) break;
      }

      return;
    }

    if (QMOD_CALLBACK == QmodCallbackType::PARABOLIC) {
      // For elliptic case, the factorization is already done.

      // To solve this quadratic modular equation we have to
      // factor the modulus and find the solution modulo the powers
      // of the prime factors. Then we combine them by using the
      // Chinese Remainder Theorem.
    }

    if (VERBOSE) {
      printf("[Call SolveEq] %s %s %s %s %s %s\n",
             coeff_quadr.ToString().c_str(),
             coeff_linear.ToString().c_str(),
             coeff_indep.ToString().c_str(),
             Modulus.ToString().c_str(),
             GcdAll.ToString().c_str(),
             ValNn.ToString().c_str());
    }

    bool interesting = false;
    SolveEquation(
        SolutionFn([&](const BigInt &Value) {
            this->SolutionX<QMOD_CALLBACK>(
                origin_translated,
                swap_xy,
                Value,
                Modulus,
                A, B, C, D, E,
                M, K,
                U, V,
                Alpha, Beta, Div, Discr);
          }),
        coeff_quadr, coeff_linear, coeff_indep,
        // XXX ValNn arg is always the same as Modulus
        Modulus, GcdAll, ValNn,
        &interesting);
    if (interesting) {
      printf("INTERESTING!\n");
      solutions.interesting_coverage = true;
    }
  }

  void DiscriminantIsZero(BigInt A, BigInt B, BigInt C,
                          BigInt D, BigInt E, BigInt F) {
    // fprintf(stderr, "disciszero coverage\n");
    // Next algorithm does not work if A = 0. In this case, exchange x and y.
    bool swap_xy = false;
    if (A == 0) {
      swap_xy = true;
      // Exchange coefficients of x^2 and y^2.
      std::swap(A, C);
      // Exchange coefficients of x and y.
      std::swap(D, E);
    }

    // ax^2 + bxy + cx^2 + dx + ey + f = 0 (1)
    // Multiplying by 4a:
    // (2ax + by)^2 + 4adx + 4aey + 4af = 0
    // Let t = 2ax + by. So (1) becomes: (t + d)^2 = uy + v.

    // Compute u <- 2(bd - 2ae)
    BigInt U = (B * D - ((A * E) << 1)) << 1;

    // Compute v <- d^2 - 4af
    BigInt V = D * D - ((A * F) << 2);

    if (U == 0) {
      // u equals zero, so (t+d)^2 = v.

      if (V < 0) {
        // There are no solutions when v is negative,
        // since a square cannot be equal to a negative number.
        return;
      }

      if (V == 0) {
        // printf("disczero_vzero coverage\n");
        // v equals zero, so (1) becomes 2ax + by + d = 0
        LinSol sol = LinearEq(A << 1, B, D);
        // Result box:
        if (swap_xy) sol.SwapXY();
        PrintLinear(sol);
        return;
      }

      // u equals zero but v does not.
      // v must be a perfect square, otherwise there are no solutions.
      BigInt G = BigInt::Sqrt(V);
      if (V != G * G) {
        // v is not perfect square, so there are no solutions.
        return;
      }

      // The original equation is now: 2ax + by + (d +/- g) = 0
      BigInt A2 = A << 1;

      // This equation represents two parallel lines.
      {
        LinSol sol = LinearEq(A2, B, D + G);
        if (swap_xy) sol.SwapXY();
        // Result box:
        PrintLinear(sol);
      }

      {
        LinSol sol = LinearEq(A2, B, D - G);
        if (swap_xy) sol.SwapXY();
        PrintLinear(sol);
      }

      return;
    }

    // At this moment u does not equal zero.
    // We have to solve the congruence
    //     T^2 = v (mod u) where T = t+d and t = 2ax+by.

    // These were actually uninitialized on this code path,
    // and are probably unused.
    const BigInt M(0);
    const BigInt K(0);

    SolveQuadModEquation<QmodCallbackType::PARABOLIC>(
        // Origin never translated on this path.
        false,
        swap_xy,
        // Coefficients and modulus
        BigInt(1), BigInt(0), -V, BigInt::Abs(U),
        // Problem state
        A, B, C, D, E,
        M, K, U, V,
        // I think these end up unused in this case, but anyway,
        // don't translate the origin.
        BigInt(0), BigInt(0), BigInt(1),
        // Discriminant is known to be zero in this function.
        BigInt(0));
  }

  // Solve ax^2+bxy+cy^2 = K
  // The quadratic modular equation algorithm requires that gcd(a, n) = 1.
  // At this point gcd(a, b, c) = 1
  // The possibilities are:
  // - gcd(a, K) = 1. There is nothing to do.
  // Otherwise perform the transformation
  //    x = PX + QY, y = RX + SY with PS - QR = 1.
  // In particular perform: x = mX + (m-1)Y, y = X + Y
  // We get:
  //    (am^2+bm+c)*X^2 + (2(am^2+bm+c) - (2am+b))*XY +
  //    ((am^2+bm+c) - (2am+b) + a)*Y^2
  // Also perform: x = X + Y, y = (m-1)X + mY
  // We get:
  //    (a+(m-1)*b+(m-1)^2*c)*X^2 + (2a + b(2m-1) + 2cm(m-1))*X*Y +
  //    (a+bm+cm^2)*Y^2
  // The discriminant of the new formula does not change.
  // Compute m=1, 2, 3,... until gcd(am^2+bm+c, K) = 1.
  // When using the second formula, change sign of m so we know the
  // formula used when undoing the unimodular transformation later.

  // Since the algorithm discovers only primitive solutions,
  // i.e. solutions (x,y) where gcd(x,y) = 1, we need to solve
  //     ax'^2+bx'y'+cy'^2 = K/R^2 where R^2 is a divisor of K.
  // Then we get x = Rx', y = Ry'.
  template<QmodCallbackType QMOD_CALLBACK>
  void NonSquareDiscriminant(bool origin_translated,
                             BigInt A, BigInt B, BigInt C,
                             BigInt K,
                             const BigInt &D,
                             BigInt Discr,
                             BigInt Alpha, BigInt Beta, const BigInt &Div) {

    // These were actually uninitialized, and probably unused?
    const BigInt U(0);
    const BigInt V(0);

    // Find GCD(a,b,c)
    BigInt GcdHomog = BigInt::GCD(BigInt::GCD(A, B), C);
    // Divide A, B, C and K by this GCD.
    if (GcdHomog != 0) {
      A = BigInt::DivExact(A, GcdHomog);
      B = BigInt::DivExact(B, GcdHomog);
      C = BigInt::DivExact(C, GcdHomog);
      K = BigInt::DivExact(K, GcdHomog);
      // Divide discriminant by the square of GCD.
      Discr /= GcdHomog;
      Discr /= GcdHomog;
    }

    if (K == 0) {
      // If k=0, the only solution is (X, Y) = (0, 0)
      (void)ShowPointOne(BigInt(0), BigInt(0), Alpha, Beta, Div);
      return;
    }

    if (VERBOSE) {
      printf("start NSD %s %s %s | %s %s | %s %s %s\n",
             A.ToString().c_str(), B.ToString().c_str(), C.ToString().c_str(),
             K.ToString().c_str(), Discr.ToString().c_str(),
             Alpha.ToString().c_str(), Beta.ToString().c_str(),
             Div.ToString().c_str());
    }

    // ughhh
    BigInt ABack = A;
    BigInt BBack = B;
    BigInt CBack = C;

    // Factor independent term.

    // Note that we modify the factors (multiplicities) in place below.
    std::vector<std::pair<BigInt, int>> factors = BigIntFactor(BigInt::Abs(K));

    if (VERBOSE) {
      for (const auto &[f, m] : factors) {
        printf("%s^%d * ", f.ToString().c_str(), m);
      }
      printf("\n");
    }

    // PERF! The whole point of this is to maintain the list of factors
    // so that we can pass them in, and not have to re-factor. See
    // alpertron64.

    // Find all indices of prime factors with even multiplicity.
    // (XXX parallel. could be pair)
    // Index of prime factors with even multiplicity
    // PORT NOTE: was 1-based in original code; now 0-based
    std::vector<int> indexEvenMultiplicity, originalMultiplicities;
    const int numFactors = factors.size();
    for (int i = 0; i < numFactors; i++) {
      const auto &[fact, multiplicity] = factors[i];
      if (multiplicity > 1) {
        // At least prime is squared.
        // Port note: The original code stored factorNbr, which was 1-based
        // because of the factor header.
        indexEvenMultiplicity.push_back(i);
        // Convert to even.
        originalMultiplicities.push_back(multiplicity & ~1);
      }
    }

    // XXX base it on the size of factors
    std::vector<int> counters(400, 0);
    std::vector<bool> is_descending(400, false);

    BigInt E = BigInt(1);
    // Loop that cycles through all square divisors of the independent term.
    BigInt M(0);
    if (BigInt::GCD(A, K) != 1) {
      // gcd(a, K) is not equal to 1.

      BigInt UU1, UU2;
      do {
        // printf("uu1uu2 loop\n");

        // Compute U1 = cm^2 + bm + a and exit loop if this
        // value is not coprime to K.

        UU2 = C * M;
        UU1 = (UU2 + B) * M + A;

        if (VERBOSE) {
          printf("%s GCD %s = %s\n",
                 UU1.ToString().c_str(),
                 K.ToString().c_str(),
                 BigInt::GCD(UU1, K).ToString().c_str());
        }

        if (BigInt::GCD(UU1, K) == 1) {
          // Increment M and change sign to indicate type.
          M = -(M + 1);
          break;
        }

        M += 1;

        // Compute U1 = am^2 + bm + c and loop while this
        // value is not coprime to K.

        UU2 = A * M;
        UU1 = (UU2 + B) * M + C;

        if (VERBOSE) {
          printf("loopy %s | %s %s | %s %s %s | %s (%s)\n",
                 M.ToString().c_str(),
                 UU1.ToString().c_str(),
                 UU2.ToString().c_str(),
                 A.ToString().c_str(),
                 B.ToString().c_str(),
                 C.ToString().c_str(),
                 K.ToString().c_str(),
                 BigInt::GCD(UU1, K).ToString().c_str());
        }

      } while (BigInt::GCD(UU1, K) != 1);

      // Compute 2am + b or 2cm + b as required.
      UU2 = (UU2 << 1) + B;

      if (M >= 0) {
        // Compute c.
        B = (UU1 - UU2);
        C = B + A;
        // Compute b.
        B += UU1;
        // Compute a.
        A = UU1;
      } else {
        // Compute c.
        B = UU1 + UU2;
        C += B;
        // Compute b.
        B += UU1;
        // Compute a.
        A = UU1;
      }
    }

    if (VERBOSE)
      printf("second NSD %s %s %s\n",
             A.ToString().c_str(),
             B.ToString().c_str(),
             C.ToString().c_str());

    // We will have to solve several quadratic modular
    // equations. To do this we have to factor the modulus and
    // find the solution modulo the powers of the prime factors.
    // Then we combine them by using the Chinese Remainder
    // Theorem. The different moduli are divisors of the
    // right hand side, so we only have to factor it once.

    CHECK(!hyperbolic_recursive_solution) << "Only set in SQME below.";

    for (;;) {

      SolveQuadModEquation<QMOD_CALLBACK>(
          origin_translated,
          // Never swapping x,y on this path.
          false,
          // Coefficients and modulus
          A, B, C, BigInt::Abs(K),
          // Problem state
          A, B, C, D, E,
          M, K, U, V,
          Alpha, Beta, Div, Discr);

      // Adjust counters.
      // This modifies the factors (multiplicities) in place.
      int index;
      CHECK(indexEvenMultiplicity.size() ==
            originalMultiplicities.size());
      if (VERBOSE) printf("factors: ");
      for (index = 0; index < (int)indexEvenMultiplicity.size(); index++) {
        if (VERBOSE) printf("%d ", index);
        // Loop that increments counters.
        if (!is_descending[index]) {
          // Ascending.

          const int fidx = indexEvenMultiplicity[index];
          if (counters[index] == originalMultiplicities[index]) {
            // Next time it will be descending.
            is_descending[index] = true;
            continue;
          } else {
            BigInt UU3 = factors[fidx].first * factors[fidx].first;
            factors[fidx].second -= 2;
            // Divide by square of prime.
            K /= UU3;
            // Multiply multiplier by prime.counters[index]++
            E *= factors[fidx].first;
            counters[index] += 2;
            break;
          }
        } else {
          // Descending.
          const int fidx = indexEvenMultiplicity[index];
          if (counters[index] <= 1) {
            // Next time it will be ascending.
            is_descending[index] = false;
            continue;
          } else {
            BigInt UU3 = factors[fidx].first * factors[fidx].first;
            factors[fidx].second += 2;
            // Multiply by square of prime.
            K *= UU3;
            // Divide multiplier by prime.counters[index]++
            E /= factors[fidx].first;
            counters[index] -= 2;
            break;
          }
        }
      }

      // Note: This seems to just be a performance hint; we've changed
      // the factors array and the modulus in tandem. But it's so messy
      // to try to keep those in sync. It seems like it'd only matter if
      // the caller runs several related queries, as we do not factor
      // in a loop within this code.
      //
      // Do not try to factor the number again.

      if (index == (int)indexEvenMultiplicity.size()) {
        // All factors have been found. Exit loop.
        break;
      }
    }
    if (VERBOSE) printf(".\n");

    if (VERBOSE) {
      printf("bottom %s %s / %s %s %s %s\n",
             K.ToString().c_str(),
             E.ToString().c_str(),
             Alpha.ToString().c_str(),
             Beta.ToString().c_str(),
             GcdHomog.ToString().c_str(),
             Discr.ToString().c_str());
    }

    if (hyperbolic_recursive_solution &&
        QMOD_CALLBACK == QmodCallbackType::HYPERBOLIC) {

      // Show recursive solution.
      GetRecursiveSolution(A, B, C,
                           ABack, BBack, CBack,
                           Alpha, Beta, GcdHomog, Discr);
    }
  }

  void NegativeDiscriminant(
      bool origin_translated,
      const BigInt &A, const BigInt &B, const BigInt &C,
      const BigInt &K,
      const BigInt &D,
      const BigInt &Discr,
      const BigInt &Alpha, const BigInt &Beta,
      const BigInt &Div) {

    NonSquareDiscriminant<QmodCallbackType::ELLIPTIC>(
        origin_translated,
        A, B, C, K, D, Discr, Alpha, Beta, Div);
  }

  void PositiveDiscriminant(
      bool origin_translated,
      const BigInt &A, const BigInt &B, const BigInt &C,
      const BigInt &K,
      const BigInt &D,
      const BigInt &Discr,
      const BigInt &Alpha, const BigInt &Beta,
      const BigInt &Div) {

    NonSquareDiscriminant<QmodCallbackType::HYPERBOLIC>(
        origin_translated,
        A, B, C, K, D, Discr, Alpha, Beta, Div);
  }

  void CallbackQuadModElliptic(
      const BigInt &A, const BigInt &B, const BigInt &C,
      const BigInt &E, const BigInt &M, const BigInt &K,
      const BigInt &Alpha, const BigInt &Beta, const BigInt &Div,
      const BigInt &Discr,
      const BigInt &Value) {

    auto pqro = PerformTransformation(A, B, C, K, Value);
    if (!pqro.has_value()) {
      // No solutions because gcd(P, Q, R) > 1.
      return;
    }

    const auto &[P, Q, R] = pqro.value();

    CHECK(Discr <= 0);

    std::optional<int64_t> plow_opt = P.ToInt();
    if (plow_opt.has_value() && plow_opt.value() >= 0) {
      int64_t plow = plow_opt.value();
      if (Discr < -4 && plow == 1) {
        // Discriminant is less than -4 and P equals 1.

        NonSquareDiscrSolutionOne(
            M, E, K,
            Alpha, Beta, Div,
            BigInt(1), BigInt(0),
            Value);

        return;
      }

      if (Discr == -4) {
        // Discriminant is equal to -4.
        BigInt G = Q >> 1;

        if (plow == 1) {
          NonSquareDiscrSolutionOne(
              M, E, K,
              Alpha, Beta, Div,
              BigInt(1), BigInt(0),
              Value);

          NonSquareDiscrSolutionOne(
              M, E, K,
              Alpha, Beta, Div,
              // (Q/2, -1)
              G, BigInt(-1),
              Value);

          return;
        } if (plow == 2) {

          NonSquareDiscrSolutionOne(
              M, E, K,
              Alpha, Beta, Div,
              // ((Q/2-1)/2, -1)
              (G - 1) >> 1, BigInt(-1),
              Value);

          NonSquareDiscrSolutionOne(
              M, E, K,
              Alpha, Beta, Div,
              // ((Q/2+1)/2, -1)
              (G + 1) >> 1, BigInt(-1),
              Value);

          return;
        }
      }

      if (Discr == -3) {
        // Discriminant is equal to -3.
        if (plow == 1) {

          NonSquareDiscrSolutionOne(
              M, E, K,
              Alpha, Beta, Div,
              BigInt(1), BigInt(0),
              Value);

          NonSquareDiscrSolutionOne(
              M, E, K,
              Alpha, Beta, Div,
              // ((Q-1)/2, -1)
              (Q - 1) >> 1, BigInt(-1),
              Value);

          NonSquareDiscrSolutionOne(
              M, E, K,
              Alpha, Beta, Div,
              // ((Q+1)/2, -1)
              (Q + 1) >> 1, BigInt(-1),
              Value);

          return;
        } else if (plow == 3) {

          // printf("plow3 coverage\n");

          NonSquareDiscrSolutionOne(
              M, E, K,
              Alpha, Beta, Div,
              // ((Q+3)/6, -1)
              (Q + 3) / 6, BigInt(-1),
              Value);

          NonSquareDiscrSolutionOne(
              M, E, K,
              Alpha, Beta, Div,
              // (Q/3, -2)
              Q / 3, BigInt(-2),
              Value);

          NonSquareDiscrSolutionOne(
              M, E, K,
              Alpha, Beta, Div,
              // ((Q-3)/6, -1)
              (Q - 3) / 6, BigInt(-1),
              Value);

          return;
        }
      }
    }

    const BigInt LL = (P << 2) / Discr;
    /*
    fprintf(stderr, "P = %s, Discr = %s, LL = %s\n",
            P.ToString().c_str(), Discr.ToString().c_str(),
            LL.ToString().c_str());
    */

    // Compute bound L = sqrt(|4P/(-D)|)
    // Port note: Original code flips the sign, but on the input
    // -10 -10 -10 -10 -8 -8, that results in sqrt(-1). Alpertron's
    // sqrt function ignores the sign.
    const BigInt L = BigInt::Sqrt(BigInt::Abs(LL));

    // Initial value of last convergent: 1/0.
    BigInt U1(1);
    BigInt V1(0);
    // Initial value of next to last convergent: 0/1.
    BigInt U2(0);
    BigInt V2(1);

    // Compute continued fraction expansion of U/V = -Q/2P.
    BigInt U = -Q;
    BigInt V = P << 1;

    while (V != 0) {
      std::tie(U, U1, U2, V, V1, V2) =
        GetNextConvergent(U, U1, U2,
                                V, V1, V2);

      // Check whether the denominator of convergent exceeds bound.
      BigInt BigTmp = L - V1;
      if (BigTmp < 0) {
        // Bound exceeded, so go out.
        break;
      }

      // Test whether P*U1^2 + Q*U1*V1 + R*V1^2 = 1.
      BigInt O = (P * U1 + Q * V1) * U1 + R * V1 * V1;

      if (O == 1) {

        // a*U1^2 + b*U1*V1 + c*V1^2 = 1.
        NonSquareDiscrSolutionOne(
            M, E, K,
            Alpha, Beta, Div,
            U1, V1,
            Value);

        std::optional<int64_t> dopt = Discr.ToInt();
        if (!dopt.has_value()) break;
        int64_t d = dopt.value();
        CHECK(d < 0) << "Original code seemed to assume this.";

        if (d < -4) {
          // Discriminant is less than -4, go out.
          break;
        }

        if (d == -3 || d == -4) {
          // Discriminant is equal to -3 or -4.
          std::tie(U, U1, U2, V, V1, V2) =
            GetNextConvergent(U, U1, U2,
                                    V, V1, V2);

          NonSquareDiscrSolutionOne(
              M, E, K,
              Alpha, Beta, Div,
              U1, V1,
              Value);

          if (d == -3) {
            std::tie(U, U1, U2, V, V1, V2) =
                GetNextConvergent(U, U1, U2, V, V1, V2);

            NonSquareDiscrSolutionOne(
                M, E, K,
                Alpha, Beta, Div,
                U1, V1,
                Value);
          }

          break;
        }
      }
    }
  }

  // Discr = G^2
  void PerfectSquareDiscriminant(
      const BigInt &A, const BigInt &B, const BigInt &C,
      const BigInt &G, const BigInt &K,
      const BigInt &Alpha, const BigInt &Beta, const BigInt &Div,
      const BigInt &Discr) {
    // only used on path where A != 0
    BigInt S(0xCAFE);
    BigInt R;
    if (A == 0) {
      // Let R = gcd(b, c)
      // (bX + cY) Y = k
      R = BigInt::GCD(B, C);
    } else {
      // Multiplying by 4a we get (2aX + (b+g)Y)(2aX + (b-g)Y) = 4ak
      // Let R = gcd(2a, b+g)
      BigInt A2 = A << 1;
      R = BigInt::GCD(A2, B + G);
      // Let S = gcd(2a, b-g)
      S = BigInt::GCD(A2, B - G);
      // Let L = 4ak
      // Port note: L is dead. It was only used in teach mode.
      // L = (A * K) << 2;
    }

    if (K == 0) {

      // k equals zero.
      if (A == 0) {
        // printf("kzeroazero coverage\n");
        // Coefficient a does equals zero.
        // Solve Dy - beta = 0

        {
          LinSol sol = LinearEq(BigInt(0), Discr, -Beta);
          // Result box:
          PrintLinear(sol);
        }

        {
          // Solve bDx + cDy - b*alpha - c*beta = 0
          const BigInt Aux0 = B * Discr;
          const BigInt Aux1 = C * Discr;
          const BigInt Aux2 = -(B * Alpha + C * Beta);

          LinSol sol = LinearEq(Aux0, Aux1, Aux2);
          // Result box:
          PrintLinear(sol);
        }

      } else {
        // printf("kzeroanzero coverage\n");
        // Coefficient a does not equal zero.

        const BigInt AAlpha2 = (A * Alpha) << 1;

        {
          // Solve 2aD x + (b+g)D y = 2a*alpha + (b+g)*beta
          const BigInt Aux0 = (A * Discr) << 1;
          const BigInt Aux1 = (B + G) * Discr;
          const BigInt Aux2 = -(AAlpha2 + (B + G) * Beta);

          LinSol sol = LinearEq(Aux0, Aux1, Aux2);
          // Result box:
          PrintLinear(sol);
        }

        {
          // Solve 2aD x + (b-g)D y = 2a*alpha + (b-g)*beta
          const BigInt Aux0 = (A * Discr) << 1;
          // Port note: At some point this was erroneously
          // multiplied by the value of aux1 above.
          const BigInt Aux1 = (B - G) * Discr;
          const BigInt Aux2 = -(AAlpha2 + (B - G) * Beta);

          LinSol sol = LinearEq(Aux0, Aux1, Aux2);
          /*
          if (sol.type == LinSolType::SOLUTION_FOUND) {
            printf("bminusg coverage\n");
          }
          */

          // Result box:
          PrintLinear(sol);
        }
      }

      return;
    }

    // k does not equal zero.
    BigInt U1, U3;
    if (A == 0) {
      // printf("knzaz coverage\n");
      // If R does not divide k, there is no solution.
      U3 = K;
      U1 = R;
    } else {
      // printf("knzanz coverage\n");
      // If R*S does not divide 4ak, there is no solution.
      U1 = R * S;
      U3 = (A * K) << 2;
    }

    if (!BigInt::DivisibleBy(U3, U1)) {
      return;
    }

    const BigInt Z = BigInt::DivExact(U3, U1);

    // We have to find all factors of the right hand side.

    // Compute all factors of Z = 4ak/RS

    // Factor positive number.
    std::vector<std::pair<BigInt, int>> factors =
      BigIntFactor(BigInt::Abs(Z));

    // Do not factor again same modulus.

    // x = (NI - JM) / D(IL - MH) and y = (JL - NH) / D(IL - MH)
    // The denominator cannot be zero here.
    // H = 2a/R, I = (b+g)/R, J = F + H * alpha + I * beta
    // L = 2a/S, M = (b-g)/S, N = Z/F + L * alpha + M * beta
    // F is any factor of Z (positive or negative).
    const int nbrFactors = factors.size();

    BigInt H, I, L, M;
    if (A == 0) {
      H = BigInt(0);
      I = BigInt(1);
      // L <- b/R
      L = B / R;
      // M <- c/R
      M = C / R;
    } else {
      // 2a
      BigInt UU3 = A << 1;
      // H <- 2a/R
      H = UU3 / R;
      // L <- 2a/S
      L = UU3 / S;
      // I <- (b+g)/R
      I = (B + G) / R;
      // M <- (b-g)/S
      M = (B - G) / S;
    }


    // Compute denominator: D(IL - MH)
    const BigInt Den = Discr * (I * L - M * H);
    // O <- L * alpha + M * beta
    const BigInt O = L * Alpha + M * Beta;
    // K <- H * alpha + I * beta
    const BigInt NewK = H * Alpha + I * Beta;

    // Loop that finds all factors of Z.
    // Use Gray code to use only one big number.
    // Gray code: 0->000, 1->001, 2->011, 3->010, 4->110, 5->111, 6->101, 7->100.
    // Change from zero to one means multiply, otherwise divide.
    std::vector<int> counters(400, 0);
    std::vector<bool> is_descending(400, false);

    BigInt CurrentFactor(1);
    for (;;) {
      // Process positive divisor.
      CheckSolutionSquareDiscr(CurrentFactor,
                                     H, I, L, M, Z,
                                     Alpha, Beta, Div);
      // Process negative divisor.
      CheckSolutionSquareDiscr(-CurrentFactor,
                                     H, I, L, M, Z,
                                     Alpha, Beta, Div);

      int fidx = 0;
      int index;
      for (index = 0; index < nbrFactors; index++) {
        // Loop that increments counters.
        if (!is_descending[index]) {
          // Ascending.
          if (counters[index] == factors[fidx].second) {
            // Next time it will be descending.
            is_descending[index] = true;
            fidx++;
            continue;
          }

          const BigInt &p = factors[fidx].first;
          CurrentFactor *= p;
          counters[index]++;
          break;
        }

        if (counters[index] == 0) {
          // Descending.
          // Next time it will be ascending.
          is_descending[index] = false;
          fidx++;
          // pstFactor++;
          continue;
        }

        // XXX same
        const BigInt &p = factors[fidx].first;
        CurrentFactor /= p;

        counters[index]--;
        break;
      }

      if (index == nbrFactors) {
        // All factors have been found. Exit loop.
        break;
      }
    }
  }

  // Used for hyperbolic curve.
  //  PQa algorithm for (P+G)/Q where G = sqrt(discriminant):
  //  Set U1 to 1 and U2 to 0.
  //  Set V1 to 0 and V2 to 1.
  //  Perform loop:
  //  Compute a as floor((U + G)/V)
  //  Set U3 to U2, U2 to U1 and U1 to a*U2 + U3
  //  Set V3 to V2, V2 to V1 and V1 <- a*V2 + V3
  //  Set U to a*V - U
  //  Set V to (D - U^2)/V
  //  Inside period when: 0 <= G - U < V
  //
  // Returns true if a solution was found.
  bool ContFrac(
      bool origin_translated,
      const BigInt &Value, enum SolutionNumber solutionNbr,
      const BigInt &A, const BigInt &B, const BigInt &E,
      const BigInt &K, const BigInt &L, const BigInt &M,
      BigInt U, BigInt V, BigInt G,
      const BigInt &Alpha, const BigInt &Beta,
      const BigInt &Div, const BigInt &Discr,
      std::optional<std::pair<BigInt, BigInt>> *sol_plus,
      std::optional<std::pair<BigInt, BigInt>> *sol_minus) {

    const bool isBeven = B.IsEven();
    // If (D-U^2) is not multiple of V, exit routine.
    if (!BigInt::DivisibleBy(L - U * U, V)) {
      return false;
    }

    int periods_to_compute = origin_translated ? 2 : 1;

    BigInt U1(1);
    BigInt U2(0);
    BigInt V1(0);
    BigInt V2(1);

    // Initialize variables.
    // Less than zero means outside period.
    // Port note: Original code left startperiodv uninitialized, though
    // it was probably only accessed when startperiodu is non-negative.
    BigInt StartPeriodU(-1);
    BigInt StartPeriodV(-1);
    int index = 0;

    if (solutionNbr == SolutionNumber::SECOND) {
      index++;
    }

    bool isIntegerPart = true;

    const bool k_neg = K < 0;
    const bool a_neg = A < 0;

    int periodIndex = 0;

    bool sol_found = false;

    for (;;) {

      const bool v_neg = V < 0;

      if (V == (isBeven ? 1 : 2) &&
          ((index & 1) == (k_neg == v_neg ? 0 : 1))) {
        // Two solutions

        // Found solution.
        if (BigInt::Abs(Discr) == 5 && (a_neg != k_neg) &&
            (solutionNbr == SolutionNumber::FIRST)) {
          // Discriminant is 5 and aK < 0.
          // Use exceptional solution (U1-U2)/(V1-V2).

          // printf("aaaaaaa coverage\n");
          sol_found =
            NonSquareDiscrSolutionTwo(
                M, E, -BigInt::Abs(K),
                Alpha, Beta, Div,
                V1 - V2, U1 - U2, Value,
                sol_plus, sol_minus);

        } else {
          // Discriminant is not 5 or aK > 0.
          // Use convergent U1/V1 as solution.

          sol_found =
            NonSquareDiscrSolutionTwo(
                M, E, -BigInt::Abs(K),
                Alpha, Beta, Div,
                V1, U1, Value,
                sol_plus, sol_minus);

        }

        if (sol_found) break;
      }

      if (StartPeriodU >= 0) {
        // Already inside period.
        periodIndex++;
        if (U == StartPeriodU &&
            V == StartPeriodV &&
            // New period started.
            (periodIndex & 1) == 0) {
          // Two periods if period length is odd, one period if even.
          periods_to_compute--;
          if (periods_to_compute == 0) {
            // Go out in this case.
            break;
          }
        }

      } else if (!isIntegerPart) {
        // Check if periodic part of continued fraction has started.

        if (CheckStartOfContinuedFractionPeriod(U, V, G)) {
          StartPeriodU = U;
          StartPeriodV = V;
        }
      }

      // Get continued fraction coefficient.
      BigInt BigTmp = U + G;
      if (V < 0) {
        // If denominator is negative, round square root upwards.
        BigTmp += 1;
      }

      // Tmp1 = Term of continued fraction.
      BigInt Tmp1 = BigInt::DivFloor(BigTmp, V);
      // Update convergents.
      // U3 <- U2, U2 <- U1, U1 <- a*U2 + U3
      BigInt U3 = U2;
      U2 = U1;
      U1 = Tmp1 * U2 + U3;

      // V3 <- V2, V2 <- V1, V1 <- a*V2 + V3
      BigInt V3 = V2;
      V2 = V1;
      V1 = Tmp1 * V2 + V3;

      // Update numerator and denominator.
      U = Tmp1 * V - U;
      V = (L - U * U) / V;

      index++;
      isIntegerPart = false;
    }

    return sol_found;
  }

  void CallbackQuadModHyperbolic(bool origin_translated,
                                 const BigInt &A,
                                 const BigInt &B,
                                 const BigInt &C,
                                 const BigInt &K,

                                 const BigInt &E,
                                 const BigInt &M,
                                 const BigInt &Alpha,
                                 const BigInt &Beta,
                                 const BigInt &Div,

                                 const BigInt &Discr,
                                 const BigInt &Value) {

    auto pqro = PerformTransformation(A, B, C, K, Value);
    if (!pqro.has_value()) {
      // No solutions because gcd(P, Q, R) > 1.
      return;
    }

    // P and Q are always overwritten below.
    // R is just dead.
    // const auto &[P_, Q_, R_] = pqro.value();

    // Expected to agree because PerformTransformation doesn't modify B?
    const bool isBeven = B.IsEven();

    // Compute P as floor((2*a*theta + b)/2)
    BigInt P = (((A << 1) * Value) + B);

    if (P.IsOdd()) P -= 1;
    P >>= 1;

    // Compute Q = a*abs(K)
    BigInt Q = BigInt::Abs(K) * A;

    // Find U, V, L so we can compute the continued fraction
    // expansion of (U+sqrt(L))/V.
    BigInt L = Discr;

    BigInt U, V;
    if (isBeven) {
      U = P;
      // Argument of square root is discriminant/4.
      L >>= 2;
      V = Q;
    } else {
      // U <- 2P+1
      U = (P << 1) + 1;
      // V <- 2Q
      V = (Q << 1);
    }

    U = -U;

    // If L-U^2 is not multiple of V, there is no solution, so go out.
    if (!BigInt::DivisibleBy(L - U * U, V)) {
      // No solutions using continued fraction.
      return;
    }

    // Set G to floor(sqrt(L))
    const BigInt G = BigInt::Sqrt(L);

    std::optional<std::pair<BigInt, BigInt>> sol_plus, sol_minus;

    // Continued fraction of (U+G)/V
    if (ContFrac(origin_translated,
                 Value, SolutionNumber::FIRST,
                 A, B, E,
                 K, L, M,
                 U, V, G,
                 Alpha, Beta, Div, Discr,
                 &sol_plus, &sol_minus))
      hyperbolic_recursive_solution = true;

    // Continued fraction of (-U+G)/(-V)
    if (ContFrac(origin_translated,
                 Value, SolutionNumber::SECOND,
                 A, B, E,
                 K, L, M,
                 -U, -V, G,
                 Alpha, Beta, Div, Discr,
                 &sol_plus, &sol_minus))
      hyperbolic_recursive_solution = true;

    if (sol_plus.has_value()) {
      const auto &[X, Y] = sol_plus.value();
      // Result box:
      ShowXYOne(X, Y);
    }

    if (sol_minus.has_value()) {
      const auto &[X, Y] = sol_minus.value();
      // Result box:
      ShowXYOne(X, Y);
    }
  }

  // Copy intentional; we modify them in place (factor out gcd).
  // PS: This is where to understand the meaning of Alpha, Beta, K, Div.
  void SolveQuadEquation(BigInt A, BigInt B, BigInt C,
                         BigInt D, BigInt E, BigInt F) {
    BigInt gcd = BigInt::GCD(BigInt::GCD(A, B),
                             BigInt::GCD(BigInt::GCD(C, D),
                                         E));

    if (gcd != 0 && !BigInt::DivisibleBy(F, gcd)) {
      // F is not multiple of GCD(A, B, C, D, E) so there are no solutions.
      return;
    }

    // Divide all coefficients by GCD(A, B, C, D, E).
    // By definition, they are known to be divisible.
    if (gcd != 0) {
      A = BigInt::DivExact(A, gcd);
      B = BigInt::DivExact(B, gcd);
      C = BigInt::DivExact(C, gcd);
      D = BigInt::DivExact(D, gcd);
      E = BigInt::DivExact(E, gcd);
      F = BigInt::DivExact(F, gcd);
    }

    if (VERBOSE)
      printf("After dividing: %s %s %s %s %s %s\n",
             A.ToString().c_str(),
             B.ToString().c_str(),
             C.ToString().c_str(),
             D.ToString().c_str(),
             E.ToString().c_str(),
             F.ToString().c_str());

    // Test whether the equation is linear. A = B = C = 0.
    if (A == 0 && B == 0 && C == 0) {
      LinSol sol = LinearEq(D, E, F);
      // Result box:
      PrintLinear(sol);
      return;
    }

    // Compute discriminant: b^2 - 4ac.
    const BigInt Discr = B * B - ((A * C) << 2);


    if (Discr == 0) {
      // Port note: This code depended on uninitialized values M, K,
      // I, Alpha, Beta, Div, but I think they end up unused inside
      // SolveQuadModEquation (or callbacks) when the discriminant is
      // zero.
      // Discriminant is zero.
      DiscriminantIsZero(A, B, C, D, E, F);
      return;
    }

    // Compute gcd(a,b,c).

    BigInt UU1 = BigInt::GCD(BigInt::GCD(A, B), C);
    BigInt Div, K, Alpha, Beta;

    bool origin_translated = false;
    // Discriminant is not zero.
    if (D == 0 && E == 0) {
      // Do not translate origin.
      Div = BigInt(1);
      K = -F;
      Alpha = BigInt(0);
      Beta = BigInt(0);
    } else {
      origin_translated = true;

      Div = Discr;
      // Translate the origin (x, y) by (alpha, beta).
      // Compute alpha = 2cd - be
      Alpha = ((C * D) << 1) - (B * E);

      // Compute beta = 2ae - bd
      Beta = ((A * E) << 1) - (B * D);

      // We get the equation ax^2 + bxy + cy^2 = k
      // where k = -D (ae^2 - bed + cd^2 + fD)

      K = (-Discr) * ((A * E * E) - (B * E * D) + (C * D * D) + (F * Discr));
    }

    // If k is not multiple of gcd(A, B, C), there are no solutions.
    if (!BigInt::DivisibleBy(K, UU1)) {
      // There are no solutions.
      return;
    }

    if (K < 0) {
      // The algorithm requires the constant coefficient
      // to be positive, so we multiply both RHS and LHS by -1.
      A = -A;
      B = -B;
      C = -C;
      K = -K;
    }

    if (Discr < 0) {
      NegativeDiscriminant(origin_translated,
                           A, B, C, K, D,
                           Discr, Alpha, Beta, Div);
      return;
    }

    const BigInt G = BigInt::Sqrt(Discr);

    if (G * G == Discr) {
      // Discriminant is a perfect square.

      PerfectSquareDiscriminant(
          A, B, C, G, K,
          Alpha, Beta, Div, Discr);

      return;
    } else {
      PositiveDiscriminant(origin_translated,
                           A, B, C, K, D,
                           Discr, Alpha, Beta, Div);
    }
  }

  void QuadBigInt(const BigInt &A, const BigInt &B, const BigInt &C,
                  const BigInt &D, const BigInt &E, const BigInt &F) {
    if (output != nullptr) {
      ShowEq(A, B, C, D, E, F, "x", "y");
      ShowText(" = 0\n");
    }

    // size_t preamble_size = (output == nullptr) ? 0 : output->size();

    SolveQuadEquation(A, B, C, D, E, F);

    if (output != nullptr &&
        !solutions.any_integers &&
        solutions.quadratic.empty() &&
        solutions.linear.empty() &&
        solutions.points.empty() &&
        solutions.recursive.empty()) {
      ShowText("The equation does not have integer solutions.");
    }
  }

  Quad() {}
};

}  // namespace

Solutions QuadBigInt(const BigInt &a, const BigInt &b, const BigInt &c,
                     const BigInt &d, const BigInt &e, const BigInt &f,
                     std::string *output) {
  std::unique_ptr<Quad> quad(new Quad);
  quad->output = output;
  quad->QuadBigInt(a, b, c, d, e, f);
  return std::move(quad->solutions);
}
