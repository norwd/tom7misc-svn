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

// Output:
// nullopt: There are no solutions because gcd(P, Q, R) > 1
// some(P, Q, R) with gcd(P, Q, R) = 1.
std::optional<std::tuple<BigInt, BigInt, BigInt>>
PerformTransformation(const BigInt &K, const BigInt &Value) {
  // These equations become simpler because of the known
  // values of A,B,C = 1,0,1.

  // Compute P as (at^2+bt+c)/K
  const BigInt P = (Value * Value + 1) / K;

  // Compute Q <- -(2at + b).
  const BigInt Q = -(Value << 1);

  // Compute R <- aK
  const BigInt &R = K;

  const BigInt I = BigInt::GCD(BigInt::GCD(P, Q), R);
  if (I == 1) {
    return {std::make_tuple(P, Q, R)};
  }

  // No solutions because gcd(P, Q, R) > 1.
  return std::nullopt;
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

  void RecordSolutionXY(const BigInt &x, const BigInt &y) {
    // Negative values are obvious, since x and y appear only under
    // squares. x and y are also interchangeable.
    if (x >= 0 && y >= 0 && x <= y) {
      solutions.points.emplace_back(PointSolution{.X = x, .Y = y});
    }
  }

  // Obtain next convergent of continued fraction of U/V
  // Previous convergents U1/V1, U2/V2, U3/V3.
  static
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

    // Undo unimodular substitution
    {
      const auto &[Temp0, Temp1] =
        UnimodularSubstitution(M, Z, O);
      RecordSolutionXY(Temp0, Temp1);
    }

    // Z: (-tu - Kv)*E
    // O: -u*E

    // Undo unimodular substitution
    {
      const auto &[Temp0, Temp1] =
        UnimodularSubstitution(M, -Z, -O);
      RecordSolutionXY(Temp0, Temp1);
    }

    return true;
  }

  void SolutionX(BigInt Value, const BigInt &Modulus,
                 const BigInt &D, const BigInt &E,
                 const BigInt &M, const BigInt &K,
                 const BigInt &U, const BigInt &V) {

    BigInt A(1);
    BigInt B(0);
    BigInt C(1);
    BigInt Discr(-4);

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
      printf("  with 1 0 1 %s %s | %s %s | %s %s\n",
             D.ToString().c_str(),
             E.ToString().c_str(),
             M.ToString().c_str(),
             K.ToString().c_str(),
             U.ToString().c_str(),
             V.ToString().c_str());
    }

    CallbackQuadModElliptic(E, M, K,
                            Value);
  }

  // Solve congruence an^2 + bn + c = 0 (mod m) where m is different from zero.
  void SolveQuadModEquation(
      BigInt Modulus,
      const BigInt &D, const BigInt &E,
      const BigInt &M, const BigInt &K, const BigInt &U, const BigInt &V,
      const BigInt &Discr) {

    const BigInt coeffQuadr(1);
    const BigInt coeffLinear(0);
    const BigInt coeffIndep(1);

    const BigInt A(1);
    const BigInt B(0);
    const BigInt C(1);

    if (Modulus == 1) {
      // Handle this case first, since various things simplify
      // below when we know the modulus is not 1.

      // Here we have 1*n^2 + 0*b + 1 = 0 mod 1.
      // i.e.           n^2 + 1 = 0 mod 1.
      // Any n satisfies this, but I guess we just want n less
      // than the modulus, which is just 0 here.

      // (In the general version of this code, this happens in the
      // ValNn == 1 case, looping up to the Gcd of 1).

      SolutionX(BigInt(0), Modulus,
                D, E,
                M, K,
                U, V);
      return;
    }

    if (VERBOSE) {
      fprintf(stderr,
              "[SQME] %s %s %s %s\n",
              coeffQuadr.ToString().c_str(),
              coeffLinear.ToString().c_str(),
              coeffIndep.ToString().c_str(),
              Modulus.ToString().c_str());
    }

    CHECK(Modulus > 1);

    // PERF: This does get called with modulus = 1.
    // We might want to shortcut that since we can skip stuff
    // like BigInt::DivisibleBy(coeff_quadr, Modulus), given
    // the known value of the coefficients.

    BigInt coeff_quadr = BigInt::CMod(coeffQuadr, Modulus);
    if (coeff_quadr < 0) coeff_quadr += Modulus;

    BigInt coeff_linear = BigInt::CMod(coeffLinear, Modulus);
    if (coeff_linear < 0) coeff_linear += Modulus;

    BigInt coeff_indep = BigInt::CMod(coeffIndep, Modulus);
    if (coeff_indep < 0) coeff_indep += Modulus;

    CHECK(coeff_quadr == 1);
    CHECK(coeff_linear == 0);
    CHECK(coeff_indep == 1);

    BigInt GcdAll = BigInt::GCD(coeff_indep,
                                BigInt::GCD(coeff_quadr, coeff_linear));

    CHECK(GcdAll == 1);

    // For a GCD of zero here, original code would cause and ignore
    // a division by zero, then read 0 from the temporary.

    // coeff_indep must be divisible by the gcd, but this is always
    // the case because the gcd is 1.

    // Divide coefficients by gcd, but this does nothing with gcd=1.

    const BigInt &ValNn = Modulus;

    // coeff_quadr (1) can't be divisible by the modulus, since the
    // modulus is greater than 1.

    if (VERBOSE) {
      printf("[Call SolveEq] %s %s %s %s %s %s\n",
             coeff_quadr.ToString().c_str(),
             coeff_linear.ToString().c_str(),
             coeff_indep.ToString().c_str(),
             Modulus.ToString().c_str(),
             GcdAll.ToString().c_str(),
             ValNn.ToString().c_str());
    }

    // PERF pass in from earlier
    std::vector<std::pair<BigInt, int>> factors =
      BigIntFactor(Modulus);

    bool interesting = false;
    SolveEquation(
        SolutionFn([&](const BigInt &Value) {
            this->SolutionX(
                Value,
                Modulus,
                D, E,
                M, K,
                U, V);
          }),
        Modulus, factors,
        &interesting);

    if (interesting) {
      printf("INTERESTING!\n");
      solutions.interesting_coverage = true;
    }
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

  // For this trimmed down version, we know the discriminant is -4.
  void NonSquareDiscriminant(BigInt K) {
    const BigInt A(1);
    const BigInt B(0);
    const BigInt C(1);
    const BigInt D(0);

    const BigInt Discr(-4);

    // These were actually uninitialized, and probably unused?
    const BigInt U(0);
    const BigInt V(0);

    // Find GCD(a,b,c)
    BigInt GcdHomog = BigInt::GCD(BigInt::GCD(A, B), C);

    CHECK(GcdHomog == 1);

    // No need to divide by gcd of 1.

    if (K == 0) {
      // If k=0, the only solution is (X, Y) = (0, 0)
      RecordSolutionXY(BigInt(0), BigInt(0));
      return;
    }

    if (VERBOSE) {
      printf("start NSD %s %s %s | %s %s | 0 0 1\n",
             A.ToString().c_str(), B.ToString().c_str(), C.ToString().c_str(),
             K.ToString().c_str(), Discr.ToString().c_str());
    }

    // Factor independent term.

    // Note that we modify the factors (multiplicities) in place below.
    std::vector<std::pair<BigInt, int>> factors =
      BigIntFactor(BigInt::Abs(K));

    if (VERBOSE) {
      for (const auto &[f, m] : factors) {
        printf("%s^%d * ", f.ToString().c_str(), m);
      }
      printf("\n");
    }

    // Find all indices of prime factors with even multiplicity.
    // (XXX parallel. could be pair)
    // Index of prime factors with even multiplicity
    // PORT NOTE: was 1-based in original code; now 0-based
    std::vector<int> indexEvenMultiplicity, originalMultiplicities;
    indexEvenMultiplicity.reserve(factors.size());
    originalMultiplicities.reserve(factors.size());
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

    std::vector<int> counters(numFactors, 0);
    std::vector<bool> is_descending(numFactors, false);

    BigInt E = BigInt(1);
    // Loop that cycles through all square divisors of the independent term.
    BigInt M(0);

    // Skip GCD(A, K) != 1 case; A is 1 so the GCD is always 1.

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

    for (;;) {

      CHECK(A == 1);
      CHECK(B == 0);
      CHECK(C == 1);

      SolveQuadModEquation(
          // Coefficients and modulus
          BigInt::Abs(K),
          // Problem state
          D, E,
          M, K, U, V,
          Discr);

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
      printf("bottom %s %s / 0 0 1 %s\n",
             K.ToString().c_str(),
             E.ToString().c_str(),
             // (alpha, beta, gcdhomog)
             Discr.ToString().c_str());
    }
  }

  void CallbackQuadModElliptic(
      const BigInt &E, const BigInt &M, const BigInt &K,
      const BigInt &Value) {

    BigInt A(1);
    BigInt B(0);
    BigInt C(1);
    const BigInt Discr(-4);

    auto pqro = PerformTransformation(K, Value);
    if (!pqro.has_value()) {
      // No solutions because gcd(P, Q, R) > 1.
      return;
    }

    const auto &[P, Q, R] = pqro.value();

    std::optional<int64_t> plow_opt = P.ToInt();
    if (plow_opt.has_value() && plow_opt.value() >= 0) {
      const int64_t plow = plow_opt.value();

      // Discriminant is equal to -4.
      BigInt G = Q >> 1;

      if (plow == 1) {

        NonSquareDiscrSolutionOne(
            M, E, K,
            BigInt(1), BigInt(0),
            Value);

        NonSquareDiscrSolutionOne(
            M, E, K,
            // (Q/2, -1)
            G, BigInt(-1),
            Value);

        return;
      } if (plow == 2) {

        NonSquareDiscrSolutionOne(
            M, E, K,
            // ((Q/2-1)/2, -1)
            (G - 1) >> 1, BigInt(-1),
            Value);

        NonSquareDiscrSolutionOne(
            M, E, K,
            // ((Q/2+1)/2, -1)
            (G + 1) >> 1, BigInt(-1),
            Value);

        return;
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
            U1, V1,
            Value);

        std::optional<int64_t> dopt = Discr.ToInt();
        if (!dopt.has_value()) break;
        int64_t d = dopt.value();
        CHECK(d < 0) << "Original code seemed to assume this.";

        CHECK(d == -4);

        // Discriminant is equal to -3 or -4.
        std::tie(U, U1, U2, V, V1, V2) =
          GetNextConvergent(U, U1, U2,
                            V, V1, V2);

        NonSquareDiscrSolutionOne(
            M, E, K,
            U1, V1,
            Value);

        break;
      }
    }
  }

  void SolveQuadEquation(const BigInt &F) {
    const BigInt gcd(1);
    // CHECK(gcd == 1) << "Expecting GCD to always be 1: " << gcd.ToString();

    // F is always divisible by gcd of 1.
    // No need to reduce coefficients by GCD of 1.

    // Not linear. Linear requires A = B = C = 0.

    // Compute discriminant: b^2 - 4ac.
    // const BigInt Discr = B * B - ((A * C) << 2);
    // 0 - (1 * 4)
    const BigInt Discr(-4);

    CHECK(Discr == -4) << "Expecting discriminant of exactly -4.";

    // Compute gcd(a,b,c).

    BigInt UU1(1);
    // BigInt::GCD(BigInt::GCD(A, B), C);
    CHECK(UU1 == 1);
    BigInt K = -F;

    // Discriminant is not zero.
    // Do not translate origin.
    // K is always divisible by the gcd of A, B, C, since that's 1.

    CHECK(K >= 0);

    NonSquareDiscriminant(K);
  }

  void QuadBigInt(const BigInt &F) {
    if (F == 0) {
      // One solution: 0^2 + 0^2.
      RecordSolutionXY(BigInt(0), BigInt(0));
    } else if (F == -1) {
      // 0^2 + 1^2
      RecordSolutionXY(BigInt(0), BigInt(1));
    } else {
      // We could support this by just returning the empty solution
      // set, but we're not trying to have this code be fully general.
      CHECK(F < -1);
      SolveQuadEquation(F);
    }
  }

  Quad() {}
};

}  // namespace

Solutions QuadBigInt(const BigInt &f) {
  std::unique_ptr<Quad> quad(new Quad);
  quad->QuadBigInt(f);
  return std::move(quad->solutions);
}
