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
#include "factorization.h"

using namespace std;

static constexpr bool SELF_CHECK = false;
static constexpr bool VERBOSE = false;

namespace {

inline int64_t DivFloor64(int64_t numer, int64_t denom) {
  // There's probably a version without %, but I verified
  // that gcc will do these both with one IDIV.
  int64_t q = numer / denom;
  int64_t r = numer % denom;
  if ((r > 0 && denom < 0) || (r < 0 && denom > 0)) {
    return q - 1;
  }
  return q;
}

struct Quad {
  // Solutions accumulated here.
  Solutions solutions;

  void RecordSolutionXY(const BigInt &X, const BigInt &Y) {
    auto xo = X.ToInt();
    auto yo = Y.ToInt();

    CHECK(xo.has_value() && yo.has_value()) << "These are squared, "
      "so for them to sum to a 64-bit number, they must be 32 bit "
      "(and here we should have room for 63!)";

    const int64_t x = xo.value();
    const int64_t y = yo.value();

    // Negative values are obvious, since x and y appear only under
    // squares. x and y are also interchangeable.
    if (x >= 0 && y >= 0 && x <= y) {
      solutions.points.emplace_back(PointSolution{
          .X = (uint64_t)x,
          .Y = (uint64_t)y
        });
    }
  }

  // Obtain next convergent of continued fraction of U/V
  // Previous convergents U1/V1, U2/V2, U3/V3.
  // For billions of 45-bit numbers, we have the maximum
  // values (log 2):
  // u: 43.84      v:  43.78     tmp:  21.35
  // u1: 22.48     v1: 21.35     tmp2: 43.78
  // u2: 21.92     v2: 20.71     tmp3: 22.48
  //
  // So we can use 64-bit numbers in here. The reason these
  // are small is that we start with the rational -Q/2P,
  // which is (2 * value) / 2 * ((value * value + 1) / K),
  // which is 1 / ((value + 1) / K),
  // which is K / (value + 1).
  // (We can't actually just pass k,value+1 though, since this
  // also uses the values of U and V directly.)

  // static
  std::tuple<int64_t, int64_t, int64_t,
             int64_t, int64_t, int64_t>
  GetNextConvergent(int64_t u, int64_t u1, int64_t u2,
                    int64_t v, int64_t v1, int64_t v2) {

    /*
    fprintf(stderr,
            "Next convergent: %lld %lld %lld | %lld %lld %lld\n",
            u, u1, u2,
            v, v1, v2);
    */

    int64_t tmp = DivFloor64(u, v);

    if (SELF_CHECK) {
      BigInt Tmp = BigInt::DivFloor(BigInt(u), BigInt(v));
      CHECK(Tmp == tmp);
    }

    // Compute new value of U and V.
    int64_t tmp2 = u - tmp * v;
    u = v;
    v = tmp2;

    // Compute new convergents: h_n = a_n*h_{n-1} + h_{n-2}
    // and also k_n = k_n*k_{n-1} + k_{n-2}
    int64_t tmp3 = tmp * u1 + u2;
    u2 = u1;
    u1 = tmp3;

    tmp *= v1;
    tmp += v2;

    v2 = v1;
    v1 = tmp;

    return std::make_tuple(u, u1, u2,
                           v, v1, v2);
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
      const BigInt &E, const BigInt &K,
      int64_t h, int64_t i,
      int64_t value) {

    BigInt H(h);
    BigInt I(i);

    /*
    fprintf(stderr, "NSDS: %s %s %s %s %lld\n",
            E.ToString().c_str(),
            K.ToString().c_str(),
            H.ToString().c_str(),
            I.ToString().c_str(),
            value);
    */

    // Port note: This used to modify the value of K based on the
    // callback type, but now we do that at the call site. (Also there
    // was something suspicious in here where it flipped the sign and
    // then set it negative.)

    // X = (tu - Kv)*E
    const BigInt Z = (H * value - K * I) * E;
    // Y = u*E
    const BigInt O = H * E;

    // (we get here with both values for two_solutions)

    // Undo unimodular substitution
    RecordSolutionXY(Z, O);
    // Z: (-tu - Kv)*E
    // O: -u*E
    RecordSolutionXY(-Z, -O);

    return true;
  }

  void SolutionX(int64_t value, uint64_t modulus,
                 uint64_t e) {
    if (VERBOSE) {
      fprintf(stderr, "SolutionX(%llu, %llu)\n", value, modulus);
    }

    // If 2*value is greater than modulus, subtract modulus.
    if (value > 0 &&
        (uint64_t)(value << 1) > modulus) {
      value -= modulus;
    }


    if (VERBOSE) {
      printf("  with 1 0 1 0 %llu | 0 %llu | 0 0\n",
             e,
             modulus);
    }

    CallbackQuadModElliptic(value, modulus, e);
  }

  // Solve congruence an^2 + bn + c = 0 (mod m) where m is different from zero.
  void SolveQuadModEquation(
      uint64_t modulus,
      // factorization of the modulus
      const std::vector<std::pair<uint64_t, int>> &factors,
      uint64_t e) {

    if (modulus == 1) {
      // Handle this case first, since various things simplify
      // below when we know the modulus is not 1.

      // Here we have 1*n^2 + 0*b + 1 = 0 mod 1.
      // i.e.           n^2 + 1 = 0 mod 1.
      // Any n satisfies this, but I guess we just want n less
      // than the modulus, which is just 0 here.

      // (In the general version of this code, this happens in the
      // ValNn == 1 case, looping up to the Gcd of 1).

      SolutionX(0, 1, e);
      return;
    }

    if (VERBOSE) {
      fprintf(stderr,
              "[SQME] 1 0 1 %llu\n",
              // (quad, linear, indep)
              modulus);
    }

    CHECK(modulus > 1);

    // This used to mod each coefficient by the modulus,
    // but this will not change the values 1,0,1.

    // For a GCD of zero here, original code would cause and ignore
    // a division by zero, then read 0 from the temporary.

    // coeff_indep must be divisible by the gcd, but this is always
    // the case because the gcd is 1.

    // Divide coefficients by gcd, but this does nothing with gcd=1.

    // coeff_quadr (1) can't be divisible by the modulus, since the
    // modulus is greater than 1.

    if (VERBOSE) {
      printf("[Call SolveEq] 1 0 1 %llu 1 %llu\n",
             // (quad, linear, indep)
             modulus,
             // (gcdall)
             // formerly Nn:
             modulus);
    }

    bool interesting = false;
    std::vector<uint64_t> values =
      SolveEquation(
          modulus, factors,
          &interesting);

    for (uint64_t value : values) {
      SolutionX(
          value,
          modulus,
          e);
    }

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

  //
  // CHECK(gcd == 1) << "Expecting GCD to always be 1: " << gcd.ToString();

  // F is always divisible by gcd of 1.
  // No need to reduce coefficients by GCD of 1.

  // Not linear. Linear requires A = B = C = 0.

  // Compute discriminant: b^2 - 4ac.
  // const BigInt Discr = B * B - ((A * C) << 2);
  // 0 - (1 * 4)
  // const BigInt Discr(-4);
  // CHECK(Discr == -4) << "Expecting discriminant of exactly -4.";

  // Compute gcd(a,b,c).

  // BigInt UU1(1);
  // BigInt::GCD(BigInt::GCD(A, B), C);
  // CHECK(UU1 == 1);
  // const BigInt K(f);

  // Discriminant is not zero.
  // Do not translate origin.
  // K is always divisible by the gcd of A, B, C, since that's 1.
  void SolveQuadEquation(uint64_t k,
                         // PERF avoid copying?
                         std::vector<std::pair<uint64_t, int>> factors) {
    /*
    const BigInt A(1);
    const BigInt B(0);
    const BigInt C(1);
    const BigInt D(0);

    const BigInt Discr(-4);

    // These were actually uninitialized, and probably unused?
    const BigInt U(0);
    const BigInt V(0);
    */

    // Gcd is always 1.

    // No need to divide by gcd of 1.

    if (k == 0) {
      // If k=0, the only solution is (X, Y) = (0, 0)
      RecordSolutionXY(BigInt(0), BigInt(0));
      return;
    }

    if (VERBOSE) {
      printf("start NSD 1 0 1 | %llu -4 | 0 0 1\n", k);
    }

    // Factor independent term.

    // Note that we modify the factors (multiplicities) in place below.
    // std::vector<std::pair<BigInt, int>> factors =
    // BigIntFactor(BigInt::Abs(K));

    CHECK(k > 1);
    // fprintf(stderr, "(outer) Factoring %llu\n", k);
    if (SELF_CHECK) {
      std::vector<std::pair<uint64_t, int>> ref_factors =
        Factorization::Factorize(k);
      CHECK(factors == ref_factors);
    }

    if (VERBOSE) {
      for (const auto &[f, m] : factors) {
        printf("%llu^%d * ", f, m);
      }
      printf("\n");
    }

    // Find all indices of prime factors with even multiplicity.

    // Index of prime factors with even multiplicity
    // Port note: was 1-based in original code; now 0-based
    struct EvenMultiplicity {
      // Index in the factors vectgor.
      int index;
      // The multiplicity of the factor, rounded down to the nearest
      // even number. This was called "originalMultiplicities" in
      // the original, confusingly.
      int even;

      int counter;
      bool is_descending;
    };

    std::vector<EvenMultiplicity> even_multiplicity;
    even_multiplicity.reserve(factors.size());
    const int numFactors = factors.size();
    for (int i = 0; i < numFactors; i++) {
      const auto &[fact, multiplicity] = factors[i];
      if (multiplicity > 1) {
        // At least prime is squared.
        // Port note: The original code stored factorNbr, which was 1-based
        // because of the factor header.
        // Convert multiplicity to even.
        even_multiplicity.push_back({
            .index = i,
            .even = multiplicity & ~1,
            .counter = 0,
            .is_descending = false,
          });
      }
    }

    // e is always some product of factors, so it fits in 64 bits
    // and is positive.
    uint64_t e = 1;
    // Loop that cycles through all square divisors of the independent term.

    // Skip GCD(A, K) != 1 case; A is 1 so the GCD is always 1.

    if (VERBOSE)
      printf("second NSD 1 0 1\n");

    // We will have to solve several quadratic modular
    // equations. To do this we have to factor the modulus and
    // find the solution modulo the powers of the prime factors.
    // Then we combine them by using the Chinese Remainder
    // Theorem. The different moduli are divisors of the
    // right hand side, so we only have to factor it once.

    for (;;) {

      // This code maintains the factor list as it computes different k.
      if (SELF_CHECK) {
        uint64_t product = 1;
        for (const auto &[p, e] : factors) {
          for (int i = 0; i < e; i++)
            product *= p;
        }
        CHECK(product == k);
      }

      SolveQuadModEquation(
          // Coefficients and modulus
          k,
          factors,
          // Problem state
          e);

      // Adjust counters.
      // This modifies the factors (multiplicities) in place.
      // PERF pass 'em!
      int index;
      if (VERBOSE) printf("factors: ");
      for (index = 0; index < (int)even_multiplicity.size(); index++) {
        EvenMultiplicity &even = even_multiplicity[index];

        if (VERBOSE) printf("%d ", index);
        // Loop that increments counters.

        const int factor_idx = even.index;
        if (!even.is_descending) {
          // Ascending.

          if (even.counter == even.even) {
            // Next time it will be descending.
            even.is_descending = true;
            continue;
          } else {
            uint64_t uu3 =
              factors[factor_idx].first * factors[factor_idx].first;
            factors[factor_idx].second -= 2;
            // Divide by square of prime.
            k /= uu3;
            // Multiply multiplier by prime.counters[index]++
            e *= factors[factor_idx].first;
            even.counter += 2;
            break;
          }
        } else {
          // Descending.

          if (even.counter <= 1) {
            // Next time it will be ascending.
            even.is_descending = false;
            continue;
          } else {
            uint64_t uu3 =
              factors[factor_idx].first * factors[factor_idx].first;
            factors[factor_idx].second += 2;
            // Multiply by square of prime.
            k *= uu3;
            // Divide multiplier by prime.counters[index]++
            e /= factors[factor_idx].first;
            even.counter -= 2;
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

      if (index == (int)even_multiplicity.size()) {
        // All factors have been found. Exit loop.
        break;
      }
    }
    if (VERBOSE) printf(".\n");

    if (VERBOSE) {
      printf("bottom %llu %llu / 0 0 1 -4\n", k, e);
      // (alpha, beta, gcdhomog, discr)
    }
  }

  void CallbackQuadModElliptic(
      int64_t value, uint64_t modulus, uint64_t e) {

    constexpr int64_t discr = -4;

    const BigInt K(modulus);
    const BigInt Value(value);
    const BigInt E(e);

    // PerformTransformation:
    // These equations become simpler because of the known
    // values of A,B,C = 1,0,1.

    // Compute P as (at^2+bt+c)/K
    const BigInt P = (Value * Value + 1) / K;
    const std::optional<int64_t> po = P.ToInt();
    // TODO: Some argument why this always fits in 64 bits?
    // I checked billions of samples, at least, and it was
    // never larger magnitude than the input value.
    CHECK(po.has_value());
    const int64_t p = po.value();

    // Compute Q <- -(2at + b).
    const BigInt Q = -(Value << 1);
    const int64_t q = -(value << 1);

    // Compute R <- aK
    const BigInt &R = K;

    if (BigInt::GCD(BigInt::GCD(P, Q), R) != 1) {
      // No solutions.
      return;
    }

    // Below we assert p >= 0...
    if (p >= 0) {

      // Discriminant is equal to -4.
      // BigInt G = Q >> 1;
      int64_t g = q >> 1;
      // XXX equal to -value, right?

      if (p == 1) {

        NonSquareDiscrSolutionOne(
            E, K,
            1, 0,
            value);

        NonSquareDiscrSolutionOne(
            E, K,
            // (Q/2, -1)
            g, -1,
            value);

        return;
      } if (p == 2) {

        NonSquareDiscrSolutionOne(
            E, K,
            // ((Q/2-1)/2, -1)
            (g - 1) >> 1, -1,
            value);

        NonSquareDiscrSolutionOne(
            E, K,
            // ((Q/2+1)/2, -1)
            (g + 1) >> 1, -1,
            value);

        return;
      }
    }

    CHECK(discr == -4);
    // const BigInt LL = -P;

    /*
    fprintf(stderr, "P = %s, discr = %lld, Q = %s\n",
            P.ToString().c_str(),
            discr,
            Q.ToString().c_str());
    */

    // Compute bound L = sqrt(|4P/(-D)|)
    // Port note: Original code flips the sign, but on the input
    // -10 -10 -10 -10 -8 -8, that results in sqrt(-1). Alpertron's
    // sqrt function ignores the sign.
    //
    // LL was P * 4 / discr, but this is just -P. Since we took
    // the absolute value before square root, we can just sqrt p.
    //
    // It was (Value^2 + 1) / K, which should be non-negative.
    CHECK(P >= 0);
    const BigInt L = BigInt::Sqrt(P);

    // Initial value of last convergent: 1/0.
    int64_t u1 = 1;
    int64_t v1 = 0;
    // Initial value of next to last convergent: 0/1.
    int64_t u2 = 0;
    int64_t v2 = 1;

    // Compute continued fraction expansion of U/V = -Q/2P.
    int64_t u = -q;
    int64_t v = p << 1;

    while (v != 0) {
      std::tie(u, u1, u2, v, v1, v2) =
        GetNextConvergent(u, u1, u2,
                          v, v1, v2);

      // Check whether the denominator of convergent exceeds bound.
      if (L < v1) {
        // Bound exceeded, so go out.
        break;
      }

      // Test whether P*U1^2 + Q*U1*V1 + R*V1^2 = 1.
      BigInt O = (P * u1 + Q * v1) * u1 + (R * v1) * v1;

      if (O == 1) {

        // a*U1^2 + b*U1*V1 + c*V1^2 = 1.
        NonSquareDiscrSolutionOne(
            E, K,
            u1, v1,
            value);

        CHECK(discr == -4);

        // Discriminant is equal to -3 or -4.
        std::tie(u, u1, u2, v, v1, v2) =
          GetNextConvergent(u, u1, u2,
                            v, v1, v2);

        NonSquareDiscrSolutionOne(
            E, K,
            u1, v1,
            value);

        break;
      }
    }
  }


  // F is non-negative (this is the reverse sign of the original
  // alpertron).
  void SolveQuad(uint64_t f,
                 const std::vector<std::pair<uint64_t, int>> &factors) {
    if (f == 0) {
      // One solution: 0^2 + 0^2.
      RecordSolutionXY(BigInt(0), BigInt(0));
    } else if (f == 1) {
      // 0^2 + 1^2
      RecordSolutionXY(BigInt(0), BigInt(1));
    } else {
      CHECK(f > 1);
      SolveQuadEquation(f, factors);
    }
  }

  Quad() {}
};

}  // namespace

Solutions SolveQuad(uint64_t f,
                    const std::vector<std::pair<uint64_t, int>> &factors) {
  std::unique_ptr<Quad> quad(new Quad);
  quad->SolveQuad(f, factors);
  return std::move(quad->solutions);
}