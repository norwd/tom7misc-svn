
#include "bignum/big.h"

#include <array>
#include <vector>
#include <cstdint>
#include <utility>

// TODO: The "factor" program from GNU coreutils would
// be a good source for algorithmic improvements to this.

// TODO: Move to big-util or something so that we don't
// need to include it by default?

using namespace std;

static constexpr array<uint16_t, 1000> PRIMES = {
  2,3,5,7,11,13,17,19,23,29,
  31,37,41,43,47,53,59,61,67,71,
  73,79,83,89,97,101,103,107,109,113,
  127,131,137,139,149,151,157,163,167,173,
  179,181,191,193,197,199,211,223,227,229,
  233,239,241,251,257,263,269,271,277,281,
  283,293,307,311,313,317,331,337,347,349,
  353,359,367,373,379,383,389,397,401,409,
  419,421,431,433,439,443,449,457,461,463,
  467,479,487,491,499,503,509,521,523,541,
  547,557,563,569,571,577,587,593,599,601,
  607,613,617,619,631,641,643,647,653,659,
  661,673,677,683,691,701,709,719,727,733,
  739,743,751,757,761,769,773,787,797,809,
  811,821,823,827,829,839,853,857,859,863,
  877,881,883,887,907,911,919,929,937,941,
  947,953,967,971,977,983,991,997,1009,1013,
  1019,1021,1031,1033,1039,1049,1051,1061,1063,1069,
  1087,1091,1093,1097,1103,1109,1117,1123,1129,1151,
  1153,1163,1171,1181,1187,1193,1201,1213,1217,1223,
  1229,1231,1237,1249,1259,1277,1279,1283,1289,1291,
  1297,1301,1303,1307,1319,1321,1327,1361,1367,1373,
  1381,1399,1409,1423,1427,1429,1433,1439,1447,1451,
  1453,1459,1471,1481,1483,1487,1489,1493,1499,1511,
  1523,1531,1543,1549,1553,1559,1567,1571,1579,1583,
  1597,1601,1607,1609,1613,1619,1621,1627,1637,1657,
  1663,1667,1669,1693,1697,1699,1709,1721,1723,1733,
  1741,1747,1753,1759,1777,1783,1787,1789,1801,1811,
  1823,1831,1847,1861,1867,1871,1873,1877,1879,1889,
  1901,1907,1913,1931,1933,1949,1951,1973,1979,1987,
  1993,1997,1999,2003,2011,2017,2027,2029,2039,2053,
  2063,2069,2081,2083,2087,2089,2099,2111,2113,2129,
  2131,2137,2141,2143,2153,2161,2179,2203,2207,2213,
  2221,2237,2239,2243,2251,2267,2269,2273,2281,2287,
  2293,2297,2309,2311,2333,2339,2341,2347,2351,2357,
  2371,2377,2381,2383,2389,2393,2399,2411,2417,2423,
  2437,2441,2447,2459,2467,2473,2477,2503,2521,2531,
  2539,2543,2549,2551,2557,2579,2591,2593,2609,2617,
  2621,2633,2647,2657,2659,2663,2671,2677,2683,2687,
  2689,2693,2699,2707,2711,2713,2719,2729,2731,2741,
  2749,2753,2767,2777,2789,2791,2797,2801,2803,2819,
  2833,2837,2843,2851,2857,2861,2879,2887,2897,2903,
  2909,2917,2927,2939,2953,2957,2963,2969,2971,2999,
  3001,3011,3019,3023,3037,3041,3049,3061,3067,3079,
  3083,3089,3109,3119,3121,3137,3163,3167,3169,3181,
  3187,3191,3203,3209,3217,3221,3229,3251,3253,3257,
  3259,3271,3299,3301,3307,3313,3319,3323,3329,3331,
  3343,3347,3359,3361,3371,3373,3389,3391,3407,3413,
  3433,3449,3457,3461,3463,3467,3469,3491,3499,3511,
  3517,3527,3529,3533,3539,3541,3547,3557,3559,3571,
  3581,3583,3593,3607,3613,3617,3623,3631,3637,3643,
  3659,3671,3673,3677,3691,3697,3701,3709,3719,3727,
  3733,3739,3761,3767,3769,3779,3793,3797,3803,3821,
  3823,3833,3847,3851,3853,3863,3877,3881,3889,3907,
  3911,3917,3919,3923,3929,3931,3943,3947,3967,3989,
  4001,4003,4007,4013,4019,4021,4027,4049,4051,4057,
  4073,4079,4091,4093,4099,4111,4127,4129,4133,4139,
  4153,4157,4159,4177,4201,4211,4217,4219,4229,4231,
  4241,4243,4253,4259,4261,4271,4273,4283,4289,4297,
  4327,4337,4339,4349,4357,4363,4373,4391,4397,4409,
  4421,4423,4441,4447,4451,4457,4463,4481,4483,4493,
  4507,4513,4517,4519,4523,4547,4549,4561,4567,4583,
  4591,4597,4603,4621,4637,4639,4643,4649,4651,4657,
  4663,4673,4679,4691,4703,4721,4723,4729,4733,4751,
  4759,4783,4787,4789,4793,4799,4801,4813,4817,4831,
  4861,4871,4877,4889,4903,4909,4919,4931,4933,4937,
  4943,4951,4957,4967,4969,4973,4987,4993,4999,5003,
  5009,5011,5021,5023,5039,5051,5059,5077,5081,5087,
  5099,5101,5107,5113,5119,5147,5153,5167,5171,5179,
  5189,5197,5209,5227,5231,5233,5237,5261,5273,5279,
  5281,5297,5303,5309,5323,5333,5347,5351,5381,5387,
  5393,5399,5407,5413,5417,5419,5431,5437,5441,5443,
  5449,5471,5477,5479,5483,5501,5503,5507,5519,5521,
  5527,5531,5557,5563,5569,5573,5581,5591,5623,5639,
  5641,5647,5651,5653,5657,5659,5669,5683,5689,5693,
  5701,5711,5717,5737,5741,5743,5749,5779,5783,5791,
  5801,5807,5813,5821,5827,5839,5843,5849,5851,5857,
  5861,5867,5869,5879,5881,5897,5903,5923,5927,5939,
  5953,5981,5987,6007,6011,6029,6037,6043,6047,6053,
  6067,6073,6079,6089,6091,6101,6113,6121,6131,6133,
  6143,6151,6163,6173,6197,6199,6203,6211,6217,6221,
  6229,6247,6257,6263,6269,6271,6277,6287,6299,6301,
  6311,6317,6323,6329,6337,6343,6353,6359,6361,6367,
  6373,6379,6389,6397,6421,6427,6449,6451,6469,6473,
  6481,6491,6521,6529,6547,6551,6553,6563,6569,6571,
  6577,6581,6599,6607,6619,6637,6653,6659,6661,6673,
  6679,6689,6691,6701,6703,6709,6719,6733,6737,6761,
  6763,6779,6781,6791,6793,6803,6823,6827,6829,6833,
  6841,6857,6863,6869,6871,6883,6899,6907,6911,6917,
  6947,6949,6959,6961,6967,6971,6977,6983,6991,6997,
  7001,7013,7019,7027,7039,7043,7057,7069,7079,7103,
  7109,7121,7127,7129,7151,7159,7177,7187,7193,7207,
  7211,7213,7219,7229,7237,7243,7247,7253,7283,7297,
  7307,7309,7321,7331,7333,7349,7351,7369,7393,7411,
  7417,7433,7451,7457,7459,7477,7481,7487,7489,7499,
  7507,7517,7523,7529,7537,7541,7547,7549,7559,7561,
  7573,7577,7583,7589,7591,7603,7607,7621,7639,7643,
  7649,7669,7673,7681,7687,7691,7699,7703,7717,7723,
  7727,7741,7753,7757,7759,7789,7793,7817,7823,7829,
  7841,7853,7867,7873,7877,7879,7883,7901,7907,7919,
};

static constexpr int FIRST_OMITTED_PRIME = 7927;

#ifndef BIG_USE_GMP

std::vector<std::pair<BigInt, int>>
BigInt::PrimeFactorization(const BigInt &x, int mf) {
  // Simple trial division.
  // It would not be hard to incorporate Fermat's method too,
  // for cases that the number has factors close to its square
  // root too (this may be common?).

  // Factors in increasing order.
  std::vector<std::pair<BigInt, int>> factors;

  BigInt cur = x;
  BigInt zero(0);
  BigInt two(2);

  // Illegal input.
  if (!BigInt::Greater(x, zero))
    return factors;

  BigInt max_factor(mf);
  // Without a max factor, use the starting number itself.
  if (mf < 0 || BigInt::Greater(max_factor, x))
    max_factor = x;

  // Add the factor, or increment its exponent if it is the
  // one already at the end. This requires that the factors
  // are added in ascending order (which they are).
  auto PushFactor = [&factors](const BigInt &b) {
      if (!factors.empty() &&
          BigInt::Eq(factors.back().first, b)) {
        factors.back().second++;
      } else {
        factors.push_back(make_pair(b, 1));
      }
    };

  // First, using the prime list.
  for (int i = 0; i < (int)PRIMES.size(); /* in loop */) {
    BigInt prime(PRIMES[i]);
    if (BigInt::Greater(prime, max_factor))
      break;

    const auto [q, r] = BigInt::QuotRem(cur, prime);
    if (BigInt::Eq(r, zero)) {
      cur = q;
      if (BigInt::Greater(max_factor, cur))
        max_factor = cur;
      PushFactor(prime);
      // But don't increment i, as it may appear as a
      // factor many times.
    } else {
      i++;
    }
  }

  // Once we exhausted the prime list, do the same
  // but with odd numbers up to the square root.
  BigInt divisor((int64_t)PRIMES.back());
  divisor = BigInt::Plus(divisor, two);
  for (;;) {
    if (mf >= 0 && BigInt::Greater(divisor, max_factor))
      break;

    // TODO: Would be faster to compute ceil(sqrt(cur)) each
    // time we have a new cur, right? We can then just have
    // the single max_factor as well.
    BigInt sq = BigInt::Times(divisor, divisor);
    if (BigInt::Greater(sq, cur))
      break;

    // TODO: Is it faster to skip ones with small factors?
    const auto [q, r] = BigInt::QuotRem(cur, divisor);
    if (BigInt::Eq(r, zero)) {
      cur = q;
      PushFactor(divisor);
      // But don't increment i, as it may appear as a
      // factor many times.
    } else {
      // At least we skip even ones.
      divisor = BigInt::Plus(divisor, two);
    }
  }

  // And the number itself, which we now know is prime
  // (unless we reached the max factor).
  if (!BigInt::Eq(cur, BigInt(1))) {
    PushFactor(cur);
  }

  return factors;
}

#endif

// Oops, there is BzSqrt! Benchmark to compare (and make sure
// it has the same rounding behavior.) Unless this is much
// faster, we should probably just use the library routine.
BigInt SqrtInternal(const BigInt &xx) {
  BigInt zero(0);
  BigInt one(1);
  BigInt two(2);
  BigInt four(4);

  BigInt cur = xx;
  std::vector<BigInt> stack;
  while (BigInt::Greater(cur, one)) {
    stack.push_back(cur);
    // PERF: with shifts
    cur = BigInt::Div(cur, four);
  }

  /*
  for (const BigInt &s : stack) {
    printf("Big Sqrt %s\n", s.ToString().c_str());
  }
  */

  BigInt ret = cur;
  while (!stack.empty()) {
    // PERF shift
    BigInt r2 = BigInt::Times(ret, two);
    BigInt r3 = BigInt::Plus(r2, one);
    /*
    printf("r2 = %s, r3 = %s\n",
           r2.ToString().c_str(),
           r3.ToString().c_str());
    */

    BigInt top = std::move(stack.back());
    stack.pop_back();

    if (BigInt::Less(top, BigInt::Times(r3, r3)))
      ret = std::move(r2);
    else
      ret = std::move(r3);
  }

  return ret;
}

BigInt BigInt::RandTo(const std::function<uint64_t()> &r,
                      const BigInt &radix) {
  if (BigInt::LessEq(radix, BigInt{1})) return BigInt{0};
  // Generate mask of all 1-bits.
  BigInt mask{1};
  int bits = 0;
  while (BigInt::Less(mask, radix)) {
    // PERF shift in place!
    mask = BigInt::Times(mask, BigInt{2});
    bits++;
  }
  mask = BigInt::Minus(mask, BigInt{1});

  int u64s = (bits + 63) / 64;

  for (;;) {
    // Generate a random bit string of the right length.
    BigInt s{0};
    for (int i = 0; i < u64s; i++) {
      uint64_t w = r();
      s = BigInt::Plus(BigInt::LeftShift(s, 64), BigInt::FromU64(w));
    }

    s = BigInt::BitwiseAnd(s, mask);
    if (BigInt::LessEq(s, radix)) return s;
  }
}

#ifdef BIG_USE_GMP

/* Factoring with Pollard's rho method.

Copyright 1995, 1997-2003, 2005, 2009, 2012, 2015 Free Software
Foundation, Inc.

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see https://www.gnu.org/licenses/.  */

// assuming sorted factors array
void BigInt::InsertFactor(std::vector<std::pair<BigInt, int>> *factors, mpz_t prime) {
  // PERF binary search
  int i;
  for (i = factors->size() - 1; i >= 0; i--) {
    int cmp = mpz_cmp ((*factors)[i].first.rep, prime);
    if (cmp == 0) {
      // Increase exponent of existing factor.
      (*factors)[i].second++;
      return;
    }
    if (cmp < 0)
      break;
  }

  // PERF: btree or something?
  // Not found. Insert new factor.
  BigInt b;
  mpz_set(b.rep, prime);
  factors->insert(factors->begin() + i + 1,
                  std::make_pair(std::move(b), 1));
}

void BigInt::InsertFactorUI(std::vector<std::pair<BigInt, int>> *factors,
                             unsigned long prime) {
  // PERF can avoid allocation when it's already present
  mpz_t pz;
  mpz_init_set_ui(pz, prime);
  InsertFactor(factors, pz);
  mpz_clear(pz);
}

void BigInt::FactorUsingDivision(mpz_t t,
                                 std::vector<std::pair<BigInt, int>> *factors) {
  unsigned long int p = mpz_scan1(t, 0);
  mpz_fdiv_q_2exp(t, t, p);
  while (p) {
    // PERF: insert multiple
    InsertFactorUI(factors, 2);
    --p;
  }

  for (int i = 1; i <= (int)PRIMES.size(); /* in loop */) {
    unsigned long int p = PRIMES[i];
    if (!mpz_divisible_ui_p(t, p)) {
      i++;
      if (mpz_cmp_ui(t, p * p) < 0)
        break;
    } else {
      mpz_tdiv_q_ui(t, t, p);
      InsertFactorUI(factors, p);
    }
  }
}

static int mp_millerrabin(mpz_srcptr n, mpz_srcptr nm1, mpz_ptr x, mpz_ptr y,
                          mpz_srcptr q, unsigned long int k) {
  unsigned long int i;

  mpz_powm (y, x, q, n);

  if (mpz_cmp_ui (y, 1) == 0 || mpz_cmp (y, nm1) == 0)
    return 1;

  for (i = 1; i < k; i++) {
    mpz_powm_ui (y, y, 2, n);
    if (mpz_cmp (y, nm1) == 0)
      return 1;
    if (mpz_cmp_ui (y, 1) == 0)
      return 0;
  }
  return 0;
}

bool BigInt::MpzIsPrime(mpz_t n) {
  int k;
  bool is_prime;
  mpz_t q, a, nm1, tmp;

  if (mpz_cmp_ui (n, 1) <= 0)
    return 0;

  /* We have already factored out small primes. */
  if (mpz_cmp_ui (n, (long) FIRST_OMITTED_PRIME * FIRST_OMITTED_PRIME) < 0)
    return 1;

  mpz_inits (q, a, nm1, tmp, NULL);

  /* Precomputation for Miller-Rabin.  */
  mpz_sub_ui (nm1, n, 1);

  /* Find q and k, where q is odd and n = 1 + 2**k * q.  */
  k = mpz_scan1 (nm1, 0);
  mpz_tdiv_q_2exp (q, nm1, k);

  mpz_set_ui (a, 2);

  /* Perform a Miller-Rabin test, which finds most composites quickly.  */
  if (!mp_millerrabin (n, nm1, a, tmp, q, k)) {
    mpz_clears (q, a, nm1, tmp, NULL);
    return false;
  }

  /* Factor n-1 for Lucas.  */
  mpz_set (tmp, nm1);

  std::vector<std::pair<BigInt, int>> factors =
    PrimeFactorizationInternal(tmp);

  /* Loop until Lucas proves our number prime, or Miller-Rabin proves our
     number composite.  */
  for (int r = 1; r < (int)PRIMES.size(); r++) {
    is_prime = true;
    for (const auto &[factor, exponent_] : factors) {
      mpz_divexact(tmp, nm1, factor.rep);
      mpz_powm(tmp, a, tmp, n);
      is_prime = mpz_cmp_ui(tmp, 1) != 0;
      if (!is_prime) break;
    }

    if (is_prime)
      goto ret1;

    mpz_set_ui (a, PRIMES[r]);

    if (!mp_millerrabin (n, nm1, a, tmp, q, k)) {
      is_prime = false;
      goto ret1;
    }
  }

  fprintf (stderr, "Lucas prime test failure.  This should not happen\n");
  abort ();

 ret1:
  mpz_clears (q, a, nm1, tmp, NULL);

  return is_prime;
}

void BigInt::FactorUsingPollardRho(mpz_t n, unsigned long a,
                                   std::vector<std::pair<BigInt, int>> *factors) {
  mpz_t x, z, y, P;
  mpz_t t, t2;
  unsigned long long k, l, i;

  mpz_inits (t, t2, NULL);
  mpz_init_set_si (y, 2);
  mpz_init_set_si (x, 2);
  mpz_init_set_si (z, 2);
  mpz_init_set_ui (P, 1);
  k = 1;
  l = 1;

  while (mpz_cmp_ui (n, 1) != 0) {
    for (;;) {
      do {
        mpz_mul (t, x, x);
        mpz_mod (x, t, n);
        mpz_add_ui (x, x, a);

        mpz_sub (t, z, x);
        mpz_mul (t2, P, t);
        mpz_mod (P, t2, n);

        if (k % 32 == 1) {
          mpz_gcd (t, P, n);
          if (mpz_cmp_ui (t, 1) != 0)
            goto factor_found;
          mpz_set (y, x);
        }
      }
      while (--k != 0);

      mpz_set (z, x);
      k = l;
      l = 2 * l;
      for (i = 0; i < k; i++) {
        mpz_mul (t, x, x);
        mpz_mod (x, t, n);
        mpz_add_ui (x, x, a);
      }
      mpz_set (y, x);
    }

  factor_found:
    do {
      mpz_mul (t, y, y);
      mpz_mod (y, t, n);
      mpz_add_ui (y, y, a);

      mpz_sub (t, z, y);
      mpz_gcd (t, t, n);
    } while (mpz_cmp_ui (t, 1) == 0);

    mpz_divexact (n, n, t); /* divide by t, before t is overwritten */

    if (!MpzIsPrime(t)) {
      FactorUsingPollardRho (t, a + 1, factors);
    } else {
      InsertFactor(factors, t);
    }

    if (MpzIsPrime(n)) {
      InsertFactor(factors, n);
      break;
    }

    mpz_mod (x, x, n);
    mpz_mod (z, z, n);
    mpz_mod (y, y, n);
  }

  mpz_clears (P, t2, t, z, x, y, nullptr);
}

std::vector<std::pair<BigInt, int>>
BigInt::PrimeFactorizationInternal(mpz_t x) {
  // Factors in increasing order.
  std::vector<std::pair<BigInt, int>> factors;

  if (mpz_sgn(x) != 0) {
    FactorUsingDivision(x, &factors);

    if (mpz_cmp_ui(x, 1) != 0) {
      if (MpzIsPrime(x))
        InsertFactor(&factors, x);
      else
        FactorUsingPollardRho(x, 1, &factors);
    }
  }

  return factors;
}


std::vector<std::pair<BigInt, int>>
BigInt::PrimeFactorization(const BigInt &x, int64_t max_factor_ignored) {
  mpz_t tmp;
  mpz_init(tmp);
  mpz_set(tmp, x.rep);
  auto ret = PrimeFactorizationInternal(tmp);
  mpz_clear(tmp);
  return ret;
}

#endif
