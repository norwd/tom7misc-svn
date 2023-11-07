
#include "bigconv.h"

#include <string>
#include <cstring>

#include "base/logging.h"
#include "bignum/big.h"
#include "bignbr.h"
#include "base/stringprintf.h"

static constexpr bool CHECK_INVARIANTS = true;
static constexpr bool VERBOSE = false;

using namespace std;

#ifndef BIG_USE_GMP
# error This program requires GMP mode for bignum.
#endif

bool ParseBigInteger(const char *str, BigInteger *big) {
  BigInt b(str);
  BigIntToBigInteger(b, big);
  return true;
}

void BigIntToBigInteger(const BigInt &b, BigInteger *g) {
  // XXX could check up front that output has enough space...
  size_t count = 0;
  mpz_export(g->Limbs.data(), &count,
             // words are little endian.
             -1,
             // word size
             sizeof(int),
             // native byte-order
             0,
             // 31 bits per word
             1,
             b.GetRep());
  // Already overwrote the end of the array in this case,
  // but maybe we could be more useful by aborting.
  CHECK(count < MAX_LEN) << count;
  if (count == 0) {
    // BigInteger wants at least one limb always.
    g->Limbs[0].x = 0;
    g->nbrLimbs = 1;
  } else {
    g->nbrLimbs = count;
  }

  if (mpz_sgn(b.GetRep()) < 0) {
    g->sign = SIGN_NEGATIVE;
  } else {
    // 0 should also be "positive."
    g->sign = SIGN_POSITIVE;
  }
  CHECK(g->nbrLimbs != 0);
}

BigInt BigIntegerToBigInt(const BigInteger *g) {
  CHECK(g->nbrLimbs > 0) << g->nbrLimbs << " limbs [" <<
    StringPrintf("%04x", g->Limbs[0].x) << "]";
  BigInt out;
  mpz_import(out.GetRep(), g->nbrLimbs,
             // words are little-endian
             -1,
             // word size
             sizeof (int),
             // native byte-order
             0,
             // "nails": high bits to skip in each word
             1,
             g->Limbs.data());
  if (g->sign == SIGN_NEGATIVE) {
    mpz_neg(out.GetRep(), out.GetRep());
  }
  return out;
}

BigInt LimbsToBigInt(const limb *limbs, int num_limbs) {
  BigInt out;
  mpz_import(out.GetRep(), num_limbs,
             // words are little-endian
             -1,
             // word size
             sizeof (int),
             // native byte-order
             0,
             // "nails": high bits to skip in each word
             1,
             limbs);
  return out;
}

int BigIntNumLimbs(const BigInt &b) {
  static constexpr int bits_per_limb = 31;

  // Number of bits in b.
  const size_t num_bits = mpz_sizeinbase(b.GetRep(), 2);

  // Round up if needed.
  return (num_bits + bits_per_limb - 1) / bits_per_limb;
}

int BigIntToLimbs(const BigInt &b, limb *limbs) {
  size_t count = 0;
  mpz_export(limbs, &count,
             // words are little endian.
             -1,
             // word size
             sizeof(int),
             // native byte-order
             0,
             // 31 bits per word
             1,
             b.GetRep());
  if (count == 0) {
    // BigInteger wants one limb always.
    limbs[0].x = 0;
    return 1;
  }
  return count;
}

void BigIntToFixedLimbs(const BigInt &b, size_t num_limbs, limb *limbs) {
  size_t count = 0;
  mpz_export(limbs, &count,
             // words are little endian.
             -1,
             // word size
             sizeof(int),
             // native byte-order
             0,
             // 31 bits per word
             1,
             b.GetRep());

  // This also handles the case where there were no limbs.
  while (count < num_limbs) {
    limbs[count].x = 0;
    count++;
  }
}

std::string BigIntegerToString(const BigInteger *g) {
  return BigIntegerToBigInt(g).ToString();
}
