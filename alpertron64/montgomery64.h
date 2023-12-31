#ifndef _MONTGOMERY64_H
#define _MONTGOMERY64_H

#include "base/int128.h"

// TODO: To cc-lib. Could also support 32 and 128 bit?

// A 64-bit number in montgomery form.
struct Montgomery64 {
  Montgomery64() = default;
  Montgomery64(const Montgomery64 &other) = default;
  Montgomery64 &operator =(const Montgomery64 &) = default;

  // As a representation invariant, this will be in [0, modulus).
  uint64_t x = 0;

private:
  // Use MontgomeryRep64::ToMontgomery to get valid Montgomery forms.
  explicit constexpr Montgomery64(uint64_t x) : x(x) {}
  friend struct MontgomeryRep64;
};

// #define FAKE_MONTGOMERY 1

#ifndef FAKE_MONTGOMERY

// Defines a montgomery form, i.e., a specific modulus.
struct MontgomeryRep64 {
  // The modulus we're working in.
  uint64_t modulus;
  uint64_t inv;
  // 2^64 mod modulus, which is the representation of 1.
  Montgomery64 r;
  // (2^64)^2 mod modulus.
  uint64_t r_squared;

  explicit constexpr MontgomeryRep64(uint64_t modulus) : modulus(modulus) {
    // PERF can be done in fewer steps with tricks
    inv = 1;
    // PERF don't need this many steps!
    for (int i = 0; i < 7; i++)
      inv *= 2 - modulus * inv;

    // Initialize r_squared = 2^128 mod modulus.
    uint64_t r2 = -modulus % modulus;
    for (int i = 0; i < 2; i++) {
      r2 <<= 1;
      if (r2 >= modulus)
        r2 -= modulus;
    }
    // r2 = r * 2^2 mod m

    for (int i = 0; i < 5; i++)
      r2 = Mult(Montgomery64(r2), Montgomery64(r2)).x;

    // r2 = r * (2^2)^(2^5) = 2^64
    r_squared = r2;

    r = Mult(Montgomery64(1), Montgomery64(r_squared));
  }

  constexpr Montgomery64 One() const { return r; }
  inline constexpr Montgomery64 ToMontgomery(uint64_t x) const {
    // PERF necessary?
    x %= modulus;
    return Mult(Montgomery64(x), Montgomery64(r_squared));
  }

  static inline constexpr bool Eq(Montgomery64 a, Montgomery64 b) {
    return a.x == b.x;
  }

  inline constexpr Montgomery64 Add(Montgomery64 a, Montgomery64 b) const {
    // PERF: Should be able to do this without int128
    uint128_t aa(a.x);
    uint128_t bb(b.x);
    uint128_t ss = aa + bb;
    if (ss >= (uint128_t)modulus) {
      ss -= (uint128_t)modulus;
    }
    return Montgomery64((uint64_t)ss);
  }

  inline constexpr Montgomery64 Sub(Montgomery64 a, Montgomery64 b) const {
    // PERF: Should be able to do this without int128
    int128_t aa(a.x);
    int128_t bb(b.x);
    int128_t ss = aa - bb;
    if (ss < 0) {
      ss += (int128_t)modulus;
    }
    return Montgomery64((uint64_t)ss);
  }

  inline constexpr Montgomery64 Mult(Montgomery64 a, Montgomery64 b) const {
    uint128_t aa(a.x);
    uint128_t bb(b.x);
    uint64_t r = Reduce(aa * bb);
    return Montgomery64(r);
  }

  // Convert from montgomery form back to integer.
  // Result will be in [0, modulus).
  inline constexpr uint64_t ToInt(Montgomery64 a) const {
    return Reduce((uint128_t)a.x);
  }

  // b^e mod m.
  // Note that the exponent e is a normal number.
  inline constexpr Montgomery64 Pow(Montgomery64 b, uint64_t e) const {
    // PowM(uint64_t b, uint64_t e, uint64_t n, uint64_t ni, uint64_t one) {
    Montgomery64 y = One();

    if (e & 1)
      y = b;

    while (e != 0) {
      // Square b
      b = Mult(b, b);
      e >>= 1;

      if (e & 1) {
        y = Mult(y, b);
      }
    }

    return y;
  }

 private:

  inline constexpr uint64_t Reduce(uint128_t x) const {
    // ((x % 2^64) * inv) % 2^64
    uint64_t q = Uint128Low64(x) * inv;
    // (x / 2^64) - (q * modulus) / 2^64
    int64_t a =
      (int64_t)(
          (int128_t)Uint128High64(x) -
          Uint128High64((int128_t)q * (int128_t)modulus));
    if (a < 0)
      a += modulus;
    return (uint64_t)a;
  }
};

#else

// Same interface as above, but representing each number as itself.
struct MontgomeryRep64 {
  // The modulus we're working in.
  uint64_t modulus = 1;
  uint64_t inv = 1;
  // 2^64 mod modulus, which is the representation of 1.
  Montgomery64 r = Montgomery64(1);
  // (2^64)^2 mod modulus.
  uint64_t r_squared = 1;

  explicit constexpr MontgomeryRep64(uint64_t modulus) : modulus(modulus) {
  }

  constexpr Montgomery64 One() const { return {.x = 1}; }
  inline constexpr Montgomery64 ToMontgomery(uint64_t x) const {
    x %= modulus;
    return Montgomery64{.x = x};
  }

  static inline constexpr bool Eq(Montgomery64 a, Montgomery64 b) {
    return a.x == b.x;
  }

  inline constexpr Montgomery64 Add(Montgomery64 a, Montgomery64 b) const {
    // PERF: Should be able to do this without int128
    uint128_t aa(a.x);
    uint128_t bb(b.x);
    uint128_t ss = aa + bb;
    if (ss >= (uint128_t)modulus) {
      ss -= (uint128_t)modulus;
    }
    return Montgomery64{.x = (uint64_t)ss};
  }

  inline constexpr Montgomery64 Sub(Montgomery64 a, Montgomery64 b) const {
    // PERF: Should be able to do this without int128
    int128_t aa(a.x);
    int128_t bb(b.x);
    int128_t ss = aa - bb;
    if (ss < 0) {
      ss += (int128_t)modulus;
    }
    return Montgomery64{.x = (uint64_t)ss};
  }

  inline Montgomery64 Mult(Montgomery64 a, Montgomery64 b) const {
    uint128_t aa(a.x);
    uint128_t bb(b.x);
    uint64_t r = (uint64_t)((aa * bb) % (uint128_t)modulus);
    return Montgomery64{.x = r};
  }

  // Convert from montgomery form back to integer.
  // Result will be in [0, modulus).
  inline constexpr uint64_t ToInt(Montgomery64 a) const {
    return a.x;
  }

  // b^e mod m.
  // Note that the exponent e is a normal number.
  inline Montgomery64 Pow(Montgomery64 b, uint64_t e) const {
    // PowM(uint64_t b, uint64_t e, uint64_t n, uint64_t ni, uint64_t one) {
    Montgomery64 y = One();

    if (e & 1)
      y = b;

    while (e != 0) {
      // Square b
      b = Mult(b, b);
      e >>= 1;

      if (e & 1) {
        y = Mult(y, b);
      }
    }

    return y;
  }

 private:
};

#endif

#endif
