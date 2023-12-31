
#include "bignbr.h"
#include "bigconv.h"

#include <cstdint>
#include <initializer_list>
#include <numeric>

#include "bignum/big.h"
#include "bignum/big-overloads.h"
#include "base/logging.h"
#include "base/stringprintf.h"
#include "base/int128.h"
#include "ansi.h"

#include "modmult.h"
#include "bigconv.h"
#include "quadmodll.h"

static void TestDivFloor() {
  for (int n = -15; n < 16; n++) {
    for (int d = -15; d < 16; d++) {
      if (d != 0) {
        int64_t q = DivFloor64(n, d);
        BigInt Q = BigInt::DivFloor(BigInt(n), BigInt(d));
        CHECK(Q == q) << n << "/" << d << " = " << q
                      << "\nWant: " << Q.ToString();
      }
    }
  }
  printf("DivFloor " AGREEN("OK") "\n");
}

static void TestJacobi() {
  for (int64_t a : std::initializer_list<int64_t>{
      -65537, -190187234, -88, -16, -15, -9, -8, -7, -6, -5, -4, -3, -2, -1,
      0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 15, 16, 17, 31337, 123874198273}) {
    for (int64_t b : std::initializer_list<int64_t>{
        1, 3, 5, 7, 9, 119, 121, 65537, 187238710001}) {

      int bigj = BigInt::Jacobi(BigInt(a), BigInt(b));
      int j = Jacobi64(a, b);
      CHECK(bigj == j) << a << "," << b <<
        "\nGot:  " << j <<
        "\nWant: " << bigj;
    }
  }
  printf("Jacobi " AGREEN("OK") "\n");
}

static void TestGCD() {
  for (int64_t a : std::initializer_list<int64_t>{
      -65537, -190187234, -88, -2, -1,
      0, 1, 2, 3, 16, 31337, 123874198273}) {
    for (int64_t b : std::initializer_list<int64_t>{
        -23897417233, -222222, -32767, -31337, -3, -1,
        0, 1, 2, 5, 6, 120, 65536, 18723871000}) {
      const auto &[gcd, x, y] = ExtendedGCD64(a, b);
      const BigInt gcd2 = BigInt::GCD(BigInt(a), BigInt(b));
      CHECK(gcd == gcd2) << a << "," << b << ": gcd="
                         << gcd << " but BigInt::GCD=" << gcd2.ToString();
      const int64_t gcd3 = std::gcd(a, b);
      CHECK(gcd == gcd3) << a << "," << b << ": gcd="
                         << gcd << " but std::gcd=" << gcd3;

      const auto &[big_gcd, big_x, big_y] =
        BigInt::ExtendedGCD(BigInt(a), BigInt(b));

      /*
        // Note: We don't necessarily return the same solution as
        // BigInt (it seems to return the *minimal* pair). I don't
        // think we need to compute the minimal.
      CHECK(big_gcd == gcd && big_x == x && big_y == y)
        << a << "," << b <<
        "\nEGCD64:  " << gcd << " " << x << " " << y <<
        "\nBigEGCD: " << big_gcd.ToString() << " " << big_x.ToString()
        << " " << big_y.ToString();
      */

      // Even though the gcd must be 64-bit, the intermediate products
      // can overflow.
      CHECK(int128_t(a) * int128_t(x) +
            int128_t(b) * int128_t(y) == int128_t(gcd))
        << "For " << a << "," << b <<
        "\nWe have " << a << " * " << x << " + " << b << " * " << y
        << " = " << (int128_t(a) * int128_t(x) +
                     int128_t(b) * int128_t(y)) <<
        "\nBut expected the gcd: " << gcd;

      if (gcd == 1 && b != 0) {
        int64_t ainv = ModularInverse64(a, b);
        std::optional<BigInt> oainv2 = BigInt::ModInverse(BigInt(a), BigInt(b));
        CHECK(oainv2.has_value()) << a << "," << b;
        if (abs(b) != 1) {
          // We don't return the same value as BigInt for a modulus of 1,
          // but we don't care (anything is an inverse as this ring is
          // degenerate). This case is not useful for alpertron, since
          // we're taking inverses mod a prime power.
          CHECK(ainv == oainv2.value())
            << a << "," << b
            << "\nGot inv: " << ainv
            << "\nBut BigInt::ModInverse: " << oainv2.value().ToString();
        }

        // Following GMP, we use |b|
        int128_t r = (int128_t(a) * int128_t(ainv)) % int128_t(abs(b));
        if (r < 0) r += int128_t(abs(b));
        // mod b again since abs(b) could be 1, and 0 is correct in this
        // case.
        int128_t onemodb = int128_t(1) % int128_t(abs(b));
        CHECK(r == onemodb)
              << "For " << a << "," << b
              << "\nhave (" << a
              << " * " << ainv << ") % " << abs(b) << " = "
              << r
              << "\nbut want " << onemodb;
      }
    }
  }
  printf("GCD " AGREEN("OK") "\n");
}

static void TestSqrtModP() {
  {
    //    int64_t s = SqrtModP(13, 17, 118587876493, 9);
  }
}

int main(int argc, char **argv) {
  ANSI::Init();
  TestDivFloor();
  TestGCD();
  TestJacobi();

  TestSqrtModP();

  printf("All explicit tests " AGREEN("OK") "\n");
  return 0;
}
