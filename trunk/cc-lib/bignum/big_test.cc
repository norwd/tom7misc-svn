
#include "bignum/big.h"
#ifndef BIG_USE_GMP
#include "bignum/bign.h"
#endif

#include <bit>
#include <cstdint>
#include <array>
#include <vector>
#include <utility>

#include "base/logging.h"
#include "base/stringprintf.h"
#include "timer.h"

using int64 = int64_t;
using namespace std;

static void TestToString() {
  for (int i = -100000; i < 100000; i++) {
    BigInt bi(i);
    string s = bi.ToString();
    CHECK(!s.empty());
    CHECK(BigInt::Eq(BigInt(s), bi));
    // We do some resizing tricks in the GMP version, which
    // can inadvertently leave nul bytes.
    for (char c : s) CHECK(c != 0) << i;
  }
}

static void CopyAndAssign() {
  BigInt a{1234};
  BigInt b{5678};
  BigInt c = a;
  a = b;
  CHECK(BigInt::Eq(a, b));
  CHECK(!BigInt::Eq(a, c));
  {
    BigInt z{8888};
    a = z;
  }
  CHECK(a.ToInt() == 8888);

  string s = "190237849028374901872390876190872349817230948719023874190827349817239048712903847190283740918273490817230948798767676767676738347482712341";
  BigInt big(s);
  CHECK(big.ToString() == s);
  big = a;
  CHECK(big.ToInt() == 8888);
}

static void LowWord() {
  BigInt a{1234};
  BigInt b{5678};
  // Technically this can be anything, but it would probably be a bug if
  // small integers don't get different values.
  CHECK(BigInt::LowWord(a) != BigInt::LowWord(b));

  string s = "190237849028374901872390876190872349817230948719023874190827349817239048712903847190283740918273490817230948798767676767676738347482712341";
  BigInt big(s);
  CHECK(BigInt::LowWord(a) != BigInt::LowWord(big));
}

static void TestMod() {
  CHECK(BigInt::Eq(BigInt::Mod(BigInt{3}, BigInt{5}), BigInt{3}));
  CHECK(BigInt::Eq(BigInt::Mod(BigInt{7}, BigInt{5}), BigInt{2}));
  CHECK(BigInt::Eq(BigInt::Mod(BigInt{10}, BigInt{5}), BigInt{0}));
  CHECK(BigInt::Eq(BigInt::Mod(BigInt{-1}, BigInt{5}), BigInt{4}));
}

static void TestEq() {
  BigInt a{1234};
  BigInt b{5678};

  CHECK(BigInt::Eq(BigInt::Times(a, b), 7006652));
}

static void TestPow() {
  BigRat q(11,15);

  BigRat qqq = BigRat::Times(q, BigRat::Times(q, q));
  BigRat qcubed = BigRat::Pow(q, 3);
  printf("%s vs %s\n", qqq.ToString().c_str(),
         qcubed.ToString().c_str());
  CHECK(BigRat::Eq(qqq, qcubed));
}

// TODO: Test/document behavior on negative inputs
static void TestQuotRem() {
  BigInt a(37);
  BigInt b(5);

  const auto [q, r] = BigInt::QuotRem(a, b);
  CHECK(BigRat::Eq(q, BigInt(7)));
  CHECK(BigRat::Eq(r, BigInt(2)));
}

static void TestPrimeFactors() {
  auto FTOS = [](const std::vector<std::pair<BigInt, int>> &fs) {
      string s;
      for (const auto &[b, i] : fs) {
        StringAppendF(&s, "%s^%d ", b.ToString().c_str(), i);
      }
      return s;
    };

  BigInt bi31337(31337);
  {
    std::vector<std::pair<BigInt, int>> factors =
      BigInt::PrimeFactorization(bi31337);

    CHECK(factors.size() == 1);
    CHECK(factors[0].second == 1);
    CHECK(BigInt::Eq(factors[0].first, bi31337));
  }

  {
    BigInt x(31337 * 71);
    std::vector<std::pair<BigInt, int>> factors =
      BigInt::PrimeFactorization(x);

    CHECK(factors.size() == 2) << FTOS(factors);
    CHECK(BigInt::Eq(factors[0].first, BigInt(71)));
    CHECK(factors[0].second == 1) << factors[0].second;
    CHECK(BigInt::Eq(factors[1].first, bi31337));
    CHECK(factors[1].second == 1) << factors[0].second;
  }

  {
    BigInt bi31337sq(31337 * 31337);
    std::vector<std::pair<BigInt, int>> factors =
      BigInt::PrimeFactorization(bi31337sq);

    CHECK(factors.size() == 1) << FTOS(factors);
    CHECK(BigInt::Eq(factors[0].first, bi31337));
    CHECK(factors[0].second == 2) << factors[0].second;
  }

  {
    BigInt x(1);
    // Must all be distinct and prime
    const array f = {
      2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 37, 41, 43, 47, 419,
      541, 547,
    };
    for (int factor : f) {
      x = BigInt::Times(x, BigInt(factor));
    }

    std::vector<std::pair<BigInt, int>> factors =
      BigInt::PrimeFactorization(x);

    CHECK(factors.size() == f.size());
    for (int i = 0; i < (int)f.size(); i++) {
      CHECK(factors[i].second == 1);
      CHECK(BigInt::Eq(factors[i].first, BigInt(f[i])));
    }
  }

  {
    // 100-digit prime; trial factoring will not succeed!
    BigInt p1("207472224677348520782169522210760858748099647"
              "472111729275299258991219668475054965831008441"
              "6732550077");
    BigInt p2("31337");

    BigInt x = BigInt::Times(p1, p2);

    // Importantly, we set a max factor (greater than p2, and
    // greater than the largest value in the primes list).
    std::vector<std::pair<BigInt, int>> factors =
      BigInt::PrimeFactorization(x, 40000);

    CHECK(factors.size() == 2);
    CHECK(factors[0].second == 1);
    CHECK(factors[1].second == 1);
    CHECK(BigInt::Eq(factors[0].first, p2));
    CHECK(BigInt::Eq(factors[1].first, p1));
  }

}

static void BenchDiv2() {
  double total_sec = 0.0;
  const int64 iters = 100;
  const BigInt zero(0);
  const BigInt two(2);
  for (int i = 0; i < iters; i++) {
    BigInt x = BigInt::Pow(two, 30000);
    Timer t;

    for (;;) {
      const auto [q, r] = BigInt::QuotRem(x, two);
      if (BigInt::Eq(r, zero)) {
        x = q;
      } else {
        break;
      }
    }

    CHECK(BigInt::Eq(x, BigInt(1)));

    total_sec += t.Seconds();
    if (i % 10 == 0) {
      printf("%d/%lld\n", i, iters);
    }
  }

  printf("%lld iters in %.5fs = %.3f/s\n",
         iters, total_sec, iters / total_sec);
}

static void TestPi() {
  printf("----\n");
  {
    BigInt i{1234567LL};
    BigInt j{33LL};
    BigInt k = BigInt::Times(i, j);
    BigInt m("102030405060708090987654321");

    printf("Integer: %s %s %s\n%s\n",
           i.ToString().c_str(),
           j.ToString().c_str(),
           k.ToString().c_str(),
           m.ToString().c_str());
    fflush(stdout);
  }

  BigRat sum;
  for (int i = 0; i < 10000; i++) {
    // + 1/1, - 1/3, + 1/5
    BigRat term{(i & 1) ? -1 : 1,
        i * 2 + 1};
    sum = BigRat::Plus(sum, term);
    if (i < 50) {
      BigRat tpi = BigRat::Times(sum, BigRat{4,1});
      printf("Approx pi: %s = %f\n",
             tpi.ToString().c_str(),
             tpi.ToDouble());
      fflush(stdout);
    } else if (i % 1000 == 0) {
      printf("%d...\n", i);
      fflush(stdout);
    }
  }

  BigRat res = BigRat::Times(sum, BigRat(4, 1));
  printf("Final approx pi: %s\n",
         res.ToString().c_str());
  fflush(stdout);


  // This sequence converges REALLY slow!
  BigRat pi_lb(314, 100);
  BigRat pi_ub(315, 100);

  CHECK(BigRat::Compare(pi_lb, pi_ub) == -1);
  CHECK(BigRat::Compare(pi_lb, res) == -1);
  CHECK(BigRat::Compare(res, pi_ub) == -1);
}

static void TestLeadingZero() {
  // TODO: Test this through BigNum interface. Ideally
  // the test should not care about the implementation.
#ifndef BIG_USE_GMP
  CHECK(BnnNumLeadingZeroBitsInDigit(1ULL) == 63);
  CHECK(BnnNumLeadingZeroBitsInDigit(0b1000ULL) == 60);
  CHECK(BnnNumLeadingZeroBitsInDigit(~0ULL) == 0);
  CHECK(BnnNumLeadingZeroBitsInDigit(0ULL) == 64);

  CHECK(std::countl_zero<BigNumDigit>(1ULL) == 63);
  CHECK(std::countl_zero<BigNumDigit>(0b1000ULL) == 60);
  CHECK(std::countl_zero<BigNumDigit>(~0ULL) == 0);
  CHECK(std::countl_zero<BigNumDigit>(0ULL) == 64);
#endif
}

static void TestToInt() {
# define ROUNDTRIP(x) do {                                              \
  BigInt bi((int64_t)(x));                                              \
  std::optional<int64_t> io = bi.ToInt();                               \
  CHECK(io.has_value()) << StringPrintf("%llx", x) << "("               \
                        << bi.ToString(16) << ")";                      \
  CHECK((x) == io.value()) << StringPrintf("%llx", x) << " vs "         \
                           << bi.ToString(16);                          \
} while (0)

  ROUNDTRIP(0);
  ROUNDTRIP(1);
  ROUNDTRIP(-1);
  ROUNDTRIP(0x7FFFFFFEL);
  ROUNDTRIP(0x7FFFFFFFL);
  ROUNDTRIP(0x80000000LL);
  ROUNDTRIP(0x80000000LL);
  ROUNDTRIP(0x7FFFFFFFFFFFFFFELL);
  ROUNDTRIP(0x7FFFFFFFFFFFFFFFLL);

# define NOROUNDTRIP(bi) do {                                     \
    std::optional<int64_t> io = (bi).ToInt();                     \
    CHECK(!io.has_value()) << #bi << " =\n" << \
      bi.ToString() << " " << io.value();     \
  } while (false)
  NOROUNDTRIP(BigInt::Plus(BigInt(0x7FFFFFFFFFFFFFFFLL), BigInt(1)));
  NOROUNDTRIP(BigInt::Minus(BigInt::Negate(BigInt(0x7FFFFFFFFFFFFFFFLL)), BigInt(1)));
  NOROUNDTRIP(BigInt::Times(BigInt(0x7FFFFFFFFFFFFFFFLL), BigInt(10000)));

# undef ROUNDTRIP
# undef NOROUNDTRIP
}

static uint64_t Sqrt64Nuprl(uint64_t xx) {
  // printf("SqrtNuprl(%llu)\n", xx);
  if (xx <= 1) return xx;
  // z = xx / 4
  uint64_t z = xx >> 2;
  uint64_t r2 = 2 * Sqrt64Nuprl(z);
  uint64_t r3 = r2 + 1;
  // printf("r2 = %llu, r3 = %llu\n", r2, r3);
  return (xx < r3 * r3) ? r2 : r3;
}

static void TestSqrt() {
  for (const uint64_t u : std::initializer_list<uint64_t>{
        1234567, 0x7FFFFFFFFFFFFFFFULL, 121,
        0, 1, 2, 3, 4, 5, 6, 9999999999999ULL, 31337 * 31337ULL}) {
    BigInt ub(u);

    // Sqare root of squares should be equal
    BigInt uu = BigInt::Times(ub, ub);
    CHECK(BigInt::Eq(BigInt::Sqrt(uu), ub)) << u;

    {
      const auto [v, vrem] = BigInt::SqrtRem(uu);
      CHECK(BigInt::Eq(v, ub));
      CHECK(BigInt::Eq(vrem, BigInt{0}));
    }

    uint64_t us64 = Sqrt64Nuprl(u);
    BigInt usb = BigInt::Sqrt(ub);
    CHECK(BigInt::Eq(usb, BigInt(us64))) << u << " " << us64;
    const auto [vsb, vrem] = BigInt::SqrtRem(ub);
    CHECK(BigInt::Eq(vsb, BigInt(us64)));
    BigInt back = BigInt::Plus(BigInt::Times(vsb, vsb), vrem);
    CHECK(BigInt::Eq(back, ub));
  }
}

static void TestGCD() {
  BigInt g = BigInt::GCD(BigInt(8), BigInt(12));
  CHECK(BigInt::Eq(g, BigInt(4)));
}

static void TestShift() {
  BigInt a{"12398471982735675717171221"};
  BigInt b = BigInt::LeftShift(a, 18);
  BigInt c = BigInt::Times(a, BigInt{262144});
  CHECK(BigInt::Eq(b, c));
}

static void TestDivExact() {
  BigInt a{"23984727341"};
  BigInt b{"12737177354116809923874293874113"};
  BigInt c = BigInt::Times(a, b);
  BigInt d = BigInt::DivExact(c, a);
  BigInt e = BigInt::DivExact(c, b);
  CHECK(BigInt::Eq(b, d));
  CHECK(BigInt::Eq(a, e));
}

int main(int argc, char **argv) {
  printf("Start.\n");
  fflush(stdout);

  CopyAndAssign();
  TestToString();

  TestEq();
  LowWord();
  TestMod();
  TestToInt();
  TestGCD();

  TestDivExact();

  TestShift();

  TestLeadingZero();

  TestPow();
  TestQuotRem();
  TestPrimeFactors();

  TestPi();

  BenchDiv2();

  TestSqrt();

  printf("OK\n");
}
