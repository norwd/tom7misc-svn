
// Careful! The builtin uint8 is a vector of 8 uints.
typedef uchar uint8_t;
typedef uint uint32_t;
typedef ulong uint64_t;
typedef long int64_t;
typedef atomic_uint atomic_uint32_t;

// XXX this was forked from factorize.cl. clean unused stuff

// These produce better code (to my eye) but aren't faster
// in benchmarks. Might be because they enable too much
// inlining or something like that.
#define PTX_SUB128 0
#define PTX_GEQ128 0

// This code is ported from cc-lib factorize.cc, which is in turn
// from gnu's factor utility.

// Subtracts 128-bit words.
// returns high, low

// The original code generates:
//  setp.lt.u64   %p84, %rd2123, %rd139;
//  selp.b64  %rd1080, -1, 0, %p84;
//  sub.s64   %rd1081, %rd2122, %rd140;
//  add.s64   %rd2122, %rd1081, %rd1080;
//  sub.s64   %rd2123, %rd2123, %rd139;
// but we can do better with inline PTX. The
// sub.cc and subc instructions subtract with
// carry, exactly for this kind of thing.

inline ulong2
Sub128(uint64_t ah, uint64_t al,
       uint64_t bh, uint64_t bl) {
#if PTX_SUB128
  uint64_t rl, rh;
  asm("/* 128-bit subtract */\n\t"
      "sub.cc.u64 %0, %2, %3;\n\t"
      "subc.u64 %1, %4, %5;"
      :
      // outputs %0, %1
      "=l"(rl), // low word
      "=l"(rh)  // high word
      :
      // inputs %2, %3
      "l"(al),
      "l"(bl),
      // inputs %4, %5
      "l"(ah),
      "l"(bh));

  ulong2 ret;
  ret.s0 = rh;
  ret.s1 = rl;

  // asm("/* End Sub128 */");
  return ret;
#else
  uint64_t carry = al < bl;
  ulong2 ret;
  ret.s0 = ah - bh - carry;
  ret.s1 = al - bl;
  return ret;
#endif
}

// Right-shifts a 128-bit quantity by count.
// returns high, low
inline ulong2
RightShift128(uint64_t ah, uint64_t al, int count) {
  asm("/* Right Shift 128 */");
  ulong2 ret;
  ret.s0 = ah >> count;
  ret.s1 = (ah << (64 - count)) | (al >> count);
  asm("/* End Right Shift 128 */");
  return ret;
}

// This is maybe slightly shorter than the C code, but the
// real benefit of doing this with subtraction is that we
// immediately do the same 128-bit subtraction in UDiv128Rem,
// so we end up reusing the calculation.
inline bool GreaterEq128(uint64_t ah, uint64_t al,
                         uint64_t bh, uint64_t bl) {
#if PTX_GEQ128
  asm("/* geq 128 */");


  uint64_t rl, rh;
  asm("/* 128-bit subtract */\n\t"
      "sub.cc.u64 %0, %2, %3;\n\t"
      "subc.u64 %1, %4, %5;"
      :
      // outputs %0, %1
      "=l"(rl), // low word
      "=l"(rh)  // high word
      :
      // inputs %2, %3
      "l"(al),
      "l"(bl),
      // inputs %4, %5
      "l"(ah),
      "l"(bh));

  // a >= b iff a - b >= 0. We only need the high word (bit!)
  // to test this.
  bool res = ((int64_t)rh) >= 0;
  asm("/* end geq 128 */");
  return res;
#else
  return ah > bh || (ah == bh && al >= bl);
#endif
}

// Divides (n1*2^64 + n0)/d, with n1 < d. Returns remainder.
inline uint64_t UDiv128Rem(uint64_t n1,
                           uint64_t n0,
                           uint64_t d) {
  // assert (n1 < d);
  uint64_t d1 = d;
  uint64_t d0 = 0;
  uint64_t r1 = n1;
  uint64_t r0 = n0;
  for (unsigned int i = 64; i > 0; i--) {
    ulong2 dd = RightShift128(d1, d0, 1);
    d1 = dd.s0;
    d0 = dd.s1;
    if (GreaterEq128(r1, r0, d1, d0)) {
      ulong2 rr = Sub128(r1, r0, d1, d0);
      r1 = rr.s0;
      r0 = rr.s1;
    }
  }
  return r0;
}


/* x B (mod n).  */
inline uint64_t Redcify(uint64_t r, uint64_t n) {
  return UDiv128Rem(r, 0, n);
}

/* Requires that a < n and b <= n */
inline uint64_t SubMod(uint64_t a, uint64_t b, uint64_t n) {
  uint64_t t = - (uint64_t) (a < b);
  return (n & t) + a - b;
}

inline uint64_t AddMod(uint64_t a, uint64_t b, uint64_t n) {
  return SubMod(a, n - b, n);
}

#define ll_B ((uint64_t) 1 << (64 / 2))
#define ll_lowpart(t)  ((uint64_t) (t) & (ll_B - 1))
#define ll_highpart(t) ((uint64_t) (t) >> (64 / 2))

// PERF: Can be done with madc, etc.
// returns (w1, w0)
inline ulong2 UMul128(uint64_t u, uint64_t v) {
  uint32_t ul = ll_lowpart(u);
  uint32_t uh = ll_highpart(u);
  uint32_t vl = ll_lowpart(v);
  uint32_t vh = ll_highpart(v);

  uint64_t x0 = (uint64_t) ul * vl;
  uint64_t x1 = (uint64_t) ul * vh;
  uint64_t x2 = (uint64_t) uh * vl;
  uint64_t x3 = (uint64_t) uh * vh;

  // This can't give carry.
  x1 += ll_highpart(x0);
  // But this indeed can.
  x1 += x2;
  // If so, add in the proper position.
  if (x1 < x2)
    x3 += ll_B;

  ulong2 ret;
  ret.s0 = x3 + ll_highpart(x1);
  ret.s1 = (x1 << 64 / 2) + ll_lowpart(x0);
  return ret;
}


#define NUM_PRIME_DELTAS 12
const uint8_t PRIME_DELTAS[NUM_PRIME_DELTAS] = {
  1,2,2,4,2,4,2,4,6,2,6,4,
};

/* Entry i contains (2i+1)^(-1) mod 2^8.  */
const unsigned char binvert_table[128] = {
  0x01, 0xAB, 0xCD, 0xB7, 0x39, 0xA3, 0xC5, 0xEF,
  0xF1, 0x1B, 0x3D, 0xA7, 0x29, 0x13, 0x35, 0xDF,
  0xE1, 0x8B, 0xAD, 0x97, 0x19, 0x83, 0xA5, 0xCF,
  0xD1, 0xFB, 0x1D, 0x87, 0x09, 0xF3, 0x15, 0xBF,
  0xC1, 0x6B, 0x8D, 0x77, 0xF9, 0x63, 0x85, 0xAF,
  0xB1, 0xDB, 0xFD, 0x67, 0xE9, 0xD3, 0xF5, 0x9F,
  0xA1, 0x4B, 0x6D, 0x57, 0xD9, 0x43, 0x65, 0x8F,
  0x91, 0xBB, 0xDD, 0x47, 0xC9, 0xB3, 0xD5, 0x7F,
  0x81, 0x2B, 0x4D, 0x37, 0xB9, 0x23, 0x45, 0x6F,
  0x71, 0x9B, 0xBD, 0x27, 0xA9, 0x93, 0xB5, 0x5F,
  0x61, 0x0B, 0x2D, 0x17, 0x99, 0x03, 0x25, 0x4F,
  0x51, 0x7B, 0x9D, 0x07, 0x89, 0x73, 0x95, 0x3F,
  0x41, 0xEB, 0x0D, 0xF7, 0x79, 0xE3, 0x05, 0x2F,
  0x31, 0x5B, 0x7D, 0xE7, 0x69, 0x53, 0x75, 0x1F,
  0x21, 0xCB, 0xED, 0xD7, 0x59, 0xC3, 0xE5, 0x0F,
  0x11, 0x3B, 0x5D, 0xC7, 0x49, 0x33, 0x55, 0xFF,
};

inline uint64_t Binv(uint64_t n) {
  uint64_t inv = binvert_table[(n / 2) & 0x7F]; /*  8 */
  inv = 2 * inv - inv * inv * n;
  inv = 2 * inv - inv * inv * n;
  inv = 2 * inv - inv * inv * n;
  return inv;
}

// PERF: This can be done with PTX that uses carries.
/* Modular two-word multiplication, r = a * b mod m, with mi = m^(-1) mod B.
   Both a and b must be in redc form, the result will be in redc form too.

   (Redc is "montgomery form". mi stands for modular inverse.
    See https://en.wikipedia.org/wiki/Montgomery_modular_multiplication )
*/
inline uint64_t
MulRedc(uint64_t a, uint64_t b, uint64_t m, uint64_t mi) {
  ulong2 r = UMul128(a, b);
  uint64_t rh = r.s0;
  uint64_t rl = r.s1;
  uint64_t q = rl * mi;
  uint64_t th = UMul128(q, m).s0;
  uint64_t xh = rh - th;
  if (rh < th)
    xh += m;

  return xh;
}

// computes b^e mod n. b in redc form. result in redc form.
uint64_t
PowM(uint64_t b, uint64_t e, uint64_t n, uint64_t ni, uint64_t one) {
  uint64_t y = one;

  if (e & 1)
    y = b;

  // PERF: for OpenCL, we might want to just run this loop a fixed
  // number (63 I think) of times?
  while (e != 0) {
    b = MulRedc(b, b, n, ni);
    e >>= 1;

    if (e & 1)
      y = MulRedc(y, b, n, ni);
  }

  return y;
}

inline uint64_t HighBitToMask(uint64_t x) {
  // This requires an arithmetic shift, which appears to be the
  // OpenCL behavior:
  // registry.khronos.org/OpenCL/specs/3.0-unified/html/OpenCL_C.html
  //    #operators-shift
  return (uint64_t)((int64_t)(x) >> (64 - 1));
}

uint64_t GCDOdd(uint64_t a, uint64_t b) {
  if ((b & 1) == 0) {
    uint64_t t = b;
    b = a;
    a = t;
  }
  if (a == 0)
    return b;

  /* Take out least significant one bit, to make room for sign */
  b >>= 1;

  for (;;) {
    while ((a & 1) == 0)
      a >>= 1;
    a >>= 1;

    uint64_t t = a - b;
    if (t == 0)
      return (a << 1) + 1;

    uint64_t bgta = HighBitToMask(t);

    /* b <-- min (a, b) */
    b += (bgta & t);

    /* a <-- |a - b| */
    a = (t ^ bgta) - bgta;
  }
}

// One Miller-Rabin test. Returns true if definitely composite;
// false if maybe prime.
// one and nm1 are just 1 and n-1 those numbers in redc form.
bool DefinitelyComposite(uint64_t n, uint64_t ni, uint64_t b, uint64_t q,
                         unsigned int k, uint64_t one, uint64_t nm1) {
  asm("/* definitelycomposite */");
  uint64_t y = PowM(b, q, n, ni, one);

  if (y == one || y == nm1)
    return false;

  // PERF idea: Could do this simultaneously for all bases?
  // No clear reason why it would be better, but it might be
  // worth testing.
  for (unsigned int i = 1; i < k; i++) {
    // y = y^2 mod n
    y = MulRedc(y, y, n, ni);

    if (y == nm1)
      return false;
    if (y == one)
      return true;
  }
  return true;
}

static const uint32_t WITNESSES[7] =
  { 2, 325, 9375, 28178, 450775, 9780504, 1795265022 };

// Are we done factoring? Yes if the final result is 1,
// or prime.
bool IsDoneOdd(uint64_t n) {
  // Since the input is odd, this is 1, 3, 5, 7. 3,5,7 are all
  // prime so we are done. 1 means we factored all the primes
  // out, so we are also done. Don't need to check 7 here for
  // correctness, but we might as well since it allows us to
  // succeed faster in some cases by just testing a different
  // constant.
  if (n <= 7) return true;

  // These are checked above.
  // if (n == 2) return true;
  // if (n == 3) return true;
  // if (n == 5) return true;

  // Need to check these for correctness.
  if (n == 13) return true;
  if (n == 19) return true;
  if (n == 73) return true;
  if (n == 193) return true;

  // PERF: Could perhaps move these into the return false cases, since
  // they should be exceptionally rare?
  if (n == 407521) return true;
  if (n == 299210837) return true;

  // Precomputation.
  uint64_t q = n - 1;
  // Count and remove trailing zeroes.
  int k = ctz(q);
  q >>= k;

  const uint64_t ni = Binv(n);                 /* ni <- 1/n mod B */
  const uint64_t one = Redcify(1, n);

  // This assumes WITNESSES[0] == 2.
  uint64_t a_prim = AddMod(one, one, n); /* i.e., redcify a = 2 */
  int a = 2;

  /* -1, but in redc representation. */
  const uint64_t nm1 = n - one;


  if (DefinitelyComposite(n, ni, a_prim, q, k, one, nm1))
    return false;

  // Skipping index 0, which we just did.
  for (int i = 1; i < 7; i++) {
    // Establish new base.
    a = WITNESSES[i];

    // compute a_prim from a
    {
      ulong2 r = UMul128(one, a);
      if (r.s0 == 0) {
        a_prim = r.s1 % n;
      } else {
        a_prim = UDiv128Rem(r.s0, r.s1, n);
      }
    }

    if (DefinitelyComposite(n, ni, a_prim, q, k, one, nm1)) {
      return false;
    }
  }

  // The test above detects all 64-bit composite numbers, so this
  // must be a prime.
  return true;
}

// PERF could do this as a 2D kernel?
// The big disadvantage is that we get much cheaper division by
// constants because they can be turned into weird multiply-shift-sub stuff.

// Idea here is to make a quick first pass without loops, to keep the
// threads as convergent as possible.
__kernel void TrialDivide(__global const uint64_t *restrict num,
                          // Just one "large factor" per num. This one
                          // may be composite if we don't succeed in
                          // completely factoring. It may be 1, in which
                          // case it should be ignored as a factor, or
                          // 0 (only for the input zero).
                          __global uint64_t *restrict large_factor,
                          // Up to MAX_FACTORS small factors.
                          // With the current implementation
                          // these could even be like uint_8, but then
                          // we'd need to copy for the second pass.
                          __global uint32_t *restrict small_factors,
                          // Number of small factors. High bit is set
                          // if we failed (and then target_out may be
                          // composite).
                          __global uint8_t *restrict num_factors) {

  const int idx = get_global_id(0);
  const uint64_t x = num[idx];
  uint32_t *small = &small_factors[idx * MAX_FACTORS];

  // Just a single prime factor.
  // We arbitrarily say 0 and 1 are "prime"; caller has to check
  // the large_factor if they want to treat these differently.
  if (x <= 3) {
    large_factor[idx] = x;
    num_factors[idx] = 0;
    return;
  }

  int nf = 0;
  uint64_t cur = x;

  // TODO PERF: If we're going to divide a lot of times (like for the
  // 3 case), it is possible to do some binary search, or at least
  // just first try dividing by 3^3 one time, or something like that.
#define TRY(p)                  \
        if (cur % p == 0) {     \
          cur /= p;             \
          small[nf++] = p;      \
        }

  // TODO: Empirically figure out how many times to run each of these.

  TRY(3); TRY(3); TRY(3); TRY(3); TRY(3); TRY(3); TRY(3); TRY(3);
  TRY(5); TRY(5); TRY(5); TRY(5); TRY(5); TRY(5); TRY(5); TRY(5);
  TRY(7); TRY(7); TRY(7); TRY(7); TRY(7); TRY(7); TRY(7);
  TRY(11); TRY(11); TRY(11); TRY(11); TRY(11); TRY(11); TRY(11);
  TRY(13); TRY(13); TRY(13); TRY(13); TRY(13); TRY(13); TRY(13);
  TRY(17); TRY(17); TRY(17); TRY(17); TRY(17); TRY(17);
  TRY(19); TRY(19); TRY(19); TRY(19); TRY(19); TRY(19);
  TRY(23); TRY(23); TRY(23); TRY(23); TRY(23);
  TRY(29); TRY(29); TRY(29); TRY(29); TRY(29);
  TRY(31); TRY(31); TRY(31); TRY(31);
  TRY(37); TRY(37); TRY(37); TRY(37);
  TRY(41); TRY(41); TRY(41);
  TRY(43); TRY(43); TRY(43);
  TRY(47); TRY(47);
  TRY(53); TRY(53);
  TRY(59);
  TRY(61);
  TRY(67);
  TRY(71);
  TRY(73);
  TRY(79);
  TRY(83);
  TRY(89);
  TRY(97);
  TRY(101);
  TRY(103);
  TRY(107);
  TRY(109);
  TRY(113);
  TRY(127);
  TRY(131);
#undef TRY

  const int twos = ctz(cur);
  if (twos) {
    for (int i = 0; i < twos; i++) {
      small[nf++] = 2;
    }
    cur >>= twos;
  }

  // Remaining number (odd); at this point can be composite, prime, or 1.
  large_factor[idx] = cur;

  if (IsDone(cur)) {
    // Success!
    num_factors[idx] = nf;
  } else {
    num_factors[idx] = 0x80 | nf;
  }
}
