
#include "fluint16.h"
#include "fluint8.h"

#include "crypt/lfsr.h"
#include "base/logging.h"

#define CHECK_EQ16(a, b) do {                               \
    uint16_t aa = (a), bb = (b);                                \
    CHECK(aa == bb) << #a << " vs " << #b << ":\n" <<       \
      StringPrintf("%04x vs %04x (%d vs %d)",               \
                   aa, bb, aa, bb);                         \
} while (false)

static constexpr uint16_t LFSR_POLY1 = 0xBDDD;
static constexpr uint16_t LFSR_POLY2 = 0xD008;

template<uint16_t POLY>
void ManyPairs(std::function<void(uint16_t, uint16_t, Fluint16, Fluint16)> f) {
  uint16_t state1 = 0x0000, state2 = 0xCAFE;

  for (int i = 0; i < 65536; i++) {
    state1++;
    state2 = LFSRNext<uint16_t, POLY>(state2);
    f(state1, state2, Fluint16(state1), Fluint16(state2));
  }

  for (int i = 0; i < 65536; i++) {
    state1++;
    state2 = LFSRNext<uint16_t, POLY>(state2);
    f(~state2, state1, Fluint16(~state2), Fluint16(state1));
  }
};

template<uint16_t POLY>
void ManyPairs8(std::function<void(uint16_t, uint8_t, Fluint16, Fluint8)> f) {
  uint16_t state1 = 0x0000, state2 = 0xCAFE;

  // Each of the words
  for (int i = 0; i < 65536; i++) {
    state1++;
    state2 = LFSRNext<uint16_t, POLY>(state2);
    f(state1, (uint8)state2, Fluint16(state1), Fluint8((uint8)state2));
  }

  // Each of the bytes
  for (int i = 0; i < 256; i++) {
    state1++;
    state2 = LFSRNext<uint16_t, POLY>(state2);
    f(~state2, (uint8)state1, Fluint16(~state2), Fluint8((uint8)state1));
  }
};


static void RightShifts() {
  for (int i = 0; i < 65536; i++) {
    uint16_t u = i;
    Fluint16 x(u);
    CHECK_EQ16(Fluint16::RightShift<1>(x).ToInt(), u >> 1);
    CHECK_EQ16(Fluint16::RightShift<2>(x).ToInt(), u >> 2);
  }
}

static void LeftShifts() {
  for (int i = 0; i < 65536; i++) {
    uint16_t u = i;
    Fluint16 x(u);
    CHECK_EQ16(Fluint16::LeftShift<1>(x).ToInt(), u << 1);
    CHECK_EQ16(Fluint16::LeftShift<2>(x).ToInt(), u << 2);
  }
}

static void BitOps() {
  ManyPairs<LFSR_POLY1>([](
      uint16_t state1, uint16_t state2,
      Fluint16 f1, Fluint16 f2) {
    const uint16_t a = state1 & state2;
    const uint16_t o = state1 | state2;
    const uint16_t x = state1 ^ state2;

    const Fluint16 fa = f1 & f2;
    const Fluint16 fo = f1 | f2;
    const Fluint16 fx = f1 ^ f2;

    CHECK_EQ16(a, fa.ToInt());
    CHECK_EQ16(o, fo.ToInt());
    CHECK_EQ16(x, fx.ToInt());
    });
}

static void Addition() {
  ManyPairs<LFSR_POLY1>([](
      uint16_t state1, uint16_t state2,
      Fluint16 f1, Fluint16 f2) {
    const uint16_t sum = state1 + state2;

    const Fluint16 fsum = f1 + f2;

    CHECK_EQ16(sum, fsum.ToInt());
    });
}

static void Add8() {
  ManyPairs8<LFSR_POLY1>([](
      uint16_t state1, uint8_t state2,
      Fluint16 f1, Fluint8 f2) {
    const uint16_t sum = state1 + state2;

    const Fluint16 fsum = f1 + f2;
    Fluint16 fsum2 = f1;
    fsum2 += f2;

    CHECK_EQ16(sum, fsum.ToInt());
    CHECK_EQ16(sum, fsum2.ToInt());
    });
}


static void Increment() {
  for (int i = 0; i < 65536; i++) {
    {
      uint16_t u = i;
      Fluint16 f(u);

      uint16_t ou = u++;
      Fluint16 of = f++;
      CHECK_EQ(ou, of.ToInt());
      CHECK_EQ(u, f.ToInt());
    }

    {
      uint16_t u = i;
      Fluint16 f(u);

      uint16_t ou = ++u;
      Fluint16 of = ++f;
      CHECK_EQ(ou, of.ToInt());
      CHECK_EQ(u, f.ToInt());
    }
  }
}

static void Subtraction() {
  ManyPairs<LFSR_POLY1>([](
      uint16_t state1, uint16_t state2,
      Fluint16 f1, Fluint16 f2) {
    const uint16_t diff = state1 - state2;

    const Fluint16 fdiff = f1 - f2;

    CHECK_EQ16(diff, fdiff.ToInt());
    });
}

static void Decrement() {
  for (int i = 0; i < 65536; i++) {
    {
      uint16_t u = i;
      Fluint16 f(u);

      uint16_t ou = u--;
      Fluint16 of = f--;
      CHECK_EQ(ou, of.ToInt());
      CHECK_EQ(u, f.ToInt());
    }

    {
      uint16_t u = i;
      Fluint16 f(u);

      uint16_t ou = --u;
      Fluint16 of = --f;
      CHECK_EQ(ou, of.ToInt());
      CHECK_EQ(u, f.ToInt());
    }
  }
}

static void Zero() {
  for (int i = 0; i < 65536; i++) {
    {
      uint16_t u = i;
      Fluint16 f(u);

      Fluint8 z = Fluint16::IsZero(f);
      Fluint8 nz = Fluint16::IsntZero(f);
      CHECK(z.ToInt() != nz.ToInt());
      CHECK(z.ToInt() <= 1);
      CHECK(nz.ToInt() <= 1);
      if (u) CHECK(nz.ToInt() == 1);
      else CHECK(z.ToInt() == 1);
    }
  }
}


int main(int argc, char **argv) {

  LeftShifts(); printf("Left Shift OK\n");
  RightShifts(); printf("Right Shift OK\n");
  BitOps(); printf("Bitops OK\n");
  Addition(); printf("Addition OK\n");
  Increment(); printf("Increment OK\n");
  Subtraction(); printf("Subtraction OK\n");
  Decrement(); printf("Decrement OK\n");
  Zero(); printf("Zero OK\n");
  Add8(); printf("Add8 OK\n");

  printf("OK\n");
  return 0;
}
