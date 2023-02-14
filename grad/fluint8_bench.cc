
#include <string>
#include <cstdint>

#include "timer.h"
#include "fluint8.h"

#include "base/stringprintf.h"
#include "ansi.h"

struct Benchmark {
  Benchmark(string name, int num_iters = 1000) :
    name(name), num_iters(num_iters) {}
  string name;

  int64 ops = 0;
  int out = 0;
  const int num_iters = 1000;

  inline void Op() {
    ops++;
  }
  inline void Observe(const Fluint8 f) {
    out += f.ToInt();
  }

  template<class F>
  void Run(const F &f) {
    printf("Bench " AWHITE("%s") ": ", name.c_str());
    Timer timer;
    for (int iters = 0; iters < num_iters; iters++) {
      f();
      if (iters % 100 == 0) printf(".");
    }
    double sec = timer.Seconds();
    printf(ANSI_RESTART_LINE);
    double os = ops / sec;
    string oss;
    if (os > 1000000.0) {
      oss = StringPrintf("%.2f " AYELLOW("M") "/sec",
                         os / 1000000.0);
    } else if (os > 1000.0) {
      oss = StringPrintf("%.2f " APURPLE("K") "/sec",
                         os / 1000.0);
    }
    printf(ABLUE("%s") ": %lld ops in %.3f" AGREY("s") " = %s\n",
           name.c_str(), ops, sec, oss.c_str());

  }
};

[[maybe_unused]]
static void BenchEq() {
  Benchmark bench("eq");

  bench.Run([&]() {
      for (int i = 0; i < 256; i++) {
        Fluint8 ii(i);
        for (int j = 0; j < 256; j++) {
          Fluint8 jj(j);

          Fluint8 c = Fluint8::Eq(ii, jj);
          bench.Op();
          bench.Observe(c);
        }
      }
    });
}

template<int N>
static void BenchRightShift() {
  Benchmark bench(StringPrintf("RightShift<%d>", N));

  bench.Run([&]() {
    for (int sub_iters = 0; sub_iters < 256; sub_iters++)
      for (int i = 0; i < 256; i++) {
        Fluint8 ii(i);

        Fluint8 c = Fluint8::RightShift<N>(ii);
        bench.Op();
        bench.Observe(c);
      }
    });
}

// Note: Can improve shifts by 5, 6
static void BenchRightShifts() {
  BenchRightShift<0>();
  BenchRightShift<1>();
  BenchRightShift<2>();
  BenchRightShift<3>();
  BenchRightShift<4>();
  BenchRightShift<5>();
  BenchRightShift<6>();
  BenchRightShift<7>();
}

static void BenchPlus() {
  Benchmark bench("Plus", 2000);

  bench.Run([&]() {
      for (int i = 0; i < 256; i++) {
        for (int j = 0; j < 256; j++) {
          Fluint8 ii(i), jj(j);

          Fluint8 c = ii + jj;
          bench.Op();
          bench.Observe(c);
        }
      }
    });
}

static void BenchAddWithCarry() {
  Benchmark bench("AddWithCarry", 2000);

  bench.Run([&]() {
      for (int i = 0; i < 256; i++) {
        for (int j = 0; j < 256; j++) {
          Fluint8 ii(i), jj(j);

          const auto &[c, d] = Fluint8::AddWithCarry(ii, jj);
          bench.Op();
          bench.Observe(c);
          bench.Observe(d);
        }
      }
    });
}

static void BenchAnd() {
  Benchmark bench("BitwiseAnd", 2000);

  bench.Run([&]() {
      for (int i = 0; i < 256; i++) {
        for (int j = 0; j < 256; j++) {
          Fluint8 ii(i), jj(j);

          Fluint8 c = ii & jj;
          bench.Op();
          bench.Observe(c);
        }
      }
    });
}

static void BenchIf() {
  Benchmark bench("If", 1000);

  bench.Run([&]() {
      for (int i = 0; i < 256; i++) {
        Fluint8 ii(i & 1);
        for (int j = 0; j < 256; j++) {
          Fluint8 jj(j);

          Fluint8 c = Fluint8::If(ii, jj);
          bench.Op();
          bench.Observe(c);
        }
      }
    });
}


int main(int argc, char **argv) {
  AnsiInit();

  /*
  BenchEq();
  BenchRightShifts();
  BenchPlus();
  BenchAddWithCarry();
  BenchAnd();
  */
  BenchIf();

  return 0;
}
