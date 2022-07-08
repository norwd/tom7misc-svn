#include "expression.h"
#include "timer.h"

#include <algorithm>
#include <functional>
#include <array>
#include <utility>

#include "half.h"

#include "grad-util.h"
#include "makefn-ops.h"
#include "arcfour.h"
#include "randutil.h"
#include "color-util.h"

using Table = Exp::Table;
using uint32 = uint32_t;
using uint8 = uint8_t;

static Table MakeTableFromFn(const std::function<half(half)> &f) {
  Table table;
  for (int i = 0; i < 65536; i++) {
    half x = Exp::GetHalf((uint16)i);
    half y = f(x);
    table[i] = Exp::GetU16(y);
  }
  return table;
}

template<size_t N>
auto Sample(ArcFour *rc,
            const std::array<std::pair<int, int>, N> &bounds) ->
  // in [0, 1], then integer in bounds
  std::pair<std::array<double, N>, std::array<int, N>> {
  std::array<double, N> ret_norm;
  std::array<int, N> ret_scaled;
  for (int i = 0; i < N; i++) {
    // Note int bounds are [low, high).
    const auto [low, high] = bounds[i];
    int offset = RandTo(rc, high - low);
    ret_scaled[i] = low + offset;
    ret_norm[i] = offset / (double)(high - low - 1);
  }
  return make_pair(ret_norm, ret_scaled);
}

template<size_t N>
auto Sample(ArcFour *rc,
            const std::array<std::pair<double, double>, N> &bounds) ->
  // in [0, 1], then scaled to bounds
  std::pair<std::array<double, N>, std::array<double, N>> {
  std::array<double, N> ret_norm;
  std::array<double, N> ret_scaled;
  for (int i = 0; i < N; i++) {
    double f = RandDouble(rc);
    ret_norm[i] = f;
    const auto [low, high] = bounds[i];
    ret_scaled[i] = low + (high - low) * f;
  }
  return make_pair(ret_norm, ret_scaled);
}

static uint32 MixRGB(float r, float g, float b, float a) {
  uint32 rr = std::clamp((int)(r * 255.0f), 0, 255);
  uint32 gg = std::clamp((int)(g * 255.0f), 0, 255);
  uint32 bb = std::clamp((int)(b * 255.0f), 0, 255);
  uint32 aa = std::clamp((int)(a * 255.0f), 0, 255);
  return (rr << 24) | (gg << 16) | (bb << 8) | aa;
}

static constexpr int IMAGE_SIZE = 1920;

static bool IsZero(const Table &table) {
  double total_sum = 0.0;
  double total_width = 2.0;
  for (half pos = (half)-1.0; pos < (half)1.0; /* in loop */) {
    half next = nextafter(pos, (half)1.0);
    uint16 upos = Exp::GetU16(pos);
    double err = Exp::GetHalf(table[upos]);
    double width = (double)next - (double)pos;
    total_sum += fabs(err) * width;
    pos = next;
  }

  return (total_sum / total_width) < (half)0.001;
}

struct Stats {
  Stats() {}

  int64 samples_in = 0;
  int64 samples_out = 0;
  int iszero = 0;
  int64 denom = 0;

  void Accumulate(const Table &result) {
    denom++;
    for (int s = 0; s < 256; s++) {
      half x = (half)( (s / 128.0) - 1.0);
      half y = Exp::GetHalf(result[Exp::GetU16(x)]);
      if (y < (half)-1.0 || y > (half)1.0) {
        samples_out++;
      } else {
        samples_in++;
      }
    }

    if (IsZero(result)) iszero++;
  }

  void Report() {
    int64 total_samples = samples_out + samples_in;
    printf("Samples in: %.3f%% out: %.3f%%\n"
           "Zero: %d/%d (%.3f%%)\n",
           (samples_in * 100.0) / total_samples,
           (samples_out * 100.0) / total_samples,
           iszero, denom, (iszero * 100.0) / denom);
  }

};

static void PlotOp2() {
  ArcFour rc("op2");
  // Not used by op2.
  Table target =
    MakeTableFromFn([](half x) {
        return sin(x * (half)3.141592653589);
      });

  ImageRGBA img(IMAGE_SIZE, IMAGE_SIZE);
  img.Clear32(0x000000FF);
  GradUtil::Grid(&img);

  static constexpr int SAMPLES = 5000;
  Stats stats;
  for (int i = 0; i < SAMPLES; i++) {
    if (i % 1000 == 0) printf("%d/%d\n", i, SAMPLES);
    const auto [norm, scaled] = Sample(&rc, Op2::DOUBLE_BOUNDS);
    auto [rf, gf, bf] = norm;

    const uint32 color =
      MixRGB(rf * 0.9 + 0.1,
             gf * 0.9 + 0.1,
             bf * 0.9 + 0.1, 0.20);

    Exp::Allocator alloc;
    const Exp *exp = Op2::GetExp(&alloc, {}, scaled, target);

    Table result = Exp::TabulateExpression(exp);
    stats.Accumulate(result);

    GradUtil::Graph(result, color, &img);
  }
  stats.Report();
  img.Save("op2.png");
}

static double StrobeOffset(std::pair<double, double> bound,
                           int s, int num_strobe, double frac) {
  auto [low, high] = bound;

  // full width of strobe
  double w = (high - low) * frac;
  // half width of strobe (+/-)
  double hw = w * 0.5;
  // magnitude per strobe sample (on each side)
  double m = hw / num_strobe;

  double r = m * (s >> 1);
  return (s & 1) ? r : -r;
}

static int StrobeOffsetI(int s) {
  return (s & 1) ? (s >> 1) : -(s >> 1);
}

static void PlotOp3() {
  ArcFour rc("op3");
  // Not used by op3.
  Table target =
    MakeTableFromFn([](half x) {
        return sin(x * (half)3.141592653589);
      });

  ImageRGBA img(IMAGE_SIZE, IMAGE_SIZE);
  img.Clear32(0x000000FF);
  GradUtil::Grid(&img);

  static constexpr int STROBE = 100;
  static constexpr double STROBE_FRAC = 0.1;
  static constexpr int SAMPLES = 50;
  static constexpr int STROBE_DIM = 0;
  Stats stats;
  for (int i = 0; i < SAMPLES; i++) {
    if (i % 1000 == 0) printf("%d/%d\n", i, SAMPLES);
    const auto [inorm, iscaled] = Sample(&rc, Op3::INT_BOUNDS);
    const auto [dnorm, dscaled] = Sample(&rc, Op3::DOUBLE_BOUNDS);
    auto [iterf] = inorm;
    auto [xf, yf, zf, wf] = dnorm;

    for (int s = 0; s < STROBE; s++) {
      auto iscaledo = iscaled;
      auto dscaledo = dscaled;
      /*
      double d = StrobeOffset(Op3::DOUBLE_BOUNDS[STROBE_DIM],
                              s, STROBE, STROBE_FRAC);
      dscaledo[STROBE_DIM] += d;
      */

      int d = StrobeOffsetI(s);
      iscaledo[STROBE_DIM] += d;
      if (iscaledo[STROBE_DIM] < Op3::INT_BOUNDS[STROBE_DIM].first ||
          iscaledo[STROBE_DIM] >= Op3::INT_BOUNDS[STROBE_DIM].second)
        continue;

      // pretty arbitrary 5D -> 3D reduction
      const auto [r1, g1, b1] =
        ColorUtil::HSVToRGB(iterf, 0.1 + xf * 0.9, 0.1 + yf * 0.9);
      const auto [r2, g2, b2] =
        ColorUtil::HSVToRGB(zf, 0.1 + wf * 0.9, 1.0);

      const float alpha = 0.05 + 0.15 * ((STROBE - s) / (float)STROBE);

      const uint32 color =
        MixRGB((r1 + r2) * 0.5, (g1 + g2) * 0.5, (b1 + b2) * 0.5,
               alpha);

      Exp::Allocator alloc;
      const Exp *exp = Op3::GetExp(&alloc, iscaledo, dscaledo, target);

      Table result = Exp::TabulateExpression(exp);
      stats.Accumulate(result);

      GradUtil::Graph(result, color, &img);
    }
  }
  stats.Report();
  img.Save("op3.png");
}

static void PlotOp4() {
  ArcFour rc("op4");

  Table target =
    MakeTableFromFn([](half x) {
        return sin(x * (half)3.141592653589);
      });

  ImageRGBA img(IMAGE_SIZE, IMAGE_SIZE);
  img.Clear32(0x000000FF);
  GradUtil::Grid(&img);

  static constexpr int STROBE = 20;
  static constexpr double STROBE_FRAC = 0.1;
  static constexpr int SAMPLES = 50;
  static constexpr int STROBE_DIM = 1;
  Stats stats;
  for (int i = 0; i < SAMPLES; i++) {
    if (i % 1000 == 0) printf("%d/%d\n", i, SAMPLES);
    const auto [inorm, iscaled] = Sample(&rc, Op4::INT_BOUNDS);
    const auto [dnorm, dscaled] = Sample(&rc, Op4::DOUBLE_BOUNDS);
    auto [a, b, c] = inorm;
    auto [d, e] = dnorm;

    for (int s = 0; s < STROBE; s++) {
      auto iscaledo = iscaled;
      auto dscaledo = dscaled;

      double d = StrobeOffset(Op4::DOUBLE_BOUNDS[STROBE_DIM],
                              s, STROBE, STROBE_FRAC);
      dscaledo[STROBE_DIM] += d;

      /*
      int d = StrobeOffsetI(s);
      iscaledo[STROBE_DIM] += d;
      if (iscaledo[STROBE_DIM] < Op4::INT_BOUNDS[STROBE_DIM].first ||
          iscaledo[STROBE_DIM] >= Op4::INT_BOUNDS[STROBE_DIM].second)
        continue;
      */

      // pretty arbitrary 5D -> 3D reduction
      const auto [r1, g1, b1] =
        ColorUtil::HSVToRGB(a, 0.1 + b * 0.9, 0.1 + c * 0.9);
      const auto [r2, g2, b2] =
        ColorUtil::HSVToRGB(d, 0.1 + e * 0.9, 1.0);

      const float alpha = 0.05 + 0.15 * ((STROBE - s) / (float)STROBE);

      const uint32 color =
        MixRGB((r1 + r2) * 0.5, (g1 + g2) * 0.5, (b1 + b2) * 0.5,
               alpha);

      Exp::Allocator alloc;
      const Exp *exp = Op4::GetExp(&alloc, iscaledo, dscaledo, target);

      Table result = Exp::TabulateExpression(exp);
      stats.Accumulate(result);

      GradUtil::Graph(result, color, &img);
    }
  }
  stats.Report();
  img.Save("op4.png");
}

static void PlotOp5() {
  ArcFour rc("op5");

  Table target =
    MakeTableFromFn([](half x) {
        return sin(x * (half)3.141592653589);
      });

  ImageRGBA img(IMAGE_SIZE, IMAGE_SIZE);
  img.Clear32(0x000000FF);
  GradUtil::Grid(&img);

  static constexpr int STROBE = 50;
  static constexpr double STROBE_FRAC = 0.25;
  static constexpr int SAMPLES = 50;
  static constexpr int STROBE_DIM = 1;
  Stats stats;
  for (int i = 0; i < SAMPLES; i++) {
    if (i % 1000 == 0) printf("%d/%d\n", i, SAMPLES);
    const auto [inorm, iscaled] = Sample(&rc, Op5::INT_BOUNDS);
    const auto [dnorm, dscaled] = Sample(&rc, Op5::DOUBLE_BOUNDS);
    auto [a, b] = inorm;
    auto [c, d, e] = dnorm;

    for (int s = 0; s < STROBE; s++) {
      auto iscaledo = iscaled;
      auto dscaledo = dscaled;

      double d = StrobeOffset(Op5::DOUBLE_BOUNDS[STROBE_DIM],
                              s, STROBE, STROBE_FRAC);
      dscaledo[STROBE_DIM] += d;

      /*
      int d = StrobeOffsetI(s);
      iscaledo[STROBE_DIM] += d;
      if (iscaledo[STROBE_DIM] < Op5::INT_BOUNDS[STROBE_DIM].first ||
          iscaledo[STROBE_DIM] >= Op5::INT_BOUNDS[STROBE_DIM].second)
        continue;
      */

      // pretty arbitrary 5D -> 3D reduction
      const auto [r1, g1, b1] =
        ColorUtil::HSVToRGB(a, 0.1 + b * 0.9, 0.1 + c * 0.9);
      const auto [r2, g2, b2] =
        ColorUtil::HSVToRGB(d, 0.1 + e * 0.9, 1.0);

      const float alpha = 0.05 + 0.15 * ((STROBE - s) / (float)STROBE);

      const uint32 color =
        MixRGB((r1 + r2) * 0.5, (g1 + g2) * 0.5, (b1 + b2) * 0.5,
               alpha);

      Exp::Allocator alloc;
      const Exp *exp = Op5::GetExp(&alloc, iscaledo, dscaledo, target);

      Table result = Exp::TabulateExpression(exp);
      stats.Accumulate(result);

      GradUtil::Graph(result, color, &img);
    }
  }
  stats.Report();
  img.Save("op5.png");
}


int main(int argc, char **argv) {

  PlotOp2();
  PlotOp3();
  PlotOp4();
  PlotOp5();
  return 0;
}
