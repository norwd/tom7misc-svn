
#include "modelinfo.h"

#include <memory>
#include <string>
#include <cstdint>
#include <cmath>

#include "image.h"
#include "network.h"
#include "base/logging.h"
#include "base/stringprintf.h"
#include "opt/opt.h"

using namespace std;
using uint32 = uint32_t;
using int64 = int64_t;

namespace {
// XXX could make sense as a standalone utility?
struct Histo {
  void Add(float f) {
    if (bound_low.has_value() && f < bound_low.value())
      f = bound_low.value();
    if (bound_high.has_value() && f > bound_high.value())
      f = bound_high.value();

    values.push_back(f);
  }

  Histo() {}
  Histo(optional<float> bound_low, optional<float> bound_high) :
    bound_low(bound_low), bound_high(bound_high) {}

  // Assumes width is the number of buckets you want.
  // If tallest_bucket is, say, 0.9, the bars are stretched to go 90%
  // of the way to the top of the image (1.0 is a sensible default but
  // can be confusing in the presence of tick marks, say).
  std::tuple<float, float, ImageA> MakeImage(
      int width, int height,
      double tallest_bucket = 1.0) const {
    ImageA img(width, height);

    float lo = std::numeric_limits<float>::infinity();
    float hi = -std::numeric_limits<float>::infinity();

    for (float v : values) {
      lo = std::min(v, lo);
      hi = std::max(v, hi);
    }

    const float ival = hi - lo;
    const float bucket_width = ival / width;
    const float oval = 1.0f / ival;
    if (ival <= 0) return make_tuple(lo, hi, img);

    vector<int64> count(width, 0);
    for (float v : values) {
      float f = (v - lo) * oval;
      int bucket = roundf(f * (width - 1));
      count[bucket]++;
    }

    int64 max_count = 0;
    int maxi = 0;
    for (int i = 0; i < count.size(); i++) {
      int64 c = count[i];
      if (c > max_count) {
        max_count = c;
        maxi = i;
      }
    }

    // Finally, fill in the image.
    for (int bucket = 0; bucket < width; bucket++) {
      double hfrac = count[bucket] / (double)max_count;
      float fh = (hfrac * tallest_bucket) * (height - 1);
      int h = fh;
      float fpart = fh - h;
      // don't allow zero pixels.
      // this is not accurate but I want to be able to see
      // non-empty buckets clearly
      if (h == 0 && count[bucket] > 0) {
        h = 1;
        fpart = 0.0f;
      }
      int nh = height - h;
      if (nh > 0) {
        uint8 v = roundf(fpart * 255);
        img.SetPixel(bucket, nh - 1, v);
      }
      for (int y = nh; y < height; y++) {
        CHECK(bucket < img.Width() && bucket >= 0 &&
              y < img.Height() && y >= 0) << bucket << " " << y;
        img.SetPixel(bucket, y, 0xFF);
      }
    }

    // Label the mode.
    float center = lo + ((maxi + 0.5f) * bucket_width);
    string label = StringPrintf("%.4f", center);
    int lw = label.size() * 9;
    // Align on left or right of the label so as not to run off the screen
    // (we could also try to avoid other buckets?)
    int x = maxi > (width / 2) ? maxi - (lw + 2) : maxi + 3;
    // Align with the peak, taking into account tallest_bucket.
    int y = (1.0 - tallest_bucket) * (height - 1);
    img.BlendText(x, y, 0xFF, label);

    return make_tuple(lo, hi, img);
  }

  // For example with tick=0.25, vertical lines at -0.25, 0, 0.25, 0.50, ...
  static ImageRGBA TickImage(int width, int height, float lo, float hi,
                             uint32 negative_tick_color,
                             uint32 zero_tick_color,
                             uint32 positive_tick_color,
                             float tick) {
    ImageRGBA img(width, height);
    const float ival = hi - lo;
    const float bucket_width = ival / width;

    for (int x = 0; x < width; x++) {
      const float bucket_lo = lo + x * bucket_width;
      const float bucket_hi = bucket_lo + bucket_width;
      // Does any tick edge reside in the bucket?
      // (Note there can be more than one...)

      // Floor here because we need rounding towards negative
      // infinity, not zero.
      const int tlo = floorf(bucket_lo / tick);
      const int thi = floorf(bucket_hi / tick);
      if (tlo != thi) {
        uint32 tick_color = 0;
        // tlo and thi are floor, so zero would fall in the bucket
        // [-1, 0]
        if (tlo == -1 && thi == 0) tick_color = zero_tick_color;
        else if (tlo < 0) tick_color = negative_tick_color;
        else tick_color = positive_tick_color;
        for (int y = 0; y < height; y++) {
          img.SetPixel32(x, y, tick_color);
        }
      }
    }
    return img;
  }

  // If present, samples are clipped to these bounds.
  optional<float> bound_low = nullopt, bound_high = nullopt;

  vector<float> values;
};
}  // namespace

static ImageRGBA ErrorImage(const string &message) {
  ImageRGBA error((message.size() + 2) * 18, 24);
  error.Clear32(0x000000FF);
  error.BlendText2x32(9, 3, 0xA00000FF, message);
  return error;
}

ImageRGBA ModelInfo::Histogram(
    const Network &net, int width, int height,
    std::optional<float> weight_bound_low,
    std::optional<float> weight_bound_high,
    std::optional<float> bias_bound_low,
    std::optional<float> bias_bound_high) {
  static constexpr int MARGIN = 8;
  // Get histogram of weights per layer.
  // Only layers with inputs (i.e. not the input layer) have weights.
  const int num_histos = net.layers.size() - 1;

  const int HISTOW = width / 2 - MARGIN;
  const int HISTOH = height / num_histos;

  ImageRGBA img{width, height};
  img.Clear32(0x000000FF);

  for (int layer_idx = 1; layer_idx < net.layers.size(); layer_idx++) {
    // Could draw each chunk as a different RGB channel, or could
    // partition the horizontal space?
    Histo bias_histo{bias_bound_low, bias_bound_high};
    Histo weight_histo{weight_bound_low, weight_bound_high};
    for (int chunk_idx = 0;
         chunk_idx < net.layers[layer_idx].chunks.size();
         chunk_idx++) {
      const Chunk &chunk = net.layers[layer_idx].chunks[chunk_idx];
      for (float f : chunk.biases) bias_histo.Add(f);
      for (float f : chunk.weights) weight_histo.Add(f);
    }

    auto DrawHisto = [&img](const Histo &histo, int x, int y, int w, int h,
                            const string &label) {
        constexpr int hmargin = 14;
        const auto [lo, hi, himg] = histo.MakeImage(w, h - hmargin, 0.95);
        ImageRGBA color(himg.Width(), himg.Height());

        // Always ticks at 0.1.
        const ImageRGBA timg = Histo::TickImage(w, h, lo, hi,
                                                0xFF777735,
                                                0xFFFFFF75,
                                                0x77FF7735,
                                                0.1f);
        color.BlendImage(0, 0, timg);

        // minor ticks if scale is very small?
        if (hi - lo <= 2.0f) {
          const ImageRGBA mtimg = Histo::TickImage(w, h, lo, hi,
                                                   0xFF777720,
                                                   0xFFFFFF40,
                                                   0x77FF7720,
                                                   0.025f);
          color.BlendImage(0, 0, mtimg);
        }

        color.BlendImage(0, 0, himg.AlphaMaskRGBA(0xFF, 0xFF, 0x00));

        img.BlendImage(x, y, color);
        const int label_y = y + (h - hmargin) + 1;
        img.BlendText32(x, label_y,
                        0xFFAAAAFF,
                        StringPrintf("%.9f", lo));
        const string his = StringPrintf("%.9f", hi);
        img.BlendText32(x + w - (his.size() * 9), label_y,
                        0xAAFFAAFF,
                        his);

        img.BlendText32(x + (w - (label.size() * 9)) / 2, label_y,
                        0xFFFFAAFF,
                        label);

      };

    DrawHisto(bias_histo,
              0, HISTOH * (layer_idx - 1), HISTOW, HISTOH,
              StringPrintf("^ Bias layer %d ^", layer_idx));

    DrawHisto(weight_histo,
              HISTOW + MARGIN, HISTOH * (layer_idx - 1), HISTOW, HISTOH,
              StringPrintf("^ Weights layer %d ^", layer_idx));
  }

  return img;
}

// Aesthetic mode
static inline uint32 GetWeightColor(float f, bool diagnostic_mode) {
  auto MapV = [diagnostic_mode](float f) -> uint8 {
      float ff = 0.0f;
      if (diagnostic_mode) {
        ff = sqrtf(f);
      } else {
        if (f > 0.25f) {
          ff = 1.0f;
        } else {
          ff = 4.0f * f;
          ff = sqrtf(sqrtf(ff));
        }
      }
      int v = roundf(255.0f * ff);
      if (v > 255) return 255;
      if (v < 0) return 0;
      return v;
    };

  // XXX from vacuum - configurable or remove?
  static constexpr float THRESHOLD = 0.00001f;
  if (diagnostic_mode && f == 0) {
    return 0xFFFF00FF;
  } else if (diagnostic_mode && fabs(f) < THRESHOLD) {
    return 0xFF00FFFF;
  } else if (f > 0) {
    uint32 v = MapV(f);
    return (v << 16) | 0xFF;
  } else {
    uint32 v = MapV(-f);
    return (v << 24) | 0xFF;
  }
}

// Find a good rectangle (returns width, height) to fit 'count'
// elements, each of which is w x h.
std::pair<int, int> MakeRectangle(int count, int w, int h) {

  const auto [bestw, bestv_] =
  Opt::Minimize1D([count, w, h](double rwidth) {
      int width = std::round(rwidth);
      CHECK(width >= 1);
      int height = count / width;
      while (width * height < count) height++;

      // Penalize empty cells.
      int empty = count - (width * height);

      // Penalize not square.
      int ww = width * w;
      int hh = height * h;
      double d = ww - hh;

      return empty * 10.0 + (d * d);
    },
  1.0,
  (double)count,
  10000);

  int resw = std::round(bestw);
  // Derive height the same way.
  int resh = count / resw;
  while (resw * resh < count) resh++;
  return make_pair(resw, resh);

#if 0
  CHECK(count >= 1);
  int w = ceilf(sqrtf(count));
  CHECK(w >= 1);
  int h = count / w;
  while (w * h < count) w++;
  return std::make_pair(w, h);
#endif
}

static ImageRGBA ChunkWeightsConvolution(
    const Chunk &chunk,
    int layer_idx,
    int chunk_idx,
    int prev_nodes,
    bool diagnostic_mode) {
  CHECK(chunk.type == CHUNK_CONVOLUTION_ARRAY);

  // For convolution layers, we show each feature in its 2D
  // orientation.
  const auto [features_across, features_down] =
    MakeRectangle(chunk.num_features,
                  chunk.pattern_width,
                  chunk.pattern_height);
  int pat_width = chunk.pattern_width;
  int pat_height = chunk.pattern_height;
  printf("For %d features of %d x %d, use %d across and %d down\n",
         chunk.num_features, pat_width, pat_height,
         features_across, features_down);

  const int ipn = pat_width * pat_height;
  CHECK(chunk.indices_per_node == ipn);

  const int padding = 2;
  int width = (pat_width + padding) * features_across;
  int height = (pat_height + padding) * features_down;

  ImageRGBA img(width, height);
  img.Clear32(0x000000FF);

  for (int fy = 0; fy < features_down; fy++) {
    for (int fx = 0; fx < features_across; fx++) {
      const int feature_num = fy * features_across + fx;
      const int weights_start = feature_num * ipn;

      if (feature_num < chunk.num_features) {
        int xpos = (pat_width + padding) * fx;
        int ypos = (pat_height + padding) * fy;

        for (int py = 0; py < pat_height; py++) {
          for (int px = 0; px < pat_width; px++) {
            const int pidx = py * pat_width + px;

            const float w =
              chunk.weights[weights_start + pidx];

            img.SetPixel32(xpos + px, ypos + py,
                           GetWeightColor(w, diagnostic_mode));
          }
        }
      }
    }
  }

  return img;
}

static ImageRGBA ChunkWeightsSparseDense(
    const Chunk &chunk,
    int layer_idx,
    int chunk_idx,
    int prev_nodes,
    bool diagnostic_mode) {

  if ((int64)chunk.num_nodes * (int64)prev_nodes > (int64)5000000) {
    return ErrorImage(StringPrintf("Too big! %d x %d",
                                   chunk.num_nodes, prev_nodes));
  }

  const int ipn = chunk.indices_per_node;

  // TODO: Consider only rendering the span. If we're using the full
  // width, we could at least highlight the span (or nodes outside the
  // span).

  //     <--- this layer's nodes --->
  // ^
  // |
  // weights
  // |
  // |
  // v
  constexpr int TOP = 12;
  constexpr int LEFT = 12;
  constexpr int RIGHT = 4;
  constexpr int BOTTOM = 4;
  ImageRGBA img(LEFT + RIGHT + chunk.num_nodes,
                TOP + BOTTOM + prev_nodes);

  const uint32 missing_weight_color =
    diagnostic_mode ? 0x000060FF : 0x000000FF;

  constexpr int WEIGHTSX = LEFT;
  constexpr int WEIGHTSY = TOP;
  img.Clear32(0x000000FF);
  img.BlendRect32(WEIGHTSX, WEIGHTSY, chunk.num_nodes, prev_nodes,
                  missing_weight_color);

  vector<bool> has_reference(prev_nodes, false);

  for (int x = 0; x < chunk.num_nodes; x++) {
    const int start = x * ipn;
    for (int i = 0; i < ipn; i++) {
      uint32 idx =
        chunk.type == CHUNK_SPARSE ?
        chunk.indices[start + i] :
        chunk.span_start + i;

      float w = chunk.weights[start + i];

      CHECK(idx < prev_nodes);
      has_reference[idx] = true;
      uint32 color = GetWeightColor(w, diagnostic_mode);
      img.SetPixel32(WEIGHTSX + x, WEIGHTSY + idx, color);
    }
  }

  if (diagnostic_mode) {
    for (int y = 0; y < prev_nodes; y++) {
      if (!has_reference[y]) {
        for (int x = 0; x < chunk.num_nodes; x++) {
          // XXX: Configurable?
          // Would be missing_weight_color if not detected.
          img.SetPixel32(WEIGHTSX + x, WEIGHTSY + y, 0xFFFF00FF);
        }
      }
    }
  }

  img.BlendText32(LEFT, 2, 0xCCCCCCFF,
                  StringPrintf("<--   Chunk %d.%d's nodes (%d)  ipn=%d  -->",
                               layer_idx, chunk_idx, chunk.num_nodes, ipn));

  return img;
}


ImageRGBA ModelInfo::ChunkWeights(const Network &net,
                                  int layer_idx,
                                  int chunk_idx,
                                  bool diagnostic_mode) {
  CHECK(layer_idx > 0 && layer_idx < net.layers.size());
  const Layer &layer = net.layers[layer_idx];
  CHECK(chunk_idx >= 0 && chunk_idx < layer.chunks.size());
  const Chunk &chunk = layer.chunks[chunk_idx];

  const int prev_nodes = net.layers[layer_idx - 1].num_nodes;

  if (chunk.type == CHUNK_CONVOLUTION_ARRAY) {
    /*
    return ErrorImage(StringPrintf("Convolution %dx%d, unimplemented",
                                   chunk.pattern_width, chunk.pattern_height));
    */
    return ChunkWeightsConvolution(
        chunk, layer_idx, chunk_idx, prev_nodes, diagnostic_mode);
  } else {
    return ChunkWeightsSparseDense(
        chunk, layer_idx, chunk_idx, prev_nodes, diagnostic_mode);
  }
}
