
#include "image.h"

#include <string>
#include <vector>
#include <cstdint>
#include <utility>
#include <cstring>
#include <tuple>
#include <algorithm>

#include "lines.h"
#include "stb_image.h"
#include "stb_image_write.h"
#include "base/logging.h"

using namespace std;
using uint8 = uint8_t;
using uint32 = uint32_t;

// Performance assumptions. You can just define this to
// do nothing (last case) if you're having trouble compiling.
// TODO: I think there is an assertion coming in the standard?
#if defined __has_builtin
# if __has_builtin (__builtin_unreachable)
# else
#  define __builtin_unreachable() do { } while (0)
# endif
#else
# define __builtin_unreachable() do { } while (0)
#endif

namespace {
// Generated by bit7/embedcc.cc.
struct EmbeddedFont {
  static constexpr int WIDTH = 9;
  static constexpr int HEIGHT = 9;

  // No bounds checking!
  static constexpr bool GetBit(int c, int x, int y) {
    int b = c * (WIDTH * HEIGHT) + y * WIDTH + x;
    return 0 != (words[b >> 6] & (1ULL << (63 - (b & 63))));
  }

  // SetPixel and ClearPixel called like SetPixel(x, y).
  template<class FS, class FC>
  static void Blit(int c, int x, int y,
                   FS SetPixel, FC ClearPixel = [](int, int){}) {
    if (c < 0 || c >= 128) return;
    for (int sy = 0; sy < HEIGHT; sy++) {
      for (int sx = 0; sx < WIDTH; sx++) {
        if (GetBit(c, sx, sy)) {
          SetPixel(x + sx, y + sy);
        } else {
          ClearPixel(x + sx, y + sy);
        }
      }
    }
  }

 private:
  static constexpr uint64_t words[162] = {
  0x0000000000000000ULL, 0x00003f3fdb6ff7fbULL, 0x6dce7e001f9fedb7ULL,
  0xfbfdcee73f000fcfULL, 0xf6dbfdcec36d9f80ULL, 0x000333fdfeff3f0fULL,
  0x030000c0f0fcff7fULL, 0x9f87818000003024ULL, 0x2110848180000000ULL,
  0x000c0f0781800000ULL, 0x007e41a0d068341bULL, 0xfcfc003f20d56935ULL,
  0x5a0dfe7e00020100ULL, 0x8043e01008040201ULL, 0x008040201f080402ULL,
  0x010000000000ff84ULL, 0x0201008040201008ULL, 0x7fc0000000002010ULL,
  0x08043fe100804020ULL, 0x1008040201008040ULL, 0x2010080402010f80ULL,
  0x0000000000000000ULL, 0x007c201008040000ULL, 0x000003e010080402ULL,
  0x01008040201f0000ULL, 0x00000000000000ffULL, 0x8000000007ffffffULL,
  0xffffc00000000000ULL, 0x000000001fffffffULL, 0xfff0f87c3e1f0f87ULL,
  0xc3e1f00783c1e0f0ULL, 0x783c1e0f07ffffffULL, 0xffffffffffffffffULL,
  0xfb6ffffedbffffb6ULL, 0xffedbdb36f6d6ddbULL, 0x5b5b76d491276d24ULL,
  0x495b493246dc9093ULL, 0x2424a48929292249ULL, 0x244992491244924cULL,
  0x9124002480000920ULL, 0x0002480000000000ULL, 0x0000000000000c06ULL,
  0x030180c000301800ULL, 0x198cc66000000000ULL, 0x00000cc667f999feULL,
  0x6633000001b3fb60ULL, 0xfc1b7f3600000702ULL, 0x999c1c1c1ccca070ULL,
  0x01f98c60773599c7ULL, 0xc0000030180c0000ULL, 0x00000000000c0c0cULL,
  0x060300c030000018ULL, 0x060180c060606000ULL, 0x0000030db7f8f0ccULL,
  0x660000000180c1f8ULL, 0xfc180c0000000000ULL, 0x0000180c06010000ULL,
  0x0000007e3f000000ULL, 0x0000000000000006ULL, 0x0300000006070707ULL,
  0x0707030000003c33ULL, 0x33db6f3330f00000ULL, 0x1c1e0b0180c060fcULL,
  0x00000f0cc060f1e0ULL, 0xc07e00000f806030ULL, 0x700c063e00000063ULL,
  0x3198fc0603018000ULL, 0x03f180f80603018fULL, 0x800000f8cc603f18ULL,
  0xcc63e00001fcc603ULL, 0x03030180c000007eULL, 0x61b0cfcc3619f800ULL,
  0x003fb0d867f0180cULL, 0x06000000000180c0ULL, 0x0030180000000000ULL,
  0xc06000180c020000ULL, 0x6060606018060180ULL, 0x00000000fc00001fULL,
  0x8000000180601806ULL, 0x0606060000007c63ULL, 0x3183818000603000ULL,
  0x7e61b6db6ce300f8ULL, 0x00001e0f0cc667fbULL, 0x0d8600003f986c37ULL,
  0xf30d86fe000007e6ULL, 0x36030180631f8000ULL, 0x0fc6330d86c3633fULL,
  0x000007fb0d80f060ULL, 0x30dfe00003fd86c0ULL, 0x7830180c0000007cULL,
  0x636031d86631f800ULL, 0x00c361b0dfec361bULL, 0x0c00003f06030180ULL,
  0xc060fc000007c0c0ULL, 0x6033199878000018ULL, 0xcc6663e19cc76180ULL,
  0x0006030180c06030ULL, 0x9fc000061b9dfedbULL, 0x6db6db6000038de6ULL,
  0xf36db3d9ec700000ULL, 0x786661b0d86661e0ULL, 0x0000fe61b0d86fe6ULL,
  0x030000001e19986cULL, 0x3679987e00003f98ULL, 0x6c361bf998c60000ULL,
  0x0fec3601f806c37fULL, 0x00000ff0c0603018ULL, 0x0c060000061b0d86ULL,
  0xc361b0cfc000030dULL, 0x8666330f07818000ULL, 0x0186c36db6cb47e3ULL,
  0x300000c3739f8307ULL, 0xe73b0c000061b0ccULL, 0xc7e0c0603000003fULL,
  0xc06061e18180ff00ULL, 0x0007c30180c06030ULL, 0x1f0000060380e038ULL,
  0x0e0380c00003e030ULL, 0x180c06030f800000ULL, 0x6078662100000000ULL,
  0x000000000000000fULL, 0xf7f80000380e0380ULL, 0x0000000000000000ULL,
  0x0f8063f318fe0000ULL, 0x180c07e3198cc67eULL, 0x0000000001f180c0ULL,
  0x601f0000003018fcULL, 0xc663318fc0000000ULL, 0x00fcc37f300f8000ULL,
  0x00786633181f0603ULL, 0x00000000001f98ccULL, 0x63f019f80060301bULL,
  0x8e6633198c00000cULL, 0x00070180c060fc00ULL, 0x00018000e030180cULL,
  0x461e000c060331b8ULL, 0xf8663300000380c0ULL, 0x6030180c1f800000ULL,
  0x00016cff6db6d860ULL, 0x00000000dc7b3998ULL, 0xcc60000000003f30ULL,
  0xd86c33f000000000ULL, 0x3f986c37f3018000ULL, 0x00000fcc6631f80cULL,
  0x070000000bc73319ULL, 0x80c00000000003f3ULL, 0x01fc067e00000180ULL,
  0xc1f830180c038000ULL, 0x0000018cc663338fULL, 0x4000000000c36199ULL,
  0x8cc3c00000000061ULL, 0xb6db67e330000000ULL, 0x00198fc183f19800ULL,
  0x0000000c6631f80cULL, 0x863e00000007e070ULL, 0x60e07e000001c180ULL,
  0xc030301807000001ULL, 0x80c06030180c0600ULL, 0x0003806030300c06ULL,
  0x0e0000000000706dULL, 0xb6c1c00000000c1eULL, 0x1d0885cee7000000ULL,
  };
};
}

inline static constexpr std::tuple<uint8, uint8, uint8, uint8>
Unpack32(uint32 color) {
  return {(uint8)((color >> 24) & 255),
          (uint8)((color >> 16) & 255),
          (uint8)((color >> 8) & 255),
          (uint8)(color & 255)};
}

inline static constexpr uint32 Pack32(uint8 r, uint8 g, uint8 b, uint8 a) {
  return
    ((uint32)r << 24) | ((uint32)g << 16) | ((uint32)b << 8) | (uint32)a;
}

// TODO: Duplicate code between the different Load routines..

ImageRGBA::ImageRGBA(const std::vector<uint32> &rgba32,
                     int width, int height) :
  width(width), height(height), rgba(rgba32) {
  CHECK((int)rgba.size() == width * height);
}

// static
ImageRGBA *ImageRGBA::Load(const string &filename) {
  vector<uint8> ret;
  int width = 0, height = 0, bpp_unused = 0;
  uint8 *stb_rgba = stbi_load(filename.c_str(),
                              &width, &height, &bpp_unused, 4);
  if (stb_rgba == nullptr) return nullptr;
  const int bytes = width * height * 4;
  ret.resize(bytes);
  // TODO: Is this portable (or even correct) wrt to endianness?
  memcpy(ret.data(), stb_rgba, bytes);
  stbi_image_free(stb_rgba);
  return new ImageRGBA(std::move(ret), width, height);
}

ImageRGBA *ImageRGBA::LoadFromMemory(const char *data, size_t size) {
  vector<uint8> ret;
  int width = 0, height = 0, bpp_unused = 0;
  uint8 *stb_rgba = stbi_load_from_memory(
      (const stbi_uc*)data, size,
      &width, &height, &bpp_unused, 4);
  if (stb_rgba == nullptr) return nullptr;
  const int bytes = width * height * 4;
  ret.resize(bytes);
  // TODO: Is this portable (or even correct) wrt to endianness?
  memcpy(ret.data(), stb_rgba, bytes);
  stbi_image_free(stb_rgba);
  return new ImageRGBA(std::move(ret), width, height);
}

ImageRGBA *ImageRGBA::LoadFromMemory(const vector<uint8> &filebytes) {
  return LoadFromMemory((const char *)filebytes.data(), filebytes.size());
}

ImageRGBA::ImageRGBA(const vector<uint8> &rgba8, int width, int height)
  : width(width), height(height) {
  CHECK((int)rgba8.size() == width * height * 4);
  rgba.resize(width * height);
  for (int i = 0; i < width * height; i++) {
    uint8 r = rgba8[i * 4 + 0];
    uint8 g = rgba8[i * 4 + 1];
    uint8 b = rgba8[i * 4 + 2];
    uint8 a = rgba8[i * 4 + 3];
    rgba[i] = Pack32(r, g, b, a);
  }
}

ImageRGBA::ImageRGBA(int width, int height)
  : width(width), height(height), rgba(width * height) {
  Clear32(0x00000000);
}

std::vector<uint8_t> ImageRGBA::ToBuffer8() const {
  std::vector<uint8_t> ret(width * height * 4);
  for (int i = 0; i < width * height; i++) {
    CHECK(i < (int)rgba.size());
    const auto &[r, g, b, a] = Unpack32(rgba[i]);
    CHECK(i * 4 + 3 < (int)ret.size());
    ret[i * 4 + 0] = r;
    ret[i * 4 + 1] = g;
    ret[i * 4 + 2] = b;
    ret[i * 4 + 3] = a;
  }
  return ret;
}

bool ImageRGBA::Save(const std::string &filename) const {
  std::vector<uint8> buffer = ToBuffer8();
  CHECK((int)buffer.size() == width * height * 4);
  return !!stbi_write_png(filename.c_str(),
                          width, height, 4, buffer.data(), 4 * width);
}

vector<uint8> ImageRGBA::SaveToVec() const {
  std::vector<uint8> buffer = ToBuffer8();
  CHECK((int)buffer.size() == width * height * 4);
  return stbi_make_png_rgba(width, height, buffer.data());
}

string ImageRGBA::SaveToString() const {
  std::vector<uint8> buffer = ToBuffer8();
  CHECK((int)buffer.size() == width * height * 4);
  const vector<uint8> v = stbi_make_png_rgba(width, height, buffer.data());
  string ret;
  ret.resize(v.size());
  memcpy(ret.data(), v.data(), v.size());
  return ret;
}

bool ImageRGBA::SaveJPG(const std::string &filename, int quality) const {
  std::vector<uint8> buffer = ToBuffer8();
  CHECK((int)buffer.size() == width * height * 4);
  CHECK(quality >= 0 && quality <= 100) << quality;
  return !!stbi_write_jpg(filename.c_str(),
                          width, height, 4, buffer.data(), quality);
}

ImageRGBA *ImageRGBA::Copy() const {
  return new ImageRGBA(rgba, width, height);
}

bool ImageRGBA::operator==(const ImageRGBA &other) const {
  return other.Width() == Width() &&
    other.Height() == Height() &&
    other.rgba == rgba;
}

std::size_t ImageRGBA::Hash() const {
  uint64_t h = width * 31337;
  for (uint32_t v : rgba) {
    h = (h << 13) | (h >> (64 - 13));
    h ^= v;
    h *= 65537;
  }
  return (std::size_t)h;
}

ImageRGBA ImageRGBA::Crop32(int x, int y, int w, int h,
                            uint32 fill_color) const {
  CHECK(w > 0 && h > 0) << w << " " << h;
  const auto [r, g, b, a] = Unpack32(fill_color);

  ImageRGBA ret{w, h};
  // xx,yy in new image's coordinates
  for (int yy = 0; yy < h; yy++) {
    const int sy = yy + y;
    for (int xx = 0; xx < w; xx++) {
      const int sx = xx + x;
      // Fill color if not in bounds
      uint8 sr = r, sg = g, sb = b, sa = a;
      if (sx >= 0 && sx < width &&
          sy >= 0 && sy < height) {
        std::tie(sr, sg, sb, sa) = GetPixel(sx, sy);
      }

      ret.SetPixel(xx, yy, sr, sg, sb, sa);
    }
  }
  return ret;
}

ImageRGBA ImageRGBA::ScaleBy(int scale) const {
  // 1 is not useful, but it does work
  CHECK(scale >= 1);
  ImageRGBA ret(width * scale, height * scale);
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      const uint32 color = GetPixel32(x, y);
      for (int yy = 0; yy < scale; yy++) {
        for (int xx = 0; xx < scale; xx++) {
          ret.SetPixel32(x * scale + xx,
                         y * scale + yy,
                         color);
        }
      }
    }
  }
  return ret;
}

ImageRGBA ImageRGBA::ScaleDownBy(int scale) const {
  // 1 is not useful, but it does work
  CHECK(scale >= 1);
  const int ww = width / scale;
  const int hh = height / scale;
  ImageRGBA ret(ww, hh);
  for (int y = 0; y < hh; y++) {
    for (int x = 0; x < ww; x++) {
      uint32 rr = 0, gg = 0, bb = 0, aa = 0;
      for (int yy = 0; yy < scale; yy++) {
        for (int xx = 0; xx < scale; xx++) {
          const auto [r, g, b, a] = GetPixel(x * scale + xx,
                                             y * scale + yy);

          // color contributions are alpha-weighted
          rr += r * a;
          gg += g * a;
          bb += b * a;
          aa += a;
        }
      }

      // Otherwise, the color can be anything, but output black.
      if (aa > 0) {
        rr /= aa;
        gg /= aa;
        bb /= aa;
        aa /= scale * scale;
      }
      ret.SetPixel(x, y, (uint8)rr, (uint8)gg, (uint8)bb, (uint8)aa);
    }
  }
  return ret;
}


void ImageRGBA::Clear32(uint32 color) {
  for (uint32 &c : rgba) c = color;
}

void ImageRGBA::Clear(uint8 r, uint8 g, uint8 b, uint8 a) {
  Clear32(Pack32(r, g, b, a));
}

// PERF: Make inline?
void ImageRGBA::BlendPixel(int x, int y,
                           uint8 r, uint8 g, uint8 b, uint8 a) {

  if ((unsigned)x >= (unsigned)width) return;
  if ((unsigned)y >= (unsigned)height) return;

  using word = uint16_t;

  const int i = (y * width + x);
  const auto &[old_r, old_g, old_b, old_a] = Unpack32(rgba[i]);
  // TODO: Figure out how to blend when dest is also transparent.

  // so a + oma = 255.
  const word oma = 0xFF - a;

  // we want (r * a/255) + (oldr * (1-a)/255),
  // which is (r * a)/255 + (oldr * (1-a))/255
  // which is (r * a + oldr * (1-a))/255
  //
  // (GCC computes the division with multiplication.)
  const word rr = (((word)r * (word)a) + (old_r * oma)) / 0xFF;
  const word gg = (((word)g * (word)a) + (old_g * oma)) / 0xFF;
  const word bb = (((word)b * (word)a) + (old_b * oma)) / 0xFF;
  // Note that the components cannot be > 0xFF.
  if (rr > 0xFF) __builtin_unreachable();
  if (gg > 0xFF) __builtin_unreachable();
  if (bb > 0xFF) __builtin_unreachable();

  rgba[i] = Pack32(rr, gg, bb, 0xFF);
}

void ImageRGBA::BlendPixel32(int x, int y, uint32 color) {
  const auto [r, g, b, a] = Unpack32(color);
  BlendPixel(x, y, r, g, b, a);
}

void ImageRGBA::BlendRect(int x, int y, int w, int h,
                          uint8 r, uint8 g, uint8 b, uint8 a) {
  // Easy to clip this to the screen. We could call an unclipped
  // version of BlendPixel here.
  if (y < 0) { w += y; y = 0; }
  if (x < 0) { h += x; x = 0; }

  const int yover = (y + h) - height;
  if (yover > 0) h -= yover;
  const int xover = (x + w) - width;
  if (xover > 0) w -= xover;

  if (w <= 0 || h <= 0) return;

  for (int yy = y; yy < y + h; yy++) {
    for (int xx = x; xx < x + w; xx++) {
      BlendPixel(xx, yy, r, g, b, a);
    }
  }
}

void ImageRGBA::BlendRect32(int x, int y, int w, int h, uint32 color) {
  const auto [r, g, b, a] = Unpack32(color);
  BlendRect(x, y, w, h, r, g, b, a);
}

void ImageRGBA::BlendBox32(int x, int y, int w, int h,
                           uint32 color, std::optional<uint32> cco) {
  const uint32 corner_color = cco.has_value() ? cco.value() : color;

  // Special cases
  if (w == 0 || h == 0) return;
  if (w == 1 && h == 1) {
    BlendPixel32(x, y, corner_color);
    return;
  } else if (w == 1) {
    BlendPixel32(x, y, corner_color);
    BlendPixel32(x, y + h - 1, corner_color);
    if (y + 1 < y + h - 1)
      BlendLine32(x, y + 1, x, y + h - 1, color);
  } else if (h == 1) {
    BlendPixel32(x, y, corner_color);
    BlendPixel32(x + w - 1, y, corner_color);
    if (x + 1 < x + w - 1)
      BlendLine32(x + 1, y, x + w - 1, y, color);
  }

  const int x1 = x + w - 1;
  const int y1 = y + h - 1;

  // PERF: straight lines can be faster by skipping bresenham

  if (x + 1 <= x1 - 1) {
    // Top
    BlendLine32(x + 1, y, x1 - 1, y, color);
    // Bottom
    BlendLine32(x + 1, y1, x1 - 1, y1, color);
  }
  if (y + 1 <= y1 - 1) {
    // Left
    BlendLine32(x, y + 1, x, y1 - 1, color);
    // Right
    BlendLine32(x1, y + 1, x1, y1 - 1, color);
  }

  BlendPixel32(x, y, corner_color);
  BlendPixel32(x1, y, corner_color);
  BlendPixel32(x, y1, corner_color);
  BlendPixel32(x1, y1, corner_color);
}

// PERF: For many of these where we call blendpixel in a loop,
// we could probably benefit by doing premultiplied alpha.

void ImageRGBA::BlendText32(int x, int y, uint32 color, const string &s) {
  auto SetPixel = [this, color](int xx, int yy) {
      this->BlendPixel32(xx, yy, color);
    };
  // SetPixel will clip, but exit early if we are totally off-screen.
  if (y >= height || y < -EmbeddedFont::HEIGHT) return;
  for (int i = 0; i < (int)s.size(); i++) {
    uint8 c = s[i];
    int xx = x + i * EmbeddedFont::WIDTH;
    if (xx >= width) return;
    EmbeddedFont::Blit(c, xx, y, SetPixel, [](int x, int y) {});
  }
}

void ImageRGBA::BlendText(int x, int y,
                          uint8 r, uint8 g, uint8 b, uint8 a,
                          const string &s) {
  BlendText32(x, y, Pack32(r, g, b, a), s);
}

void ImageRGBA::BlendText2x32(int x, int y, uint32 color, const string &s) {
  for (int i = 0; i < (int)s.size(); i++) {
    // Here we draw to "0,0", and then this function scales and translates.
    auto SetPixel = [x, y, i, this, color](int px, int py) {
        const int xx = x + (i * EmbeddedFont::WIDTH * 2) + px * 2;
        const int yy = y + py * 2;
        this->BlendPixel32(xx, yy, color);
        this->BlendPixel32(xx + 1, yy, color);
        this->BlendPixel32(xx, yy + 1, color);
        this->BlendPixel32(xx + 1, yy + 1, color);
      };

    const uint8 c = s[i];
    EmbeddedFont::Blit(c, 0, 0, SetPixel, [](int x, int y) {});
  }
}

void ImageRGBA::BlendText2x(int x, int y,
                            uint8 r, uint8 g, uint8 b, uint8 a,
                            const string &s) {
  BlendText2x32(x, y, Pack32(r, g, b, a), s);
}

void ImageRGBA::BlendLine32(int x1, int y1, int x2, int y2,
                            uint32 color) {
  const auto [r, g, b, a] = Unpack32(color);
  for (const auto [x, y] : Line<int>{x1, y1, x2, y2}) {
    BlendPixel(x, y, r, g, b, a);
  }
}

void ImageRGBA::BlendLine(int x1, int y1, int x2, int y2,
                          uint8 r, uint8 g, uint8 b, uint8 a) {
  for (const auto [x, y] : Line<int>{x1, y1, x2, y2}) {
    BlendPixel(x, y, r, g, b, a);
  }
}

void ImageRGBA::BlendLineAA32(float x1, float y1, float x2, float y2,
                              uint32 color) {
  uint8 r, g, b, a;
  std::tie(r, g, b, a) = Unpack32(color);
  auto Plot = [this, r, g, b, a](int x, int y, float f) {
      uint8 aa = f * a;
      this->BlendPixel(x, y, r, g, b, aa);
    };
  LineAA::Draw<int>(x1, y1, x2, y2, Plot);
}

void ImageRGBA::BlendLineAA(float x1, float y1, float x2, float y2,
                            uint8 r, uint8 g, uint8 b, uint8 a) {
  auto Plot = [this, r, g, b, a](int x, int y, float f) {
      uint8 aa = f * a;
      this->BlendPixel(x, y, r, g, b, aa);
    };
  LineAA::Draw<int>(x1, y1, x2, y2, Plot);
}

void ImageRGBA::BlendFilledCircle32(int x, int y, int r, uint32 color) {
  const int xmin = std::max(0, x - r);
  const int ymin = std::max(0, y - r);
  const int xmax = std::min(width - 1, x + r);
  const int ymax = std::min(height - 1, y + r);
  const int rr = r * r;
  for (int yy = ymin; yy <= ymax; yy++) {
    const int dy = yy - y;
    for (int xx = xmin; xx <= xmax; xx++) {
      const int dx = xx - x;
      const int dd = (dy * dy) + (dx * dx);
      if (dd <= rr) {
        BlendPixel32(xx, yy, color);
      }
    }
  }
}

// TODO: Does not handle overlap correctly.
void ImageRGBA::BlendImage(int x, int y, const ImageRGBA &other) {
  // PERF can factor out the pixel clipping here, supposing the
  // compiler cannot.
  for (int yy = 0; yy < other.height; yy++) {
    int yyy = y + yy;
    // Exit early if off-screen.
    if (yyy >= height) break;
    for (int xx = 0; xx < other.width; xx++) {
      int xxx = x + xx;
      if (xxx >= width) break;
      BlendPixel32(xxx, yyy, other.GetPixel32(xx, yy));
    }
  }
}

// TODO: Does not handle overlap correctly.
void ImageRGBA::BlendImageRect(int dstx, int dsty, const ImageRGBA &other,
                               int srcx, int srcy, int srcw, int srch) {
  for (int yy = 0; yy < srch; yy++) {
    const int syy = srcy + yy;
    const int dyy = dsty + yy;
    // Exit early if outside dstination.
    if (dyy >= height) break;
    // Exit early if outside source.
    if (syy >= other.height) break;

    if (syy >= 0 && dyy >= 0) {
      for (int xx = 0; xx < srcw; xx++) {
        const int sxx = srcx + xx;
        const int dxx = dstx + xx;
        if (dxx >= width) break;
        if (sxx >= other.width) break;

        if (sxx >= 0 && dxx >= 0) {
          BlendPixel32(dxx, dyy, other.GetPixel32(sxx, syy));
        }
      }
    }
  }
}

// This version assumes every pixel read/write is in bounds.
// Assumes the two images do not alias.
static void CopyImageRectUnclipped(
    uint32_t *dest, int destwidth,
    int dstx, int dsty,
    const uint32_t *source, int sourcewidth,
    int srcx, int srcy, int srcw, int srch) {
  for (int yy = 0; yy < srch; yy++) {
    const int syy = srcy + yy;
    const int dyy = dsty + yy;

    memcpy(dest + (dyy * destwidth) + dstx,
           source + (syy * sourcewidth) + srcx,
           srcw * sizeof (uint32_t));
  }
}

void ImageRGBA::CopyImage(int x, int y, const ImageRGBA &other) {
  // Compute clipped rectangle.
  int srcx = 0, srcy = 0;
  int srcw = other.Width(), srch = other.Height();
  // Clip left and top edges.
  if (x < 0) { srcx -= x; srcw += x; x = 0; }
  if (y < 0) { srcy -= y; srch += y; y = 0; }

  // Clip right and bottom.
  int right = (x + srcw) - width;
  if (right > 0) srcw -= right;
  int bottom = (y + srch) - height;
  if (bottom > 0) srch -= bottom;

  // Completely clipped; no-op.
  if (srcw <= 0 || srch <= 0) return;

  if (this == &other) {
    // Overlapping copies don't work. (PERF: If the src/dest
    // rectangles don't overlap, we can still use memcpy.)
    // If they come from the same image, always use slow path
    // via temporary:
    ImageRGBA tmp(srcw, srch);
    CopyImageRectUnclipped(
        tmp.rgba.data(), tmp.Width(),
        0, 0,
        other.rgba.data(), other.Width(),
        srcx, srcy, srcw, srch);
    CopyImageRectUnclipped(
        rgba.data(), width,
        x, y,
        tmp.rgba.data(), tmp.Width(),
        0, 0, srcw, srch);
  } else {
    CopyImageRectUnclipped(
        rgba.data(), width,
        x, y,
        other.rgba.data(), other.width,
        srcx, srcy, srcw, srch);
  }
}

// PERF merge with CopyImage (that approach is like 7x faster and
// handles overlap)
void ImageRGBA::CopyImageRect(int dstx, int dsty, const ImageRGBA &other,
                              int srcx, int srcy, int srcw, int srch) {
  for (int yy = 0; yy < srch; yy++) {
    const int syy = srcy + yy;
    const int dyy = dsty + yy;
    // Exit early if outside destination.
    if (dyy >= height) break;
    // Exit early if outside source.
    if (syy >= other.height) break;

    if (syy >= 0 && dyy >= 0) {
      for (int xx = 0; xx < srcw; xx++) {
        const int sxx = srcx + xx;
        const int dxx = dstx + xx;
        if (dxx >= width) break;
        if (sxx >= other.width) break;

        if (sxx >= 0 && dxx >= 0) {
          SetPixel32(dxx, dyy, other.GetPixel32(sxx, syy));
        }
      }
    }
  }
}


template<class F>
inline static ImageA Extract(const ImageRGBA &img, const F &f) {
  ImageA ret(img.Width(), img.Height());
  for (int y = 0; y < img.Height(); y++) {
    for (int x = 0; x < img.Width(); x++) {
      const auto [r, g, b, a] = img.GetPixel(x, y);
      ret.SetPixel(x, y, f(r, g, b, a));
    }
  }
  return ret;
}

ImageA ImageRGBA::Red() const {
  return Extract(*this, [](uint8 r, uint8 g, uint8 b, uint8 a) { return r; });
}

ImageA ImageRGBA::Green() const {
  return Extract(*this, [](uint8 r, uint8 g, uint8 b, uint8 a) { return g; });
}

ImageA ImageRGBA::Blue() const {
  return Extract(*this, [](uint8 r, uint8 g, uint8 b, uint8 a) { return b; });
}

ImageA ImageRGBA::Alpha() const {
  return Extract(*this, [](uint8 r, uint8 g, uint8 b, uint8 a) { return a; });
}

ImageRGBA ImageRGBA::FromChannels(const ImageA &red,
                                  const ImageA &green,
                                  const ImageA &blue,
                                  const ImageA &alpha) {
  const int width = red.Width();
  const int height = red.Height();
  CHECK(green.Width() == width);
  CHECK(green.Height() == height);
  CHECK(blue.Width() == width);
  CHECK(blue.Height() == height);
  CHECK(alpha.Width() == width);
  CHECK(alpha.Height() == height);
  ImageRGBA out(width, height);
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      const uint8 r = red.GetPixel(x, y);
      const uint8 g = green.GetPixel(x, y);
      const uint8 b = blue.GetPixel(x, y);
      const uint8 a = alpha.GetPixel(x, y);
      out.SetPixel(x, y, r, g, b, a);
    }
  }
  return out;
}

std::tuple<float, float, float, float>
ImageRGBA::SampleBilinear(float x, float y) const {
  // Truncate to integer pixels.
  int ix = x;
  int iy = y;

  // subpixel values give us the interpolants
  float fx = x - ix;
  float fy = y - iy;

  // Get these four values.
  //
  //  v00 ----- v10
  //   |   :fy   |
  //   |...*     | 1.0
  //   | fx      |
  //  v01 ----- v11
  //       1.0

  auto BClipPixel = [this](int x, int y) {
      if (x < 0) x = 0;
      if (y < 0) y = 0;
      if (x >= width) x = width - 1;
      if (y >= height) y = height - 1;
      return GetPixel(x, y);
    };

  using rgba8 = std::tuple<uint8, uint8, uint8, uint8>;

  rgba8 v00 = BClipPixel(ix, iy);
  rgba8 v10 = BClipPixel(ix + 1, iy);
  rgba8 v01 = BClipPixel(ix, iy + 1);
  rgba8 v11 = BClipPixel(ix + 1, iy + 1);

  auto Component = [fx, fy](
      uint8 c00, uint8 c10, uint8 c01, uint8 c11) -> float {
      // c0 interpolates between c00 and c10 at fx.
      float c0 = (float)c00 + (float)(c10 - c00) * fx;
      float c1 = (float)c01 + (float)(c11 - c01) * fx;

      float c = c0 + (c1 - c0) * fy;
      return std::clamp(c, 0.0f, 255.0f);
    };

  // TODO: Don't just average alpha; the average of #FF0000FF and
  // #00000000 should not be dark red.
  return std::make_tuple(
      // R
      Component(std::get<0>(v00), std::get<0>(v10),
                std::get<0>(v01), std::get<0>(v11)),
      // G
      Component(std::get<1>(v00), std::get<1>(v10),
                std::get<1>(v01), std::get<1>(v11)),
      // B
      Component(std::get<2>(v00), std::get<2>(v10),
                std::get<2>(v01), std::get<2>(v11)),
      // A
      Component(std::get<3>(v00), std::get<3>(v10),
                std::get<3>(v01), std::get<3>(v11)));
}


ImageA::ImageA(const vector<uint8> &alpha, int width, int height)
    : width(width), height(height), alpha(alpha) {
  CHECK((int)alpha.size() == width * height);
}

ImageA::ImageA(int width, int height) : width(width), height(height),
                                        alpha(width * height, 0) {
}

ImageA *ImageA::Copy() const {
  return new ImageA(alpha, width, height);
}

bool ImageA::operator==(const ImageA &other) const {
  return other.Width() == Width() &&
    other.Height() == Height() &&
    other.alpha == alpha;
}

std::size_t ImageA::Hash() const {
  uint64_t h = width * 31337;
  for (uint8 v : alpha) {
    // PERF: Work on 64 bits at a time...
    h = (h << 13) | (h >> (64 - 13));
    h ^= v;
    h *= 65537;
  }
  return (std::size_t)h;
}

ImageA ImageA::ScaleBy(int scale) const {
  // 1 is not useful, but it does work
  CHECK(scale >= 1);
  ImageA ret(width * scale, height * scale);
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      const uint8 color = GetPixel(x, y);
      for (int yy = 0; yy < scale; yy++) {
        for (int xx = 0; xx < scale; xx++) {
          ret.SetPixel(x * scale + xx,
                       y * scale + yy,
                       color);
        }
      }
    }
  }
  return ret;
}

void ImageA::Clear(uint8 value) {
  for (int i = 0; i < (int)alpha.size(); i++) alpha[i] = value;
}

float ImageA::SampleBilinear(float x, float y) const {
  // Truncate to integer pixels.
  int ix = x;
  int iy = y;

  // subpixel values give us the interpolants
  float fx = x - ix;
  float fy = y - iy;

  // Get these four values.
  //
  //  v00 ----- v10
  //   |   :fy   |
  //   |...*     | 1.0
  //   | fx      |
  //  v01 ----- v11
  //       1.0

  auto BClipPixel = [this](int x, int y) {
      if (x < 0) x = 0;
      if (y < 0) y = 0;
      if (x >= width) x = width - 1;
      if (y >= height) y = height - 1;
      return GetPixel(x, y);
    };

  uint8 v00 = BClipPixel(ix, iy);
  uint8 v10 = BClipPixel(ix + 1, iy);
  uint8 v01 = BClipPixel(ix, iy + 1);
  uint8 v11 = BClipPixel(ix + 1, iy + 1);

  // v0 interpolates between v00 and v10 at fx.
  float v0 = (float)v00 + (float)(v10 - v00) * fx;
  float v1 = (float)v01 + (float)(v11 - v01) * fx;

  float v = v0 + (v1 - v0) * fy;
  return std::clamp(v, 0.0f, 255.0f);
}

ImageA ImageA::ResizeBilinear(int nwidth, int nheight) const {
  ImageA ret{nwidth, nheight};
  // XXX Sampling is probably a little off wrt width-1 stuff?
  const float scale_width = width / (float)nwidth;
  const float scale_height = height / (float)nheight;
  for (int y = 0; y < nheight; y++) {
    const float sy = y * scale_height;
    for (int x = 0; x < nwidth; x++) {
      const float sx = x * scale_width;

      float fv = SampleBilinear(sx, sy);
      ret.SetPixel(x, y, std::roundf(fv));
    }
  }
  return ret;
}

ImageA ImageA::ResizeNearest(int nwidth, int nheight) const {
  ImageA ret{nwidth, nheight};
  // XXX Sampling is probably a little off wrt width-1 stuff?
  for (int y = 0; y < nheight; y++) {
    const float fy = y / (float)(nheight - 1);
    const int sy = std::roundf(fy * (height - 1));
    for (int x = 0; x < nwidth; x++) {
      const float fx = x / (float)(nwidth - 1);
      const int sx = std::roundf(fx * (width - 1));

      ret.SetPixel(x, y, GetPixel(sx, sy));
    }
  }
  return ret;
}

void ImageA::BlendText(int x, int y, uint8 v, const string &s) {
  auto SetPixel = [this, v](int xx, int yy) {
      this->BlendPixel(xx, yy, v);
    };
  // SetPixel will clip, but exit early if we are totally off-screen.
  if (y >= height || y < -EmbeddedFont::HEIGHT) return;
  for (int i = 0; i < (int)s.size(); i++) {
    uint8 c = s[i];
    int xx = x + i * EmbeddedFont::WIDTH;
    if (xx >= width) return;
    EmbeddedFont::Blit(c, xx, y, SetPixel, [](int x, int y) {});
  }
}

void ImageA::BlendImage(int x, int y, const ImageA &other) {
  // PERF can factor out the pixel clipping here, supposing the
  // compiler cannot.
  for (int yy = 0; yy < other.Height(); yy++) {
    int yyy = y + yy;
    // Exit early if we're off-screen.
    if (yyy >= height) break;
    for (int xx = 0; xx < other.Width(); xx++) {
      int xxx = x + xx;
      if (xxx >= width) break;
      BlendPixel(xxx, yyy, other.GetPixel(xx, yy));
    }
  }
}


ImageRGBA ImageA::GreyscaleRGBA() const {
  ImageRGBA rgba(width, height);
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      const uint8 v = GetPixel(x, y);
      // PERF if compiler can't optimize out the bounds checks here,
      // we can do it
      rgba.SetPixel(x, y, v, v, v, 0xFF);
    }
  }
  return rgba;
}

ImageRGBA ImageA::AlphaMaskRGBA(uint8 r, uint8 g, uint8 b) const {
  ImageRGBA rgba(width, height);
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      const uint8 v = GetPixel(x, y);
      // PERF if compiler can't optimize out the bounds checks here,
      // we can do it
      rgba.SetPixel(x, y, r, g, b, v);
    }
  }
  return rgba;
}


// --- ImageF ---

ImageF::ImageF(const vector<float> &alpha, int width, int height)
    : width(width), height(height), alpha(alpha) {
  CHECK((int)alpha.size() == width * height);
}

ImageF::ImageF(int width, int height) : width(width), height(height),
                                        alpha(width * height, 0.0f) {
}

ImageF::ImageF(const ImageA &other) : width(other.Width()),
                                      height(other.Height()) {
  alpha.resize(width * height);
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      alpha[y * width + x] =
        other.GetPixel(x, y) * (1.0f / 255.0f);
    }
  }
}

void ImageF::Clear(float value) {
  for (int i = 0; i < (int)alpha.size(); i++) alpha[i] = value;
}


// Shared by multiple versions of SampleBilinear. clip_pixel
// takes integer coordinates and returns the float value
// at that coordinate, handling clipping.
template<class ClipPixel>
static inline float SampleBilinearInternal(const ClipPixel &clip_pixel,
                                           float x, float y) {
  // Truncate to integer pixels.
  int ix = x;
  int iy = y;

  // subpixel values give us the interpolants
  float fx = x - ix;
  float fy = y - iy;

  // Get these four values.
  //
  //  v00 ----- v10
  //   |   :fy   |
  //   |...*     | 1.0
  //   | fx      |
  //  v01 ----- v11
  //       1.0

  float v00 = clip_pixel(ix, iy);
  float v10 = clip_pixel(ix + 1, iy);
  float v01 = clip_pixel(ix, iy + 1);
  float v11 = clip_pixel(ix + 1, iy + 1);

  // v0 interpolates between v00 and v10 at fx.
  float v0 = v00 + (v10 - v00) * fx;
  float v1 = v01 + (v11 - v01) * fx;

  float v = v0 + (v1 - v0) * fy;
  return v;
}

float ImageF::SampleBilinear(float x, float y) const {
  return SampleBilinearInternal(
      [this](int x, int y) -> float {
        if (x < 0) x = 0;
        if (y < 0) y = 0;
        if (x >= width) x = width - 1;
        if (y >= height) y = height - 1;
        return GetPixel(x, y);
      }, x, y);
}

float ImageF::SampleBilinear(float x, float y, float outside_value) const {
  return SampleBilinearInternal(
      [this, outside_value](int x, int y) -> float {
        if (x < 0 || y < 0 ||
            x >= width || y >= height)
          return outside_value;
        return GetPixel(x, y);
      }, x, y);
}

ImageF ImageF::ResizeBilinear(int nwidth, int nheight) const {
  ImageF ret{nwidth, nheight};
  // XXX Sampling is probably a little off wrt width-1 stuff?
  const float scale_width = width / (float)nwidth;
  const float scale_height = height / (float)nheight;
  for (int y = 0; y < nheight; y++) {
    const float sy = y * scale_height;
    for (int x = 0; x < nwidth; x++) {
      const float sx = x * scale_width;

      float fv = SampleBilinear(sx, sy);
      ret.SetPixel(x, y, fv);
    }
  }
  return ret;
}

void ImageF::BlendText(int x, int y, float v, const string &s) {
  auto SetPixel = [this, v](int xx, int yy) {
      this->BlendPixel(xx, yy, v);
    };
  // SetPixel will clip, but exit early if we are totally off-screen.
  if (y >= height || y < -EmbeddedFont::HEIGHT) return;
  for (int i = 0; i < (int)s.size(); i++) {
    uint8 c = s[i];
    int xx = x + i * EmbeddedFont::WIDTH;
    if (xx >= width) return;
    EmbeddedFont::Blit(c, xx, y, SetPixel, [](int x, int y) {});
  }
}


ImageA ImageF::Make8Bit() const {
  ImageA out(width, height);
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      const uint8 v = std::roundf(GetPixel(x, y) * 255.0f);
      out.SetPixel(x, y, v);
    }
  }
  return out;
}
