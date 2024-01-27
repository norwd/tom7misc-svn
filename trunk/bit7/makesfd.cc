// Generate a text SFD file from a PNG containing a proportional font.
// The characters are arranged in a grid a la makegrid.cc. White
// (#fff) pixels give the character shapes, with a solid vertical
// black (#000) line gives (one past) the width.

// TODO: I changed this encoding and I need to update the DFX fonts!

#include <cstdint>
#include <string>
#include <vector>
#include <set>
#include <map>

#include "util.h"
#include "image.h"
#include "base/logging.h"
#include "fonts/island-finder.h"
#include "fonts/ttf.h"
#include "font-image.h"

using namespace std;
using uint8 = uint8_t;
using uint32 = uint32_t;
using uint64 = uint64_t;

using Glyph = FontImage::Glyph;

// TODO: Would be nice for this to be configurable!
// Standard layout of input file, based on characters
// that were useful for DestroyFX. If -1, the spot
// is unclaimed. Should be fine to extend past this many
// characters by increasing CHARS_DOWN too.
constexpr int MAPPED_CHARS_ACROSS = 16;
constexpr int MAPPED_CHARS_DOWN = 24;

static constexpr array<int, MAPPED_CHARS_ACROSS * MAPPED_CHARS_DOWN>
CODEPOINTS = {
  // First line
  // BLACK HEART SUIT
  0x2665,
  // BEAMED EIGHTH NOTES
  0x266B,
  // INFINITY
  0x221E,
  // SQUARE ROOT
  0x221A,
  // LESS THAN OR EQUAL TO
  0x2264,
  // GREATER THAN OR EQUAL TO
  0x2265,
  // APPROXIMATELY EQUAL
  0x2248,
  // EURO SIGN
  0x20AC,
  // ARROWS: LEFT, UP, RIGHT, DOWN
  0x2190, 0x2191, 0x2192, 0x2193,

  // EN DASH, EM DASH
  0x2013, 0x2014,

  // LEFT SINGLE QUOTE, RIGHT SINGLE QUOTE
  0x2018, 0x2019,
  // Second line

  // LEFT DOUBLE QUOTE, RIGHT DOUBLE QUOTE
  0x201C, 0x201D,

  // BULLET
  0x2022,
  // HORIZONTAL ELLIPSIS
  0x2026,
  // EMOJI: CLOUD
  0x2601,
  // EMOJI: ROCKET
  0x1F680,
  // EMOJI: NO ENTRY
  0x26D4,

  // dagger, double-dagger
  0x2020, 0x2021,

  // checkmark, heavy checkmark,
  0x2713, 0x2714,
  // ballot x, heavy ballot x,
  0x2717, 0x2718,

  // Trade Mark Sign
  0x2122,

  // Ideographic full stop (big japanese period)
  0x3002,
  // turnstile (a.k.a. right tack)
  0x22A2,

  // space for emoji
  // EMOJI: LIGHT BULB
  0x1F4A1,
  // EMOJI: BEER MUG
  0x1F37A,
  // EMOJI: WASTEBASKET
  0x1F5D1,
  // EMOJI: MOAI HEAD
  0x1F5FF,
  // EMOJI: HIGH VOLTAGE
  0x26A1,
  // EMOJI: MAGNET
  0x1F9F2,
  // EMOJI: SKULL
  0x1F480,
  // EMOJI: SKULL AND CROSSBONES
  0x2620,
  // EMOJI: DROPLET
  0x1F4A7,
  // EMOJI: HUNDRED POINTS
  0x1F4AF,
  // EMOJI: ANGER SYMBOL
  0x1F4A2,
  // EMOJI: ZZZ
  0x1F4A4,
  // EMOJI: PAGE FACING UP
  0x1F4C4,
  // EMOJI: BOMB
  0x1F4A3,
  // EMOJI: GLOBE WITH MERIDIANS
  0x1F310,
  // EMOJI: EYES
  0x1F440,

  // Emoji line 2.

  // EMOJI: TOOTHBRUSH
  0x1FAA5,
  // EMOJI: HEADSTONE
  0x1FAA6,
  // EMOJI: PLACARD (Signpost)
  0x1FAA7,
  // EMOJI: ROCK
  0x1FAA8,
  // EMJOI: FLY
  0x1FAB0,

  // EMOJI: MAGIC WAND
  0x1FA84,
  // EMOJI: COIN
  0x1FA99,
  // EMOJI: LADDER
  0x1FA9C,

  // EMOJI: HOT PEPPER
  0x1F336,

  // EMOJI: GHOST
  0x1F47B,

  // EMOJI: KEY
  0x1F511,

  // EMOJI: LOCK (LOCKED)
  0x1F512,
  // EMOJI: OPEN LOCK
  0x1F513,

  // EMOJI: HEAVY DOLLAR SIGN
  0x1F4B2,

  -1, -1,


  // ASCII, in order
  0x20, 0x21, 0x22, 0x23, 0x24, 0x25, 0x26, 0x27, 0x28, 0x29, 0x2A, 0x2B, 0x2C, 0x2D, 0x2E, 0x2F,
  0x30, 0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37, 0x38, 0x39, 0x3A, 0x3B, 0x3C, 0x3D, 0x3E, 0x3F,
  0x40, 0x41, 0x42, 0x43, 0x44, 0x45, 0x46, 0x47, 0x48, 0x49, 0x4A, 0x4B, 0x4C, 0x4D, 0x4E, 0x4F,
  0x50, 0x51, 0x52, 0x53, 0x54, 0x55, 0x56, 0x57, 0x58, 0x59, 0x5A, 0x5B, 0x5C, 0x5D, 0x5E, 0x5F,
  0x60, 0x61, 0x62, 0x63, 0x64, 0x65, 0x66, 0x67, 0x68, 0x69, 0x6A, 0x6B, 0x6C, 0x6D, 0x6E, 0x6F,
  0x70, 0x71, 0x72, 0x73, 0x74, 0x75, 0x76, 0x77, 0x78, 0x79, 0x7A, 0x7B, 0x7C, 0x7D, 0x7E, -1,

  // white king, queen, rook, bishop, knight, pawn
  0x2654, 0x2655, 0x2656, 0x2657, 0x2658, 0x2659,
  // black
  0x265A, 0x265B, 0x265C, 0x265D, 0x265E, 0x265F,

  // Three free before replacement char
  -1, -1, -1,
  // <?> replacement char
  0xFFFD,

  // Black circle, black square
  0x25CF, 0x25A0,
  // geometric shapes line, unclaimed
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,

  // Unicode Latin-1 Supplement, mapped to itself.
  // See https://en.wikibooks.org/wiki/Unicode/Character_reference/0000-0FFF
  0xA0, 0xA1, 0xA2, 0xA3, 0xA4, 0xA5, 0xA6, 0xA7, 0xA8, 0xA9, 0xAA, 0xAB, 0xAC, 0xAD, 0xAE, 0xAF,
  0xB0, 0xB1, 0xB2, 0xB3, 0xB4, 0xB5, 0xB6, 0xB7, 0xB8, 0xB9, 0xBA, 0xBB, 0xBC, 0xBD, 0xBE, 0xBF,
  0xC0, 0xC1, 0xC2, 0xC3, 0xC4, 0xC5, 0xC6, 0xC7, 0xC8, 0xC9, 0xCA, 0xCB, 0xCC, 0xCD, 0xCE, 0xCF,
  0xD0, 0xD1, 0xD2, 0xD3, 0xD4, 0xD5, 0xD6, 0xD7, 0xD8, 0xD9, 0xDA, 0xDB, 0xDC, 0xDD, 0xDE, 0xDF,
  0xE0, 0xE1, 0xE2, 0xE3, 0xE4, 0xE5, 0xE6, 0xE7, 0xE8, 0xE9, 0xEA, 0xEB, 0xEC, 0xED, 0xEE, 0xEF,
  0xF0, 0xF1, 0xF2, 0xF3, 0xF4, 0xF5, 0xF6, 0xF7, 0xF8, 0xF9, 0xFA, 0xFB, 0xFC, 0xFD, 0xFE, 0xFF,

  // unclaimed
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,

  // Greek. We skip the characters that look the same in the Latin
  // alphabet: A B E Z H I K M N O P T Y X v u x.

  // Gamma, Delta, Theta, Lambda, Xi
  0x0393, 0x0394, 0x0398, 0x039B, 0x039E,
  // Pi, Sigma, Phi, Psi, Omega,
  0x03A0, 0x03A3, 0x03A6, 0x03A8, 0x03A9,
  // alpha, beta, gamma, delta, epsilon, zeta
  0x03B1, 0x03B2, 0x03B3, 0x03B4, 0x03B5, 0x03B6,

  // Line 2:
  // eta, theta, iota, kappa, lambda, mu, (no nu), xi, (no omicron)
  0x03B7, 0x03B8, 0x03B9, 0x03BA, 0x03BB, 0x03BC, 0x03BE,
  // pi, rho, (final) sigma, sigma, tau, (no upsilon), phi, (no chi), omega
  0x03C0, 0x03C1, 0x03C2, 0x03C3, 0x03C4, 0x03C6, 0x03C8, 0x03C9,
  // one unclaimed spot at the end of greek
  -1,

  // math
  // exists, forall
  0x2203, 0x2200,
  // rest of math, unclaimed
  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,

  // Block Elements, in unicode order
  0x2580, 0x2581, 0x2582, 0x2583, 0x2584, 0x2585, 0x2586, 0x2587,
  0x2588, 0x2589, 0x258A, 0x258B, 0x258C, 0x258D, 0x258E, 0x258F,
  0x2590, 0x2591, 0x2592, 0x2593, 0x2594, 0x2595, 0x2596, 0x2597,
  0x2598, 0x2599, 0x259A, 0x259B, 0x259C, 0x259D, 0x259E, 0x259F,
};

// e.g. use the glyph for hyphen (0x2D) to render U+2212 (minus).
static constexpr std::initializer_list<std::pair<int, int>>
REUSE_FOR = {
  // hyphen used as minus
  {0x2D, 0x2212},
  // ascii -> cyrillic
  {'S', 0x0405},
  {'J', 0x0408},
  {'A', 0x0410},
  {'B', 0x0412},
  {'E', 0x0415},
  {'M', 0x041C},
  {'H', 0x041D},
  {'O', 0x041E},
  {'P', 0x0420},
  {'C', 0x0421},
  {'T', 0x0422},
  {'X', 0x0425},
  {'a', 0x0430},
  {'e', 0x0435},
  {'o', 0x043E},
  {'p', 0x0440},
  {'c', 0x0441},
  {'x', 0x0445},
  {'s', 0x0455},
  {'i', 0x0456},
  {'j', 0x0458},
  // TODO: More cyrillic can be copied from Latin-1, Greek.

  // ascii -> greek
  {'J', 0x037F},
  {'A', 0x0391},
  {'B', 0x0392},
  {'E', 0x0395},
  {'Z', 0x0396},
  {'H', 0x0397},
  {'I', 0x0399},
  {'K', 0x039A},
  {'M', 0x039C},
  {'N', 0x039D},
  {'O', 0x039F},
  {'P', 0x03A1},
  {'T', 0x03A4},
  {'Y', 0x03A5},
  {'X', 0x03A7},
  {'v', 0x03BD},
  {'o', 0x03BF},
  {'u', 0x03C5},
  {'x', 0x03C7},

  // Full-width comma
  {',', 0xFF0C},
  // Ideographic space
  {' ', 0x3000},
  // Fullwidth parentheses
  // Actually we can just map these all from ascii?
  // en.wikipedia.org/wiki/Halfwidth_and_Fullwidth_Forms_(Unicode_block)
  {'(', 0xFF08},
  {')', 0xFF09},

  // bullet -> katakana middle dot
  {0x2022, 0x30FB},

  // Black circle -> black circle for record
  {0x25CF, 0x23FA},
  // Same for square
  {0x25A0, 0x23F9},

};


// XXX remove
static constexpr bool VERBOSE = false;

template<class C, class K>
static bool ContainsKey(const C &c, const K &k) {
  return c.find(k) != c.end();
}

static Config ParseAndCheckConfig(const std::string &cfgfile) {
  Config config = Config::ParseConfig(cfgfile);
  CHECK(!config.pngfile.empty()) << "Required config line: pngfile";
  CHECK(!config.name.empty()) << "Required config line: name";

  CHECK(config.charbox_width > 0) << "Config line charbox-width must be >0";
  CHECK(config.charbox_height > 0) << "Config line charbox-height must be >0";

  CHECK(config.descent >= 0) << "Config line charbox-height must be >= 0";

  CHECK(config.chars_across > 0);
  CHECK(config.chars_down > 0);

  return config;
}

// Given a series of points on the grid that trace a proper outline
// (e.g., it has nonzero area, consecutive points are in different
// locations), generate an equivalent but more efficient outline
// by skipping points that are colinear with their neighbors. (The
// routine below generates one on every pixel corner, even when
// unnecessary.) Note that this does not handle a case like
// 1 ----- 3 ----- 2
// where a line doubles back on itself.
static vector<pair<int, int>> RemoveColinearPoints(
    const vector<pair<int, int>> &points) {
  CHECK(points.size() >= 3) << "Degenerate contour; too small!";
  // Start on a corner so that we don't need to think about that
  // edge case.
  const int corner_idx = [&points](){
      for (int idx = 0; idx < (int)points.size(); idx++) {
        int prev_idx = idx == 0 ? points.size() - 1 : (idx - 1);
        int next_idx = idx == ((int)points.size() - 1) ? 0 : idx + 1;

        // If these three points are colinear, then they will
        // share an x coordinate or y coordinate. Otherwise,
        // the center one is a corner.
        const auto [px, py] = points[prev_idx];
        const auto [x, y] = points[idx];
        const auto [nx, ny] = points[next_idx];
        if ((px == x && x == nx) ||
            (py == y && y == ny)) {
          // colinear
        } else {
          return idx;
        }
      }
      LOG(FATAL) << "Degenerate contour; no area!";
    }();

  // We definitely keep the corner index. Now loop over all the
  // points starting there, and emit points if they
  vector<pair<int, int>> out;
  out.reserve(points.size());
  out.push_back(points[corner_idx]);
  auto Observe = [&points, &out](int idx) {
      CHECK(idx >= 0 && idx < (int)points.size());
      int next_idx = idx == ((int)points.size() - 1) ? 0 : idx + 1;
      auto [px, py] = out.back();
      auto [x, y] = points[idx];
      auto [nx, ny] = points[next_idx];
      if ((px == x && x == nx) ||
          (py == y && y == ny)) {
        // colinear. skip it.
      } else {
        out.emplace_back(x, y);
      }
    };
  for (int i = corner_idx + 1; i < (int)points.size(); i++)
    Observe(i);
  for (int i = 0; i < corner_idx; i++)
    Observe(i);
  return out;
}

// Scale these coordinates, probably?
static TTF::Contour MakeContour(const vector<pair<int, int>> &points) {
  // Just return straight lines between these edge points.
  TTF::Contour ret;
  for (const auto &[ex, ey] : points) {
    ret.paths.emplace_back((float)ex, (float)ey);
  }
  return ret;
}

// Trace a single pixel blob in a bitmap, producing a single clockwise
// contour. This returns a series of points, on each of the pixel
// corners, which should be further simplified to remove colinear
// points. (0, 0) is the top-left corner of the top-left pixel, with
// y increasing downward.
//
// The bitmap must have oine contiguous non-empty region with values >0,
// which is the shape to trace.
//
// This is based on code that was tracing SDFs (from ../lowercase) which
// might account for some overkill therein?
static vector<pair<int, int>> VectorizeOne(const ImageA &bitmap) {
  auto InBlob = [&bitmap](int x, int y) -> bool {
      if (x < 0 || y <0 || x >= bitmap.Width() || y >= bitmap.Height())
        return false;
      return bitmap.GetPixel(x, y) > 0;
    };

  // First, find a pixel inside the blob.
  // This pixel has the property that there is no pixel with
  // a smaller y coordinate, which is also in the blob.
  const auto [startpx, startpy] = [&bitmap, InBlob]() ->
    std::pair<int, int> {
      for (int y = 0; y < bitmap.Height(); y++) {
        for (int x = 0; x < bitmap.Width(); x++) {
          if (InBlob(x, y)) return make_pair(x, y);
        }
      }
      LOG(FATAL) << "VectorizeOne requires a non-empty bitmap!";
  }();

  // We wind around the pixel blob's exterior, always maintaining a
  // direction and a pair of pixels, one in, and one out. The starting
  // pixel we just found is such an example we scanned from top to
  // bottom. We'll be done when we return to the start pixel.
  CHECK(!InBlob(startpx, startpy - 1)) << "Need the uppermost pixel "
    "in this column.";

  // Discrete direction. The code below is written for a pattern
  // where we are moving right, with the blob down, and the
  // exterior up (which is the start condition), but it naturally
  // rotates to the other directions.
  enum Dir {
    UP,
    DOWN,
    LEFT,
    RIGHT,
  };

  // Get the orthogonal "normal" direction, which is up for Right.
  auto Normal = [](Dir d) {
      switch (d) {
      case RIGHT: return UP;
      case LEFT: return DOWN;
      case DOWN: return RIGHT;
      case UP: return LEFT;
      }
      CHECK(false) << "Bad dir";
    };

  auto TurnCCW = Normal;

  auto TurnCW = [](Dir d) {
      switch (d) {
      case RIGHT: return DOWN;
      case DOWN: return LEFT;
      case LEFT: return UP;
      case UP: return RIGHT;
      }
      CHECK(false) << "Bad dir";
    };

  auto Move = [](int x, int y, Dir d) -> pair<int, int> {
    switch (d) {
    case RIGHT: return make_pair(x + 1, y);
    case LEFT: return make_pair(x - 1, y);
    case DOWN: return make_pair(x, y + 1);
    case UP: return make_pair(x, y - 1);
    }
    CHECK(false) << "Bad dir";
  };

  // Return the starting corner when we are visiting the pixel at px,py
  // in the given direction. For example, if we are moving right, then
  // the top-left corner of the pixel (which is px,py) is the start of
  // a clockwise path around it. If we are moving up, then it is the
  // bottom-left corner (px, py+1), etc.
  auto SourceCorner = [](int px, int py, Dir dir) {
      switch (dir) {
      case RIGHT: return make_pair(px, py);
      case LEFT: return make_pair(px + 1, py + 1);
      case DOWN: return make_pair(px + 1, py);
      case UP: return make_pair(px, py + 1);
      }
      CHECK(false) << "Bad dir";
    };

  // Pixel we're currently looking at.
  int px = startpx;
  int py = startpy;

  // Direction we're currently heading. We are at the top of the
  // blob, so go right for clockwise. (It seems any local top
  // would work; compare for example the inner top edges of an 's'
  // shape.)
  Dir right = RIGHT;


  vector<std::pair<int, int>> edge_points;
  for (;;) {
    Dir up = Normal(right);
    const auto [upx, upy] = Move(px, py, up);
    // Invariant is that we are on the edge, so px,py is
    // in the blob and the pixel "above" it is not.
    CHECK(InBlob(px, py));
    CHECK(!InBlob(upx, upy));

    edge_points.push_back(SourceCorner(px, py, right));

    // We're in a situation like this (perhaps under some rotation),
    // traveling along the top edge of the filled pixel at px,py.
    // We'll proceed by case analysis on the two pixels ahead of us.
    //
    // source corner
    //    |  +--+
    //    v  |a?|
    //    +->+--+
    //    |##|b?|
    //    +--+--+

    const auto [ax, ay] = Move(upx, upy, right);
    const auto [bx, by] = Move(px, py, right);
    const bool a = InBlob(ax, ay);
    const bool b = InBlob(bx, by);

    if (!a && b) {
      // +--+--+
      // |  |a |
      // +--+--+
      // |##|b#|
      // +--+--+
      // Just continue in the same direction.
      px = bx;
      py = by;
    } else if (!a && !b) {
      // +--+--+
      // |  |a |
      // +--+--+
      // |##|b |
      // +--+--+
      // Make a 90 degree turn around this pixel, but
      // stay on it.
      right = TurnCW(right);
    } else {
      CHECK(a);
      // +--+--+
      // |  |a#|
      // +--+--+
      // |##|b?|
      // +--+--+
      // Don't care what b is (we are using 4-connectivity);
      // if it's open we'll get there separately.

      px = ax;
      py = ay;
      right = TurnCCW(right);
    }

    // Consider the case of a single pixel. We should
    // only end when we approach it with right = RIGHT, right?
    if (px == startpx &&
        py == startpy && right == RIGHT)  {
      // printf("Loop finished!\n");
      break;
    }
  }

  return edge_points;
}

// Returns coordinates measured in pixels, not in the normalized [0,1]
// representation. Caller should normalize.
static TTF::Char Vectorize(const Glyph &glyph) {
  const IslandFinder::Maps maps = IslandFinder::Find(glyph.pic);

  // Tracing follows the same recursive approach as in SDF-based code
  // (../lowercase/font-problem) but is much simpler here since we
  // have exact data.

  std::function<vector<TTF::Contour>(int, uint8)> VectorizeRec =
    [&maps, &VectorizeRec](int d, uint8 parent) -> vector<TTF::Contour> {
      if (VERBOSE) printf("DEPTH %d/%d\n", d, maps.max_depth);
      if (d > maps.max_depth) return {};

      std::vector<TTF::Contour> contours;

      // Get the equivalence classes at this depth, paired with the set
      // of (strict) descendants of that class (we want to remove
      // these holes when tracing the contour to simplify our lives).
      // Only consider classes that have 'parent' as an ancestor.
      // We'll run the routine below on each one.
      std::map<uint8, std::set<uint8>> eqclasses;
      for (int y = 0; y < maps.eqclass.Height(); y++) {
        for (int x = 0; x < maps.eqclass.Width(); x++) {
          if ((int)maps.depth.GetPixel(x, y) >= d) {
            uint8 eqc = maps.eqclass.GetPixel(x, y);
            if (maps.HasAncestor(eqc, parent)) {
              // This must exist because the depth is at least d, and
              // d > 0.
              uint8 ancestor = maps.GetAncestorAtDepth(d, eqc);
              eqclasses[ancestor].insert(eqc);
            }
          }
        }
      }

      // Now, for each component...
      for (const auto &[this_eqc, descendants] : eqclasses) {
        if (VERBOSE) {
          printf("Tracing eqc %d (descendants:", this_eqc);
          for (uint8 d : descendants) printf(" %d", d);
          printf(")\n");
        }
        // Generate a simplified bitmap.
        ImageA bitmap(maps.eqclass.Width(), maps.eqclass.Height());
        for (int y = 0; y < maps.eqclass.Height(); y++) {
          for (int x = 0; x < maps.eqclass.Width(); x++) {
            uint8 eqc = maps.eqclass.GetPixel(x, y);
            bool inside = eqc == this_eqc || ContainsKey(descendants, eqc);
            bitmap.SetPixel(x, y, inside ? 0xFF : 0x00);
          }
        }

        vector<pair<int, int>> points =
          RemoveColinearPoints(VectorizeOne(bitmap));
        contours.push_back(MakeContour(points));

        // Now recurse on descendants, if any.
        if (!descendants.empty()) {
          const auto child_contours = VectorizeRec(d + 1, this_eqc);
          // ... but reverse the winding order so that these cut out
          // (or maybe get reversed again).
          for (TTF::Contour cc : child_contours)
            contours.push_back(TTF::ReverseContour(cc));
        }
      }

      return contours;
    };

  // Start at depth 1, since we do not want any outline for the
  // surrounding "sea".
  vector<TTF::Contour> contours = VectorizeRec(1, 0);

  TTF::Char ttf_char;
  // TODO: Honor left_edge
  ttf_char.contours = std::move(contours);
  ttf_char.width = (float)glyph.pic.Width();
  return ttf_char;
}

int main(int argc, char **argv) {
  CHECK(argc == 3 || argc == 4) <<
    "Usage: ./makesfd.exe config.cfg out.sfd [testpattern.png]\n";

  const Config config = ParseAndCheckConfig(argv[1]);
  if (VERBOSE)
    printf("Converting from %s\n", argv[1]);
  const string out_sfd = argv[2];
  const string out_test_png = (argc > 3) ? argv[3] : "";

  FontImage font_image(config);
  // XXX
  auto &font = font_image.glyphs;

  if (config.no_lowercase) {
    for (int c = 'A'; c <= 'Z'; c++) {
      int lc = c | 32;
      bool lc_missing = font.find(lc) == font.end() ||
        (config.fixed_width && FontImage::EmptyGlyph(font[lc]));
      if (font.find(c) != font.end() && lc_missing) {
        font[lc] = font[c];
      }
    }
  }

  if (!out_test_png.empty()) {
    const int BORDER = 2;
    const int output_height = config.charbox_height + config.extra_linespacing;

    // Output test pattern PNG.
    // Heart can't actually go here because it is \0.
    #define HEART "<3"
    #define INFTY "\x02"
    #define NOTES "\x01"
    #define PLUSMINUS "\x03"
    #define DEG    "\x04"
    vector<string> testpattern = {
     "  Welcome to my font!  it is cozy here " HEART "  (ok) ",
     "  Now is the FALL-TIME of our DISCONTENT !!|1Il ",
     "",
     "",
     "  " NOTES " Enable hyper-drive      for (;;) {",
     "  " NOTES " Enable ultra-disc         printf(\"hi?\\n\"); ",
     "  " NOTES " Disable introspection   }",
     "",
     "  Mr. Jock, TV Quiz Ph.D., bags few lynx!  ",
     "  (glib jocks quiz nymph to vex dwarf) ",
     "  (SYMPATHIZING WOULD FIX QUAKER OBJECTIVES.) ",
     "  XW!@#$%^&*()-=_+{}[]\\|:\";'<>?,./ZXCVB~` ",
     "",
     "  " PLUSMINUS "123,456 * 7,890 = 974,067,840" DEG,
     "",
     "  jungle quip, " INFTY " If you knew where you'd fall,",
     "  TTTTTT QQQQ` " INFTY  " you'd put a pillow!",
     "  http://.com/ " INFTY  " (watch--said I--beloved)",
    };

    ImageRGBA test(config.charbox_width * 48 + 2 * BORDER,
                   output_height * testpattern.size() + 2 * BORDER);
    test.Clear32(0x000033FF);

    for (int lidx = 0; lidx < (int)testpattern.size(); lidx++) {
      const string &line = testpattern[lidx];
      const int ypos = BORDER + lidx * output_height;
      int xpos = BORDER;
      for (int cidx = 0; cidx < (int)line.size(); cidx++) {
        const int codepoint = line[cidx];
        auto it = font.find(codepoint);
        if (it != font.end()) {
          // Converting to imagergba each time obviously wasteful...!
          ImageRGBA glyph = it->second.pic.AlphaMaskRGBA(0xEE, 0xEE, 0xFF);
          test.BlendImage(xpos, ypos, glyph);
          xpos += glyph.Width();
        }
      }
    }

    test.ScaleBy(3).Save(out_test_png);
  }

  const double one_pixel = 1.0 / config.charbox_height;

  const int chars_across = config.chars_across;

  TTF::Font ttf_font;
  for (const auto &[index, glyph] : font) {
    // Don't warn about empty glyphs if the font is fixed-width, as
    // there is no other way to indicate a missing glyph.
    bool ok_missing = config.fixed_width && FontImage::EmptyGlyph(glyph);

    if (index >= (int)CODEPOINTS.size()) {
      if (!ok_missing) {
        printf("Skipping glyph at %d,%d because it is outside the codepoint "
               "array!\n", index % chars_across,
               index / chars_across);
        printf("%s", FontImage::GlyphString(glyph).c_str());
      }
      continue;
    }
    CHECK(index >= 0 && index < (int)CODEPOINTS.size());
    const int codepoint = CODEPOINTS[index];
    if (codepoint < 0) {
      if (!ok_missing) {
        printf("Skipping glyph at %d,%d because the codepoint is not "
               "configured!\n", index % chars_across, index / chars_across);
        printf("%s", FontImage::GlyphString(glyph).c_str());
      }
    } else {
      TTF::Char ch = Vectorize(glyph);

      ch.width *= one_pixel;
      TTF::MapCoords([one_pixel](float x, float y) {
          return make_pair(x * one_pixel, y * one_pixel);
        }, &ch);

      ttf_font.chars[codepoint] = std::move(ch);
    }
  }

  for (const auto &[src, dst] : REUSE_FOR) {
    // If we do have the source, but don't have the dest, copy.
    // (PERF: I think it's possible for a single outline to be used
    // for multiple characters, so we should do that instead! Or
    // dedupe as a separate matter.)
    if (ttf_font.chars.find(src) != ttf_font.chars.end() &&
        ttf_font.chars.find(dst) == ttf_font.chars.end()) {
      printf("Copy %04x to %04x\n", src, dst);
      ttf_font.chars[dst] = ttf_font.chars[src];
    }
  }

  ttf_font.baseline = 1.0 - (config.descent * one_pixel);
  ttf_font.linegap = config.extra_linespacing * one_pixel;
  // Might only affect FontForge, but it at least looks better in the
  // editor without anti-aliasing.
  ttf_font.antialias = false;
  ttf_font.bitmap_grid_height = config.charbox_height;
  // Reserved for Tom 7!
  ttf_font.vendor = {'F', 'r', 'o', 'g'};
  ttf_font.copyright = config.copyright;

  const string sfd = ttf_font.ToSFD(config.name);
  Util::WriteFile(out_sfd, sfd);

  return 0;
}
