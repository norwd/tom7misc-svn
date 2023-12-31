
#include "pactom.h"

#include <string>
#include <vector>
#include <memory>
#include <cstdint>

#include "base/logging.h"
#include "base/stringprintf.h"
#include "geom/latlon.h"
#include "bounds.h"
#include "image.h"
#include "lines.h"
#include "arcfour.h"
#include "randutil.h"
#include "color-util.h"
#include "threadutil.h"
#include "pactom-util.h"

using namespace std;


using int64 = int64_t;

static constexpr int WIDTH = 1920;
static constexpr int HEIGHT = 1080;
static constexpr int SCALE = 4;
// Additional pixels to draw for line (0 = 1 pixel thick)
static constexpr int RADIUS = 4;
// circle at end while in motion
static constexpr int DOT_RADIUS = 18;

static constexpr int NUM_FRAMES = 512 + 256;

int main(int argc, char **argv) {
  ArcFour rc("pactom");
  unique_ptr<PacTom> pactom = PacTomUtil::Load(false);
  CHECK(pactom.get() != nullptr);

  std::vector<uint32_t> colors;
  for (const auto &r : pactom->runs) {
    uint32_t color =
      PacTomUtil::RandomBrightColor(&rc) & 0xFFFFFF77; // XXX
    colors.emplace_back(color);
  }

  const LatLon home = LatLon::FromDegs(40.452911, -79.936313);
  LatLon::Projection Project = LatLon::Gnomonic(home);
  const std::pair<double, double> home_pt = Project(home);

  // Find the extrema.
  Bounds bounds;
  for (const auto &r : pactom->runs) {
    for (const auto &[latlon, elev] : r.path) {
      auto [x, y] = Project(latlon);
      bounds.Bound(x, y);
    }
  }
  bounds.AddMarginFrac(0.05);

  // The bounds for the full data are the endpoint.
  Bounds::Scaler scaler_end = bounds.ScaleToFit(WIDTH * SCALE,
                                                HEIGHT * SCALE).FlipY();

  // Screen location of home.
  auto home_screen = scaler_end.Scale(home_pt);

  // Start scale; centered on home.
  Bounds::Scaler scaler_start = scaler_end.
    // put home at 0,0
    PanScreen(-home_screen.first, -home_screen.second).
    // zoom
    Zoom(6.0, 6.0).
    // center screen on 0, 0
    PanScreen(WIDTH * SCALE * 0.5, HEIGHT * SCALE * 0.5);

  auto ScaleInterp = [&](double f, std::pair<double, double> pt) ->
    std::pair<double, double> {
      auto [x0, y0] = scaler_start.Scale(pt);
      auto [x1, y1] = scaler_end.Scale(pt);

      return std::make_pair(f * x1 + (1.0 - f) * x0,
                            f * y1 + (1.0 - f) * y0);
    };

  auto MakeFrame = [&](int64_t frame) {
      const double frame_frac = frame / (double)(NUM_FRAMES - 1);

      ImageRGBA image(WIDTH * SCALE, HEIGHT * SCALE);
      image.Clear32(0x000000FF);

      for (const auto &[name, path] : pactom->hoods) {
        constexpr uint32 color = 0x909090FF;
        for (int i = 0; i < path.size() - 1; i++) {
          const LatLon latlon0 = path[i];
          const LatLon latlon1 = path[i + 1];
          auto [x0, y0] = ScaleInterp(frame_frac, Project(latlon0));
          auto [x1, y1] = ScaleInterp(frame_frac, Project(latlon1));

          PacTomUtil::DrawThickLine<RADIUS>(&image, x0, y0, x1, y1, color);
        }
      }

      for (int idx = 0; idx < pactom->runs.size(); idx++) {
        const uint32_t color = colors[idx];
        const auto &p = pactom->runs[idx].path;
        const int last_pt =
          std::clamp((int)std::round(p.size() * frame_frac), 0, (int)p.size());
        for (int i = 0; i < last_pt - 1; i++) {
          const auto &[latlon0, elev0] = p[i];
          const auto &[latlon1, elev1] = p[i + 1];
          auto [x0, y0] = ScaleInterp(frame_frac, Project(latlon0));
          auto [x1, y1] = ScaleInterp(frame_frac, Project(latlon1));

          PacTomUtil::DrawThickLine<RADIUS>(&image, x0, y0, x1, y1, color);

          if (i == last_pt - 2 && last_pt != p.size()) {
            uint32_t dot_color = color | 0x77;
            image.BlendFilledCircle32(x1, y1, DOT_RADIUS, dot_color);
          }
        }
      }

      ImageRGBA out = image.ScaleDownBy(SCALE);

      out.Save(StringPrintf("zoomout/zoom%04d.png", (int)frame));
      printf(".");
    };

  ParallelComp(NUM_FRAMES, MakeFrame, 8);

  return 0;
}
