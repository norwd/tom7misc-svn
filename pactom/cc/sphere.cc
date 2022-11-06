
#include "yocto_geometry.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <variant>
#include <utility>
#include <optional>
#include <functional>
#include <cmath>
#include <numbers>

#include "threadutil.h"
#include "image.h"
#include "color-util.h"
#include "base/logging.h"
#include "base/stringprintf.h"
#include "timer.h"
#include "randutil.h"
#include "arcfour.h"
#include "opt/opt.h"

static constexpr double PI = std::numbers::pi;

static constexpr int THREADS = 24;

static constexpr double EARTH_RADIUS = 2.8f;
// static constexpr double SUN_DISTANCE = 23214.0f * EARTH_RADIUS;
static constexpr double SUN_DISTANCE = 1000.0f;
static constexpr double SUN_RADIUS = 108.15f * EARTH_RADIUS;

static constexpr double FAR_DIST = 100.0f;
static constexpr double FAR_WIDTH = 100.0f;
// static constexpr double NEAR_WIDTH = 0.01f;
static constexpr double NEAR_WIDTH = 0.005f;

static constexpr float TOTAL_MILES = 800 + 1400 + 1800 + 25;
// Normalized to [0, 1].
static constexpr ColorUtil::Gradient INNER_EARTH {
  GradRGB(0.0f / TOTAL_MILES, 0xEEEECC),
  GradRGB(790.0f / TOTAL_MILES, 0xEEEECC),
  GradRGB(801.0f / TOTAL_MILES, 0xf7dd24),
  GradRGB(2200.0f / TOTAL_MILES, 0xf7dd24),
  GradRGB(2201.0f / TOTAL_MILES, 0xf7dd24),
  GradRGB(4000.0f / TOTAL_MILES, 0xcf1857),
  GradRGB(4001.0f / TOTAL_MILES, 0x724719),
  GradRGB(4024.0f / TOTAL_MILES, 0x724719),
  GradRGB(1.0f, 0x046374),
};

using namespace yocto;

struct Convert {
  // uv coordinates here are given in terms of the earth
  // texture map image (already flipped horizontally).

  // Wellsville
  static constexpr double xu0 = 11921.0 / (21600.0 * 2.0);
  static constexpr double yu0 = 5928.0 / 21600.0;

  static constexpr double lat0 = 40.600658; // 40.5815756;
  static constexpr double lon0 = -80.646501; // -80.5735614;

  // Confluence
  static constexpr double xu1 = 11999.0 / (21600.0 * 2.0);
  static constexpr double yu1 = 5947.0 / 21600.0;
  static constexpr double lat1 = 40.442487; // 40.4452293;
  static constexpr double lon1 = -80.015946; // -80.140803;

  static constexpr double lattoy = (yu1 - yu0) / (lat1 - lat0);
  static constexpr double lontox = (xu1 - xu0) / (lon1 - lon0);

  static constexpr double ytolat = (lat1 - lat0) / (yu1 - yu0);
  static constexpr double xtolon = (lon1 - lon0) / (xu1 - xu0);

  // Near Pittsburgh. Assumes world is flat, which is of
  // course not true. Remember lat,lon is y,x.
  static constexpr std::pair<double, double> ToUV(double lat, double lon) {
    return make_pair((lon - lon0) * lontox + xu0,
                     (lat - lat0) * lattoy + yu0);
  }

  static constexpr std::pair<double, double> ToLatLon(double u, double v) {
    return make_pair((u - xu0) * xtolon + lon0,
                     (v - yu0) * ytolat + lat0);
  }
};

static ImageRGBA *bluemarble = nullptr;

mat3f RotYaw(float a) {
  const float cosa = cos(a);
  const float sina = sin(a);
  return mat3f
    {cosa, -sina, 0.0f,
     sina, cosa,  0.0f,
     0.0f, 0.0f,  1.0f};
}

mat3f RotPitch(float a) {
  const float cosa = cos(a);
  const float sina = sin(a);

  return mat3f
    {cosa,  0.0f, sina,
     0.0f,  1.0f, 0.0f,
     -sina, 0.0f, cosa};
}

mat3f RotRoll(float a) {
  const float cosa = cos(a);
  const float sina = sin(a);

  return mat3f
    {1.0f, 0.0f, 0.0f,
     0.0f, cosa, -sina,
     0.0f, sina, cosa};
}

mat3f Rot(double yaw, double pitch, double roll) {
  mat3f mr = RotRoll(roll);
  mat3f mp = RotPitch(pitch);
  mat3f my = RotYaw(yaw);
  mat3f m = mp * my;
  mat3f n = mr * m;
  return n;
}

struct Tetrahedron {
  vec3f p0, p1, p2, p3;
};

// TODO: Doesn't tell us which face it hit; same UV coordinates no matter
// what.
inline prim_intersection intersect_tetrahedron(
    const ray3f& ray, const Tetrahedron &tet) {

  //    1---3
  //   / \ /
  //  0---2

  // TODO: These all have the same winding order, but probably we want
  // to pick them so that they have consistent UV coordinates or something.
  prim_intersection p0 = intersect_triangle(ray, tet.p0, tet.p1, tet.p2);
  prim_intersection p1 = intersect_triangle(ray, tet.p1, tet.p3, tet.p2);
  prim_intersection p2 = intersect_triangle(ray, tet.p0, tet.p3, tet.p1);
  prim_intersection p3 = intersect_triangle(ray, tet.p2, tet.p3, tet.p0);

  // take the closest one
  for (const auto &p : {p1, p2, p3}) {
    if (p.hit && p.distance < p0.distance) p0 = p;
  }

  return p0;
}

struct Sphere {
  vec3f origin = {0.0f, 0.0f, 0.0f};
  float radius = 0.0f;
};

struct Prim {
  std::variant<Sphere, Tetrahedron> v;
};

// From yocto_geometry.h, but fixing a bug (?) where the UV coordinates
// are always from the back of the sphere
inline prim_intersection intersect_sphere_front(
    const ray3f& ray, const vec3f& p, float r) {
  // compute parameters
  auto a = dot(ray.d, ray.d);
  auto b = 2 * dot(ray.o - p, ray.d);
  auto c = dot(ray.o - p, ray.o - p) - r * r;

  // check discriminant
  auto dis = b * b - 4 * a * c;
  if (dis < 0) return {};

  auto MaybeIntersect = [&](float t) -> prim_intersection {
      prim_intersection isect;
      isect.hit = false;
      if (t < ray.tmin || t > ray.tmax) return isect;

      // compute local point for uvs
      auto plocal = ((ray.o + ray.d * t) - p) / r;
      auto u      = atan2(plocal.y, plocal.x) / (2 * pif);
      if (u < 0) u += 1;
      auto v = acos(clamp(plocal.z, -1.0f, 1.0f)) / pif;

      // intersect front
      isect.hit = true;
      isect.uv.x = u;
      isect.uv.y = v;
      isect.distance = t;
      return isect;
    };

  prim_intersection front = MaybeIntersect((-b - sqrt(dis)) / (2 * a));
  if (front.hit) return front;

  return MaybeIntersect((-b + sqrt(dis)) / (2 * a));
}

double OptimizeMe(double a1, double a2, double a3, double distance) {
  const double frame_width = 1920;
  const double frame_height = 1080;
  const double aspect = frame_width / (double)frame_height;

  const double near_height = NEAR_WIDTH / aspect;
  const double far_height = FAR_WIDTH / aspect;

  Sphere earth;
  earth.radius = EARTH_RADIUS;
  earth.origin = {0.0f, 0.0f, 0.0f};

  double TEX_WIDTH = 21600 * 2;
  double TEX_HEIGHT = 21600;
  double ux = 11973.0 / TEX_WIDTH;
  double uy = 5930.0 / TEX_HEIGHT;
  double uw = 65.0 / TEX_WIDTH;
  double uh = 37.0 / TEX_HEIGHT;

  auto GetUV = [&](float xf, float yf) ->
    std::optional<std::pair<float, float>> {
      float far_dist = FAR_DIST;
      float near_dist = -distance;
      float near_x = std::lerp(-NEAR_WIDTH * 0.5, NEAR_WIDTH * 0.5, xf);
      float near_y = std::lerp(-near_height * 0.5, near_height * 0.5, yf);

      float far_x = std::lerp(-FAR_WIDTH * 0.5, FAR_WIDTH * 0.5, xf);
      float far_y = std::lerp(-far_height * 0.5, far_height * 0.5, yf);

      vec3f near_pt = {near_x, near_y, -near_dist};
      vec3f far_pt = {far_x, far_y, -far_dist};

      ray3f ray;
      ray.o = near_pt;
      ray.d = far_pt - near_pt;
      ray.tmin = 0.000000001f;
      ray.tmax = 10.0f;

      mat3f rot = Rot(a1, a2, a3);
      frame3f rot_frame = make_frame(rot, {0.0, 0.0, 0.0});
      ray = transform_ray(rot_frame, ray);

      prim_intersection pi = intersect_sphere_front(
          ray, earth.origin, earth.radius);

      if (!pi.hit) return std::nullopt;
      return {make_pair(1.0f - pi.uv.x, pi.uv.y)};
    };

  auto pi00 = GetUV(0.0f, 0.0f);
  if (!pi00.has_value()) return 9999999.0f;

  auto pi10 = GetUV(1.0f, 0.0f);
  if (!pi10.has_value()) return 9999999.0f;

  auto pi01 = GetUV(0.0f, 1.0f);
  if (!pi01.has_value()) return 9999999.0f;

  auto pi11 = GetUV(1.0f, 1.0f);
  if (!pi11.has_value()) return 9999999.0f;

  auto [ux00, uy00] = pi00.value();
  auto [ux10, uy10] = pi10.value();
  auto [ux01, uy01] = pi01.value();
  auto [ux11, uy11] = pi11.value();

  auto SqDist = [](double x0, double y0,
                   double x1, double y1) {
      double xx = x0 - x1;
      double yy = y0 - y1;
      return xx * xx + yy * yy;
    };

  return
    SqDist(ux00, uy00, ux, uy) +
    SqDist(ux10, uy10, ux + uw, uy) +
    SqDist(ux01, uy01, ux, uy + uh) +
    SqDist(ux11, uy11, ux + uw, uy + uh);
}

static std::tuple<double, double, double, double> Optimize() {
  static constexpr double d = 2.801;
  auto [a1, a2, a3] =
    Opt::Minimize3D([](double a1, double a2, double a3) {
        return OptimizeMe(a1, a2, a3, d);
      },
    std::make_tuple(0.0, 0.0, 0.0),
    std::make_tuple(2 * PI, 2 * PI, 2 * PI),
    10000, 1, 1000).first;
  printf("Best:\n"
         "const double a1 = %.11f;\n"
         "const double a2 = %.11f;\n"
         "const double a3 = %.11f;\n"
         "const double distance = %.11f\n",
         a1, a2, a3, d);

  return std::make_tuple(a1, a2, a3, d);
}

struct Scene {
  std::vector<Prim> prims;


  // Could do an actual 3D rotation here, but the distance is made up
  // anyway.
# define EARTH_TILT_RAD 0.41f
# define SUN_ANGLE (2.0)
  Sphere sun = {.origin = {float(SUN_DISTANCE * sin(EARTH_TILT_RAD)),
                           float(SUN_DISTANCE * sin(SUN_ANGLE)),
                           float(SUN_DISTANCE * cos(SUN_ANGLE))},
                .radius = SUN_RADIUS};

  std::vector<std::pair<int, prim_intersection>> AllIntersections(
      const ray3f &ray_in) {
    ray3f ray = ray_in;
    std::vector<std::pair<int, prim_intersection>> hits;
    for (;;) {
      auto p = NextIntersection(&ray);
      if (!p.second.hit) return hits;
      else hits.push_back(p);
    }
  }

  // Updates ray with distance of next intersection (if any).
  std::pair<int, prim_intersection> NextIntersection(
      ray3f *ray) {
    // Obviously this should use spatial data structures if the
    // scene is big!
    int isect_idx = -1;
    prim_intersection isect;
    isect.distance = flt_max;
    isect.hit = false;
    for (int idx = 0; idx < prims.size(); idx++) {
      const Prim &p = prims[idx];
      if (const Sphere *sphere = std::get_if<Sphere>(&p.v)) {
        prim_intersection pi =
          intersect_sphere_front(*ray, sphere->origin, sphere->radius);
        if (pi.hit && pi.distance < isect.distance) {
          isect = pi;
          isect_idx = idx;
        }

      } else if (const Tetrahedron *tet = std::get_if<Tetrahedron>(&p.v)) {
        prim_intersection pi =
          intersect_tetrahedron(*ray, *tet);
        if (pi.hit && pi.distance < isect.distance) {
          isect = pi;
          isect_idx = idx;
        }

      } else {
        CHECK(false) << "Unknown prim??";
      }
    }

    // XXX some more principled epsilon; nextafter?
    if (isect.hit) ray->tmin = isect.distance + 0.00001;
    return make_pair(isect_idx, isect);
  }
};

ImageRGBA *tile = nullptr;

uint32_t EarthColor(float ux, float uy) {
  static constexpr auto tile_top = Convert::ToUV(
      40.576817, -80.179213);
  static constexpr auto tile_bot = Convert::ToUV(
      40.293648, -79.515399);

  auto In = [ux, uy](std::pair<double, double> a,
                     std::pair<double, double> b) ->
    std::optional<std::pair<double, double>> {

    auto InOrder = [](double a, double b, double c) {
        return a < c ? (a <= b && b < c) : (c <= b && b < a);
      };

    if (InOrder(a.first, ux, b.first) &&
        InOrder(a.second, uy, b.second)) {
        return std::make_pair((ux - a.first) / (b.first - a.first),
                              (uy - a.second) / (b.second - a.second));
      } else {
        return std::nullopt;
      }
    };

  /*
  static constexpr auto usa1 = Convert::ToUV(48.348262, -124.675907);
  static constexpr auto usa2 = Convert::ToUV(32.703243, -79.675909);
  if (auto pou = In(usa1, usa2)) {
    return 0xFFFF77FF;

  } else
  */

  if (auto po = In(tile_top, tile_bot)) {
    auto [x, y] = po.value();
    const auto [r, g, b, a_] = tile->SampleBilinear(
        // XXX maybe don't bother scaling down?
        x * tile->Width(),
        y * tile->Height());
    return ColorUtil::FloatsTo32(r / 255.0f, g / 255.0f, b / 255.0f, 1.0f);
  } else {
    const auto [r, g, b, a_] = bluemarble->SampleBilinear(
        ux * bluemarble->Width(),
        uy * bluemarble->Height());

    return ColorUtil::FloatsTo32(r / 255.0f, g / 255.0f, b / 255.0f, 1.0f);
  }
}


// A good distance is 20.0
// A good angle is 3.14159 + 0.25
static ImageRGBA RenderFrame(
    int frame_width,
    int frame_height,
    int oversample,
    // Distance from Earth
    double distance,
    // Rotation of Earth along its axis
    double a1, double a2, double a3) {

  ImageRGBA img(frame_width * oversample, frame_height * oversample);

  // Right-handed:
  //
  //   ^ +y   -z (into scene)
  //   |    /
  //   |  /
  //   |/
  //   *-------> +x
  //  /
  // +z

  Scene scene;
  Sphere earth;
  earth.radius = EARTH_RADIUS;
  earth.origin = {0.0f, 0.0f, 0.0f};
  scene.prims.emplace_back(earth);
  vec3f tp1 = {0, -5, 0};
  vec3f tp2 = {0, +5, 0};
  vec3f tp3 = {5, 0, +5};
  vec3f tp4 = {5, 0, -5};
  Tetrahedron mouth(tp1, tp2, tp3, tp4);
  // scene.prims.emplace_back(mouth);

  const double aspect = frame_width / (double)frame_height;

  const double near_height = NEAR_WIDTH / aspect;
  const double far_height = FAR_WIDTH / aspect;


  ParallelComp(
      /* img.Width(), */ img.Height(),
      [distance, a1, a2, a3, near_height, far_height, &img, &scene](int py) {
      double yf = py / (double)img.Height();
      for (int px = 0; px < img.Width(); px++) {
        double xf = px / (double)img.Width();

        // image plane and far plane are parallel to the
        // xy plane.
        float far_dist = FAR_DIST;
        // (maybe we want much less to zoom to a
        // region on the planet?)
        float near_dist = -distance;

        float near_x = std::lerp(-NEAR_WIDTH * 0.5, NEAR_WIDTH * 0.5, xf);
        float near_y = std::lerp(-near_height * 0.5, near_height * 0.5, yf);

        float far_x = std::lerp(-FAR_WIDTH * 0.5, FAR_WIDTH * 0.5, xf);
        float far_y = std::lerp(-far_height * 0.5, far_height * 0.5, yf);

        vec3f near_pt = {near_x, near_y, -near_dist};
        vec3f far_pt = {far_x, far_y, -far_dist};
        // vec3f far_pt = {0, 0, 0};

        ray3f ray;
        ray.o = near_pt;
        ray.d = far_pt - near_pt;
        ray.tmin = 0.000000001f;
        ray.tmax = 10.0f;

        // Rotate...
        // double axis = 0.0;
        // double lat = (40.4314779 / 90) * 3.1415926535;
        // double lon = (-80.0507117 / 90) * 3.1415926535;
        // mat3f rot = Rot(0.0, angle, 3.14159 * 0.5);
        mat3f rot = Rot(a1, a2, a3);
        // Rotate around origin.
        frame3f rot_frame = make_frame(rot, {0.0, 0.0, 0.0});
        ray = transform_ray(rot_frame, ray);


        const std::vector<std::pair<int, prim_intersection>> hits =
          scene.AllIntersections(ray);

        int num_tet = 0;
        for (const auto &[idx, pi] : hits)
          if (idx == 1) num_tet++;

        // Assume that if there are odd tetrahedron intersections,
        // it's because we started inside.
        bool in_mouth = !!(num_tet & 1);

        // 0.0f = fully shadow
        auto Illumination = [&scene](vec3f start) -> float {
            // XXX
            return 1.0f;

            // TODO: Sample sun as disc
            ray3f shadow_ray;
            vec3f dir = scene.sun.origin - start;
            shadow_ray.o = start;
            shadow_ray.d = dir;
            shadow_ray.tmin = 0.001f;
            shadow_ray.tmax = SUN_DISTANCE * 2.0f;

            const std::vector<std::pair<int, prim_intersection>> shadow_hits =
              scene.AllIntersections(shadow_ray);
            /*
            if (!shadow_hits.empty())
              printf("Shadow hits: %d\n", shadow_hits.size());
            */

            int s_num_tet = 0;
            for (const auto &[idx, pi] : shadow_hits)
              if (idx == 1) s_num_tet++;
            bool s_in_mouth = !!(s_num_tet & 1);

            for (const auto &[idx, pi] : shadow_hits) {
              if (idx == 0) {
                // sphere
                if (!s_in_mouth) return 0.0f;
              } else if (idx == 1) {
                if (s_in_mouth) {
                  // Exit mouth.
                  // Is this the inner surface of the sphere?
                  vec3f pt = shadow_ray.o + pi.distance * shadow_ray.d;
                  // XXX hard coded location of earth
                  float r = length(pt);
                  if (r <= EARTH_RADIUS) {
                    return 0.0f;
                  } else {
                    // Exiting to free space
                    s_in_mouth = false;
                  }

                } else {
                  // Enter mouth.
                  s_in_mouth = true;
                }
              } else {
                CHECK(false) << "unknown prim";
              }
            }

            // reaches sun
            return 1.0f;
          };

        auto Shade = [&](vec3f pt, uint32_t color) {
            float light = std::clamp(Illumination(pt), 0.25f, 1.0f);
            auto [r, g, b, a_] = ColorUtil::U32ToFloats(color);
            return ColorUtil::FloatsTo32(light * r,
                                         light * g,
                                         light * b,
                                         1.0f);
          };

        auto Trace = [&](){
          for (const auto &[idx, pi] : hits) {
            CHECK(pi.hit);

            if (idx == 0) {
              // sphere
              if (in_mouth) {
                // Ignore the surface of the sphere while
                // inside the mouth cutout.
              } else {
                // Normal hit on sphere.
                uint32_t color = EarthColor(1.0f - pi.uv.x, pi.uv.y);

                vec3f pt = ray.o + pi.distance * ray.d;
                color = Shade(pt, color);

                img.SetPixel32(px, py, color);
                return;
              }
            } else if (idx == 1) {

              if (in_mouth) {
                // Exit mouth.
                // Is this the inner surface of the sphere?
                vec3f pt = ray.o + pi.distance * ray.d;
                // XXX hard coded location of earth
                float r = length(pt);
                if (r <= EARTH_RADIUS) {
                  uint32_t color = ColorUtil::LinearGradient32(
                      INNER_EARTH, r / EARTH_RADIUS);

                  // But make it darker
                  // img.BlendPixel32(px, py, 0x48270677);
                  color = Shade(pt, color);
                  img.SetPixel32(px, py, color);
                  in_mouth = false;
                  return;
                } else {
                  // Exiting to free space
                  in_mouth = false;
                }

              } else {
                // Enter mouth.
                in_mouth = true;
              }
            } else {
              CHECK(false) << "unknown prim";
            }
          }

          // did not hit solid, so render space
          img.SetPixel32(px, py, 0x111100FF);
        };

        Trace();

      }
      }, THREADS);

  ImageRGBA down = img.ScaleDownBy(oversample);
  down.BlendText32(5, 5, 0xFFFF22FF,
                   StringPrintf("%.9f,%.9f,%.9f", a1, a2, a3));
  return down;
}

int main(int argc, char **argv) {

  {
    static constexpr auto usa1 = Convert::ToUV(48.348262, -124.675907);
    static constexpr auto usa2 = Convert::ToUV(32.703243, -79.675909);

    printf("%.6f,%.6f .. %.6f,%.6f\n",
           usa1.first, usa1.second,
           usa2.first, usa2.second);
  }

  //  const auto [a1, a2, a3, d_] = Optimize();
  const double a1 = 0.117112122;
  const double a2 = 3.008712215;
  const double a3 = 2.284319047;
  const double distance = 2.80100000000;

  tile = ImageRGBA::Load("tile.png");
  CHECK(tile != nullptr);

  printf("Loading marble...\n");
  // bluemarble = ImageRGBA::Load("world_shaded_43k.jpg");
  // bluemarble = ImageRGBA::Load("bluemarble.png");

  // stb_image can't decode the original (integer overflow) but two
  // hemispheres do fit. So load them individually and blit them into
  // the full sphere texture.
  string west_file = "land_shallow_topo_west.jpg";
  string east_file = "land_shallow_topo_east.jpg";

  std::unique_ptr<ImageRGBA> west, east;

  #if 1
  InParallel(
      [&](){ west.reset(ImageRGBA::Load(west_file)); },
      [&](){ east.reset(ImageRGBA::Load(east_file)); });
  // std::unique_ptr<ImageRGBA> west(ImageRGBA::Load(west_file));
  CHECK(west.get() != nullptr);
  // std::unique_ptr<ImageRGBA> east(ImageRGBA::Load(east_file));
  CHECK(east.get() != nullptr);

  CHECK(west->Height() == east->Height());
  CHECK(west->Width() == east->Width());
  Timer alloc_timer;
  bluemarble = new ImageRGBA(west->Width() * 2, west->Height());
  printf("Alloc in %.4fs\n", alloc_timer.Seconds());
  CHECK(bluemarble != nullptr);

  Timer copy_timer;
  bluemarble->CopyImage(0, 0, *west);
  bluemarble->CopyImage(west->Width(), 0, *east);
  printf("Copy in %.4fs\n", copy_timer.Seconds());

  west.reset();
  east.reset();
  printf("Done. Earth texture %d x %d\n",
         bluemarble->Width(), bluemarble->Height());
  #else
  bluemarble = ImageRGBA::Load("bluemarble.png");
  #endif

  auto FrameDistance = [](float f) {
      return std::lerp(2.80000428964, 6.0, f);
      // return std::lerp(3, 20, f);
    };
  auto FrameAngle = [a1, a2, a3](float f) {
      return std::make_tuple(a1, a2 - 0.0001 * f, a3 + 0.0001 * f);
    };

  // return std::lerp(0, 0.25, f); };

  #if 1
  Timer run_timer;
  const int NUM_FRAMES = 120;
  Asynchronously async(8);
  // static constexpr int FRAME_WIDTH = 1920;
  // static constexpr int FRAME_HEIGHT = 1080;
  static constexpr int FRAME_WIDTH = 2880;
  static constexpr int FRAME_HEIGHT = 1620;
  static constexpr int OVERSAMPLE = 2;
  for (int i = 0; i < NUM_FRAMES; i++) {
    double f = i / (double)(NUM_FRAMES - 1);
    double distance = FrameDistance(f);
    const auto [fa1, fa2, fa3] = FrameAngle(f);
    ImageRGBA img = RenderFrame(FRAME_WIDTH, FRAME_HEIGHT, OVERSAMPLE,
                                distance, fa1, fa2, fa3);
    async.Run([i, img = std::move(img)]() {
        img.Save(StringPrintf("sphere%d.png", i));
        printf("%d\n", i);
      });
  }
  double sec = run_timer.Seconds();
  printf("Wrote %d frames in %.1f sec (%.3f sec/frame).\n",
         NUM_FRAMES, sec, sec / NUM_FRAMES);
  #else
  ArcFour rc(StringPrintf("%lld", time(nullptr)));
  RandomGaussian gauss(&rc);
  static constexpr int ACROSS = 5;
  static constexpr int DOWN = 5;
  static constexpr int ONEW = 288 * 2;
  static constexpr int ONEH = 162 * 2;
  static constexpr int OVERSAMPLE = 2;

  // best
  /*
  const double a1 = 0.580167;
  const double a2 = 2.986757;
  const double a3 = 2.287507;
  */

  /*
  const double a1 = 0.486626907;
  const double a2 = 3.019623116;
  const double a3 = 2.282994048;
  */

  /*
  const double a1 = 0.50481783779;
  const double a2 = 2.99730660274;
  const double a3 = 2.28557750908;
  */
  /*
  const double a1 = 0.49672948489;
  const double a2 = 3.00654714633;
  const double a3 = 2.28445155746;
  */

  ImageRGBA mosaic(ONEW * ACROSS, ONEH * DOWN);
  for (int y = 0; y < DOWN; y++) {
    for (int x = 0; x < ACROSS; x++) {
      int i = y * DOWN + x;
      double f = i / (double)(ACROSS * DOWN - 1);
      double distance = 2.80000428964; // FrameDistance(f);
      // double angle = FrameAngle(f);

      double scale = 0.0001;
      double aa1 = i == 0 ? a1 : a1 + gauss.Next() * scale;
      double aa2 = i == 0 ? a2 : a2 + gauss.Next() * scale;
      double aa3 = i == 0 ? a3 : a3 + gauss.Next() * scale;

      ImageRGBA img = RenderFrame(ONEW, ONEH, OVERSAMPLE,
                                  distance,
                                  aa1, aa2, aa3);
      mosaic.CopyImage(x * ONEW, y * ONEH, img);
    }
  }
  mosaic.Save("mosaic.png");
  #endif


  delete bluemarble;
  return 0;
}
