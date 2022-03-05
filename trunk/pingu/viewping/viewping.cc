
#include <cstdint>
#include <string>
#include <mutex>

#include "SDL.h"

#include "sdl/sdlutil.h"
#include "sdl/font.h"
#include "sdl/chars.h"
#include "sdl/cursor.h"

#include "image64.h"
#include "base/logging.h"
#include "util.h"
#include "threadutil.h"
#include "geom/hilbert-curve.h"
#include "base/stringprintf.h"
#include "timer.h"


using namespace std;
using int64 = int64_t;
static constexpr uint8_t TIMEOUT = 0;
static constexpr uint8_t WRONG_DATA = 255;

static constexpr int SCREENW = 1920;
static constexpr int SCREENH = 1080;

std::mutex screen_mutex;
static SDL_Surface *screen = nullptr;
static Font *font = nullptr;
static SDL_Cursor *cursor_arrow = nullptr, *cursor_bucket = nullptr;
static SDL_Cursor *cursor_hand = nullptr, *cursor_hand_closed = nullptr;
static SDL_Cursor *cursor_eraser = nullptr;

// First image is the full resolution (64k x 64k) and each
// successive one is one quarter size (half on each dimension).
static std::vector<Image64RGBA *> mipmaps;



static void ReadRaw(int octet_c, Image64RGBA *huge) {
  /*
  {
    std::unique_lock<std::mutex> ul(screen_mutex);
    sdlutil::drawbox(screen, octet_c * 7, SCREENH / 2 - 10, 7, 20,
                     octet_c, 0xFF, 0x33);
    SDL_Flip(screen);
  }
  */
  
  string filename = StringPrintf("../oping/ping%d.dat", (int)octet_c);
  std::vector<uint8_t> pings = Util::ReadFileBytes(filename);
  CHECK(!pings.empty()) << filename;
  CHECK(pings.size() == 256 * 256 * 256) << "Incomplete/bad file "
                                         << filename;

  for (int a = 0; a < 256; a++) {
    for (int b = 0; b < 256; b++) {
      for (int d = 0; d < 256; d++) {
        uint64_t pos = ((uint32)a << 24) | ((uint32)b << 16) |
          ((uint32)octet_c << 8) | d;
        const auto [x32, y32] = HilbertCurve::To2D(16, pos);
        const int x64 = x32;
        const int y64 = y32;
          
        const int64 pos_in_file = (a << 16) | (b << 8) | d;
        const uint8_t b = pings[pos_in_file];
        const uint32_t color = (b == TIMEOUT || b == WRONG_DATA) ?
          0x000000FF : 0xFFFFFFFF;
        huge->SetPixel32(x64, y64, color);
      }
    }
  }

  printf(".");
  /*
  {
    std::unique_lock<std::mutex> ul(screen_mutex);
    sdlutil::FillRectRGB(screen, octet_c * 7, SCREENH / 2 - 10, 7, 20,
                         octet_c, 0xFF, 0x33);
    SDL_Flip(screen);
  }
  */
}

static void Load() {
  Timer load_timer;
  Image64RGBA *huge = new Image64RGBA(65536, 65536);
  CHECK(huge != nullptr);
  printf("Allocated image of size %lld x %lld\n",
         huge->Width(), huge->Height());
  huge->Clear32(0xFF00FFFF);

  printf("Cleared. Loading data:\n");
  ParallelComp(256, [huge](int c) { ReadRaw(c, huge); }, 8);

  mipmaps.push_back(huge);

  for (int i = 0; i < 15; i++) {
    printf("Mipmap %d/15:\n", i);
    Image64RGBA *last = mipmaps.back();
    // PERF can avoid this copy, which is significant (gigabytes)
    mipmaps.push_back(last->ScaleDownBy(2).Copy());
  }

  // mipmaps[5]->Save("mipmap5.png");
  printf("All loaded in %.1fs\n", load_timer.Seconds());
}

// Blit a rectangle from an Image64 to the screen.
static void BlitImagePart(const Image64RGBA &img,
                          int64 srcx, int64 srcy,
                          int srcw, int srch,
                          int destx, int desty) {
  for (int y = 0; y < srch; y++) {
    for (int x = 0; x < srcw; x++) {
      int sx = srcx + x;
      int sy = srcy + y;
      int dx = destx + x;
      int dy = desty + y;
      if (sx >= 0 && sx < img.Width() &&
          sy >= 0 && sy < img.Height() &&
          dy >= 0 && dy < SCREENH &&
          dx >= 0 && dx < SCREENW) {
        // auto [r, g, b, _] = img.GetPixel32(sx, sy);
        uint32 rgba = img.GetPixel32(sx, sy);
        // PERF maybe much faster to inline and skip bpp test
        // Uint32 color = SDL_MapRGB(screen->format, r, g, b);
        // Uint32 color = 0xFF | (r << 16) | (g << 8) | b;
        uint32 color = 0xFF000000 | (rgba >> 8);
        Uint32 *bufp = (Uint32 *)screen->pixels + dy*screen->pitch/4 + dx;
        *bufp = color;

        // sdlutil::drawpixel(screen, dx, dy, r, g, b);
      }
    }
  }
}

struct UI {
  int64 scrollx = 0, scrolly = 0;

  int mousex = 0, mousey = 0;

  int current_zoom = 3;
  
  // if dragging, the original scrollx,scrolly and start mouse
  // position.
  std::optional<std::tuple<int, int, int, int>> drag_origin;
  bool ui_dirty = true;

  void DrawPings() {
    // Draw from corresponding mipmap to the screen, using the
    // current scrollx,scrolly.
    BlitImagePart(*mipmaps[current_zoom],
                  scrollx, scrolly, SCREENW, SCREENH, 0, 0);
  }

  double total_draw = 0.0;
  int num_draws = 0;
  void Draw() {
    Timer draw_timer;
    sdlutil::clearsurface(screen, 0xFF001100);
    DrawPings();
    total_draw += draw_timer.Seconds();
    num_draws++;
    if (num_draws % 100 == 0) {
      printf("%d draws. Average draw time %.1fms\n",
             num_draws, (total_draw * 1000.0) / num_draws);
    }
  }

  void SetScroll(int sx, int sy) {
    scrollx = sx;
    scrolly = sy;
    
    int64 mmw = mipmaps[current_zoom]->Width();
    int64 mmh = mipmaps[current_zoom]->Height();    
    if (mmw < SCREENW) {
      scrollx = - ((SCREENW - mmw) >> 1);
    } else {
      scrollx = std::clamp(scrollx, int64(0), mmw - SCREENW);
    }

    if (mmh < SCREENH) {
      scrolly = - ((SCREENH - mmh) >> 1);
    } else {
      scrolly = std::clamp(scrolly, int64(0), mmh - SCREENH);
    }
  }

  // Supposing the mouse is at mx,my, get the normalized (in [0,1]) position
  // that the mouse is pointing at for the current zoom.
  std::pair<double, double> GetNormalizedPos(int mx, int my) {
    int64 cx = scrollx + mx;
    int64 cy = scrolly + my;
    return std::make_pair((double)(cx) / mipmaps[current_zoom]->Width(),
                          (double)(cy) / mipmaps[current_zoom]->Height());
  }

  std::pair<int64, int64> GetScrollToNormalized(int mx, int my,
                                                double fx, double fy) {
    // Return scroll such that the normalized pos fx,fy is at the
    // screen position mx,my.
    int64 cx = fx * mipmaps[current_zoom]->Width();
    int64 cy = fy * mipmaps[current_zoom]->Height();
    return std::make_pair(cx - mx, cy - my);
  }
  
  void Loop() {
    for (;;) {

      SDL_Event event;
      if (SDL_PollEvent(&event)) {
        switch (event.type) {
        case SDL_QUIT:
          printf("QUIT.\n");
          return;

        case SDL_MOUSEMOTION: {
          SDL_MouseMotionEvent *e = (SDL_MouseMotionEvent*)&event;

          mousex = e->x;
          mousey = e->y;

          if (drag_origin.has_value()) {
            // In the midst of a drag.
            auto [osx, osy, omx, omy] = drag_origin.value();
            SetScroll(osx - (mousex - omx), osy - (mousey - omy));
            
            ui_dirty = true;
          }
          break;
        }

        case SDL_KEYDOWN: {
          switch (event.key.keysym.sym) {
          case SDLK_ESCAPE:
            printf("ESCAPE.\n");
            return;

            // TODO zoom:
#if 0
          case SDLK_KP_PLUS:
          case SDLK_EQUALS:
          case SDLK_PLUS:
            if (current_value < 0xF0) current_value += 0x10;
            else current_value = 0xFF;
            ui_dirty = true;
            break;

          case SDLK_KP_MINUS:
          case SDLK_MINUS:
            if (current_value >= 0x10) current_value -= 0x10;
            else current_value = 0x00;
            ui_dirty = true;
            break;
#endif
            
          default:;
          }
          break;
        }

        case SDL_MOUSEBUTTONDOWN: {
          // LMB/RMB, drag, etc.
          SDL_MouseButtonEvent *e = (SDL_MouseButtonEvent*)&event;
          printf("Mouse button %d\n", e->button);
          mousex = e->x;
          mousey = e->y;

          if (e->button == SDL_BUTTON_LEFT) {
            drag_origin = make_optional(make_tuple(
                                            scrollx, scrolly,
                                            mousex, mousey));
            SDL_SetCursor(cursor_hand_closed);
            ui_dirty = true;
          } else if (e->button == SDL_BUTTON_WHEELUP) {
            auto [fx, fy] = GetNormalizedPos(mousex, mousey);
            
            // scroll up = zoom in
            if (current_zoom > 0) {
              current_zoom--;
              ui_dirty = true;
            }
            printf("Old f %.2f,%.2f now zoom %d\n", fx, fy, current_zoom);

            auto [sx, sy] = GetScrollToNormalized(mousex, mousey, fx, fy);
            SetScroll(sx, sy);
          } else if (e->button == SDL_BUTTON_WHEELDOWN) {
            auto [fx, fy] = GetNormalizedPos(mousex, mousey);
            if (current_zoom < mipmaps.size() - 1) {
              current_zoom++;
              ui_dirty = true;
            }

            printf("Old f %.2f,%.2f now zoom %d\n", fx, fy, current_zoom);

            auto [sx, sy] = GetScrollToNormalized(mousex, mousey, fx, fy);
            SetScroll(sx, sy);
          }

          break;
        }

        case SDL_MOUSEBUTTONUP: {
          SDL_MouseButtonEvent *e = (SDL_MouseButtonEvent*)&event;
          printf("Mouse button %d\n", e->button);
          // LMB/RMB, drag, etc.
          if (e->button == SDL_BUTTON_LEFT) {
            drag_origin.reset();
            SDL_SetCursor(cursor_arrow);
          }
          break;
        }

        default:;
        }
      }

      if (ui_dirty) {
        Draw();
        SDL_Flip(screen);
        ui_dirty = false;
      }
    }
    
  }
  
};

int main(int argc, char **argv) {

  CHECK(SDL_Init(SDL_INIT_VIDEO) >= 0);
  fprintf(stderr, "SDL initialized OK.\n");

  screen = sdlutil::makescreen(SCREENW, SCREENH);
  CHECK(screen != nullptr);

  font = Font::Create(screen,
                      "..\\..\\cc-lib\\sdl\\font.png",
                      FONTCHARS,
                      FONTWIDTH, FONTHEIGHT, FONTSTYLES, 1, 3);
  CHECK(font != nullptr) << "Couldn't load font.";
  
  CHECK((cursor_arrow = Cursor::MakeArrow()));
  CHECK((cursor_bucket = Cursor::MakeBucket()));
  CHECK((cursor_hand = Cursor::MakeHand()));
  CHECK((cursor_hand_closed = Cursor::MakeHandClosed()));
  CHECK((cursor_eraser = Cursor::MakeEraser()));
  
  Load();
  UI ui;
  ui.Loop();

  for (Image64RGBA *img : mipmaps) delete img;
  mipmaps.clear();
  
  SDL_Quit();
  return 0;
}  
