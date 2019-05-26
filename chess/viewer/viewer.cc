
#include <string>
#include <vector>
#include <shared_mutex>
#include <cstdint>
#include <deque>
#include <unordered_map>
#include <unistd.h>
#include <cmath>

#include "../../cc-lib/threadutil.h"
#include "../../cc-lib/randutil.h"
#include "../../cc-lib/arcfour.h"
#include "../../cc-lib/base/logging.h"
#include "../../cc-lib/base/stringprintf.h"

#include "../chess.h"
#include "../pgn.h"
#include "../bigchess.h"
#include "../subprocess.h"
#include "../stockfish.h"
#include "../player-util.h"
#include "timer.h"

// #include "unblinder.h"
// #include "unblinder-mk0.h"

#include "SDL.h"
#include "SDL_main.h"
#include "../cc-lib/sdl/sdlutil.h"
#include "../cc-lib/sdl/font.h"
#include "../cc-lib/sdl/cursor.h"
#include "../cc-lib/lines.h"
#include "../cc-lib/util.h"

/* there are some non-ascii symbols in the font */
#define CHECKMARK "\xF2"
#define ESC "\xF3"
#define HEART "\xF4"
/* here L means "long" */
#define LCMARK1 "\xF5"
#define LCMARK2 "\xF6"
#define LCHECKMARK LCMARK1 LCMARK2
#define LRARROW1 "\xF7"
#define LRARROW2 "\xF8"
#define LRARROW LRARROW1 LRARROW2
#define LLARROW1 "\xF9"
#define LLARROW2 "\xFA"
#define LLARROW LLARROW1 LLARROW2

/* BAR_0 ... BAR_10 are guaranteed to be consecutive */
#define BAR_0 "\xE0"
#define BAR_1 "\xE1"
#define BAR_2 "\xE2"
#define BAR_3 "\xE3"
#define BAR_4 "\xE4"
#define BAR_5 "\xE5"
#define BAR_6 "\xE6"
#define BAR_7 "\xE7"
#define BAR_8 "\xE8"
#define BAR_9 "\xE9"
#define BAR_10 "\xEA"
#define BARSTART "\xEB"

#define FONTCHARS " ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789`-=[]\\;',./~!@#$%^&*()_+{}|:\"<>?" CHECKMARK ESC HEART LCMARK1 LCMARK2 BAR_0 BAR_1 BAR_2 BAR_3 BAR_4 BAR_5 BAR_6 BAR_7 BAR_8 BAR_9 BAR_10 BARSTART LRARROW LLARROW

#define FONTSTYLES 7

#define STOCKFISH_COLOR 0x339900CC

using namespace std;

using int64 = int64_t;
using uint32 = uint32_t;

template<class T>
using PositionMap = std::unordered_map<Position, T, PositionHash, PositionEq>;
using Move = Position::Move;
#define EVALUATION_THREADS 8
// #define EVALUATION_THREADS 1

#define FONTWIDTH 9
#define FONTHEIGHT 16
static Font *font = nullptr, *font2x = nullptr, *font4x = nullptr;
static Font *chessfont = nullptr, *chessfont3x = nullptr;

static SDL_Cursor *cursor_arrow = nullptr, *cursor_bucket = nullptr;
static SDL_Cursor *cursor_hand = nullptr, *cursor_hand_closed = nullptr;
static SDL_Cursor *cursor_eraser = nullptr;
#define VIDEOH 1080
#define STATUSH 128
#define SCREENW 1920
#define SCREENH (VIDEOH + STATUSH)

static SDL_Surface *screen = nullptr;

static constexpr uint32 COLORS[] = {
  // black, white
  /* 1 */ 0xFF000000,
  /* 2 */ 0xFFFFFFFF,
  // 'nice' RGB
  /* 3 */ 0xFFbd3838, // R
  /* 4 */ 0xFF0ba112, // G
  /* 5 */ 0xFF1664CE, // B
  // intense
  /* 6 */ 0xFFf943ea, // magenta
  /* 7 */ 0xFFeef943, // yellow
  /* 8 */ 0xFF29edef, // cyan
  // misc
  /* 9 */ 0xFFc77e00, // orange/brown
  // Note zero is actually like "10"
  /* 0 */ 0xFF7500c7, // purple
};

// Mode basically controls what happens when we use the mouse.
enum class Mode {
  ERASING,
  DRAWING,
  FILLING,
  CHESS,
};

namespace {
struct Evaluator {
  enum Status {
    CHECKMATE,
    DRAW,
    PLAYING,
  };
  
  struct Evaluation {
    Status status;
    // Best move, if not game over.
    Move move;
    Stockfish::Score score;
  };
  
  Evaluator() {}

  void Enqueue(const Position &position, bool at_front = false) {
    WriteMutexLock ml(&m);

    // PERF: Could move to front if at_front is true, and already
    // added but not yet running.
    bool &b = added[position];
    if (b)
      return;
    b = true;

    // fprintf(stderr, "Enqueue %s first\n", position.ToFEN(1, 1).c_str());
    // fflush(stderr);
    
    if (at_front) {
      todo.push_front(position);
    } else {
      todo.push_back(position);
    }
  };

  int64 QueueSize() {
    ReadMutexLock ml(&m);
    return todo.size();
  }

  void AddThread() {
    WriteMutexLock ml(&m);
    threads.emplace_back([this]() { WorkThread(); });
  }

  const Evaluation *GetEval(const Position &pos) {
    ReadMutexLock ml(&m);
    auto it = cache.find(pos);
    if (it == cache.end()) return nullptr;
    return it->second;
  }
  
  ~Evaluator() {
    fprintf(stderr, "Waiting for evaluators to die...\n");
    fflush(stderr);
    WriteWithLock(&m, &should_die, true);
    
    for (std::thread &t : threads)
      t.join();
    fprintf(stderr, "Done.\n");
    fflush(stderr);
  }
  
private:
  // OK to have many of these. Each creates its own private stockfish
  // process.
  void WorkThread() {
    std::unique_ptr<Stockfish> stockfish{new Stockfish(20, 500'000 /* ' */)};
    CHECK(stockfish.get() != nullptr);
    
    for (;;) {
      /*
      {
	ReadMutexLock ml(&m);
	if (should_die)
	  return;

	if (!todo.empty())
	  break;
      }
      */

      Position p;
      {
	WriteMutexLock ml(&m);
	if (should_die)
	  return;
	
	// Lost race?
	if (todo.empty()) {
	  usleep(1);
	  continue;
	}

	p = std::move(todo.front());
	todo.pop_front();
      }

      // fprintf(stderr, "process: %s\n", p.BoardString().c_str());
      
      // Do work, not holding lock.
      Evaluation *evaluation = new Evaluation;
      if (p.HasLegalMoves()) {
	string fen = p.ToFEN(1, 1);
	// fprintf(stderr, "FEN: %s\n", fen.c_str());
	// fflush(stderr);
	string movestring;
	stockfish->GetMove(fen, &movestring, &evaluation->score);
	CHECK(PlayerUtil::ParseLongMove(movestring,
					p.BlackMove(),
					&evaluation->move)) <<
	  movestring << "\n" << p.BoardString();
	// fprintf(stderr, "Move: %s\n", movestring.c_str());
	// fflush(stderr);
	evaluation->status = PLAYING;
      } else {
	evaluation->status = p.IsInCheck() ? CHECKMATE : DRAW;
      }

      {
	WriteMutexLock ml(&m);
	cache[p] = evaluation;
      }
    }
  }
  
  std::shared_mutex m;
  bool should_die = false;
  std::deque<Position> todo;
  PositionMap<bool> added;
  PositionMap<const Evaluation *> cache;

  std::vector<thread> threads;
};

struct UI {
  Mode mode = Mode::CHESS;
  bool ui_dirty = true;

  // If interp_frames > 0, then count is
  // some value in [0, frames), and incremented
  // each frame until then.
  int interp_frames = 0;
  int interp_count = 0;
  
  UI();
  void LoadMoves(const string &movedata);
  void Loop();
  void DrawStatus();
  void Draw();
  bool InSquare(int screenx, int screeny, std::pair<int, int> *square);

  bool MakeMove(Move m);
  
  void SaveUndo() {
    undo_buffer.push_back(sdlutil::duplicate(drawing));
    while (undo_buffer.size() > MAX_UNDO) {
      SDL_FreeSurface(undo_buffer.front());
      undo_buffer.pop_front();
    }

    while (!redo_buffer.empty()) {
      SDL_FreeSurface(redo_buffer.front());
      redo_buffer.pop_front();
    }
  }

  uint64 current_bitmap = 0xFFFF00000000FFFFULL;
  Position position;
  // Position is after the move. String is short move name.
  vector<std::tuple<Position, Move, string>> movie;
  // Invariant is that all moves up to this point are reflected
  // in the position.
  int movie_idx = 0;
  
  uint32 current_color = 0x77AA0000;
  
  bool draw_only_bits = false;
  bool draw_stockfish = true;
  bool draw_meter = true;

  float old_meter_value = 0.0f / 0.0f;
  
  SDL_Surface *drawing = nullptr;
  int mousex = 0, mousey = 0;
  bool dragging = false;
  // If both are non-negative, represents a currently grabbed piece on
  // the board. As row, col.
  std::pair<int, int> drag_source = {-1, -1};
  int drag_handlex = 0, drag_handley = 0;
  static constexpr int CHESSX = 64, CHESSY = 64, CHESSSCALE = 32 * 3;

  static constexpr int MAX_UNDO = 32;
  deque<SDL_Surface *> undo_buffer;
  deque<SDL_Surface *> redo_buffer;
  
  void Interpolate(int frames) {
    // If already in progress, piggy back on the current interpolation.
    if (interp_frames > 0)
      return;
    interp_frames = frames;
    interp_count = 0;
  }
  
  bool UpdateEval() {
    if (!PositionEq()(position, current_eval_position)) {
      current_eval = nullptr;
      current_eval_position = position;
    }
    
    if (current_eval == nullptr) {
      current_eval = evaluator->GetEval(current_eval_position);
      if (current_eval == nullptr) {
	evaluator->Enqueue(current_eval_position, true);
      }

      return !!current_eval;
    }

    return false;
  }

  uint64 PositionBits() const;
  
  Position current_eval_position;
  const Evaluator::Evaluation *current_eval = nullptr;
  
  std::unique_ptr<Evaluator> evaluator;
};
}  // namespace

UI::UI() {
  drawing = sdlutil::makesurface(SCREENW, SCREENH, true);
  sdlutil::ClearSurface(drawing, 0, 0, 0, 0);
  CHECK(drawing != nullptr);
  Position::ParseFEN(
      "rnbqkbnr/pppppppp/8/8/8/8/PPPPPPPP/RNBQKBNR w KQkq - 0 1",
      &position);

  evaluator.reset(new Evaluator);
  for (int i = 0; i < EVALUATION_THREADS; i++) {
    evaluator->AddThread();
  }
}

void UI::LoadMoves(const string &movestring) {
  std::vector<PGN::Move> moves;
  CHECK(PGN::ParseMoves(movestring, &moves)) << movestring;
  position = Position();
  movie_idx = 0;
  movie.clear();
  for (const PGN::Move move : moves) {
    Position::Move m;
    position.ParseMove(move.move.c_str(), &m);
    CHECK(MakeMove(m));
    evaluator->Enqueue(position);
  }
}

static void DrawThick(SDL_Surface *surf, int x0, int y0,
		      int x1, int y1, 
		      Uint32 color) {  
  static constexpr int THICKNESS = 3;
  Line<int> l{x0, y0, x1, y1};

  const int w = surf->w, h = surf->h;
  
  Uint32 *bufp = (Uint32 *)surf->pixels;
  int stride = surf->pitch >> 2;
  auto SetPixel = [color, w, h, bufp, stride](int x, int y) {
      if (x >= 0 && y >= 0 &&
	  x < w && y < h) {
	bufp[y * stride + x] = color;
      }
    };

  auto ThickPixel = [&SetPixel](int x, int y) {
      static constexpr int LO = THICKNESS >> 1;
      // static constexpr int RO = THICKNESS - LO;

      for (int xx = x - LO; xx < x - LO + THICKNESS; xx++) {
	for (int yy = y - LO; yy < y - LO + THICKNESS; yy++) {
	  SetPixel(xx, yy);
	}
      }
    };

  ThickPixel(x0, y0);
  
  for (const std::pair<int, int> point : Line<int>{x0, y0, x1, y1}) {
    const int x = point.first, y = point.second;
    ThickPixel(x, y);
  }
}

static void FloodFill(SDL_Surface *surf, int x, int y,
		      Uint32 color) {
  const int w = surf->w, h = surf->h;
  
  Uint32 *bufp = (Uint32 *)surf->pixels;
  int stride = surf->pitch >> 2;

  const Uint32 replace_color = bufp[y * stride + x];
  
  auto GetPixel = [w, h, bufp, stride, replace_color](int x, int y) {
      if (x >= 0 && y >= 0 &&
	  x < w && y < h) {
	return bufp[y * stride + x];
      } else {
	return ~replace_color;
      }
    };

  auto SetPixel = [color, w, h, bufp, stride](int x, int y) {
      if (x >= 0 && y >= 0 &&
	  x < w && y < h) {
	bufp[y * stride + x] = color;
      }
    };

  std::vector<std::pair<int, int>> todo;
  if (color != replace_color)
    todo.emplace_back(x, y);

  while (!todo.empty()) {
    int xx, yy;
    std::tie(xx, yy) = todo.back();
    todo.pop_back();

    Uint32 c = GetPixel(xx, yy);
    if (c == replace_color) {
      SetPixel(xx, yy);
      todo.emplace_back(xx - 1, yy);
      todo.emplace_back(xx + 1, yy);
      todo.emplace_back(xx, yy - 1);
      todo.emplace_back(xx, yy + 1);
    }
  }
}

// x, y in screen space
bool UI::InSquare(int x, int y, std::pair<int, int> *square) {
  if (x >= CHESSX && y >= CHESSY) {
    int bx = (x - CHESSX) / CHESSSCALE;
    int by = (y - CHESSY) / CHESSSCALE;
    if (bx < 8 && by < 8) {
      // Note: row, col is y, x
      *square = {by, bx};
      return true;
    }
  }
  return false;
}

void UI::Loop() {
  for (;;) {

    SDL_Event event;
    if (SDL_PollEvent(&event)) {
      switch (event.type) {
      case SDL_QUIT:
	printf("QUIT.\n");
	return;

      case SDL_MOUSEMOTION: {
	SDL_MouseMotionEvent *e = (SDL_MouseMotionEvent*)&event;

	const int oldx = mousex, oldy = mousey;
	
	mousex = e->x;
	mousey = e->y;

	if (dragging) {
	  switch (mode) {
	  case Mode::DRAWING:
	    DrawThick(drawing, oldx, oldy, mousex, mousey, current_color);
	    break;
	  case Mode::ERASING:
	    DrawThick(drawing, oldx, oldy, mousex, mousey, 0x00000000);
	    break;
	  default:;
	  }
	  ui_dirty = true;
	}
	break;
      }

      case SDL_KEYDOWN: {
	switch (event.key.keysym.sym) {
	case SDLK_ESCAPE:
	  printf("ESCAPE.\n");
	  return;

	case SDLK_HOME: {
	  movie_idx = 0;
	  position = Position();
	  ui_dirty = true;
	  break;
	}
	  
	case SDLK_LEFT: {
	  if (movie_idx > 0) {
	    movie_idx--;
	    if (movie_idx > 0) {
	      position = std::get<0>(movie[movie_idx - 1]);
	    } else {
	      position = Position();
	    }
	    ui_dirty = true;
	  }
	  break;
	}

	case SDLK_RIGHT: {
	  if (movie_idx < movie.size()) {
	    MakeMove(std::get<1>(movie[movie_idx]));
	    ui_dirty = true;
	  }
	  break;
	}

	case SDLK_e: {
	  mode = Mode::ERASING;
	  SDL_SetCursor(cursor_eraser);
	  ui_dirty = true;
	  break;
	}
	  
	case SDLK_d: {
	  mode = Mode::DRAWING;
	  SDL_SetCursor(cursor_arrow);
	  ui_dirty = true;
	  break;
	}

	case SDLK_f: {
	  mode = Mode::FILLING;
	  SDL_SetCursor(cursor_bucket);
	  ui_dirty = true;
	  break;
	}
	  
	case SDLK_c: {
	  mode = Mode::CHESS;
	  SDL_SetCursor(cursor_hand);
	  ui_dirty = true;
	  break;
	}

	case SDLK_q: {
	  SaveUndo();
	  sdlutil::clearsurface(drawing, 0x0);
	  ui_dirty = true;
	  break;
	}

	case SDLK_z: {
	  // XXX check ctrl?
	  if (!undo_buffer.empty()) {
	    redo_buffer.push_front(drawing);
	    drawing = undo_buffer.back();
	    undo_buffer.pop_back();
	    ui_dirty = true;
	    printf("Undo size %d\n", undo_buffer.size());
	    fflush(stdout);
	  }
	  break;
	}

	case SDLK_y: {
	  // XXX check ctrl?
	  if (!redo_buffer.empty()) {
	    undo_buffer.push_back(drawing);
	    drawing = redo_buffer.front();
	    redo_buffer.pop_front();
	    ui_dirty = true;
	  }
	  break;
	}

	case SDLK_b: {
	  draw_only_bits = !draw_only_bits;
	  ui_dirty = true;
	  break;
	}

	case SDLK_s: {
	  draw_stockfish = !draw_stockfish;
	  ui_dirty = true;
	  break;
	}
	  
	case SDLK_0:
	case SDLK_1:
	case SDLK_2:
	case SDLK_3:
	case SDLK_4:
	case SDLK_5:
	case SDLK_6:
	case SDLK_7:
	case SDLK_8:
	case SDLK_9: {
	  // Not only are these in order, but they map to their
	  // ASCII Values!
	  const int n = ((event.key.keysym.sym - SDLK_0) + 9) % 10;
	  current_color = COLORS[n];
	  ui_dirty = true;
	  break;
	}
	  
	default:;
	}
	break;
      }

      case SDL_MOUSEBUTTONDOWN: {
	// LMB/RMB, drag, etc.
	SDL_MouseMotionEvent *e = (SDL_MouseMotionEvent*)&event;
	mousex = e->x;
	mousey = e->y;
	
	dragging = true;
	if (mode == Mode::DRAWING) {
	  SaveUndo();
	  // Make sure that a click also makes a pixel.
	  DrawThick(drawing, mousex, mousey, mousex, mousey, current_color);
	  ui_dirty = true;
	} else if (mode == Mode::ERASING) {
	  SaveUndo();
	  DrawThick(drawing, mousex, mousey, mousex, mousey, 0x00000000);
	  ui_dirty = true;
	  
	} else if (mode == Mode::FILLING) {
	  SaveUndo();

	  FloodFill(drawing, mousex, mousey, current_color);
	  
	  ui_dirty = true;
	} else if (mode == Mode::CHESS) {

	  SDL_SetCursor(cursor_hand_closed);
	  ui_dirty = true;
	  
	  // Get indicated square.
	  std::pair<int, int> square;
	  if (InSquare(mousex, mousey, &square)) {
	    uint8 p = position.PieceAt(square.first, square.second);
	    if (p != Position::EMPTY &&
		(p & Position::COLOR_MASK) ==
		(position.BlackMove() ? Position::BLACK : Position::WHITE)) {
	      drag_source = square;
	      drag_handlex = mousex - (CHESSX + square.second * CHESSSCALE);
	      drag_handley = mousey - (CHESSY + square.first * CHESSSCALE);
	      break;
	    }
	  }
	  drag_source = {-1, -1};
	}
	
	break;
      }

      case SDL_MOUSEBUTTONUP: {
	// LMB/RMB, drag, etc.
	dragging = false;	

	if (mode == Mode::CHESS) {
	  std::pair<int, int> square;
	  if (InSquare(mousex, mousey, &square)) {
	    Move m;
	    m.src_row = drag_source.first;
	    m.src_col = drag_source.second;

	    m.dst_row = square.first;
	    m.dst_col = square.second;

	    // (No way to pick this, currently, but it must be filled
	    // in correctly!);
	    m.promote_to = 0;
	    if ((position.PieceAt(m.src_row, m.src_col) &
		 Position::TYPE_MASK) ==
		Position::PAWN) {
	      if (m.dst_row == 0)
		m.promote_to = Position::WHITE | Position::QUEEN;
	      else if (m.dst_row == 7)
		m.promote_to = Position::BLACK | Position::QUEEN;
	    }

	    if (!MakeMove(m)) {
	      // visual feedback?
	      fprintf(stderr, "Illegal move %s.\n",
		      Position::DebugMoveString(m).c_str());
	      fflush(stderr);
	    }
	  }

	  ui_dirty = true;
	  drag_source = {-1, -1};
	  SDL_SetCursor(cursor_hand);
	}
	break;
      }
	
      default:;
      }
    }

    if (UpdateEval()) {
      Interpolate(15);
      ui_dirty = true;
    }

    if (interp_frames > 0) {
      // If in interpolation, always draw.
      ui_dirty = true;
      interp_count++;
      if (interp_count == interp_frames) {
	// Interpolation is over.
	interp_frames = interp_count = 0;
      }
    }
    
    if (ui_dirty) {
      sdlutil::clearsurface(screen, 0xFFFFFFFF);
      Draw();
      // font->draw(180, 5, StringPrintf("%d/%d", interp_count, interp_frames));
      SDL_Flip(screen);
      ui_dirty = false;
    }
  }
  
}

bool UI::MakeMove(Move m) {
  if (!position.IsLegal(m))
    return false;

  string s = position.ShortMoveString(m);

  // Position orig_pos = position;
  position.ApplyMove(m);
  // Could add check, checkmate here.

  CHECK(movie_idx <= movie.size()) << movie_idx << " " << movie.size();
  if (movie_idx == movie.size()) {
    movie.emplace_back(position, m, s);
    movie_idx++;
  } else if (Position::MoveEq(std::get<1>(movie[movie_idx]), m)) {
    movie_idx++;
  } else {
    movie.resize(movie_idx);
    movie.emplace_back(position, m, s);
    movie_idx++;
  }
  return true;
}

void UI::DrawStatus() {
  int erasecolor = 1;
  int drawcolor = 1;
  int fillcolor = 1;
  int chesscolor = 1;
  
  switch (mode) {
  case Mode::ERASING:
    erasecolor = 2;
    break;

  case Mode::DRAWING:
    drawcolor = 2;
    break;

  case Mode::FILLING:
    fillcolor = 2;
    break;
    
  case Mode::CHESS:
    chesscolor = 2;
    break;
  }
  const string modestring =
    StringPrintf("[^3E^<]^%drase^<  "
		 "[^3D^<]^%draw^<  "
		 "[^3F^<]^%dill^<  "
		 "[^3C^<]^%dhess^<  ",
		 erasecolor,
		 drawcolor, fillcolor, chesscolor);
  font2x->draw(5, SCREENH - (FONTHEIGHT * 2) - 1, modestring);

  // Color swatches.
  switch (mode) {
  case Mode::ERASING:
  case Mode::DRAWING:
  case Mode::FILLING: {
    const int yy = SCREENH - (FONTHEIGHT * 4) - 1;
    for (int i = 0; i < 10; i++) {
      static constexpr int SWATCHWIDTH = 64;
      sdlutil::fillrect(screen, COLORS[i],
			SWATCHWIDTH * i, yy,
			SWATCHWIDTH, FONTHEIGHT * 2);
      font2x->draw(SWATCHWIDTH * i + (SWATCHWIDTH >> 1) - FONTWIDTH, yy + 1,
		   StringPrintf("%d", (i + 1) % 10));
      if ((COLORS[i] & 0x00FFFFFF) ==
	  (current_color & 0x00FFFFFF)) {
	uint32 fake_brightness = (current_color & 0xFF) +
	  ((current_color >> 8) & 0xFF) +
	  ((current_color >> 16) & 0xFF);

	uint32 outline = fake_brightness > 0x77 ? 0xFF000000 : 0xFFFFFF00;

	for (int w = 0; w < 3; w++) {
	  sdlutil::DrawBox32(screen,
			     SWATCHWIDTH * i + w, yy + w,
			     SWATCHWIDTH - w * 2, FONTHEIGHT * 2 - w * 2,
			     outline);
	}
      }
    }
    break;
  }

  case Mode::CHESS:
    font2x->draw(5, SCREENH - (FONTHEIGHT * 4) - 1,
		 "[home] [<-] [->] navigate movie, ...");
    break;
  }
  
}

struct Typewriter {
  Typewriter(Font *fon, int x, int y, int w, int h) :
    fon(fon),
    startx(x), starty(y),
    width(w), height(h) { }

  // Assumes 
  void Write(const string &s) {
    int w = fon->sizex(s);
    if (current_line_width > 0 &&
	current_line_width + w > width) {
      current_line++;
      current_line_width = 0;
    }
    fon->draw(startx + current_line_width,
	      current_line * fon->height,
	      s);
    current_line_width += w;
  }
  
  int current_line = 0;
  int current_line_width = 0;
  Font *fon = nullptr;
  const int startx = 0, starty = 0, width = 0, height = 0;
};

// XXX note no support for legal
static uint32 SquareColor(bool black, bool inmove, bool legal) {
  if (black) {
    if (inmove) return 0xFF8fb43c;
    else return 0xFF86a666;
  } else {
    if (inmove) return 0xFFd6e882;
    else return 0xFFffffdd;
  }
}

static void DrawArrow(SDL_Surface *surf,
		      int x0, int y0, int x1, int y1,
		      Uint32 color) {
  // Main stem.
  DrawThick(surf, x0, y0, x1, y1, color);
  sdlutil::DrawCircle32(surf, x1, y1, 6, color);
}

void UI::Draw() {
  // Status stuff
  DrawStatus();
  
  // On-screen stuff

  // Board first.
  // TODO: When dragging, draw legal destinations
  bool have_move = false;
  Move last_move;
  if (movie_idx - 1 >= 0) {
    last_move = std::get<1>(movie[movie_idx - 1]);
    have_move = true;
  }
  
  for (int r = 0; r < 8; r++) {
    const int yy = CHESSY + r * CHESSSCALE;
    for (int c = 0; c < 8; c++) {
      const int xx = CHESSX + c * CHESSSCALE;
      const bool black = (r + c) & 1;
      const bool inmove =
	have_move && ((r == last_move.src_row &&
		       c == last_move.src_col) ||
		      (r == last_move.dst_row &&
		       c == last_move.dst_col));
      const uint32 color = SquareColor(black, inmove, false);
      sdlutil::fillrect(screen, color, xx, yy, CHESSSCALE, CHESSSCALE);
    }
  }

  // Draw evaluation.
  if (draw_stockfish && current_eval != nullptr) {
    switch (current_eval->status) {
    case Evaluator::CHECKMATE:
      font2x->draw(5, 16, "CHECKMATE.");
      break;
    case Evaluator::DRAW:
      font2x->draw(5, 16, "DRAW.");
      break;
    case Evaluator::PLAYING: {
      Move m = current_eval->move;
      DrawArrow(screen,
		CHESSX + m.src_col * CHESSSCALE + (CHESSSCALE >> 1),
		CHESSY + m.src_row * CHESSSCALE + (CHESSSCALE >> 1),
		CHESSX + m.dst_col * CHESSSCALE + (CHESSSCALE >> 1),
		CHESSY + m.dst_row * CHESSSCALE + (CHESSSCALE >> 1),
		STOCKFISH_COLOR);
      
      font2x->draw(100, 16, Position::DebugMoveString(m));
      int value = current_eval->score.value;
      // Evaluation is always from engine perspective, but since we
      // are evaluating both sides, negate for black to instead always
      // make it white's perspective.
      if (position.BlackMove()) value = -value;
      font2x->draw(5, 16,
		   StringPrintf("%s%d",
				(current_eval->score.is_mate ? "#" : ""),
				value));
	
      break;
    }
    }
  }

  // Can draw meter if we have an up to date evaluation or
  // the old one hasn't been invalidated (nan).
  if (draw_meter &&
      (current_eval != nullptr ||
       !std::isnan(old_meter_value))) {
    
      // Draw "meter"
      const float new_meter_value = [&]{
	  if (current_eval == nullptr)
	    return old_meter_value;

	  int value = current_eval->score.value;
	  if (position.BlackMove()) value = -value;
	  
	  if (current_eval->score.is_mate)
	    return value < 0 ? 0.0f : 1.0f;
	  // Otherwise, sigmoid
	  static constexpr float stretch = -1.0f / 500.0f;
	  return 1.0f / (1.0f + expf(value * stretch));
	}();

      const float interp = interp_count / (float)interp_frames;
      const float f =
	interp_frames > 0 ? 
	((1.0f - interp) * old_meter_value) +
	(interp * new_meter_value) :
	new_meter_value;

      if (interp_frames == 0)
	old_meter_value = new_meter_value;

      constexpr int xx = CHESSX + CHESSSCALE * 8 + 3;
      constexpr int yy = CHESSY + 1;
      constexpr int ww = 6;
      constexpr int hh = CHESSSCALE * 8 - 2;
      // White is on bottom, so bigger f advantages white
      const int botheight = f * hh;
      const int topheight = hh - botheight;

      sdlutil::fillrect(screen, 0xFF000033,
			xx,
			yy,
			ww, topheight);
      sdlutil::fillrect(screen, 0xFFFFFFCC,
			xx,
			yy + topheight,
			ww, botheight);

      sdlutil::DrawBox32(screen,
			 xx - 1,
			 yy - 1,
			 ww + 2,
			 hh + 2,
			 0xFF444444);
  };
      
  auto DrawPieceAt = [this](int x, int y, uint8 piece) {
      if (draw_only_bits) {
	font4x->draw(x, y + 8, "^9" BARSTART BAR_10);
      } else {
	uint8 typ = piece & Position::TYPE_MASK;
	if (typ == Position::C_ROOK) typ = Position::ROOK;
	string str = " ";
	if ((piece & Position::COLOR_MASK) == Position::BLACK) {
	  str[0] = "?pnbrqk"[typ];
	} else {
	  str[0] = "?PNBRQK"[typ];
	}

	chessfont3x->draw(x, y, str);
      } 
    };

  for (int r = 0; r < 8; r++) {
    const int yy = CHESSY + r * CHESSSCALE;
    for (int c = 0; c < 8; c++) {
      const int xx = CHESSX + c * CHESSSCALE;
      uint8 piece = position.PieceAt(r, c);
      if ((piece & Position::TYPE_MASK) != Position::EMPTY) {
	if (drag_source.first != r || drag_source.second != c) {
	  DrawPieceAt(xx, yy + 8, piece);
	}
      }
    }
  }
  
  if (drag_source.first >= 0 && drag_source.second >= 0) {
    uint8 p = position.PieceAt(drag_source.first, drag_source.second);
    CHECK((p & Position::TYPE_MASK) != Position::EMPTY);
    DrawPieceAt(mousex - drag_handlex, mousey - drag_handley, p);
  }

  if (draw_only_bits) {
    font2x->draw(CHESSX, CHESSY + 8 * CHESSSCALE + 8,
		 StringPrintf("0x%llx", PositionBits()));
  }
  
  {
    Typewriter t{font2x, 900, 16, 1900 - 900, 1000 - 16};
    // t.Write(StringPrintf("^3Move %d. ", movie_idx));

    /*
    for (int i = 0; i < movie.size(); i++) {
      t.Write(StringPrintf("[%d] ^2%s ", i, std::get<2>(movie[i]).c_str()));
    }
    */
    
    for (int i = 0; i < movie.size(); i++) {
      int color = i < movie_idx ? (i == movie_idx - 1 ? 6 : 0) : 4;
      const string &ms = std::get<2>(movie[i]);
      if (i & 1) {
	// Black's move.
	t.Write(StringPrintf("^%d%s  ", color, ms.c_str()));
      } else {
	t.Write(StringPrintf("^%d%d. %s ", color, (i >> 1) + 1,
			     ms.c_str()));
      }
    }
  }
  
  sdlutil::blitall(drawing, screen, 0, 0);
}

uint64 UI::PositionBits() const {
  uint64 out = 0ULL;
  for (int r = 0; r < 8; r++) {
    for (int c = 0; c < 8; c++) {
      out <<= 1;
      if (position.PieceAt(r, c) != Position::EMPTY) {
	out |= 1ULL;
      }
    }
  }
  return out;
}
      
int main(int argc, char **argv) {
  /* Initialize SDL and network, if we're using it. */
  CHECK(SDL_Init(SDL_INIT_VIDEO |
		 SDL_INIT_TIMER |
		 SDL_INIT_AUDIO) >= 0);
  fprintf(stderr, "SDL initialized OK.\n");

  SDL_EnableKeyRepeat(SDL_DEFAULT_REPEAT_DELAY,
		      SDL_DEFAULT_REPEAT_INTERVAL);

  SDL_EnableUNICODE(1);

  SDL_Surface *icon = SDL_LoadBMP("../blind/unblind-icon.bmp");
  if (icon != nullptr) {
    SDL_WM_SetIcon(icon, nullptr);
  }
  
  screen = sdlutil::makescreen(SCREENW, SCREENH);
  CHECK(screen);

  font = Font::create(screen,
		      "../blind/font.png",
		      FONTCHARS,
		      FONTWIDTH, FONTHEIGHT, FONTSTYLES, 1, 3);
  CHECK(font != nullptr) << "Couldn't load font.";

  font2x = Font::CreateX(2,
			 screen,
			 "../blind/font.png",
			 FONTCHARS,
			 FONTWIDTH, FONTHEIGHT, FONTSTYLES, 1, 3);
  CHECK(font2x != nullptr) << "Couldn't load font.";

  font4x = Font::CreateX(4,
			 screen,
			 "../blind/font.png",
			 FONTCHARS,
			 FONTWIDTH, FONTHEIGHT, FONTSTYLES, 1, 3);
  CHECK(font4x != nullptr) << "Couldn't load font.";

  //# define CHESSFONT "../blind/chessfont.png"
# define CHESSFONT "chessfont-cylon.png"
  
  chessfont = Font::create(screen,
			   CHESSFONT,
			   " PNBRQKpnbrqk",
			   32, 32, 1, 0, 1);
  CHECK(chessfont != nullptr) << "Couldn't load chessfont.";

  chessfont3x = Font::CreateX(3,
			      screen,
			      CHESSFONT,
			      " PNBRQKpnbrqk",
			      32, 32, 1, 0, 1);
  CHECK(chessfont3x != nullptr) << "Couldn't load chessfont3x.";

  CHECK((cursor_arrow = Cursor::MakeArrow()));
  CHECK((cursor_bucket = Cursor::MakeBucket()));
  CHECK((cursor_hand = Cursor::MakeHand()));
  CHECK((cursor_hand_closed = Cursor::MakeHandClosed()));
  CHECK((cursor_eraser = Cursor::MakeEraser()));
  
  SDL_SetCursor(cursor_arrow);
  SDL_ShowCursor(SDL_ENABLE);
  
  UI ui;
  if (argc > 1) {
    ui.LoadMoves(argv[1]);
  }
  ui.Loop();
  
  SDL_Quit();
  return 0;
}

