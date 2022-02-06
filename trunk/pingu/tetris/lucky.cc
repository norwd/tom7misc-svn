
#include <string>
#include <vector>
#include <memory>
#include <cstdint>
#include <cmath>
#include <unordered_set>

#include "../fceulib/emulator.h"
#include "../fceulib/simplefm2.h"
#include "../fceulib/simplefm7.h"
#include "../fceulib/x6502.h"

#include "base/logging.h"
#include "base/stringprintf.h"
#include "arcfour.h"
#include "image.h"

#include "tetris.h"
#include "nes-tetris.h"
#include "encoding.h"

using namespace std;
using uint8 = uint8_t;
using uint16 = uint16_t;

[[maybe_unused]]
static void HistoMenu(int iters) {
  ArcFour rc("tetris");

  std::unique_ptr<Emulator> emu;
  emu.reset(Emulator::Create("tetris.nes"));
  CHECK(emu.get() != nullptr);

  vector<uint8> startmovie =
    SimpleFM7::ParseString(
        "!" // 1 player
        "554_7t" // press start on main menu
        "142_"
        "7c" // signals a place where we can wait to affect RNG (maybe)
        "45_5t" // press start on game type menu
        "82_7c" // another place to potentially wait (likely?)
                           );

  vector<uint8> entergame =
    SimpleFM7::ParseString(
        "!" // 1 player
        "8_7t" // press start on A-type menu to begin game
        "68_" // game started by now
                           );

  for (uint8 c : startmovie) emu->Step(c, 0);
  for (uint8 c : entergame) emu->Step(c, 0);
  // SaveScreenshot("start.png", emu.get());

  // const std::vector<uint8> menu = emu->SaveUncompressed();

  std::vector<int> histo(NUM_PIECES * NUM_PIECES, 0);
  std::vector<int> example(NUM_PIECES * NUM_PIECES, -1);
  for (int wait = 0; wait < iters; wait++) {
    // for (int f = 0; f < wait; f++) emu->Step(0, 0);

    // (could do this incrementally instead of starting from
    // the beginning each time..)

    // save on menu
    const std::vector<uint8> menustill = emu->SaveUncompressed();

    // XXX
    // blank out what I think is the RNG state
    // (unexplained why this produces both O|T and T|L
    // multiple times? Demo mode does not increment the drop
    // counter nor modify the current drop piece, so I would
    // expect this to be totally deterministic.
    emu->SetRAM(MEM_RNG1, 0x00);
    emu->SetRAM(MEM_RNG2, 0x00);

    // And see what pieces we get if we start.
    for (uint8 c : entergame) emu->Step(c, 0);
    Piece cur = DecodePiece((Shape)emu->ReadRAM(MEM_CURRENT_PIECE));
    Piece next = DecodePiece((Shape)emu->ReadRAM(MEM_NEXT_PIECE));

    histo[next * NUM_PIECES + cur]++;
    if (example[next * NUM_PIECES + cur] == -1)
      example[next * NUM_PIECES + cur] = wait;

    if (wait % 100 == 0) printf("%d/%d\n", wait, iters);

    emu->LoadUncompressed(menustill);

    // wait one more frame.
    uint8 b = rc.Byte();
    // None of these inputs do anything on this menu
    emu->Step(b & (INPUT_L | INPUT_A | INPUT_U | INPUT_S), 0);
  }

  for (int next = 0; next < NUM_PIECES; next++) {
    for (int cur = 0; cur < NUM_PIECES; cur++) {
      const int idx = next * NUM_PIECES + cur;
      printf("%c|%c|%d times|ex %d\n",
             PieceChar((Piece)next), PieceChar((Piece)cur),
             histo[idx], example[idx]);
    }
  }
}

[[maybe_unused]]
static void PieceDrop() {

  ArcFour rc("tetris");

  std::unique_ptr<Emulator> emu;
  emu.reset(Emulator::Create("tetris.nes"));
  CHECK(emu.get() != nullptr);

  vector<uint8> begin =
    SimpleFM7::ParseString(
        "!"
        "547_6t9_5t10_7t34_"
        "5d5_5d3_5d7_5d5_5d5_4d6_4d8_4d18_4d14_4d18_3d19_4d24_"
        "6t4_" // wait here, paused
                           );
  vector<uint8> drop_piece =
    SimpleFM7::ParseString(
        "!"
        "6t" // unpause
        "134_" // wait for piece to land
                           );

  for (uint8 c : begin) emu->Step(c, 0);

  std::vector<int> histo(NUM_PIECES * NUM_PIECES, 0);
  std::vector<int> example(NUM_PIECES * NUM_PIECES, -1);
  for (int wait = 0; wait < 65536; wait++) {
    // for (int f = 0; f < wait; f++) emu->Step(0, 0);

    // (could do this incrementally instead of starting from
    // the beginning each time..)

    const std::vector<uint8> paused = emu->SaveUncompressed();

    // And see what pieces we get if we start.
    for (uint8 c : drop_piece) emu->Step(c, 0);
    Piece cur = DecodePiece((Shape)emu->ReadRAM(MEM_CURRENT_PIECE));
    Piece next = DecodePiece((Shape)emu->ReadRAM(MEM_NEXT_PIECE));

    histo[next * NUM_PIECES + cur]++;
    if (example[next * NUM_PIECES + cur] == -1)
      example[next * NUM_PIECES + cur] = wait;

    if (wait % 100 == 0) printf("%d/%d\n", wait, 65536);
    if (wait == 0) SaveScreenshot("wait0.png", emu.get());

    emu->LoadUncompressed(paused);

    // wait one more frame.
    // uint8 b = rc.Byte();
    emu->Step(0, 0);
  }

  for (int next = 0; next < NUM_PIECES; next++) {
    for (int cur = 0; cur < NUM_PIECES; cur++) {
      const int idx = next * NUM_PIECES + cur;
      printf("%c|%c|%d times|ex %d\n",
             PieceChar((Piece)next), PieceChar((Piece)cur),
             histo[idx], example[idx]);
    }
  }
}

// It doesn't work... we get stuck!
[[maybe_unused]]
static void MakeStreak1() {

  ArcFour rc("tetris");

  std::unique_ptr<Emulator> emu;
  emu.reset(Emulator::Create("tetris.nes"));
  CHECK(emu.get() != nullptr);

  vector<uint8> startmovie =
    SimpleFM7::ParseString(
        "!" // 1 player
        "554_7t" // press start on main menu
        "142_"
        "7c" // signals a place where we can wait to affect RNG (maybe)
        "45_5t" // press start on game type menu
        "82_7c" // another place to potentially wait (likely?)
        "81_" // computed wait
        "8_7t" // press start on A-type menu to begin game
        "68_" // game started by now
                           );

  for (uint8 c : startmovie) emu->Step(c, 0);
  SaveScreenshot("start.png", emu.get());

  vector<uint8> outmovie = startmovie;

  vector<uint8> savestate;
  int savelen = 0;
  auto Save = [&emu, &outmovie, &savestate, &savelen]() {
      savelen = outmovie.size();
      savestate = emu->SaveUncompressed();
    };

  auto Restore = [&emu, &savestate, &savelen, &outmovie]() {
      outmovie.resize(savelen);
      emu->LoadUncompressed(savestate);
    };

  Save();

  // Now, repeatedly...

  int target_column = 0;
  uint8 prev_counter = 0x7F;
  int retry_count = 0;
  // just one of these.
  bool square_ok = true;
  for (;;) {
    Shape cur = (Shape)emu->ReadRAM(MEM_CURRENT_PIECE);
    Shape next_byte = (Shape)emu->ReadRAM(MEM_NEXT_PIECE);
    Piece next = DecodePiece(next_byte);
    uint8 current_counter = emu->ReadRAM(MEM_DROP_COUNT);

    #if 0
    printf("%d frames. prev %d cur %d. cur %02x next %02x\n",
           (int)outmovie.size(),
           prev_counter, current_counter, cur, next_byte);
    #endif

    if (prev_counter != current_counter) {
      // piece has landed.
      const bool next_ok = (square_ok && next == PIECE_O) ||
        next == PIECE_I;
      if (!next_ok) {
        // We didn't get I as expected.
        // Go back to previous piece. Maybe should increment
        // some failure count and eventually start pausing?
        retry_count++;
        if (retry_count > 64) {
          printf("%d retries on piece %d. got next=%02x :/\n",
                 retry_count, current_counter, next_byte);
        }
        Restore();
        continue;
      }

      // Otherwise, new checkpoint.
      Save();
      retry_count = 0;
      square_ok = false;

      prev_counter = current_counter;
      if (true || prev_counter % 8 == 0) {
        SaveScreenshot(StringPrintf("lucky%d.png", prev_counter),
                       emu.get());

        SimpleFM2::WriteInputs("lucky.fm2",
                               "tetris.nes",
                               "base64:Ww5XFVjIx5aTe5avRpVhxg==",
                               outmovie);
        printf("Saved after dropping %d.\n", prev_counter);
      }

      // This also works for the two squares we place at the
      // beginning, which end up being equivalent to two I.
      target_column++;
      target_column %= 10;
    }

    uint8 cur_x = emu->ReadRAM(MEM_CURRENT_X);
    uint8 input = 0;
    switch (cur) {
    case 0x0A:
      // We start with two squares, which go in the
      // leftmost column. But this is column "1" for
      // squares.
      if (cur_x > 1) input = INPUT_L;
      break;
    case 0x12:
      // Bar but in wrong orientation.
      input = INPUT_B;
      break;
    case 0x11:
      if (cur_x > target_column) input = INPUT_L;
      else if (cur_x < target_column) input = INPUT_R;
      break;
    default:
      SaveScreenshot("error.png", emu.get());
      SimpleFM2::WriteInputs("error.fm2",
                             "tetris.nes",
                             "base64:Ww5XFVjIx5aTe5avRpVhxg==",
                             outmovie);
      LOG(FATAL) << StringPrintf(
          "counter %d retries %d target %d cur %02x next %02x\n",
          current_counter, retry_count, target_column,
          cur, next_byte);
    }

    if (input == 0) {
      input = rc.Byte() & INPUT_D;
      // If we're having trouble, even try pausing
      if (retry_count > 64) input |= rc.Byte() & INPUT_T;
    }

    emu->Step(input, 0);
    outmovie.push_back(input);
  }
}

[[maybe_unused]]
static void MakeStreak2() {

  std::vector<Move> schedule = {
    // piece type, desired orientation, column
    Move(SQUARE, 0),
    Move(L_UP, 2),
    Move(SQUARE, 0),
    Move(J_UP, 4),
    Move(L_DOWN, 2),
    Move(J_DOWN, 4),
    Move(L_UP, 6),
    Move(J_UP, 8),
    Move(L_DOWN, 6),
    Move(J_DOWN, 8),
  };

  ArcFour rc("tetris");

  std::unique_ptr<Emulator> emu;
  emu.reset(Emulator::Create("tetris.nes"));
  CHECK(emu.get() != nullptr);

  vector<uint8> startmovie =
    SimpleFM7::ParseString(
        "!" // 1 player
        "554_7t" // press start on main menu
        "142_"
        "7c" // signals a place where we can wait to affect RNG (maybe)
        "45_5t" // press start on game type menu
        "82_7c" // another place to potentially wait (likely?)
        // "81_" // computed wait
        "13_" // computed wait for L|O, to match schedule below
        "8_7t" // press start on A-type menu to begin game
        "68_" // game started by now
                           );

  for (uint8 c : startmovie) emu->Step(c, 0);
  SaveScreenshot("start.png", emu.get());

  vector<uint8> outmovie = startmovie;

  vector<uint8> savestate;
  int savelen = 0;
  auto Save = [&emu, &outmovie, &savestate, &savelen]() {
      savelen = outmovie.size();
      savestate = emu->SaveUncompressed();
    };

  auto Restore = [&emu, &savestate, &savelen, &outmovie]() {
      outmovie.resize(savelen);
      emu->LoadUncompressed(savestate);
    };

  Save();

  // Now, repeatedly...

  // This starts at 2 because the current and next count as drops.
  uint8 prev_counter = 0x02;
  int retry_count = 0;
  int schedule_idx = 0;
  int64 frame = 0;
  int pieces = 0;
  for (;;) {
    // (note this is not really a frame counter currently, as it
    // does not get saved/restored)
    frame++;

    Shape cur = (Shape)emu->ReadRAM(MEM_CURRENT_PIECE);
    Shape next_byte = (Shape)emu->ReadRAM(MEM_NEXT_PIECE);
    Piece next = DecodePiece(next_byte);
    uint8 current_counter = emu->ReadRAM(MEM_DROP_COUNT);

    #if 0
    printf("%d frames sched %d. prev %d cur %d. cur %02x next %02x\n",
           (int)outmovie.size(), schedule_idx,
           prev_counter, current_counter, cur, next_byte);
    #endif

    if (cur == 0x13) {
      // Wait for lines to finish. Not a valid piece so we don't
      // want to enter the code below.
      // TODO: tetrises?
      emu->Step(0, 0);
      outmovie.push_back(0);
      continue;
    }

    if (prev_counter != current_counter) {
      // A piece just landed.
      // This means that 'cur' should just now have the next piece that
      // we previously established as correct.
      CHECK(DecodePiece(cur) ==
            DecodePiece(schedule[(schedule_idx + 1) %
                                 schedule.size()].shape));;

      // Now check that the NEW next piece is what we expect.
      const Piece expected =
        DecodePiece(schedule[(schedule_idx + 2) % schedule.size()].shape);
      if (next != expected) {
        // We didn't get the expected piece.
        // Go back to previous piece. Also keep track of how many
        // times this has recently happened.
        retry_count++;
        if (retry_count > 64) {
          static bool first = true;
          if (first) SaveScreenshot("stuck.png", emu.get());
          first = false;

          const uint8 rng1 = emu->ReadRAM(MEM_RNG1);
          const uint8 rng2 = emu->ReadRAM(MEM_RNG2);
          printf("%d tries on piece %d (sched %d). cur: %02x=%c got next=%c want %c RNG %02x%02x fs %d\n",
                 retry_count, current_counter, schedule_idx,
                 cur, PieceChar(DecodePiece(cur)),
                 PieceChar(DecodePiece(next_byte)),
                 PieceChar(expected),
                 rng1, rng2,
                 (int)outmovie.size());
        }
        Restore();
        continue;
      }

      // Otherwise, new checkpoint.
      Save();
      retry_count = 0;
      schedule_idx++;
      schedule_idx %= schedule.size();
      pieces++;

      prev_counter = current_counter;
      if (pieces % 25 == 0) {
        SaveScreenshot(StringPrintf("lucky%d.png", pieces),
                       emu.get());

        SimpleFM2::WriteInputs("lucky.fm2",
                               "tetris.nes",
                               "base64:Ww5XFVjIx5aTe5avRpVhxg==",
                               outmovie);
        printf("Saved after dropping %d pieces.\n", pieces);
      }
    }

    uint8 cur_x = emu->ReadRAM(MEM_CURRENT_X);
    Shape cur_shape = (Shape)emu->ReadRAM(MEM_CURRENT_PIECE);
    uint8 input = 0;

    const Move move = schedule[schedule_idx];
    // The shape should only mismatch due to orientation.
    CHECK(DecodePiece(cur_shape) == DecodePiece(move.shape)) <<
      StringPrintf("Expect %c but have %c=%02x?",
                   PieceChar(DecodePiece(move.shape)),
                   PieceChar(DecodePiece(cur_shape)), cur_shape);

    const uint8 target_nes_x =
      move.col + ShapeOffset(move.shape);

    if ((frame % 2) == 0) {
      // PERF: Can rotate in the most efficient direction.
      if (move.shape != cur_shape)
        input |= INPUT_A;

      if (target_nes_x < cur_x) input |= INPUT_L;
      else if (target_nes_x > cur_x) input |= INPUT_R;
    }

    if (move.shape == cur_shape &&
        target_nes_x == cur_x) {
      // printf("OK %d %d\n", orientation, column);
      input = rc.Byte() & INPUT_D;
      // If we're having trouble, even try pausing
      if (retry_count > 64) input |= rc.Byte() & INPUT_T;
    }

    emu->Step(input, 0);
    outmovie.push_back(input);
  }
}


[[maybe_unused]]
static void HistoSimulate(int iters) {
  std::vector<int> histo(NUM_PIECES * NUM_PIECES, 0);
  std::vector<int> example(NUM_PIECES * NUM_PIECES, -1);

  RNGState state;
  for (int wait = 0; wait < iters; wait++) {


    RNGState p1 = NextPiece(state);
    RNGState p2 = NextPiece(p1);

    uint8 cur = DecodePiece((Shape)p1.last_drop);
    uint8 next = DecodePiece((Shape)p2.last_drop);

    histo[next * NUM_PIECES + cur]++;
    if (example[next * NUM_PIECES + cur] == -1)
      example[next * NUM_PIECES + cur] = wait;

    state = NextRNG(state);
  }

  for (int next = 0; next < NUM_PIECES; next++) {
    for (int cur = 0; cur < NUM_PIECES; cur++) {
      const int idx = next * NUM_PIECES + cur;
      printf("%c|%c|%d times|ex %d\n",
             PieceChar((Piece)next), PieceChar((Piece)cur),
             histo[idx], example[idx]);
    }
  }
}

struct RngHisto {
  RngHisto() : count8(256, 0), counts(65536, 0) {}
  std::vector<int> count8;
  std::vector<int> counts;

  void Increment(uint8 rng1, uint8 rng2) {
    uint16 rng_state = ((uint16)rng1 << 8) | rng2;
    counts[rng_state]++;
    count8[rng1]++;
  }

  void Save(const string &filename) {
    int mosti = 0, most = 0;
    for (int i = 0; i < 65536; i++) {
      if (counts[i]) {
        if (counts[i] > most) {
          most = counts[i];
          mosti = i;
        }
        printf("%02x: %d\n", i, counts[i]);
      }
    }
    printf("Most at idx %02x: %d\n", mosti, most);

    ImageRGBA img(256, 257);
    for (int y = 0; y < 256; y++) {
      for (int x = 0; x < 256; x++) {
        int idx = y * 256 + x;
        if (counts[idx]) {
          uint8 c = std::clamp(counts[idx] * 16, 0x60, 0xFF);
          img.SetPixel(x, y, c, c, c, 0xFF);
        } else {
          img.SetPixel32(x, y, 0x000000FF);
        }
      }
    }

    int max8 = 0;
    for (int c : count8) max8 = std::max(c, max8);
    for (int x = 0; x < 256; x++) {
      if (count8[x]) {
        double f = count8[x] / (double)max8;
        uint8 r = std::clamp((int)round(f * (255 - 60)), 60, 255);
        img.SetPixel(x, 256, 0, r, 0, 0xFF);
      } else {
        img.SetPixel32(x, 256, 0x000000FF);
      }
    }
    img.ScaleBy(6).Save(filename);
  }
};

[[maybe_unused]]
static void HistoRNG() {

  ArcFour rc("tetris");

  std::unique_ptr<Emulator> emu;
  emu.reset(Emulator::Create("tetris.nes"));
  CHECK(emu.get() != nullptr);

  vector<uint8> startmovie =
    SimpleFM7::ParseString(
        "!" // 1 player
        "554_7t" // press start on main menu
        "142_"
        "7c" // signals a place where we can wait to affect RNG (maybe)
        "45_5t" // press start on game type menu
        "82_"
        "7t82_"); // and enter game...

  for (uint8 c : startmovie) emu->Step(c, 0);

  RngHisto histo;
  for (int iters = 0; iters < 131072; iters++) {
    histo.Increment(emu->ReadRAM(MEM_RNG1),
                    emu->ReadRAM(MEM_RNG2));
    emu->Step(0, 0);
    if (iters % 1000 == 0) printf("%d/65536\n", iters);
  }
  SaveScreenshot("actualnes-screen.png", emu.get());

  histo.Save("actualnes.png");
}

[[maybe_unused]]
static void HistoSimulatedRNG() {

  std::vector<int> counts(65536, 0);
  RNGState state;
  for (int warmup = 0; warmup < 1024; warmup++) {
    state = NextRNG(state);
  }

  RngHisto histo;
  for (int iters = 0; iters < 131072; iters++) {
    histo.Increment(state.rng1, state.rng2);
    state = NextRNG(state);
    if (iters % 1000 == 0) printf("%d/65536\n", iters);
  }

  histo.Save("simulated.png");
}

static const vector<std::tuple<float, float, float, float>> rgb_ramp = {
  { 0.0f,  1.0f, 0.0f, 0.0f},
  { 0.25f, 1.0f, 1.0f, 0.0f},
  { 0.5f,  0.0f, 1.0f, 0.0f},
  { 0.75f, 0.0f, 1.0f, 1.0f},
  { 1.0f,  0.0f, 0.0f, 1.0f},
};

static std::tuple<float, float, float>
LinearRamp(float t,
           const vector<std::tuple<float, float, float, float>> &ramp) {
  CHECK(!ramp.empty());
  auto prev = ramp[0];

  {
    const auto [x, r, g, b] = prev;
    if (t < x) {
      return make_tuple(r, g, b);
    }
  }

  for (int i = 1; i < (int)ramp.size(); i++) {
    const auto now = ramp[i];
    const auto [px, pr, pg, pb] = prev;
    const auto [x, r, g, b] = now;
    if (t < x) {
     // linear interpolation
      const float w = x - px;
      const float f = (t - px) / w;
      const float omf = 1.0f - f;
      return make_tuple(f * r + omf * pr,
                        f * g + omf * pg,
                        f * b + omf * pb);
    }
    prev = now;
  }

  {
    const auto [x, r, g, b] = prev;
    return make_tuple(r, g, b);
  }
}

// hard to read plot of the cycle
[[maybe_unused]]
static void PlotSimulatedRNG() {
  constexpr int PERIOD = 32767;
  constexpr int SCALE = 7;
  ImageRGBA img(SCALE * 256, SCALE * 256);
  img.Clear32(0x000000FF);
  RNGState state;

  auto C = [](float f) -> uint8 {
      return std::clamp(std::roundf(f * 255.0f), 0.0f, 255.0f);
    };

  int prev_x = state.rng1 * SCALE + SCALE / 2;
  int prev_y = state.rng2 * SCALE + SCALE / 2;
  for (int iters = 0; iters < PERIOD; iters++) {
    state = NextRNG(state);
    double f = iters / (double)PERIOD;
    const auto [rf, gf, bf] = LinearRamp(f, rgb_ramp);
    int x = state.rng1 * SCALE + SCALE / 2;
    int y = state.rng2 * SCALE + SCALE / 2;

    img.BlendLineAA(prev_x, prev_y, x, y,
                    C(rf), C(gf), C(bf), 0x33);

    prev_x = x;
    prev_y = y;
  }

  const RNGState start_state;
  CHECK(state.rng1 == start_state.rng1 &&
        state.rng2 == start_state.rng2);

  img.Save("simulated-plot.png");
}

[[maybe_unused]]
static void HeatSimulatedRNG() {
  constexpr int PERIOD = 32767;
  constexpr int SCALE = 7;
  ImageRGBA img(256, 256);
  img.Clear32(0x000000FF);
  RNGState state;

  auto C = [](float f) -> uint8 {
      return std::clamp(std::roundf(f * 255.0f), 0.0f, 255.0f);
    };

  for (int iters = 0; iters < PERIOD; iters++) {
    state = NextRNG(state);
    double f = iters / (double)PERIOD;
    const auto [rf, gf, bf] = LinearRamp(f, rgb_ramp);
    int x = state.rng1;
    int y = state.rng2;

    img.SetPixel(x, y, C(rf), C(gf), C(bf), 0xFF);
  }

  const RNGState start_state;
  CHECK(state.rng1 == start_state.rng1 &&
        state.rng2 == start_state.rng2);

  img.ScaleBy(SCALE).Save("simulated-heat.png");
}

[[maybe_unused]]
static void PrintSequence() {
  constexpr int PERIOD = 32767;
  RNGState state;

  for (int iters = 0; iters < PERIOD; iters++) {
    printf("%d. %d %d\n", iters, state.rng1, state.rng2);
    state = NextRNG(state);
  }

  const RNGState start_state;
  CHECK(state.rng1 == start_state.rng1 &&
        state.rng2 == start_state.rng2);
}

[[maybe_unused]]
static void HeatPair() {
  constexpr int PERIOD = 32767;
  constexpr int SCALE = 7;
  ImageRGBA img(256 * SCALE, 256 * SCALE);
  img.Clear32(0x000000FF);
  RNGState state;

  int prev = state.rng1;
  for (int iters = 0; iters < PERIOD; iters++) {
    state = NextRNG(state);
    int cur = state.rng1;

    img.BlendBox32(cur * SCALE + 1, prev * SCALE + 1,
                   SCALE - 2, SCALE - 2,
                   0xFFFFFFFF, 0x777777FF);
    img.BlendRect32(cur * SCALE + 2, prev * SCALE + 2,
                    SCALE - 4, SCALE -4, 0xFFFFFFFF);
    prev = cur;
  }

  const RNGState start_state;
  CHECK(state.rng1 == start_state.rng1 &&
        state.rng2 == start_state.rng2);

  img.Save("heat-pair.png");
}


// Why do I get a checkerboard pattern for the NES
// version if I don't play any inputs at start, and
// the same for my transcribed code -- but get a
// more degenerate looking pattern (more like what I
// expected for this routine) if I first play
// the movie to go to the game type menu??
// - could be that the flags entering the routine
//   make a difference? On the title screen I see
//   UV set, CZ clear.
// - same on game type menu
// - ok but on A-type menu, the carry flag is
//   definitely set sometimes!
// - buuuut, this routine always does a CLC so it
//   should not matter.
// - another possibility is that it's calling
//   the RNG many times on the frame, dependent
//   on the rng state. so this only leaves certain
//   RNG states.
// - hmm yes that seems to be true, as this phenomenon
//   goes away during gameplay.

// Even though the first byte takes on all possible values and "looks
// random", the dual-roll approach in NextPiece is not well
// distributed:
//
// Given an rng byte, only two values are possible for the next byte,
// (rng >> 1) and (128 + rng >> 1).
//
// In the re-roll scenario, we have that (rng + drop_count) & 0b111 is
// either 0b111 or the last drop index. The next roll is one of
// ((rng >> 1) + drop_count) & 0b111  or
// (128 + (rng >> 1) + drop_count) & 0b111. These expressions are
// congruent mod 8, so we really just have
// ((rng >> 1) + drop_count) & 0b111.
//  (could complete this analytical approach but it's also easy
//   to just sample all states!)
//
// As a result (see AllRerolls), even if you have complete control
// over the RNG state (but not prev piece and drop count), you
// have at most 3 different outcomes possible in the reroll case.
// For some setups, there is only 1 (rare) or 2 possible rerolls.
// If one of those isn't the drop_piece, too bad! This happens in
// 14.29% of cases for most pieces, except Z where it is 10.27%.
//
// So our easiest strategies would be ones where we alternate pieces,
// since in any state we can get a piece different from the last drop
// on the first roll. Another possibility would be to detect the 14%
// of cases where it is impossible to continue the streak of our
// desired piece (this can be a table or even analytical?) and instead
// manage to roll a garbage piece in those cases that we have some
// strategy for dealing with.

// The only place that writes these bytes outside of
// initialization really appears to be the RNG
// routine. (Use "Forbid" breakpoint range on the
// RNG code, in addition to a "write" breakpoint
// on the rng state.)


struct TinySet {
  // Perhaps a faster data structure since max 3 items!
  std::unordered_set<int> items;
  void Add(int item) {
    items.insert(item);
    if (items.size() > 3) {
      printf("Too many items for %d,%d:\n", z, last_piece);
      for (int item : items)
        printf("%d, ", item);
      LOG(FATAL) << "just added " << item;
    }
  }

  TinySet(int z, int last_piece) : z(z), last_piece(last_piece) {}
  const int z, last_piece;
};

// TODO: Worst case "droughts" for getting a certain piece in a state?
// It could be that rng1 stays even for half the period, as one simple
// example. This would be 8 minutes of real time at 60hz!
[[maybe_unused]]
static void AllRerolls() {
  constexpr std::array<uint8, 8> PIECES = {
    0x02, 0x07, 0x08, 0x0A, 0x0B, 0x0E, 0x12,
    // not used ?
    0x02,
  };

  // for any counter / drop
  int histo_states[4] = {};
  int max_piece[7] = {};
  int min_piece[7] = {9, 9, 9, 9, 9, 9, 9};
  int cant_dup[7] = {};
  for (uint8 last_drop = 0; last_drop < 7; last_drop++) {
    for (int z = 0; z < 256; z++) {
      const uint8 tmp_drop_count = z;

      // show that the random state doesn't even matter

      TinySet tes(z, last_drop);

      for (int x = 0; x < 256; x++) {
        const uint8 tmp_rng1 = x;
        const uint8 first_roll = (tmp_rng1 + (tmp_drop_count + 1)) & 0x7;

        if (first_roll == last_drop || first_roll == 0x7) {
          // These are the re-roll scenarios.
          // Try with all possible RNG states.
          for (int y = 0; y < 256; y++) {
            if (x == 0 && y == 0) continue;

            RNGState state;
            state.rng1 = x;
            state.rng2 = y;
            state.drop_count = z;
            state.last_drop = last_drop;
            RNGState after = NextPiece(state);

            tes.Add(after.last_drop);
          }
        }
      }

      histo_states[tes.items.size()]++;
      max_piece[last_drop] = std::max(max_piece[last_drop],
                                      (int)tes.items.size());
      min_piece[last_drop] = std::min(min_piece[last_drop],
                                      (int)tes.items.size());

      if (tes.items.find(last_drop) == tes.items.end())
        cant_dup[last_drop]++;
    }
  }
  printf("OK\n");
  printf("Count with 0: %d, 1: %d, 2: %d, 3: %d\n",
         histo_states[0], histo_states[1], histo_states[2], histo_states[3]);

  for (int p = 0; p < 7; p++) {
    char pc = PieceChar(DecodePiece((Shape)PIECES[p]));
    printf("Best case for piece %d=%c: %d, worst: %d, can't dup: %d/%d = %.02f%%\n",
           p, pc, max_piece[p], min_piece[p], cant_dup[p],
           256 * 7, (100.0 * cant_dup[p]) / (256 * 7));
  }
}

int main(int argc, char **argv) {
  MakeStreak2();

  return 0;
}
