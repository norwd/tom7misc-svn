
// Takes the best solutions from solutions.txt and verifies that
// they are correct; writes best-solutions.txt and play.png.
//
// With an argument, executes just one of the solutions and
// displays the steps in ascii.

#include <set>
#include <map>
#include <string>
#include <cstdint>

#include "tetris.h"
#include "encoding.h"

#include "util.h"
#include "image.h"

using namespace std;
using uint8 = uint8_t;
using uint16 = uint16_t;

using Tetris = TetrisDepth<6>;

int main(int argc, char **argv) {
  static constexpr const char *solfile = "solutions.txt";
 
  std::map<uint8_t, std::vector<Move>> sols =
    Encoding::ParseSolutions(solfile);

  if (argc == 2) {
    string v = argv[1];
    CHECK(v.size() == 2 && Util::IsHexDigit(v[0]) &&
          Util::IsHexDigit(v[1])) <<
      "Want a two-digit hex byte: " << v;
    uint8 t = Util::HexDigitValue(v[0]) * 16 +
      Util::HexDigitValue(v[1]);

    auto it = sols.find(t);
    CHECK(it != sols.end()) << "No solution for " << v <<
      " in " << solfile;

    const std::vector<Move> &movie = it->second;
    Tetris tetris;
    const uint16 full_target = Encoding::FullTarget(t);
    for (Move m : movie) {
      printf("%s", Encoding::GraphicalMoveString(m).c_str());
      printf("%s %s <- target\n\n", tetris.BoardString().c_str(),
             RowString(full_target).c_str());
      CHECK(tetris.Place(m.shape, m.col));
    }

    printf("Final:\n"
           "%s %s <- target\n\n", tetris.BoardString().c_str(),
           RowString(full_target).c_str());

  } else {
    int total_moves = 0;
    int best_moves = 99999, worst_moves = 0;
    for (const auto &[b, movie] : sols) {
      CHECK(b >= 0 && b < 256);

      Tetris tetris;
      for (Move m : movie) {
        CHECK(tetris.Place(m.shape, m.col));
      }

      int blank_prefix = 0;
      for (int r = 0; r < Tetris::MAX_DEPTH; r++) {
        if (tetris.rows[r] != 0) {
          break;
        }
        blank_prefix++;
      }

      const uint16 full_target = Encoding::FullTarget(b);

      const uint16 last_line1 = tetris.rows[Tetris::MAX_DEPTH - 4];
      const uint16 last_line2 = tetris.rows[Tetris::MAX_DEPTH - 3];
      const uint16 last_line3 = tetris.rows[Tetris::MAX_DEPTH - 2];

      CHECK(blank_prefix == Tetris::MAX_DEPTH - 4 &&
            tetris.rows[Tetris::MAX_DEPTH - 1] == full_target &&
            last_line1 == Encoding::STDPOS1 &&
            last_line2 == Encoding::STDPOS2 &&
            last_line3 == Encoding::STDPOS3) << "Supposed solution "
        "for " << (int)b << " actually made board:\n" <<
        tetris.BoardString() <<
        " " << RowString(full_target) << " <- target";


      if (!movie.empty() && Encoding::IsOkToEnd(DecodePiece(movie[0].shape))) {
        printf("Sol for %02x starts with an ending piece: %c.\n",
               b, PieceChar(DecodePiece(movie[0].shape)));
      }

      if (!Encoding::IsOkToEnd(tetris.GetLastPiece())) {
        printf("Sol for %02x doesn't end with an allowed ending piece: %c.\n",
               b, PieceChar(tetris.GetLastPiece()));
      }

      total_moves += movie.size();
      best_moves = std::min(best_moves, (int)movie.size());
      worst_moves = std::max(worst_moves, (int)movie.size());
    }
    
    printf("Total moves: %d\n"
           "Best: %d\n"
           "Worst: %d\n"
           "Average : %.2f\n",
           total_moves, best_moves, worst_moves,
           total_moves / 256.0);

    const int width = worst_moves - best_moves;
    if (width == 0) {
      printf("Can't generate image because all solutions are same length!\n");
    } else {
      ImageRGBA image(16, 16);
      for (const auto &[b, movie] : sols) {
        int size = (int)movie.size();
        float f = (size - best_moves) / (float)width;
        int x = b & 0xF;
        int y = (b >> 4) & 0xF;
        uint8 v = std::clamp((int)roundf(f * 255.0f), 0, 255);
        image.SetPixel(x, y, v, v, v, 0xFF);
      }
      image.ScaleBy(4).Save("play.png");
      printf("Wrote play.png\n");
    }

    CHECK(sols.size() == 256) << sols.size();

    Encoding::SaveSolutions("best-solutions.txt", sols);
    printf("Wrote best-solutions.txt\n");
  }

  printf("OK!\n");
  return 0;
}
