
#include "tetris.h"
#include "encoding.h"

#include "util.h"

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
    for (const auto &[b, movie] : sols) {
      CHECK(b >= 0 && b < 256);

      Tetris tetris;
      for (Move m : movie) {
        CHECK(tetris.Place(m.shape, m.col));
      }

      const uint16 full_target = Encoding::FullTarget(b);

      const uint16 last_line1 = tetris.rows[Tetris::MAX_DEPTH - 4];
      const uint16 last_line2 = tetris.rows[Tetris::MAX_DEPTH - 3];
      const uint16 last_line3 = tetris.rows[Tetris::MAX_DEPTH - 2];

      CHECK(tetris.rows[Tetris::MAX_DEPTH - 1] == full_target &&
            last_line1 == Encoding::STDPOS1 &&
            last_line2 == Encoding::STDPOS2 &&
            last_line3 == Encoding::STDPOS3) << "Supposed solution "
        "for " << (int)b << " actually made board:\n" <<
        tetris.BoardString() <<
        " " << RowString(full_target) << " <- target";
    }

    CHECK(sols.size() == 256) << sols.size();
  }

  printf("OK!\n");
  return 0;
}
