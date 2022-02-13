
#include "nes-tetris.h"

#include <string>
#include <vector>
#include <memory>
#include <cstdint>
#include <cmath>

#include "base/logging.h"
#include "base/stringprintf.h"

#include "tetris.h"

using namespace std;

// EF = empty
// 7B = white square w/shine, like in T, square, line
// 7D = blue, like in J
// 7C = cyan, like Z

static std::vector<uint8_t> StringBoard(const std::string &s) {
  CHECK(s.size() == 20 * 10);
  std::vector<uint8_t> ret(20 * 10, 0);
  for (int i = 0; i < 20 * 10; i++) {
	switch (s[i]) {
	case ' ': ret[i] = 0xEF; break;
	case '*': ret[i] = 0x7B; break;
	case '#': ret[i] = 0x7C; break;
	case '@': ret[i] = 0x7D; break;
	case '%': ret[i] = 0x77; break;
	default:
	  LOG(FATAL) << "Unknown board char " << s[i];
	}
  }
  return ret;
}

static std::string BoardString(const std::vector<uint8_t> &b) {
  CHECK(b.size() == 20 * 10);
  std::string ret(20 * 10, '?');
  for (int i = 0; i < 20 * 10; i++) {
	switch (b[i]) {
	case 0xEF: ret[i] = ' '; break;
	case 0x7B: ret[i] = '*'; break;
	case 0x7C: ret[i] = '#'; break;
	case 0x7D: ret[i] = '@'; break;
	case 0x77: ret[i] = '%'; break;
	default:
	  LOG(FATAL) << "Unknown board byte " << b[i];
	}
  }
  return ret;
}


static void TestDraw() {
  string s =
	"          "
	"          "
	"          "
	"          "
	"          "
	"          "
	"          "
	"          "
	"          "
	"          "
	"          "
	"          "
	"          "
	"          "
	"          "
	"          "
	"          "
	"   #  ****"
	"   # @@   "
	"#### @@   ";
  vector<uint8_t> board = StringBoard(s);
  CHECK(s == BoardString(board));

  DrawShapeOnBoard(0x77, I_VERT, 0, 3, &board);
  DrawShapeOnBoard(0x77, SQUARE, 8, 19, &board);
  DrawShapeOnBoard(0x7D, S_HORIZ, 3, 0, &board);
  
  // printf("%s\n", BoardString(board).c_str());
  CHECK(BoardString(board) ==
		"%  @@     "
		"%         "
		"%         "
		"%         "
		"          "
		"          "
		"          "
		"          "
		"          "
		"          "
		"          "
		"          "
		"          "
		"          "
		"          "
		"          "
		"          "
		"   #  ****"
		"   # @@ %%"
		"#### @@ %%");

}

int main(int argc, char **argv) {
  TestDraw();
  printf("OK\n");
  return 0;
}
