
#include "simplefm2.h"

#include "util.h"
#include "base/stringprintf.h"

using namespace std;

vector<uint8> SimpleFM2::ReadInputs(const string &filename) {
  vector<pair<uint8, uint8>> twop = ReadInputs2P(filename);
  vector<uint8> ret;
  ret.reserve(twop.size());
  for (const auto &p : twop) ret.push_back(p.first);
  return ret;
}

vector<pair<uint8, uint8>> SimpleFM2::ReadInputsEx(
    const string &filename,
    vector<pair<int, string>> *subtitles) {
  vector<string> contents = Util::ReadFileToLines(filename);
  vector<pair<uint8, uint8>> out;
  for (string line : contents) {
    if (subtitles != nullptr &&
        Util::StartsWith(line, "subtitle ")) {
      string rest = line.substr(9, string::npos);
      string framenum = Util::chop(rest);
      subtitles->emplace_back(atoi(framenum.c_str()),
                              Util::LoseWhiteL(rest));
      continue;
    }

    if (line.empty() || line[0] != '|')
      continue;

    // Parse the line.
    if (line.size() < 12) {
      fprintf(stderr, "Illegal line: [%s]\n", line.c_str());
      abort();
    }

    if (!(line[1] == '0' ||
          (line[1] == '2' && out.empty()))) {
      fprintf(stderr, "Command must be zero except hard "
              "reset in first input: [%s]\n", line.c_str());
      abort();
    }

    /*
      |2|........|........||
      |0|....T...|........||
    */

    // For one-player FM2 files, sometimes the second player inputs
    // are just blank. If so, expand.
    if (line.size() > 3 &&
        line[line.size() - 1] == '|' &&
        line[line.size() - 2] == '|' &&
        line[line.size() - 3] == '|') {
      line.resize(line.size() - 2);
      line.reserve(line.size() + 10);
      line += "........||";
    }

    const string player1 = line.substr(3, 8);
    const string player2 = line.substr(12, 8);
    auto Command = [](const string &s) -> uint8 {
      uint8 command = 0;
      for (int j = 0; j < 8; j++) {
        if (s[j] != '.') {
          command |= (1 << (7 - j));
        }
      }
      return command;
    };
    out.push_back({Command(player1), Command(player2)});
  }

  return out;
}

vector<pair<uint8, uint8>> SimpleFM2::ReadInputs2P(const string &filename) {
  return ReadInputsEx(filename, nullptr);
}

vector<pair<uint8, uint8>> SimpleFM2::ExpandTo2P(const vector<uint8> &inputs) {
  vector<pair<uint8, uint8>> out;
  out.reserve(inputs.size());
  for (uint8 i : inputs) {
    out.push_back({i, 0});
  }
  return out;
}

void SimpleFM2::WriteInputs(const string &outputfile,
                            const string &romfilename,
                            const string &romchecksum,
                            const vector<uint8> &inputs) {
  WriteInputsWithSubtitles2P(outputfile, romfilename, romchecksum,
                             ExpandTo2P(inputs), {});
}

void SimpleFM2::WriteInputs2P(const string &outputfile,
                              const string &romfilename,
                              const string &romchecksum,
                              const vector<pair<uint8, uint8>> &inputs) {
  WriteInputsWithSubtitles2P(outputfile, romfilename, romchecksum,
                             inputs, {});
}

void SimpleFM2::WriteInputsWithSubtitles(
    const string &outputfile,
    const string &romfilename,
    const string &romchecksum,
    const vector<uint8> &inputs,
    const vector<pair<int, string>> &subtitles) {
  return WriteInputsWithSubtitles2P(outputfile, romfilename, romchecksum,
                                    ExpandTo2P(inputs), subtitles);
}

void SimpleFM2::WriteInputsWithSubtitles2P(
    const string &outputfile,
    const string &romfilename,
    const string &romchecksum,
    const vector<pair<uint8, uint8>> &inputs,
    const vector<pair<int, string>> &subtitles) {
  // XXX Create one of these by hashing inputs.
  string fakeguid = "FDAEE33C-B32D-B38C-765C-FADEFACE0000";
  FILE *f = fopen(outputfile.c_str(), "wb");
  fprintf(f,
          "version 3\n"
          "emuversion 9815\n"
          "romFilename %s\n"
          "romChecksum %s\n"
          "guid %s\n"
          // Read from settings?
          "palFlag 0\n"
          "NewPPU 0\n"
          "fourscore 0\n"
          "microphone 0\n"
          "port0 1\n"
          "port1 1\n"
          "port2 0\n"
          // ?
          "FDS 1\n"
          "comment author tasbot-simplefm2\n",
          romfilename.c_str(),
          romchecksum.c_str(),
          fakeguid.c_str());

  // XXX could check that subtitles are in range for the movie,
  // and sorted, and unique?
  for (const pair<int, string> &sub : subtitles)
    fprintf(f, "subtitle %d %s\n", sub.first, sub.second.c_str());

  for (size_t i = 0; i < inputs.size(); i++) {
    fprintf(f, "|%c|", (i == 0) ? '2' : '0');
    auto Controller = [f](uint8 input) {
      static constexpr char gamepad[] = "RLDUTSBA";
      for (int j = 0; j < 8; j++) {
        fprintf(f, "%c",
                (input & (1 << (7 - j))) ? gamepad[j] : '.');
      }
    };

    // 1p
    Controller(inputs[i].first);
    fprintf(f, "|");
    // 2p
    Controller(inputs[i].second);
    fprintf(f, "||\n");
  }
  fclose(f);
}

vector<pair<int, string>> SimpleFM2::MakeSparseSubtitles(
    const vector<string> &dense_subtitles) {
  vector<pair<int, string>> out;
  const string *last = nullptr;
  for (size_t i = 0; i < dense_subtitles.size(); i++) {
    if (last == nullptr || *last != dense_subtitles[i]) {
      out.emplace_back(i, dense_subtitles[i]);
    }
    last = &dense_subtitles[i];
  }
  return out;
}

string SimpleFM2::InputToString(uint8 input) {
  char f[9] = {0};
  static constexpr char gamepad[] = "RLDUTSBA";
  for (int j = 0; j < 8; j++) {
    f[j] = (input & (1 << (7 - j))) ? gamepad[j] : '.';
  }
  return (string)f;
}

string SimpleFM2::InputToColorString(uint8 input) {
  string color = "";
  string out;
  static constexpr char DOTCOLOR[] = "#999";
  static constexpr char gamepad[] = "RLDUTSBA";
  static constexpr char const *colors[] = {
    "#000",
    "#000",
    "#000",
    "#000",
    "#009",
    "#009",
    "#900",
    "#900",
  };
  for (int j = 0; j < 8; j++) {
    bool button_down = input & (1 << (7 - j));
    string this_color = button_down ? colors[j] : DOTCOLOR;
    char c = button_down ? gamepad[j] : '.';
    if (color != this_color) {
      if (color != "") out += "</span>";
      out += "<span style=\"color:" + this_color + "\">";
      color = this_color;
    }
    out += StringPrintf("%c", c);
  }
  if (color != "") out += "</span>";
  return out;
}
