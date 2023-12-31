
/* Writes a WAV file containing NES sound. Only supports
   16-bit mono audio at a sample rate of XXX. This code
   is mostly the same as fceu/wave.* but wrapped in a
   class. */

// XXX should just merge this with cc-lib/wavesave.h

#ifndef _SMEIGHT_WAVE_H
#define _SMEIGHT_WAVE_H

#include <string>
#include <vector>
#include <cstdint>
#include <cstdio>

struct WaveFile {
  // Opens the file.
  WaveFile(const std::string &filename, int sample_rate);
  ~WaveFile();

  // Also closed on destruction, but this makes it explicit.
  void Close();

  void Write(const std::vector<int16_t> &vec);

 private:
  FILE *outfile;
  long wsize;
};

#endif
