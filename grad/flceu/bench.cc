
#include "emulator.h"

#include <bit>
#include <cmath>
#include <string>
#include <vector>
#include <memory>
#include <sys/time.h>
#include <sstream>
#include <unistd.h>
#include <cstdio>

#include "base/logging.h"
#include "base/stringprintf.h"
#include "test-util.h"
#include "rle.h"
#include "simplefm2.h"
#include "simplefm7.h"
#include "timer.h"
#include "image.h"
#include "periodically.h"
#include "ansi.h"

#include "opcodes.h"
#include "x6502.h"

#define ROMDIR "../../fceulib/roms/"
static constexpr const char *ROMFILE = ROMDIR "mario.nes";

// Beats Mario 1-1.
static constexpr const char *MOVIE =
  "!14_6t5_4t3_3t6_4t2_3t4_3t11_,b7rb4b228rb21+a21rb3+a37rb45+a20rb"
  "29+a20rb7+a5rb2+a9rb20+a2rb27+a84rb,+a2rb11+a45rb15+a35rb14+a49rb"
  "28+a13rb13+a52rb13+a46rb26+a44rb25+a14rb11+a44rb20+a15rb19+a45rb"
  "30+a36rb16+a39rb25+a2rb,b11lb,b8ba13+r21rb16+a23rb6b3lb5b15lb10b"
  "669_";
static constexpr uint64 expected_nes = 0x8eed7b59ec6a1163ULL;
static constexpr uint64 expected_img = 0xad90a1030e98e8e5ULL;
static constexpr uint64 expected_inst = 0xb106a13670bbd088ULL;
static constexpr uint64 expected_mem = 0x70db77bbffa5a4c6ULL;

static uint64 InstHistoChecksum(Emulator *emu) {
  uint64 out = 0xCAFE1234C0FFEEEEULL;
  const X6502 *x = emu->GetFC()->X;
  for (int i = 0; i < 256; i++) {
    // printf("%02x: %lld (%llx)\n", i, x->inst_histo[i], out);
    out += x->inst_histo[i];
    // Random number
    out *= 0x45c1d4edd839ff81ULL;
    out = std::rotr(out, 17);
  }
  return out;
}

static uint64 MemTrace(Emulator *emu) {
  return emu->GetFC()->X->mem_trace;
}

std::tuple<uint64, uint64, uint64, uint64, double> RunBenchmark(
    Emulator *emu,
    const vector<uint8> &start,
    const vector<uint8> &movie) {
  emu->LoadUncompressed(start);
  emu->GetFC()->X->ClearInstHisto();
  emu->GetFC()->X->ClearMemTrace();
  Timer exec_timer;
  Periodically slow_per(10.0);

  // Only the last step needs to be full, so that we render the image.
  for (int i = 0; i < (int)movie.size() - 1; i++) {
    emu->Step(movie[i], 0);
    /*
    if (i % 600 == 0) {
      ImageRGBA img(emu->GetImage(), 256, 256);
      img.Save(StringPrintf("mario-%d.png", i));
      printf("Wrote %d\n", i);
    }
    */
    if (i % 50 == 0) {
      if (slow_per.ShouldRun()) {
        double sec = exec_timer.Seconds();
        double fps = i / sec;
        double left = ((int)movie.size() - i) / fps;
        printf("%d/%d frames in %.2fs (" ACYAN("%.4f") " fps) eta "
               ARED("%.1f") "s\n",
               i, (int)movie.size(),
               sec, i / sec, left);
      }
    }
  }
  emu->StepFull(movie[movie.size() - 1], 0);

  const double exec_seconds = exec_timer.Seconds();
  return make_tuple(emu->MachineChecksum(), emu->ImageChecksum(),
                    InstHistoChecksum(emu),
                    MemTrace(emu),
                    exec_seconds);
}

int main(int argc, char **argv) {
  AnsiInit();
  char date[256] = {};
  int64 now = time(nullptr);
  strftime(date, 254, "%d-%b-%Y %H:%M:%S", localtime(&now));

  string desc = "unknown";
  if (argc > 1) {
    desc = argv[1];
  }
  printf("Starting benchmark at [%s] called [%s]\n",
         date, desc.c_str());

  Timer warm_timer;

  Timer startup_timer;
  // TODO: This is not really fair since it counts all the IO.
  std::unique_ptr<Emulator> emu(Emulator::Create(ROMFILE));
  CHECK(emu.get() != nullptr);
  const double startup_seconds = startup_timer.Seconds();
  vector<uint8> start = emu->SaveUncompressed();

  vector<uint8> movie = SimpleFM7::ParseString(MOVIE);
  printf("Movie: %d frames\n", (int)movie.size());
  CHECK(!movie.empty());

  printf("Benchmarking %s. Startup in %0.4f sec...\n",
         ROMFILE, startup_seconds);

  double exec_seconds = -1.0;

  // This is simplified to just use a time budget, since each
  // execution is so slow (and the number of seconds so high)
  // that convergence otherwise takes way too long.
  static constexpr int MIN_EXECUTIONS = 1;
  static constexpr int MAX_EXECUTIONS = 10;
  static constexpr int MAX_SECONDS = 60.0 * 30.0;
  static constexpr int MIN_SECONDS = 60.0;
  int executions = 0;
  double total_time = 0.0;
  vector<int> last_means;
  for (int i = 0; /* exit upon convergence */; i++) {
    const auto [ram, img, inst, mem, sec] =
      RunBenchmark(emu.get(), start, movie);

    CHECK(ram == expected_nes &&
          img == expected_img &&
          inst == expected_inst) <<
      StringPrintf("Failed checksums on round %d:\n"
                   "RAM:  Got %016llx wanted %016llx\n"
                   "IMG:  Got %016llx wanted %016llx\n"
                   "INST: Got %016llx wanted %016llx\n"
                   "MEM:  Got %016llx wanted %016llx\n",
                   i,
                   ram, expected_nes,
                   img, expected_img,
                   inst, expected_inst,
                   mem, expected_mem);

    executions++;
    total_time += sec;
    double mean = total_time / (double)executions;

    if ((executions >= MIN_EXECUTIONS &&
         total_time >= MIN_SECONDS) &&
        (executions >= MAX_EXECUTIONS ||
         total_time > MAX_SECONDS)) {
      // Done!
      exec_seconds = mean;
      break;
    }

    printf("Round %4d in %.4f, mean %.4f, %lld cheats\n",
           executions, sec, mean, hfluint8::NumCheats());
    fflush(stdout);
  }

  const uint64 nes_checksum = emu->MachineChecksum();
  const uint64 img_checksum = emu->ImageChecksum();
  fprintf(stderr,
          "NES checksum: %016llx\n"
          "Img checksum: %016llx\n",
          nes_checksum,
          img_checksum);

  fprintf(stderr,
          "Startup time: %.4fs\n"
          "Exec time:    %.4fs\n",
          startup_seconds, exec_seconds);

  const int64 cheats = hfluint8::NumCheats();
  const double cpr = (double)cheats / (double)executions;
  fprintf(stderr,
          "Cheats:       %lld\n"
          "Cheats/round: %.4f\n",
          cheats, cpr);

  int status = 0;
  if (nes_checksum != expected_nes) {
    fprintf(stderr, "*** NES checksum mismatch. "
            "Got %016llx; wanted %016llx!\n",
            nes_checksum, expected_nes);
    status = -1;
  }
  if (img_checksum != expected_img) {
    fprintf(stderr, "*** Img checksum mismatch. "
            "Got %016llx; wanted %016llx!\n",
            img_checksum, expected_img);
    status = -1;
  }

  {
    FILE *f = fopen("bench.csv", "a");
    CHECK(f);
    fprintf(f, "%s,%s,%.4f,%.4f,%.4f,%lld,%.4f,%s\n",
            date,
            status ? "FAIL" : "OK",
            // was warm_time
            0.0,
            startup_seconds, exec_seconds,
            cheats, cpr, desc.c_str());
    fclose(f);
  }

  {
    string inst;
    int64 total = 0;
    const X6502 *X = emu->GetFC()->X;
    for (int i = 0; i < 256; i++) {
      int64 count = X->inst_histo[i];
      total += count;
      StringAppendF(&inst, "%02x:\t%s\t\t%lld\n",
                    i, fceulib_2a03_opcode_name[i],
                    count);
    }
    Util::WriteFile("instruction-counts-bench.txt", inst);
    printf("Total instructions: %lld\n", total);
    printf("Wrote histogram to instruction-counts-bench.txt\n");
  }

  return status;
}
