
#ifndef _LOWERCASE_LOADFONTS_H
#define _LOWERCASE_LOADFONTS_H

#include <functional>
#include <string>
#include <shared_mutex>
#include <thread>
#include <cstdint>
#include <memory>

#include "fontdb.h"
#include "fonts/ttf.h"
#include "image.h"
#include "font-problem.h"

// Load fonts in the background.
struct VectorLoadFonts {
  using int64 = int64_t;

  // Starts loading fonts into memory immediately, in separate
  // threads. Note that we reject the entire font if any letter
  // exceeds row_max_points, because we don't want to create bias
  // in the training data (e.g. to make Q less common).
  VectorLoadFonts(
      // If this returns true (e.g. because the startup process is
      // aborted), just stop loading.
      std::function<bool()> ExitEarly,
      const std::vector<int> &row_max_points,
      int max_parallelism,
      int64 max_fonts);

  ~VectorLoadFonts();

  // Call once. Waits for font vector to be complete. After this
  // returns, safe to access the vector (from a single thread) without
  // taking the mutex.
  void Sync();

  std::shared_mutex fonts_m;
  // Protected by fonts_m.
  // The font pointers are owned by this object.
  std::vector<TTF *> fonts;

private:
  void Init();

  const int max_parallelism;
  const int64 max_fonts;
  const std::function<bool()> ExitEarly;
  const std::vector<int> row_max_points;

  std::unique_ptr<FontDB> font_db;
  std::unique_ptr<std::thread> init_thread;
};

struct SDFLoadFonts {
  using int64 = int64_t;

  // The SDFs are somewhat expensive to compute, so they are generated up front
  // and kept in memory according to the SDFConfig.
  // At 64x64, 10,000 fonts take about 2 GB.
  SDFLoadFonts(
      std::function<bool()> ExitEarly,
      FontProblem::SDFConfig config,
      int max_parallelism,
      int64 max_fonts);

  ~SDFLoadFonts();

  // As above.
  void Sync();

  struct Font {
    TTF *ttf = nullptr;
    // 52 entries. a-z then A-Z.
    std::vector<ImageA> sdfs;
  };
  std::shared_mutex fonts_m;
  // Protected by fonts_m.
  // The font pointers are owned by this object.
  std::vector<Font> fonts;

  // Sync() must have been called or lock held.
  int64 NumFailed() const {
    return num_failed;
  }

  void Init();

private:
  const int max_parallelism;
  const int64 max_fonts;
  const std::function<bool()> ExitEarly;
  const FontProblem::SDFConfig config;

  int64 num_failed = 0;

  std::unique_ptr<FontDB> font_db;
  std::unique_ptr<std::thread> init_thread;
};


#endif
