
#ifndef _LOWERCASE_LOADFONTS_H
#define _LOWERCASE_LOADFONTS_H

#include <functional>
#include <string>
#include <shared_mutex>
#include <thread>
#include <cstdint>
#include <memory>

#include "fontdb.h"
#include "ttf.h"

// Load fonts in the background.
struct LoadFonts {
  using int64 = int64_t;

  // Starts loading fonts into memory immediately, in separate
  // threads. Note that we reject the entire font if any letter
  // exceeds row_max_points, because we don't want to create bias
  // in the training data (e.g. to make Q less common).
  LoadFonts(
      // If this returns true (e.g. because the startup process is
      // aborted), just stop loading.
      std::function<bool()> ExitEarly,
      const vector<int> &row_max_points,
      int max_parallelism,
      int64 max_fonts);

  ~LoadFonts();
  
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


#endif
