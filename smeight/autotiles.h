
#ifndef _SMEIGHT_AUTOTILES_H
#define _SMEIGHT_AUTOTILES_H

#include "smeight.h"

#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>

#include "../fceulib/emulator.h"
#include "autocamera.h"

struct AutoTiles {
  static constexpr int NUM_EMULATORS = 16;

  // Like autocamera, we keep a bunch of emulators around so that we
  // can do stuff in parallel.
  explicit AutoTiles(const std::string &game) {
    printf("Creating %d emulators for AutoTiles...\n", NUM_EMULATORS);
    for (int i = 0; i < NUM_EMULATORS; i++) {
      emus.push_back(Emulator::Create(game));
    }
  }

  ~AutoTiles() {
    for (Emulator *emu : emus) delete emu;
  }

  // Here we try to compute a simple mapping from tile id to bool,
  // where true indicates solid.
  static uint64_t TilesCRC(const Emulator *emu);

  struct Tileset {
    std::unordered_map<uint32_t, bool> is_solid;
  };

  enum Solidity { SOLID, OPEN, UNKNOWN, };

  struct Tile {
    uint32_t fourtiles = 0;
    Solidity solidity = UNKNOWN;
  };

  // Get the solidity of the 16x16 tiles in the emulator state.
  // Supports scroll position, but treats it coarsely (16x16 pixel
  // blocks). If there are unknown tiles, pauses to experiment on them
  // (but may not succeed!).
  std::vector<Tile> GetTileInfo(
      Emulator *emu,
      // Need to know how to make the player face
      // left and right.
      const AngleRule &left, const AngleRule &right,
      bool is_top,
      const std::vector<AutoCamera::XYSprite> &cams);

  // All of the tilesets we've seen, keyed by CRC.
  std::unordered_map<uint64_t, Tileset*> tilesets;

 private:
  std::vector<Emulator *> emus;
};

#endif
