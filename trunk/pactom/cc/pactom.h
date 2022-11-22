
#ifndef _PACTOM_H
#define _PACTOM_H

#include <memory>
#include <string>
#include <vector>
#include <optional>
#include <utility>
#include <map>

#include "geom/latlon.h"
#include "bounds.h"

struct PacTom {

  static std::unique_ptr<PacTom> FromFiles(
      const std::vector<std::string> &files,
      const std::optional<std::string> &neighborhoods);

  // Parallel to above.
  struct Run {
    std::string name;
    // Might not have a date :/
    int year = 0, month = 0, day = 0;
    // Position and elevation, but no timing info.
    std::vector<std::pair<LatLon, double>> path;
  };
  std::vector<Run> runs;

  // Borders of neighborhoods, if loaded.
  std::map<std::string, std::vector<LatLon>> hoods;

  std::vector<std::string> neighborhood_names;

  // Return -1 if not in any neighborhood; otherwise neighborhood id.
  int InNeighborhood(LatLon pos) const;

 private:
  PacTom();

  std::vector<std::pair<Bounds, const std::vector<LatLon> *>> hood_boxes;
};


#endif
