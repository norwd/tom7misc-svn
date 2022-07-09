
#include <optional>
#include <array>

#include "base/logging.h"
#include "base/stringprintf.h"
#include "image.h"
#include "expression.h"
#include "half.h"
#include "hashing.h"

#include "choppy.h"
#include "grad-util.h"
#include "color-util.h"
#include "arcfour.h"

// Makes a database of "choppy" functions.

using DB = Choppy::DB;
using Allocator = Exp::Allocator;
using Table = Exp::Table;

static const Exp *GetOp5(Exp::Allocator *alloc,
                         int iters1, int iters2,
                         int off1, int off2,
                         double prescale, double postscale,
                         bool do_postscale = true) {
  const Exp *fa =
    alloc->TimesC(
        // x - 4
        alloc->PlusC(alloc->TimesC(alloc->Var(),
                                   Exp::GetU16((half)prescale)), off1),
        // * 0.999 ...
        0x3bffu,
        iters1);

  const Exp *fb =
    alloc->TimesC(
        // x - 4
        alloc->PlusC(alloc->TimesC(alloc->Var(),
                                   Exp::GetU16((half)prescale)), off2),
        // * 0.999 ...
        0x3bffu,
        iters2);

  const Exp *f0 =
    alloc->PlusE(fa,
                 alloc->Neg(fb));

  return do_postscale ?
    alloc->TimesC(f0, Exp::GetU16((half)postscale)) : f0;
}


static void Explore5(DB *db) {
  // Original int bounds: {49758, 49152};
  const std::array<double, 3> BASE_DOUBLES = {0.0039, -3.9544};

  ArcFour rc(StringPrintf("explore.%lld", time(nullptr)));
  RandomGaussian gauss(&rc);

  Exp::Allocator *alloc = &db->alloc;

  // Can return nullptr if not feasible.
  auto MakeExp = [&](int iters1, int iters2,
                     int a, int b,
                     double x, double y,
                     bool do_postscale) -> const Exp * {
      // Derive the scale.
      const Exp *exp0 = GetOp5(alloc, iters1, iters2, a, b, x, y,
                               do_postscale);

      const half mid_x = (half)(1.0/(Choppy::GRID * 2.0));

      // First put x=0 (really, the midpoint of zero) at y=0.
      uint16 offset = Exp::EvaluateOn(exp0, mid_x);

      // flip sign bit so that we are subtracting the computed
      // offset.
      const Exp *exp1 = alloc->PlusC(exp0, 0x8000 ^ offset);

      // The result should be exactly zero (or negative zero).
      const uint16 new_offset = Exp::EvaluateOn(exp1, mid_x) & 0x7FFF;
      if (new_offset != 0x000) {
        printf("After offsetting, unexpectedly got %04x??\n",
               new_offset);
        return nullptr;
      }

      // XXX this loop does not necessarily find a scale for a function
      // that can in fact be scaled!
      // Now find any nonzero point.

      // If zero, we will fail.
      half multiplier = (half)0.0f;

      for (int i = 0; i < Choppy::GRID; i++) {
        half x = (half)((i / (double)(Choppy::GRID/2)) - 1.0);
        x += (half)(1.0/(Choppy::GRID * 2.0));

        // Ignore the sign so that scaling preserves it.
        half y = fabs(Exp::GetHalf(Exp::EvaluateOn(exp1, Exp::GetU16(x))));
        if (y != (half)0.0) {
          // ... and just scale it to be 1.0.
          half yinv = (half)(1.0 / (Choppy::GRID/2.0)) / y;

          /*
          const Exp *exp2 = alloc->TimesC(exp1, Exp::GetU16(yinv));
          const half newy =
            Exp::GetHalf(Exp::EvaluateOn(exp2, Exp::GetU16(x)));
            printf("Scale y = %.5f by %.5f, expect %.5f got %.5f\n",
            (float)y, (float)yinv, (float)(y * yinv), (float)newy);
          */

          if (yinv > multiplier) multiplier = yinv;
        }
      }

      if (multiplier > (half)0.0) {
        return alloc->TimesC(exp1, Exp::GetU16(multiplier));
      }

      // Was all zero; this is common but not feasible!
      return nullptr;
    };

  int64 infeasible = 0;
  int64 added_new = 0, added_smaller = 0;
  int64 not_choppy = 0, outside_grid = 0, not_new = 0;
  int64 done = 0;

  // run a full grid. bc00 = -1, c600 = -6.
  const int64 SIZE = 0xc600 - 0xbc00;
  const int64 TOTAL = SIZE * SIZE;

  // Each loop takes about 1m22s.
  for (int loop = 0; loop < 10; loop++)
  for (int i1 = 0xbc00; i1 < 0xc600; i1++) {
    for (int i2 = 0xbc00; i2 < 0xc600; i2++) {
  // for (int i1 : {49758}) {
  // for (int i2 : {49152}) {
      done++;
      if (done % 25000 == 0) {
        fprintf(stderr,
                "Ran %lldk/%dk (%.1f%%). %lld new %lld opt, "
                "%lld infeasible (%.1f%%)\n",
                done / 1024, TOTAL / 1024,
                (done * 100.0) / TOTAL,
                added_new, added_smaller,
                infeasible, (infeasible * 100.0) / done);
      }

      auto [d1, d2, d3_] = BASE_DOUBLES;

      double r1 = gauss.Next() * 0.5 - 0.25;
      double r2 = gauss.Next() * 0.5 - 0.25;

      int io1 = gauss.Next() * 100 - 50;
      int io2 = gauss.Next() * 100 - 50;

      int iters1 = std::clamp(200 + io1, 10, 1000);
      int iters2 = std::clamp(300 + io2, 10, 1000);

      bool do_postscale = RandFloat(&rc) < 0.9;

      const Exp *exp = MakeExp(
          iters1, iters2,
          i1 + io1, i2 + io2,
          d1 + r1,
          d2 + r2,
          do_postscale);
      if (!exp) {
        infeasible++;
        continue;
      }

      #if 0
      ImageRGBA img(1920, 1920);
      img.Clear32(0x000000FF);
      GradUtil::Grid(&img);
      Table table = Exp::TabulateExpression(exp);
      GradUtil::Graph(table, 0xFFFFFF88, &img);
      img.Save("findchop.png");
      #endif

      DB::AddResult ar = db->Add(exp);
      switch (ar) {
      case DB::AddResult::SUCCESS_NEW:
        added_new++;
        break;
      case DB::AddResult::SUCCESS_SMALLER:
        added_smaller++;
        break;
      case DB::AddResult::NOT_CHOPPY:
        not_choppy++;
        break;
      case DB::AddResult::OUTSIDE_GRID:
        outside_grid++;
        break;
      case DB::AddResult::NOT_NEW:
        not_new++;
        break;
      }
    }
  }

  printf("Total done: %lld\n"
         "Infeasible: %lld\n"
         "Added new: %lld\n"
         "Added smaller: %lld\n"
         "Not choppy: %lld\n"
         "Outside grid: %lld\n"
         "Not new: %lld\n",
         done, infeasible,
         added_new, added_smaller,
         not_choppy, outside_grid, not_new);
}

static void Expand(DB *db) {
  // Copy so that we don't have to deal with iterator
  // invalidation.
  std::vector<const Exp *> all;
  for (const auto &[k, v] : db->fns) all.push_back(v);

  const Exp *var = db->alloc.Var();

  int64 added_new = 0, added_smaller = 0;
  int64 not_choppy = 0, outside_grid = 0, not_new = 0;
  int64 done = 0;

  for (const Exp *e : all) {
    for (int xo = -8; xo <= 8; xo++) {
      if (xo != 0) {
        done++;
        half dx = (half)(xo/(Choppy::GRID * 2.0));
        const Exp *shifted =
          Exp::Subst(&db->alloc,
                     e,
                     db->alloc.PlusC(var, Exp::GetU16(dx)));

        DB::AddResult ar = db->Add(shifted);
        switch (ar) {
        case DB::AddResult::SUCCESS_NEW:
          added_new++;
          break;
        case DB::AddResult::SUCCESS_SMALLER:
          added_smaller++;
          break;
        case DB::AddResult::NOT_CHOPPY:
          not_choppy++;
          break;
        case DB::AddResult::OUTSIDE_GRID:
          outside_grid++;
          break;
        case DB::AddResult::NOT_NEW:
          not_new++;
          break;
        }

      }

    }
  }

  printf("Total done: %lld\n"
         "Added new: %lld\n"
         "Added smaller: %lld\n"
         "Not choppy: %lld\n"
         "Outside grid: %lld\n"
         "Not new: %lld\n",
         done,
         added_new, added_smaller,
         not_choppy, outside_grid, not_new);
}

#if 0
[[maybe_unused]]
static void SeedDBFromHeader(DB *db) {
  Allocator *alloc = &db->alloc;
  #include "seed-db.h"
  for (const char *s : FNS) {
    string error;
    const Exp *e = Exp::Deserialize(alloc, s, &error);
    CHECK(e != nullptr) << error << "\n" << s;
    CHECK(db->Add(e));
  }
}
#endif

static void LoadDB(DB *db) {
  std::vector<string> lines =
    Util::ReadFileToLines("chopdb.txt");
  int64 rejected = 0;
  Allocator *alloc = &db->alloc;
  for (const string &line : lines) {
    if (line.empty()) continue;
    if (line[0] == '/') continue;
    string error;
    const Exp *e = Exp::Deserialize(alloc, line, &error);
    CHECK(e != nullptr) << error << "\n" << line;
    // Can be rejected by new stricter criteria (e.g. outside_grid)
    // but nothing in the database should be invalid.
    auto res = db->Add(e);
    if (res == DB::AddResult::NOT_CHOPPY)
      rejected++;
  }
  if (rejected) {
    printf("Warning: Rejected %lld as not choppy.\n", rejected);
  }
  printf("Chop DB size: %d\n", (int)db->fns.size());
}


int main(int argc, char **argv) {
  DB db;
  LoadDB(&db);

  Explore5(&db);
  // Expand(&db);

  fprintf(stderr, "\n\ndb now has %d distinct fns\n", (int)db.fns.size());

  Util::WriteFile("chopdb.txt", db.Dump());
  return 0;
}
