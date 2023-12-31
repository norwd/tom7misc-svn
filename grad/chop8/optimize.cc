
#include <optional>
#include <array>
#include <utility>

#include "base/logging.h"
#include "base/stringprintf.h"
#include "image.h"
#include "expression.h"
#include "half.h"
#include "hashing.h"

#include "choppy8.h"
#include "grad-util.h"
#include "color-util.h"
#include "arcfour.h"
#include "ansi.h"
#include "threadutil.h"
#include "timer.h"
#include "periodically.h"

#include "state.h"

using DB = Choppy::DB;
using Allocator = Exp::Allocator;
using Table = Exp::Table;

constexpr int MAX_THREADS = 4;

static void MaybeSaveDB(DB *db) {
  static std::mutex *m = new std::mutex;
  static Periodically *save_per = new Periodically(60.0);

  {
    std::unique_lock<std::mutex> ml(*m);
    if (save_per->ShouldRun()) {
      Util::WriteFile("optimized-checkpoint.txt", db->Dump());
      printf("Wrote optimized checkpoint.\n");
    }
  }
}


// Reduce iterations multiplicatively.
// Sometimes there are way too many iterations, and it would
// be faster to use binary search to find the minimal number.
// This is hard to arrange with the "depth" concept. So this
// is basically just like the DropNode pass, but it reduces
// the number of iterations by a multiplicative factor,
// rather than one at a time.
struct ReduceIters {
  const std::string Name() const { return "ReduceIters"; }

  const Exp *Run(Allocator *alloc_, const Exp *e) {
    alloc = alloc_;
    did_drop = false;
    depth = 0;
    // idea is to keep some secondary part of the cursor,
    // like the lower/upper bound?
    return Rec(e);
  }

  bool Done() const { return !did_drop; }

  string CurrentDepth() const {
    return StringPrintf("%lld", stop_depth);
  }
  void NextDepth() { stop_depth++; }

  const Exp *Rec(const Exp *e) {
    if (depth == stop_depth) {
      // Drop the current node.
      did_drop = true;
      switch (e->type) {
      case VAR:
        return e;
      case PLUS_C:
        if (e->iters > 1) {
          int new_iters = std::max(1, (int)std::round(e->iters * 0.90));
          return alloc->PlusC(e->a, e->c, new_iters);
        } else {
          return e;
        }
        break;
      case TIMES_C:
        CHECK(e->iters != 0);
        if (e->iters > 1) {
          int new_iters = std::max(1, (int)std::round(e->iters * 0.90));
          return alloc->TimesC(e->a, e->c, new_iters);
        } else {
          return e;
        }
        break;
      case PLUS_E:
        return e;
      default:
        CHECK(false);
        return nullptr;
      }

    } else {
      // Go deeper.
      depth++;
      switch (e->type) {
      case VAR:
        return e;
      case PLUS_C:
        return alloc->PlusC(Rec(e->a), e->c, e->iters);
      case TIMES_C:
        return alloc->TimesC(Rec(e->a), e->c, e->iters);
      case PLUS_E: {
        const Exp *ea = Rec(e->a);
        const Exp *eb = Rec(e->b);
        return alloc->PlusE(ea, eb);
      }
      default:
        CHECK(false);
        return nullptr;
      }
    }
  }

private:
  Allocator *alloc = nullptr;
  int64 depth = 0;
  int64 stop_depth = 0;
  bool did_drop = false;
};


struct DropNode {
  const std::string Name() const { return "DropNode"; }

  const Exp *Run(Allocator *alloc_, const Exp *e) {
    alloc = alloc_;
    did_drop = false;
    depth = 0;
    return Rec(e);
  }

  bool Done() const { return !did_drop; }

  string CurrentDepth() const {
    return StringPrintf("%lld", stop_depth);
  }
  void NextDepth() { stop_depth++; }

  const Exp *Rec(const Exp *e) {
    if (depth == stop_depth) {
      // Drop the current node.
      did_drop = true;
      switch (e->type) {
      case VAR:
        return e;
      case PLUS_C:
        CHECK(e->iters != 0);
        if (e->iters == 1) {
          return e->a;
        } else {
          return alloc->PlusC(e->a, e->c, e->iters - 1);
        }
        break;
      case TIMES_C:
        CHECK(e->iters != 0);
        if (e->iters == 1) {
          return e->a;
        } else {
          return alloc->TimesC(e->a, e->c, e->iters - 1);
        }
        break;
      case PLUS_E:
        // We don't expect it to be useful to drop a whole
        // arm (nor do we use this in the expressions that need
        // most to be optimized), so we don't worry about trying
        // both a and b here.
        return e->a;
      default:
        CHECK(false);
        return nullptr;
      }

    } else {
      // Go deeper.
      depth++;
      switch (e->type) {
      case VAR:
        return e;
      case PLUS_C:
        return alloc->PlusC(Rec(e->a), e->c, e->iters);
      case TIMES_C:
        return alloc->TimesC(Rec(e->a), e->c, e->iters);
      case PLUS_E: {
        const Exp *ea = Rec(e->a);
        const Exp *eb = Rec(e->b);
        return alloc->PlusE(ea, eb);
      }
      default:
        CHECK(false);
        return nullptr;
      }
    }
  }

private:
  Allocator *alloc = nullptr;
  int64 depth = 0;
  int64 stop_depth = 0;
  bool did_drop = false;
};

static int64 TotalDepth(const Exp *e) {
  int64 depth = 0;
  std::function<void(const Exp *)> Rec = [&Rec, &depth](const Exp *e) {
      depth++;
      switch (e->type) {
      case VAR:
        return;
      case PLUS_C:
        Rec(e->a);
        return;
      case TIMES_C:
        Rec(e->a);
        return;
      case PLUS_E: {
        Rec(e->a);
        Rec(e->b);
        return;
      }
      default:
        CHECK(false);
        return;
      }
    };

  Rec(e);
  return depth;
}

static const Exp *CleanRec(Allocator *alloc, const Exp *exp) {
  switch (exp->type) {
  case VAR:
    return exp;
  case PLUS_C: {
    // Adding zero or negative zero (usually) does nothing.
    if (exp->c == 0x0000 || exp->c == 0x8000)
      return CleanRec(alloc, exp->a);

    const Exp *ea = CleanRec(alloc, exp->a);
    if (ea->type == PLUS_C &&
        ea->c == exp->c) {
      int32 new_iters = (int32)exp->iters + (int32)ea->iters;
      if (new_iters <= 65535) {
        return alloc->PlusC(ea->a, exp->c, new_iters);
      }
    }

    return alloc->PlusC(ea, exp->c, exp->iters);
  }
  case TIMES_C: {
    // Multiplying by one does nothing.
    if (exp->c == 0x3c00)
      return CleanRec(alloc, exp->a);

    // TODO: squash multiply by -1 with another multiplication.

    const Exp *ea = CleanRec(alloc, exp->a);
    if (ea->type == TIMES_C &&
        ea->c == exp->c) {
      int32 new_iters = (int32)exp->iters + (int32)ea->iters;
      if (new_iters <= 65535) {
        return alloc->TimesC(ea->a, exp->c, new_iters);
      }
    }

    return alloc->TimesC(ea, exp->c, exp->iters);
  }
  case PLUS_E:
    return alloc->PlusE(CleanRec(alloc, exp->a),
                        CleanRec(alloc, exp->b));
  default:
    CHECK(false);
    return nullptr;
  }
}

static void OptimizeOne(DB *db,
                        const DB::key_type &key,
                        const Exp *exp) {
  Allocator *alloc = &db->alloc;

  static constexpr double SAVE_EVERY = 60.0 * 5.0;

  auto StillWorks = [&](const Exp *exp) {
      auto chopo = Choppy::GetChoppy(exp);
      if (!chopo.has_value()) return false;
      if (chopo.value() != key) return false;

      return true;
    };

  auto StillWorksLinear = [&](const vector<Step> &steps) {
      Allocator alloc;
      const Exp *exp = State::GetExpressionFromSteps(&alloc, steps);
      return StillWorks(exp);
    };

  const int start_size = Exp::ExpSize(exp);

  // Make sure the entry has the key it purports to.
  CHECK(StillWorks(exp));

  // First, delete expressions that do nothing.
  const Exp *clean_exp = CleanRec(alloc, exp);
  if (StillWorks(exp)) {
    exp = clean_exp;
  } else {
    CPrintf(ANSI_RED "Clean did not preserve behavior!"
            ANSI_RESET "\n");
  }


  if (State::CanBeLinearized(exp)) {
    CPrintf(ANSI_GREEN "Linearizable." ANSI_RESET "\n");
    vector<Step> steps = State::Linearize(exp);


    {
      CPrintf("Running ChopSection on expression of size " ANSI_YELLOW
              "%d" ANSI_RESET "\n", (int)steps.size());

      int64 tries = 0;
      Timer loop_timer;
      // average one per second across all threads.
      Periodically status_per((double)MAX_THREADS);
      // Try to skip the first report, though.
      (void)status_per.ShouldRun();
      // Save occasionally so that we don't lose too much
      // progress if we stop early.
      Periodically checkpoint_per(SAVE_EVERY);
      (void)checkpoint_per.ShouldRun();

      constexpr int MAX_CHOP = 16;

      for (int start_idx = 0; start_idx < steps.size(); start_idx++) {
        // PERF: At the end, this tries (harmlessly) chopping non-existent
        // steps.
        for (int chop_size = MAX_CHOP; chop_size > 0; chop_size--) {
          if (status_per.ShouldRun()) {
            int64 step_size = 0;
            for (const Step &step : steps)
              step_size += step.iters;

            CPrintf(ANSI_GREY "[%s] " ANSI_RESET
                    "Tries " ANSI_YELLOW "%lld" ANSI_RESET " ("
                    ANSI_CYAN "%.3f " ANSI_RESET "/s) size "
                    ANSI_PURPLE "%lld" ANSI_RESET
                    " depth " ANSI_RED "%d" ANSI_RESET "/"
                    ANSI_YELLOW "%d" ANSI_RESET "\n",
                    "ChopSection",
                    tries, (tries / loop_timer.Seconds()), step_size,
                    start_idx, (int)steps.size());
          }

          // Try chopping.
          vector<Step> chopped;
          chopped.reserve(steps.size());
          for (int i = 0; i < steps.size(); i++) {
            if (i >= start_idx && i < start_idx + chop_size) {
              // skip it.
            } else {
              chopped.push_back(steps[i]);
            }
          }

          tries++;
          if (StillWorksLinear(chopped)) {
            steps = std::move(chopped);
            CPrintf(ANSI_GREY "[%s] " ANSI_RESET
                    "Chopped " ANSI_BLUE "%d" ANSI_RESET " steps at "
                    ANSI_PURPLE "%d" ANSI_RESET "! Now "
                    ANSI_YELLOW "%d" ANSI_RESET " steps.\n",
                    "ChopSection", chop_size, start_idx,
                    (int)steps.size());

            exp = State::GetExpressionFromSteps(alloc, steps);
            if (checkpoint_per.ShouldRun()) {
              db->Add(exp);
              MaybeSaveDB(db);
            }

            // Reset chop size so that we keep trying to chop at
            // this position (the steps have been replaced).
            chop_size = MAX_CHOP + 1;
          }
        }
      }
    }

    {
      CPrintf("Running Meld on expression of size " ANSI_YELLOW
              "%d" ANSI_RESET "\n", (int)steps.size());

      int64 tries = 0;
      Timer loop_timer;
      // average one per second across all threads.
      Periodically status_per((double)MAX_THREADS);
      // Try to skip the first report, though.
      (void)status_per.ShouldRun();
      // Save occasionally so that we don't lose too much
      // progress if we stop early.
      Periodically checkpoint_per(SAVE_EVERY);
      (void)checkpoint_per.ShouldRun();

      constexpr int MAX_CHOP = 16;

      for (int start_idx = 0; start_idx < steps.size(); start_idx++) {
        // PERF: At the end, this tries (harmlessly) melding non-existent
        // steps.
        for (int meld_size = MAX_CHOP; meld_size > 1; meld_size--) {
          if (status_per.ShouldRun()) {
            int64 step_size = 0;
            for (const Step &step : steps)
              step_size += step.iters;

            CPrintf(ANSI_GREY "[%s] " ANSI_RESET
                    "Tries " ANSI_YELLOW "%lld" ANSI_RESET " ("
                    ANSI_CYAN "%.3f " ANSI_RESET "/s) size "
                    ANSI_PURPLE "%lld" ANSI_RESET
                    " depth " ANSI_RED "%d" ANSI_RESET "/"
                    ANSI_YELLOW "%d" ANSI_RESET "\n",
                    "Meld",
                    tries, (tries / loop_timer.Seconds()), step_size,
                    start_idx, (int)steps.size());
          }

          // Try chopping.
          vector<Step> chopped;
          chopped.reserve(steps.size());
          double sum = 0.0;
          double product = 1.0;
          for (int i = 0; i < steps.size(); i++) {
            if (i >= start_idx && i < start_idx + meld_size) {
              const Step &step = steps[i];
              double c = (double)Exp::GetHalf(step.c);

              // This kind of doesn't make sense if there are
              // a mix of plus/times, but we try anyway.
              if (step.type == STEP_PLUS) {
                // Or just multiply?
                for (int z = 0; z < step.iters; z++)
                  sum += c;
              } else {
                CHECK(step.type == STEP_TIMES);
                for (int z = 0; z < step.iters; z++)
                  product *= c;
              }

            } else  {
              if (i == start_idx + meld_size) {
                if (sum != 0.0)
                  chopped.push_back(
                      Step(STEP_PLUS, Exp::GetU16((half)sum), 1));
                if (product != 1.0)
                  chopped.push_back(
                      Step(STEP_TIMES, Exp::GetU16((half)product), 1));
              }
              chopped.push_back(steps[i]);
            }
          }

          tries++;
          if (chopped.size() < steps.size() &&
              StillWorksLinear(chopped)) {
            steps = std::move(chopped);
            CPrintf(ANSI_GREY "[%s] " ANSI_RESET
                    "Melded " ANSI_BLUE "%d" ANSI_RESET " steps at "
                    ANSI_PURPLE "%d" ANSI_RESET "! Now "
                    ANSI_YELLOW "%d" ANSI_RESET " steps.\n",
                    "Meld", meld_size, start_idx,
                    (int)steps.size());

            exp = State::GetExpressionFromSteps(alloc, steps);
            if (checkpoint_per.ShouldRun()) {
              db->Add(exp);
              MaybeSaveDB(db);
            }

            // Reset chop size so that we keep trying to chop at
            // this position (the steps have been replaced).
            meld_size = MAX_CHOP + 1;
          }
        }
      }
    }


    {
      CPrintf("Running BinaryIters on expression of size " ANSI_YELLOW
              "%d" ANSI_RESET "\n", (int)steps.size());

      int64 tries = 0;
      Timer loop_timer;
      // average one per second across all threads.
      Periodically status_per((double)MAX_THREADS);
      // Try to skip the first report, though.
      (void)status_per.ShouldRun();
      // Save occasionally so that we don't lose too much
      // progress if we stop early.
      Periodically checkpoint_per(SAVE_EVERY);
      (void)checkpoint_per.ShouldRun();

      for (int start_idx = 0; start_idx < steps.size(); start_idx++) {
        if (status_per.ShouldRun()) {
          int64 step_size = 0;
          for (const Step &step : steps)
            step_size += step.iters;

          CPrintf(ANSI_GREY "[%s] " ANSI_RESET
                  "Tries " ANSI_YELLOW "%lld" ANSI_RESET " ("
                  ANSI_CYAN "%.3f " ANSI_RESET "/s) size "
                  ANSI_PURPLE "%lld" ANSI_RESET
                  " depth " ANSI_RED "%d" ANSI_RESET "/"
                  ANSI_YELLOW "%d" ANSI_RESET "\n",
                  "BinaryIters",
                  tries, (tries / loop_timer.Seconds()), step_size,
                  start_idx, (int)steps.size());
        }

        // Binary search for minimal iters.
        {
          CHECK(start_idx >= 0 && start_idx < steps.size());
          const int original_iters = steps[start_idx].iters;
          int search_steps = 0;

          // iters does not work for values < lower_bound,
          // so minimal iters is >= this value.
          int lower_bound = 1;
          // iters does work at this value, so minimal iters is
          // <= this value.
          int upper_bound = steps[start_idx].iters;

          // Note that in the common case that iters == 1, we
          // are already done.
          while (lower_bound != upper_bound) {
            CHECK(lower_bound < upper_bound);
            // Rounding down, since we already know the result for upper_bound.
            const int next_iters =
              lower_bound + ((upper_bound - lower_bound) >> 1);

            search_steps++;
            steps[start_idx].iters = next_iters;
            tries++;
            if (StillWorksLinear(steps)) {
              CHECK(next_iters < upper_bound);
              upper_bound = next_iters;
              CPrintf(ANSI_GREY "[%s] " ANSI_RESET
                      "lb " ANSI_BLUE "%d" ANSI_RESET " ub "
                      ANSI_PURPLE "%d" ANSI_RESET " ("
                      ANSI_YELLOW "%d" ANSI_RESET " steps).\n",
                      "BinaryIters",
                      lower_bound, upper_bound, search_steps);
            } else {
              lower_bound = next_iters + 1;
            }
          }

          steps[start_idx].iters = upper_bound;
          if (upper_bound < original_iters) {

            CPrintf(ANSI_GREY "[%s] " ANSI_RESET
                    ANSI_GREEN "Reduced!" ANSI_RESET
                    " Iters " ANSI_BLUE "%d" ANSI_RESET " -> "
                    ANSI_PURPLE "%d" ANSI_RESET
                    " in " ANSI_CYAN "%d" ANSI_RESET " steps.\n",
                    "BinaryIters",
                    original_iters, upper_bound,
                    search_steps);
            exp = State::GetExpressionFromSteps(alloc, steps);
            if (checkpoint_per.ShouldRun()) {
              db->Add(exp);
              MaybeSaveDB(db);
            }
          }
        }
      }
    }
  }

  // For each expression, see if we can remove it (or reduce
  // its iterations) but retain the property.

  auto DoPhase = [db, alloc, &exp, start_size, &StillWorks]<typename Phase>(
      Phase phase) {
    const int total_depth = TotalDepth(exp);
    int64 tries = 0;
    Timer loop_timer;
    // average one per second across all threads.
    Periodically status_per((double)MAX_THREADS);
    // Try to skip the first report, though.
    (void)status_per.ShouldRun();
    // Save occasionally so that we don't lose too much
    // progress if we stop early.
    Periodically checkpoint_per(SAVE_EVERY);
    (void)checkpoint_per.ShouldRun();

    int exp_size = Exp::ExpSize(exp);
    for (;;) {
      Allocator local_alloc;
      tries++;
      if (status_per.ShouldRun()) {
        CPrintf(ANSI_GREY "[%s] " ANSI_RESET
                "Tries " ANSI_YELLOW "%lld" ANSI_RESET " ("
                ANSI_CYAN "%.3f " ANSI_RESET "/s) size "
                ANSI_PURPLE "%d" ANSI_RESET
                " depth " ANSI_RED "%s" ANSI_RESET "/"
                ANSI_YELLOW "%d" ANSI_RESET "\n",
                phase.Name().c_str(),
                tries, (tries / loop_timer.Seconds()), exp_size,
                phase.CurrentDepth().c_str(), total_depth);
      }

      const Exp *e = phase.Run(&local_alloc, exp);

      // Nothing changed.
      if (phase.Done())
        break;

      // Does the trimmed expression still work?
      if (StillWorks(e)) {
        int new_size = Exp::ExpSize(e);
        if (new_size < exp_size) {
          CPrintf(ANSI_GREY "[%s] " ANSI_RESET
                  ANSI_GREEN "Reduced!" ANSI_RESET
                  " Start " ANSI_BLUE "%d" ANSI_RESET " now "
                  ANSI_PURPLE "%d" ANSI_RESET "\n",
                  phase.Name().c_str(),
                  start_size, new_size);
          exp = alloc->Copy(e);
          exp_size = new_size;
          if (checkpoint_per.ShouldRun()) {
            db->Add(exp);
            MaybeSaveDB(db);
          }
          // Keep the same stop depth, as we may be able to apply
          // another improvement at the same position.
        } else {
          // This can happen because we don't actually drop VAR nodes,
          // but we pretend we did (for uniformity).
          phase.NextDepth();
        }
      } else {
        // Eventually this gets larger than the expression and
        // we will be Done().
        phase.NextDepth();
      }
    }
  };

  if (!State::CanBeLinearized(exp)) {
    if (start_size > 1000) {
      CPrintf(ANSI_RED "Slow iter reduction: can't be linearized."
              ANSI_RESET "\n");
    }
    DoPhase(ReduceIters());
  }

  if (Exp::ExpSize(exp) < start_size) {
    db->Add(exp);
    MaybeSaveDB(db);
  }

  DoPhase(DropNode());

  const int end_size = Exp::ExpSize(exp);
  if (end_size < start_size) {
    CPrintf("Reduced from " ANSI_BLUE "%d" ANSI_RESET " to "
            ANSI_PURPLE "%d" ANSI_RESET "\n", start_size, end_size);
    db->Add(exp);
    MaybeSaveDB(db);
  } else {
    CPrintf(ANSI_GREY "No reduction (still %d)" ANSI_RESET "\n",
            start_size);
  }
}

int main(int argc, char **argv) {
  AnsiInit();
  CHECK(argc == 2) << "Give a database file on the command line.";

  DB db;
  printf("Load database:\n");
  db.LoadFile(argv[1]);

  std::vector<std::pair<const DB::key_type &,
                        const Exp *>> all;
  for (const auto &[k, v] : db.fns)
    all.emplace_back(k, v);

  ParallelApp(all,
              [&](const std::pair<
                  const DB::key_type &, const Exp *> &arg) {
                OptimizeOne(&db, arg.first, arg.second);
              },
              MAX_THREADS);

  Util::WriteFile("optimized.txt", db.Dump());

  printf("OK\n");
  return 0;
}
