
#ifndef _CC_LIB_ATOMIC_UTIL_H
#define _CC_LIB_ATOMIC_UTIL_H

#include <cstdint>
#include <atomic>

// For performance reasons, eight counters are defined in
// batch. Declare some counters (typically at file scope in
// a .cc file) like this:
//
// DECLARE_COUNTERS(bytes, lines, errors, u1_, u2_, u3_, u4_, u5_);
//
// Then you have
//   bytes++;
//   printf("%llu\n", bytes.Read());
//   bytes.Reset();
//
// Only incrementing is efficient here, but it should be a lot
// faster than std::atomic<uint64_t> (or using a mutex) when there
// is a lot of contention.

#define DECLARE_COUNTERS(a, b, c, d, e, f, g, h) \
  [[maybe_unused]] static internal::EightCounters ec_ ## a;        \
  [[maybe_unused]] static AtomicCounter a(&ec_ ## a, 0);           \
  [[maybe_unused]] static AtomicCounter b(&ec_ ## a, 1);           \
  [[maybe_unused]] static AtomicCounter c(&ec_ ## a, 2);           \
  [[maybe_unused]] static AtomicCounter d(&ec_ ## a, 3);           \
  [[maybe_unused]] static AtomicCounter e(&ec_ ## a, 4);           \
  [[maybe_unused]] static AtomicCounter f(&ec_ ## a, 5);           \
  [[maybe_unused]] static AtomicCounter g(&ec_ ## a, 6);           \
  [[maybe_unused]] static AtomicCounter h(&ec_ ## a, 7)

namespace internal {
// Based on a great article by Travis Downs.
// https://travisdowns.github.io/blog/2020/07/06/concurrency-costs.html
class EightCounters {
public:
  EightCounters() {
    for (size_t off = 0; off < 8; off++) {
      Reset(off);
    }
  }

  // Increment the logical counter value.
  uint64_t Increment(size_t off) {
    // must be one of the 8 counters
    if (off != (off & 7)) __builtin_unreachable();

    for (;;) {
      std::atomic<uint64_t> &counter = buckets[idx].counters[off];

      // Try storing without a lock.
      uint64_t cur = counter.load();
      if (counter.compare_exchange_strong(cur, cur + 1)) {
        return cur;
      }

      // CAS failure indicates contention,
      // so try again at a different index.
      // PERF: Would be best if different threads used
      // different strategies (i.e. generators) here.
      idx = (idx + 1) % NUM_BUCKETS;
    }
  }

  // Get the counter's value. Has to sum up all the buckets.
  uint64_t Read(size_t off) const {
    // must be one of the 8 counters
    if (off != (off & 7)) __builtin_unreachable();

    uint64_t sum = 0ULL;
    for (const Cacheline &line : buckets) {
      sum += line.counters[off].load();
    }
    return sum;
  }

  // Reset the counter's value to zero.
  void Reset(size_t off) {
    // must be one of the 8 counters
    if (off != (off & 7)) __builtin_unreachable();

    for (size_t i = 0; i < NUM_BUCKETS; i++) {
      buckets[i].counters[off].store(0ULL);
    }
  }

 private:
  // Essentially the maximum number of threads that can be
  // concurrently accessing the counter without contending. OK to
  // raise this, although it consumes nontrivial memory (since there
  // are 64 bytes per bucket). It's also okay to lower it; we
  // just get more contention when that happens.
  static constexpr size_t NUM_BUCKETS = 32;

  // Each bucket is a single cacheline, so that we don't get
  // false contention across threads. But we can use up the
  // 64 bytes for eight 64-bit counters.
  struct Cacheline {
    alignas(64) std::atomic<uint64_t> counters[8];
  };

  // Each thread gets its own index; these start at zero but move
  // when contention is observed.
  static thread_local size_t idx;
  Cacheline buckets[NUM_BUCKETS];
};
}  // namespace internal

// Represents one of the 8 slots in an EightCounters instance.
class AtomicCounter {
 public:
  inline uint64_t operator++(int suffix_) {
    return ec->Increment(offset);
  }

  inline uint64_t Read() const {
    return ec->Read(offset);
  }

  inline void Reset() {
    ec->Reset(offset);
  }

  // Use the macro.
  AtomicCounter(internal::EightCounters *ec, size_t offset) :
    ec(ec), offset(offset) {}

 private:
  internal::EightCounters *ec = nullptr;
  // 0-7
  size_t offset = 0;
};

#endif
