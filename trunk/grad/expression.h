
#ifndef _GRAD_EXPRESSION_H
#define _GRAD_EXPRESSION_H

#include <cstdint>
#include <vector>
#include <array>
#include <string>

#include "base/logging.h"
#include "base/stringprintf.h"
#include "half.h"

using half_float::half;

// Expression of one variable.

enum ExpType {
  VAR,
  PLUS_C,
  TIMES_C,
  PLUS_E,
};

struct Exp {
  ExpType type;
  uint16_t c = 0x0000;
  // Plus_c and Times_c can be iterated.
  uint16_t iters = 1;
  const Exp *a = nullptr, *b = nullptr;

  struct Allocator {

    // TODO: Some way to release/copy an expression!
    ~Allocator() {
      for (Exp *e : allocations) delete e;
      allocations.clear();
    }

    // TODO: Could index variable (using c, perhaps) to allow for
    // linear functions of multiple variables.
    static const Exp *Var() {
      return new Exp(VAR);
    }

    static const Exp *PlusC(const Exp *e, uint16_t c, uint16_t iters = 1) {
      Exp *ret = new Exp(PLUS_C);
      ret->a = e;
      ret->c = c;
      ret->iters = iters;
      return ret;
    }

    static const Exp *TimesC(const Exp *e, uint16_t c, uint16_t iters = 1) {
      Exp *ret = new Exp(TIMES_C);
      ret->a = e;
      ret->c = c;
      ret->iters = iters;
      return ret;
    }

    // TODO: Verify that this is equivalent to unary negation.
    static const Exp *Neg(const Exp *e) {
      Exp *ret = new Exp(TIMES_C);
      ret->a = e;
      ret->c = 0xbc00;  // -1.0
      return ret;
    }

    static const Exp *PlusE(const Exp *a, const Exp *b) {
      Exp *ret = new Exp(PLUS_E);
      ret->a = a;
      ret->b = b;
      return ret;
    }

  private:
    inline Exp *New(ExpType t) {
      Exp *e = new Exp(t);
      allocations.push_back(e);
      return e;
    }
    std::vector<Exp *> allocations;
  };

  static inline half GetHalf(uint16_t u) {
    half h;
    static_assert(sizeof (h) == sizeof (u));
    memcpy((void*)&h, (void*)&u, sizeof (u));
    return h;
  }

  static inline uint16_t GetU16(half h) {
    uint16_t u;
    static_assert(sizeof (h) == sizeof (u));
    memcpy((void*)&u, (void*)&h, sizeof (u));
    return u;
  }

  using Table = std::array<uint16_t, 65536>;

  struct TimesTable {
    TimesTable(uint16_t cu, int iters) {
      // PERF this can be computed by squaring, right?
      // (And we can certainly compute multiple tables at once.)
      half c = GetHalf(cu);
      for (int i = 0; i < 65536; i++) {
        uint16_t xu = i;
        half y = GetHalf(xu);
        for (int z = 0; z < iters; z++)
          y *= c;
        table[xu] = GetU16(y);
      }
    }
    Table table;
  };

  // Fast iteration of * 0x3bff.
  static uint16_t FastTimes0x3bff(uint16_t lhs, int iters) {
    static TimesTable *times_table100 =
      new TimesTable(0x3bff, 100);

    static TimesTable *times_table10 =
      new TimesTable(0x3bff, 10);


    while (iters > 100) {
      lhs = times_table100->table[lhs];
      iters -= 100;
    }

    while (iters > 10) {
      lhs = times_table10->table[lhs];
      iters -= 10;
    }

    half y = GetHalf(lhs);
    const half c = GetHalf(0x3bff);
    for (int z = 0; z < iters; z++)
      y *= c;
    return GetU16(y);
  }

  // PERF: For iterated plus/times, we can have precomputed
  // tables that apply the function e.g. 100 times.
  static uint16_t EvaluateOn(const Exp *e, uint16_t x) {
    switch (e->type) {
    case VAR: return x;
    case PLUS_C: {
      half res = GetHalf(EvaluateOn(e->a, x));
      half rhs = GetHalf(e->c);
      for (int i = 0; i < e->iters; i++)
        res += rhs;
      return GetU16(res);
    }
    case TIMES_C: {
      uint16_t lhs = EvaluateOn(e->a, x);

      if (e->c == 0x3bffu) {
        return FastTimes0x3bff(lhs, e->iters);
      }

      half res = GetHalf(lhs);
      half rhs = GetHalf(e->c);
      for (int i = 0; i < e->iters; i++)
        res *= rhs;
      return GetU16(res);
    }
    case PLUS_E:
      return GetU16(GetHalf(EvaluateOn(e->a, x)) +
                    GetHalf(EvaluateOn(e->b, x)));
    default:
      CHECK(false) << "Unknown expression type";
      return 0;
    }
  }

  static Table TabulateExpression(const Exp *e) {
    Table ret;
    // PERF: Can skip running on nans.
    for (int x = 0; x < 65536; x++) {
      uint16_t y = EvaluateOn(e, x);
      ret[x] = y;
    }
    return ret;
  }

  // Only fills the table in the range [low, high].
  static Table TabulateExpressionIn(const Exp *e,
                                    half low, half high) {
    Table ret;
    CHECK(low < high);
    for (half pos = low; pos <= high; /* in loop */) {
      half next = nextafter(pos, HUGE_VALH);
      uint16_t upos = GetU16(pos);

      ret[upos] = EvaluateOn(e, upos);

      pos = next;
    }
    return ret;
  }

  static std::string ExpString(const Exp *e) {
    switch (e->type) {
    case VAR: return "V";
    case PLUS_C: {
      std::string lhs = ExpString(e->a);
      return StringPrintf("P(%s,0x%04x,%d)",
                          lhs.c_str(), e->c, e->iters);
    }
    case TIMES_C: {
      std::string lhs = ExpString(e->a);
      return StringPrintf("T(%s,0x%04x,%d)",
                          lhs.c_str(), e->c, e->iters);
    }
    case PLUS_E: {
      std::string lhs = ExpString(e->a);
      std::string rhs = ExpString(e->b);
      return StringPrintf("E(%s,%s)", lhs.c_str(), rhs.c_str());
    }
    default:
      CHECK(false) << "Unknown expression type";
      return "";
    }
  }

private:
  Exp(ExpType t) : type(t) {}
};


#endif
