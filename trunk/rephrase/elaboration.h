
#ifndef _REPHRASE_ELABORATION_H
#define _REPHRASE_ELABORATION_H

#include "el.h"
#include "il.h"

#include "context.h"
#include "initial.h"

// Elaboration is a recursive transformation from EL to IL.
struct Elaboration {

  Elaboration(il::AstPool *pool) : pool(pool), init(pool) {}

  void SetVerbose(int v) { verbose = v; }

  const il::Exp *Elaborate(const el::Exp *el_exp);

private:

  const std::pair<const il::Exp *, const il::Type *> Elab(
      const il::Context &G,
      const el::Exp *el_exp);

  const il::Type *ElabType(
      const il::Context &G,
      const el::Type *el_type);

  const il::Type *NewEVar();

  const std::pair<const il::Exp *, const il::Type *> ElabPat(
      const il::Exp *target,
      const el::Pat *pat,
      const il::Type *type,
      std::function<
        std::pair<const il::Exp *, const il::Type *>(const il::Context &G)>
      cont
      // TODO: failure continuation?
                                                             );

  int verbose = 0;
  il::AstPool *pool;
  il::Initial init;
};

#endif

