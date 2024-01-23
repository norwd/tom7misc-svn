
#ifndef _REPHRASE_INITIAL_H
#define _REPHRASE_INITIAL_H

#include "context.h"
#include "il.h"

struct Initial {
  explicit Initial(il::AstPool *pool);

  const Context &InitialContext() const;

  // Types needed by elaboration.

private:
  Context ctx;
};

#endif
