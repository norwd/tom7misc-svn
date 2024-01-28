
#include "initial.h"

#include <vector>
#include <string>
#include <utility>

namespace il {

Initial::Initial(AstPool *pool) {

  auto PairType = [&](const Type *a, const Type *b) {
      return pool->Product({a, b});
    };

  auto BinOpType = [&](const Type *a, const Type *b, const Type *ret) ->
    const Type * {
      return pool->Arrow(PairType(a, b), ret);
    };

  auto Mono = [&](const Type *t) -> PolyType {
      return PolyType{.tyvars = {}, .type = t};
    };

  const il::Type *Unit = pool->RecordType({});
  const il::Type *Alpha = pool->VarType("a");
  const il::Type *Int = pool->IntType();
  auto Ref = [&](const Type *a) { return pool->RefType(a); };
  // This is probably wrong: We need to expand the type of list,
  // or better
  auto List = [&](const Type *a) { return pool->VarType("list", {a}); };

  const std::vector<std::pair<std::string, PolyType>> exp_vars = {
    {"+", Mono(BinOpType(Int, Int, Int))},
    {"-", Mono(BinOpType(Int, Int, Int))},
    {":=", PolyType{
        .tyvars = {"a"},
        .type = BinOpType(Ref(Alpha), Alpha, Unit)}},
    {"!", PolyType{
        .tyvars = {"a"},
        .type = pool->Arrow(Ref(Alpha), Alpha)}},
    {"ref", PolyType{
        .tyvars = {"a"},
        .type = pool->Arrow(Alpha, Ref(Alpha))}},

    /*
    {"SOME", PolyType{
        .tyvars = {"a"},
        .type = pool->Arrow(Alpha, Option(Alpha))}},
    */

    // TODO: List stuff can just be introduced by a preamble.
    {"::", PolyType{
        .tyvars = {"a"},
        .type = BinOpType(Alpha, List(Alpha), List(Alpha))}},

  };

  const il::Type *String = pool->StringType();
  auto Kind0 = [&](const Type *t) {
      return SingletonKind{.tyvars = {}, .type = t};
    };

  const std::vector<std::pair<std::string, SingletonKind>> type_vars = {
    // These need datatype declarations.
    //    {"list", 1},
    //    {"option", 1},
    //    {"bool", 0},
    {"int", Kind0(Int)},
    {"string", Kind0(String)},
    {"ref", SingletonKind{.tyvars = {"a"}, .type = Ref(Alpha)}},
  };

  ctx = Context(exp_vars, type_vars);
}

const Context &Initial::InitialContext() const { return ctx; }

}  // il
