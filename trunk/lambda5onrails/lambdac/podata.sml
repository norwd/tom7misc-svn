(* primop data *)
structure Podata =
struct

  exception Podata of string

  (* base types that describe arguments to and return values from primops *)
  datatype potype =
    PT_INT | PT_STRING | PT_REF of potype | PT_UNITCONT | PT_BOOL | PT_VAR of Variable.var
  | PT_UNIT

  local open Primop
        fun mono (dom, cod) = { worlds = nil, tys = nil, dom = dom, cod = cod }
  in

    fun potype (PJointext i) = { worlds = nil, tys = nil,
                                 dom = List.tabulate (i, Util.K PT_STRING),
                                 cod = PT_STRING }
      | potype PHalt =
          let val a = Variable.namedvar "a"
          in { worlds = nil, tys = [a], dom = nil, cod = PT_VAR a }
          end
      | potype (B (PCmp _)) = mono ([PT_INT, PT_INT], PT_BOOL)
      | potype (B PTimes) = mono ([PT_INT, PT_INT], PT_INT)
      | potype (B PPlus) = mono ([PT_INT, PT_INT], PT_INT)
      | potype (B PMinus) = mono ([PT_INT, PT_INT], PT_INT)
      | potype PEqs = mono ([PT_STRING, PT_STRING], PT_BOOL)

      | potype PSet = 
          let val a = Variable.namedvar "a"
          in { worlds = nil, tys = [a], dom = [PT_REF (PT_VAR a), PT_VAR a], 
               cod = PT_UNIT }
          end

      | potype PGet =
          let val a = Variable.namedvar "a"
          in { worlds = nil, tys = [a], dom = [PT_REF (PT_VAR a)], cod = PT_VAR a }
          end

      | potype PRef =
          let val a = Variable.namedvar "a"
          in { worlds = nil, tys = [a], dom = [PT_VAR a], cod = PT_REF (PT_VAR a) }
          end

      | potype p = raise Podata ("unimplemented potype " ^ tostring p)

  end (* local *)
end