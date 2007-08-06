(* Sometimes it is useful for debugging to be able to
   have more memorable names than "x$2679". This structure
   facilitates that. *)

structure Nonce =
struct
  
  val words = ref
   ["cat", "dog", "toy", "robot", "shirt", "house", "helicopter",
    "comb", "skeletor", "kimono", "beard", "police", "olestra",
    "umbra", "cedar", "crow", "chagrin", "fourier", "kilt",
    "dinosaur", "erasmus", "monopoly", "liver", "lion", "politics",
    "lemur", "lepus", "crayon", "velvet", "languid", "quiche",
    "photon", "graviton", "muon", "lepton", "electron", "shock",
    "uracil", "methane", "pentane", "hexane", "butane", "kinesin",
    "opalescent", "cloud", "child", "icon", "kitten", "lilac",
    "topaz", "ruby", "sapphire", "hubris", "fluorite", "myriad",
    "cheese", "onions", "phlegm", "melisma", "jurist", "hacienda"]

  val ctr = ref 0

  fun nonce () =
    case !words of
      nil => (ctr := !ctr + 1; "nonce$" ^ Int.toString (!ctr))
    | h :: t => (words := t; h)

end