signature PACTOM =
sig

  exception PacTom of string

  (* Pactom data *)
  type pactom
  val paths : pactom -> (LatLon.pos * real) Vector.vector Vector.vector
  val overlays : pactom -> { href : string,
                             rotation : real,
                             alpha : Word8.word,
                             north : real, south : real,
                             east : real, west : real } Vector.vector

  val projectpaths : LatLon.projection -> pactom ->
      { (* x, y, elevation *)
        paths : (real * real * real) Vector.vector Vector.vector,
        (* bounding rectangle in projected space *)
        bounds : Bounds.bounds }


  val home : LatLon.pos

  val fromkmlfile : string -> pactom
  val fromkmlfiles : string list -> pactom

  type neighborhoods = (string * LatLon.pos Vector.vector) Vector.vector

  (* Loads neighborhoods kml, normalizes, and returns the data.
     The kml must consist only of <LinearRing> geometric elements,
     because it assumes that all <coordinates> tags contain polygons. *)
  val neighborhoodsfromkml : string -> neighborhoods

  val projecthoods : LatLon.projection -> neighborhoods ->
    { borders : (string * (real * real) Vector.vector) Vector.vector,
      bounds : Bounds.bounds }

  (* GPS trails are not usually efficient drawing paths, in that they
     record points even when traveling in a straight line. Removes from
     points from paths when doing so introduces no more than the given
     number of square meters of error (approximate; assumes flat earth).
     Ignores elevations. For the removal candidate X:

                               *
                    .X.       /
                 .-`   `.    /
              .-`  error `. /
             *-------------*
            /
           /
          *
     *)
  val simplify_paths : real -> pactom -> pactom

  (* Distances are accurate great-circle distances, ignoring elevation. *)
  structure G : UNDIRECTEDGRAPH where type weight = real

  (* As indices into the vector of paths, and into those paths. *)
  type waypoint = { path : int, pt : int }
  (* Put every point into a data structure for querying by location.
     Each is associated with the path it's from, by index into the vector. *)
  val latlontree : pactom -> waypoint LatLonTree.tree

  (* Compute the graph of reachability, which is all of the points in the
     paths, connected when they are close enough to imply connectedness. *)
  val graph : pactom -> { graph : waypoint G.graph,
                          (* Same indices as in original paths. *)
                          promote : waypoint -> waypoint G.node,
                          (* Same as above. *)
                          latlontree : waypoint LatLonTree.tree }

  (* Like graph, but adds and moves points to simplify it.
     Currently, nearly equivalent edges are merged, preserving
     connectedness. XXX: This should add points at intersections,
     then also merge edges that are nearly overlapping and colinear.

     Because this changes the set of points, the  XXX

  val simplified_graph : pactom -> { graph : waypoint G.graph,
                                     promote : waypoint -> waypoint G.node,
                                     latlontree : waypoint LatLonTree.tree }
     *)

  (* A minimal spanning tree, rooted at the node closest to the given
     position.
     For each node, it has its distance back to the root (or NONE, if
     unreachable), and if it is reachable, then the parent node to
     take to get there (except for the root, which is NONE). *)
  val spanning_graph : pactom -> LatLon.pos ->
                         { graph : waypoint G.span G.graph,
                           promote : waypoint -> waypoint G.span G.node }

  val rand : unit -> Word32.word
  val randf : unit -> real
  (* As hex string (rrggbb) *)
  val randomanycolor : unit -> string
  (* High saturation and value, random hue. *)
  val randombrightcolor : unit -> string
  (* Low saturation, high value, random hue. *)
  val randompalecolor : unit -> string

  (* Suitable for XML output. Avoids expontential notation
     and tilde for negative. (XXX should also cap out or fail on inf, nan) *)
  val rtos : real -> string

end