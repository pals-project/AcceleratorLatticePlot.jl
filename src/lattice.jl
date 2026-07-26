# Extraction of a flat, plotting-friendly view of an expanded PALS lattice.
#
# Plotting reads `lat.full_expanded`, not `lat.expanded`: the two hold the same
# lattice, but `expanded` is pruned back to what the author wrote, and placement
# is never written by hand. `FloorP`, `s_position` and the derived members of
# `BendP` exist only in `full_expanded`.
#
# That tree nests as
#
#   <root>
#     <lattice-name>:            # kind: Lattice
#       branches:
#         - <branch-name>:
#             line:
#               - <ele-name>:    # kind: Quadrupole, Bend, Drift, ...
#                   length: ...
#                   FloorP: {x, y, z, theta, phi, psi}   # UPSTREAM (entrance) end
#                   BendP:  {angle_ref, ...}             # bends only
#                   s_position: ...
#               - branch_end:    # kind: Placeholder, appended by the bookkeeper
#                   FloorP: ...  # the downstream end of the last element
#                   s_position: ...
#
# `ElementTable` flattens every element of every branch into parallel arrays
# (struct-of-arrays) so the geometry and render stages can work on contiguous
# `Vector`s rather than chasing pointers through the tree for each of what may be
# tens of thousands of elements. Each element keeps a `YAMLNode` handle so the
# GUI can list its full parameter set on demand.

using PALSJulia
import PALSJulia as pj

"""
    ElementTable

Flat, struct-of-arrays view of every lattice element in an expanded lattice.
Field `i` of each vector describes the same element. `node[i]` is the element's
`YAMLNode`, kept so the GUI can display its parameters.

`FloorP` values are the element's **upstream (entrance)** end in the Bmad global
reference system `(x, y, z)` with heading `theta` (rotation about the vertical
`y` axis). `angle` is the total bend angle in radians (0 for straight elements).
"""
struct ElementTable
  name::Vector{String}
  kind::Vector{String}
  length::Vector{Float64}
  x::Vector{Float64}       # global floor coordinates at the entrance end
  y::Vector{Float64}
  z::Vector{Float64}
  theta::Vector{Float64}   # heading: rotation about the vertical (y) axis
  angle::Vector{Float64}   # total bend angle (rad); 0 when straight
  s::Vector{Float64}       # longitudinal s position
  branch::Vector{Int}      # 1-based index into `branch_names`
  node::Vector{pj.YAMLNode}
  branch_names::Vector{String}
end

Base.length(t::ElementTable) = length(t.name)

# Value of a direct scalar child as Float64, or `default` if absent/unparseable.
# Note: the wrapper's `is_scalar` can report false for a keyed numeric leaf whose
# text `String` still returns correctly (e.g. "1e+02"), so gate only on the node
# not being a map/sequence and let `String` + `tryparse` do the rest.
function _num(node::pj.YAMLNode, key::String, default::Float64=0.0)
  pj.haskey(node, key) || return default
  child = node[key]
  (pj.is_map(child) || pj.is_sequence(child)) && return default
  s = try
    String(child)
  catch
    return default
  end
  v = tryparse(Float64, s)
  return v === nothing ? default : v
end

# The value node of a single-key wrapper map `{name: value}`, plus its key.
# Sequence entries in `branches`/`line` are exactly these one-key maps.
function _unwrap(entry::pj.YAMLNode)
  ks = pj.keys(entry)
  isempty(ks) && return (nothing, nothing)
  k = ks[1]
  return (k, entry[k])
end

# Collect the Lattice nodes under `root`, whether `root` is itself a Lattice or a
# container map whose children are lattices.
function _lattices(root::pj.YAMLNode)
  out = pj.YAMLNode[]
  pj.is_map(root) || return out
  _is_lat(n) = pj.is_map(n) && pj.haskey(n, "kind") && String(n["kind"]) == "Lattice"
  if _is_lat(root)
    push!(out, root)
  else
    for (_, v) in root
      _is_lat(v) && push!(out, v)
    end
  end
  return out
end

"""
    element_table(lat::Lattices) -> ElementTable

Flatten every element of every branch of `lat.full_expanded` into an
`ElementTable`. Elements are appended branch by branch in lattice order.

`full_expanded` rather than `expanded` because only it carries the parameters
the bookkeeper computed — `FloorP`, `s_position` and the bend geometry — which
is everything the drawing is placed from.
"""
function element_table(lat::pj.Lattices)
  name = String[]; kind = String[]; len = Float64[]
  x = Float64[]; y = Float64[]; z = Float64[]
  theta = Float64[]; angle = Float64[]; s = Float64[]
  branch = Int[]; node = pj.YAMLNode[]; bnames = String[]

  for latnode in _lattices(lat.full_expanded)
    pj.haskey(latnode, "branches") || continue
    for bentry in latnode["branches"]
      bname, bmap = _unwrap(bentry)
      (bmap === nothing || !pj.is_map(bmap) || !pj.haskey(bmap, "line")) && continue
      push!(bnames, bname)
      bidx = length(bnames)
      for eentry in bmap["line"]
        ename, emap = _unwrap(eentry)
        (emap === nothing || !pj.is_map(emap)) && continue

        ekind = pj.haskey(emap, "kind") ? String(emap["kind"]) : ""
        L = _num(emap, "length")

        # Floor coordinates (entrance end). Skip elements that carry none; they
        # cannot be placed.
        pj.haskey(emap, "FloorP") || continue
        fp = emap["FloorP"]
        push!(name, ename); push!(kind, ekind); push!(len, L)
        push!(x, _num(fp, "x")); push!(y, _num(fp, "y")); push!(z, _num(fp, "z"))
        push!(theta, _num(fp, "theta"))
        push!(s, _num(emap, "s_position"))
        push!(branch, bidx); push!(node, emap)

        # Bend angle, read the way the expander's own floor propagation reads it
        # (`element_LS` in pals_expand.cpp): `angle_ref`, which the bookkeeper
        # derives from whichever pair of curvature/length/chord the author gave,
        # falling back to g_ref * length for a bend it could not reduce.
        a = 0.0
        if pj.haskey(emap, "BendP")
          bp = emap["BendP"]
          a = pj.haskey(bp, "angle_ref") ? _num(bp, "angle_ref") : _num(bp, "g_ref") * L
        end
        push!(angle, a)
      end
    end
  end

  return ElementTable(name, kind, len, x, y, z, theta, angle, s, branch, node, bnames)
end
