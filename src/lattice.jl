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
#                   BendP:  {angle_ref, tilt_ref, ...}   # bends only
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

using PALSParserJ
import PALSParserJ as pj

"""
    ElementTable

Flat, struct-of-arrays view of every lattice element in an expanded lattice.
Field `i` of each vector describes the same element. `node[i]` is the element's
`YAMLNode`, kept so the GUI can display its parameters.

`FloorP` values are the element's **upstream (entrance)** end in the global
reference system: position `(x, y, z)` and the three orientation angles `theta`
(azimuth), `phi` (pitch) and `psi` (roll), which together give the orientation
matrix `W = R_y(theta) R_x(phi) R_z(psi)` (see [`w_matrix`](@ref)). `angle` is
the total reference bend angle in radians (0 for straight elements) and
`tilt_ref` the bend's reference tilt, which rolls its bend plane.

All six placement quantities are needed to draw a lattice that leaves the
horizontal plane; a lattice that stays in it has `phi = psi = tilt_ref = 0`
throughout.
"""
struct ElementTable
  name::Vector{String}
  kind::Vector{String}
  length::Vector{Float64}
  x::Vector{Float64}         # global floor coordinates at the entrance end
  y::Vector{Float64}
  z::Vector{Float64}
  theta::Vector{Float64}     # orientation angles at the entrance end: azimuth,
  phi::Vector{Float64}       # pitch,
  psi::Vector{Float64}       # roll
  angle::Vector{Float64}     # total reference bend angle (rad); 0 when straight
  tilt_ref::Vector{Float64}  # bend reference tilt (rad); rolls the bend plane
  s::Vector{Float64}         # longitudinal s position
  branch::Vector{Int}        # 1-based index into `branch_names`
  node::Vector{pj.YAMLNode}
  branch_names::Vector{String}
end

"""
    ElementTable(name, kind, length, x, y, z, theta, angle, s, branch, node,
                 branch_names)

Build a table for a lattice that lies in the horizontal plane, taking
`phi = psi = tilt_ref = 0`. Convenient for synthesizing a table by hand; the
extraction path fills in all of them.
"""
function ElementTable(name, kind, len, x, y, z, theta, angle, s, branch, node,
                      branch_names)
  n = Base.length(name)
  return ElementTable(name, kind, len, x, y, z, theta, zeros(n), zeros(n),
                      angle, zeros(n), s, branch, node, branch_names)
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
  theta = Float64[]; phi = Float64[]; psi = Float64[]
  angle = Float64[]; tilt = Float64[]; s = Float64[]
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
        # All three orientation angles: `theta` alone places only a lattice that
        # stays in the horizontal plane.
        push!(theta, _num(fp, "theta"))
        push!(phi, _num(fp, "phi"))
        push!(psi, _num(fp, "psi"))
        push!(s, _num(emap, "s_position"))
        push!(branch, bidx); push!(node, emap)

        # Bend angle and tilt, read the way the expander's own floor propagation
        # reads them (`element_LS` in pals_expand.cpp): `angle_ref`, which the
        # bookkeeper derives from whichever pair of curvature/length/chord the
        # author gave, falling back to g_ref * length for a bend it could not
        # reduce, and `tilt_ref`, which rolls the plane the arc bends in.
        a = 0.0; tl = 0.0
        if pj.haskey(emap, "BendP")
          bp = emap["BendP"]
          a = pj.haskey(bp, "angle_ref") ? _num(bp, "angle_ref") : _num(bp, "g_ref") * L
          tl = _num(bp, "tilt_ref")
        end
        push!(angle, a); push!(tilt, tl)
      end
    end
  end

  return ElementTable(name, kind, len, x, y, z, theta, phi, psi, angle, tilt, s,
                      branch, node, bnames)
end
