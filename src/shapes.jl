# Mapping lattice elements to the shape, color, size and label used to draw them.
#
# This follows the Tao `floor_plan_drawing` model (see the Tao manual, "Lat_layout
# and Floor_plan Drawings Shape Definition"): an ordered list of rules, each of
# the form `<class>::<name-glob> -> shape, color, size, label`. The first rule an
# element matches wins. A built-in default table is used for any element not
# matched by a user rule.

using Colors

# The shape primitives we know how to draw. Curved (bend) elements follow their
# arc regardless of which of these is chosen.
const SHAPE_NAMES = (:box, :xbox, :x, :bow_tie, :rbow_tie, :diamond,
                     :u_triangle, :d_triangle, :l_triangle, :r_triangle,
                     :circle, :none)

"""
How to draw one element.

`size` is the horizontal transverse half-height in meters and `size2` the
vertical one: the 2D drawing is a profile in the plane of `size`, and the 3D
drawing extrudes that same profile by `size2` either side of the centerline, so
one table serves both views. `size2` defaults to `size`, making an element square
in cross-section.
"""
struct ShapeSpec
  shape::Symbol
  color::Colorant
  size::Float64
  size2::Float64
  label::Symbol       # :name, :none, or :s
  draw::Bool
  line_width::Float64
end

ShapeSpec(shape, color; size=0.5, size2=nothing, label=:name, draw=true,
          line_width=1.5) =
  ShapeSpec(shape, parse(Colorant, string(color)), size,
            size2 === nothing ? Float64(size) : Float64(size2),
            Symbol(label), draw, Float64(line_width))

"One match rule: element class (or `*`) plus a name glob, mapped to a `ShapeSpec`."
struct ShapeRule
  class::String        # element `kind`, matched case-insensitively; "*" = any
  glob::String         # element name glob (`*`, `?`); matched case-insensitively
  spec::ShapeSpec
end

"""
    ele_shape(ele_id, shape, color; size, size2, label, draw, line_width) -> ShapeRule

Build a shape rule. `ele_id` is `"<class>::<name-glob>"` (e.g. `"Quadrupole::q*"`)
or just `"<class>"` (implying `::*`). `"*"` for the class matches any class.
`size2` sets the vertical half-height used by the 3D drawing, and defaults to
`size`.
"""
function ele_shape(ele_id::AbstractString, shape, color; kwargs...)
  class, glob = occursin("::", ele_id) ? split(ele_id, "::"; limit=2) : (ele_id, "*")
  return ShapeRule(String(class), String(glob), ShapeSpec(Symbol(shape), color; kwargs...))
end

# Translate a `*`/`?` glob to an anchored, case-insensitive regex, once per rule.
function _glob_regex(glob::AbstractString)
  glob == "*" && return r""  # sentinel: matches everything (handled in _match)
  buf = IOBuffer(); print(buf, "^")
  for c in glob
    if c == '*';     print(buf, ".*")
    elseif c == '?'; print(buf, ".")
    else             print(buf, "\\Q", c, "\\E")
    end
  end
  print(buf, "\$")
  return Regex(String(take!(buf)), "i")
end

_class_matches(rule_class, kind) =
  rule_class == "*" || lowercase(rule_class) == lowercase(kind)

"""
    ShapeMap(rules; defaults=true)

An ordered collection of `ShapeRule`s. Call `mapshape(map, kind, name)` to get
the `ShapeSpec` for an element. User `rules` are tried first; if `defaults` is
true the built-in default rules are appended.
"""
struct ShapeMap
  rules::Vector{ShapeRule}
  regexes::Vector{Regex}
end

function ShapeMap(rules::Vector{ShapeRule}=ShapeRule[]; defaults::Bool=true)
  all = defaults ? vcat(rules, default_shapes()) : copy(rules)
  return ShapeMap(all, [_glob_regex(r.glob) for r in all])
end

"""
    mapshape(m::ShapeMap, kind, name) -> ShapeSpec

The `ShapeSpec` for the first rule matching an element of class `kind` and name
`name`. Falls back to a neutral gray box if nothing matches.
"""
function mapshape(m::ShapeMap, kind::AbstractString, name::AbstractString)
  for (r, rx) in zip(m.rules, m.regexes)
    _class_matches(r.class, kind) || continue
    (r.glob == "*" || occursin(rx, name)) && return r.spec
  end
  return ShapeSpec(:box, colorant"gray50", 0.3, 0.3, :none, true, 1.0)
end

"""
    default_shapes() -> Vector{ShapeRule}

Built-in element-class → shape table covering the element kinds PALS defines
(`lattice-element-kinds.md`), loosely following Tao's floor-plan conventions:
bends blue, quadrupoles black boxed-X, sextupoles green, cavities red, and so
on. Drifts and the zero-length bookkeeping elements draw only the centerline
(`:none`).

An element of a kind not listed here — `Feedback` and `Girder`, which are not
part of a beam line, or anything the standard adds later — falls through to the
`"*"` rule and is drawn as a small unlabeled box.
"""
function default_shapes()
  R(class, shape, color; kw...) = ele_shape(class, shape, color; kw...)
  return ShapeRule[
    # Magnets and RF cavities
    R("Bend",        :box,     :blue,     size=0.5),
    R("Quadrupole",  :xbox,    :black,    size=0.5),
    R("Sextupole",   :box,     :green,    size=0.4),
    R("Octupole",    :box,     :green,    size=0.4),
    R("Solenoid",    :box,     :cyan,     size=0.5),
    R("Wiggler",     :box,     :green,    size=0.5),
    R("Multipole",   :diamond, :purple,   size=0.3),
    R("Kicker",      :box,     :magenta,  size=0.3),
    R("ACKicker",    :box,     :magenta,  size=0.3),
    R("RFCavity",    :box,     :red,      size=0.5),
    R("CrabCavity",  :box,     :red,      size=0.4),
    # Beam and plasma
    R("BeamBeam",    :diamond, :orange,   size=0.4),
    # Sources and collimation
    R("EGun",        :r_triangle, :red,   size=0.5),
    R("Converter",   :r_triangle, :sienna, size=0.4),
    R("Foil",        :box,     :sienna,   size=0.2),
    R("Mask",        :bow_tie, :sienna,   size=0.4),
    # Instrumentation
    R("Instrument",  :box,     :darkgreen, size=0.3),
    # Map elements
    R("Match",       :box,     :gray40,   size=0.3),
    R("Taylor",      :box,     :gray40,   size=0.3),
    # Bookkeeping: the ones that redirect the reference curve are worth seeing,
    # the rest only mark a point on it.
    R("Fork",        :u_triangle, :orange, size=0.5),
    R("Patch",       :rbow_tie, :gray40,  size=0.3),
    R("FloorShift",  :x,       :gray40,   size=0.3),
    R("ReferenceChange", :x,   :gray40,   size=0.3),
    R("Fiducial",    :x,       :gray40,   size=0.3),
    R("BeginningEle", :none,   :gray60,   size=0.0, label=:none),
    R("Marker",      :none,    :gray60,   size=0.0, label=:none),
    R("Placeholder", :none,    :gray60,   size=0.0, label=:none),
    R("Drift",       :none,    :gray60,   size=0.0, label=:none),
    # Grouping
    R("UnionEle",    :box,     :gray40,   size=0.4),
    R("*",           :box,     :black,    size=0.3, label=:none),
  ]
end
