# Turning an ElementTable + ShapeMap into batched 3D draw data.
#
# Same idea as geometry.jl, and the same shape table: an element's profile in
# normalized `(f, w)` coordinates is swept along its reference curve. Where the
# 2D drawing projects that sweep onto a plane, this one extrudes it by the
# shape's `size2` either side of the centerline, giving a solid whose silhouette
# from above is the floor plan.
#
# Everything lands in three buffers so the whole lattice is a fixed number of
# Makie calls regardless of element count:
#
#   * one triangle mesh    -> one `mesh!`          (every element solid)
#   * cap outlines/strokes -> one `linesegments!`  (the edges over the shading)
#   * reference curve      -> one `lines!`         (NaN-separated per branch)
#
# Triangles carry their own vertices and an explicit face normal rather than
# sharing vertices between faces. That costs memory but buys flat shading with no
# ambiguity about which way a face points -- a magnet reads as a solid object,
# and nothing is smoothed across the edge between a cap and a side wall. It also
# gives picking a vertex-to-element map for free (`vertex_ele`), which is what
# turns a GPU pick index back into an element.

using GeometryBasics: Point3f, Vec3f, Vec3d, GLTriangleFace
using Colors

const _NAN3 = Point3f(NaN32, NaN32, NaN32)

@inline dot3(u, v) = u[1] * v[1] + u[2] * v[2] + u[3] * v[3]
@inline cross3(u, v) = Vec3f(u[2] * v[3] - u[3] * v[2],
                             u[3] * v[1] - u[1] * v[3],
                             u[1] * v[2] - u[2] * v[1])

"""
    FloorGeometry3

Batched 3D draw data for a whole lattice, plus the per-element data selection
needs. Positions are in *drawn* coordinates -- the global axes permuted by the
`view` string that built it -- not in global coordinates.
"""
struct FloorGeometry3
  mesh_pts::Vector{Point3f}
  mesh_nrm::Vector{Vec3f}
  mesh_faces::Vector{GLTriangleFace}
  mesh_col::Vector{RGBA{Float32}}
  vertex_ele::Vector{Int32}          # mesh vertex -> element index, for picking

  edge_pts::Vector{Point3f}          # outline edges and interior strokes; pairs
  edge_col::Vector{RGBA{Float32}}

  ref_pts::Vector{Point3f}           # reference curve, NaN-separated per branch
  ref_col::RGBA{Float32}

  label_pos::Vector{Point3f}
  label_str::Vector{String}
  label_ele::Vector{Int}

  # Element centerline, sampled coarsely: entrance, midpoint and exit. Picking
  # measures a mouse ray against those two segments, which is exact for a
  # straight element and close enough to the arc for a bend.
  ele_entrance::Vector{Point3f}
  ele_center::Vector{Point3f}
  ele_exit::Vector{Point3f}
end

# ── convex clipping in f ──────────────────────────────────────────────────────
# A curved element's cap cannot be one flat polygon: it has to follow the arc.
# The cap polygons out of `_template` are convex, and a convex polygon clipped to
# a slab `fa <= f <= fb` is still convex, so slicing one into a strip along `f`
# is two Sutherland-Hodgman passes -- after which each slice is short enough in
# `f` to be treated as flat.

function _clip_half(poly::Vector{Tuple{Float64,Float64}}, bound::Float64, keep_ge::Bool)
  out = Tuple{Float64,Float64}[]
  isempty(poly) && return out
  inside(p) = keep_ge ? p[1] >= bound - 1.0e-12 : p[1] <= bound + 1.0e-12
  n = length(poly)
  for k in 1:n
    p = poly[k]; q = poly[k % n + 1]
    pin = inside(p); qin = inside(q)
    pin && push!(out, p)
    if pin != qin && q[1] != p[1]
      t = (bound - p[1]) / (q[1] - p[1])
      push!(out, (bound, p[2] + t * (q[2] - p[2])))
    end
  end
  return out
end

_clip_f(poly, fa, fb) = _clip_half(_clip_half(poly, fa, true), fb, false)

# Centroid of a template's profile, used to turn its side walls outward.
function _profile_centroid(tpl::ShapeTemplate)
  sf = 0.0; sw = 0.0; n = 0
  for face in tpl.faces, p in face
    sf += p[1]; sw += p[2]; n += 1
  end
  n == 0 && return (0.5, 0.0)
  return (sf / n, sw / n)
end

# ── main builder ──────────────────────────────────────────────────────────────

"""
    build_geometry3(tab::ElementTable, smap::ShapeMap; view="zxy", arc_tol=0.08,
                    circle_sides=16, edges=true) -> FloorGeometry3

Compute batched 3D draw data for the whole lattice.

`view` is a three-character permutation of the global axes `x`, `y`, `z`, one per
drawn axis, the last being the drawn vertical. The default `"zxy"` puts global
`z` and `x` in the horizontal plane with the global vertical `y` up, so looking
straight down at the result gives the default `"zx"` floor plan.

`arc_tol` is the maximum angular step (rad) a bend's arc is tessellated at;
`circle_sides` sets the facet count of `:circle` elements, which are drawn as
spheres. `edges=false` drops the outline strokes, leaving bare shaded solids.
`label_sep` works as it does in [`build_geometry`](@ref) -- how close two label
anchors must be to count as colliding, in units of the elements' half-height --
but the labels it separates are stacked vertically rather than outward, since a
billboarded label has no outward direction of its own.
"""
function build_geometry3(tab::ElementTable, smap::ShapeMap;
                         view::AbstractString="zxy", arc_tol::Real=0.08,
                         circle_sides::Int=16, edges::Bool=true,
                         label_sep::Real=1.0)
  length(view) == 3 ||
    throw(ArgumentError("a 3D `view` names three axes, got \"$view\""))
  a, b, c = view[1], view[2], view[3]

  mesh_pts = Point3f[]; mesh_nrm = Vec3f[]; mesh_faces = GLTriangleFace[]
  mesh_col = RGBA{Float32}[]; vertex_ele = Int32[]
  edge_pts = Point3f[]; edge_col = RGBA{Float32}[]
  ref_pts = Point3f[]
  label_str = String[]; label_ele = Int[]
  label_base = Point3f[]      # anchor before the collision stack is applied
  label_up = Vec3f[]          # the direction that stack runs in
  label_anchor = Point3f[]    # element midpoint, what collisions are judged on
  label_h = Float32[]
  ele_entrance = Vector{Point3f}(undef, length(tab))
  ele_center = Vector{Point3f}(undef, length(tab))
  ele_exit = Vector{Point3f}(undef, length(tab))

  # One triangle, with its own three vertices and one shared face normal. Slivers
  # are dropped: a profile that pinches to a point (a diamond's tip, a bow tie's
  # crossing) tiles into a few zero-area triangles, which draw nothing and would
  # only give a renderer something degenerate to chew on.
  function tri!(p1, p2, p3, nrm, col, i)
    hypot(cross3(p2 - p1, p3 - p1)...) < 1.0e-12 && return
    base = length(mesh_pts)
    push!(mesh_pts, p1, p2, p3)
    push!(mesh_nrm, nrm, nrm, nrm)
    push!(mesh_col, col, col, col)
    push!(vertex_ele, i, i, i)
    push!(mesh_faces, GLTriangleFace(base + 1, base + 2, base + 3))
  end

  edge!(p1, p2, col) = (push!(edge_pts, p1, p2); push!(edge_col, col, col))

  for i in 1:length(tab)
    L = tab.length[i]; ang = tab.angle[i]; tilt = tab.tilt_ref[i]
    e = entrance(tab, i)

    # Centerline point and the local axes there, in drawn coordinates. The axes
    # go through the same permutation as the point, so they stay the true branch
    # axes rather than something rebuilt in the drawn frame.
    function at(f)
      pl = place(e.r, e.W, L, ang, tilt, f)
      return (proj3(a, b, c, pl.r), _vec3(a, b, c, xaxis(pl.W)),
              _vec3(a, b, c, yaxis(pl.W)), _vec3(a, b, c, zaxis(pl.W)))
    end

    ele_entrance[i] = at(0.0)[1]
    pmid, midx, midy, _ = at(0.5)
    ele_center[i] = pmid
    ele_exit[i] = at(1.0)[1]

    # Reference curve, broken between branches.
    i > 1 && tab.branch[i] != tab.branch[i - 1] && push!(ref_pts, _NAN3)
    nsamp = max(2, ceil(Int, abs(ang) / arc_tol) + 1)
    for k in 0:(nsamp - 1)
      push!(ref_pts, at(k / (nsamp - 1))[1])
    end

    spec = mapshape(smap, tab.kind[i], tab.name[i])
    if spec.draw && spec.shape !== :none
      col = RGBA{Float32}(spec.color)
      if spec.shape === :circle
        _sphere!(tri!, edge!, pmid, spec.size, col, i, circle_sides, edges)
      else
        _solid!(tri!, edge!, at, _template(spec.shape), L, ang, spec.size,
                spec.size2, col, i, arc_tol, edges)
      end
    end

    # Label anchor: out past the shape on the +x side of the centerline, the same
    # side the floor plan puts it on, and lifted a little so it clears the solid
    # when the machine is seen edge-on.
    if spec.label !== :none
      hmax = Float32(max(spec.size, spec.size2))
      off = 1.6f0 * hmax
      push!(label_base, pmid + off * midx + 0.4f0 * off * midy)
      push!(label_up, midy)
      push!(label_anchor, pmid)
      push!(label_h, hmax)
      push!(label_str,
            spec.label === :s ? string(round(tab.s[i]; digits=3)) : tab.name[i])
      push!(label_ele, i)
    end
  end

  label_pos = _stack_labels3(label_base, label_up, label_anchor, label_h, label_sep)

  return FloorGeometry3(mesh_pts, mesh_nrm, mesh_faces, mesh_col, vertex_ele,
                        edge_pts, edge_col, ref_pts, RGBA{Float32}(0, 0, 0, 1),
                        label_pos, label_str, label_ele,
                        ele_entrance, ele_center, ele_exit)
end

@inline _vec3(a::Char, b::Char, c::Char, v::Vec3d) =
  Vec3f(_axis(a, v), _axis(b, v), _axis(c, v))

# Elements that share a piece of the machine -- a pickup with the two correctors
# wound around it -- share a label anchor, and the floor plan deals with that by
# stacking their labels outward along the ray they have in common (see
# `_stack_labels!`). The same problem exists here, but the 2D answer does not
# transfer: these labels billboard, so there is no fixed ray to stack along and
# no reading direction to measure a stack in characters of.
#
# They are stacked *upward* instead, along the element's own vertical, by a step
# set by the element's size. That direction survives the camera moving: it stays
# a separation on screen from every angle except looking straight down the
# vertical, and it needs no idea of how wide the text will turn out to be.
#
# The run-detection is the 2D rule: labels arrive in branch order, so a group of
# coincident elements is a run of consecutive entries, each tested against the
# one before it.
function _stack_labels3(base, up, anchor, h, sep)
  pos = copy(base)
  isempty(pos) && return pos
  step = 0.0f0
  for k in 2:length(pos)
    d = hypot((anchor[k] - anchor[k - 1])...)
    step = d < sep * 0.5 * (h[k] + h[k - 1]) ? step + 2.2f0 * h[k] : 0.0f0
    pos[k] = base[k] + step * up[k]
  end
  return pos
end

# Extrude one shape template along an element's reference curve.
#
# `at(f)` gives the centerline point and the local axes there, all in drawn
# coordinates; a profile point `(f, w)` at height `v` sits at
# `at(f).p + w*h*x̂(f) + v*h2*ŷ(f)`. Caps come from the template's convex faces,
# sliced along `f` so a bend's arc is followed rather than chorded; side walls
# come from its outline loops, whose edges are subdivided the same way.
function _solid!(tri!, edge!, at, tpl::ShapeTemplate, L, ang, h, h2, col, i,
                 arc_tol, edges)
  curved = abs(ang) > 1e-9
  nslice = curved ? max(1, ceil(Int, abs(ang) / arc_tol)) : 1

  # A profile point, at height v = +-1 (top or bottom cap).
  function pt(f, w, v)
    p, x̂, ŷ, _ = at(f)
    return p + Float32(w * h) * x̂ + Float32(v * h2) * ŷ
  end

  # ── caps ──
  for face in tpl.faces, k in 0:(nslice - 1)
    poly = nslice == 1 ? face : _clip_f(face, k / nslice, (k + 1) / nslice)
    length(poly) < 3 && continue
    fbar = sum(p[1] for p in poly) / length(poly)
    _, _, ŷ, _ = at(fbar)
    for v in (1.0, -1.0)
      # The cap's outward normal is the local vertical, up on top and down
      # below; the winding is then made to agree so both faces light correctly.
      nrm = Float32(v) * ŷ
      for t in 2:(length(poly) - 1)
        p1 = pt(poly[1]..., v); p2 = pt(poly[t]..., v); p3 = pt(poly[t + 1]..., v)
        dot3(cross3(p2 - p1, p3 - p1), nrm) < 0 && ((p2, p3) = (p3, p2))
        tri!(p1, p2, p3, nrm, col, i)
      end
    end
  end

  # ── side walls ──
  # Raised along every outline edge. The wall's outward direction is the edge
  # normal in the profile plane turned away from the profile's interior, which
  # is taken to be the centroid of the cap faces -- inside every one of the
  # templates, including the two self-intersecting ones, whose crossing point is
  # exactly that centroid.
  cf, cw = _profile_centroid(tpl)
  for loop in tpl.loops
    for k in 1:(length(loop) - 1)
      f0, w0 = loop[k]; f1, w1 = loop[k + 1]
      nsub = curved ? max(1, ceil(Int, abs(ang) * abs(f1 - f0) / arc_tol)) : 1
      for j in 0:(nsub - 1)
        t0 = j / nsub; t1 = (j + 1) / nsub
        fa = f0 + (f1 - f0) * t0; wa = w0 + (w1 - w0) * t0
        fb = f0 + (f1 - f0) * t1; wb = w0 + (w1 - w0) * t1

        # Worked out at physical scale, so a long thin element gets it right,
        # then carried into 3D on the local axes at the segment midpoint.
        dp = (fb - fa) * L; dq = (wb - wa) * h
        np, nq = dq, -dp
        n = hypot(np, nq)
        n <= 1.0e-12 && continue      # degenerate edge (zero-length element end)
        np /= n; nq /= n
        if np * (0.5 * (fa + fb) - cf) * L + nq * (0.5 * (wa + wb) - cw) * h < 0
          np = -np; nq = -nq          # turn it away from the interior
        end
        _, x̂, _, ẑ = at(0.5 * (fa + fb))
        nrm = Float32(nq) * x̂ + Float32(np) * ẑ

        A = pt(fa, wa, -1.0); B = pt(fb, wb, -1.0)
        C = pt(fb, wb, 1.0);  D = pt(fa, wa, 1.0)
        for (q1, q2, q3) in ((A, B, C), (A, C, D))
          if dot3(cross3(q2 - q1, q3 - q1), nrm) < 0
            tri!(q1, q3, q2, nrm, col, i)
          else
            tri!(q1, q2, q3, nrm, col, i)
          end
        end
        edges && (edge!(A, B, col); edge!(D, C, col))
      end
    end
  end

  # ── interior strokes, on both caps ──
  if edges
    for ((f0, w0), (f1, w1)) in tpl.segs
      nsub = curved ? max(1, ceil(Int, abs(ang) * abs(f1 - f0) / arc_tol)) : 1
      for j in 0:(nsub - 1), v in (1.0, -1.0)
        t0 = j / nsub; t1 = (j + 1) / nsub
        edge!(pt(f0 + (f1 - f0) * t0, w0 + (w1 - w0) * t0, v),
              pt(f0 + (f1 - f0) * t1, w0 + (w1 - w0) * t1, v), col)
      end
    end
  end
end

# A `:circle` element in 3D. The 2D drawing puts a disc at the element midpoint
# rather than a shape swept along it, and the solid whose silhouette is that disc
# from every direction -- not only from above -- is a sphere.
function _sphere!(tri!, edge!, center, r, col, i, nseg, edges)
  nlat = max(3, nseg ÷ 2)
  vert(u, v) = (θ = 2π * u / nseg; φ = π * v / nlat;
                Vec3f(sin(φ) * cos(θ), cos(φ), sin(φ) * sin(θ)))
  p(n) = center + Float32(r) * n
  # On a sphere the outward normal is the vertex direction itself; a facet takes
  # the mean of its three, renormalized so flat shading still gets a unit normal.
  unit(n) = (m = hypot(n...); m > 0 ? n / m : Vec3f(0, 1, 0))
  for v in 0:(nlat - 1), u in 0:(nseg - 1)
    n1 = vert(u, v);         n2 = vert(u + 1, v)
    n3 = vert(u + 1, v + 1); n4 = vert(u, v + 1)
    v > 0 && tri!(p(n1), p(n2), p(n3), unit(n1 + n2 + n3), col, i)
    v < nlat - 1 && tri!(p(n1), p(n3), p(n4), unit(n1 + n3 + n4), col, i)
    edges && v > 0 && edge!(p(n1), p(n2), col)
  end
end

"""
    element_outline3(tab, i; view="zxy", pad=1.15, minsize=0.3, arc_tol=0.08)
        -> Vector{Point3f}

Wireframe box hugging element `i`, as points to be drawn as line segments
(consecutive pairs). Used to highlight the selected element: the 3D counterpart
of [`element_outline`](@ref).
"""
function element_outline3(tab::ElementTable, i::Int; view::AbstractString="zxy",
                          pad::Real=1.15, minsize::Real=0.3, arc_tol::Real=0.08)
  a, b, c = view[1], view[2], view[3]
  L = tab.length[i]; ang = tab.angle[i]; tilt = tab.tilt_ref[i]
  e = entrance(tab, i)
  h = Float32(max(minsize, 0.3) * pad)

  nsamp = max(2, ceil(Int, abs(ang) / arc_tol) + 1)
  corner = [Point3f[] for _ in 1:4]
  for k in 0:(nsamp - 1)
    pl = place(e.r, e.W, L, ang, tilt, k / (nsamp - 1))
    p = proj3(a, b, c, pl.r)
    x̂ = _vec3(a, b, c, xaxis(pl.W)); ŷ = _vec3(a, b, c, yaxis(pl.W))
    push!(corner[1], p + h * (x̂ + ŷ)); push!(corner[2], p + h * (x̂ - ŷ))
    push!(corner[3], p + h * (-x̂ - ŷ)); push!(corner[4], p + h * (-x̂ + ŷ))
  end

  pts = Point3f[]
  for cn in corner, k in 1:(nsamp - 1)      # the four long edges
    push!(pts, cn[k], cn[k + 1])
  end
  for k in (1, nsamp), j in 1:4             # the two end rings
    push!(pts, corner[j][k], corner[j % 4 + 1][k])
  end
  return pts
end
