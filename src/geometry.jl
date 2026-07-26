# Turning an ElementTable + ShapeMap into batched 2D draw data.
#
# Everything is accumulated into a handful of flat vectors so the whole lattice
# can be drawn with only a few Makie calls regardless of element count:
#   * element outlines  -> one `lines!`        (NaN-separated closed loops)
#   * "X" strokes etc.  -> one `linesegments!`
#   * reference orbit   -> one `lines!`        (NaN-separated per branch)
#   * labels            -> one `text!`
#
# Each element's shape is defined by a template in normalized coordinates
# `(f, w)` with `f in [0,1]` along the element length and `w in [-1,1]`
# transverse. The template is mapped through the element's (possibly curved)
# centerline, so bends automatically follow their arc.

using GeometryBasics: Point2f
using Colors

const _NAN2 = Point2f(NaN32, NaN32)

# Geometry ready to hand to Makie, plus per-element data for picking/highlight.
struct FloorGeometry
  outline_pts::Vector{Point2f}
  outline_col::Vector{RGBA{Float32}}
  outline_wid::Vector{Float32}

  seg_pts::Vector{Point2f}          # X strokes; consecutive pairs are segments
  seg_col::Vector{RGBA{Float32}}

  ref_pts::Vector{Point2f}          # reference-orbit centerline, per branch
  ref_col::RGBA{Float32}

  label_pos::Vector{Point2f}
  label_str::Vector{String}
  label_ele::Vector{Int}

  ele_center::Vector{Point2f}       # projected element midpoint (for picking)
end

# ── projection ────────────────────────────────────────────────────────────────
# `view` is a two-character plane like "zx": first char -> horizontal axis,
# second -> vertical, each one of 'x','y','z' in the global reference system.

_axis(c::Char, x, y, z) = c == 'x' ? x : c == 'y' ? y : z

# Projected element frame: the forward direction ĝ and the transverse direction
# n̂, for an element whose entrance heading is `th`.
#
# Both are the *branch axes* projected through the view, not one derived from the
# other by rotating 90° in the picture plane: which way that rotation points
# depends on the handedness of the axis pair, so "zx" and "xz" -- the same plane
# seen with the screen axes swapped -- would get opposite transverse directions,
# and a bend would curve the wrong way in one of them.
#
# ĝ is the branch z axis (sinθ, 0, cosθ) and n̂ the branch x axis
# (cosθ, 0, -sinθ). Only the heading θ is used; `phi` and `psi` are ignored, so a
# reference curve pitched or rolled out of the picture plane is drawn in it (see
# the known gap in the README).
@inline function _frame(a::Char, b::Char, th)
  s, c = sin(th), cos(th)
  gx = _axis(a, s, 0.0, c);  gy = _axis(b, s, 0.0, c)
  nx = _axis(a, c, 0.0, -s); ny = _axis(b, c, 0.0, -s)
  gn = hypot(gx, gy); gn > 0 && (gx /= gn; gy /= gn)
  nn = hypot(nx, ny); nn > 0 && (nx /= nn; ny /= nn)
  return gx, gy, nx, ny
end

# ── shape templates ───────────────────────────────────────────────────────────
# Each returns (loops, segs) in normalized (f, w) coordinates. `loops` are closed
# outlines; `segs` are independent strokes drawn inside/around the element.

function _template(shape::Symbol)
  box = [(0.0, -1.0), (1.0, -1.0), (1.0, 1.0), (0.0, 1.0), (0.0, -1.0)]
  X   = [((0.0, -1.0), (1.0, 1.0)), ((0.0, 1.0), (1.0, -1.0))]
  if shape === :box
    return ([box], NTuple{2,NTuple{2,Float64}}[])
  elseif shape === :xbox
    return ([box], X)
  elseif shape === :x
    return (Vector{Tuple{Float64,Float64}}[], X)
  elseif shape === :diamond
    return ([[(0.0, 0.0), (0.5, -1.0), (1.0, 0.0), (0.5, 1.0), (0.0, 0.0)]],
            NTuple{2,NTuple{2,Float64}}[])
  elseif shape === :bow_tie
    return ([[(0.0, -1.0), (1.0, 1.0), (1.0, -1.0), (0.0, 1.0), (0.0, -1.0)]],
            NTuple{2,NTuple{2,Float64}}[])
  elseif shape === :rbow_tie
    return ([[(0.0, -1.0), (1.0, -1.0), (0.0, 1.0), (1.0, 1.0), (0.0, -1.0)]],
            NTuple{2,NTuple{2,Float64}}[])
  elseif shape === :u_triangle
    return ([[(0.0, -1.0), (1.0, -1.0), (0.5, 1.0), (0.0, -1.0)]],
            NTuple{2,NTuple{2,Float64}}[])
  elseif shape === :d_triangle
    return ([[(0.0, 1.0), (1.0, 1.0), (0.5, -1.0), (0.0, 1.0)]],
            NTuple{2,NTuple{2,Float64}}[])
  elseif shape === :r_triangle
    return ([[(0.0, -1.0), (1.0, 0.0), (0.0, 1.0), (0.0, -1.0)]],
            NTuple{2,NTuple{2,Float64}}[])
  elseif shape === :l_triangle
    return ([[(1.0, -1.0), (0.0, 0.0), (1.0, 1.0), (1.0, -1.0)]],
            NTuple{2,NTuple{2,Float64}}[])
  else  # :circle and :none handled specially by the caller
    return (Vector{Tuple{Float64,Float64}}[], NTuple{2,NTuple{2,Float64}}[])
  end
end

# ── centerline of one element ─────────────────────────────────────────────────
# Given the projected entrance point, projected heading unit vectors ĝ (forward)
# and n̂ (transverse), length L and bend angle, return position + transverse unit
# normal at length fraction `f`. Straight elements have a constant normal; bends
# rotate the normal along the arc.

@inline function _centerline_at(P0::Point2f, gx, gy, nx, ny, L, ang, f)
  if abs(ang) < 1e-9
    p = Point2f(P0[1] + f * L * gx, P0[2] + f * L * gy)
    return p, (nx, ny)
  end
  rho = L / ang
  t = f * ang
  st, ct = sin(t), cos(t)
  lo = rho * st            # longitudinal offset along ĝ
  tr = rho * (ct - 1)      # transverse offset along n̂
  p = Point2f(P0[1] + lo * gx + tr * nx, P0[2] + lo * gy + tr * ny)
  # normal rotates by t: n(f) = sin(t)·ĝ + cos(t)·n̂
  nnx = st * gx + ct * nx
  nny = st * gy + ct * ny
  return p, (nnx, nny)
end

# ── main builder ──────────────────────────────────────────────────────────────

"""
    build_geometry(tab::ElementTable, smap::ShapeMap; view="zx",
                   arc_tol=0.08, circle_sides=24) -> FloorGeometry

Compute batched draw data for the whole lattice. `view` selects the projection
plane. `arc_tol` is the max angular step (rad) used to tessellate bend arcs.
"""
function build_geometry(tab::ElementTable, smap::ShapeMap;
                        view::AbstractString="zx", arc_tol::Real=0.08,
                        circle_sides::Int=24)
  a, b = view[1], view[2]

  outline_pts = Point2f[]; outline_col = RGBA{Float32}[]; outline_wid = Float32[]
  seg_pts = Point2f[];     seg_col = RGBA{Float32}[]
  ref_pts = Point2f[]
  label_pos = Point2f[];   label_str = String[]; label_ele = Int[]
  ele_center = Vector{Point2f}(undef, length(tab))

  # Emit a normalized closed loop mapped through element i's centerline. Edges
  # that span a range of `f` on a curved element are subdivided so the outline
  # follows the arc rather than cutting across it as a chord.
  function emit_loop!(loop, P0, gx, gy, nx, ny, L, ang, h, col, wid)
    isempty(loop) && return
    curved = abs(ang) > 1e-9
    push_pt(f, w) = begin
      c, (nnx, nny) = _centerline_at(P0, gx, gy, nx, ny, L, ang, f)
      push!(outline_pts, Point2f(c[1] + w * h * nnx, c[2] + w * h * nny))
      push!(outline_col, col); push!(outline_wid, wid)
    end
    m = length(loop)
    for k in 1:(m - 1)
      f0, w0 = loop[k]; f1, w1 = loop[k + 1]
      nsub = curved ? max(1, ceil(Int, abs(ang) * abs(f1 - f0) / arc_tol)) : 1
      for j in 0:(nsub - 1)          # [f0,f1): the closing vertex is emitted below
        t = j / nsub
        push_pt(f0 + (f1 - f0) * t, w0 + (w1 - w0) * t)
      end
    end
    push_pt(loop[end]...)             # closing vertex (templates repeat the first)
    push!(outline_pts, _NAN2); push!(outline_col, col); push!(outline_wid, wid)
  end

  function emit_seg!(seg, P0, gx, gy, nx, ny, L, ang, h, col)
    for pt in seg
      f, w = pt
      c, (nnx, nny) = _centerline_at(P0, gx, gy, nx, ny, L, ang, f)
      push!(seg_pts, Point2f(c[1] + w * h * nnx, c[2] + w * h * nny))
      push!(seg_col, col)
    end
  end

  prev_branch = 0
  for i in 1:length(tab)
    xi, yi, zi = tab.x[i], tab.y[i], tab.z[i]
    th = tab.theta[i]; L = tab.length[i]; ang = tab.angle[i]

    P0 = Point2f(_axis(a, xi, yi, zi), _axis(b, xi, yi, zi))
    gx, gy, nx, ny = _frame(a, b, th)

    # centerline midpoint (for picking + label anchor)
    cmid, (mnx, mny) = _centerline_at(P0, gx, gy, nx, ny, L, ang, 0.5)
    ele_center[i] = cmid

    # reference orbit: break the polyline between branches
    if tab.branch[i] != prev_branch && !isempty(ref_pts)
      push!(ref_pts, _NAN2)
    end
    prev_branch = tab.branch[i]
    nsamp = max(2, ceil(Int, abs(ang) / arc_tol) + 1)
    for k in 0:(nsamp - 1)
      c, _ = _centerline_at(P0, gx, gy, nx, ny, L, ang, k / (nsamp - 1))
      push!(ref_pts, c)
    end

    spec = mapshape(smap, tab.kind[i], tab.name[i])
    (!spec.draw || spec.shape === :none) && continue
    col = RGBA{Float32}(spec.color)
    h = Float32(spec.size); wid = Float32(spec.line_width)

    if spec.shape === :circle
      loop = Vector{Tuple{Float64,Float64}}(undef, circle_sides + 1)
      for k in 0:circle_sides
        φ = 2π * k / circle_sides
        push!(outline_pts, Point2f(cmid[1] + h * cos(φ), cmid[2] + h * sin(φ)))
        push!(outline_col, col); push!(outline_wid, wid)
      end
      push!(outline_pts, _NAN2); push!(outline_col, col); push!(outline_wid, wid)
    else
      loops, segs = _template(spec.shape)
      for lp in loops
        emit_loop!(lp, P0, gx, gy, nx, ny, L, ang, h, col, wid)
      end
      for sg in segs
        emit_seg!(sg, P0, gx, gy, nx, ny, L, ang, h, col)
      end
    end

    # label, offset transversely just past the shape
    if spec.label !== :none
      off = 1.6f0 * h
      push!(label_pos, Point2f(cmid[1] + off * mnx, cmid[2] + off * mny))
      push!(label_str, spec.label === :s ? string(round(tab.s[i]; digits=3)) : tab.name[i])
      push!(label_ele, i)
    end
  end

  return FloorGeometry(outline_pts, outline_col, outline_wid,
                       seg_pts, seg_col, ref_pts, RGBA{Float32}(0, 0, 0, 1),
                       label_pos, label_str, label_ele, ele_center)
end

"""
    element_outline(tab, i; view="zx", pad=1.0, arc_tol=0.08) -> Vector{Point2f}

Closed rectangular outline hugging element `i` (its box template inflated by
`pad`), used to highlight the selected element. Returns points in projected
coordinates.
"""
function element_outline(tab::ElementTable, i::Int; view::AbstractString="zx",
                         pad::Real=1.15, arc_tol::Real=0.08, minsize::Real=0.3)
  a, b = view[1], view[2]
  xi, yi, zi = tab.x[i], tab.y[i], tab.z[i]
  th = tab.theta[i]; L = tab.length[i]; ang = tab.angle[i]
  P0 = Point2f(_axis(a, xi, yi, zi), _axis(b, xi, yi, zi))
  gx, gy, nx, ny = _frame(a, b, th)
  h = Float32(max(minsize, 0.3) * pad)

  nsamp = max(2, ceil(Int, abs(ang) / arc_tol) + 1)
  top = Point2f[]; bot = Point2f[]
  for k in 0:(nsamp - 1)
    f = k / (nsamp - 1)
    c, (nnx, nny) = _centerline_at(P0, gx, gy, nx, ny, L, ang, f)
    push!(top, Point2f(c[1] + h * nnx, c[2] + h * nny))
    push!(bot, Point2f(c[1] - h * nnx, c[2] - h * nny))
  end
  return vcat(top, reverse(bot), top[1:1])
end
