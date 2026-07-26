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
# transverse. The template is swept along the element's reference curve using the
# placement math in frame.jl, so bends automatically follow their arc -- and a
# lattice that leaves the horizontal plane is drawn as its true projection onto
# the chosen plane rather than as though it stayed in it.

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
  label_rot::Vector{Float32}                  # radians, in [-pi/2, pi/2]
  label_align::Vector{Tuple{Symbol,Symbol}}
  label_stack::Vector{Float32}                # collision shift, in font sizes
  label_ele::Vector{Int}

  ele_center::Vector{Point2f}       # projected element midpoint (for picking)
end

# ── shape templates ───────────────────────────────────────────────────────────
# A template is a profile in normalized (f, w) coordinates: `f in [0,1]` along
# the element length, `w in [-1,1]` transverse. It carries three things:
#
#   * `loops` -- closed outlines, what the 2D drawing strokes and what the 3D
#     drawing raises side walls along;
#   * `segs`  -- independent interior strokes (the X of an `xbox`);
#   * `faces` -- *convex* polygons tiling the same profile, which the 3D drawing
#     uses to cap the extruded solid.
#
# `faces` are given per shape rather than triangulated on the fly because two of
# the profiles (`bow_tie` and `rbow_tie`) are self-intersecting, so no fan from a
# centroid tiles them: each is two triangles meeting at the crossing point, which
# is what is written out below.

struct ShapeTemplate
  loops::Vector{Vector{Tuple{Float64,Float64}}}
  segs::Vector{NTuple{2,NTuple{2,Float64}}}
  faces::Vector{Vector{Tuple{Float64,Float64}}}
end

const _NOSEG = NTuple{2,NTuple{2,Float64}}[]
const _NOFACE = Vector{Tuple{Float64,Float64}}[]

function _template(shape::Symbol)
  box = [(0.0, -1.0), (1.0, -1.0), (1.0, 1.0), (0.0, 1.0), (0.0, -1.0)]
  boxface = [[(0.0, -1.0), (1.0, -1.0), (1.0, 1.0), (0.0, 1.0)]]
  X = [((0.0, -1.0), (1.0, 1.0)), ((0.0, 1.0), (1.0, -1.0))]
  if shape === :box
    return ShapeTemplate([box], _NOSEG, boxface)
  elseif shape === :xbox
    return ShapeTemplate([box], X, boxface)
  elseif shape === :x
    return ShapeTemplate(Vector{Tuple{Float64,Float64}}[], X, _NOFACE)
  elseif shape === :diamond
    return ShapeTemplate([[(0.0, 0.0), (0.5, -1.0), (1.0, 0.0), (0.5, 1.0), (0.0, 0.0)]],
                         _NOSEG,
                         [[(0.0, 0.0), (0.5, -1.0), (1.0, 0.0), (0.5, 1.0)]])
  elseif shape === :bow_tie
    # The two diagonals cross at (0.5, 0), leaving a triangle on each end.
    return ShapeTemplate([[(0.0, -1.0), (1.0, 1.0), (1.0, -1.0), (0.0, 1.0), (0.0, -1.0)]],
                         _NOSEG,
                         [[(0.0, -1.0), (0.5, 0.0), (0.0, 1.0)],
                          [(1.0, -1.0), (0.5, 0.0), (1.0, 1.0)]])
  elseif shape === :rbow_tie
    # Same crossing, but the untouched pair of edges is top and bottom.
    return ShapeTemplate([[(0.0, -1.0), (1.0, -1.0), (0.0, 1.0), (1.0, 1.0), (0.0, -1.0)]],
                         _NOSEG,
                         [[(0.0, -1.0), (1.0, -1.0), (0.5, 0.0)],
                          [(0.0, 1.0), (1.0, 1.0), (0.5, 0.0)]])
  elseif shape === :u_triangle
    return ShapeTemplate([[(0.0, -1.0), (1.0, -1.0), (0.5, 1.0), (0.0, -1.0)]],
                         _NOSEG, [[(0.0, -1.0), (1.0, -1.0), (0.5, 1.0)]])
  elseif shape === :d_triangle
    return ShapeTemplate([[(0.0, 1.0), (1.0, 1.0), (0.5, -1.0), (0.0, 1.0)]],
                         _NOSEG, [[(0.0, 1.0), (1.0, 1.0), (0.5, -1.0)]])
  elseif shape === :r_triangle
    return ShapeTemplate([[(0.0, -1.0), (1.0, 0.0), (0.0, 1.0), (0.0, -1.0)]],
                         _NOSEG, [[(0.0, -1.0), (1.0, 0.0), (0.0, 1.0)]])
  elseif shape === :l_triangle
    return ShapeTemplate([[(1.0, -1.0), (0.0, 0.0), (1.0, 1.0), (1.0, -1.0)]],
                         _NOSEG, [[(1.0, -1.0), (0.0, 0.0), (1.0, 1.0)]])
  else  # :circle and :none handled specially by the caller
    return ShapeTemplate(Vector{Tuple{Float64,Float64}}[], _NOSEG, _NOFACE)
  end
end

# ── projected centerline of one element ───────────────────────────────────────
# Position and transverse direction at length fraction `f`, both projected onto
# the view plane. The transverse direction is deliberately *not* renormalized
# after projection: an element tilted out of the picture plane must come out
# foreshortened, which is exactly the length the projection of its own frame has.

@inline function _at2(a::Char, b::Char, e::Placement, L, ang, tilt, f)
  pl = place(e.r, e.W, L, ang, tilt, f)
  return proj2(a, b, pl.r), proj2(a, b, xaxis(pl.W))
end

# ── label collisions ──────────────────────────────────────────────────────────
# Elements that occupy the same piece of the machine -- a pickup and the two
# correctors wound around it, say -- have their midpoints, and so their label
# anchors, on top of each other. Perpendicular labels do not help there: they
# are all on the same ray and print over one another. Such labels get pushed
# out along that ray, one after the next, so they read as a stack.
#
# The shift is measured in font sizes rather than meters because that is what
# the text is: a fixed size on screen, unaffected by zoom. The caller turns it
# into a pixel offset (see `render.jl`).

const _CHAR_EM = 0.62f0    # mean glyph advance of the label font, in font sizes
const _GAP_EM = 0.9f0      # blank left between one stacked label and the next

# Labels arrive in element order, which is branch order, so a run of colliding
# labels is a run of consecutive entries. Each is tested against the one before
# it -- chaining, so three coincident elements form one stack of three -- and
# the first of the run keeps the anchor, leaving it closest to the centerline.
function _stack_labels!(stack, pos, str, hs, sep)
  isempty(stack) && return stack
  stack[1] = 0
  for k in 2:length(stack)
    d = hypot(pos[k][1] - pos[k - 1][1], pos[k][2] - pos[k - 1][2])
    stack[k] = if d < sep * 0.5 * (hs[k] + hs[k - 1])
      stack[k - 1] + _CHAR_EM * length(str[k - 1]) + _GAP_EM
    else
      0.0f0
    end
  end
  return stack
end

# ── main builder ──────────────────────────────────────────────────────────────

"""
    build_geometry(tab::ElementTable, smap::ShapeMap; view="zx",
                   arc_tol=0.08, circle_sides=24, label_sep=1.0) -> FloorGeometry

Compute batched draw data for the whole lattice. `view` selects the projection
plane. `arc_tol` is the max angular step (rad) used to tessellate bend arcs.

`label_sep` sets when two labels count as colliding: they do when their anchors
are closer than `label_sep` times the mean of the two elements' half-heights,
which are the scale the drawing (and so the readable font size) is built around.
Colliding labels are stacked outward along the centerline normal. Raise it to
stack more eagerly, or set it to 0 to switch the stacking off.
"""
function build_geometry(tab::ElementTable, smap::ShapeMap;
                        view::AbstractString="zx", arc_tol::Real=0.08,
                        circle_sides::Int=24, label_sep::Real=1.0)
  a, b = view[1], view[2]

  outline_pts = Point2f[]; outline_col = RGBA{Float32}[]; outline_wid = Float32[]
  seg_pts = Point2f[];     seg_col = RGBA{Float32}[]
  ref_pts = Point2f[]
  label_pos = Point2f[];   label_str = String[]; label_ele = Int[]
  label_rot = Float32[];   label_align = Tuple{Symbol,Symbol}[]
  label_h = Float32[]      # each label's element half-height; collision scale
  ele_center = Vector{Point2f}(undef, length(tab))

  # Emit a normalized closed loop swept along element i's reference curve. Edges
  # that span a range of `f` on a curved element are subdivided so the outline
  # follows the arc rather than cutting across it as a chord.
  function emit_loop!(loop, e, L, ang, tilt, h, col, wid)
    isempty(loop) && return
    curved = abs(ang) > 1e-9
    push_pt(f, w) = begin
      c, n = _at2(a, b, e, L, ang, tilt, f)
      push!(outline_pts, Point2f(c[1] + w * h * n[1], c[2] + w * h * n[2]))
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

  function emit_seg!(seg, e, L, ang, tilt, h, col)
    for (f, w) in seg
      c, n = _at2(a, b, e, L, ang, tilt, f)
      push!(seg_pts, Point2f(c[1] + w * h * n[1], c[2] + w * h * n[2]))
      push!(seg_col, col)
    end
  end

  prev_branch = 0
  for i in 1:length(tab)
    L = tab.length[i]; ang = tab.angle[i]; tilt = tab.tilt_ref[i]
    e = entrance(tab, i)

    # Centerline midpoint (for picking + label anchor), and the *unit* projected
    # transverse direction there, which is what a label is laid out along.
    mid = place(e.r, e.W, L, ang, tilt, 0.5)
    cmid = proj2(a, b, mid.r)
    _, _, mnx, mny = proj_axes(a, b, mid.W)
    ele_center[i] = cmid

    # reference orbit: break the polyline between branches
    if tab.branch[i] != prev_branch && !isempty(ref_pts)
      push!(ref_pts, _NAN2)
    end
    prev_branch = tab.branch[i]
    nsamp = max(2, ceil(Int, abs(ang) / arc_tol) + 1)
    for k in 0:(nsamp - 1)
      c, _ = _at2(a, b, e, L, ang, tilt, k / (nsamp - 1))
      push!(ref_pts, c)
    end

    spec = mapshape(smap, tab.kind[i], tab.name[i])
    (!spec.draw || spec.shape === :none) && continue
    col = RGBA{Float32}(spec.color)
    h = Float32(spec.size); wid = Float32(spec.line_width)

    if spec.shape === :circle
      for k in 0:circle_sides
        φ = 2π * k / circle_sides
        push!(outline_pts, Point2f(cmid[1] + h * cos(φ), cmid[2] + h * sin(φ)))
        push!(outline_col, col); push!(outline_wid, wid)
      end
      push!(outline_pts, _NAN2); push!(outline_col, col); push!(outline_wid, wid)
    else
      tpl = _template(spec.shape)
      for lp in tpl.loops
        emit_loop!(lp, e, L, ang, tilt, h, col, wid)
      end
      for sg in tpl.segs
        emit_seg!(sg, e, L, ang, tilt, h, col)
      end
    end

    # Label, offset transversely just past the shape and set running along the
    # transverse normal, i.e. perpendicular to the centerline. Neighbouring
    # elements are strung out *along* the line, so labels laid across it stay
    # clear of each other where horizontal ones would collide.
    if spec.label !== :none
      off = 1.6f0 * h
      push!(label_pos, Point2f(cmid[1] + off * mnx, cmid[2] + off * mny))
      push!(label_str, spec.label === :s ? string(round(tab.s[i]; digits=3)) : tab.name[i])
      # Reading direction is n̂ folded into the right half-plane, so text is
      # never upside down and the angle lands in [-pi/2, pi/2]. Where that
      # folding reverses n̂ the anchor flips with it, so the label always grows
      # away from the centerline rather than back across the machine.
      flip = mnx < 0 || (mnx == 0 && mny < 0)
      rx, ry = flip ? (-mnx, -mny) : (mnx, mny)
      push!(label_rot, Float32(atan(ry, rx)))
      push!(label_align, flip ? (:right, :center) : (:left, :center))
      push!(label_h, h)
      push!(label_ele, i)
    end
  end

  label_stack = Vector{Float32}(undef, length(label_pos))
  _stack_labels!(label_stack, label_pos, label_str, label_h, label_sep)

  return FloorGeometry(outline_pts, outline_col, outline_wid,
                       seg_pts, seg_col, ref_pts, RGBA{Float32}(0, 0, 0, 1),
                       label_pos, label_str, label_rot, label_align, label_stack,
                       label_ele, ele_center)
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
  L = tab.length[i]; ang = tab.angle[i]; tilt = tab.tilt_ref[i]
  e = entrance(tab, i)
  h = Float32(max(minsize, 0.3) * pad)

  nsamp = max(2, ceil(Int, abs(ang) / arc_tol) + 1)
  top = Point2f[]; bot = Point2f[]
  for k in 0:(nsamp - 1)
    c, n = _at2(a, b, e, L, ang, tilt, k / (nsamp - 1))
    push!(top, Point2f(c[1] + h * n[1], c[2] + h * n[2]))
    push!(bot, Point2f(c[1] - h * n[1], c[2] - h * n[2]))
  end
  return vcat(top, reverse(bot), top[1:1])
end
