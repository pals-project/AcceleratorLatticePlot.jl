# Makie rendering and interaction for 3D drawings.
#
# The counterpart of render.jl. The whole lattice is three plot objects -- one
# mesh, one linesegments, one lines -- for the same reason the floor plan is
# four: a plot object per element does not survive a machine with tens of
# thousands of them.
#
# Two things do not carry over from 2D and are done differently here:
#
#   * Picking. The floor plan finds the element nearest the cursor in data
#     coordinates, which a 3D axis has no equivalent of: a cursor is a ray, not a
#     point. `_pick3` below asks the backend to pick first -- exact, and it
#     respects occlusion -- and falls back to measuring the ray against element
#     centerlines, which needs no GPU and so also works in CairoMakie and in the
#     tests.
#
#   * Labels. Text billboards toward the camera: it is drawn at a fixed size in
#     pixels, facing the viewer, wherever the camera happens to be. So how much
#     room a label needs, and which way "clear of the element" or "up" points for
#     it, are pixel quantities that only exist once there is a camera -- none of
#     them can be settled in meters while building the geometry, the way the 2D
#     drawing can settle them in a plane it knows will not rotate. The geometry
#     therefore fixes only the anchor, and `_layout_labels3` below places, spaces
#     and thins the labels in screen space on every camera change.

using GeometryBasics: Point2f, Point3f, Vec2f, Vec4f, widths
using Colors
using Makie
import GeometryBasics
import PALSJulia as pj

"""
Handle to a live 3D drawing: the figure, its axis, and the source data. Mirrors
[`FloorPlot`](@ref) field for field, so anything scripted against one -- setting
`selected`, saving `figure`, adding overlays -- works on the other.
"""
struct FloorPlot3
  figure::Figure
  axis::Axis3
  table::ElementTable
  geometry::FloorGeometry3
  view::String
  selected::Observable{Int}
end

Base.display(fp::FloorPlot3) = display(fp.figure)

"""
    floor_plot3(lat; view="zxy", shapes=ShapeMap(), title="Floor Plan (3D)",
                label_min_px=4, size=(1500, 850), edges=true) -> FloorPlot3

Build a 3D drawing of the expanded lattice `lat`: every element a solid, placed
and oriented from the floor coordinates the expander computed. Returns a
`FloorPlot3`, which like [`FloorPlot`](@ref) is just a handle -- `using GLMakie`
and `display(fp)` opens an interactive window, `using CairoMakie` and
`save("floor.png", fp.figure)` writes a file.

GLMakie is the backend this is meant for. CairoMakie will render it, but it sorts
whole primitives rather than pixels, so a large mesh comes out with sorting
artifacts; it is fine for a quick still and for the tests, not for looking at a
machine.

Controls (Makie's default `Axis3` bindings, plus selection): **left-drag**
orbits, **right-drag** pans, **scroll** zooms, **ctrl-left-click** resets the
view, and a plain **left-click** selects an element, listing its parameters in
the side panel.

Keyword arguments:
  * `view` — a permutation of `"xyz"`, one global axis per drawn axis, the last
    being the drawn vertical. The default `"zxy"` stands the machine up the way
    the floor plan lies flat: look straight down and you see `floor_plot`'s
    default view.
  * `shapes` — a `ShapeMap`. The same one the floor plan uses; its `size2`
    controls the vertical half-height of the extruded solids.
  * `label_min_px` — the smallest an element may appear on screen and still be
    labelled, measured across its own projected silhouette. This only skips
    elements too small to see; what decides how many labels a crowded view
    carries is the collision handling below, not this.
  * `label_sep` — padding left around a label, in units of the font size, when
    deciding whether it collides with one already placed. Labels are placed
    nearest-camera first; one whose space is taken is bumped up a row and given a
    leader line back to its element, or dropped if there is still no room after a
    few rows. So a zoomed-out view shows the labels that fit rather than all of
    them on top of each other, and the rest appear as you zoom in. `0` disables
    this, leaving labels free to overlap.
  * `edges` — draw element outlines over the shading.

A machine that lies in a plane is a hundred metres of one-metre magnets, and with
the honest `aspect = :data` its 3D box comes out a pancake with most of the frame
empty. Scroll to zoom and drag to orbit, or set `fp.axis.aspect = (1, 1, 0.4)` to
stretch the vertical if what you are after is a diagram rather than the shape of
the real thing.
"""
function floor_plot3(lat::pj.Lattices; view::AbstractString="zxy",
                     shapes::ShapeMap=ShapeMap(),
                     title::AbstractString="Floor Plan (3D)",
                     label_min_px::Real=4, label_sep::Real=1.0,
                     size=(1500, 850), edges::Bool=true, arc_tol::Real=0.08)
  return floor_plot3(element_table(lat); view=view, shapes=shapes, title=title,
                     label_min_px=label_min_px, label_sep=label_sep, size=size,
                     edges=edges, arc_tol=arc_tol)
end

"""
    floor_plot3(tab::ElementTable; kwargs...) -> FloorPlot3

Draw from an already-extracted `ElementTable` rather than extracting one from a
lattice.
"""
function floor_plot3(tab::ElementTable; view::AbstractString="zxy",
                     shapes::ShapeMap=ShapeMap(),
                     title::AbstractString="Floor Plan (3D)",
                     label_min_px::Real=4, label_sep::Real=1.0,
                     size=(1500, 850), edges::Bool=true, arc_tol::Real=0.08)
  geom = build_geometry3(tab, shapes; view=view, edges=edges, arc_tol=arc_tol,
                         label_sep=label_sep)
  return _floor_plot3(tab, geom; view=view, title=title,
                      label_min_px=label_min_px, size=size)
end

function _floor_plot3(tab::ElementTable, geom::FloorGeometry3; view, title,
                      label_min_px, size)
  fig = Figure(; size=size)
  # `aspect = :data` is not optional: a machine is a hundred metres long and a
  # metre wide, and Axis3's default normalizes each axis to the same box, which
  # would inflate the transverse directions by that whole ratio.
  ax = Axis3(fig[1, 1]; aspect=:data, title=title,
             xlabel=_axis_label(view[1]), ylabel=_axis_label(view[2]),
             zlabel=_axis_label(view[3]))

  # Element solids: one mesh for the whole lattice. Vertices are not shared
  # between faces and each carries its face's normal, so this shades flat.
  if !isempty(geom.mesh_faces)
    mesh!(ax, GeometryBasics.Mesh(geom.mesh_pts, geom.mesh_faces;
                                  normal=geom.mesh_nrm);
          color=geom.mesh_col, transparency=false)
  end

  # Outlines over the shading, and the reference curve through every element.
  isempty(geom.edge_pts) ||
    linesegments!(ax, geom.edge_pts; color=geom.edge_col, linewidth=1)
  isempty(geom.ref_pts) ||
    lines!(ax, geom.ref_pts; color=geom.ref_col, linewidth=1)

  # Labels, laid out in screen space on every camera change (see
  # `_layout_labels3`). Culling a label is done by moving it to NaN rather than
  # by a `visible` flag, which is per-plot rather than per-label.
  if !isempty(geom.label_pos)
    layout = lift(ax.scene.camera.projectionview, ax.scene.viewport,
                  ax.finallimits) do _, _, _
      _layout_labels3(geom, ax.scene, label_min_px, _LABEL_FS)
    end
    text!(ax, lift(l -> l.pos, layout); text=geom.label_str, fontsize=_LABEL_FS,
          align=lift(l -> l.align, layout), offset=lift(l -> l.offset, layout),
          color=:gray25, depth_shift=-0.01f0)

    # Leader lines, tying a label that had to be bumped back to the element it
    # names. One end of each is a position in pixels, so they are drawn in
    # `space = :pixel` -- which on this scene is the same coordinate space
    # `Makie.project` reports, so the endpoints land where the layout put them.
    linesegments!(ax.scene, lift(l -> l.leader, layout); space=:pixel,
                  color=(:gray50, 0.8), linewidth=0.8)
  end

  # Selection highlight: a wireframe cage around the element, drawn with
  # `overdraw` so it stays visible when the element is behind something else --
  # the whole point of highlighting it is to find it.
  selected = Observable(0)
  hi_pts = lift(selected) do i
    i == 0 ? Point3f[] : element_outline3(tab, i; view=view)
  end
  linesegments!(ax, hi_pts; color=:red, linewidth=3, overdraw=true)

  info = _param_panel!(fig)

  # Click-to-select. As in 2D this goes through the interaction system, so
  # Makie's mouse state machine has already told a click from the drag that
  # orbits the camera, and ctrl-left-click is left to Makie's own view reset.
  register_interaction!(ax, :selectelement) do event::MouseEvent, axis
    event.type === MouseEventTypes.leftclick || return Consume(false)
    ispressed(axis.scene, Keyboard.left_control) && return Consume(false)
    i = _pick3(geom, axis, event.px)
    i == 0 && return Consume(false)
    selected[] = i
    info[] = _info_text(tab, i)
    return Consume(true)
  end

  return FloorPlot3(fig, ax, tab, geom, String(view), selected)
end

# ── label layout ──────────────────────────────────────────────────────────────

const _LABEL_FS = 11          # label font size, points
const _GAP_PX = 5.0f0         # readable gap between an element and its label
const _ROW_TRIES = 4          # how far up a colliding label may be bumped

# Depth of a data point in the axis's clip space. Only the ordering is used --
# smaller is nearer the camera -- so the exact convention does not matter, just
# that it is monotonic. `Makie.project` computes this on the way to a screen
# position and then discards it, and the part it discards is short enough to redo.
@inline function _ndc_depth(model, pv, p)
  c = pv * (model * Vec4f(p[1], p[2], p[3], 1))
  return Float32(c[4] == 0 ? Inf32 : c[3] / c[4])
end

# Grid cells a pixel-space box covers, clipped to the grid rather than clamped,
# so a label hanging off the edge occupies only the part that is on screen.
@inline function _cells(bx0, by0, bx1, by1, cell, nx, ny)
  i0 = max(1, floor(Int, bx0 / cell) + 1); i1 = min(nx, floor(Int, bx1 / cell) + 1)
  j0 = max(1, floor(Int, by0 / cell) + 1); j1 = min(ny, floor(Int, by1 / cell) + 1)
  return i0:i1, j0:j1
end

"""
    _layout_labels3(geom, scene, label_min_px, fontsize) -> (; pos, offset, align, leader)

Place the labels for one camera. Returns, per label, the data-space position to
draw it at (`_NAN3` for one that is not drawn), the pixel `offset` from that
position, the text alignment, and pixel-space leader-line endpoints for those
that had to be moved far enough to need one.

Three things happen here that cannot happen while the geometry is being built:

  * **Level of detail, per element.** An element is labelled when *it* is big
    enough on screen, measured across its own projected silhouette. The threshold
    it replaced was one number for the whole view, so on a long machine nothing
    was labelled until it was well zoomed in, and then everything was at once.

  * **Placing the label clear of its element.** How far out that is cannot be
    fixed in meters: what a label has to clear is the element's *silhouette*,
    and a box a meter wide covers anything from two pixels to half the window
    depending on where the camera is and which way the box is turned. So the
    element's three half-extents -- the two transverse probes from the geometry,
    and half the centerline from entrance to exit -- are projected and measured
    along the direction the label goes in, and their sum is how far it is
    pushed. That direction is also what the alignment follows, so the text grows
    away from the centerline rather than back across it.

  * **Collisions.** What overlaps on screen is whatever *projects* close
    together, which is not the same set as what is close together in meters --
    looking down a beamline, elements metres apart land a few pixels apart, and
    a label separation measured in meters shrinks to nothing at the same time.
    So labels are placed greedily against a pixel occupancy grid: bumped up a
    row if the space is taken, dropped if it is still taken after `_ROW_TRIES`.

Labels are placed nearest-camera first. When two want the same pixels the one on
the element you are looking at should win, and the one behind it is the one that
can be dropped without the picture losing anything -- which is also why the text
can keep drawing over the shading rather than depth-testing against it: what
survives this pass is, by construction, not sitting on top of a nearer label.
"""
function _layout_labels3(geom::FloorGeometry3, scene, label_min_px::Real,
                         fontsize::Real)
  n = length(geom.label_pos)
  pos = fill(_NAN3, n)
  offset = fill(Vec2f(0, 0), n)
  align = fill((:left, :center), n)
  leader = Point2f[]
  out = (; pos, offset, align, leader)
  n == 0 && return out

  vp = Makie.viewport(scene)[]
  W, H = widths(vp)
  (W > 0 && H > 0) || return out
  model = scene.transformation.model[]
  pv = scene.camera.projectionview[]

  # Pass 1: project, and drop what cannot be labelled at all.
  scr = Vector{Point2f}(undef, n)     # the anchor, on screen
  dirv = Vector{Vec2f}(undef, n)      # unit direction the label goes out in
  clear = Vector{Float32}(undef, n)   # how far out that has to be, in pixels
  depth = Vector{Float32}(undef, n)
  ord = Int[]
  for k in 1:n
    i = geom.label_ele[k]
    S = Makie.project(scene, geom.label_pos[k])
    (isfinite(S[1]) && isfinite(S[2])) || continue
    (0 <= S[1] <= W && 0 <= S[2] <= H) || continue
    M = Makie.project(scene, geom.ele_center[i])
    A = Makie.project(scene, geom.ele_entrance[i])
    B = Makie.project(scene, geom.ele_exit[i])
    P = Makie.project(scene, geom.label_probe[k])

    # The element's three projected half-extents. Measuring all three rather
    # than the centerline alone is what lets a zero-length element -- a pickup,
    # a marker, most of what a user actually wants named -- be labelled at all.
    d1 = Vec2f(S[1] - M[1], S[2] - M[2])
    d2 = Vec2f(P[1] - M[1], P[2] - M[2])
    d3 = Vec2f(0.5f0 * (B[1] - A[1]), 0.5f0 * (B[2] - A[2]))

    # Level of detail: is the element visible at all? This is only a cheap way
    # of skipping the ones that are not; the greedy pass below is what actually
    # decides how many labels a view can carry.
    2 * max(hypot(d1[1], d1[2]), hypot(d2[1], d2[2]), hypot(d3[1], d3[2])) <
      label_min_px && continue

    # Which way is "away from the centerline" on screen. It collapses when the
    # camera looks straight down the element's own transverse axis, and then any
    # fixed direction will do as well as another.
    m = hypot(d1[1], d1[2])
    u = m < 1 ? Vec2f(1, 0) : Vec2f(d1[1] / m, d1[2] / m)

    # How far out the label has to sit: the silhouette's reach along `u`, being
    # the three half-extents each measured along `u` and summed. An axis that
    # projects to nothing contributes nothing, which is exactly right.
    reach = abs(d1[1] * u[1] + d1[2] * u[2]) +
            abs(d2[1] * u[1] + d2[2] * u[2]) +
            abs(d3[1] * u[1] + d3[2] * u[2])

    scr[k] = Point2f(S)
    dirv[k] = u
    # Offsets are measured from the anchor, which is already `m` pixels out.
    clear[k] = reach - m + _GAP_PX
    depth[k] = _ndc_depth(model, pv, geom.label_pos[k])
    push!(ord, k)
  end
  isempty(ord) && return out
  sort!(ord; by = k -> depth[k])

  # Pass 2: place greedily against a pixel occupancy grid.
  fs = Float32(fontsize)
  # Cells well under a line of text: the grid rounds every box outward to a cell
  # boundary, so a coarse one would have neighbouring labels colliding over
  # pixels neither of them actually covers.
  cell = max(4.0f0, 0.5f0 * fs)
  nx = max(1, ceil(Int, W / cell)); ny = max(1, ceil(Int, H / cell))
  taken = falses(nx, ny)
  pad = Float32(geom.label_sep) * fs * 0.5f0
  th = 1.2f0 * fs
  # A whole padded box, so a bumped label clears the one it was bumped over
  # rather than landing back on it -- plus a cell, because two boxes that merely
  # touch still share the grid cell their common edge falls in, and would read as
  # a collision for as long as they abut exactly.
  rowstep = th + 2 * pad + cell
  tries = geom.label_sep <= 0 ? 1 : _ROW_TRIES

  for k in ord
    tw = _CHAR_EM * fs * length(geom.label_str[k])
    for r in 0:(tries - 1)
      off = Vec2f(dirv[k][1] * clear[k], dirv[k][2] * clear[k] + r * rowstep)
      x0 = scr[k][1] + off[1]; y0 = scr[k][2] + off[2]
      left = dirv[k][1] >= 0
      bx0 = (left ? x0 : x0 - tw) - pad
      bx1 = (left ? x0 + tw : x0) + pad
      by0 = y0 - 0.5f0 * th - pad; by1 = y0 + 0.5f0 * th + pad
      xs, ys = _cells(bx0, by0, bx1, by1, cell, nx, ny)
      any(@view taken[xs, ys]) && continue
      taken[xs, ys] .= true
      pos[k] = geom.label_pos[k]
      offset[k] = off
      align[k] = left ? (:left, :center) : (:right, :center)
      # A label bumped clear of its neighbours is no longer obviously attached to
      # anything, so it gets a leader back to the element it names.
      r > 0 && push!(leader, scr[k], Point2f(x0, y0))
      break
    end
  end
  return out
end

# Which element is under the cursor. Returns 0 for a click on nothing.
#
# The backend is asked first. A GPU pick is exact and, more to the point, knows
# what is in front of what -- which in 3D is most of the problem -- and because
# the mesh's vertices are per-element (`vertex_ele`) its answer maps straight back
# to an element with no search at all.
#
# Backends that cannot pick fall through to a screen-space search: project each
# element's centerline and take the one passing nearest the click. That happens
# in pixels rather than by unprojecting the cursor into a ray, because `project`
# already accounts for everything the axis does to the data on the way to the
# screen (its own scaling, the Float32 conversion, the model matrix) and
# reproducing all of that to cast a ray the other way would only be a way of
# getting it subtly wrong. The trade is that this fallback has no notion of
# depth: where two elements overlap on screen it takes the nearer one in the
# picture, not the nearer one to the camera.
function _pick3(geom::FloorGeometry3, ax::Axis3, px)
  # `px` arrives relative to the scene's origin; `pick` works in the figure's.
  vp = Makie.viewport(ax.scene)[]
  plt, idx = try
    Makie.pick(ax.scene, px .+ minimum(vp))
  catch
    (nothing, 0)
  end
  if plt isa Makie.Mesh && 1 <= idx <= length(geom.vertex_ele)
    return Int(geom.vertex_ele[idx])
  end

  # As generous as the floor plan's own tolerance: a click near the machine
  # should hit something rather than nothing.
  rmax = 0.1 * maximum(widths(vp))
  best = 0; bestd = Float64(rmax)^2
  q = Point2f(px)
  @inbounds for i in eachindex(geom.ele_center)
    A = Makie.project(ax.scene, geom.ele_entrance[i])
    M = Makie.project(ax.scene, geom.ele_center[i])
    B = Makie.project(ax.scene, geom.ele_exit[i])
    d = min(_seg_dist2(q, A, M), _seg_dist2(q, M, B))
    if d < bestd
      bestd = d; best = i
    end
  end
  return best
end

# Squared distance from point `p` to the segment `a`--`b`, in the plane.
function _seg_dist2(p, a, b)
  dx = b[1] - a[1]; dy = b[2] - a[2]
  len2 = dx * dx + dy * dy
  t = len2 <= 0 ? 0.0 : clamp(((p[1] - a[1]) * dx + (p[2] - a[2]) * dy) / len2, 0.0, 1.0)
  ex = a[1] + t * dx - p[1]; ey = a[2] + t * dy - p[2]
  return Float64(ex * ex + ey * ey)
end
