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
#   * Labels. Text billboards toward the camera, so the floor plan's trick of
#     laying labels across the centerline has no meaning here; coincident ones
#     are stacked vertically instead (see `_stack_labels3`). They are also culled
#     harder than in 2D -- drawn only when the view is zoomed in far enough *and*
#     the element is inside the current limits -- because a 3D view crowds much
#     faster than a top-down one, and showing fewer labels is what keeps it
#     readable.

using GeometryBasics: Point2f, Point3f, Vec2f, Vec3f, widths
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
                label_min_px=14, size=(1500, 850), edges=true) -> FloorPlot3

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
  * `label_min_px` — labels appear only when an element's on-screen size reaches
    this many pixels. Higher than the 2D default, because a 3D view crowds
    sooner.
  * `label_sep` — how close two labels must be to count as colliding, in units of
    the elements' half-height. Colliding labels are stacked vertically; `0` turns
    that off.
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
                     label_min_px::Real=14, label_sep::Real=1.0,
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
                     label_min_px::Real=14, label_sep::Real=1.0,
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

  # Labels. Level-of-detail as in 2D -- hide them until elements are visually
  # resolvable -- and on top of that, cull the ones outside the current view box:
  # zooming into part of a machine should not leave the labels of the rest of it
  # floating over the picture. Culling is done by moving a label to NaN rather
  # than by a `visible` flag, which is per-plot rather than per-label.
  if !isempty(geom.label_pos)
    pos = lift(ax.finallimits, ax.scene.viewport) do lims, vp
      mpp = Float32(maximum(lims.widths) / max(1, vp.widths[1]))  # metres/pixel
      if 1.0f0 / mpp < label_min_px
        return fill(_NAN3, length(geom.label_pos))
      end
      lo = minimum(lims); hi = maximum(lims)
      return map(geom.label_pos) do p
        all(lo .<= p .<= hi) ? p : _NAN3
      end
    end
    text!(ax, pos; text=geom.label_str, fontsize=11, align=(:left, :center),
          offset=Vec2f(4, 0), color=:gray25, depth_shift=-0.01f0)
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
