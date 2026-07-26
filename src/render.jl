# Makie rendering and interaction for floor-plan drawings.
#
# The whole lattice is drawn with a few batched calls (see geometry.jl). Pan and
# zoom come from Makie's default Axis interactions. Clicking an element selects
# it: the element is highlighted and its full parameter tree is listed in a side
# panel. Selection uses a nearest-centroid search in data coordinates rather than
# GPU picking, which keeps it independent of how the shapes are drawn.
#
# This builds on `Makie` itself rather than on a particular backend: every call
# here is Makie core, and a backend is needed only to put the figure somewhere --
# GLMakie for an interactive window, CairoMakie to rasterize it to a file. So the
# caller loads whichever they want and PALSPlot has no opinion, which is also
# what lets the tests run on a machine with no OpenGL at all.

using GeometryBasics: Point2f
using Colors
using Makie
import PALSJulia as pj

"Handle to a live floor-plan window: the figure, its axis, and the source data."
struct FloorPlot
  figure::Figure
  axis::Axis
  table::ElementTable
  geometry::FloorGeometry
  selected::Observable{Int}
end

Base.display(fp::FloorPlot) = display(fp.figure)

# Axis labels chosen from the projection, echoing Tao's "SMART LABEL".
_axis_label(c::Char) = string(uppercase(c), " (m)")

# Index of the element whose centroid is nearest `p`, or 0 if none within `rmax`.
function _nearest(centers::Vector{Point2f}, p, rmax::Float64)
  best = 0; bestd = rmax^2
  @inbounds for i in eachindex(centers)
    dx = centers[i][1] - p[1]; dy = centers[i][2] - p[2]
    d = dx * dx + dy * dy
    if d < bestd
      bestd = d; best = i
    end
  end
  return best
end

# Multiline "key: value" dump of an element node, for the parameter panel.
function _node_text(node::pj.YAMLNode)
  io = IOBuffer()
  _walk(io, node, 0)
  return String(take!(io))
end

# A leaf's text. The wrapper's `is_scalar` can be false for keyed numeric leaves,
# so recursion is decided by is_map/is_sequence and everything else is a leaf.
_scalar_str(node) = try String(node) catch; "" end

function _walk(io, node, ind)
  pad = "  "^ind
  if pj.is_map(node)
    for (k, v) in node
      if pj.is_map(v) || pj.is_sequence(v)
        println(io, pad, k, ":")
        _walk(io, v, ind + 1)
      else
        println(io, pad, k, ": ", _scalar_str(v))
      end
    end
  elseif pj.is_sequence(node)
    for child in node
      _walk(io, child, ind)
    end
  else
    println(io, pad, _scalar_str(node))
  end
end

# Pixel offsets for the labels that collided and were stacked (see geometry.jl,
# which measures the stack in font sizes).
#
# The separation has to be in *pixels*: what a stacked label must clear is the
# text below it, which is a fixed size on screen, so a data-space gap would come
# out right at one zoom only. Makie's `offset` is exactly that -- a pixel shift,
# applied without the text's own rotation -- so the direction is resolved here.
# It runs along the label's reading direction, negated for the labels whose
# anchor the readability fold flipped, which points away from the centerline
# either way. The axis is a `DataAspect`, so the data-space angle the geometry
# computed is also the on-screen angle.
function _label_offsets(geom::FloorGeometry, fontsize::Real)
  return map(geom.label_stack, geom.label_rot, geom.label_align) do st, rot, al
    st == 0 && return Vec2f(0, 0)
    d = (al[1] === :left ? 1 : -1) * st * fontsize
    return Vec2f(d * cos(rot), d * sin(rot))
  end
end

"""
    floor_plot(lat; view="zx", shapes=ShapeMap(), title="Floor Plan",
               label_min_px=6, label_sep=1.0, size=(1500, 850)) -> FloorPlot

Build a floor-plan of the expanded lattice `lat`. Returns a `FloorPlot`; what you
do with it is the backend's business — `using GLMakie` and `display(fp)` opens an
interactive window (a script should `wait` on it to keep it open), `using
CairoMakie` and `save("floor.pdf", fp.figure)` writes a file. The figure itself
is built with Makie alone, so no backend is needed to construct one.

Controls, once displayed in an interactive backend (Makie's default `Axis`
bindings, plus selection): **scroll** zooms at the cursor, **right-drag** pans,
**left-drag** rubber-band-zooms to a rectangle, **ctrl-left-click** resets to the
full view, and a plain **left-click** selects an element (listing its parameters
in the side panel). Note that a double-click does nothing — Makie's `Axis` has no
double-click binding.

Keyword arguments:
  * `view`  — projection plane, e.g. `"zx"` (default, horizontal plane from above).
  * `shapes` — a `ShapeMap` controlling element shapes/colors/labels.
  * `label_min_px` — element labels hide when the element's on-screen size falls
    below this many pixels (level-of-detail, so labels don't swamp a zoomed-out
    view of a large machine).
  * `label_sep` — how close two label anchors must be to count as colliding, in
    units of the elements' half-height; colliding labels are stacked outward
    rather than printed on top of each other. `0` disables the stacking.
"""
function floor_plot(lat::pj.Lattices; view::AbstractString="zx",
                    shapes::ShapeMap=ShapeMap(), title::AbstractString="Floor Plan",
                    label_min_px::Real=6, label_sep::Real=1.0, size=(1500, 850))
  tab = element_table(lat)
  return floor_plot(tab; view=view, shapes=shapes, title=title,
                    label_min_px=label_min_px, label_sep=label_sep, size=size)
end

"""
    floor_plot(tab::ElementTable; kwargs...) -> FloorPlot

Draw a floor-plan from an already-extracted `ElementTable`, rather than
extracting one from a lattice. Useful when the table is synthesized directly, or
was built once and is being drawn in several projections.
"""
function floor_plot(tab::ElementTable; view::AbstractString="zx",
                    shapes::ShapeMap=ShapeMap(), title::AbstractString="Floor Plan",
                    label_min_px::Real=6, label_sep::Real=1.0, size=(1500, 850))
  geom = build_geometry(tab, shapes; view=view, label_sep=label_sep)
  return _floor_plot(tab, geom; view=view, title=title,
                     label_min_px=label_min_px, size=size)
end

function _floor_plot(tab::ElementTable, geom::FloorGeometry; view, title,
                     label_min_px, size)
  fig = Figure(; size=size)
  ax = Axis(fig[1, 1]; aspect=DataAspect(), title=title,
            xlabel=_axis_label(view[1]), ylabel=_axis_label(view[2]))
  # Interactions are Makie's Axis defaults -- scroll = zoom, right-drag = pan,
  # left-drag = rubber-band zoom, ctrl-left-click = reset (`LimitReset`), and
  # ctrl-shift-left-click = autolimits. Plain left-click selects (below).

  # Reference orbit (thin black centerline through every element).
  isempty(geom.ref_pts) || lines!(ax, geom.ref_pts; color=geom.ref_col, linewidth=1)

  # Element outlines and interior strokes (one batched call each).
  isempty(geom.outline_pts) ||
    lines!(ax, geom.outline_pts; color=geom.outline_col, linewidth=geom.outline_wid)
  isempty(geom.seg_pts) ||
    linesegments!(ax, geom.seg_pts; color=geom.seg_col, linewidth=1.2)

  # Labels, with level-of-detail: shown only when the view is zoomed in enough
  # that elements are visually resolvable.
  if !isempty(geom.label_pos)
    vis = lift(ax.finallimits, ax.scene.viewport) do lims, vp
      mpp = Float32(lims.widths[1] / vp.widths[1])   # meters per horizontal pixel
      typ_half = 0.5f0                               # typical element half-height (m)
      (2 * typ_half) / mpp >= label_min_px
    end

    fontsize = 11
    text!(ax, geom.label_pos; text=geom.label_str, fontsize=fontsize,
          rotation=geom.label_rot, align=geom.label_align,
          offset=_label_offsets(geom, fontsize),
          visible=vis, color=:gray25)
  end

  # Selection highlight overlay, driven by an Observable element index.
  selected = Observable(0)
  hi_pts = lift(selected) do i
    i == 0 ? Point2f[] : element_outline(tab, i; view=view)
  end
  lines!(ax, hi_pts; color=:red, linewidth=3)

  # Parameter side-panel.
  panel = fig[1, 2] = GridLayout(; tellheight=false)
  colsize!(fig.layout, 2, Fixed(340))
  Label(panel[1, 1], "Element parameters"; fontsize=14, font=:bold,
        halign=:left, tellwidth=false)
  info = Observable("Click an element to see its parameters.")
  Label(panel[2, 1], info; halign=:left, valign=:top, justification=:left,
        fontsize=11, tellheight=false, tellwidth=false, word_wrap=true)

  # Click-to-select. Going through the interaction system rather than raw
  # `mousebutton` events means Makie's mouse state machine has already sorted a
  # plain click from a drag (its threshold is 2 px, so zoom-independent).
  #
  # Ctrl-left-click is Makie's own view reset (its `LimitReset`), so it is not a
  # selection. Deferring to it explicitly matters: `ax.interactions` is an
  # unordered Dict, so which of the two sees the click first is not something
  # this can rely on.
  register_interaction!(ax, :selectelement) do event::MouseEvent, axis
    event.type === MouseEventTypes.leftclick || return Consume(false)
    ispressed(axis.scene, Keyboard.left_control) && return Consume(false)
    lims = axis.finallimits[]
    i = _nearest(geom.ele_center, event.data, 0.5 * lims.widths[1])
    i == 0 && return Consume(false)
    selected[] = i
    info[] = string(tab.name[i], "  [", tab.kind[i], "]\n",
                    "branch: ", tab.branch_names[tab.branch[i]], "\n\n",
                    _node_text(tab.node[i]))
    return Consume(true)
  end

  return FloorPlot(fig, ax, tab, geom, selected)
end
