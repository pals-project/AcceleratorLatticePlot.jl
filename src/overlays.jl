# Overlays: things drawn on a floor plan that are not lattice elements.
#
# Two kinds, both taking their input in *global* coordinates and projecting it
# the way the plot it is added to was built, so the same call works on a 2D and a
# 3D drawing:
#
#   * `add_curve!` -- an arbitrary polyline. What a reference-orbit or
#     measured-orbit overlay is: a curve through the machine that did not come
#     out of the shape table.
#   * `add_wall!`  -- a building wall, given as its outline on the floor and a
#     height. Flat in the floor plan, extruded in 3D.
#
# Both are ordinary Makie calls on the existing axis, so they pan, zoom and
# rotate with everything else, and both return the plot object for the caller to
# tweak or `delete!`.

using GeometryBasics: Point2f, Point3f, Vec3d, GLTriangleFace
using Colors
using Makie

# Accept any reasonable spelling of a point: a Vec/Point/tuple of 3, or of 2 in
# which case the missing global coordinate is the vertical (`y`), taken as 0.
_as_vec3(p::Vec3d) = p
_as_vec3(p) = length(p) >= 3 ? Vec3d(p[1], p[2], p[3]) : Vec3d(p[1], 0.0, p[2])

"""
    add_curve!(fp, points; color=:red, linewidth=2, kwargs...)

Draw a polyline through `points`, given in global coordinates, on the floor plan
or 3D drawing `fp`. Points may be 3-vectors `(x, y, z)`, or 2-vectors taken as
`(x, z)` on the `y = 0` plane. Extra keyword arguments go to Makie's `lines!`.

This is the overlay a reference or measured orbit goes on: AcceleratorLatticePlot
has no opinion about where the curve comes from, only about placing it in the
same coordinates as the machine.

```julia
fp = floor_plot(lat)
add_curve!(fp, orbit_points; color=:orange, linewidth=3)
```
"""
function add_curve!(fp::FloorPlot, points; color=:red, linewidth=2, kwargs...)
  a, b = fp.view[1], fp.view[2]
  pts = [proj2(a, b, _as_vec3(p)) for p in points]
  return lines!(fp.axis, pts; color=color, linewidth=linewidth, kwargs...)
end

function add_curve!(fp::FloorPlot3, points; color=:red, linewidth=2, kwargs...)
  a, b, c = fp.view[1], fp.view[2], fp.view[3]
  pts = [proj3(a, b, c, _as_vec3(p)) for p in points]
  return lines!(fp.axis, pts; color=color, linewidth=linewidth, kwargs...)
end

"""
    add_wall!(fp, outline; base=0.0, height=3.0, color=:gray70, alpha=0.25,
              kwargs...)

Draw a building wall whose footprint on the floor is the polyline `outline`,
given in global coordinates and read in the horizontal plane. `base` and `height`
are the wall's bottom and top on the global vertical axis.

On a floor plan this is the footprint, a plain polyline: the floor plan is a plan,
and a wall's height does not show in it. On a 3D drawing it is the wall itself,
extruded from `base` to `base + height` and drawn translucent so the machine
inside stays visible.

```julia
add_wall!(fp3, [(0, 0, -5), (0, 0, 120), (30, 0, 120)]; height=4)
```
"""
function add_wall!(fp::FloorPlot, outline; base=0.0, height=3.0, color=:gray70,
                   alpha=0.25, linewidth=2, kwargs...)
  a, b = fp.view[1], fp.view[2]
  pts = [proj2(a, b, _as_vec3(p)) for p in outline]
  return lines!(fp.axis, pts; color=(color, max(alpha, 0.6)),
                linewidth=linewidth, kwargs...)
end

function add_wall!(fp::FloorPlot3, outline; base=0.0, height=3.0, color=:gray70,
                   alpha=0.25, kwargs...)
  a, b, c = fp.view[1], fp.view[2], fp.view[3]
  up = c            # the global axis the drawn vertical shows

  # Two triangles per segment, from the footprint at `base` to the footprint at
  # `base + height`.
  pts = Point3f[]; faces = GLTriangleFace[]
  raise(p, h) = proj3(a, b, c, _with_axis(_as_vec3(p), up, h))
  for k in 1:(length(outline) - 1)
    n = length(pts)
    push!(pts, raise(outline[k], base), raise(outline[k + 1], base),
               raise(outline[k + 1], base + height), raise(outline[k], base + height))
    push!(faces, GLTriangleFace(n + 1, n + 2, n + 3),
                 GLTriangleFace(n + 1, n + 3, n + 4))
  end
  isempty(faces) && return nothing
  # No normals: a wall is a flat slab seen from both sides, and letting Makie
  # light it either way looks better than picking a side to face.
  return mesh!(fp.axis, pts, faces; color=(color, alpha), transparency=true,
               shading=Makie.NoShading, kwargs...)
end

# A copy of `v` with the global axis named by `c` set to `h`.
@inline _with_axis(v::Vec3d, c::Char, h) =
  c == 'x' ? Vec3d(h, v[2], v[3]) :
  c == 'y' ? Vec3d(v[1], h, v[3]) : Vec3d(v[1], v[2], h)
