# Floor-coordinate frames: the placement math every drawing is carried on.
#
# The PALS standard describes the orientation of the branch coordinate system at
# a point on the reference curve by the 3x3 rotation matrix W (coordinates.md,
# Eq. www), whose columns are the branch x, y and z axes in global coordinates.
# An element's `FloorP` gives the placement at its *upstream* end; the drawing has
# to carry that placement across the element to put the far end -- and every
# shape vertex in between -- where the expander put it.
#
# That carry is exactly what the expander itself does (`floor_propagate` with the
# (L, S) pair from `straight_LS`/`bend_LS` in pals-cpp's pals_floor.cpp), so it
# is done here the same way rather than approximated: `place` below is those two
# functions composed, evaluated at an arbitrary fraction along the element.
# pals-cpp carries the orientation as a quaternion; a matrix is used here because
# what the drawing wants at each point is the branch axes, which are its columns.
#
# Both the 2D and the 3D drawing sit on this: 2D projects the result onto a
# plane, 3D sweeps the shape templates through it directly. Doing the geometry in
# 3D and projecting afterwards is what lets a lattice that leaves the horizontal
# plane be drawn correctly in either -- the earlier 2D-only code built its frame
# from the heading `theta` alone, which drew such a lattice as though it stayed
# in the plane.

using GeometryBasics: Point2f, Point3f, Vec3d

"""
    Mat3

A 3x3 rotation matrix in global (floor) coordinates, stored row-major. Used for
the standard's orientation matrix `W`, whose columns are the branch x, y, z axes
([`xaxis`](@ref), [`yaxis`](@ref), [`zaxis`](@ref)).
"""
struct Mat3
  m11::Float64; m12::Float64; m13::Float64
  m21::Float64; m22::Float64; m23::Float64
  m31::Float64; m32::Float64; m33::Float64
end

"The identity rotation."
const I3 = Mat3(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)

@inline function Base.:*(A::Mat3, B::Mat3)
  return Mat3(A.m11 * B.m11 + A.m12 * B.m21 + A.m13 * B.m31,
              A.m11 * B.m12 + A.m12 * B.m22 + A.m13 * B.m32,
              A.m11 * B.m13 + A.m12 * B.m23 + A.m13 * B.m33,
              A.m21 * B.m11 + A.m22 * B.m21 + A.m23 * B.m31,
              A.m21 * B.m12 + A.m22 * B.m22 + A.m23 * B.m32,
              A.m21 * B.m13 + A.m22 * B.m23 + A.m23 * B.m33,
              A.m31 * B.m11 + A.m32 * B.m21 + A.m33 * B.m31,
              A.m31 * B.m12 + A.m32 * B.m22 + A.m33 * B.m32,
              A.m31 * B.m13 + A.m32 * B.m23 + A.m33 * B.m33)
end

@inline function Base.:*(A::Mat3, v::Vec3d)
  return Vec3d(A.m11 * v[1] + A.m12 * v[2] + A.m13 * v[3],
               A.m21 * v[1] + A.m22 * v[2] + A.m23 * v[3],
               A.m31 * v[1] + A.m32 * v[2] + A.m33 * v[3])
end

"Branch x axis (the horizontal transverse direction) in global coordinates."
@inline xaxis(W::Mat3) = Vec3d(W.m11, W.m21, W.m31)
"Branch y axis (the vertical transverse direction) in global coordinates."
@inline yaxis(W::Mat3) = Vec3d(W.m12, W.m22, W.m32)
"Branch z axis (the direction of travel) in global coordinates."
@inline zaxis(W::Mat3) = Vec3d(W.m13, W.m23, W.m33)

# The elementary rotations of the standard (Eq. wtt0t), about the branch x, y and
# z axes.
@inline rot_x(a) = (s = sin(a); c = cos(a); Mat3(1.0, 0.0, 0.0, 0.0, c, -s, 0.0, s, c))
@inline rot_y(a) = (s = sin(a); c = cos(a); Mat3(c, 0.0, s, 0.0, 1.0, 0.0, -s, 0.0, c))
@inline rot_z(a) = (s = sin(a); c = cos(a); Mat3(c, -s, 0.0, s, c, 0.0, 0.0, 0.0, 1.0))

# Rotation by `angle` about a unit `axis` (Rodrigues). `bend_LS` needs it about
# an axis that tilt_ref swings out of the vertical.
function axis_angle(axis::Vec3d, angle)
  n = sqrt(axis[1]^2 + axis[2]^2 + axis[3]^2)
  n <= 0 && return I3
  ux, uy, uz = axis[1] / n, axis[2] / n, axis[3] / n
  s, c = sin(angle), cos(angle)
  t = 1 - c
  return Mat3(t * ux * ux + c,      t * ux * uy - s * uz, t * ux * uz + s * uy,
              t * ux * uy + s * uz, t * uy * uy + c,      t * uy * uz - s * ux,
              t * ux * uz - s * uy, t * uy * uz + s * ux, t * uz * uz + c)
end

"""
    w_matrix(theta, phi, psi) -> Mat3

The standard's orientation matrix `W = R_y(theta) * R_x(phi) * R_z(psi)`
(Eq. www) built from the three `FloorP` angles. Its third column, the direction
of travel, is `(sinθcosφ, -sinφ, cosθcosφ)`; with `phi = psi = 0` it reduces to a
heading in the horizontal plane.
"""
@inline w_matrix(theta, phi, psi) = rot_y(theta) * (rot_x(phi) * rot_z(psi))

"""
    Placement

The position `r` and orientation `W` of the branch coordinate system at one point
on the reference curve -- the standard's floor placement, and what pals-cpp
carries as a `FloorState`.
"""
struct Placement
  r::Vec3d
  W::Mat3
end

"""
    place(r0, W0, len, angle, tilt, f) -> Placement

Placement a fraction `f` of the way along an element whose upstream end is at
`r0` with orientation `W0`, of arc length `len`, reference bend angle `angle`
(0 when straight) and reference tilt `tilt`.

This is pals-cpp's `floor_propagate` applied to the `(L, S)` pair that
`straight_LS`/`bend_LS` build, evaluated at `f * angle` rather than only at the
downstream end. Following the expander's own construction is what makes the
drawing close on its floor coordinates: at `f = 1` this reproduces the next
element's `FloorP` rather than merely coming near it.
"""
@inline function place(r0::Vec3d, W0::Mat3, len, angle, tilt, f)
  if abs(angle) < 1e-12
    # Straight element (Eq. l00l): L = (0, 0, len), S = identity.
    return Placement(r0 + W0 * Vec3d(0.0, 0.0, f * len), W0)
  end
  # Bend (Eqs. lrztt, ustt). The displacement in the untilted bend plane, rotated
  # by tilt_ref about z; the coordinate rotation is about an axis the same tilt
  # swings out of the vertical. A positive angle_ref moves the exit toward
  # negative branch x, which is why the transverse term is rho*(cos - 1).
  #
  # The standard writes the bend's coordinate rotation twice: as an axis and
  # angle (Eq. ustt, u = (sin θ_tr, -cos θ_tr, 0)) and as a conjugation
  # (Eq. srrr, S = R_z(θ_tr) R_y(-α_b) R_z(-θ_tr)). They are the same rotation,
  # and both agree with the displacement in Eq. lrztt -- which is R_z(θ_tr)·L̃,
  # the untilted bend turned bodily about z, so its coordinate rotation is the
  # untilted one conjugated the same way. The axis is just the untilted axis
  # (0, -1, 0) carried around by that same R_z(θ_tr).
  #
  # This used to disagree: Eq. ustt read u = (-sin θ_tr, -cos θ_tr, 0), which is
  # a *different* rotation for a bend rolled out of its branch's x-z plane, and
  # left the frame's z axis off the tangent to the arc it is travelling along by
  # about 2 sin(α_b) sin(θ_tr). pals-cpp's `bend_LS` implemented that sign and
  # this function followed it to keep the drawing flush with the expander's
  # `FloorP`. Both have since been corrected; the `place` tests check the two
  # forms against each other so a regression on either side fails here rather
  # than turning into a mystery about the drawing.
  rho = len / angle
  a = f * angle
  Lt = Vec3d(rho * (cos(a) - 1.0), 0.0, rho * sin(a))
  L = rot_z(tilt) * Lt
  S = axis_angle(Vec3d(sin(tilt), -cos(tilt), 0.0), a)
  return Placement(r0 + W0 * L, W0 * S)
end

"""
    entrance(tab, i) -> Placement

The placement of element `i` at its upstream end, straight from its `FloorP`.
"""
@inline entrance(tab, i::Int) =
  Placement(Vec3d(tab.x[i], tab.y[i], tab.z[i]),
            w_matrix(tab.theta[i], tab.phi[i], tab.psi[i]))

"""
    place(tab, i, f) -> Placement

Placement a fraction `f` along element `i` of an [`ElementTable`](@ref).
"""
@inline function place(tab, i::Int, f)
  e = entrance(tab, i)
  return place(e.r, e.W, tab.length[i], tab.angle[i], tab.tilt_ref[i], f)
end

# ── projection ────────────────────────────────────────────────────────────────
# A `view` string names global axes, one per drawn axis: two characters for a
# plane ("zx" -> horizontal screen axis is global z, vertical is global x), three
# for a 3D drawing ("xzy" -> the drawn vertical is global y). Vectors and points
# go through the same permutation, so a projected direction stays the projection
# of the true direction rather than something reconstructed in the picture plane.

@inline _axis(c::Char, v::Vec3d) = c == 'x' ? v[1] : c == 'y' ? v[2] : v[3]

@inline proj2(a::Char, b::Char, v::Vec3d) = Point2f(_axis(a, v), _axis(b, v))
@inline proj3(a::Char, b::Char, c::Char, v::Vec3d) =
  Point3f(_axis(a, v), _axis(b, v), _axis(c, v))

# Projected branch axes as *unit screen directions*: the direction of travel and
# the horizontal transverse direction, each normalized in the picture plane.
#
# Both are projections of the branch axes rather than one derived from the other
# by rotating 90 degrees in the picture plane: the sign of that rotation depends
# on the handedness of the axis pair, so "zx" and "xz" -- the same plane with the
# screen axes swapped -- would get opposite transverse directions, and a bend
# would curve the wrong way in one of them.
#
# Normalizing is right for *directions* (a label has to run somewhere even when
# the element it belongs to is pointed at the viewer) but wrong for the geometry,
# which uses the unnormalized projection so that an element tilted out of the
# picture plane is drawn foreshortened, as it should be.
@inline function proj_axes(a::Char, b::Char, W::Mat3)
  g = proj2(a, b, zaxis(W)); n = proj2(a, b, xaxis(W))
  gn = hypot(g[1], g[2]); g = gn > 1.0e-12 ? g ./ gn : Point2f(1, 0)
  nn = hypot(n[1], n[2]); n = nn > 1.0e-12 ? n ./ nn : Point2f(-g[2], g[1])
  return g[1], g[2], n[1], n[2]
end
