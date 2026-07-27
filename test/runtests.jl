using Test
using PALSPlot
using PALSJulia          # for parse_and_expand_pals in the extraction tests
import PALSPlot as pp
import PALSJulia as pj

# CairoMakie is the backend under test. PALSPlot itself depends only on Makie, so
# a backend is needed here purely to prove the figures render -- and Cairo is
# pure software, which means these tests run anywhere, including the CI runners
# that cannot give GLMakie an OpenGL context.
#
# Loading a Makie backend and then calling into the pals-cpp parser is also the
# pairing that aborts the process when the library was built against the wrong
# C++ runtime (see the README), so running the two together is deliberate.
using CairoMakie

# Triangle vertex indices as plain Ints. Faces are GeometryBasics `OffsetInteger`
# triples, which index correctly but do not convert to Int on their own.
faceverts(f) = Int.(pp.GeometryBasics.value.(Tuple(f)))

# A throwaway node to fill ElementTable.node in synthetic tables.
const NODE = pj.parse_string("kind: Drift\n")

# Build a one-branch ElementTable by hand. `phi`, `psi` and `tilt` default to
# zero, which is the horizontal-plane case; pass them to build a lattice that
# leaves the plane.
function synth_table(; names, kinds, lengths, x, z, theta, angle,
                     y=nothing, phi=nothing, psi=nothing, tilt=nothing)
    n = length(names)
    zed(v) = v === nothing ? zeros(n) : collect(v)
    return pp.ElementTable(collect(names), collect(kinds), collect(lengths),
                           collect(x), zed(y), collect(z), collect(theta),
                           zed(phi), zed(psi), collect(angle), zed(tilt),
                           collect(z), ones(Int, n), fill(NODE, n), ["b"])
end

@testset "shape rules" begin
    r = ele_shape("Quadrupole::q*", :xbox, :black; size=0.5, label=:name)
    @test r.class == "Quadrupole"
    @test r.glob == "q*"
    @test r.spec.shape == :xbox
    @test r.spec.size == 0.5

    r2 = ele_shape("Drift", :none, :gray60)   # no "::" -> glob "*"
    @test r2.glob == "*"
end

@testset "default shape map" begin
    m = ShapeMap()
    @test mapshape(m, "Quadrupole", "q1").shape == :xbox
    @test mapshape(m, "Bend", "b1").shape == :box
    @test mapshape(m, "Sextupole", "s1").shape == :box
    @test mapshape(m, "Drift", "d1").shape == :none
    # unknown kind falls through to the "*" default
    @test mapshape(m, "Wibbler", "w1").shape == :box
    @test mapshape(m, "Wibbler", "w1").label == :none
end

# Every kind a beam line can hold should be in the table rather than reaching the
# "*" fallback, so a kind renamed in the standard shows up here as a failure.
@testset "default table covers the PALS element kinds" begin
    m = ShapeMap()
    kinds = ["ACKicker", "BeamBeam", "BeginningEle", "Bend", "Converter",
             "CrabCavity", "Drift", "EGun", "Fiducial", "FloorShift", "Foil",
             "Fork", "Instrument", "Kicker", "Marker", "Mask", "Match",
             "Multipole", "Octupole", "Patch", "Placeholder", "Quadrupole",
             "ReferenceChange", "RFCavity", "Sextupole", "Solenoid", "Taylor",
             "UnionEle", "Wiggler"]
    fallback = mapshape(m, "Wibbler", "w1")
    for k in kinds
        @test mapshape(m, k, "e1") !== fallback
    end
end

# The orientation matrix of the standard (Eq. www). Everything drawn is carried
# on this, so it is checked against the standard's own statement of it rather
# than against the code's.
@testset "W matrix" begin
    for (th, ph, ps) in ((0.0, 0.0, 0.0), (0.3, 0.0, 0.0), (0.3, -0.2, 0.7),
                         (-1.1, 0.4, -0.9), (2.8, -0.63, -2.17))
        W = pp.w_matrix(th, ph, ps)

        # The direction of travel is the third column, which the standard gives
        # in closed form.
        z = pp.zaxis(W)
        @test isapprox(z[1], sin(th) * cos(ph); atol=1e-12)
        @test isapprox(z[2], -sin(ph); atol=1e-12)
        @test isapprox(z[3], cos(th) * cos(ph); atol=1e-12)

        # ...and W is a rotation: its columns are orthonormal and right-handed.
        cols = (pp.xaxis(W), pp.yaxis(W), pp.zaxis(W))
        for (i, u) in enumerate(cols), (j, v) in enumerate(cols)
            @test isapprox(pp.dot3(u, v), i == j ? 1.0 : 0.0; atol=1e-12)
        end
        @test isapprox(pp.dot3(pp.cross3(cols[1], cols[2]), cols[3]), 1.0; atol=1e-6)
    end

    # With no pitch or roll it is a heading in the horizontal plane, which is all
    # the drawing used to know about.
    W = pp.w_matrix(0.7, 0.0, 0.0)
    @test isapprox(pp.zaxis(W), pp.Vec3d(sin(0.7), 0, cos(0.7)); atol=1e-12)
    @test isapprox(pp.xaxis(W), pp.Vec3d(cos(0.7), 0, -sin(0.7)); atol=1e-12)
end

# `place` is pals-cpp's floor_propagate applied to straight_LS/bend_LS, evaluated
# partway along the element. Each of those pieces is checked separately, and then
# the property that matters most: propagating in steps must agree with
# propagating in one go, since that is what carrying a shape along an arc does.
@testset "reference-curve propagation" begin
    r0 = pp.Vec3d(1.0, 2.0, 3.0)
    W0 = pp.w_matrix(0.4, -0.2, 0.9)

    # f = 0 is the element's own upstream placement, untouched.
    p = pp.place(r0, W0, 5.0, 0.3, 0.1, 0.0)
    @test isapprox(p.r, r0; atol=1e-12)
    @test isapprox(pp.zaxis(p.W), pp.zaxis(W0); atol=1e-12)

    # Straight: L = (0, 0, len) in the branch frame, orientation unchanged.
    p = pp.place(r0, W0, 5.0, 0.0, 0.0, 1.0)
    @test isapprox(p.r, r0 + 5.0 * pp.zaxis(W0); atol=1e-12)
    @test isapprox(pp.xaxis(p.W), pp.xaxis(W0); atol=1e-12)

    # Bend, untilted (Eq. lrztt): rho*sin(angle) along z, rho*(cos(angle)-1)
    # along x -- a positive angle moving the exit toward negative branch x.
    len, ang = 4.0, 0.7
    rho = len / ang
    p = pp.place(r0, W0, len, ang, 0.0, 1.0)
    expect = r0 + rho * sin(ang) * pp.zaxis(W0) + rho * (cos(ang) - 1) * pp.xaxis(W0)
    @test isapprox(p.r, expect; atol=1e-12)

    # A tilt of pi/2 rolls the same arc into the vertical plane: what was a
    # displacement along branch x becomes one along branch y.
    p = pp.place(r0, W0, len, ang, π / 2, 1.0)
    expect = r0 + rho * sin(ang) * pp.zaxis(W0) + rho * (cos(ang) - 1) * pp.yaxis(W0)
    @test isapprox(p.r, expect; atol=1e-9)

    # Propagation composes: half an element, then half again from there, is the
    # whole element. A drawing that got this wrong would still start and end in
    # the right places while bulging in between.
    for ang in (0.0, 0.7, -0.4)
        whole = pp.place(r0, W0, len, ang, 0.0, 1.0)
        half = pp.place(r0, W0, len, ang, 0.0, 0.5)
        rest = pp.place(half.r, half.W, len / 2, ang / 2, 0.0, 1.0)
        @test isapprox(whole.r, rest.r; atol=1e-9)
        @test isapprox(pp.zaxis(whole.W), pp.zaxis(rest.W); atol=1e-9)
        @test isapprox(pp.xaxis(whole.W), pp.xaxis(rest.W); atol=1e-9)
    end

    # A bend of vanishing angle is a straight element, not a division by zero.
    p = pp.place(r0, W0, 2.0, 0.0, 0.4, 1.0)
    @test isapprox(p.r, r0 + 2.0 * pp.zaxis(W0); atol=1e-12)
end

# A tilted bend follows pals-cpp's `bend_LS`, which implements the standard's
# Eq. ustt. That is a deliberate choice and not an obvious one, because the
# standard's own alternative form for the same rotation, Eq. srrr, is a *different
# rotation*: it is the untilted bend conjugated by R_z(tilt_ref), which is what
# the displacement in Eq. lrztt already is, whereas Eq. ustt flips the sign of the
# rotation axis's first component. Under Eq. ustt the frame's z axis leaves the
# tangent to its own arc.
#
# PALSPlot follows the expander so that the drawing meets the `FloorP` the
# expander wrote at the bend's exit face (see the closure tests below). These
# tests exist to say which convention that is out loud: if pals-cpp switches to
# Eq. srrr, they fail here and point at the reason.
@testset "tilted bends follow the expander's convention" begin
    ustt(a, t) = pp.axis_angle(pp.Vec3d(-sin(t), -cos(t), 0.0), a)
    srrr(a, t) = pp.rot_z(t) * (pp.rot_y(-a) * pp.rot_z(-t))
    r0 = pp.Vec3d(0.0, 0.0, 0.0)

    for (a, t) in ((0.4, 0.3), (0.4, π / 2), (-0.9, 1.3))
        p = pp.place(r0, pp.I3, 2.0, a, t, 1.0)
        @test isapprox(pp.zaxis(p.W), pp.zaxis(ustt(a, t)); atol=1e-9)

        # The two forms really are different rotations, so following one rather
        # than the other is a choice with consequences.
        @test !isapprox(pp.zaxis(ustt(a, t)), pp.zaxis(srrr(a, t)); atol=1e-3)

        # ...and it is Eq. srrr whose frame stays tangent to the arc the
        # displacement of Eq. lrztt traces out.
        tangent = pp.rot_z(t) * pp.Vec3d(-sin(a), 0.0, cos(a))
        @test isapprox(pp.zaxis(srrr(a, t)), tangent; atol=1e-9)
        @test !isapprox(pp.zaxis(ustt(a, t)), tangent; atol=1e-3)
    end

    # With no tilt the question does not arise: the two forms coincide.
    for a in (0.4, -0.9)
        @test isapprox(pp.zaxis(ustt(a, 0.0)), pp.zaxis(srrr(a, 0.0)); atol=1e-12)
    end
end

@testset "rule precedence and globbing" begin
    m = ShapeMap([ele_shape("Quadrupole::qs*", :circle, :red)])
    @test mapshape(m, "Quadrupole", "qs3").shape == :circle   # user rule
    @test mapshape(m, "Quadrupole", "q9").shape == :xbox      # default
    @test mapshape(m, "Quadrupole", "QS3").shape == :circle   # case-insensitive

    m2 = ShapeMap([ele_shape("Quadrupole::qs*", :circle, :red)]; defaults=false)
    # with defaults off, non-matching elements get the neutral fallback
    @test mapshape(m2, "Quadrupole", "q9").shape == :box
end

@testset "geometry: straight box placement" begin
    tab = synth_table(names=["s1"], kinds=["Sextupole"], lengths=[0.4],
                      x=[0.0], z=[0.0], theta=[0.0], angle=[0.0])
    g = pp.build_geometry(tab, ShapeMap(); view="zx")
    @test length(g.ele_center) == 1
    # midpoint at f=0.5 -> (z, x) = (0.2, 0)
    @test isapprox(g.ele_center[1][1], 0.2; atol=1e-5)
    @test isapprox(g.ele_center[1][2], 0.0; atol=1e-5)
    # first outline vertex is box corner (f=0, w=-1) -> (0, -0.4)
    @test isapprox(g.outline_pts[1][1], 0.0; atol=1e-5)
    @test isapprox(g.outline_pts[1][2], -0.4; atol=1e-5)
    # a plain box has no interior strokes
    @test isempty(g.seg_pts)
end

@testset "geometry: xbox has interior X" begin
    tab = synth_table(names=["q1"], kinds=["Quadrupole"], lengths=[0.3],
                      x=[0.0], z=[0.0], theta=[0.0], angle=[0.0])
    g = pp.build_geometry(tab, ShapeMap(); view="zx")
    @test length(g.seg_pts) == 4        # two diagonal segments
end

@testset "geometry: bend follows an arc" begin
    straight = pp.build_geometry(
        synth_table(names=["b"], kinds=["Bend"], lengths=[1.0],
                    x=[0.0], z=[0.0], theta=[0.0], angle=[0.0]),
        ShapeMap(); view="zx")
    bent = pp.build_geometry(
        synth_table(names=["b"], kinds=["Bend"], lengths=[1.0],
                    x=[0.0], z=[0.0], theta=[0.0], angle=[0.6]),
        ShapeMap(); view="zx")
    # subdividing the arc yields more outline vertices than the straight chord
    @test length(bent.outline_pts) > length(straight.outline_pts)
end

@testset "geometry: view selects projection plane" begin
    tab = synth_table(names=["d"], kinds=["Sextupole"], lengths=[1.0],
                      x=[0.0], z=[0.0], theta=[0.0], angle=[0.0])
    gzx = pp.build_geometry(tab, ShapeMap(); view="zx")  # horiz=z
    gxz = pp.build_geometry(tab, ShapeMap(); view="xz")  # horiz=x
    @test isapprox(gzx.ele_center[1][1], 0.5; atol=1e-5)  # advances along z
    @test isapprox(gxz.ele_center[1][2], 0.5; atol=1e-5)  # ...now on vertical axis
end

# "zx" and "xz" are the same plane with the screen axes swapped, so they must be
# the same drawing transposed -- including which way a bend curves. Deriving the
# transverse direction by rotating the heading 90° in the picture plane got this
# wrong, because the sign of that rotation depends on the handedness of the axis
# pair; the frame is built from the projected branch axes instead.
@testset "geometry: swapping the view axes transposes the drawing" begin
    tab = synth_table(names=["b"], kinds=["Bend"], lengths=[2.0],
                      x=[0.0], z=[0.0], theta=[0.3], angle=[0.7])
    gzx = pp.build_geometry(tab, ShapeMap(); view="zx")
    gxz = pp.build_geometry(tab, ShapeMap(); view="xz")
    @test length(gzx.outline_pts) == length(gxz.outline_pts)
    for (p, q) in zip(gzx.outline_pts, gxz.outline_pts)
        if isnan(p[1])
            @test isnan(q[1])
        else
            @test isapprox(p[1], q[2]; atol=1e-5)
            @test isapprox(p[2], q[1]; atol=1e-5)
        end
    end
end

# Sign convention, fixed by pals-cpp's bend_LS (Eq. lrztt): the displacement in
# the bend plane is (rho*(cos(angle)-1), 0, rho*sin(angle)) along the branch
# axes, so a positive angle moves the exit toward *negative* branch x. Drawn in
# "zx" from theta = 0, branch x is the vertical screen axis.
@testset "geometry: a positive bend angle curves toward -x" begin
    tab = synth_table(names=["b"], kinds=["Bend"], lengths=[1.0],
                      x=[0.0], z=[0.0], theta=[0.0], angle=[0.5])
    g = pp.build_geometry(tab, ShapeMap(); view="zx")
    @test g.ele_center[1][2] < 0            # midpoint has swung to -x
    @test g.ele_center[1][1] > 0            # ...while still advancing along z
    # and the mirror image for a negative angle
    tabm = synth_table(names=["b"], kinds=["Bend"], lengths=[1.0],
                       x=[0.0], z=[0.0], theta=[0.0], angle=[-0.5])
    gm = pp.build_geometry(tabm, ShapeMap(); view="zx")
    @test isapprox(gm.ele_center[1][2], -g.ele_center[1][2]; atol=1e-5)
    @test isapprox(gm.ele_center[1][1], g.ele_center[1][1]; atol=1e-5)
end

# Every shape must produce drawable geometry: no error, no NaN except the loop
# separators, and something actually emitted for the shapes that draw.
@testset "geometry: every shape renders" begin
    for shape in pp.SHAPE_NAMES
        tab = synth_table(names=["e"], kinds=["Widget"], lengths=[1.0],
                          x=[0.0], z=[0.0], theta=[0.4], angle=[0.3])
        m = ShapeMap([ele_shape("Widget", shape, :black; size=0.5)]; defaults=false)
        g = pp.build_geometry(tab, m; view="zx")
        pts = vcat(g.outline_pts, g.seg_pts)
        finite = filter(p -> !isnan(p[1]), pts)
        @test all(p -> isfinite(p[1]) && isfinite(p[2]), finite)
        shape === :none || @test !isempty(finite)
        @test length(g.outline_col) == length(g.outline_pts)
        @test length(g.outline_wid) == length(g.outline_pts)
        @test length(g.seg_col) == length(g.seg_pts)
    end
end

@testset "geometry: label content follows the rule" begin
    tab = synth_table(names=["q1"], kinds=["Quadrupole"], lengths=[1.0],
                      x=[0.0], z=[7.0], theta=[0.0], angle=[0.0])
    byname = pp.build_geometry(tab, ShapeMap([ele_shape("Quadrupole", :box, :black;
                                                        label=:name)]); view="zx")
    @test byname.label_str == ["q1"]
    bys = pp.build_geometry(tab, ShapeMap([ele_shape("Quadrupole", :box, :black;
                                                     label=:s)]); view="zx")
    @test bys.label_str == ["7.0"]          # synth_table uses z as s
    none = pp.build_geometry(tab, ShapeMap([ele_shape("Quadrupole", :box, :black;
                                                      label=:none)]); view="zx")
    @test isempty(none.label_str)
end

# Labels are laid across the line rather than along it, so that neighbours --
# which are strung out along the line -- do not run into each other.
@testset "geometry: labels sit perpendicular to the centerline" begin
    smap = ShapeMap([ele_shape("Quadrupole", :box, :black; label=:name)])
    # Headings all round the compass, one element each, so the label rotation
    # has to track the centerline rather than a fixed direction.
    for (k, th) in enumerate(range(0, 2π, length=17)[1:16])
        tab = synth_table(names=["q$k"], kinds=["Quadrupole"], lengths=[1.0],
                          x=[0.0], z=[0.0], theta=[th], angle=[0.0])
        g = pp.build_geometry(tab, smap; view="zx")
        rot = only(g.label_rot)

        # Always readable: never upside down, never past the vertical.
        @test -π / 2 - 1e-6 <= rot <= π / 2 + 1e-6

        # And perpendicular to the element: the centerline runs along ĝ, the
        # text along the rotation, so the two are a right angle apart (mod π,
        # since the readability fold can reverse the text direction).
        gx, gy, = pp.proj_axes('z', 'x', pp.w_matrix(th, 0, 0))
        along = atan(gy, gx)
        @test isapprox(mod(rot - along, π), π / 2; atol = 1e-5)
    end

    # A label starts at its anchor and runs outward, away from the centerline,
    # whichever way the readability fold left it pointing.
    for th in (0.0, π)          # normal points +x on one, -x on the other
        tab = synth_table(names=["q1"], kinds=["Quadrupole"], lengths=[1.0],
                          x=[0.0], z=[0.0], theta=[th], angle=[0.0])
        g = pp.build_geometry(tab, smap; view="zx")
        _, _, nx, ny = pp.proj_axes('z', 'x', pp.w_matrix(th, 0, 0))
        rot = only(g.label_rot)
        reading = (cos(rot), sin(rot))          # direction the glyphs run
        outward = (nx, ny)                      # direction the anchor is offset
        toward_outside = reading[1] * outward[1] + reading[2] * outward[2]
        # :left grows along the reading direction, :right against it; either way
        # the growth must have a positive component along the outward normal.
        grow = only(g.label_align)[1] === :left ? toward_outside : -toward_outside
        @test grow > 0
    end
end

# Overlapping elements -- a pickup with correctors wound around it, say -- put
# their labels on the same anchor, and being perpendicular does not separate
# them: they share a ray. Those get stacked out along it.
@testset "geometry: colliding labels stack outward in branch order" begin
    smap = ShapeMap([ele_shape("Instrument", :box, :black; label=:name),
                     ele_shape("Kicker", :box, :black; label=:name)])
    # Three zero-length elements at one spot, the shape in the bta example.
    coincident = synth_table(names=["pue_a12", "dhca12", "dvca12"],
                             kinds=["Instrument", "Kicker", "Kicker"],
                             lengths=[0.0, 0.0, 0.0], x=[0.0, 0.0, 0.0],
                             z=[10.0, 10.0, 10.0], theta=zeros(3), angle=zeros(3))
    g = pp.build_geometry(coincident, smap; view="zx")
    @test g.label_str == ["pue_a12", "dhca12", "dvca12"]

    # All three keep the one anchor: the stack is a shift off it, not a move.
    @test allequal(g.label_pos)

    # First in the branch stays on the centerline, the rest step out past it, in
    # branch order.
    @test g.label_stack[1] == 0
    @test issorted(g.label_stack)
    @test allunique(g.label_stack)

    # Each step clears the label it has to get past: at least its length in
    # characters, times the mean glyph advance, plus a gap.
    for k in 2:3
        step = g.label_stack[k] - g.label_stack[k - 1]
        @test step >= length(g.label_str[k - 1]) * pp._CHAR_EM
    end

    # Elements far enough apart are left where they are.
    apart = synth_table(names=["pue_a12", "dhca12", "dvca12"],
                        kinds=["Instrument", "Kicker", "Kicker"],
                        lengths=[0.0, 0.0, 0.0], x=[0.0, 0.0, 0.0],
                        z=[10.0, 20.0, 30.0], theta=zeros(3), angle=zeros(3))
    @test all(iszero, pp.build_geometry(apart, smap; view="zx").label_stack)

    # ...and `label_sep=0` turns the stacking off entirely.
    off = pp.build_geometry(coincident, smap; view="zx", label_sep=0)
    @test all(iszero, off.label_stack)

    # The shift reaches the plot as a pixel offset pointing away from the line.
    offs = pp._label_offsets(g, 11)
    @test offs[1] == Makie.Vec2f(0, 0)
    _, _, nx, ny = pp.proj_axes('z', 'x', pp.w_matrix(0.0, 0, 0))
    for k in 2:3
        @test offs[k][1] * nx + offs[k][2] * ny > 0      # outward, not inward
        @test hypot(offs[k]...) ≈ 11 * g.label_stack[k] rtol = 1e-5
    end
end

# ── 3D geometry ───────────────────────────────────────────────────────────────

# The mesh is handed to a renderer whole, so what matters is that it is *valid*:
# indices in range, no NaN, unit normals, and a vertex-to-element map that lines
# up with the vertices, since picking reads elements out of it.
@testset "3D geometry: the mesh is well formed" begin
    for shape in pp.SHAPE_NAMES
        tab = synth_table(names=["e"], kinds=["Widget"], lengths=[1.0],
                          x=[0.0], z=[0.0], theta=[0.4], angle=[0.3])
        m = ShapeMap([ele_shape("Widget", shape, :black; size=0.5, size2=0.25)];
                     defaults=false)
        g = build_geometry3(tab, m)

        @test length(g.mesh_nrm) == length(g.mesh_pts)
        @test length(g.mesh_col) == length(g.mesh_pts)
        @test length(g.vertex_ele) == length(g.mesh_pts)
        @test all(p -> all(isfinite, p), g.mesh_pts)
        @test all(n -> isapprox(hypot(n...), 1; atol=1e-4), g.mesh_nrm)
        @test all(f -> all(k -> 1 <= k <= length(g.mesh_pts), faceverts(f)), g.mesh_faces)
        @test all(==(1), g.vertex_ele)
        # Every shape but the two that draw nothing solid makes a solid.
        shape in (:x, :none) || @test !isempty(g.mesh_faces)
    end
end

# Faces must be wound to agree with the normal they carry, or half a magnet
# lights as though it were inside out.
@testset "3D geometry: faces are wound to match their normals" begin
    for shape in (:box, :xbox, :diamond, :bow_tie, :rbow_tie, :u_triangle,
                  :r_triangle, :circle)
        for ang in (0.0, 0.9)      # straight, and curved so the caps get sliced
            tab = synth_table(names=["e"], kinds=["Widget"], lengths=[2.0],
                              x=[0.0], z=[0.0], theta=[0.4], angle=[ang])
            m = ShapeMap([ele_shape("Widget", shape, :black; size=0.5, size2=0.25)];
                         defaults=false)
            g = build_geometry3(tab, m)
            for f in g.mesh_faces
                i1, i2, i3 = faceverts(f)
                gn = pp.cross3(g.mesh_pts[i2] - g.mesh_pts[i1],
                               g.mesh_pts[i3] - g.mesh_pts[i1])
                @test hypot(gn...) > 1e-9              # no slivers were emitted
                @test pp.dot3(gn, g.mesh_nrm[i1]) > 0  # ...and it faces outward
            end
        end
    end
end

# A straight box is the case with an answer that can be written down: a
# rectangular prism of the element's length by 2*size by 2*size2, and no more
# than eight distinct corners.
@testset "3D geometry: a straight box is a box" begin
    tab = synth_table(names=["e"], kinds=["Widget"], lengths=[2.0],
                      x=[0.0], z=[0.0], theta=[0.0], angle=[0.0])
    m = ShapeMap([ele_shape("Widget", :box, :black; size=0.5, size2=0.25)];
                 defaults=false)
    g = build_geometry3(tab, m; view="zxy")     # drawn x = global z, along it
    @test length(unique(g.mesh_pts)) == 8
    @test length(g.mesh_faces) == 12            # 6 quads
    @test extrema(p -> p[1], g.mesh_pts) == (0.0f0, 2.0f0)
    @test extrema(p -> p[2], g.mesh_pts) == (-0.5f0, 0.5f0)
    @test extrema(p -> p[3], g.mesh_pts) == (-0.25f0, 0.25f0)

    # size2 defaults to size, which makes it square in cross-section.
    g2 = build_geometry3(tab, ShapeMap([ele_shape("Widget", :box, :black; size=0.4)];
                                       defaults=false); view="zxy")
    @test extrema(p -> p[3], g2.mesh_pts) == (-0.4f0, 0.4f0)
end

# The point of extruding the 2D profile rather than inventing a cross-section:
# the 3D drawing seen from above *is* the floor plan. Checked as the silhouette,
# since that is the part the two have to agree on.
@testset "3D geometry: the silhouette from above is the floor plan" begin
    lat_tab = synth_table(names=["b", "q", "s"], kinds=["Bend", "Quadrupole", "Sextupole"],
                          lengths=[2.0, 0.5, 0.3], x=[0.0, 0.0, 0.0],
                          z=[0.0, 2.0, 2.5], theta=[0.0, 0.3, 0.3],
                          angle=[0.3, 0.0, 0.0])
    g2 = pp.build_geometry(lat_tab, ShapeMap(); view="zx")
    g3 = build_geometry3(lat_tab, ShapeMap(); view="zxy")
    flat = filter(p -> !isnan(p[1]), g2.outline_pts)
    for k in 1:2
        lo2, hi2 = extrema(p -> p[k], flat)
        lo3, hi3 = extrema(p -> p[k], g3.mesh_pts)
        @test isapprox(lo2, lo3; atol=1e-4)
        @test isapprox(hi2, hi3; atol=1e-4)
    end
end

@testset "3D geometry: view permutes the drawn axes" begin
    tab = synth_table(names=["e"], kinds=["Widget"], lengths=[2.0],
                      x=[0.0], z=[0.0], theta=[0.0], angle=[0.0])
    m = ShapeMap([ele_shape("Widget", :box, :black; size=0.5)]; defaults=false)
    zxy = build_geometry3(tab, m; view="zxy")
    xzy = build_geometry3(tab, m; view="xzy")
    @test extrema(p -> p[1], zxy.mesh_pts) == extrema(p -> p[2], xzy.mesh_pts)
    @test extrema(p -> p[2], zxy.mesh_pts) == extrema(p -> p[1], xzy.mesh_pts)
    @test_throws ArgumentError build_geometry3(tab, m; view="zx")
end

# The geometry fixes only where a label is anchored -- one half-height off the
# centerline on the element's +x, with the matching +y probe that lets the
# renderer measure the element on screen. Nothing about spacing happens here:
# in 3D that is a screen-space question and belongs to `_layout_labels3`.
@testset "3D geometry: label anchors straddle the centerline" begin
    smap = ShapeMap([ele_shape("Widget", :box, :black; size=0.5, size2=0.2,
                               label=:name)]; defaults=false)
    tab = synth_table(names=["w"], kinds=["Widget"], lengths=[2.0], x=[0.0],
                      z=[0.0], theta=[0.0], angle=[0.0])
    g = build_geometry3(tab, smap; view="zxy")
    # Drawn axes are (global z, global x, global y); the element runs along +z
    # unrotated, so its local x is drawn y and its local y is drawn z.
    @test only(g.label_pos) ≈ Point3f(1.0, 0.52, 0.0)
    @test only(g.label_probe) ≈ Point3f(1.0, 0.0, 0.22)
    @test only(g.label_ele) == 1
    @test g.label_sep == 1.0

    # Coincident elements come out with coincident anchors: the geometry does
    # not try to separate them, because it cannot know what "apart" means on a
    # screen it has not got a camera for.
    co = synth_table(names=["pue_a12", "dhca12", "dvca12"],
                     kinds=["Widget", "Widget", "Widget"],
                     lengths=[0.0, 0.0, 0.0], x=[0.0, 0.0, 0.0],
                     z=[10.0, 10.0, 10.0], theta=zeros(3), angle=zeros(3))
    gc = build_geometry3(co, smap; view="zxy")
    @test gc.label_str == ["pue_a12", "dhca12", "dvca12"]
    @test allequal(gc.label_pos)
end

# The layout pass, which is the whole point of the 3D label handling: a label's
# size, the room it needs and the direction anything moves it in are all pixel
# quantities, so they are resolved against a live camera rather than baked into
# the geometry in meters.
@testset "3D labels are laid out in screen space" begin
    smap = ShapeMap([ele_shape("Widget", :box, :black; size=0.3, label=:name)];
                    defaults=false)
    co = synth_table(names=["pue_a12", "dhca12", "dvca12"],
                     kinds=fill("Widget", 3), lengths=zeros(3), x=zeros(3),
                     z=fill(10.0, 3), theta=zeros(3), angle=zeros(3))

    # Where each label actually lands on screen: its anchor, projected, plus the
    # pixel offset the layout gave it.
    function screen_at(fp, L)
        [Makie.project(fp.axis.scene, fp.geometry.label_pos[k]) .+ L.offset[k]
         for k in eachindex(L.pos) if all(isfinite, L.pos[k])]
    end

    fp = floor_plot3(co; shapes=smap, size=(900, 700))
    Makie.update_state_before_display!(fp.figure)
    L = pp._layout_labels3(fp.geometry, fp.axis.scene, 4, 11)

    # All three are drawn, none on top of another, and the two that had to be
    # moved to manage it are joined back to their element by a leader.
    @test count(p -> all(isfinite, p), L.pos) == 3
    pts = screen_at(fp, L)
    @test allunique(pts)
    @test minimum(abs(a[2] - b[2]) for a in pts, b in pts if a !== b) > 11
    @test length(L.leader) == 2 * 2      # two leaders, two endpoints each

    # Each label clears its own element rather than printing over it, and reads
    # away from the centerline: left-aligned text sits to the right of the
    # anchor and vice versa.
    for k in eachindex(L.pos)
        all(isfinite, L.pos[k]) || continue
        @test hypot(L.offset[k]...) > 0
        @test (L.align[k][1] === :left) == (L.offset[k][1] >= 0)
    end

    # The property the old data-space stacking did not have: it separated labels
    # along the element's own vertical, which projects to nothing when you look
    # straight down it -- and straight down is the default view's own axis. In
    # screen space there is no camera angle that collapses the separation.
    for (el, az) in ((pi / 2 - 1.0f-3, 0.0), (0.0, 0.0), (0.3, 1.9), (-0.7, -2.5))
        fp.axis.elevation[] = el
        fp.axis.azimuth[] = az
        Makie.update_state_before_display!(fp.figure)
        Lr = pp._layout_labels3(fp.geometry, fp.axis.scene, 4, 11)
        @test count(p -> all(isfinite, p), Lr.pos) == 3
        @test allunique(screen_at(fp, Lr))
    end
end

@testset "3D labels: level of detail and label_sep" begin
    smap = ShapeMap([ele_shape("Widget", :box, :black; size=0.3, label=:name)];
                    defaults=false)
    co = synth_table(names=["pue_a12", "dhca12", "dvca12"],
                     kinds=fill("Widget", 3), lengths=zeros(3), x=zeros(3),
                     z=fill(10.0, 3), theta=zeros(3), angle=zeros(3))

    fp = floor_plot3(co; shapes=smap, size=(900, 700))
    Makie.update_state_before_display!(fp.figure)
    drawn(L) = count(p -> all(isfinite, p), L.pos)

    # An element smaller on screen than `label_min_px` is not labelled at all.
    @test drawn(pp._layout_labels3(fp.geometry, fp.axis.scene, 10_000, 11)) == 0

    # `label_sep = 0` turns the collision handling off: the labels are no longer
    # bumped apart, so only the first to claim the space is drawn.
    flat = floor_plot3(co; shapes=smap, size=(900, 700), label_sep=0)
    Makie.update_state_before_display!(flat.figure)
    Lf = pp._layout_labels3(flat.geometry, flat.axis.scene, 4, 11)
    @test flat.geometry.label_sep == 0
    @test drawn(Lf) == 1
    @test isempty(Lf.leader)

    # And a lattice whose elements are far apart needs no bumping at all.
    apart = synth_table(names=["a", "b"], kinds=["Widget", "Widget"],
                        lengths=[1.0, 1.0], x=[0.0, 0.0], z=[0.0, 40.0],
                        theta=zeros(2), angle=zeros(2))
    fa = floor_plot3(apart; shapes=smap, size=(900, 700))
    Makie.update_state_before_display!(fa.figure)
    La = pp._layout_labels3(fa.geometry, fa.axis.scene, 4, 11)
    @test drawn(La) == 2
    @test isempty(La.leader)
end

@testset "3D geometry: element_outline3 hugs the element" begin
    tab = synth_table(names=["e"], kinds=["Widget"], lengths=[2.0],
                      x=[0.0], z=[0.0], theta=[0.0], angle=[0.0])
    pts = element_outline3(tab, 1; view="zxy")
    @test !isempty(pts)
    @test iseven(length(pts))              # drawn as segment pairs
    @test all(p -> all(isfinite, p), pts)
    @test extrema(p -> p[1], pts) == (0.0f0, 2.0f0)
end

# Extraction from a real expanded lattice, if the sibling PALSJulia checkout with
# its example files is present (the same side-by-side layout the parser needs).
_lattice_file(name) = normpath(joinpath(@__DIR__, "..", "..", "PALSJulia",
                                        "lattice_files", name))
const CONVERT = _lattice_file("convert.pals.yaml")
const BTA = _lattice_file("bta.pals.yaml")
const FORK = _lattice_file("fork.pals.yaml")

# The centerline ends of element `i` as actually drawn, recovered from the
# public `element_outline`: it walks f = 0 -> 1 along the top edge and back along
# the bottom, both at the same half-height, so opposite pairs average to the
# centerline. Returns (entrance, exit).
function drawn_ends(tab, i; view="zx")
    p = element_outline(tab, i; view=view)
    n = (length(p) - 1) ÷ 2            # p = [top(1:n); reverse(bot); top[1]]
    return ((p[1] .+ p[2n]) ./ 2, (p[n] .+ p[n+1]) ./ 2)
end

# A reference curve that leaves the picture plane used to be drawn as though it
# did not: the drawing frame was built from the heading `theta` alone, so a
# pitched element came out at its full length instead of its projected length --
# overshooting by L(1 - cos(phi)), which on the 1.3 m cavity at phi = -0.63 in
# convert.pals.yaml was a quarter of a metre. The frame is now the standard's W
# matrix, so the projection is the real one.
@testset "an element pitched out of the plane is drawn foreshortened" begin
    L = 1.3
    for phi in (0.0, -0.63, 0.9)
        tab = synth_table(names=["c"], kinds=["RFCavity"], lengths=[L],
                          x=[0.0], z=[0.0], theta=[0.0], angle=[0.0], phi=[phi])
        entrance, exit = drawn_ends(tab, 1; view="zx")
        # Heading is along z, pitched by phi, so the "zx" projection of the
        # element's length is L*cos(phi) and none of it lands on the x axis.
        @test isapprox(exit[1] - entrance[1], L * cos(phi); atol=1e-5)
        @test isapprox(exit[2] - entrance[2], 0.0; atol=1e-5)

        # The part that left the plane has to show up in a view that contains
        # the vertical: in "zy" the same element rises by -L*sin(phi).
        _, exit_zy = drawn_ends(tab, 1; view="zy")
        entrance_zy, _ = drawn_ends(tab, 1; view="zy")
        @test isapprox(exit_zy[2] - entrance_zy[2], -L * sin(phi); atol=1e-5)
    end
end

# Roll shows up as the *width* of an element rather than its length: a magnet
# rolled about its own axis presents a narrower face to a viewer above it.
@testset "a rolled element is drawn narrower" begin
    smap = ShapeMap([ele_shape("Quadrupole", :box, :black; size=0.5)])
    for psi in (0.0, 0.5, π / 3)
        tab = synth_table(names=["q"], kinds=["Quadrupole"], lengths=[1.0],
                          x=[0.0], z=[0.0], theta=[0.0], angle=[0.0], psi=[psi])
        g = pp.build_geometry(tab, smap; view="zx")
        pts = filter(p -> !isnan(p[1]), g.outline_pts)
        halfwidth = maximum(p -> abs(p[2]), pts)
        @test isapprox(halfwidth, 0.5 * cos(psi); atol=1e-5)
    end
end

if isfile(CONVERT)
    @testset "extraction from expanded lattice" begin
        lat = parse_and_expand_pals(CONVERT; problems=:none)
        tab = element_table(lat)
        @test length(tab) > 0
        @test "drift1" in tab.name
        i = findfirst(==("drift1"), tab.name)
        @test isapprox(tab.length[i], 100.0; atol=1e-6)   # length: 1e+02
        @test tab.z[end] > tab.z[1]                        # positions advance
        @test "Bend" in tab.kind

        # The bookkeeper caps each branch with a zero-length `branch_end`
        # Placeholder carrying the downstream end of the last real element.
        @test tab.name[end] == "branch_end"
        @test tab.kind[end] == "Placeholder"
        @test tab.s[end] > tab.s[1]

        # The bend angle comes from BendP.angle_ref, which the bookkeeper
        # derives; it must agree with g_ref * length.
        b = findfirst(==("Bend"), tab.kind)
        @test tab.angle[b] != 0.0
        bend = tab.node[b]["BendP"]
        @test isapprox(tab.angle[b], Float64(bend["g_ref"]) * tab.length[b];
                       rtol=1e-9)
    end

    # Placement lives only in `full_expanded`; `expanded` is pruned back to what
    # the author wrote, so reading it would yield an empty table.
    @testset "placement comes from full_expanded" begin
        lat = parse_and_expand_pals(CONVERT; problems=:none)
        ele = lat.full_expanded["lat"]["branches"][1]["ring"]["line"][1]["beg"]
        @test pj.haskey(ele, "FloorP")
        @test pj.haskey(ele, "s_position")
        pruned = lat.expanded["lat"]["branches"][1]["ring"]["line"][1]["beg"]
        @test !pj.haskey(pruned, "FloorP")
    end
else
    @info "Skipping extraction tests: $CONVERT not found (needs sibling PALSJulia checkout)"
end

# The drawing has to agree with the floor coordinates the expander computed, not
# merely be self-consistent: each element is *placed* from its own FloorP, but
# its length, heading and bend angle are what carry the drawing across it. If any
# of the three is read wrong, the far end of one element stops meeting the near
# end of the next -- which is the one thing FloorP pins down independently.
#
# bta.pals.yaml is a real machine (159 elements, 20 bends) and lies in the
# horizontal plane, so "zx" and "xz" show it without foreshortening.
if isfile(BTA)
    @testset "drawn geometry closes on the expander's floor coordinates" begin
        lat = parse_and_expand_pals(BTA; problems=:none)
        tab = element_table(lat)
        @test length(tab) > 100
        @test count(==("Bend"), tab.kind) > 0

        for view in ("zx", "xz")
            worst = 0.0
            for i in 1:(length(tab) - 1)
                tab.branch[i] == tab.branch[i + 1] || continue
                _, exit = drawn_ends(tab, i; view=view)
                entrance, _ = drawn_ends(tab, i + 1; view=view)
                worst = max(worst, hypot(exit[1] - entrance[1],
                                         exit[2] - entrance[2]))
            end
            # Point2f is Float32; over a ~113 m lattice that is ~1e-5 of slack.
            @test worst < 1e-4
        end

        # The first element starts where FloorP says, in either view.
        e_zx, _ = drawn_ends(tab, 1; view="zx")
        @test isapprox(e_zx[1], tab.z[1]; atol=1e-4)
        @test isapprox(e_zx[2], tab.x[1]; atol=1e-4)
        e_xz, _ = drawn_ends(tab, 1; view="xz")
        @test isapprox(e_xz[1], tab.x[1]; atol=1e-4)
        @test isapprox(e_xz[2], tab.z[1]; atol=1e-4)
    end
end

# convert.pals.yaml is the lattice that does *not* stay in the horizontal plane:
# nine of its elements carry a pitch and a roll, and one bend has a tilt_ref. It
# is the case the drawing used to get wrong, so it gets checked the hard way --
# against the expander's own floor coordinates rather than for self-consistency.
if isfile(CONVERT)
    @testset "the drawing closes on floor coordinates out of the plane" begin
        lat = parse_and_expand_pals(CONVERT; problems=:none)
        tab = element_table(lat)

        # If this lattice ever becomes planar the rest of the testset stops
        # meaning anything, so say so here rather than passing vacuously.
        @test count(!=(0.0), tab.phi) > 0
        @test count(!=(0.0), tab.psi) > 0
        @test count(!=(0.0), tab.tilt_ref) > 0

        # Carrying an element's placement across itself must land on the next
        # element's FloorP -- position *and* orientation. Position alone would
        # pass on a drawing whose elements were all correctly placed but rolled.
        worst_r = 0.0; worst_W = 0.0
        for i in 1:(length(tab) - 1)
            tab.branch[i] == tab.branch[i + 1] || continue
            p = pp.place(tab, i, 1.0)
            nxt = pp.Vec3d(tab.x[i + 1], tab.y[i + 1], tab.z[i + 1])
            worst_r = max(worst_r, hypot((p.r - nxt)...))
            Wn = pp.w_matrix(tab.theta[i + 1], tab.phi[i + 1], tab.psi[i + 1])
            for f in fieldnames(pp.Mat3)
                worst_W = max(worst_W, abs(getfield(p.W, f) - getfield(Wn, f)))
            end
        end
        @test worst_r < 1e-9
        @test worst_W < 1e-9

        # ...and the drawing itself, not just the math behind it: the drawn exit
        # of each element is the drawn entrance of the next.
        g = build_geometry3(tab, ShapeMap(); view="zxy")
        worst = 0.0
        for i in 1:(length(tab) - 1)
            tab.branch[i] == tab.branch[i + 1] || continue
            worst = max(worst, hypot((g.ele_exit[i] - g.ele_entrance[i + 1])...))
        end
        @test worst < 1e-3          # Point3f over a ~140 m lattice

        # Every 2D projection closes too, including the two that show the
        # vertical -- which is where drawing from `theta` alone came apart.
        for view in ("zx", "zy", "xy")
            worst2 = 0.0
            for i in 1:(length(tab) - 1)
                tab.branch[i] == tab.branch[i + 1] || continue
                _, exit = drawn_ends(tab, i; view=view)
                entrance, _ = drawn_ends(tab, i + 1; view=view)
                worst2 = max(worst2, hypot(exit[1] - entrance[1],
                                           exit[2] - entrance[2]))
            end
            @test worst2 < 1e-3
        end
    end
end

# A fork spawns a new branch, so the table spans several lines rather than one.
if isfile(FORK)
    @testset "extraction across forked branches" begin
        lat = parse_and_expand_pals(FORK; problems=:none)
        tab = element_table(lat)
        @test length(tab.branch_names) > 1
        @test "Fork" in tab.kind
        @test allunique(tab.branch_names)

        # Every element names a branch that exists, and branches are appended
        # whole rather than interleaved.
        @test all(1 .<= tab.branch .<= length(tab.branch_names))
        @test issorted(tab.branch)

        # Each branch runs from its own beginning and is capped by the
        # bookkeeper's branch_end.
        for b in 1:length(tab.branch_names)
            idx = findall(==(b), tab.branch)
            @test !isempty(idx)
            @test tab.name[last(idx)] == "branch_end"
            @test issorted(tab.s[idx])
        end

        # The reference orbit is broken between branches so the polyline does
        # not jump from the end of one line to the start of the next.
        g = pp.build_geometry(tab, ShapeMap(); view="zx")
        @test count(p -> isnan(p[1]), g.ref_pts) == length(tab.branch_names) - 1
    end
end

# The render stage, up to but not including opening a window: building the
# figure exercises every Makie call the plotter makes, and is where a Makie or
# GLMakie change that PALSPlot has not followed would surface.
@testset "floor_plot builds a figure" begin
    tab = synth_table(names=["q1", "b1"], kinds=["Quadrupole", "Bend"],
                      lengths=[0.5, 1.0], x=[0.0, 0.0], z=[0.0, 0.5],
                      theta=[0.0, 0.0], angle=[0.0, 0.2])
    fp = floor_plot(tab; view="zx", title="test")
    @test fp isa FloorPlot
    @test fp.axis.title[] == "test"
    @test fp.selected[] == 0
    @test length(fp.table) == 2
    @test !isempty(fp.geometry.outline_pts)

    # Selecting an element drives the highlight overlay.
    fp.selected[] = 2
    @test !isempty(element_outline(fp.table, 2; view="zx"))
end

# The mouse bindings the docs promise. These are registered `Axis` interactions,
# so they can be driven headlessly by handing `process_interaction` a synthetic
# MouseEvent -- no window, no backend, no real cursor. Without this the click and
# double-click paths are only ever exercised by hand.
@testset "mouse bindings" begin
    tab = synth_table(names=["q1", "b1"], kinds=["Quadrupole", "Bend"],
                      lengths=[0.5, 1.0], x=[0.0, 0.0], z=[0.0, 0.5],
                      theta=[0.0, 0.0], angle=[0.0, 0.2])
    fp = floor_plot(tab; view="zx")
    ax = fp.axis

    # Selection is registered under the name the source registers it by, next to
    # Makie's own bindings (the reset among them).
    @test haskey(Makie.interactions(ax), :selectelement)
    @test haskey(Makie.interactions(ax), :limitreset)

    mev(type, pos) = Makie.MouseEvent(type, 0.0, Point2d(pos), Point2f(0, 0),
                                      0.0, Point2d(pos), Point2f(0, 0))
    fire(type, pos) = Makie.process_axis_event(ax, mev(type, pos))
    ctrl!(down) = (down ? push! : delete!)(events(ax.scene).keyboardstate,
                                           Keyboard.left_control)

    # A left-click near an element selects it; the second element sits near
    # (z, x) = (1.0, 0.0), the first near (0.25, 0.0).
    fire(Makie.MouseEventTypes.leftclick, fp.geometry.ele_center[2])
    @test fp.selected[] == 2
    fire(Makie.MouseEventTypes.leftclick, fp.geometry.ele_center[1])
    @test fp.selected[] == 1

    # Zoom the way the rubber-band and scroll interactions do, by moving
    # `targetlimits`. (Not `limits!`, which sets the axis's *stored* limits --
    # what a reset resets back *to*, so a reset would then be a no-op.)
    zoom_in!() = (ax.targetlimits[] = Rect2(0.4, -0.1, 0.2, 0.2))

    # Ctrl-left-click is Makie's reset. It must not also select: `ax.interactions`
    # is an unordered Dict, so this only holds because selection checks the
    # modifier itself rather than trusting that the reset ran first.
    zoom_in!()
    zoomed = ax.finallimits[]
    ctrl!(true)
    fire(Makie.MouseEventTypes.leftclick, fp.geometry.ele_center[2])
    ctrl!(false)
    @test ax.finallimits[] != zoomed
    @test fp.selected[] == 1

    # Without the modifier the same click selects and leaves the view alone.
    zoom_in!()
    zoomed = ax.finallimits[]
    fire(Makie.MouseEventTypes.leftclick, fp.geometry.ele_center[2])
    @test ax.finallimits[] == zoomed
    @test fp.selected[] == 2
end

# Building a figure is not the same as being able to draw one: rasterizing it
# runs every plot object through a backend, which is where a Makie change that
# PALSPlot has not followed would actually bite.
@testset "the figure renders through a backend" begin
    tab = synth_table(names=["q1", "b1"], kinds=["Quadrupole", "Bend"],
                      lengths=[0.5, 1.0], x=[0.0, 0.0], z=[0.0, 0.5],
                      theta=[0.0, 0.0], angle=[0.0, 0.2])
    fp = floor_plot(tab; view="zx", size=(400, 300))

    img = Makie.colorbuffer(fp.figure)
    @test size(img, 1) > 0 && size(img, 2) > 0
    @test length(unique(img)) > 1        # something was actually drawn, not a blank

    mktempdir() do dir
        for ext in ("png", "pdf", "svg")
            path = joinpath(dir, "floor.$ext")
            save(path, fp.figure)
            @test isfile(path)
            @test filesize(path) > 0
        end
    end
end

if isfile(CONVERT)
    @testset "floor_plot from a Lattices, and the parameter panel" begin
        lat = parse_and_expand_pals(CONVERT; problems=:none)
        fp = floor_plot(lat; view="zx")
        @test length(fp.table) > 0
        @test length(fp.geometry.ele_center) == length(fp.table)

        # The side panel lists an element's parameters by walking its node.
        b = findfirst(==("Bend"), fp.table.kind)
        text = pp._node_text(fp.table.node[b])
        @test occursin("kind: Bend", text)
        @test occursin("BendP:", text)
        @test occursin("angle_ref:", text)   # a derived value, so full_expanded
    end
end

# ── the 3D plotter ────────────────────────────────────────────────────────────

@testset "floor_plot3 builds a figure" begin
    tab = synth_table(names=["q1", "b1"], kinds=["Quadrupole", "Bend"],
                      lengths=[0.5, 1.0], x=[0.0, 0.0], z=[0.0, 0.5],
                      theta=[0.0, 0.0], angle=[0.0, 0.2])
    fp = floor_plot3(tab; title="test3")
    @test fp isa FloorPlot3
    @test fp.axis isa Makie.Axis3
    @test fp.axis.title[] == "test3"
    @test fp.selected[] == 0
    @test length(fp.table) == 2
    @test !isempty(fp.geometry.mesh_faces)

    # A machine is long and thin, so an axis that normalized its three
    # dimensions into one box would fatten it by that whole ratio.
    @test fp.axis.aspect[] === :data

    # The handle carries the same contract as the 2D one, which is what lets a
    # script drive either.
    @test fp.view == "zxy"
    fp.selected[] = 2
    @test !isempty(element_outline3(fp.table, 2; view=fp.view))
end

# Building a figure is not the same as being able to draw one: rasterizing runs
# every plot object through a backend. CairoMakie sorts whole primitives rather
# than pixels, so a big 3D mesh comes out with artifacts -- that is a reason not
# to look at machines this way, not a reason to leave the path untested.
@testset "the 3D figure renders through a backend" begin
    tab = synth_table(names=["q1", "b1"], kinds=["Quadrupole", "Bend"],
                      lengths=[0.5, 1.0], x=[0.0, 0.0], z=[0.0, 0.5],
                      theta=[0.0, 0.0], angle=[0.0, 0.2])
    fp = floor_plot3(tab; size=(400, 300))

    img = Makie.colorbuffer(fp.figure)
    @test size(img, 1) > 0 && size(img, 2) > 0
    @test length(unique(img)) > 1        # something was drawn, not a blank

    mktempdir() do dir
        for ext in ("png", "pdf", "svg")
            path = joinpath(dir, "floor3.$ext")
            save(path, fp.figure)
            @test isfile(path)
            @test filesize(path) > 0
        end
    end
end

# Click-to-inspect in 3D. A cursor in a 3D axis is a ray, not a point, so this
# does not work the way the floor plan's picking does and needs its own test.
# There is no window here: the click is a synthetic `MouseEvent` at the pixel
# `project` says the element sits at, which is exactly what a real click carries.
@testset "3D mouse bindings" begin
    tab = synth_table(names=["q1", "b1", "s1"],
                      kinds=["Quadrupole", "Bend", "Sextupole"],
                      lengths=[0.5, 1.0, 0.4], x=[0.0, 0.0, 0.0],
                      z=[0.0, 2.0, 6.0], theta=zeros(3), angle=[0.0, 0.0, 0.0])
    fp = floor_plot3(tab; size=(800, 600))
    ax = fp.axis
    Makie.update_state_before_display!(fp.figure)

    # Registered next to Makie's own Axis3 bindings, the reset among them.
    @test haskey(Makie.interactions(ax), :selectelement)
    @test haskey(Makie.interactions(ax), :limitreset)
    @test haskey(Makie.interactions(ax), :dragrotate)

    mev(type, px) = Makie.MouseEvent(type, 0.0, Point2d(0, 0), Point2f(px),
                                     0.0, Point2d(0, 0), Point2f(px))
    click(px) = Makie.process_interaction(Makie.interactions(ax)[:selectelement][2],
                                          mev(Makie.MouseEventTypes.leftclick, px), ax)
    at(i) = Makie.project(ax.scene, fp.geometry.ele_center[i])
    ctrl!(down) = (down ? push! : delete!)(events(ax.scene).keyboardstate,
                                           Keyboard.left_control)

    for i in (2, 1, 3)
        click(at(i))
        @test fp.selected[] == i
    end

    # Ctrl-left-click is Makie's view reset, not a selection. `ax.interactions`
    # is an unordered Dict, so this holds only because selection checks the
    # modifier itself rather than trusting that the reset ran first.
    click(at(1))
    ctrl!(true)
    click(at(3))
    ctrl!(false)
    @test fp.selected[] == 1

    # A click far from the machine selects nothing rather than the least-distant
    # element in the lattice.
    fp.selected[] = 0
    click(Point2f(-4000, -4000))
    @test fp.selected[] == 0
end

# Overlays are ordinary plots on the existing axis, and take their input in
# global coordinates whichever drawing they go on.
@testset "overlays" begin
    tab = synth_table(names=["q1", "b1"], kinds=["Quadrupole", "Bend"],
                      lengths=[0.5, 1.0], x=[0.0, 0.0], z=[0.0, 0.5],
                      theta=[0.0, 0.0], angle=[0.0, 0.2])
    curve = [(0.0, 0.0, 0.0), (0.1, 0.2, 1.0), (0.0, 0.0, 2.0)]
    wall = [(-2.0, 0.0, -1.0), (-2.0, 0.0, 3.0), (2.0, 0.0, 3.0)]

    fp = floor_plot(tab)
    n0 = length(fp.axis.scene.plots)
    c = add_curve!(fp, curve; color=:orange)
    w = add_wall!(fp, wall)
    @test length(fp.axis.scene.plots) == n0 + 2
    # Projected the way the plot was built: "zx" puts global z on the horizontal.
    @test c[1][][2] ≈ Point2f(1.0, 0.1)

    fp3 = floor_plot3(tab)
    n3 = length(fp3.axis.scene.plots)
    c3 = add_curve!(fp3, curve; color=:orange)
    w3 = add_wall!(fp3, wall; base=-1.0, height=4.0)
    @test length(fp3.axis.scene.plots) == n3 + 2
    # "zxy" puts global z on the drawn x, global x on the drawn y, y up.
    @test c3[1][][2] ≈ Point3f(1.0, 0.1, 0.2)
    # The wall is a solid between base and base + height on the drawn vertical.
    @test extrema(p -> p[3], w3[1][].position) == (-1.0f0, 3.0f0)

    # Two-vectors are read as (x, z) on the y = 0 plane.
    c2 = add_curve!(fp3, [(0.0, 0.0), (1.0, 5.0)])
    @test c2[1][][2] ≈ Point3f(5.0, 1.0, 0.0)

    # Both figures still render with the overlays on them.
    @test length(unique(Makie.colorbuffer(fp.figure))) > 1
    @test length(unique(Makie.colorbuffer(fp3.figure))) > 1
end

if isfile(BTA)
    @testset "floor_plot3 from a real lattice" begin
        lat = parse_and_expand_pals(BTA; problems=:none)
        fp = floor_plot3(lat)
        @test length(fp.table) > 100
        @test length(fp.geometry.ele_center) == length(fp.table)
        @test !isempty(fp.geometry.mesh_faces)

        # The whole machine is a fixed number of plot objects, not one per
        # element: that is what keeps a big lattice drawable.
        @test length(fp.axis.scene.plots) < 10

        # Picking maps a mesh vertex back to the element that owns it.
        @test length(fp.geometry.vertex_ele) == length(fp.geometry.mesh_pts)
        @test extrema(fp.geometry.vertex_ele) ⊆ 1:length(fp.table)
    end
end
