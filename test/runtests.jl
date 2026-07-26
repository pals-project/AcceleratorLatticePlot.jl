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

# A throwaway node to fill ElementTable.node in synthetic tables.
const NODE = pj.parse_string("kind: Drift\n")

# Build a one-branch ElementTable by hand.
function synth_table(; names, kinds, lengths, x, z, theta, angle)
    n = length(names)
    return pp.ElementTable(collect(names), collect(kinds), collect(lengths),
                           collect(x), zeros(n), collect(z), collect(theta),
                           collect(angle), collect(z), ones(Int, n),
                           fill(NODE, n), ["b"])
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
        gx, gy, = pp._frame('z', 'x', th)
        along = atan(gy, gx)
        @test isapprox(mod(rot - along, π), π / 2; atol = 1e-5)
    end

    # A label starts at its anchor and runs outward, away from the centerline,
    # whichever way the readability fold left it pointing.
    for th in (0.0, π)          # normal points +x on one, -x on the other
        tab = synth_table(names=["q1"], kinds=["Quadrupole"], lengths=[1.0],
                          x=[0.0], z=[0.0], theta=[th], angle=[0.0])
        g = pp.build_geometry(tab, smap; view="zx")
        _, _, nx, ny = pp._frame('z', 'x', th)
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
    _, _, nx, ny = pp._frame('z', 'x', 0.0)
    for k in 2:3
        @test offs[k][1] * nx + offs[k][2] * ny > 0      # outward, not inward
        @test hypot(offs[k]...) ≈ 11 * g.label_stack[k] rtol = 1e-5
    end
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
