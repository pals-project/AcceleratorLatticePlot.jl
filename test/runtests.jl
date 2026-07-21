using Test
using PALSPlot
using PALSJulia          # for parse_and_expand_pals in the extraction tests
import PALSPlot as pp
import PALSJulia as pj

# These tests never open a window, so GLMakie is never loaded and the pals-cpp
# parser stays safe to call (see the note in render.jl / README).

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
    @test mapshape(m, "SBend", "b1").shape == :box
    @test mapshape(m, "Sextupole", "s1").shape == :box
    @test mapshape(m, "Drift", "d1").shape == :none
    # unknown kind falls through to the "*" default
    @test mapshape(m, "Wibbler", "w1").shape == :box
    @test mapshape(m, "Wibbler", "w1").label == :none
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
        synth_table(names=["b"], kinds=["SBend"], lengths=[1.0],
                    x=[0.0], z=[0.0], theta=[0.0], angle=[0.0]),
        ShapeMap(); view="zx")
    bent = pp.build_geometry(
        synth_table(names=["b"], kinds=["SBend"], lengths=[1.0],
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

# Extraction from a real expanded lattice, if the sibling PALSJulia checkout with
# its example files is present (the same side-by-side layout the parser needs).
const CONVERT = normpath(joinpath(@__DIR__, "..", "..", "PALSJulia",
                                  "lattice_files", "convert.pals.yaml"))
if isfile(CONVERT)
    @testset "extraction from expanded lattice" begin
        lat = parse_and_expand_pals(CONVERT; problems=:none)
        tab = element_table(lat)
        @test length(tab) > 0
        @test "drift1" in tab.name
        i = findfirst(==("drift1"), tab.name)
        @test isapprox(tab.length[i], 100.0; atol=1e-6)   # length: 1e+02
        @test tab.z[end] > tab.z[1]                        # positions advance
        @test "SBend" in tab.kind
    end
else
    @info "Skipping extraction tests: $CONVERT not found (needs sibling PALSJulia checkout)"
end
