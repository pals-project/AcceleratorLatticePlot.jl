# The GLMakie path, which the main suite deliberately does not cover.
#
# AcceleratorLatticePlot depends on Makie rather than on a backend, and
# `runtests.jl` renders through CairoMakie because that works on any machine.
# GLMakie is what an actual user opens a window with, though, so it needs
# checking somewhere -- but it needs a real OpenGL context, which a headless
# macOS CI runner cannot provide at all (GLMakie will not even precompile there).
# So this lives outside the suite and is run only where a context is available: a
# developer's machine, or Linux CI under xvfb.
#
# Run with GLMakie available:
#     julia --project=. test/glmakie_smoke.jl

using Test
using GLMakie
using AcceleratorLatticePlot
import AcceleratorLatticePlot as alp
import PALSParserJ as pj

const NODE = pj.parse_string("kind: Drift\n")

@testset "GLMakie renders a floor plan" begin
    tab = alp.ElementTable(["q1", "b1"], ["Quadrupole", "Bend"], [0.5, 1.0],
                          [0.0, 0.0], [0.0, 0.0], [0.0, 0.5], [0.0, 0.0],
                          [0.0, 0.2], [0.0, 0.5], [1, 1], fill(NODE, 2), ["b"])

    @test Makie.current_backend() === GLMakie

    fp = floor_plot(tab; view="zx", title="smoke", size=(400, 300))
    @test fp isa FloorPlot

    # Rasterizing goes all the way through the GL pipeline, which is the part
    # CairoMakie cannot stand in for.
    img = Makie.colorbuffer(fp.figure)
    @test size(img, 1) > 0 && size(img, 2) > 0
    @test length(unique(img)) > 1

    # And the interaction wiring the window relies on is present.
    @test fp.selected[] == 0
    fp.selected[] = 2
    @test fp.selected[] == 2
end

# The 3D drawing is really a GLMakie feature: CairoMakie will render one, but it
# sorts whole primitives rather than pixels, so a machine-sized mesh comes out
# wrong. This is also the only place GPU picking can be exercised -- it reads a
# framebuffer, which is the one thing a software backend has not got.
@testset "GLMakie renders a 3D drawing, and picks from it" begin
    #                     name          kind          length  x  y  z  theta
    #                     angle  tilt   s  branch  node  branch_names
    tab = alp.ElementTable(["q1", "b1", "s1"], ["Quadrupole", "Bend", "Sextupole"],
                          [0.5, 1.0, 0.4], zeros(3), zeros(3), [0.0, 2.0, 6.0],
                          zeros(3), zeros(3), [0.0, 2.0, 6.0], ones(Int, 3),
                          fill(NODE, 3), ["b"])

    fp = floor_plot3(tab; title="smoke 3D", size=(600, 450))
    @test fp isa FloorPlot3
    @test !isempty(fp.geometry.mesh_faces)

    # A pick reads a rendered framebuffer, so this needs a real screen, and the
    # figure has to be drawn on it before anything can be asked of it. Rendering
    # goes through that same screen rather than `colorbuffer(figure)`, which
    # would open a second one -- GLMakie allows a figure only one.
    screen = GLMakie.Screen(; visible=false)
    display(screen, fp.figure)

    img = Makie.colorbuffer(screen)
    @test size(img, 1) > 0 && size(img, 2) > 0
    @test length(unique(img)) > 1

    for i in 1:3
        px = Makie.project(fp.axis.scene, fp.geometry.ele_center[i])
        @test alp._pick3(fp.geometry, fp.axis, px) == i
    end
end

# Labels are laid out against the live camera, and their leader lines are drawn
# in `space = :pixel`. Both go through the backend's projection rather than
# through anything AcceleratorLatticePlot computes itself, so a backend is the
# only place they can be checked -- and GLMakie is the one this view is meant
# for.
@testset "GLMakie lays out 3D labels against a live camera" begin
    # Three elements on the same point, the case that forces the layout to
    # separate labels and to draw leaders back to what they name.
    tab = alp.ElementTable(["pue_a12", "dhca12", "dvca12"],
                          ["Instrument", "Kicker", "Kicker"], zeros(3),
                          zeros(3), zeros(3), fill(10.0, 3), zeros(3), zeros(3),
                          fill(10.0, 3), ones(Int, 3), fill(NODE, 3), ["b"])

    fp = floor_plot3(tab; title="smoke labels", size=(900, 700))
    screen = GLMakie.Screen(; visible=false)
    display(screen, fp.figure)
    Makie.colorbuffer(screen)

    L = alp._layout_labels3(fp.geometry, fp.axis.scene, 4, 11)
    @test count(p -> all(isfinite, p), L.pos) == 3
    @test length(L.leader) == 4          # two leaders, two endpoints each

    # The leaders reach the rasterizer, rather than merely being computed:
    # hiding them has to change the picture. `space = :pixel` is the one thing
    # here that CairoMakie cannot vouch for on GLMakie's behalf.
    leaders = only(p for p in fp.axis.scene.plots
                     if p isa Makie.LineSegments && to_value(p.space) === :pixel)
    @test length(leaders[1][]) == 4
    lit = copy(Makie.colorbuffer(screen))
    leaders.visible[] = false
    @test count(lit .!= Makie.colorbuffer(screen)) > 0
    leaders.visible[] = true

    # No two land on the same pixels, from any camera angle.
    for (el, az) in ((pi / 2 - 1.0f-3, 0.0), (0.3, 1.9), (-0.7, -2.5))
        fp.axis.elevation[] = el
        fp.axis.azimuth[] = az
        Makie.colorbuffer(screen)
        Lr = alp._layout_labels3(fp.geometry, fp.axis.scene, 4, 11)
        at = [Makie.project(fp.axis.scene, fp.geometry.label_pos[k]) .+ Lr.offset[k]
              for k in eachindex(Lr.pos) if all(isfinite, Lr.pos[k])]
        @test length(at) == 3
        @test allunique(at)
    end
end
