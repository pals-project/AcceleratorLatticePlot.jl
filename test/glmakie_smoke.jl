# The GLMakie path, which the main suite deliberately does not cover.
#
# PALSPlot depends on Makie rather than on a backend, and `runtests.jl` renders
# through CairoMakie because that works on any machine. GLMakie is what an actual
# user opens a window with, though, so it needs checking somewhere -- but it
# needs a real OpenGL context, which a headless macOS CI runner cannot provide at
# all (GLMakie will not even precompile there). So this lives outside the suite
# and is run only where a context is available: a developer's machine, or Linux
# CI under xvfb.
#
# Run with GLMakie available:
#     julia --project=. test/glmakie_smoke.jl

using Test
using GLMakie
using PALSPlot
import PALSPlot as pp
import PALSJulia as pj

const NODE = pj.parse_string("kind: Drift\n")

@testset "GLMakie renders a floor plan" begin
    tab = pp.ElementTable(["q1", "b1"], ["Quadrupole", "Bend"], [0.5, 1.0],
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
    tab = pp.ElementTable(["q1", "b1", "s1"], ["Quadrupole", "Bend", "Sextupole"],
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
        @test pp._pick3(fp.geometry, fp.axis, px) == i
    end
end
