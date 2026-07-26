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
