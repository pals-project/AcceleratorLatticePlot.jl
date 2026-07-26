# Draw a PALS lattice in three dimensions.
#
# Run with:
#     julia --project=. examples/floor_plan_3d.jl [lattice-file]
#
# in an environment with GLMakie added. A window opens: left-drag to orbit,
# right-drag to pan, scroll to zoom, left-click an element to list its parameters
# in the side panel.
#
# The 2D floor plan is examples/floor_plan.jl; the two share their shape table
# and their placement math, so the same machine looks the same in both -- look
# straight down at this one and you have the floor plan.

using PALSJulia
using PALSPlot
using GLMakie      # the backend; PALSPlot itself depends only on Makie

# Default to the bundled convert example in the sibling PALSJulia checkout: it is
# the one that does not stay in the horizontal plane, so it has something to show
# in 3D that a floor plan cannot.
default_file = normpath(joinpath(@__DIR__, "..", "..", "PALSJulia",
                                 "lattice_files", "convert.pals.yaml"))
file = isempty(ARGS) ? default_file : ARGS[1]

println("Parsing $file ...")
lat = parse_and_expand_pals(file; problems=:none)
isempty(lat.problems) || println("(", length(lat.problems), " expansion problems ignored)")

fp = floor_plot3(lat; title=basename(file))

# A machine that lies flat is a hundred metres of one-metre magnets, and at the
# honest `aspect = :data` its box comes out a pancake. Uncomment to stretch the
# vertical for a diagram rather than a picture of the real thing:
#     fp.axis.aspect[] = (1, 1, 0.4)

screen = display(fp)      # opens the GLMakie window

println("""
Controls: left-drag = orbit,  right-drag = pan,  scroll = zoom,
          ctrl-left-click = reset view,  left-click = select an element.""")

# Keep the process alive while the window is open when run as a script.
if !isinteractive()
    println("Close the window to exit.")
    try
        wait(screen)
    catch
        println("Press Enter to exit."); readline()
    end
end
