# Draw an interactive floor plan of a PALS lattice.
#
# Run with:
#     julia --project=. examples/floor_plan.jl [lattice-file]
#
# in an environment with GLMakie added -- AcceleratorLatticePlot depends on Makie
# and leaves the backend to the caller. A window opens: scroll to zoom,
# right-drag to pan, left-click an element to list its parameters in the side
# panel.

using PALSParserJ
using AcceleratorLatticePlot
using GLMakie      # the backend; this package depends only on Makie

# Default to the bundled convert example in the sibling PALSParserJ checkout.
default_file = normpath(joinpath(@__DIR__, "..", "..", "PALSParserJ",
                                 "lattice_files", "convert.pals.yaml"))
file = isempty(ARGS) ? default_file : ARGS[1]

println("Parsing $file ...")
lat = parse_and_expand_pals(file; problems=:none)
isempty(lat.problems) || println("(", length(lat.problems), " expansion problems ignored)")

fp = floor_plot(lat; view="zx", title=basename(file), label_min_px=6)
screen = display(fp)      # opens the GLMakie window

println("""
Controls: scroll = zoom,  right-drag = pan,  left-drag = rubber-band zoom,
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
