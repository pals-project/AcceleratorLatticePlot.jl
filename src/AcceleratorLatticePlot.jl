"""
    AcceleratorLatticePlot

Floor-plan plotting for expanded PALS lattices (see
[`PALSParserJ`](https://github.com/pals-project/PALSParserJ.jl)).

Given the `Lattices` object returned by `PALSParserJ.parse_and_expand_pals`,
draw the machine projected onto a plane ([`floor_plot`](@ref)) or as solids in
three dimensions ([`floor_plot3`](@ref)), with each element rendered as a shape
sized, colored and labeled by a Tao-style rule table. Both windows support pan,
zoom, and click-to-inspect: clicking an element lists its full parameter set.

Quick start:

```julia
using PALSParserJ, AcceleratorLatticePlot
using GLMakie                 # or any other Makie backend
lat = parse_and_expand_pals("machine.pals.yaml")
fp = floor_plot(lat)          # 2D floor plan
display(fp)                   # opens an interactive window

fp3 = floor_plot3(lat)        # ...the same machine in 3D
display(fp3)
```

AcceleratorLatticePlot builds on `Makie` rather than on a particular backend, so
the caller chooses one: GLMakie for an interactive window, CairoMakie to write a
file. None is needed to build the figure, or to use the extraction and geometry
stages.

The geometry, shape mapping and rendering stages are separate so the geometry can
also be used headless (see [`element_table`](@ref) and [`build_geometry`](@ref)).
"""
module AcceleratorLatticePlot

include("lattice.jl")
include("frame.jl")
include("shapes.jl")
include("geometry.jl")
include("geometry3.jl")
include("render.jl")
include("render3.jl")
include("overlays.jl")

export floor_plot, FloorPlot
export floor_plot3, FloorPlot3
export element_table, ElementTable
export ShapeMap, ShapeRule, ShapeSpec, ele_shape, mapshape, default_shapes
export build_geometry, FloorGeometry, element_outline
export build_geometry3, FloorGeometry3, element_outline3
export w_matrix, place, Placement
export add_curve!, add_wall!

end # module
