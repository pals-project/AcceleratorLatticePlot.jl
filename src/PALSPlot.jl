"""
    PALSPlot

Floor-plan plotting for expanded PALS lattices (see
[`PALSJulia`](https://github.com/pals-project/PALSJulia.jl)).

Given the `Lattices` object returned by `PALSJulia.parse_and_expand_pals`, draw
the machine projected onto a plane, with each element rendered as a shape sized,
colored and labeled by a Tao-style rule table. The window supports pan, zoom, and
click-to-inspect: clicking an element lists its full parameter set.

Quick start:

```julia
using PALSJulia, PALSPlot
using GLMakie                 # or any other Makie backend
lat = parse_and_expand_pals("machine.pals.yaml")
fp = floor_plot(lat)
display(fp)                   # opens an interactive window
```

PALSPlot builds on `Makie` rather than on a particular backend, so the caller
chooses one: GLMakie for an interactive window, CairoMakie to write a file. None
is needed to build the figure, or to use the extraction and geometry stages.

The geometry, shape mapping and rendering stages are separate so the geometry can
also be used headless (see [`element_table`](@ref) and [`build_geometry`](@ref)).
"""
module PALSPlot

include("lattice.jl")
include("shapes.jl")
include("geometry.jl")
include("render.jl")

export floor_plot, FloorPlot
export element_table, ElementTable
export ShapeMap, ShapeRule, ShapeSpec, ele_shape, mapshape, default_shapes
export build_geometry, FloorGeometry, element_outline

end # module
