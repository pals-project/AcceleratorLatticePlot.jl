# PALSPlot.jl

Floor-plan plotting for [PALS](https://github.com/campa-consortium/pals) lattices,
built on [PALSJulia](https://github.com/pals-project/PALSJulia.jl) and
[GLMakie](https://docs.makie.org/stable/).

Given the expanded lattice from `PALSJulia.parse_and_expand_pals`, PALSPlot draws
the machine projected onto a plane. Each element is rendered as a shape — sized,
colored and labeled by a [Tao](https://www.classe.cornell.edu/bmad/)-style rule
table — placed and oriented from the floor coordinates the expander computes.
Bends follow their true arc. The window supports **pan**, **zoom**, and
**click-to-inspect**: click an element and its full parameter set is listed in a
side panel.

It is designed to stay responsive on machines with tens of thousands of elements:
the entire lattice is drawn with a handful of batched Makie calls (one for all
element outlines, one for all interior strokes, one for the reference orbit, one
for the labels) rather than one plot object per element, and element labels use
level-of-detail so they appear only when the view is zoomed in far enough.

## Installation

PALSPlot depends on PALSJulia, which is a wrapper around the `yaml_c_wrapper` C
library shipped with [pals-cpp](https://github.com/pals-project/pals-cpp). Clone
all three side by side and build the C library first:

```
some-dir/
  pals-cpp/       # build this first: cmake -S . -B build && cmake --build build
  PALSJulia/
  PALSPlot/
```

Then, from the `PALSPlot` directory:

```console
julia --project=. -e 'using Pkg; Pkg.develop(path="../PALSJulia"); Pkg.instantiate()'
```

## Quick start

```julia
using PALSJulia
using PALSPlot

lat = parse_and_expand_pals("machine.pals.yaml")
fp  = floor_plot(lat)    # opens an interactive GLMakie window
display(fp)
```

or run the bundled example:

```console
julia --project=. examples/floor_plan.jl path/to/machine.pals.yaml
```

**Controls** (Makie's default `Axis` bindings)

| gesture | action |
|---|---|
| scroll wheel | zoom in/out at the cursor |
| **right-drag** | pan |
| left-drag | rubber-band zoom to a rectangle |
| double-click | reset to the full view |
| left-click | select an element (lists its parameters in the side panel) |

### Build pals-cpp against libc++

PALSPlot needs the pals-cpp C library to be built with the **same C++ runtime as
Julia and GLMakie — LLVM libc++ (Apple clang)**. If pals-cpp is instead built
with GCC/libstdc++ (e.g. because `CC`/`CXX` point at MacPorts GCC), then loading
GLMakie and then parsing a lattice file **aborts the process** (signal 6): the two
C++ exception runtimes clash in the library's file reader.

Build pals-cpp with clang:

```console
cd pals-cpp && rm -rf build
env -u CC -u CXX cmake -S . -B build \
    -DCMAKE_C_COMPILER=/usr/bin/clang -DCMAKE_CXX_COMPILER=/usr/bin/clang++
cmake --build build
```

Verify with `otool -L build/libyaml_c_wrapper.dylib`: it should show
`/usr/lib/libc++.1.dylib` and **no** `libstdc++`.

## Choosing shapes

Element appearance is controlled by an ordered list of rules, each mapping an
`"<class>::<name-glob>"` to a shape, color, size and label — the first rule an
element matches wins, with a built-in default table appended:

```julia
shapes = ShapeMap([
    ele_shape("Quadrupole::q*", :xbox,   :black; size=0.6, label=:name),
    ele_shape("SBend",          :box,    :blue;  size=0.5),
    ele_shape("Wiggler",        :box,    :green; size=0.5),
])
floor_plot(lat; shapes=shapes)
```

Available shapes: `:box`, `:xbox`, `:x`, `:bow_tie`, `:rbow_tie`, `:diamond`,
`:u_triangle`, `:d_triangle`, `:l_triangle`, `:r_triangle`, `:circle`, `:none`.
`size` is the transverse half-height in meters; `label` is `:name`, `:s`, or
`:none`. Pass `defaults=false` to `ShapeMap` to use only your own rules.

## Projection

`floor_plot(lat; view="zx")` selects the projection plane. `view` is a
two-character string of global axes (`x`, `y`, `z`): the first maps to the
horizontal screen axis, the second to the vertical. The default `"zx"` looks at
the horizontal plane from above. The axis aspect ratio is kept 1:1 so the drawing
is not distorted.

## Using from Python

The Julia API is thin and callable from Python via
[`juliacall`](https://juliapy.github.io/PythonCall.jl/stable/):

```python
from juliacall import Main as jl
jl.seval("using PALSJulia, PALSPlot")
lat = jl.parse_and_expand_pals("machine.pals.yaml")
fp  = jl.floor_plot(lat)
jl.display(fp)
```

## Headless use

The extraction and geometry stages are independent of GLMakie and can be used
without opening a window:

```julia
tab  = element_table(lat)                 # struct-of-arrays over every element
geom = build_geometry(tab, ShapeMap())    # batched 2D draw data
```

## Status

Initial development. 2D floor plans with pan/zoom, click-to-inspect, Tao-style
shape mapping, and curved bends are working. Planned: 3D views, building-wall
overlays, and reference-orbit overlays.
