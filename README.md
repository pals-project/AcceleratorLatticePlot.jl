# PALSPlot.jl

[![Julia Tests](https://github.com/pals-project/PALSPlot.jl/actions/workflows/test.yaml/badge.svg)](https://github.com/pals-project/PALSPlot.jl/actions/workflows/test.yaml)

Floor-plan plotting for [PALS](https://github.com/campa-consortium/pals) lattices,
built on [PALSJulia](https://github.com/pals-project/PALSJulia.jl) and
[Makie](https://docs.makie.org/stable/).

Given the expanded lattice from `PALSJulia.parse_and_expand_pals`, PALSPlot draws
the machine projected onto a plane. Each element is rendered as a shape — sized,
colored and labeled by a [Tao](https://www.classe.cornell.edu/bmad/)-style rule
table — placed and oriented from the floor coordinates the expander computes
(the `full_expanded` view; see [Which expanded view](#which-expanded-view)).
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

### Choosing a backend

PALSPlot depends on **Makie**, not on any one Makie backend: every drawing call
it makes is Makie core, and a backend is needed only to put the figure
somewhere. So you add the one you want.

```console
julia --project=. -e 'using Pkg; Pkg.add("GLMakie")'      # interactive window
julia --project=. -e 'using Pkg; Pkg.add("CairoMakie")'   # write PDF/PNG/SVG
```

Loading a backend activates it, so `using GLMakie` is all it takes. The
extraction and geometry stages, and building a `Figure`, need no backend at all.

## Quick start

```julia
using PALSJulia
using PALSPlot
using GLMakie

lat = parse_and_expand_pals("machine.pals.yaml")
fp  = floor_plot(lat)
display(fp)              # opens an interactive window
```

or, with no display available, render straight to a file:

```julia
using CairoMakie
save("floor.pdf", floor_plot(lat).figure)
```

or run the bundled example:

```console
julia --project=. examples/floor_plan.jl path/to/machine.pals.yaml
```

**Controls**

| gesture | action |
|---|---|
| scroll wheel | zoom in/out at the cursor |
| **right-drag** | pan |
| left-drag | rubber-band zoom to a rectangle |
| `ctrl` + **left**-click | reset to the full view |
| `ctrl` + `shift` + **left**-click | reset to limits recomputed from the data |
| left-click (no modifier) | select an element (lists its parameters in the side panel) |

Everything but selection is a Makie default `Axis` binding. Makie has no
double-click binding, so **double-clicking does nothing** — `ctrl`-left-click is
the reset.

### Which expanded view

`parse_and_expand_pals` returns the expanded lattice twice. PALSPlot reads
**`lat.full_expanded`**, the view in which the expander has computed every
dependent parameter: each element's `FloorP` (position and orientation at its
upstream end), its `s_position`, and — for a bend — the `BendP.angle_ref` the
arc is drawn from. `lat.expanded` holds the same lattice pruned back to what the
author wrote, so none of that is in it and nothing could be placed.

`full_expanded` also caps each branch with a zero-length `branch_end`
`Placeholder` holding the downstream end of the last element. It appears in the
element table like any other element; the default shape table draws it as bare
centerline.

### Build pals-cpp against libc++

PALSPlot needs the pals-cpp C library to be built with the **same C++ runtime as
Julia and Makie — LLVM libc++ (Apple clang)**. If pals-cpp is instead built with
GCC/libstdc++ (e.g. because `CC`/`CXX` point at MacPorts GCC), then loading a
Makie backend and then parsing a lattice file **aborts the process** (signal 6):
the two C++ exception runtimes clash in the library's file reader.

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
    ele_shape("Bend",           :box,    :blue;  size=0.5),
    ele_shape("Wiggler",        :box,    :green; size=0.5),
])
floor_plot(lat; shapes=shapes)
```

The class is a PALS element kind (`Bend`, `Quadrupole`, `RFCavity`, `Fork`, …
from the standard's `lattice-element-kinds.md`), matched case-insensitively;
`"*"` matches any kind. The built-in table covers every kind a beam line can
hold, and anything else falls through to a small unlabeled box.

Available shapes: `:box`, `:xbox`, `:x`, `:bow_tie`, `:rbow_tie`, `:diamond`,
`:u_triangle`, `:d_triangle`, `:l_triangle`, `:r_triangle`, `:circle`, `:none`.
`size` is the transverse half-height in meters; `label` is `:name`, `:s`, or
`:none`. Pass `defaults=false` to `ShapeMap` to use only your own rules.

Labels are drawn perpendicular to the centerline, running outward from the
element, and tilted only within ±90° so none of them come out upside down.
Elements are strung out *along* the line, so laying their labels across it is
what keeps neighbouring names from running into each other.

Elements that sit on top of each other — a pickup with its two correctors, say —
defeat that, since their labels share one anchor and one ray. Those are detected
and stacked along the ray instead, the element earliest in the branch keeping the
place nearest the centerline. `floor_plot(...; label_sep=…)` sets how close two
anchors have to be to count as colliding, in units of the elements' half-height;
`label_sep=0` turns the stacking off.

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
jl.seval("using PALSJulia, PALSPlot, GLMakie")
lat = jl.parse_and_expand_pals("machine.pals.yaml")
fp  = jl.floor_plot(lat)
jl.display(fp)
```

## Headless use

No backend is needed to extract a lattice, compute its geometry, or build the
figure — only to put that figure somewhere:

```julia
tab  = element_table(lat)                 # struct-of-arrays over every element
geom = build_geometry(tab, ShapeMap())    # batched 2D draw data
fig  = floor_plot(tab).figure             # a Makie Figure, undisplayed
```

With CairoMakie loaded that figure goes straight to a file, no display or OpenGL
involved anywhere:

```julia
using CairoMakie
save("floor.pdf", fig)
```

## Status

Initial development. 2D floor plans with pan/zoom, click-to-inspect, Tao-style
shape mapping, and curved bends are working. Planned: 3D views, building-wall
overlays, and reference-orbit overlays.

Known gap: **a reference curve that leaves the projection plane is drawn as
though it did not.** Every element is placed exactly, from its own `FloorP`, but
the frame the drawing is carried across it on is built from the heading `theta`
alone — `FloorP.phi` and `FloorP.psi` are ignored, as is `BendP.tilt_ref`. So on
a lattice that stays in the plane the drawing closes on the expander's floor
coordinates to rounding (checked in the tests against `bta.pals.yaml`), while on
one that does not:

* a straight element is drawn at its full length rather than its projected
  length, overshooting by `L·(1 − cos φ)` — 0.25 m on a 1.3 m cavity at
  `phi = −0.63` in `convert.pals.yaml`;
* a bend with a non-zero `tilt_ref` gets the right endpoints and the wrong arc
  between them.

Both fall out of building the frame from the standard's W matrix
(`pals_floor.h`, Eq. www) rather than from `theta`.

## Testing

```console
julia --project=. -e 'using Pkg; Pkg.develop(path="../PALSJulia"); Pkg.instantiate(); Pkg.test()'
```

CI runs the suite on Julia 1.11 and 1.12, on Linux and macOS, against **`main`
of PALSJulia and pals-cpp** — both are checked out fresh and pals-cpp is built
from source on every run. PALSJulia binds the pals-cpp C structs by layout and
PALSPlot draws from the floor coordinates pals-cpp computes, so a change to
either that PALSPlot has not followed shows up as a failure here rather than as
an empty window later.

The suite renders through **CairoMakie**, which is pure software and so needs no
display, no OpenGL and no X server. That is what lets it run everywhere: a
headless macOS runner cannot give GLMakie an OpenGL context at all — GLMakie will
not even precompile there, since its precompile workload opens a GLFW window.

GLMakie is what you actually open a window with, so it is covered separately by
`test/glmakie_smoke.jl`, which CI runs on Linux under `xvfb-run`. To run it
yourself, in an environment with GLMakie added:

```console
julia --project=. test/glmakie_smoke.jl
```
