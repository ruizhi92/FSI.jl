# FSI.jl

_A framework for simulating rigid body systems dynamically interacting with viscous incompressible flows_


[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://ruizhi92.github.io/FSI.jl/latest)
[![Build Status](https://travis-ci.org/ruizhi92/FSI.jl.png?branch=master)](https://travis-ci.org/ruizhi92/FSI.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/o7du221qb5fqa5s1/branch/master?svg=true)](https://ci.appveyor.com/project/ruizhi92/fsi-jl/branch/master)

## About the package

This is a framework for simulating rigid body systems interacting with two-dimensional flow, i.e. Fluid-Structure Interaction (FSI).
Considering the components, rigid body solver take advantage of Dyn3d.jl and incompressible flow solver ViscousFlow.jl.
The monolithic and versatile algorithm for this solver refers to https://escholarship.org/uc/item/0nq1t5zw

The package is currently stable for Julia 1.3 with some possible warnings, and it
- allows for both passive and active rigid body systems in an incompressible flow
- allows for both infinitely thin body (1d body) and body with finite area (2d body)
- allows for arbitrarily small fluid-body density ratios, including zero mass and neutrally buoyant cases


**FSI.jl** is registered in the general Julia registry. To install, type
e.g.,
```julia
] add FSI
```

Then, in any version, type
```julia
julia> using FSI
```
See the example Jupyter notebooks in the examples folder.

![](https://github.com/ruizhi92/FSI.jl/raw/master/example_gif.gif)
