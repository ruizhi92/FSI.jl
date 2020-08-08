# FSI.jl

*A framework for simulating rigid body systems dynamically interacting with viscous incompressible flows*

## About the package

This is a framework for simulating rigid body systems interacting with two-dimensional flow, i.e. Fluid-Structure Interaction (FSI).
Considering the components, rigid body solver take advantage of Dyn3d.jl and incompressible flow solver ViscousFlow.jl.
The FSI solver is monolithic (fully-coupled), i.e. all unknowns in the fluid-body system are solved simultaneously.

The package is currently stable for Julia 1.3 with some possible warnings, and it
- allows for both passive and active rigid body systems in an incompressible flow
- allows for both infinitely thin body (1d body) and body with finite area (2d body)
- allows for arbitrarily small fluid-body density ratios, including zero mass and neutrally buoyant cases


![](https://github.com/ruizhi92/FSI.jl/raw/master/example_gif.gif)

## Installation

This package supports *Julia* 1.3 version for now. Possible warnings may show up but
**FSI.jl** is registered in the general Julia registry. To install, type
e.g.,
```julia
] add FSI
```

Then, in any version, type
```julia
julia> using FSI
```
## References

[^1]: Ruizhi Yang and Jeff Eldredge https://escholarship.org/uc/item/0nq1t5zw
