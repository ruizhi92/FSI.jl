# Problem set up
In order to create a fluid-structure interaction problem, we need to specify three things:
- rigid body systems properties
- flow domain set up
- create fluid-body interface

Rigid body system properties are specified using functions imported from Dyn3d.jl. User should
refer to documentation of *Dyn3d.jl*. For fluid domain set up, user should refer to *ViscousFlow.jl*
or *CartesianGrids.jl*. Multiple FSI examples are provided in examples folder.


## Methods
```@autodocs
Modules = [TimeMarching, SaddlePointSystems, Systems, Tools]
Order   = [:type, :function]
```

## Index
```@index
Pages = ["problemsetup.md"]
```
