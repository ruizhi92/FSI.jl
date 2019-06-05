module FSI

using Reexport

# modules from ViscousFlow
@reexport using ViscousFlow.Utils
@reexport using ViscousFlow.Fields
@reexport using ViscousFlow.RigidBodyMotions
@reexport using ViscousFlow.SaddlePointSystems


# use modules defined in FSI
include("TimeMarching.jl")
@reexport using .TimeMarching

include("Systems.jl")
@reexport using .Systems

end
