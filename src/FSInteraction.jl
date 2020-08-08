module FSInteraction

using Reexport

# modules from Dyn3d
@reexport using Dyn3d.ConfigDataType
@reexport using Dyn3d.SpatialAlgebra
@reexport using Dyn3d.ConstructSystem
@reexport using Dyn3d.UpdateSystem
@reexport using Dyn3d.RigidBodyDynamics
@reexport using Dyn3d.FluidInteraction
@reexport using Dyn3d.ConfigDataType

# modules from ViscousFlow
@reexport using ViscousFlow.Utils
@reexport using ViscousFlow.Fields
@reexport using ViscousFlow.RigidBodyMotions

# modules from FSInteraction
include("SaddlePointSystems.jl")
@reexport using .SaddlePointSystems

include("Systems.jl")
@reexport using .Systems

include("TimeMarching.jl")
@reexport using .TimeMarching

include("Tools.jl")
@reexport using .Tools

end
