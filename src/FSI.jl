module FSI

using Reexport

# modules from Dyn3d
@reexport using Dyn3d:ConfigDataType, SpatialAlgebra, ConstructSystem, UpdateSystem
@reexport using Dyn3d:RigidBodyDynamics, FluidInteraction

# modules from Whirl
@reexport using Whirl:Utils, Fields, RigidBodyMotions

# modules from FSI
include("SaddlePointSystems.jl")
# import SaddlePointSystems
@reexport using .SaddlePointSystems

include("Systems.jl")
# import Systems
@reexport using .Systems

include("TimeMarching.jl")
# import TimeMarching
@reexport using .TimeMarching

end
