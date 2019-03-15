module Systems

using Reexport
using ..TimeMarching

import Whirl:@get
import FSI:RKParams,TimeMarching.RK31
@reexport using Whirl:Fields, RigidBodyMotions
@reexport using Dyn3d:ConstructSystem, RigidBodyDynamics, FluidInteraction, SpatialAlgebra, UpdateSystem

include("systems/fluidstruct.jl")

end
