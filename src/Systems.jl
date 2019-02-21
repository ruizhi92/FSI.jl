module Systems

using Reexport
using ..TimeMarching

import FSI:RKParams,TimeMarching.RK31
@reexport using Whirl:Fields, RigidBodyMotions
@reexport using Dyn3d:ConstructSystem, RigidBodyDynamics, FluidInteraction

include("systems/fluidstruct.jl")

end
