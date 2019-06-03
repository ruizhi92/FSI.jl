module Systems

using Reexport

export FluidStruct

import ViscousFlow:@get
import Dyn3d:RKParams, TimeMarching.RK31
@reexport using ViscousFlow.Fields
@reexport using ViscousFlow.RigidBodyMotions

include("systems/fluidstruct.jl")

end
