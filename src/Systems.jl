module Systems

using Reexport

export FluidStruct

import Whirl:@get
import Dyn3d:RKParams, TimeMarching.RK31
@reexport using Whirl:Fields, RigidBodyMotions

include("systems/fluidstruct.jl")

end
