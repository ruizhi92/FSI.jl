module Systems

using Reexport

export NavierStokes

import ViscousFlow:@get

using ..Fields
using ..TimeMarching
using ..RigidBodyMotions

include("systems/navier_stokes.jl")

end
