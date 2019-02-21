module TimeMarching

using Reexport

import Whirl:@get, Nodes
import Dyn3d:RKParams, TimeMarching.RK31, TimeMarching.Euler

export RKParams, RK31, Euler
export IFHERKCoupled, r₁, r₂, B₂, B₁ᵀ, plan_constraints, F, G₁ᵀ, G₂, M, gti

# Functions that get extended by fluid systems
function r₁ end
function r₂ end
function B₂ end
function B₁ᵀ end
function plan_constraints end

# function scaffold for rigid body systems
function M end
function F end
function G₁ᵀ end
function G₂ end
function gti end

using ..SaddlePointSystems

include("timemarching/ifherk_coupled.jl")

end
