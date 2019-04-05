module TimeMarching

using Reexport

import Whirl:@get
import Dyn3d:RKParams, TimeMarching.RK31, TimeMarching.Euler

@reexport using Whirl:Fields, RigidBodyMotions
@reexport using Dyn3d:ConstructSystem, RigidBodyDynamics, FluidInteraction, SpatialAlgebra, UpdateSystem

export RKParams, RK31, Euler
export IFHERK_coupled, r₁, B₂, B₁ᵀ, plan_constraints
export F, G₁ᵀ, G₂, M⁻¹, gti, UpP, UpV, T₁ᵀ, T₂, getX̃

# Functions that get extended by fluid systems
function r₁ end
function B₂ end
function B₁ᵀ end
function plan_constraints end

# function scaffold for rigid body systems
function M⁻¹ end
function F end
function G₁ᵀ end
function G₂ end
function gti end
function UpP end
function UpV end

# functions for fluid-structure coupling
function T₁ᵀ end
function T₂ end
function getX̃ end

using ..SaddlePointSystems
using ..Systems

include("timemarching/timemarching_ops.jl")
include("timemarching/ifherk_coupled.jl")

end
