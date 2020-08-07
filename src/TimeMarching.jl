module TimeMarching

using Reexport
using LinearAlgebra

import ViscousFlow:@get
import Dyn3d.RKParams
import Dyn3d.TimeMarching.RK31, Dyn3d.TimeMarching.Euler, Dyn3d.TimeMarching.RK4

@reexport using ViscousFlow.Fields
@reexport using ViscousFlow.RigidBodyMotions
@reexport using Dyn3d.ConstructSystem
@reexport using Dyn3d.RigidBodyDynamics
@reexport using Dyn3d.FluidInteraction
@reexport using Dyn3d.SpatialAlgebra
@reexport using Dyn3d.UpdateSystem


export RKParams, RK31, Euler, RK32, RK4
export IFHERK_coupled, r₁, U_inf, B₂, B₁ᵀ
export F, G₁ᵀ, G₂, M, gti, UpP, UpV, T₁ᵀ, T₂, getX̃, BodyGridToVectorData

# Functions that get extended by fluid systems
function r₁ end
function U_inf end
function B₂ end
function B₁ᵀ end

# function scaffold for rigid body systems
function M end
function Mf end
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
