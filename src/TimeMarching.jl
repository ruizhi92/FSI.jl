module TimeMarching

using Reexport
using LinearAlgebra

import ViscousFlow:@get
import Base: +,*

using ..RigidBodyMotions
using ..SaddlePointSystems
using ..Fields:VectorData, ScalarData

export RKParams, RK31, Euler
export IFHERK_sc2d, r₁, r₂, B₂, B₁ᵀ, plan_constraints


# Functions that get extended by individual systems
function r₁ end
function r₂ end
function B₂ end
function B₁ᵀ end
function plan_constraints end

struct RKParams{N}
  c::Vector{Float64}
  a::Matrix{Float64}
end

const RK31 = RKParams{3}([0.5, 1.0, 1.0],
                      [1/2        0        0
                       √3/3 (3-√3)/3        0
                       (3+√3)/6    -√3/3 (3+√3)/6])
const Euler = RKParams{1}([1.0],ones(1,1))


# extend method from ViscousFlow
function VectorData(a::Array{Float64,2})
    return VectorData(a[:,1],a[:,2])
end
function VectorData(a::Vector{Float64})
    N = length(a)÷2
    return VectorData(a[1:N],a[N+1:end])
end
function (+)(a::VectorData, b::VectorData)
    c = VectorData(a)
    c.u .= a.u .+ b.u
    c.v .= a.v .+ b.v
    return c
end
function (+)(a::ScalarData, b::ScalarData)
    c = ScalarData(a)
    c.data .= a.data .+ b.data
    return c
end
function (*)(ω::T, a::VectorData) where T<: Real
    c = VectorData(a)
    c.u .= ω.*a.u
    c.v .= ω.*a.v
    return c
end
function (*)(a::VectorData,ω::T) where T<: Real
    c = VectorData(a)
    c.u .= ω.*a.u
    c.v .= ω.*a.v
    return c
end
function (*)(no::Nothing,a::VectorData)
    return VectorData(a)
end


include("timemarching/ifherk_sc2d.jl")

end
