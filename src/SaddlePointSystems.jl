module SaddlePointSystems

using LinearMaps
using LinearAlgebra

import Base: \,+,*
import LinearAlgebra: ldiv!, inv

using ..Fields:VectorData, ScalarData
export SaddleSystem1d,SaddleSystem2d

# import ViscousFlow.SaddlePointSystems with a different name
using ImportMacros
@import ViscousFlow.SaddlePointSystems as Sad

# extend method from FSI
function VectorData(a::Array{T,2}) where T<: Real
    return VectorData(a[:,1],a[:,2])
end
function VectorData(a::Vector{T}) where T<: Real
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

function(*)(S::Matrix{T},a::VectorData) where T<: Real
    return VectorData(S*[a.u;a.v])
end

include("saddlepointsystem/saddlepointsystem_1dbody.jl")
include("saddlepointsystem/saddlepointsystem_2dbody.jl")

end
