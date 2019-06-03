using Pkg
Pkg.activate("..")

using Compat
using Compat.Test
using LinearMaps
using LinearAlgebra
using FSI

include("SaddlePointSystems.jl")
include("FSIOperators.jl")
