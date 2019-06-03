module SaddlePointSystems

using LinearMaps
using LinearAlgebra

import Base: \,+,*
import LinearAlgebra: ldiv!

using ..Fields:VectorData, ScalarData
export SaddleSystem

# extend method from FSI
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
"""
    SaddleSystem((ċ, u̇, f, λ), (A⁻¹,B₁ᵀ,B₂), (M⁻¹,G₁ᵀ,G₂), (T₁ᵀ, T₂); [tol=1e-4])

Construct the computational operators for a saddle-point system of the form
\$[A 0 B₁ᵀ 0; 0 M T₁ᵀ G₁ᵀ; B₂ -T₂ 0 0; 0 G₂ 0 0][ċ; u̇; f; λ]\$.
Note that the constituent operators are passed in as a
tuple in the order seen here. Each of these operators could act on its corresponding
data type in a function-like way, e.g. `A⁻¹(u)`, or in a matrix-like way, e.g.,
`A⁻¹*u`.

The optional argument `tol` sets the tolerance for iterative solution (if
  applicable). Its default is 1e-4.

# Arguments

- `u` : example of state vector data.
- `f` : example of constraint force vector data. This data must be of
        AbstractVector supertype.
- `A⁻¹` : operator evaluating the inverse of `A` on data of type `u`, return type `u`
- `B₁ᵀ` : operator evaluating the influence of constraint force,
            acting on `f` and returning type `u`
- `B₂` : operator evaluating the influence of state vector on constraints,
            acting on `u` and returning type `f`
"""
struct SaddleSystem{TC,TU,TF,Tλ,FA,FB2,FT2,FM,FG2,FAB,FMT,FMG,Nf,Nλ}
    A⁻¹rċ :: TC
    M⁻¹ru̇ :: TU
    tmpvec :: Vector{Float64}

    A⁻¹ :: FA
    B₂ :: FB2
    T₂ :: FT2
    M⁻¹ :: FM
    G₂ :: FG2
    A⁻¹B₁ᵀ :: FAB
    M⁻¹T₁ᵀ :: FMT
    M⁻¹G₁ᵀ :: FMG
    S  :: LinearMap
    S⁻¹ :: Factorization{Float64}

    tol :: Float64
end


function (::Type{SaddleSystem})(state::Tuple{TC,TU,TF,Tλ},
                                fluidop::Tuple{FA,FB1,FB2},
                                bodyop::Tuple{FM,FG1,FG2},
                                fsiop::Tuple{FT1,FT2};
                                tol::Float64=1e-4) where {TC,TU,TF,Tλ,FA,FB1,FB2,FM,FG1,FG2,FT1,FT2}
    ċ, u̇, f, λ = state

    A⁻¹, B₁ᵀ, B₂ = fluidop
    M⁻¹, G₁ᵀ, G₂ = bodyop
    T₁ᵀ, T₂ = fsiop

    # check for fluid methods
    fsys = (A⁻¹, B₁ᵀ, B₂)
    foptypes = (TC,TF,TC)
    fopnames = ("A⁻¹","B₁ᵀ","B₂")
    fops = []

    for (i,typ) in enumerate(foptypes)
      if hasmethod(fsys[i],Tuple{typ})
        push!(fops,fsys[i])
    elseif hasmethod(*,Tuple{typeof(fsys[i]),typ})
        push!(fops,x->fsys[i]*x)
      else
        error("No valid operator for $(fopnames[i]) supplied")
      end
    end
    A⁻¹, B₁ᵀ, B₂ = fops

    # check for body methods
    bsys = (M⁻¹, G₁ᵀ, G₂)
    boptypes = (TU,Tλ,TU)
    bopnames = ("M⁻¹","G₁ᵀ","G₂")
    bops = []

    for (i,typ) in enumerate(boptypes)
        push!(bops,x->bsys[i]*x)
    end
    M⁻¹, G₁ᵀ, G₂ = bops

    fbuffer = deepcopy(f)
    λbuffer = deepcopy(λ)
    Nf = length(f)
    Nλ = length(λ)
    N = Nf + Nλ
    tmpvec = zeros(N)
    fλbuffer = zeros(N)

    # functions in SaddleSystem attributes
    A⁻¹B₁ᵀ(f::TF) = (A⁻¹∘B₁ᵀ)(f)
    M⁻¹T₁ᵀ(f::TF) = (M⁻¹∘T₁ᵀ)(f)
    M⁻¹G₁ᵀ(λ::Tλ) = (M⁻¹∘G₁ᵀ)(λ)

    # functions used to create LinearMap lhsmat
    B₂A⁻¹B₁ᵀ(f::TF) = (B₂∘A⁻¹∘B₁ᵀ)(f)
    T₂M⁻¹T₁ᵀ(f::TF) = (T₂∘M⁻¹∘T₁ᵀ)(f)
    B₂A⁻¹B₁ᵀT₂M⁻¹T₁ᵀ(f::TF) = B₂A⁻¹B₁ᵀ(f) - T₂M⁻¹T₁ᵀ(f)
    T₂M⁻¹G₁ᵀ(λ::Tλ) = (T₂∘M⁻¹∘G₁ᵀ)(λ)
    G₂M⁻¹T₁ᵀ(f::TF) = (G₂∘M⁻¹∘T₁ᵀ)(f)
    G₂M⁻¹G₁ᵀ(λ::Tλ) = (G₂∘M⁻¹∘G₁ᵀ)(λ)

    # create the left hand side matrix
    out = zeros(N)
    function lhsmat(x::AbstractVector{Float64})
        fbuffer .= x[1:Nf]
        λbuffer .= x[Nf+1:end]

        out[1:Nf] .= B₂A⁻¹B₁ᵀT₂M⁻¹T₁ᵀ(fbuffer)
        out[1:Nf] .-= T₂M⁻¹G₁ᵀ(λbuffer)
        out[Nf+1:end] .= G₂M⁻¹T₁ᵀ(fbuffer)
        out[Nf+1:end] .+= G₂M⁻¹G₁ᵀ(λbuffer)
        return out
    end

    S = LinearMap(lhsmat,N;ismutating=false,issymmetric=false,isposdef=false)

    # construct S into real matrix
    Smat = zeros(N,N)
    for i = 1:N
      tmpvec[i] = 1.0
      fλbuffer .= S*tmpvec
      Smat[1:N,i] .= fλbuffer
      tmpvec[i] = 0.0
    end
    S⁻¹ = factorize(Smat)
    tmpvec .= 0.0

    saddlesys = SaddleSystem{TC,TU,TF,Tλ,typeof(A⁻¹),typeof(B₂),typeof(T₂),typeof(M⁻¹),typeof(G₂),
                            typeof(A⁻¹B₁ᵀ),typeof(M⁻¹T₁ᵀ),typeof(M⁻¹G₁ᵀ),Nf,Nλ}(
                                ċ,u̇,tmpvec,
                                A⁻¹,B₂,T₂,M⁻¹,G₂,A⁻¹B₁ᵀ,M⁻¹T₁ᵀ,M⁻¹G₁ᵀ,S,S⁻¹,tol)

    return saddlesys
end


function Base.show(io::IO, S::SaddleSystem{TC,TU,TF,Tλ,FA,FB2,FT2,FM,FG2,FAB,FMT,FMG,Nf,Nλ}) where {TC,TU,TF,Tλ,FA,FB2,FT2,FM,FG2,FAB,FMT,FMG,Nf,Nλ}
    println(io, "Saddle system with $Nf constraints on fluid and $Nλ constraints on body")
    println(io, "   Fluid state of type $TC")
    println(io, "   Fluid force of type $TF")
    println(io, "   Body state of type $TU")
    println(io, "   Joint force of type $Tλ")
end

function ldiv!(state::Tuple{TC,TU,TF,Tλ},
                    sys::T,
                    rhs::Tuple{TC,TU,TF,Tλ}) where {TC,TU,TF,Tλ,T<:SaddleSystem}
    rċ, ru̇, rf, rλ = rhs
    ċ, u̇, f, λ = state
    Nf = length(f)

    # create right hand side for the apply constraints step
    sys.A⁻¹rċ .= sys.A⁻¹(rċ)
    sys.M⁻¹ru̇ .= sys.M⁻¹(ru̇)
    rf .= -rf
    rf .+= sys.B₂(sys.A⁻¹rċ)
    rf .-= sys.T₂(sys.M⁻¹ru̇)
    rλ .= -rλ
    rλ .+= sys.G₂(sys.M⁻¹ru̇)

    # solve for forcing terms
    # tmpvec = gmres(sys.S, [rf;rλ], tol=1e-4)
    tmpvec = ldiv!(sys.S⁻¹, [rf;rλ])
    f .= VectorData(tmpvec[1:Nf])
    λ .= tmpvec[Nf+1:end]

    # correction step
    ċ .= sys.A⁻¹rċ
    ċ .-= sys.A⁻¹B₁ᵀ(f)

    u̇ .= sys.M⁻¹ru̇
    u̇ .-= sys.M⁻¹T₁ᵀ(f)
    u̇ .-= sys.M⁻¹G₁ᵀ(λ)

    state = ċ, u̇, f, λ
end


\(sys::T,rhs::Tuple{TC,TU,TF,Tλ}) where {TC,TU,TF,Tλ,T<:SaddleSystem} =
    ldiv!(deepcopy.(rhs),sys,rhs)

end
