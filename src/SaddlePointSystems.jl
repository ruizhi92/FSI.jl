module SaddlePointSystems
import Base: \, A_ldiv_B!

using LinearMaps
using IterativeSolvers

export SaddleSystem


"""
    SaddleSystem((ċ, u̇, f, λ), (A⁻¹,B₁ᵀ,B₂), (M⁻¹,G₁ᵀ,G₂), (T₁ᵀ, T₂); [tol=1e-3])

Construct the computational operators for a saddle-point system of the form
\$[A 0 B₁ᵀ 0; 0 M T₁ᵀ G₁ᵀ; B₂ -T₂ 0 0; 0 G₂ 0 0][ċ; u̇; f; λ]\$.
Note that the constituent operators are passed in as a
tuple in the order seen here. Each of these operators could act on its corresponding
data type in a function-like way, e.g. `A⁻¹(u)`, or in a matrix-like way, e.g.,
`A⁻¹*u`.

The optional argument `tol` sets the tolerance for iterative solution (if
  applicable). Its default is 1e-3.

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
struct SaddleSystem{TU,TF,FA,FB2,FT2,FM,FG2,FAB,FMT,FMG,Nf,Nλ}
    A⁻¹rċ :: TU
    M⁻¹ru̇ :: TU
    tmpvec :: Vector{Float64}

    A⁻¹ :: FA
    B₂ :: FB2
    T₂ :: FB2
    M⁻¹ :: FM
    G₂ :: FG2
    A⁻¹B₁ᵀ :: FAB
    M⁻¹T₁ᵀ :: FMT
    M⁻¹G₁ᵀ :: FMG
    S  :: LinearMap

    tol :: Float64
end


function (::Type{SaddleSystem})(state::Tuple{TU,TU,TF,TF},
                                fluidop::Tuple{FH,FB1,FB2},
                                bodyop::Tuple{FM,FG1,FG2},
                                fsiop::Tuple{FT1,FT2};
                                tol::Float64=1e-3) where {TU,TF,FH,FB1,FB2,FM,FG1,FG2,FT1,FT2}
    ċ, u̇, f, λ = state

    A⁻¹, B₁ᵀ, B₂ = fluidop
    M⁻¹, G₁ᵀ, G₂ = bodyop
    T₁ᵀ, T₂ = fsiop
    sys = (A⁻¹, B₁ᵀ, B₂, M⁻¹, G₁ᵀ, G₂, T₁ᵀ, T₂)

    optypes = (TU,TF,TU,TU,TF,TU,TF,TU)
    opnames = ("A⁻¹","B₁ᵀ","B₂","M⁻¹","G₁ᵀ","G₂","T₁ᵀ","T₂")
    ops = []

    # check for methods
    for (i,typ) in enumerate(optypes)
      if method_exists(sys[i],Tuple{typ})
        push!(ops,sys[i])
    elseif method_exists(*,Tuple{typeof(sys[i]),typ})
        # generate a method that acts on TU
        push!(ops,x->sys[i]*x)
      else
        error("No valid operator for $(opnames[i]) supplied")
      end
    end

    A⁻¹, B₁ᵀ, B₂, M⁻¹, G₁ᵀ, G₂, T₁ᵀ, T₂ = ops

    ċbuffer = deepcopy(ċ)
    u̇buffer = deepcopy(u̇)
    Nf = length(f)
    Nλ = length(λ)
    tmpvec = zeros(Nf+Nλ)

    # functions in SaddleSystem attributes
    A⁻¹B₁ᵀ(f::TF) = (A⁻¹∘B₁ᵀ)(f)
    M⁻¹T₁ᵀ(f::TF) = (M⁻¹∘T₁ᵀ)(f)
    M⁻¹G₁ᵀ(λ::TF) = (M⁻¹∘G₁ᵀ)(λ)

    # functions used to create LinearMap here
    B₂A⁻¹B₁ᵀ(f::TF) = (B₂∘A⁻¹∘B₁ᵀ)(f)
    T₂M⁻¹T₁ᵀ(f::TF) = (T₂∘M⁻¹∘T₁ᵀ)(f)
    B₂A⁻¹B₁ᵀT₂M⁻¹T₁ᵀ(f::TF) = B₂A⁻¹B₁ᵀ(f) - T₂M⁻¹T₁ᵀ(f)
    T₂M⁻¹G₁ᵀ(λ::TF) = (T₂∘M⁻¹∘G₁ᵀ)(λ)
    G₂M⁻¹T₁ᵀ(f::TF) = (G₂∘M⁻¹∘T₁ᵀ)(f)
    G₂M⁻¹G₁ᵀ(λ::TF) = (G₂∘M⁻¹∘G₁ᵀ)(λ)

    function lhsmat(x::AbstractVector)
        out = zeros(Nf+Nλ)
        out[1:Nf] .= B₂A⁻¹B₁ᵀT₂M⁻¹T₁ᵀ(x[1:Nf])
        out[1:Nf] .-= T₂M⁻¹G₁ᵀ(x[Nf+1:end])
        out[Nf+1:end] .= G₂M⁻¹T₁ᵀ(x[1:Nf])
        out[Nf+1:end] .+= G₂M⁻¹G₁ᵀ(x[Nf+1:end])
        return out
    end

    S = LinearMap(lhsmat,Nf+Nλ;ismutating=false,issymmetric=false,isposdef=false)

    saddlesys = SaddleSystem{TU,TF,typeof(A⁻¹),typeof(B₂),typeof(T₂),typeof(M⁻¹),typeof(G₂),
                            typeof(A⁻¹B₁ᵀ),typeof(M⁻¹T₁ᵀ),typeof(M⁻¹G₁ᵀ),Nf,Nλ}(
                                ċbuffer,u̇buffer,tmpvec,
                                A⁻¹,B₂,T₂,M⁻¹,G₂,A⁻¹B₁ᵀ,M⁻¹T₁ᵀ,M⁻¹G₁ᵀ,S,tol)

    return saddlesys
end


function Base.show(io::IO, S::SaddleSystem{TU,TF,FA,FB2,FT2,FM,FG2,FAB,FMT,FMG,Nf,Nλ}) where {TU,TF,FA,FB2,FT2,FM,FG2,FAB,FMT,FMG,Nf,Nλ}
    println(io, "Saddle system with $Nf constraints on fluid and $Nλ constraints on body")
    println(io, "   State of type $TU")
    println(io, "   Force of type $TF")
end

# This form has error message said type mismatch
# function A_ldiv_B!(state::Tuple{TU,TU,TF,TF},
#                     sys::SaddleSystem{TU,TF,FA,FB2,FT2,FM,FG2,FAB,FMT,FMG,Nf,Nλ},
#                     rhs::Tuple{TU,TU,TF,TF}) where {TU,TF,FA,FB2,FT2,FM,FG2,FAB,FMT,FMG,Nf,Nλ}

function A_ldiv_B!(state::Tuple{TU,TU,TF,TF},
                    sys::T,
                    rhs::Tuple{TU,TU,TF,TF}) where {TU,TF,T<:SaddleSystem}
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
    tmpvec = gmres(sys.S, [rf;rλ], tol=1e-4)
    f .= tmpvec[1:Nf]
    λ .= tmpvec[Nf+1:end]

    # correction step
    ċ .= sys.A⁻¹rċ
    ċ .-= sys.A⁻¹B₁ᵀ(f)

    u̇ .= sys.M⁻¹ru̇
    u̇ .-= sys.M⁻¹T₁ᵀ(f)
    u̇ .-= sys.M⁻¹G₁ᵀ(λ)

    state = ċ, u̇, f, λ  
end


# \(sys::SaddleSystem{TU,TF,FA,FB2,FT2,FM,FG2,FAB,FMT,FMG,Nf,Nλ},rhs::Tuple{TU,TU,TF,TF}) where {TU,TF,FA,FB2,FT2,FM,FG2,FAB,FMT,FMG,Nf,Nλ} =
#     A_ldiv_B!(similar.(rhs),sys,rhs)
\(sys::T,rhs::Tuple{TU,TU,TF,TF}) where {TU,TF,T<:SaddleSystem} =
    A_ldiv_B!(similar.(rhs),sys,rhs)

end
