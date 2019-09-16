"""
    SaddleSystem1d((ċ, f, u̇, λ), (A⁻¹,B₁ᵀ,B₂), (M⁻¹,G₁ᵀ,G₂), (T₁ᵀ, T₂); [tol=1e-3])

Construct the computational operators for a saddle-point system of the form
\$[A B₁ᵀ 0 0; B₂ 0 -T₂  0; 0 T₁ᵀ M G₁ᵀ; 0 0 G₂ 0][ċ; f; u̇; λ]\$.
Note that the constituent operators are passed in as a
tuple in the order seen here. Each of these operators could act on its corresponding
data type in a function-like way, e.g. `A⁻¹(u)`, or in a matrix-like way, e.g.,
`A⁻¹*u`.

The optional argument `tol` sets the tolerance for iterative solution (if
  applicable). Its default is 1e-4.

# Arguments

- `ċ` : change of fluid vorticity in Node form
- `f` : constraint fluid force in VectorData form
- `u̇` : change of joint velocity
- `λ` : constraint joint force
- `A⁻¹` : operator evaluating the inverse of `A` on data of type `ċ`, return type `ċ`,
            here represents the diffusion operator
- `B₁ᵀ` : operator evaluating the influence of fluid constraint force on fluid state,
            acting on `f` and returning type `ċ`
- `B₂` : operator evaluating the influence of fluid state on fluid constraints,
            acting on `ċ` and returning type `f`
- `M⁻¹` : operator evaluating the inverse of body chain inertia matrix
- `G₁ᵀ` : operator evaluating the influence of joint constraint force on body state,
            acting on `λ` and returning type `u̇`
- `G₂` : operator evaluating the influence of body state on joint constraints,
            acting on `u̇` and returning type `λ`
- `T₁ᵀ` : operator evaluating the influence of fluid constraint force on body state,
            acting on `f` and returning type `u̇`
- `T₂` : operator evaluating the influence of body state on fluid constraints,
            acting on `u̇` and returning type `f`
"""
struct SaddleSystem1d{TC,TF,TU,Tλ,FAB,FMMOD,FT1,FT2,FSF,Nf,Nλ}

    # basic operators needed
    A⁻¹B₁ᵀ :: FAB
    A⁻¹B₁ᵀf :: TC
    G₁ᵀ :: Matrix{Float64}
    G₂ :: Matrix{Float64}
    Mmod :: FMMOD
    T₁ᵀ :: FT1
    T₂ :: FT2

    # fluid and body Schur complement
    Sf :: FSF   # B₂Hᵢ₋₁,ᵢB₁ᵀ
    Sf⁻¹mat :: Matrix{Float64}
    # Sb :: LinearMap   # G₂M⁻¹G₁ᵀ

    # scratch space
    fbuffer :: TF
    tol :: Float64
end

function (::Type{SaddleSystem1d})(state::Tuple{TC,TF,TU,Tλ},
                                fluidop::Tuple{FA,FB1,FB2},
                                bodyop::Tuple{FM,FG1,FG2},
                                fsiop::Tuple{FT1,FT2};
                                tol::Float64=1e-6) where {TC,TF,TU,Tλ,FA,FB1,FB2,FM,FG1,FG2,FT1,FT2}
    ċ, f, u̇, λ = state

    A⁻¹, B₁ᵀ, B₂ = fluidop
    M, G₁ᵀ, G₂ = bodyop
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
    # bsys = (M, G₁ᵀ, G₂)
    # boptypes = (TU,Tλ,TU)
    # bopnames = ("M","G₁ᵀ","G₂")
    bsys = (M,)
    boptypes = (TU,)
    bopnames = ("M",)
    bops = []

    for (i,typ) in enumerate(boptypes)
        push!(bops,x->bsys[i]*x)
    end
    M, = bops

    Nċ = size(ċ,1)
    ċbuffer = deepcopy(ċ)
    fbuffer = deepcopy(f)
    u̇buffer = deepcopy(u̇)
    λbuffer = deepcopy(λ)
    Nu̇ = length(u̇)
    Nf = length(f)
    Nλ = length(λ)

    # fluid saddlesystem using ViscousFlow.SaddleSystem
    Sf = Sad.SaddleSystem((ċ,f),(A⁻¹,B₁ᵀ,B₂),tol=tol,
                issymmetric=false,isposdef=true,store=true,precompile=true)
    Sf⁻¹mat = inv(Sf.S⁻¹)

    # body saddlesystem using gmres
    Mmod(u̇::TU) = M(u̇) + T₁ᵀ(Sf⁻¹mat*T₂(u̇))

    # out = zeros(Nu̇+Nλ)
    # function lhsmat(x::AbstractVector{Float64})
    #     u̇buffer .= x[1:Nu̇]
    #     λbuffer .= x[Nu̇+1:end]
    #     out[1:Nf] .= Mmod(u̇buffer)
    #     out[1:Nf] .+= G₁ᵀ(λbuffer)
    #     out[Nf+1:end] .= G₂(u̇buffer)
    #     return out
    # end
    # Sb = LinearMap(lhsmat,Nu̇+Nλ;ismutating=false,issymmetric=false,isposdef=true)

    # Mmod⁻¹mat(u̇) = inv(Mmod(u̇))
    # Sb = Sad.SaddleSystem((u̇,λ),(Mmod⁻¹mat,G₁ᵀ,G₂),tol=tol,
    #             issymmetric=false,isposdef=true,store=true,precompile=true)


    # functions for correction step
    A⁻¹B₁ᵀ(f::TF) = (A⁻¹∘B₁ᵀ)(f)

    saddlesys1d = SaddleSystem1d{TC,TF,TU,Tλ,typeof(A⁻¹B₁ᵀ),typeof(Mmod),typeof(T₁ᵀ),typeof(T₂),typeof(Sf),Nf,Nλ}(
                                A⁻¹B₁ᵀ,ċbuffer,G₁ᵀ,G₂,Mmod,T₁ᵀ,T₂,Sf,Sf⁻¹mat,fbuffer,tol)

    return saddlesys1d
end


function Base.show(io::IO, S::SaddleSystem1d{TC,TF,TU,Tλ,FAB,FMMOD,FT1,FT2,FSF,Nf,Nλ}) where {TC,TF,TU,Tλ,FAB,FMMOD,FT1,FT2,FSF,Nf,Nλ}
    println(io, "Saddle system with $Nf constraints on fluid and $Nλ constraints on body")
    println(io, "   Fluid state of type $TC")
    println(io, "   Fluid force of type $TF")
    println(io, "   Body state of type $TU")
    println(io, "   Joint force of type $Tλ")
end

function ldiv!(state::Tuple{TC,TF,TU,Tλ},
                sys::TSYS,
                rhs::Tuple{TC,TF,TU,Tλ}) where {TC,TF,TU,Tλ,TSYS<:SaddleSystem1d}

    # retrive states and rhs
    rċ, rf, ru̇, rλ = rhs
    ċ, f, u̇, λ = state

    # solve for fluid with "stationary " body only
    ċ, f = sys.Sf\(rċ, rf)
# println("step1 ċ: ",ċ)
# println("step1 f: ",f)

    # solve for body with fluid added mass
    u̇ .= sys.G₂\rλ
    ru̇ .-= sys.T₁ᵀ(f)
    ru̇ .-= sys.Mmod(u̇)
    λ .= sys.G₁ᵀ\ru̇
    # tmp = gmres(sys.Sb, [ru̇;rλ], tol=sys.tol)
    # Nu̇ = length(u̇)
    # u̇ .= VectorData(tmp[1:Nu̇])
    # λ .= tmp[Nu̇+1:end]
    # u̇, λ = sys.Sb\(ru̇, rλ)
    sys.fbuffer .= sys.T₂(u̇)
# println("Sf⁻¹mat: ",sys.Sf⁻¹mat)
    sys.fbuffer .= sys.Sf⁻¹mat*sys.fbuffer

# println("step2 u̇: ",u̇)
# println("step2 λ: ",λ)

    # correct fluid state and force with moving body effect
    sys.A⁻¹B₁ᵀf .= sys.A⁻¹B₁ᵀ(sys.fbuffer)
    ċ .-= sys.A⁻¹B₁ᵀf
    f .+= sys.fbuffer

# println("step3 ċ: ",ċ)
# println("step3 f: ",f)

    state = ċ, f, u̇, λ
end


\(sys::TSYS,rhs::Tuple{TC,TF,TU,Tλ}) where {TC,TF,TU,Tλ,TSYS<:SaddleSystem1d} =
    ldiv!(deepcopy.(rhs),sys,rhs)
