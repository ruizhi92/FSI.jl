"""
    IFHERK(rk::u,f,Δt,plan_intfact,B₁ᵀ,B₂,r₁,r₂;[tol=1e-3])

Construct an integrator to advance a system of the form

du/dt - Au = -B₁ᵀf + r₁(u,t)
B₂u = r₂(u,t)

The resulting integrator will advance the system `(u,f)` by one time step, `Δt`.
The optional argument `tol` sets the tolerance of iterative saddle-point solution,
if applicable.

# Arguments

- `u` : example of state vector data
- `f` : example of constraint force vector data
- `Δt` : time-step size
- `plan_intfact` : constructor to set up integrating factor operator for `A` that
              will act on type `u` (by left multiplication) and return same type as `u`
- `plan_constraints` : constructor to set up the
- `B₁ᵀ` : operator acting on type `f` and returning type `u`
- `B₂` : operator acting on type `u` and returning type `f`
- `r₁` : operator acting on type `u` and `t` and returning `u`
- `r₂` : operator acting on type `u` and `t` and returning type `f`
"""
struct IFHERK_coupled{FH,FB1,FB2,FM,FG1,FG2,FUPP,FUPV,FT1,FT2,FX,FR11,FR12,FR22,FS,TW,TF}

    # time step size
    Δt :: Float64

    # rk parameters
    rk :: RKParams
    rkdt :: RKParams

    # Fluid operators
    H :: Vector{FH}
    B₁ᵀ :: FB1
    B₂ :: FB2

    # Body operators
    M⁻¹ :: FM
    G₁ᵀ :: FG1
    G₂ :: FG2
    UpP :: FUPP
    UpV :: FUPV

    # FSI operators
    T₁ᵀ :: FT1
    T₂ :: FT2
    getX̃ :: FX

    # right hand side operators
    r₁₁ :: FR11
    r₁₂ :: FR12
    r₂₂ :: FR22

    # saddle-point system
    S :: Vector{FS}

    # body grid information
    bgs :: Vector{BodyGrid}

    # scratch space
    w₀ :: TW
    qJ₀ :: Vector{Float64}
    v₀ :: Vector{Float64}
    wbuffer :: TW
    vbuffer :: Vector{Float64}
    fbuffer :: TF
    λbuffer :: Vector{Float64}
    ċ :: Vector{TW}
    vJ :: Vector{Vector{Float64}}
    v̇ :: Vector{Vector{Float64}}

    # iterative solution tolerance
    tol :: Float64

end

function (::Type{IFHERK_coupled})(Δt::Float64, bd::BodyDyn, bgs::Vector{BodyGrid},
                state::Tuple{TW,Vector{Float64},Vector{Float64},TF,Vector{Float64}},
                fluidop::Tuple{FI,FB1,FB2},
                bodyop::Tuple{FM,FG1,FG2,FUPP,FUPV},
                fsiop::Tuple{FT1,FT2,FX},
                rhs::Tuple{FR11,FR12,FR22},
                fsys::FluidStruct{NX,NY,N};
                tol::Float64=1e-3, rk::RKParams=RK31) where {TW,TF,FI,FB1,FB2,FM,FG1,FG2,FUPP,FUPV,FT1,FT2,FX,FR11,FR12,FR22,NX,NY,N}

    w, qJ, v, f, λ = state
    plan_intfact, B₁ᵀ, B₂ = fluidop
    M⁻¹, G₁ᵀ, G₂, UpP, UpV = bodyop
    T₁ᵀ, T₂, getX̃ = fsiop
    r₁₁,r₁₂,r₂₂ = rhs

    # scratch space
    w₀ = deepcopy(w)
    qJ₀ = deepcopy(qJ)
    v₀ = deepcopy(v)
    wbuffer = deepcopy(w)
    vbuffer = deepcopy(v)
    fbuffer = deepcopy(f)
    λbuffer = deepcopy(λ)
    ċ = [deepcopy(w) for i = 1:rk.st]
    v̇ = [deepcopy(v) for i = 1:rk.st]
    vJ = [deepcopy(qJ) for i = 1:rk.st+1]

    # construct an array of operators for the integrating factor
    # H[i-1] corresponds to H((cᵢ - cᵢ₋₁)Δt)
    dclist = diff(rk.c)
    Hlist = [plan_intfact(dc*Δt,w) for dc in unique(dclist)]
    H = [Hlist[i] for i in indexin(dclist,unique(dclist))]

    # fuse the time step size into the coefficients for some cost savings
    rkdt = deepcopy(rk)
    rkdt.a .*= Δt
    rkdt.c .*= Δt

    # preform the saddle-point systems, they are overwritten with time
    Slist = [construct_saddlesys(0.0,1,rkdt.a,bd,bgs,qJ,vJ,(w,v,f,λ),(H[i],B₁ᵀ,B₂),
        (M⁻¹,G₁ᵀ,G₂,UpP,r₁₂,r₂₂),(T₁ᵀ,T₂,getX̃),fsys)[1][1] for i=1:rk.st]
    S = [Slist[i] for i in indexin(dclist,unique(dclist))]

    htype,_ = typeof(H).parameters
    stype,_ = typeof(S).parameters

    # actually construct ifherk_coupled object
    ifherksys = IFHERK_coupled{htype,typeof(B₁ᵀ),typeof(B₂),typeof(M⁻¹),typeof(G₁ᵀ),typeof(G₂),
                    typeof(UpP),typeof(UpV),typeof(T₁ᵀ),typeof(T₂),typeof(getX̃),
                    typeof(r₁₁),typeof(r₁₂),typeof(r₂₂),stype,TW,TF}(Δt,rk,rkdt,
                                H,B₁ᵀ,B₂,
                                M⁻¹,G₁ᵀ,G₂,UpP,UpV,
                                T₁ᵀ,T₂,getX̃,
                                r₁₁,r₁₂,r₂₂,
                                S,bgs,
                                w₀,qJ₀,v₀,wbuffer,vbuffer,fbuffer,λbuffer,
                                ċ,vJ,v̇,tol)

    return ifherksys
end

function Base.show(io::IO, scheme::IFHERK_coupled)
    println(io, "Stage-$(scheme.rk.st)+ IF-HERK integrator with")
    println(io, "   Time step size $(scheme.Δt)")
end

"""
    construct_saddlesys()
takes in qJ,vJ to construct the lhs saddle point system, also return operators used
in ifherk_coupled.
"""
function construct_saddlesys(t::Float64, stage::Int64, rkdt_a::Matrix{Float64},
                             bd::BodyDyn, bgs::Vector{BodyGrid},
                             qJ::Vector{Float64}, vJ::Vector{Vector{Float64}},
                             state::Tuple{TW,Vector{Float64},TF,Vector{Float64}},
                             fluidop::Tuple{FH,FB1,FB2},
                             bodyop::Tuple{FM,FG1,FG2,FUPP,FF,FGTI},
                             fsiop::Tuple{FT1,FT2,FX},fsys::FluidStruct{NX,NY,N};
                             tol::Float64=1e-3) where {TW,TF,FH,FB1,FB2,FM,FG1,FG2,FUPP,FF,FGTI,FT1,FT2,FX,NX,NY,N}

    w, v, f, λ = state
    H, B₁ᵀ, B₂ = fluidop
    M⁻¹, G₁ᵀ, G₂, UpP, F, gti = bodyop
    T₁ᵀ, T₂, getX̃ = fsiop

    #----------------------- Operators at time tᵢ₋₁ ----------------------------
    # Body operators
    Mᵢ₋₁ = M⁻¹(bd)
    Fᵢ₋₁ = F(bd)
    G₁ᵀᵢ₋₁ = G₁ᵀ(bd)

    # get current body points coordinates
    fsys.X̃ = getX̃(bd, bgs)
    regop = Regularize(fsys.X̃,fsys.Δx;issymmetric=true)
    fsys.Hmat, _ = RegularizationMatrix(regop,VectorData{N}(),Edges{Primal,NX,NY}())

    # Fluid operators
    B₁ᵀᵢ₋₁ = f -> B₁ᵀ(f, fsys)

    # FSI operators
    T₁ᵀᵢ₋₁ = f -> T₁ᵀ(bd, bgs, f)

    # Integrate joints qJ and update bs and js in timemarching
    for k = 1:stage-1
        qJ += rkdt_a[stage,k]*vJ[k]
    end
    bd = UpP(bd, qJ)

    #---------- Operators at time tᵢ with updated body coordinates -------------
    # Body operators
    G₂ᵢ = G₂(bd)
    gtiᵢ = gti(bd, t)

    # get current body points coordinates
    fsys.X̃ = getX̃(bd, bgs)
    regop = Regularize(fsys.X̃,fsys.Δx;issymmetric=true)
    _, fsys.Emat = RegularizationMatrix(regop,VectorData{N}(),Edges{Primal,NX,NY}())

    # Fluid operators
    B₂ᵢ = w -> B₂(w, fsys)        # Fluid operators

    T₂ᵢ = v -> T₂(bd, bgs, v)     # FSI operators

    # Actually call SaddleSystem
    S = SaddleSystem((w, v, f, λ),
                     (H, B₁ᵀᵢ₋₁, B₂ᵢ),
                     (Mᵢ₋₁, G₁ᵀᵢ₋₁, G₂ᵢ),
                     (T₁ᵀᵢ₋₁, T₂ᵢ))

    return (S, B₂ᵢ, G₂ᵢ, T₂ᵢ, Fᵢ₋₁, gtiᵢ), qJ, bd

end


function (scheme::IFHERK_coupled{FH,FB1,FB2,FM,FG1,FG2,FUPP,FUPV,FT1,FT2,FX,FR11,FR12,FR22,FS,TW,TF})(t::Float64,
    w::TW,qJ::Vector{Float64},v::Vector{Float64},bd::BodyDyn,fsys::FluidStruct{NX,NY,N}) where {FH,FB1,
    FB2,FM,FG1,FG2,FUPP,FUPV,FT1,FT2,FX,FR11,FR12,FR22,FS,TW,TF,NX,NY,N}

    @get scheme (rk,rkdt,H,B₁ᵀ,B₂,M⁻¹,G₁ᵀ,G₂,UpP,UpV,T₁ᵀ,T₂,getX̃,r₁₁,r₁₂,r₂₂,S,bgs,tol)
    @get scheme (w₀,qJ₀,v₀,wbuffer,vbuffer,fbuffer,λbuffer,ċ,vJ,v̇)
    @get bd (bs,js,sys)

    f = deepcopy(fbuffer)
    λ = deepcopy(λbuffer)

    #---------------------------- first stage, i = 1 ---------------------------
    i = 1
    tᵢ = t
    w₀ .= w
    qJ₀ .= qJ
    v₀ .= v
    # update vJ using v
    bd, vJ[1] = UpV(bd, v₀)

    #----------------------------- stage 2 to st+1 -----------------------------
    for i = 2:rk.st+1
        # set time
        tᵢ = t + rkdt.c[i]

        # construct saddlesys
        (S, B₂ᵢ, G₂ᵢ, T₂ᵢ, Fᵢ₋₁, gtiᵢ), qJ, bd = construct_saddlesys(tᵢ,i,rkdt.a,
            bd,bgs,qJ₀,vJ,(w,v,f,λ),(H[i-1],B₁ᵀ,B₂),(M⁻¹,G₁ᵀ,G₂,UpP,r₁₂,r₂₂),(T₁ᵀ,T₂,getX̃),fsys)

        # forward w₀ by recursion
        w₀ .= H[i-1]*w₀

        # temporarily use vbuffer to store accumulated v
        vbuffer .= v₀
        for j = 1:i-2
            vbuffer .+= rkdt.a[i,j] .* v̇[j]
        end

        # construct r₂₁
        fbuffer .= -T₂ᵢ(vbuffer)
        wbuffer.data .= w₀
        for j = 1:i-2
            wbuffer .+= rkdt.a[i,j] .* ċ[j]
        end
        fbuffer .+= B₂ᵢ(wbuffer)
        fbuffer .*= -1.0/rkdt.a[i,i-1]

        # construct r₂₂
        λbuffer .= gtiᵢ
        λbuffer .+= G₂ᵢ*vbuffer
        λbuffer .*= -1.0/rkdt.a[i,i-1]

        # construct r₁₁
        wbuffer .= r₁₁(w,tᵢ)

        # construct r₁₂
        vbuffer .= Fᵢ₋₁

        # solve the linear system
        ċ[i-1], v̇[i-1], f, λ = S\(wbuffer,vbuffer,fbuffer,λbuffer)

        # forward ċ[j] by recursion
        for j = 1:i-2
            ċ[j] .= H[i-1]*ċ[j]
        end

        # accumulate fluid vorticity up to the current stage
        w .= w₀
        for j = 1:i-1
            w .+=  rkdt.a[i,j] .* ċ[j]
        end

        # accumulate body velocity up to the current stage
        v .= v₀
        for j = 1:i-1
            v .+= rkdt.a[i,j] .* v̇[j]
        end

        # update vJ using updated v
        bd, vJ[i] = UpV(bd, v)

    end

    return tᵢ, (w, f), (qJ, v, λ), bd
end
