# Basic operators for any FluidStruct system

# Integrating factor -- rescale the time-step size
Fields.plan_intfact(Δt,w,sys::FluidStruct{NX,NY}) where {NX,NY} =
        Fields.plan_intfact(Δt/(sys.Re*sys.Δx^2),w)

# RHS of Navier-Stokes (non-linear convective term)
function TimeMarching.r₁(w::Nodes{Dual,NX,NY},t,sys::FluidStruct{NX,NY}) where {NX,NY}

  Ww = sys.Ww
  Qq = sys.Qq
  L = sys.L
  Δx⁻¹ = 1/sys.Δx

  shift!(Qq,curl(L\w)) # -velocity, on dual edges
  Qq.u .-= sys.U∞[1]
  Qq.v .-= sys.U∞[2]

  return scale!(divergence(Qq∘shift!(Ww,w)),Δx⁻¹) # -∇⋅(wu)

end

# RHS of Navier-Stokes (non-linear convective term)
# U∞ represents free stream velocity, which is subtracted from the b.c.s from rigid body
function TimeMarching.r₁(w::Nodes{Dual,NX,NY},t,sys::FluidStruct{NX,NY},U∞::RigidBodyMotions.RigidBodyMotion) where {NX,NY}

  Ww = sys.Ww
  Qq = sys.Qq
  L = sys.L
  Δx⁻¹ = 1/sys.Δx

  shift!(Qq,curl(L\w)) # -velocity, on dual edges
  _,ċ,_,_,_,_ = U∞(t)
  Qq.u .-= real(ċ)
  Qq.v .-= imag(ċ)

  return scale!(divergence(Qq∘shift!(Ww,w)),Δx⁻¹) # -∇⋅(wu)

end

# uniform flow velocity
function TimeMarching.U_inf(w::Nodes{Dual,NX,NY},sys::FluidStruct{NX,NY,N}) where {NX,NY,N}
    ΔV = VectorData(sys.X̃) # create a struct ΔV with same shape of sys.X̃ and initialize it with 0
    ΔV.u .= sys.U∞[1]
    ΔV.v .= sys.U∞[2]
    return ΔV
end

# Constraint operators, using stored regularization and interpolation operators
# B₁ᵀ = CᵀEᵀ, B₂ = -ECL⁻¹
function TimeMarching.B₁ᵀ(f,sys::FluidStruct{NX,NY,N}) where {NX,NY,N}
    return Curl()*(get(sys.Hmat)*f)
end

function TimeMarching.B₂(w,sys::FluidStruct{NX,NY,N}) where {NX,NY,N}
    return -(get(sys.Emat)*(Curl()*(sys.L\w)))
end


# wrap functions from Dyn3d
function TimeMarching.M⁻¹(bd::BodyDyn)
    return HERKFuncM⁻¹(bd.sys)
end

function TimeMarching.F(bd::BodyDyn)
    f_exi = zeros(Float64,bd.sys.nbody,6)
    return HERKFuncf(bd.bs, bd.js, bd.sys, f_exi)
end

function TimeMarching.G₁ᵀ(bd::BodyDyn)
    return HERKFuncGT(bd.bs, bd.sys)
end

function TimeMarching.G₂(bd::BodyDyn)
    return HERKFuncG(bd.bs, bd.sys)
end

function TimeMarching.gti(bd::BodyDyn, t::Float64)
    return HERKFuncgti(bd.js, bd.sys, t)
end

function TimeMarching.UpP(bd::BodyDyn, qJ::Vector{Float64})
    bd.bs, bd.js, bd.sys = UpdatePosition!(bd.bs, bd.js, bd.sys, qJ)
    return bd
end

function TimeMarching.UpV(bd::BodyDyn, v::Vector{Float64})
    bd.bs, bd.js, bd.sys, vJ = UpdateVelocity!(bd.bs, bd.js, bd.sys, v)
    return bd, vJ
end

# define functions for fsi operators

# T₁ᵀ takes in 2d x-y plane fluid force f of all body grid points in VectorData form
# and calculate integrated force on each body in 1d Array form(line up dimension of nbody*6_dof)
function TimeMarching.T₁ᵀ(bd::BodyDyn, bgs::Vector{BodyGrid}, f::VectorData, Δx::Float64)
    # Note that force from fluid solver need to be multiplied by Δx^2 before going into body solver
    fbuffer = deepcopy(f)
    fbuffer .*= Δx^2

    # Assign f on body grid points to BodyGrid structure of each body
    b_cnt = 1
    ref = 0
    f_exis = zeros(bd.sys.nbody,6)
    np_total = round(Int,length(fbuffer)/2)
    for i = 1:np_total
        # move to the next bgs if i exceed bgs[b_cnt].np
        if i > ref + bgs[b_cnt].np
            ref += bgs[b_cnt].np
            b_cnt += 1
        end
        bgs[b_cnt].f_ex3d[i-ref][1] = fbuffer.u[i]
        bgs[b_cnt].f_ex3d[i-ref][2] = fbuffer.v[i]
    end

    # Integrate total forces from all body points on a body
    # then transform the external forces from inertial frame to body frame
    # Mind Newton's 3rd law(the negative sign)
    bgs = IntegrateBodyGridDynamics(bd,bgs)
    for i = 1:bd.sys.nbody
        f_exis[i,:] = -bgs[i].f_ex6d
        f_exis[i,:] = bd.bs[i].Xb_to_i'*f_exis[i,:]
    end
    
    return (f_exis')[:]
end

# T₂ takes in the lined-up 6d spatial velocity of all body in 1d Array form
# and acquire body grid points's x-y plane 2d velocity to return VectorData
function TimeMarching.T₂(bd::BodyDyn, bgs::Vector{BodyGrid}, u::Array{Float64,1})
    @get bd (bs, sys)

    count = 0
    for i = 1:sys.nbody
        bs[i].v = u[count+1:count+6]
        count += 6
    end

    # the j-th v_i in body points of a body is calculated by transferring to
    # a coordinate that sits at the beginning point of the first body but with zero angle.
    X_ref = zeros(Float64,6)
    for i = 1:length(bgs)
        b = bs[bgs[i].bid]
        if b.bid == 1
            X_ref = TransMatrix([zeros(Float64,3);b.x_i])
        end
    end

    for i = 1:length(bgs)
        b = bs[bgs[i].bid]
        for j = 1:bgs[i].np
            v_temp = bs[i].v + [zeros(Float64, 3); cross(bs[i].v[1:3],bgs[i].points[j])]
            bgs[i].v_i[j] = (X_ref*b.Xb_to_i*v_temp)[4:6]
        end
    end

    motion = hcat(bgs[1].v_i...)'[:,[1,2]]
    for i = 2:length(bgs)
        motion = [motion[1:end-1,:]; hcat(bgs[i].v_i...)'[:,[1,2]]]
    end
    return VectorData(motion)
end

# getX̃ use body bs[i].x_i and acquire body grid points's coordinates by
# the j-th q_i in body points of a body = bs[i].x_i + Xb_to_i*points[j], return VectorData
function TimeMarching.getX̃(bd::BodyDyn, bgs::Vector{BodyGrid})
    @get bd (bs, )

    for i = 1:length(bgs)
        b = bs[bgs[i].bid]
        for j = 1:bgs[i].np
            q_temp = [zeros(Float64, 3); bgs[i].points[j]]
            q_temp = [zeros(Float64, 3); b.x_i] + b.Xb_to_i*q_temp
            bgs[i].q_i[j] = q_temp[4:6]
        end
    end

    coord = hcat(bgs[1].q_i...)'[:,[1,2]]
    for i = 2:length(bgs)
        coord = [coord[1:end-1,:]; hcat(bgs[i].q_i...)'[:,[1,2]]]
    end
    return VectorData(coord)
end
