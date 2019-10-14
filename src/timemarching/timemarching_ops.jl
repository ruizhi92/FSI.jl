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

  interpolate!(Qq,curl(L\w)) # -velocity, on dual edges
  Qq.u .-= sys.U∞[1]
  Qq.v .-= sys.U∞[2]

  return rmul!(divergence(Qq∘interpolate!(Ww,w)),Δx⁻¹) # -∇⋅(wu)

end

# RHS of Navier-Stokes (non-linear convective term)
# U∞ represents free stream velocity, which is subtracted from the b.c.s from rigid body
function TimeMarching.r₁(w::Nodes{Dual,NX,NY},t,sys::FluidStruct{NX,NY},U∞::RigidBodyMotions.RigidBodyMotion) where {NX,NY}

  Ww = sys.Ww
  Qq = sys.Qq
  L = sys.L
  Δx⁻¹ = 1/sys.Δx

  interpolate!(Qq,curl(L\w)) # -velocity, on dual edges
  _,ċ,_,_,_,_ = U∞(t)
  Qq.u .-= real(ċ)
  Qq.v .-= imag(ċ)

  return rmul!(divergence(Qq∘interpolate!(Ww,w)),Δx⁻¹) # -∇⋅(wu)

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
    return Curl()*(sys.Hmat*f)
end

function TimeMarching.B₂(w,sys::FluidStruct{NX,NY,N}) where {NX,NY,N}
    return -(sys.Emat*(Curl()*(sys.L\w)))
end


# wrap functions from Dyn3d
function TimeMarching.M(bd::BodyDyn)
    return HERKFuncM(bd.sys)
end

# Whenever buoyancy needs to be accounted for, it means that the body has finite
# volume in x-y space. Thus the body must account for fictitious fluid effect as
# well. These two effects got accounted for by setting ρ in theory to ρ-1 in code
function TimeMarching.F(bd::BodyDyn)
    f_exi = zeros(Float64,bd.sys.nbody,6)
    return HERKFuncf(bd.bs, bd.js, bd.sys, f_exi; influid=true, bodydim=1)
end

function TimeMarching.F(bd::BodyDyn,ρb::Float64)
    f_exi = zeros(Float64,bd.sys.nbody,6)
    return HERKFuncf(bd.bs, bd.js, bd.sys, f_exi; influid=true, bodydim=2, ρb=ρb)
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


"""
    BodyGridToVectorData(bgs::Vector{BodyGrid},datatype::String;plane::Vector{Int}=[1,2])

Given BodyGrid structure, extract motion or coordinates and return as VectorData,
arranged by Lagrangian points order used in ViscousFlow. By default the body
constructed in Dyn3d is in x-y plane(1d body like a flat plate). If we want to have
a 2d body for FSI such as a cylinder, set plane=[1,3].
"""
function BodyGridToVectorData(bgs::Vector{BodyGrid},datatype::String;plane::Vector{Int}=[1,2])
    if datatype == "coord"
        coord = hcat(bgs[1].q_i...)'[:,plane]
        for i = 2:length(bgs)
            coord = [coord; hcat(bgs[i].q_i...)'[:,plane]]
        end
        return VectorData(coord)
    elseif datatype == "motion"
        motion = hcat(bgs[1].v_i...)'[:,plane]
        for i = 2:length(bgs)
            motion = [motion; hcat(bgs[i].v_i...)'[:,plane]]
        end
        return VectorData(motion)
    end
end


"""
    T₁ᵀ(bd::BodyDyn,bgs::Vector{BodyGrid},f::VectorData,Δx::Float64;plane::Vector{Int}=[1,2])

T₁ᵀ takes in 2d x-y plane fluid force per unit are f of all body grid points in VectorData form
and calculate integrated force on each body in 1d Vector form(line up dimension of nbody*6_dof).

Note that force from fluid solver need to be multiplied by Δx^2 before going into body solver.
"""
function TimeMarching.T₁ᵀ(bd::BodyDyn,bgs::Vector{BodyGrid},f::VectorData,Δx::Float64;plane::Vector{Int}=[1,2])

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
        bgs[b_cnt].f_ex3d[i-ref][plane[1]] = fbuffer.u[i]
        bgs[b_cnt].f_ex3d[i-ref][plane[2]] = fbuffer.v[i]
        bgs[b_cnt].m_ex3d[i-ref] .= cross(bgs[b_cnt].q_i[i-ref],bgs[b_cnt].f_ex3d[i-ref])
    end

    # Integrate total forces from all body points on a body
    # then transform the external forces from inertial frame to body frame
    bgs = IntegrateBodyGridDynamics(bd,bgs)
    for i = 1:bd.sys.nbody
        f_exis[i,:] = bgs[i].f_ex6d
        f_exis[i,:] = bd.bs[i].Xb_to_i'*f_exis[i,:]
    end

    return (f_exis')[:]
end


"""
    T₂(bd::BodyDyn,bgs::Vector{BodyGrid},u::Array{Float64,1};plane::Vector{Int}=[1,2])

T₂ takes in the lined-up 6d spatial velocity of all body in 1d Vector form and
acquire body grid points's x-y plane 2d velocity to return VectorData
"""
function TimeMarching.T₂(bd::BodyDyn,bgs::Vector{BodyGrid},u::Array{Float64,1};plane::Vector{Int}=[1,2])
    @get bd (bs, sys)
    @get sys.pre_array (la_tmp1,la_tmp2)

    v_temp = zeros(Float64,6)
    la_tmp1 .= 0.0
    la_tmp2 .= 0.0

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
            X_ref = TransMatrix([zeros(Float64,3);b.x_i],la_tmp1,la_tmp2)
        end
    end

    for i = 1:length(bgs)
        b = bs[bgs[i].bid]
        for j = 1:bgs[i].np
            v_temp = bs[i].v + [zeros(Float64, 3); cross(bs[i].v[1:3],bgs[i].points[j])]
            bgs[i].v_i[j] = (X_ref*b.Xb_to_i*v_temp)[4:6]
        end
    end

    motion = BodyGridToVectorData(bgs,"motion";plane=plane)
    return motion
end


"""
    getX̃(bd::BodyDyn,bgs::Vector{BodyGrid};plane::Vector{Int}=[1,2])

getX̃ use body bs[i].x_i and acquire body grid points's coordinates by the j-th
q_i in body points of a body = bs[i].x_i + Xb_to_i*points[j], return VectorData.
"""
function TimeMarching.getX̃(bd::BodyDyn,bgs::Vector{BodyGrid};plane::Vector{Int}=[1,2])
    @get bd (bs, )

    for i = 1:length(bgs)
        b = bs[bgs[i].bid]
        for j = 1:bgs[i].np
            q_temp = [zeros(Float64, 3); bgs[i].points[j]]
            q_temp = [zeros(Float64, 3); b.x_i] + b.Xb_to_i*q_temp
            bgs[i].q_i[j] = q_temp[4:6]
        end
    end

    coord = BodyGridToVectorData(bgs,"coord";plane=plane)
    return coord
end
