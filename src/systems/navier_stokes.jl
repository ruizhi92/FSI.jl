mutable struct NavierStokes{NX, NY, N, isstatic}  #<: System{Unconstrained}
    # Physical Parameters
    "Reynolds number"
    Re::Float64
    "Free stream velocities"
    U∞::Tuple{Float64, Float64}

    # Discretization
    "Grid spacing"
    Δx::Float64
    "Time step (used to determine integrating factor diffusion rate)"
    Δt::Float64
    "Runge-Kutta method"
    rk::TimeMarching.RKParams

    # Operators
    "Laplacian operator"
    L::Fields.Laplacian{NX,NY}

    # Body coordinate data, if present
    # if a static problem, these coordinates are in inertial coordinates
    # if a non-static problem, in their own coordinate systems
    X̃::VectorData{N}

    # Pre-stored regularization and interpolation matrices (if present)
    Hmat::RegularizationMatrix
    Emat::InterpolationMatrix

    # Scratch space

    ## Pre-allocated space for intermediate values
    Vb::VectorData{N}
    Fq::Edges{Primal, NX, NY}
    Ww::Edges{Dual, NX, NY}
    Qq::Edges{Dual, NX, NY}

    # Flags
    _isstore :: Bool

end

function NavierStokes(dims::Tuple{Int, Int}, Re, Δx, Δt;
                       U∞ = (0.0, 0.0), X̃ = VectorData{0}(),
                       isstore = false,
                       isstatic = true,
                       rk::TimeMarching.RKParams=TimeMarching.RK31)
    NX, NY = dims

    α = Δt/(Re*Δx^2)

    L = plan_laplacian((NX,NY),with_inverse=true)

    Vb = VectorData(X̃)
    Fq = Edges{Primal,NX,NY}()
    Ww   = Edges{Dual, NX, NY}()
    Qq  = Edges{Dual, NX, NY}()
    N = length(X̃)÷2

    if length(N) > 0 && isstore && isstatic
      # in this case, X̃ is assumed to be in inertial coordinates
      regop = Regularize(X̃,Δx;issymmetric=true)
      Hmat, Emat = RegularizationMatrix(regop,VectorData{N}(),Edges{Primal,NX,NY}())
    else
        regop = Regularize(X̃,Δx,issymmetric=true)
        Hmat, Emat = RegularizationMatrix(regop,Vb,Fq)
    end

    # should be able to set up time marching operator here...

    NavierStokes{NX, NY, N, isstatic}(Re, U∞, Δx, Δt, rk, L, X̃, Hmat, Emat, Vb, Fq, Ww, Qq, isstore)
end

function Base.show(io::IO, sys::NavierStokes{NX,NY,N,isstatic}) where {NX,NY,N,isstatic}
    print(io, "Navier-Stokes system on a grid of size $NX x $NY")
end


# Integrating factor -- rescale the time-step size
Fields.plan_intfact(Δt,w,sys::NavierStokes{NX,NY}) where {NX,NY} =
        Fields.plan_intfact(Δt/(sys.Re*sys.Δx^2),w)

# RHS of Navier-Stokes (non-linear convective term)
function TimeMarching.r₁(w::Nodes{Dual,NX,NY},t,sys::NavierStokes{NX,NY}) where {NX,NY}

  Ww = sys.Ww
  Qq = sys.Qq
  L = sys.L
  Δx⁻¹ = 1/sys.Δx

  cellshift!(Qq,curl(L\w)) # -velocity, on dual edges
  Qq.u .-= sys.U∞[1]
  Qq.v .-= sys.U∞[2]

  return Fields.rmul!(divergence(Qq∘cellshift!(Ww,w)),Δx⁻¹) # -∇⋅(wu)

end

# RHS of Navier-Stokes (non-linear convective term)
# U∞ represents free stream velocity, which is subtracted from the b.c.s from rigid body
function TimeMarching.r₁(w::Nodes{Dual,NX,NY},t,sys::NavierStokes{NX,NY},U∞::RigidBodyMotions.RigidBodyMotion) where {NX,NY}

  Ww = sys.Ww
  Qq = sys.Qq
  L = sys.L
  Δx⁻¹ = 1/sys.Δx

  cellshift!(Qq,curl(L\w)) # -velocity, on dual edges
  _,ċ,_,_,_,_ = U∞(t)
  Qq.u .-= real(ċ)
  Qq.v .-= imag(ċ)

  return Fields.rmul!(divergence(Qq∘cellshift!(Ww,w)),Δx⁻¹) # -∇⋅(wu)

end


# Constraint operators, using stored regularization and interpolation operators
# B₁ᵀ = CᵀEᵀ, B₂ = -ECL⁻¹
TimeMarching.B₁ᵀ(f,sys::NavierStokes{NX,NY,N,C}) where {NX,NY,N,C} = Curl()*(sys.Hmat*f)
TimeMarching.B₂(w,sys::NavierStokes{NX,NY,N,C}) where {NX,NY,N,C} = -(sys.Emat*(Curl()*(sys.L\w)))


# RHS of Navier-Stokes (non-linear convective term) for constant U∞
"""
Boundary condition at body points for fluids with uniform constrant free stream U∞
r₂ takes in motion as linear velocity of body points in the inertial frame and substract
the free stream velocity.
"""
function TimeMarching.r₂(u::Nodes{Dual,NX,NY},t,sys::NavierStokes{NX,NY,N,false},
            motion::Array{Float64,2}) where {NX,NY,N}
    ΔV = VectorData(motion)
    ΔV.u .-= sys.U∞[1]
    ΔV.v .-= sys.U∞[2]
# println("motion ", ΔV.v, "\n")
    return ΔV
end

"""
Boundary condition at body points for fluids with uniform time-varying free stream U∞(t)
r₂ takes in motion as linear velocity of body points in the inertial frame and substract
the free stream velocity.
"""
function TimeMarching.r₂(u::Nodes{Dual,NX,NY},t,sys::NavierStokes{NX,NY,N,false},
            motion::Array{Float64,2},U∞::RigidBodyMotions.RigidBodyMotion) where {NX,NY,N}
    ΔV = VectorData(motion)
    _,ċ,_,_,_,_ = U∞(t)
    ΔV.u .-= real(ċ)
    ΔV.v .-= imag(ċ)
    return ΔV
end

"""
plan_constraints takes in coordinates of body points in inertial frame and construct
B₁ᵀ, B₂ operator from it.
"""
function TimeMarching.plan_constraints(u::Nodes{Dual,NX,NY},t,sys::NavierStokes{NX,NY,N,false},
            coord::Array{Float64,2}) where {NX,NY,N}
# println("coord ", coord,"\n")
    X = VectorData(coord)
    regop = Regularize(X,sys.Δx;issymmetric=true)
    if sys._isstore
      Hmat, Emat = RegularizationMatrix(regop,VectorData{N}(),Edges{Primal,NX,NY}())
      sys.Hmat = Hmat
      sys.Emat = Emat
      return f->TimeMarching.B₁ᵀ(f,sys), w->TimeMarching.B₂(w,sys)
    else
      return f->TimeMarching.B₁ᵀ(f,regop,sys), w->TimeMarching.B₂(w,regop,sys)
    end
end
