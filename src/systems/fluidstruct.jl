mutable struct FluidStruct{NX, NY, N}
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
    rk::RKParams

    # Operators
    "Laplacian operator"
    L::Fields.Laplacian{NX,NY}

    # Body coordinate data
    X̃::VectorData{N}

    # Pre-stored regularization and interpolation matrices (if present)
    Hmat::RegularizationMatrix
    Emat::InterpolationMatrix

    # Pre-allocated space for intermediate values
    Vb::VectorData{N}
    Fq::Edges{Primal, NX, NY}
    Ww::Edges{Dual, NX, NY}
    Qq::Edges{Dual, NX, NY}

end

function FluidStruct(dims::Tuple{Int, Int}, Re, Δx, Δt;
                       U∞ = (0.0, 0.0), X̃ = VectorData{0}(),
                       rk::RKParams=RK31)
    NX, NY = dims

    L = plan_laplacian((NX,NY),with_inverse=true)

    Vb = VectorData(X̃)
    Fq = Edges{Primal,NX,NY}()
    Ww   = Edges{Dual, NX, NY}()
    Qq  = Edges{Dual, NX, NY}()
    N = length(X̃)÷2

    # set interpolation and regularization matrices
    regop = Regularize(X̃,Δx,issymmetric=true)
    Hmat, Emat = RegularizationMatrix(regop,Vb,Fq)

    FluidStruct{NX, NY, N}(Re, U∞, Δx, Δt, rk, L, X̃, Hmat, Emat, Vb, Fq, Ww, Qq)
end

function Base.show(io::IO, sys::FluidStruct{NX,NY,N}) where {NX,NY,N}
    print(io, "Fluid-Structure Interaction system on a grid of size $NX x $NY")
end
