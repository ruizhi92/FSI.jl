# This test only checks for function syntax mistakes using a falling plate example
@testset "Operators in FSI system" begin
    # problem dimension
    ndim = 2
    # numerical params
    tf = 6
    dt = 1e-3
    scheme = "Liska"
    st = 3
    tol = 1e-4
    num_params = NumParams(tf, dt, scheme, st, tol)
    # gravity
    gravity = [0., -1.0, 0.]
    # set up system config info
    config_system = ConfigSystem(ndim, gravity, num_params)
    # set up bodys
    nbody = 1
    config_body = ConfigBody(nbody, 4,
       [0. 0.; 1. 0.; 1. 1.0/nbody; 0. 1.0/nbody], 1.0)
    config_bodys = fill(config_body, nbody)
    # set up joints
    njoint = nbody
    config_joints = Vector{ConfigJoint}(undef,njoint)
    # set the first passive joint with no stiff and damp
    dof_1 = Dof(5, "passive", 0., 0., Motions())
    config_joints[1] = ConfigJoint(njoint, "custom",
        [0.,0.,0.,1.0,1.0,0.], zeros(Float64,6), 0, [dof_1], [0.])

    bs, js, bsys = BuildChain(config_bodys, config_joints, config_system)
    bd = BodyDyn(bs, js, bsys)

    bd, soln₀ = InitSystem!(bd)

    @get bd (bs, js, sys)
    qJ_dim = sys.ndof
    λ_dim = sys.ncdof_HERK
    u = zeros(qJ_dim)
    λ = zeros(λ_dim);

    Re = 200 # Reynolds number
    U = 1.0 # Free stream velocity
    U∞ = (U, 0.0)

    nx = 152; ny = 102;
    Ly = 2.0;
    Δx = Ly/(ny-2);
    w₀ = Nodes(Dual,(nx,ny))
    xg, yg = coordinates(w₀,dx=Δx)

    w₀ .= 0.0;

    # bgs short for body grid system
    bgs = GenerateBodyGrid(bd; np=DetermineNP(nbody, Δx))
    bgs = CutOut2d(bd,bgs);

    bgs = AcquireBodyGridKinematics(bd,bgs)
    coord_init = hcat(bgs[1].q_i...)'[:,[1,2]]
    for i = 2:length(bgs)
        coord_init = [coord_init[1:end-1,:]; hcat(bgs[i].q_i...)'[:,[1,2]]]
    end

    X̃ = VectorData(coord_init)
    f = VectorData(X̃);

    Δt = 0.01
    t = Δt
    fsys = Systems.FluidStruct((nx,ny),Re,Δx,Δt,U∞ = U∞, X̃ = X̃, rk=RK31)
    N = length(X̃)÷2

    # Integrating factor not depending on time
    dc = 1.0
    H = plan_intfact(dc*Δt,w₀)

    # Operators at time t_i-1
    function TimeMarching.F(bd::BodyDyn)
        f_exi = zeros(Float64,bd.sys.nbody,6)
        return HERKFuncf(bd.bs, bd.js, bd.sys, f_exi)
    end
    # Body operators
    Mᵢ₋₁ = M⁻¹(bd)
    Fᵢ₋₁ = F(bd)
    G₁ᵀᵢ₋₁ = G₁ᵀ(bd)
    # Fluid operators
    B₁ᵀᵢ₋₁ = f -> B₁ᵀ(f,fsys)
    # FSI operators
    T₁ᵀᵢ₋₁ = f -> T₁ᵀ(bd,bgs,f,Δx);

    # Integrate joints qJ and update bs and js in time marching
    qJ = deepcopy(soln₀.qJ)
    vJ = deepcopy(qJ)
    # get vJ from v to update qJ
    bd, vJ = UpV(bd, soln₀.v)
    dc = 0.5
    qJ += Δt*0.5*vJ
    # use new qJ to update system position b.x_i
    bd = UpP(bd, qJ);

    # Operators at time t_i with updated body coordinates
    # Body operators
    G₂ᵢ = G₂(bd)
    gtiᵢ = gti(bd, t)
    # Update current body points coordinates and fill in fsys.X̃
    X̃ = getX̃(bd, bgs)
    # Fluid operators
    B₂ᵢ = w -> B₂(w,fsys)
    # FSI operators
    T₂ᵢ = u -> T₂(bd,bgs,u);

    # Fill operators into saddle point system
    S = SaddleSystem((w₀,u,f,λ),
                 (H, B₁ᵀᵢ₋₁, B₂ᵢ),
                 (Mᵢ₋₁, G₁ᵀᵢ₋₁, G₂ᵢ),
                 (T₁ᵀᵢ₋₁, T₂ᵢ))

    # Test all operators
    S.A⁻¹(w₀)
    S.B₂(w₀)
    S.T₂(u)
    S.M⁻¹(u)
    S.G₂(u)
    S.A⁻¹B₁ᵀ(f)
    S.M⁻¹T₁ᵀ(f)
    S.M⁻¹G₁ᵀ(λ);

    # Assign right hand side
    rw₀ = deepcopy(w₀)
    rw₀.data .= 0.0
    ru = ones(Float64,6*bd.sys.nbody)
    rf = deepcopy(f)
    rλ =zeros(size(λ))
    rhs = (rw₀, ru, rf, rλ);

    # Solve the system
    w₀, u, f, λ = S\rhs;

    @test isapprox(sum(f.u),0.0,atol=1e-10)
end
