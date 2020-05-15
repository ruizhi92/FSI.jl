@testset "Construct Block system with LinearMaps" begin
    using IterativeSolvers

    # using Julia's default \ operator
    global m1 = [ 0.533442  0.367231  0.217876;
     0.937095  0.838996  0.38914;
     0.399906  0.409811  0.571862]
    global m2 = [ 0.515455  0.912171; 0.197663  0.165777; 0.3626 0.735575]
    global m3 = [0.763043  0.36013   0.779911; 0.676889  0.262959  0.0735239]
    global m4 = [0.981832  0.695422; 0.906051  0.432737]
    b = [0.71537;0.37217;0.728864;0.72988;0.881618]
    result1 = [m1 m2;m3 m4]\b

    # construct block system using LinearMap and solve with IterativeSolvers
    function al(x::AbstractVector)
        f(x::AbstractVector) =  m1*x
        g(x::AbstractVector) = m2*x
        m(x::AbstractVector) = m3*x
        n(x::AbstractVector) = m4*x
        # has to do the matrix mul manually to return 1d output
        out = zeros(5)
        out[1:3] = f(x[1:3]) + g(x[4:5])
        out[4:5] = m(x[1:3]) + n(x[4:5])
        return out
    end
    S = LinearMap(al,5;ismutating=false)
    result2 = gmres(S,b)

    @test isapprox(norm(result1-result2),0.0,atol=1e-10)
end


@testset "SaddlePointSystems for 1d body" begin
    # construct a random block system
    # block size
    m = 20
    n = 8
    p = 10
    q = 6

    # operators
    A = rand(m,m); A⁻¹ = inv(A)
    B₁ᵀ = rand(m,n)
    B₂ = rand(n,m)
    M = rand(p,p)
    G₁ᵀ = rand(p,q)
    G₂ = rand(q,p)
    T₁ᵀ = rand(p,n)
    T₂ = rand(n,p)

    # zero matrix
    Omp = zeros(m,p)
    Omq = zeros(m,q)
    Onn = zeros(n,n)
    Onq = zeros(n,q)
    Opm = zeros(p,m)
    Oqm = zeros(q,m)
    Oqn = zeros(q,n)
    Oqq = zeros(q,q)

    # saddle point system
    S = [A B₁ᵀ Omp Omq; B₂ Onn -T₂ Onq; Opm T₁ᵀ M G₁ᵀ; Oqm Oqn G₂ Oqq]

    # rhs
    rċ = rand(m)
    rf = rand(n)
    ru̇ = rand(p)
    rλ = rand(q)
    b = [rċ;rf;ru̇;rλ];

    # test using Julia's default \ solver
    @test isapprox(norm(S*(S\b)-b),0.0,atol=1e-10)
    @test rank(S) == m+n+p+q

    # construct saddle point system
    ċ = zeros(m)
    f = zeros(n)
    u̇ = zeros(p)
    λ = zeros(q)
    St = SaddleSystem1d((ċ, f, u̇, λ), (A⁻¹, B₁ᵀ, B₂),
                      (M, G₁ᵀ, G₂), (x->T₁ᵀ*x, x->T₂*x))
    bt = (rċ, ru̇, rf, rλ);

    # test using SaddlePointSystems
    aa,bb,cc,dd = St\bt
    x = [aa;bb;cc;dd]
    @test isapprox(norm(S*x-b),0.0,atol=1e-10)

end

@testset "SaddlePointSystems for 2d body" begin
    # construct a random block system
    # block size
    m = 10
    n = 8
    p = 6
    q = 4

    # operators
    A = rand(m,m); A⁻¹ = inv(A)
    B₁ᵀ = rand(m,n)
    B₂ = rand(n,m)
    M = rand(p,p)
    G₁ᵀ = rand(p,q)
    G₂ = rand(q,p)
    T₁ᵀ = rand(p,n)
    T₂ = rand(n,p)

    # zero matrix
    Omp = zeros(m,p)
    Omq = zeros(m,q)
    Onn = zeros(n,n)
    Onq = zeros(n,q)
    Opm = zeros(p,m)
    Oqm = zeros(q,m)
    Oqn = zeros(q,n)
    Oqq = zeros(q,q)

    # fictitious fluid
    ρb = 2.0
    Mf = 1.0/ρb*M

    # saddle point system
    S = [A B₁ᵀ -B₁ᵀ*T₂*Mf Omq; B₂ Onn -T₂ Onq; Opm T₁ᵀ M G₁ᵀ; Oqm Oqn G₂ Oqq]

    # rhs
    rċ = rand(m)
    rf = rand(n)
    ru̇ = rand(p)
    rλ = rand(q)
    b = [rċ;rf;ru̇;rλ];

    # test using Julia's default \ solver
    @test isapprox(norm(S*(S\b)-b),0.0,atol=1e-10)
    @test rank(S) == m+n+p+q

    # construct saddle point system
    ċ = zeros(m)
    f = zeros(n)
    u̇ = zeros(p)
    λ = zeros(q)
    St = SaddleSystem2d((ċ, f, u̇, λ), (A⁻¹, B₁ᵀ, B₂),
                      (M, G₁ᵀ, G₂), (x->T₁ᵀ*x, x->T₂*x);
                      ρb=ρb)
    bt = (rċ, ru̇, rf, rλ);

    # test using SaddlePointSystems
    aa,bb,cc,dd = St\bt
    x = [aa;bb;cc;dd]
    @test isapprox(norm(S*x-b),0.0,atol=1e-10)

end
