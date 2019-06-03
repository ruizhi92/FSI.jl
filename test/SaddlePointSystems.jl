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


@testset "SaddlePointSystems" begin
    using ViscousFlow.Fields
    import Base: +,*

    # construct a random block system
    A = rand(80,80)
    A⁻¹ = inv(A)
    B₁ᵀ = rand(80,12)
    B₂ = rand(12,80)
    M = rand(12,12)
    M⁻¹ = inv(M)
    G₁ᵀ = rand(12,12)
    G₂ = rand(12,12)
    T₁ᵀ = rand(12,12)
    T₂ = rand(12,12)
    O32 = zeros(80,12)
    O23 = zeros(12,80)
    O22 = zeros(12,12)
    S = [A O32 B₁ᵀ O32; O23 M T₁ᵀ G₁ᵀ; B₂ -T₂ O22 O22; O23 G₂ O22 O22]
    rċ = rand(80)
    ru̇ = rand(12)
    rf = rand(12)
    rλ = rand(12)
    b = [rċ;ru̇;rf;rλ];

    # test using Julia's default \ solver
    @test isapprox(norm(S*(S\b)-b),0.0,atol=1e-10)
    @test rank(S) == 116

    # construct saddle point system
    ċ = zeros(80)
    u̇ = zeros(12)
    f = zeros(12)
    λ = zeros(12)
    St = SaddleSystem((ċ, u̇, f, λ), (A⁻¹, B₁ᵀ, B₂),
                      (M⁻¹, G₁ᵀ, G₂), (x->T₁ᵀ*x, x->T₂*x))
    bt = (rċ, ru̇, rf, rλ);

    # test using SaddlePointSystems
    aa,bb,cc,dd = St\bt
    x = [aa;bb;cc;dd]
    @test isapprox(norm(S*x-b),0.0,atol=1e-10)

end
