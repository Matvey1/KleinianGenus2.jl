using Test, KleinianGenus2, LinearAlgebra, Random, Polynomials, QuadGK
const K = KleinianGenus2
Random.seed!(20260923)
const R64 = [-2.,-1.,0.,1.,3.]
period_error(C,z) = norm(K.normalizeToPeriods(C.Periods,z))
curve_error(C,P) = abs(P[2]^2-C.F(P[1]))/max(1,abs(C.F(P[1])))

@testset "API and input regressions" begin
    C = Genus2WCurve(R64,source="Roots",is_real=true)
    @test fieldnames(typeof(C)) == (:F,:Roots,:Periods,:Eta,:KleinianFunc,
        :KleinianFuncSigmaSquared,:KleinianFuncWeight2,:JacobiInversion,:AbelMap)
    @test all(f -> getfield(C,f) isa Function,fieldnames(typeof(C))[5:end])
    cfs = [-1,0,1,2,3,4]; snapshot=copy(cfs)
    CI=Genus2WCurve(cfs)
    @test cfs == snapshot
    @test CI.F == Polynomial(complex.(float.(cfs)))
    @test Genus2WCurve(Polynomial(cfs),source="Polynomial").F == CI.F
    cfloat=[-1.,0.,1.,2.,3.,4.]; snap=copy(cfloat)
    Genus2WCurve(cfloat)
    @test isequal(cfloat,snap)
    @test eltype(Genus2WCurveDivisor((1,0)).P1) == ComplexF64
    @test eltype(Genus2WCurveDivisor((1.,0.),(big"2",big"3")).P1) == Complex{BigFloat}
    @test isempty(Genus2WCurveDivisor((),(1.,0.)).P2)
    @test occursin("0",sprint(show,Genus2WCurveDivisor()))
    @test C.AbelMap(Genus2WCurveDivisor()) == zeros(2)
    @test isempty(C.JacobiInversion(zeros(2)).P1)
    @test !isempty(C.JacobiInversion([1e-8,5e-9]).P1)
    @test isempty(C.JacobiInversion(C.Periods[:,1]).P1)
    @test all(isfinite,Genus2WCurve(R64,source="Roots",is_real=true,e=10).KleinianFunc([.1,.2]))
    for bad in (0,-1,Inf,NaN)
        @test_throws ArgumentError Genus2WCurve(R64,source="Roots",e=bad)
    end
    @test_throws ArgumentError Genus2WCurve(R64,source="Roots",n=0)
    @test_throws ArgumentError Genus2WCurve(R64,source="Roots",n=1.5)
    @test_throws ArgumentError Genus2WCurve(R64 .- im,source="Roots",is_real=true)
    @test_throws ArgumentError Genus2WCurve([0,0,1,2,3],source="Roots")
    @test_throws ArgumentError Genus2WCurve([0,1,2,3,Inf],source="Roots")
    @test_throws ArgumentError C.AbelMap(Genus2WCurveDivisor((4.,1.)))
    @test_throws ArgumentError C.KleinianFuncWeight2([1.])
    @test_throws ArgumentError C.KleinianFuncWeight2([NaN,1.])
    @test_throws ArgumentError Genus2WCurve(R64,source="Roots",triple_of_discs=((0,1),(0,2)))
    @test_logs (:warn,r"Richelot iteration did not reach") Genus2WCurve(R64,source="Roots",n=1,is_real=false)
end

@testset "Disc construction and membership" begin
    for _ in 1:1000
        R=10 .* (rand(ComplexF64,5) .- (.5+.5im))
        discs=K.generateTripleOfDiscs(R)
        @test K.verifySortingMethod(R,discs)==0
        sm=K.GenerateSortingMethod(discs)(R)
        @test sort(length.(sm))==[1,2,2]
    end
    @test !K.belongs_to_circle(1,K.Circle(0,1))
    @test !K.belongs_to_circle(1,K.Circle(0,-1))
    @test K.check_separation([K.Circle(0,-10),K.Circle(0,-20),K.Circle(0,1)],R64)!=0
    @test K.check_separation([K.Circle(0,NaN),K.Circle(3,1),K.Circle(0,-10)],R64)!=0
    setprecision(256) do
        R=Complex{BigFloat}.([1+2im,-9+13im,3-7im,-11+4im,42+42im])
        @test K.verifySortingMethod(R,K.generateTripleOfDiscs(R))==0
    end
end

@testset "Independent quadratures and identities ($bits bits)" for bits in (53,128,256,512)
    setprecision(bits) do
        T=bits==53 ? Float64 : BigFloat
        R=T.([-2,-1,0,1,3]); e=100*eps(T); tol=100000*eps(T)
        C=Genus2WCurve(R,source="Roots",is_real=true,e=e,n=40)
        quad(f,a,b)=quadgk(f,a,b;atol=e,rtol=zero(T),order=max(7,bits÷8))[1]
        x=T(4); y=sqrt(C.F(x)); D=Genus2WCurveDivisor((x,y))
        integral=[quad(u->u^j/sqrt(prod(1-r*u^2 for r in R)),zero(T),inv(sqrt(x))) for j in (2,0)]
        @test period_error(C,C.AbelMap(D)+integral)<tol
        @test period_error(C,C.AbelMap(D,true)+integral)<tol
        @test period_error(C,C.AbelMap(Genus2WCurveDivisor((x,y),(x,y)))+2integral)<tol
        @test iszero(norm(C.AbelMap(Genus2WCurveDivisor((x,y),(x,-y)))))
        x2=T(5);DD=Genus2WCurveDivisor((x,y),(x2,-sqrt(C.F(x2))))
        interval=[quad(t->t^j/sqrt(C.F(t)),x,x2) for j in (0,1)]
        @test period_error(C,C.AbelMap(DD)+interval)<tol
        # Remove square-root endpoint singularities analytically.
        a,b=R[2],R[3]; rest=R[[1,4,5]]
        W=[quad(t->begin xx=(a+b)/2+(b-a)/2*cos(t);xx^j/sqrt(complex(-prod(xx-r for r in rest))) end,zero(T),T(pi)) for j in (0,1)]
        @test min(norm(W-C.Periods[:,1]),norm(W+C.Periods[:,1]))<tol
        z=Complex{T}.([T(12)/100+im*T(2)/100,T(5)/100-im*T(3)/100])
        DD=C.JacobiInversion(z)
        @test maximum(curve_error(C,P) for P in (DD.P1,DD.P2))<tol
        @test period_error(C,C.AbelMap(DD)-z)<tol
        A=C.KleinianFunc(z);B=C.KleinianFuncSigmaSquared(z);S=C.KleinianFuncWeight2(z)
        @test abs(A[1]^2-B[1])<tol
        @test norm(collect(A[2:end]).-collect(B[2:end]))<tol
        @test norm(C.KleinianFuncWeight2(-z)[1]-S[1])<tol
        @test norm(C.KleinianFuncWeight2(-z)[2]+S[2])<tol
        @test abs(C.KleinianFunc(-z)[1]+A[1])<tol
        @test abs(K.KummerEquationCheck(S[1]/norm(S[1]),C.F))<tol
        for j in 1:4
            shifted=C.KleinianFuncSigmaSquared(z+C.Periods[:,j])
            @test norm(collect(shifted[4:6])-collect(B[4:6]))<100tol
            @test norm(collect(shifted[2:3])-collect(B[2:3])-C.Eta[:,j])<100tol
        end
        h=bits==53 ? T(1e-5) : sqrt(sqrt(e))
        for j in 1:2
            v=h .* [j==1,j==2]
            fd=(C.KleinianFuncWeight2(z+v)[1]-C.KleinianFuncWeight2(z-v)[1])/(2h)
            @test norm(fd-S[j+1])<100*h^2
        end
    end
end

@testset "Special divisors" begin
    C=Genus2WCurve(R64,source="Roots",is_real=true)
    for x in R64
        z=C.AbelMap(Genus2WCurveDivisor((x,0.)))
        @test all(isfinite,z)
        @test period_error(C,2z)<1e-12
    end
    C=Genus2WCurve([-3.,-2.,-1.,1.,2.],source="Roots",is_real=true)
    D=Genus2WCurveDivisor((0.,sqrt(C.F(0.))))
    z=C.AbelMap(D)
    S=C.KleinianFuncWeight2(z)
    @test abs(S[1][1])/norm(S[1])<1e-10
    @test abs(S[1][3]/S[1][2])<1e-9
    @test abs(-(S[2][3])/S[1][2]-D.P1[2])<1e-8
end

@testset "Complex curves and local Abel integrals" begin
    for _ in 1:12
        R=4 .* (rand(ComplexF64,5) .- (.5+.5im))
        C=Genus2WCurve(R,source="Roots",e=1e-10,n=30)
        z=[.12+.02im,.05-.03im];D=C.JacobiInversion(z)
        @test maximum(curve_error(C,P) for P in (D.P1,D.P2))<1e-8
        # For complex curves only two periods are exposed: compare divisors,
        # since Abel coordinates are defined modulo the full rank-four lattice.
        DDback=C.JacobiInversion(C.AbelMap(D))
        dist(P,Q)=norm(P-Q)/max(1,norm(P),norm(Q))
        @test min(max(dist(D.P1,DDback.P1),dist(D.P2,DDback.P2)),
                  max(dist(D.P1,DDback.P2),dist(D.P2,DDback.P1)))<1e-7
        a=maximum(abs,R)+2;b=a+.5
        DD=Genus2WCurveDivisor((a,sqrt(C.F(a))),(b,-sqrt(C.F(b))))
        integral=[quadgk(x->x^j/sqrt(C.F(x)),a,b;atol=1e-13,rtol=0)[1] for j in (0,1)]
        @test norm(C.AbelMap(DD)+integral)<1e-8
    end
end

@testset "Experiment reproducibility and precision scope" begin
    include(joinpath(@__DIR__,"..","Numerical experiments","Experiments_genus2.jl"))
    original_precision=precision(BigFloat)
    res=testWithIntegration(2,2)
    @test length(res)==6
    @test maximum(res[5:6])<1e-9
    res=testWithIntegrationHighPrec(2,2,128;order=16)
    @test maximum(res[5:6])<big"1e-28"
    @test precision(BigFloat)==original_precision
    res=testWithMultiprecisionG2(1,2)
    @test length(res)==7 && all(length(r)==4 for r in res)
    @test all(isfinite,vcat(res...))
    @test res[4][end]<big"1e-105"
    @test res[7][end]<big"1e-105"
    @test precision(BigFloat)==original_precision
end

@testset "Close branch points at sufficient precision" begin
    setprecision(256) do
        for delta in (big"1e-6",big"1e-9")
            R=[big"0",delta,big"2",big"3",big"4"]
            C=Genus2WCurve(R,source="Roots",is_real=true,n=30)
            z=Complex{BigFloat}.([.1+.02im,.05-.03im])
            @test period_error(C,C.AbelMap(C.JacobiInversion(z))-z)<big"1e-65"
        end
    end
end

include("richelot_stability.jl")

@testset "Randomized Abel inversion regression" begin
    Random.seed!(2)
    original_precision=precision(BigFloat)
    result=testWithMultiprecisionG2(3,4;randomize=true,rng=MersenneTwister(2))
    @test all(isfinite,vcat(result...))
    @test all(result[7] .< [big"1e-20",big"1e-45",big"1e-90",big"1e-130"])
    @test precision(BigFloat)==original_precision
end
