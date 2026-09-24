using KleinianGenus2
using LinearAlgebra
using Random
using QuadGK

# All timings exclude compilation of the functions being measured.
# Return order is retained for compatibility with the original scripts.
function integrationExperiment(m,n,T; rng=MersenneTwister(20260923), order=7,
                               randomize=false)
    m>0 && n>0 || throw(ArgumentError("m and n must be positive"))
    tol=100*eps(T)
    shifts=T.([-2,-1,0,1,2])./10
    roots=[sort(4 .* (rand(rng,T,5) .- T(1)/2)) .+ shifts for _ in 1:m]
    construct(r)=Genus2WCurve(r,source="Roots",e=tol,n=80,is_real=true)
    foreach(construct,roots)
    start=time_ns();curves=construct.(roots);tcon=time_ns()-start
    lb=[rand(rng,T,n) .+ T(5)/2 for _ in 1:m]
    rb=[lb[j] .+ 1 .+ rand(rng,T,n) for j in 1:m]
    f1=[x->inv(2sqrt(prod(x-r for r in roots[j]))) for j in 1:m]
    f2=[x->x*f1[j](x) for j in 1:m]
    q(f,a,b)=quadgk(f,a,b;atol=tol,rtol=zero(T),order=order)
    # Warm up at the SAME arithmetic precision and tolerance as the measurement.
    for j in 1:m
        q(f1[j],lb[j][1],rb[j][1]);q(f2[j],lb[j][1],rb[j][1])
    end
    start=time_ns();Q1=[[q(f1[j],lb[j][k],rb[j][k]) for k in 1:n] for j in 1:m];ti1=time_ns()-start
    start=time_ns();Q2=[[q(f2[j],lb[j][k],rb[j][k]) for k in 1:n] for j in 1:m];ti2=time_ns()-start
    all(v[2]<=tol for row in Q1 for v in row) || error("QuadGK did not meet atol for dx/y")
    all(v[2]<=tol for row in Q2 for v in row) || error("QuadGK did not meet atol for x dx/y")
    divisors=[[Genus2WCurveDivisor((lb[j][k],sqrt(curves[j].F(lb[j][k]))),
                 (rb[j][k],-sqrt(curves[j].F(rb[j][k])))) for k in 1:n] for j in 1:m]
    for j in 1:m
        curves[j].AbelMap(divisors[j][1],randomize)
    end
    start=time_ns();Z=[[curves[j].AbelMap(D,randomize) for D in divisors[j]] for j in 1:m];ta=time_ns()-start
    # A(P_a + conjugate(P_b) - 2infinity) = -integral_a^b, modulo periods.
    residual=[[PeriodReduction(Z[j][k]+[Q1[j][k][1],Q2[j][k][1]],curves[j].Periods) for k in 1:n] for j in 1:m]
    err1=sum(abs(v[1]) for row in residual for v in row)/(m*n)
    err2=sum(abs(v[2]) for row in residual for v in row)/(m*n)
    return (tcon*1e-6/m,ti1*1e-6/(m*n),ti2*1e-6/(m*n),ta*1e-6/(m*n),err1,err2)
end

function testWithIntegration(m,n;kwargs...)
    integrationExperiment(m,n,Float64;kwargs...)
end

function testWithIntegrationHighPrec(m,n,prec;kwargs...)
    setprecision(BigFloat,prec) do
        integrationExperiment(m,n,BigFloat;kwargs...)
    end
end

function testWithMultiprecisionG2(m,n; rng=MersenneTwister(20260923),randomize=false)
    m>0 && n>0 || throw(ArgumentError("m and n must be positive"))
    setprecision(BigFloat,512) do
        shifts=BigFloat.([-2,-1,0,1,2])./10
        roots=[sort(4 .* (rand(rng,BigFloat,5) .- BigFloat(1)/2)) .+ shifts for _ in 1:m]
        construct(r,tol)=Genus2WCurve(r,source="Roots",e=tol,n=80,is_real=true)
        reference=[construct(r,big"1e-140") for r in roots]
        periods=[C.Periods for C in reference]
        z=[[periods[j]*(rand(rng,BigFloat,4) .- BigFloat(1)/2) for _ in 1:n] for j in 1:m]
        Sref=[reference[j].KleinianFuncWeight2.(z[j]) for j in 1:m]
        Dref=[reference[j].JacobiInversion.(z[j]) for j in 1:m]
        tc=zeros(4);tf=zeros(4);ta=zeros(4)
        errors=[zeros(BigFloat,4) for _ in 1:4]
        tolerances=[big"1e-15",big"1e-30",big"1e-60",big"1e-120"]
        for (k,tol) in enumerate(tolerances)
            foreach(r->construct(r,tol),roots)
            start=time_ns()
            curves=[construct(r,tol) for r in roots]
            tc[k]=(time_ns()-start)*1e-6/m
            for j in 1:m
                curves[j].KleinianFunc(z[j][1]);curves[j].AbelMap(Dref[j][1],randomize)
            end
            start=time_ns()
            for j in 1:m, w in z[j];curves[j].KleinianFunc(w);end
            tf[k]=(time_ns()-start)*1e-6/(m*n)
            S=[curves[j].KleinianFuncWeight2.(z[j]) for j in 1:m]
            start=time_ns()
            A=[[curves[j].AbelMap(D,randomize) for D in Dref[j]] for j in 1:m]
            ta[k]=(time_ns()-start)*1e-6/(m*n)
            for q in 1:3
                errors[q][k]=sum(norm(S[j][l][q]-Sref[j][l][q]) for j in 1:m for l in 1:n)/(m*n)
            end
            errors[4][k]=sum(norm(PeriodReduction(z[j][l]-A[j][l],periods[j])) for j in 1:m for l in 1:n)/(m*n)
        end
        return (tc,tf,ta,errors...)
    end
end

# Keep the period reduction's arithmetic type; never demote BigFloat to Float64.
PeriodReduction(z,P)=KleinianGenus2.normalizeToPeriods(P,z)
function PeriodReductionAlt(z,P)
    M=KleinianGenus2.RealPeriodMatrix(P)
    c=M \ [real(z[1]),imag(z[1]),real(z[2]),imag(z[2])]
    return P*(c .- floor.(c))
end
