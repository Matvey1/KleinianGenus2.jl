include("richelot_failure.jl")

@testset "Richelot branch retention and continuation" begin
    setprecision(512) do
        c=richelot_failure_case()
        image(w)=K.richelot_image(c.A,w)
        ws=K.RichelotInversion(c.data,c.M,c.F,c.V,c.e)
        @test length(ws)==8
        @test minimum(K.ProjectiveDistance(image(w),c.v) for w in ws)<big"1e-145"
        w=K.richelot_newton(c.A,c.v,c.v)
        @test w !== nothing
        if w !== nothing
            @test K.ProjectiveDistance(image(w),c.v)<big"1e-145"
            @test K.kummer_residual(w,c.F)<big"1e-145"
            @test K.ProjectiveDistance(w,c.v)<big"1e-30"
        end
    end
end

# Independent formula for the canonical quadratic map (before either linear
# coordinate change). Exercise 0, 1, 2 and 3 small/vanishing discriminants.
function canonical_richelot(d,x,u)
    a=x.^2
    [2d^2*u[1]*u[4]+2d^4*a[1]*u[2]*u[3],
     2d^2*u[2]*u[4]+2d^4*a[2]*u[1]*u[3],
     2d^2*u[3]*u[4]+2d^4*a[3]*u[1]*u[2],
     u[4]^2+d^4*(a[2]*a[3]*u[1]^2+a[1]*a[3]*u[2]^2+a[1]*a[2]*u[3]^2)]
end
@testset "Canonical formulas with small discriminants" begin
    setprecision(256) do
        d=complex(big"1.2",big"0.1")
        u=Complex{BigFloat}.([1,2,3,4]) .+ im/10
        M=Matrix{Complex{BigFloat}}(I,4,4)
        for small in (big"0",big"1e-40"), nsmall in 0:3
            x=Complex{BigFloat}.([big"0.2",big"0.3",big"0.4"])
            x[1:nsmall].=small
            target=canonical_richelot(d,x,u)
            ws=K.RichelotInversion([d;x],M,Polynomial([one(d)]),target,big"1e-60")
            @test minimum(K.ProjectiveDistance(w,u) for w in ws)<big"1e-65"
            @test all(K.ProjectiveDistance(canonical_richelot(d,x,w),target)<big"1e-65" for w in ws)
        end
    end
end

@testset "Projective Newton charts and scaling" begin
    setprecision(256) do
        # Q(w)=w1*w is the projective identity on w1 != 0. The solution has
        # w4=0, where a fixed fourth-coordinate chart cannot be used.
        A=zeros(Complex{BigFloat},4,4,4)
        for k in 1:4
            A[k,1,k]+=big"0.5"; A[k,k,1]+=big"0.5"
        end
        target=Complex{BigFloat}.([1,big"0.2",big"0.3",0])
        for factor in (big"1e-200",big"1",big"1e200")
            w=K.richelot_newton(factor*A,(2+im)*target,(3-im)*(target .+ big"0.01"))
            @test w !== nothing
            if w !== nothing
                @test K.ProjectiveDistance(w,target)<big"1e-70"
            end
        end
        @test K.richelot_newton(zero(A),target,target) === nothing
    end
end
