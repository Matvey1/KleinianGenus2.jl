struct Genus2WCurve{T<:AbstractFloat,KF<:Function,SKF<:Function,KF2<:Function,JI<:Function,AM<:Function}
	F::Polynomial{Complex{T}}
	Roots::Array{Complex{T}}
	Periods::Matrix{Complex{T}}
	Eta::Matrix{Complex{T}}
	KleinianFunc::KF
	KleinianFuncSigmaSquared::SKF
	KleinianFuncWeight2::KF2
	JacobiInversion::JI
	AbelMap::AM
	
	function Genus2WCurve(data;source = "Coefficients", n = 15, e = nothing, is_real = false, triple_of_discs = nothing)
        n isa Integer && n > 0 || throw(ArgumentError("n must be a positive integer"))
        is_real isa Bool || throw(ArgumentError("is_real must be Boolean"))
        source in ("Coefficients", "Roots", "Polynomial") || throw(ArgumentError("Invalid source"))
        raw = source == "Polynomial" ? coeffs(data) : collect(data)
        length(raw) == (source == "Roots" ? 5 : 6) || throw(ArgumentError("Invalid number of roots or coefficients"))
        values = complex.(float.(raw))
        T = promote_type(Float64, typeof(real(sum(values))))
        all(isfinite, values) || throw(ArgumentError("Roots and coefficients must be finite"))
        if source == "Roots"
            Roots = Complex{T}.(values)
            F = 4*fromroots(Roots)
        else
            abs(values[6] - 4) <= eps(T) || throw(ArgumentError("Leading coefficient is not equal to 4"))
            values[6] = 4
            F = Polynomial(Complex{T}.(values))
            Roots = PolynomialRoots.roots(coeffs(F))
        end
        if is_real
            maximum(abs, imag.(Roots)) <= 100*eps(T)*max(one(T), maximum(abs, Roots)) ||
                throw(ArgumentError("The roots are not real"))
            Roots = complex.(real.(Roots))
            F = 4*fromroots(Roots)
        end
        all(isfinite, Roots) || throw(ArgumentError("Root calculation failed"))
        length(unique(Roots)) == 5 || throw(ArgumentError("The curve must have five distinct roots"))
        e = isnothing(e) ? 10*eps(T) : T(e)
        isfinite(e) && e > 0 || throw(ArgumentError("e must be finite and positive"))
		ArrA = undef
		ArrChi = undef
		Lcoeff = undef
		p = undef
		W = undef
		E = undef
		AA = undef
		ChiChi = undef
		m = undef
		AMM1 = undef
		AMM2 = undef
		AMV = undef
		HatF = undef
		if(is_real)
			ArrAleft,ArrChileft,Lcoeffleft,pleft,AAleft,ChiChileft,mleft,AMM1left,AMM2left,AMVleft,HatFleft = RichelotPrecomputation(Roots,n,e,LeftRealSortMethod,T)
			ArrAright,ArrChiright,Lcoeffright,pright,AAright,ChiChiright,mright,AMM1right,AMM2right,AMVright,HatFright = RichelotPrecomputation(Roots,n,e,RightRealSortMethod,T)
			Wleft,Eleft = CalculatePeriods(ArrChileft,Lcoeffleft,pleft)
			Wright,Eright = CalculatePeriods(ArrChiright,Lcoeffright,pright,mode="right")
			W = hcat(Wleft,Wright)
			E = hcat(Eleft,Eright)
			if(length(ArrAleft) < length(ArrAright) || (length(ArrAleft) == length(ArrAright) && mleft >= mright))
				ArrA,ArrChi,Lcoeff,p,AA,ChiChi,m,AMM1,AMM2,AMV,HatF = ArrAleft,ArrChileft,Lcoeffleft,pleft,AAleft,ChiChileft,mleft,AMM1left,AMM2left,AMVleft,HatFleft
			else
				ArrA,ArrChi,Lcoeff,p,AA,ChiChi,m,AMM1,AMM2,AMV,HatF = ArrAright,ArrChiright,Lcoeffright,pright,AAright,ChiChiright,mright,AMM1right,AMM2right,AMVright,HatFright
			end
		else
			isn = isnothing(triple_of_discs)
			if(isn)
				triple_of_discs = generateTripleOfDiscs(Roots)
			end
			ind = verifySortingMethod(Roots, triple_of_discs)
			if(isn && ind != 0)
				throw(ErrorException("Cannot construct a valid triple of discs; try increasing precision or specifying triple_of_discs"))
			end
			if(!isn && ind == 1)
				throw(ArgumentError("The specified triple of discs is not disjoint"))
			end
			if(!isn && ind == 2)
				throw(ArgumentError("The specified triple of discs does not split the Weierstrass points into three pairs"))
			end
			sorting_method = GenerateSortingMethod(triple_of_discs)
			ArrA,ArrChi,Lcoeff,p,AA,ChiChi,m,AMM1,AMM2,AMV,HatF = RichelotPrecomputation(Roots,n,e,sorting_method,T)
			W,E = CalculatePeriods(ArrChi,Lcoeff,p)
		end
		function KlFWeight2(w)
			S = SFuncGeneric(ArrA,ArrChi,Lcoeff,p,w)
			return S
		end
		function KlSigmaSq(w)
			return GenericKleinianSigmaSquared(SFuncGeneric(ArrA,ArrChi,Lcoeff,p,w))
		end
		function Kl(w)
			return GenericKleinianDuplication(ArrA,ArrChi,Lcoeff,p,AA,ChiChi,w)
		end
		function JacInv(z)
            if is_real
                original_scale = max(one(T),norm(z))
                z = reduce_periods(z)
                norm(z) <= 8*eps(T)*original_scale && return Genus2WCurveDivisor()
            end
            all(iszero,z) && return Genus2WCurveDivisor()
            S = KlFWeight2(z)
            # Compare projective coordinates, not their absolute magnitudes.
            if abs(S[1][1]) <= 100*e*max(abs(S[1][2]),abs(S[1][3]))
                iszero(S[1][2]) && throw(ErrorException("Jacobi inversion is ill-conditioned; increase precision"))
				x = -S[1][3]/S[1][2]
				y = -(S[2][2]/S[1][2])*x  - S[2][3]/S[1][2]#(S[2][2] * S[1][3] - S[2][3]*S[1][2]) / (S[1][2])^2
				return Genus2WCurveDivisor((x,y))
			end
			p22 = S[1][2]/S[1][1]
			p12 = S[1][3]/S[1][1]
			xvec = PolynomialRoots.roots([-p12, -p22, 1])
			p222 = (S[3][2] - S[3][1]*p22)/S[1][1]
			p221 = (S[3][3] - S[3][1]*p12)/S[1][1]
			yvec = xvec.*p222 .+ p221
			return Genus2WCurveDivisor((xvec[1], yvec[1]), (xvec[2], yvec[2]))
		end
        # Factor the real period matrix once, not at every Abel evaluation.
        period_factor = is_real ? factorize(RealPeriodMatrix(W)) : nothing
        reduce_periods = z -> isnothing(period_factor) ? z :
            z - W*round.(period_factor \ [real(z[1]), imag(z[1]), real(z[2]), imag(z[2])])

        function abel_generic(D)
            target = KummerCoord(D, F)
            v = target
            v = RichelotSequenceInversion(ArrA, AMM1, AMM2, AMV, HatF, v/norm(v), e)
            x = PolynomialRoots.roots([-v[3], -v[2], v[1]])
            length(x) == 2 || throw(ErrorException("Degenerate Abel inversion; use randomize=true"))
            z = DegenerateAbel(Lcoeff, p, x)
            S = KlFWeight2(z)
            ProjectiveDistance(S[1],target) <= max(100*e,1000*eps(T)) ||
                throw(ErrorException("Abel inversion failed the Kummer-coordinate check; increase precision"))
            if !isempty(D.P2)
                slope = (D.P2[2]-D.P1[2])/(D.P2[1]-D.P1[1])
                intercept = D.P1[2]-slope*D.P1[1]
                expected = [slope, intercept]
                actual = ((S[3] - S[1]*(S[3][1]/S[1][1]))/S[1][1])[2:3]
            else
                # On the sigma divisor y = -(x*d_1 S22 + d_1 S12)/S22.
                # This expression also works at x=0 (no division by x or S12).
                expected = D.P1[2]
                actual = -(D.P1[1]*S[2][2]+S[2][3])/S[1][2]
            end
            all(isfinite, z) && all(isfinite, actual) ||
                throw(ErrorException("Non-finite Abel inversion; use randomize=true or increase precision"))
            return norm(actual+expected) < norm(actual-expected) ? -z : z
        end
        function Ab(DD::Genus2WCurveDivisor, randomize = false)
            for P in (DD.P1, DD.P2)
                isempty(P) && continue
                all(isfinite, P) || throw(ArgumentError("Divisor coordinates must be finite"))
                scale = max(one(T), abs(P[2])^2, sum(abs(F[k])*abs(P[1])^k for k in 0:5))
                abs(P[2]^2-F(P[1])) <= 1000*max(e,eps(T))*scale ||
                    throw(ArgumentError("Divisor point is not on the curve"))
            end
            isempty(DD.P1) && return zeros(Complex{T}, 2)
            if !isempty(DD.P2) && DD.P1[1] == DD.P2[1]
                DD.P1[2] == -DD.P2[2] && return zeros(Complex{T}, 2)
                randomize = true # repeated point: avoid the secant's 0/0
            end
            if randomize
                aux = isempty(DD.P2) ? [DD.P1[1]] : [DD.P1[1],DD.P2[1]]
                for attempt in 1:8
                    a = GenerateRandomPoint(Roots, aux)
                    b = sqrt(F(a))
                    first = Genus2WCurveDivisor(DD.P1, (a,b))
                    second = isempty(DD.P2) ? Genus2WCurveDivisor((a,-b)) : Genus2WCurveDivisor(DD.P2,(a,-b))
                    try
                        return reduce_periods(abel_generic(first)+abel_generic(second))
                    catch err
                        if !(err isa ErrorException || err isa LinearAlgebra.SingularException) || attempt == 8
                            rethrow()
                        end
                    end
                end
            end
            z = try
                abel_generic(DD)
            catch err
                if err isa ErrorException || err isa LinearAlgebra.SingularException
                    return Ab(DD, true)
                end
                rethrow()
            end
            if is_real && isempty(DD.P2) && iszero(DD.P1[2])
                # A branch point is exactly a point of order two in the Jacobian.
                coordinates = period_factor \ [real(z[1]),imag(z[1]),real(z[2]),imag(z[2])]
                z = W*(round.(2 .* coordinates)./2)
            end
            return reduce_periods(z)
        end
		return new{T, typeof(Kl),typeof(KlSigmaSq),typeof(KlFWeight2),typeof(JacInv),typeof(Ab)}(F, Roots, W, E, Kl, KlSigmaSq, KlFWeight2, JacInv, Ab)	
	end
	
end

Base.show(io::IO,C::Genus2WCurve) = print(io, Genus2WCurveString(C))

function Genus2WCurveString(C::Genus2WCurve)
	s = "2WCurve{" * string(typeof(real(C.Roots[1]))) * "}\n" * "equation: y^2 = " * string(C.F) * "\n" *
			"roots = " * string(C.Roots[1]) * ", " * string(C.Roots[2]) * ", " * string(C.Roots[3])* ", " * string(C.Roots[4]) * ", " * string(C.Roots[5]) * "\n"
	s = s * "Periods = "
	for i in 1:size(C.Periods)[2]
		s = s * string(C.Periods[:,i]) 
		if(i < size(C.Periods)[2])
			s = s * ", "
		end
	end
	s = s * "\n"
	s = s * "Eta-periods = "
	for i in 1:size(C.Eta)[2]
		s = s * string(C.Eta[:,i]) 
		if(i < size(C.Eta)[2])
			s = s * ", "
		end
	end
	s = s * "\nFunctions: KleinianFunc, KleinianFuncSigmaSquared, KleinianFuncWeight2, JacobiInversion, AbelMap"
	return s
end

function RichelotPrecomputation(Roots, n, e, sortmethod, T)
	ArrA = Array{Array{Complex{T}, 3}}([])
	ArrChi = Array{Matrix{Complex{T}}}([])
	Lcoeff = Complex{T}(4)
	AA = undef
	ChiChi = undef
	AMM1 = Array{Matrix{Complex{T}}}([])
	AMM2 = Array{Matrix{Complex{T}}}([])
	AMV = Array{Array{Complex{T}}}([])
	HatF = Array{Polynomial{Complex{T}}}([])
	R = sortmethod(Roots)
	m = 1
	for j in 1:n
		p = Lcoeff*fromroots(R[1])
		q = fromroots(R[2])
		r = fromroots(R[3])
		d = DDD(p,q,r)
		xp = Lcoeff * AuxDiff(R[1])
		xq = AuxDiff(R[2])
		xr = AuxDiff(R[3])
		A,Chi,Lcoeff,D,C,R = RichelotData(p,q,r,T)
		if(j == 1)
			RR = [R[1:2], R[3:4], R[5:length(R)]]
			pp = Lcoeff*fromroots(RR[1])
			qq = fromroots(RR[2])
			rr = fromroots(RR[3])
			AA,ChiChi, LLcoeffcoeff, DD, CC, RR = RichelotData(pp,qq,rr,T)
		end
		R = sortmethod(R)
		push!(ArrA, copy(A))
		push!(ArrChi, copy(Chi))
		push!(AMM1, D)
		push!(AMM2, C)
		push!(AMV, [d, xp, xq, xr])
		push!(HatF, Lcoeff*fromroots(vcat(R[1], R[2], R[3])))
        if minimum(abs(R[k][1]-R[k][2])^2 for k in 1:3) < e
            m = j
        end
        if maximum(abs(R[k][1]-R[k][2]) for k in 1:3) < e
            break
        end
        j == n && @warn "Richelot iteration did not reach e; increase n or working precision" n e
	end
	p = SplitToLimitRoots(R)
	return ArrA, ArrChi, Lcoeff, p, AA, ChiChi, m, AMM1, AMM2, AMV,HatF
end

function RichelotData(p,q,r,T)
	(pp,qq,rr) = RichelotTriple(p,q,r)
	R = (a,b,c,d,e)->Rexpr(p,q,r,pp,qq,rr,a,b,c,d,e)
	d = DDD(p,q,r)
	D,C = AbelMapMatrices(d,p,q,r,pp,qq,rr)
	cc = 1/(4*d)
	Dp = Dis(p)
	Dq = Dis(q)
	Dr = Dis(r)
	Vp = [-p[0], -p[1]/2, p[2], 0]
	Vq = [-q[0], -q[1]/2, q[2], 0]
	Vr = [-r[0], -r[1]/2, r[2], 0]
	Vinf = [4*R(0,0,0,2,2) + R(2,1,1,0,0) + R(2,0,0,1,1), -(d*p[1]*q[1]*r[1]/2 + R(1,0,2,1,1)), 
				 -(4*R(2,2,2,0,0) + R(0,1,1,2,2) + R(0,2,2,1,1)), -d*16]
	A = zeros(Complex{T}, 4, 4, 4)
	for a in 1:3
		A[a,1,1] = R(3-a,1,1,0,0) - 4*R(3-a,0,2,0,0)
		A[a,1,2] = (R(3-a,1,1,0,1) - 4*R(3-a,0,2,0,1))/2
		A[a,1,3] = -R(3-a,1,1,0,2) + 4*R(3-a,0,2,0,2)
		A[a,2,2] = (R(3-a,1,1,1,1) - 4*R(3-a,0,2,1,1))/4
		A[a,2,3] = -(R(3-a,1,1,1,2) - 4*R(3-a,0,2,1,2))/2
		A[a,3,3] = R(3-a,1,1,2,2) - 4*R(3-a,0,2,2,2)
		A[a,2,1] = A[a,1,2]
		A[a,3,1] = A[a,1,3]
		A[a,3,2] = A[a,2,3]
		if(a > 1)
			A[a,:,:] = -A[a,:,:]
		end
		A[a,a,:] = A[a,a,:] .+ Vinf./2
		A[a,:,a] = A[a,:,a] .+ Vinf./2
	end
	A = A*(2*d^2)
	A[4,:,:] = (-(p[0]*q[1]*r[1] + p[1]*q[0]*r[1] + p[1]*q[1]*r[0])*A[1,:,:] .- 0.5*p[1]*q[1]*r[1]*A[2,:,:] 
				.+ (p[2]*q[1]*r[1] + p[1]*q[2]*r[1] + p[1]*q[1]*r[2])*A[3,:,:] 
				-d*Vinf*transpose(Vinf) .- d^3*(Dq*Dr*Vp*transpose(Vp) + Dp*Dr*Vq*transpose(Vq) + Dp*Dq * Vr*transpose(Vr)))/8
	Chi = [ -A[4,1,4]   -A[4,2,4];
			-A[4,2,4]   A[4,3,4] ] / (A[4,4,4]/2)
	return A, Chi, cc*pp[2]*qq[2]*rr[2], D,C, vcat(PolynomialRoots.roots(coeffs(pp)), PolynomialRoots.roots(coeffs(qq)), PolynomialRoots.roots(coeffs(rr)))
end

function AbelMapMatrices(d,p,q,r,pp,qq,rr)
	D = [pp[2] qq[2] rr[2] 0;
		 -pp[1] -qq[1] -rr[1] 0;
		 -pp[0] -qq[0] -rr[0] 0;
		 (d - p[1]*pp[1])*q[1]*r[1]/8 (d - q[1]*qq[1])*p[1]*r[1]/8 (d - r[1]*rr[1])*p[1]*q[1]/8 -d/8]
	C = hcat(HalfPeriodCoords(1/(4*d), pp, qq*rr), HalfPeriodCoords(1/(4*d), qq, pp*rr), HalfPeriodCoords(1/(4*d), rr, pp*qq), [0,0,0,-1/(16*d)])
	return inv(D), C
end
