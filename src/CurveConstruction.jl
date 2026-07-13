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
		F = undef
		Roots = undef
		if(source == "Coefficients")
			if(length(data) != 6)
				throw(ArgumentError("Invalid number of coefficients of the polynomial"))
			end
			if(abs(data[6] - 4.0) > eps(typeof(real(data[6]))))
				throw(ArgumentError("Leading coefficient is not equal to 4"))
			end
			data[6] = 4.0
			F = Polynomial(data .* 1.0 .+ 0.0im)
			Roots = PolynomialRoots.roots(coeffs(F))
			if(is_real == true && maximum(imag.(Roots)) > eps(typeof(real(Roots[1]))))
				throw(ArgumentError("The roots are not real"))
			end
		elseif(source == "Roots")
			if(length(data) != 5)
				throw(ArgumentError("Invalid number of roots"))
			end
			Roots = data .* 1.0 .+ 0.0im
			if(is_real == true && maximum(imag.(Roots)) > eps(typeof(real(Roots[1]))))
				throw(ArgumentError("The roots are not real"))
			end
			if(is_real == true)
				Roots = real.(Roots) .+ 0.0im
			end
			F = fromroots(Roots)*4
		elseif(source == "Polynomial")
			if(degree(data) != 5)
				throw(ArgumentError("Invalid degree of the polynomial"))
			end
			cfs = coeffs(data .* 1.0 .+ 0.0im)
			if(abs(cfs[6] - 4.0) > eps(typeof(real(cfs[6]))))
				throw(ArgumentError("Leading coefficient is not equal to 4"))
			end
			cfs[6] = 4.0
			F = Polynomial(cfs)
			Roots = PolynomialRoots.roots(coeffs(F))
			if(is_real == true && maximum(imag.(Roots)) > eps(typeof(real(Roots[1]))))
				throw(ArgumentError("The roots are not real"))
			end
		else 
			throw(ArgumentError("Invalid source: please ensure that source is equal to either \"Coefficients\", or \"Roots\", or \"Polynomial\";\nby default source is equal to \"Coefficients\""))
		end
		T = typeof(real(sum(Roots)))
		F = Polynomial(Complex{T}.(coeffs(F)))
		if(isnothing(e))
			e = 10*eps(T)
		end
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
				throw(InexactError("Cannot construct a valid triple of discs; try increasing precision"))
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
			S = KlFWeight2(z)
			if(abs(S[1][1]) < 100*e)
				if(abs(S[1][2]) < 100*e)
					return Genus2WCurveDivisor()
				end
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
		function Ab(DD::Genus2WCurveDivisor)
			if(length(DD.P2) > 0)
				a = GenerateRandomPoint(Roots)
				b = sqrt(F(a))
				D = [Genus2WCurveDivisor(DD.P1, (a,b)), Genus2WCurveDivisor(DD.P2, (a,-b))]
				v = [KummerCoord(D[j], F) for j in 1:2]
				v = [RichelotSequenceInversion(AMM1, AMM2, AMV, HatF, v[j]/norm(v[j]), e) for j in 1:2]
				x = [PolynomialRoots.roots([-v[j][3], -v[j][2], v[j][1]]) for j in 1:2]
				z = [DegenerateAbel(Lcoeff, p, x[j]) for j in 1:2]
				for j in 1:2
					S = KlFWeight2(z[j])
					diff2Pz = ((S[3] - S[1]*(S[3][1]/S[1][1]))/S[1][1])[2:3]
					diff2Ptrue = [(D[j].P2[2] - D[j].P1[2])/(D[j].P2[1] - D[j].P1[1]), (D[j].P2[1]*D[j].P1[2] - D[j].P1[1]*D[j].P2[2])/(D[j].P2[1] - D[j].P1[1])]
					if( norm(diff2Pz + diff2Ptrue) < norm(diff2Pz - diff2Ptrue) )
						z[j] = -z[j]
					end
				end
				if(is_real)
					return normalizeToPeriods(W, z[1] + z[2])
				end
				return z[1] + z[2]
			end
			if(length(DD.P1) > 0)
				a = GenerateRandomPoint(Roots)
				b = sqrt(F(a))
				D = [Genus2WCurveDivisor(DD.P1, (a,b)), Genus2WCurveDivisor((a,-b))]
				v = [KummerCoord(D[j], F) for j in 1:2]
				v = [RichelotSequenceInversion(AMM1, AMM2, AMV, HatF, v[j]/norm(v[j]), e) for j in 1:2]
				x = [PolynomialRoots.roots([-v[j][3], -v[j][2], v[j][1]]) for j in 1:2]
				z = [DegenerateAbel(Lcoeff, p, x[j]) for j in 1:2]
				S = KlFWeight2(z[1])
				diff2Pz = ((S[3] - S[1]*(S[3][1]/S[1][1]))/S[1][1])[2:3]
				diff2Ptrue = [(D[1].P2[2] - D[1].P1[2])/(D[1].P2[1] - D[1].P1[1]), (D[1].P2[1]*D[1].P1[2] - D[1].P1[1]*D[1].P2[2])/(D[1].P2[1] - D[1].P1[1])]
				if( norm(diff2Pz + diff2Ptrue) < norm(diff2Pz - diff2Ptrue) )
					z[1] = -z[1]
				end
				S = KlFWeight2(z[2])
				diff1Pz = [S[3][2]/S[1][2] - S[3][3]/S[1][3], S[2][2]/S[1][2] - S[2][3]/S[1][3]]
				diff1Ptrue = [0, b/a]
				println(diff1Pz, "\n", diff1Ptrue)
				if( norm(diff1Pz + diff1Ptrue) < norm(diff1Pz - diff1Ptrue) )
					z[2] = -z[2]
				end
				if(is_real)
					return normalizeToPeriods(W, z[1] + z[2])
				end
				return z[1] + z[2]
			end
			return zeros(Complex{T}, 2)
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
		push!(ArrA, copy(A))
		push!(ArrChi, copy(Chi))
		push!(AMM1, D)
		push!(AMM2, C)
		push!(AMV, [d, xp, xq, xr])
		push!(HatF, Lcoeff*fromroots(R))
		R = sortmethod(R)
		if(maximum([abs(R[j][1] - R[j][2]) for j in 1:3]) < e)
			break
		end
		if(minimum([abs(R[j][1] - R[j][2])^2 for j in 1:3]) < e)
			m = j
		end
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
