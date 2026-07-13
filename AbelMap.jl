function DegenerateAbel(lcoeff, p, x)
	q = sqrt(lcoeff)
	l1 = log((x[2] - p[1])/(x[1] - p[1]))/((p[1] - p[2])*(p[1] - p[3]))
	l2 = log((x[2] - p[2])/(x[1] - p[2]))/((p[2] - p[1])*(p[2] - p[3]))
	l3 = log((x[2] - p[3])/(x[1] - p[3]))/((p[3] - p[1])*(p[3] - p[2]))
	return [(l1 + l2 + l3)/q, (p[1]*l1 + p[2]*l2 + p[3]*l3)/q]
end

function RichelotSequenceInversion(AMM1, AMM2, AMV, HatF, v, e)
	for i in 1:length(AMM1)
		V = copy(v)
		V = AMM1[i]*V
		Varr = RichelotInversion(AMV[i], AMM2[i], HatF[i], V, e)
		j = findmin([truedist(Varr[j], v) for j in 1:length(Varr)])[2]
		v = Varr[j]
	end
	return v
end

function RichelotInversion(DataVec, M, F, v, e)
	d = DataVec[1]
	x = DataVec[2:4]
	discr = x.^2
	ind = (Int64).(abs.(discr) .> 10*e)
	J = z->(   [2*d^2*z[4]						2*d^4*discr[1]*z[3]				2*d^4*discr[1]*z[2]				2*d^2*z[1];
				2*d^4*discr[2]*z[3]				2*d^2*z[4]						2*d^4*discr[2]*z[1]				2*d^2*z[2];
				2*d^4*discr[3]*z[2]				2*d^4*discr[3]*z[1]				2*d^2*z[4]						2*d^2*z[3];
				2*d^4*discr[2]*discr[3]*z[1]	2*d^4*discr[1]*discr[3]*z[2]	2*d^4*discr[1]*discr[2]*z[3]	2*z[4]])		
	IM = y->([2*d^4*discr[1]*y[2]*y[3] + 2*d^2*y[1]*y[4], 
		      2*d^4*discr[2]*y[1]*y[3] + 2*d^2*y[2]*y[4],
			  2*d^4*discr[3]*y[1]*y[2] + 2*d^2*y[3]*y[4],
			  d^4*(discr[2]*discr[3]*y[1]^2 + discr[1]*discr[3]*y[2]^2 + discr[1]*discr[2]*y[3]^2) + y[4]^2] .- v)
		
	function NewtonCorrection(w)
		u = copy(w)
		i = 0
		while(norm(IM(u)) > 10*e || i < 30)
			u = u - inv(J(u))*IM(u)
			i = i+1
		end
		return u
	end
	
	if sum(ind) == 3
		setprecision(BigFloat, 2*precision(real(d))) do
			CD = Complex{BigFloat}
			xd = CD.(x)
			vd = CD.(v)
			dd = CD(d)

			T = [xd[2]*xd[3], xd[1]*xd[3], xd[1]*xd[2], one(CD)]

			a = T[1]*vd[1]; b = T[2]*vd[2]; c = T[3]*vd[3]; e = T[4]*vd[4]
			u = a + b;  p = a - b
			q = c + e;  r = c - e

			qv = [ u + q,
				   q - u,
				  -(p + r),
				   p - r ]

			w = sqrt.(qv)

			all_signs = ((1,1,1,1), (-1,1,1,1), (1,-1,1,1), (1,1,-1,1),
						 (-1,-1,1,1), (-1,1,-1,1), (1,-1,-1,1), (-1,-1,-1,1))
			Warr = [w .* collect(s) for s in all_signs]

			pr = prod(sort(xd, by = abs)) * dd^2

			P = [xd[1] -xd[1] -xd[1] xd[1];
				 xd[2] -xd[2] xd[2] -xd[2];
				 xd[3] xd[3] -xd[3] -xd[3];
				 pr pr pr pr]

			Warr = [P*Warr[j] for j in 1:8]
		end 
		WarrT = [Complex{typeof(real(d))}.(w) for w in Warr]

		WarrT = sort(WarrT, by = w -> abs(KummerEquationCheck(M*w/norm(M*w), F)))[1:4]
		WarrT = [M*WarrT[j] for j in 1:4]
		WarrT = [WarrT[j]/norm(WarrT[j]) for j in 1:4]
		return WarrT
	end
	
	
	Warr = undef
	if(sum(ind) == 2)
		j = findmin(ind)[2]
		j1 = mod(j,3) + 1
		j2 = mod(j1,3) + 1
		t = x[j1]*x[j2]
		o = [(v[4] + t*v[j])^(0.5), (v[4] - t*v[j])^(0.5)]
		O = inv([1 t*d^2; 1 -t*d^2])
		Warrtmp = [O*o, O*(o.*[-1,1])]
		Warrtmp2 = [ inv([2*d^2*Warrtmp[j][1] 2*d^4*Warrtmp[j][2]*discr[j1]; 2*d^4*Warrtmp[j][2]*discr[j2] 2*d^2*Warrtmp[j][1]]) * [v[j1], v[j2]] for j in 1:2]
		Warr = [zeros(typeof(v[1]), 4), zeros(typeof(v[1]), 4)]
		for k in 1:2
			Warr[k][j] = Warrtmp[k][2]
			Warr[k][4] = Warrtmp[k][1]
			Warr[k][j1] = Warrtmp2[k][1]
			Warr[k][j2] = Warrtmp2[k][2]
		end
	end
	if(sum(ind) == 1)
		j = findmax(ind)[2]
		j1 = mod(j,3) + 1
		j2 = mod(j1,3) + 1
		w = [0,0,0,sqrt(v[4])]
		w[j1] = v[j1]/(2*d^2*w[4])
		w[j2] = v[j2]/(2*d^2*w[4])
		w[j] = (v[j] - 2*d^4*w[j1]*w[j2]*discr[j])/(2*d^2*w[4])
		Warr = [w]
	end
	if(sum(ind) == 0)
		w = [0,0,0,sqrt(v[4])]
		for j in 1:3
			w[j] = v[j]/(2*d^2*w[4])
		end
		Warr = [w]
	end
	res = NewtonCorrection.(Warr)
	Warr = [M*res[j] for j in 1:length(res)]
	Warr = [Warr[j]/norm(Warr[j]) for j in 1:length(res)]
	return Warr
end
