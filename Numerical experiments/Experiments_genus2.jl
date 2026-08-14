using KleinianGenus2
using LinearAlgebra
using Random
using QuadGK

function testWithIntegration(m,n)
	e = [4*(rand(Float64, 5) .- 0.5) for j in 1:m]
	e = sort!.(e)
	e = [e[j] .+ [-0.2, -0.1, 0, 0.1, 0.2] for j in 1:m]
	for j in 1:m
		Genus2WCurve(e[j], source = "Roots", e = 1e-16, n = 80, is_real = true)
	end
	
	t_constr = time_ns()
	C = [Genus2WCurve(e[j], source = "Roots", e = 1e-16, n = 80, is_real = true) for j in 1:m]
	t_constr = time_ns() - t_constr
	
	lb = [rand(Float64, n) .+ 2.5 for j in 1:m]
	rb = [(rand(Float64, n) .+ 1) + lb[j] for j in 1:m]
	
	f1 = [x -> 0.5 / sqrt(prod([x - e[j][k] for k in 1:5])) for j in 1:m]
	f2 = [x -> 0.5 * x / sqrt(prod([x - e[j][k] for k in 1:5])) for j in 1:m]
	for j in 1:m
		f1[j](3), f2[j](3)
		quadgk(f1[j], 3, 4, rtol = 1e-4), quadgk(f2[j], 3, 4, rtol = 1e-4)
	end
	
	t_int1 = time_ns()
	I1 = [[quadgk(f1[j], lb[j][k], rb[j][k], atol = 1e-15)[1] for k in 1:n] for j in 1:m]
	t_int1 = time_ns() - t_int1
	
	t_int2 = time_ns()
	I2 = [[quadgk(f2[j], lb[j][k], rb[j][k], atol = 1e-15)[1] for k in 1:n] for j in 1:m]
	t_int2 = time_ns() - t_int2
	
	DD = [[Genus2WCurveDivisor((lb[j][k], 2 * sqrt(prod([lb[j][k] - e[j][l] for l in 1:5]))), (rb[j][k], -2 * sqrt(prod([rb[j][k] - e[j][l] for l in 1:5]))))
					for k in 1:n] for j in 1:m]
	for j in 1:m
		for k in 1:n
			C[j].AbelMap(DD[j][k])
		end
	end
	
	t_ab = time_ns()
	Z = [[C[j].AbelMap(DD[j][k]) for k in 1:n] for j in 1:m]
	t_ab = time_ns() - t_ab
	
	#Z = [[PeriodReductionAlt(Z[j][k], C[j].Periods) for k in 1:n] for j in 1:m]
	
	
	I1A = [[Z[j][k][1] for k in 1:n] for j in 1:m]
	I2A = [[Z[j][k][2] for k in 1:n] for j in 1:m]
	
	err_1 = sum([sum(abs.(I1[j] .+ I1A[j]))/n for j in 1:m])/m
	err_2 = sum([sum(abs.(I2[j] .+ I2A[j]))/n for j in 1:m])/m
	return (t_constr * 1e-6 / m, t_int1 * 1e-6 / (n * m), t_int2* 1e-6 / (n * m), t_ab* 1e-6 / (n * m), err_1, err_2)
end

function testWithIntegrationHighPrec(m,n,prec)
	setprecision(prec)
	e = [4*(rand(BigFloat, 5) .- 0.5) for j in 1:m]
	e = sort!.(e)
	e = [e[j] .+ [-0.2, -0.1, 0, 0.1, 0.2] for j in 1:m]
	for j in 1:m
		Genus2WCurve(e[j], source = "Roots", e = 100*eps(BigFloat), n = 80, is_real = true)
	end
	
	t_constr = time_ns()
	C = [Genus2WCurve(e[j], source = "Roots", e = 100*eps(BigFloat), n = 80, is_real = true) for j in 1:m]
	t_constr = time_ns() - t_constr
	
	lb = [rand(BigFloat, n) .+ 2.5 for j in 1:m]
	rb = [(rand(BigFloat, n) .+ 1) + lb[j] for j in 1:m]
	
	f1 = [x -> 0.5 / sqrt(prod([x - e[j][k] for k in 1:5])) for j in 1:m]
	f2 = [x -> 0.5 * x / sqrt(prod([x - e[j][k] for k in 1:5])) for j in 1:m]
	for j in 1:m
		f1[j](3), f2[j](3)
		quadgk(f1[j], 3, 4, rtol = 1e-4), quadgk(f2[j], 3, 4, rtol = 1e-4)
	end
	
	t_int1 = time_ns()
	I1 = [[quadgk(f1[j], lb[j][k], rb[j][k], atol = 100*eps(BigFloat))[1] for k in 1:n] for j in 1:m]
	t_int1 = time_ns() - t_int1
	
	t_int2 = time_ns()
	I2 = [[quadgk(f2[j], lb[j][k], rb[j][k], atol = 100*eps(BigFloat))[1] for k in 1:n] for j in 1:m]
	t_int2 = time_ns() - t_int2
	
	DD = [[Genus2WCurveDivisor((lb[j][k], 2 * sqrt(prod([lb[j][k] - e[j][l] for l in 1:5]))), (rb[j][k], -2 * sqrt(prod([rb[j][k] - e[j][l] for l in 1:5]))))
					for k in 1:n] for j in 1:m]
	for j in 1:m
		for k in 1:n
			C[j].AbelMap(DD[j][k])
		end
	end
	
	t_ab = time_ns()
	Z = [[C[j].AbelMap(DD[j][k]) for k in 1:n] for j in 1:m]
	t_ab = time_ns() - t_ab
	
	#Z = [[PeriodReductionAlt(Z[j][k], C[j].Periods) for k in 1:n] for j in 1:m]
	
	
	I1A = [[Z[j][k][1] for k in 1:n] for j in 1:m]
	I2A = [[Z[j][k][2] for k in 1:n] for j in 1:m]
	
	err_1 = sum([sum(abs.(I1[j] .+ I1A[j]))/n for j in 1:m])/m
	err_2 = sum([sum(abs.(I2[j] .+ I2A[j]))/n for j in 1:m])/m
	return (t_constr * 1e-6 / m, t_int1 * 1e-6 / (n * m), t_int2* 1e-6 / (n * m), t_ab* 1e-6 / (n * m), err_1, err_2)
end

function testWithMultiprecisionG2(M,n)
	m = M
	setprecision(512)
	e = [4*(rand(BigFloat, 5) .- 0.5) for j in 1:m]
	e = sort!.(e)
	e = [e[j] .+ [-0.2, -0.1, 0, 0.1, 0.2] for j in 1:m]
	C = [Genus2WCurve(e[j], source = "Roots", e = 1e-140, n = 80, is_real = true) for j in 1:m]
	P = [C[j].Periods for j in 1:m]
	er = BigFloat(1e-15) #1e-15
	z = [[(rand(BigFloat) .- 0.5)*P[j][:,1] + (rand(BigFloat) .- 0.5)*P[j][:,2] + 
				(rand(BigFloat) .- 0.5)*P[j][:,3] + (rand(BigFloat) .- 0.5)*P[j][:,4] for k in 1:n] for j in 1:m]
	SS_max_prec = [C[j].KleinianFuncWeight2.(z[j]) for j in 1:m]
	S_max_prec = [[SS_max_prec[j][k][1] for k in 1:n] for j in 1:m]
	S1_max_prec = [[SS_max_prec[j][k][2] for k in 1:n] for j in 1:m]
	S2_max_prec = [[SS_max_prec[j][k][3] for k in 1:n] for j in 1:m]
	
	DD_max_prec = [C[j].JacobiInversion.(z[j]) for j in 1:m]
	
	t_constr = zeros(UInt64, 4)
	t_func = zeros(UInt64, 4)
	t_ab = zeros(UInt64, 4)
	
	err_s = zeros(BigFloat, 4)
	err_s1 = zeros(BigFloat, 4)
	err_s2 = zeros(BigFloat, 4)
	err_ab = zeros(BigFloat, 4)
	
	println("progress: high precision precomputed")
	
	for k = 1:4
		for j in 1:m
		Genus2WCurve(e[j], source = "Roots", e = er, n = 80, is_real = true)
		end
		
		t_constr[k] = time_ns()
		for p in 1:4
		Ctmp = [Genus2WCurve(e[j], source = "Roots", e = er, n = 80, is_real = true) for j in 1:m]
		end
		Ctmp = [Genus2WCurve(e[j], source = "Roots", e = er, n = 80, is_real = true) for j in 1:m]
		t_constr[k] = time_ns() - t_constr[k]
		
		for j in 1:m
			Ctmp[j].KleinianFunc(z[j][1])
			Ctmp[j].AbelMap(DD_max_prec[j][1])
		end
		
		t_func[k] = time_ns()
		for p in 1:4
		K = [Ctmp[j].KleinianFunc.(z[j]) for j in 1:m]
		end
		K = [Ctmp[j].KleinianFunc.(z[j]) for j in 1:m]
		t_func[k] = time_ns() - t_func[k]
		SS_tmp = [Ctmp[j].KleinianFuncWeight2.(z[j]) for j in 1:m]
		S_tmp = [[SS_tmp[j][k][1] for k in 1:n] for j in 1:m]
		S1_tmp = [[SS_tmp[j][k][2] for k in 1:n] for j in 1:m]
		S2_tmp = [[SS_tmp[j][k][3] for k in 1:n] for j in 1:m]
		
		println("progress: before Abel, prec = ", k)
		
		t_ab[k] = time_ns()
		abtmp = [Ctmp[j].AbelMap.(DD_max_prec[j]) for j in 1:m]
		t_ab[k] = time_ns() - t_ab[k]
		
		println("progress: after Abel, prec = ", k)
		
		err_s[k] = sum([sum(norm.(S_max_prec[j] .- S_tmp[j]))/n for j in 1:m])/m
		err_s1[k] = sum([sum(norm.(S1_max_prec[j] .- S1_tmp[j]))/n for j in 1:m])/m
		err_s2[k] = sum([sum(norm.(S2_max_prec[j] .- S2_tmp[j]))/n for j in 1:m])/m
		
		diff_ab = [[PeriodReduction(z[j][k] - abtmp[j][k], P[j]) for k in 1:n] for j in 1:m]
		err_ab[k] = sum([sum(norm.(diff_ab[j]))/n for j in 1:m])/m
			
		er = er^2
	end
	return (t_constr .* 1e-6 ./ (5 * m), t_func .* 1e-6 ./ (5 * m * n), t_ab .* 1e-6 ./ (m * n), err_s, err_s1, err_s2, err_ab)
end

function PeriodReduction(z, P)
	a = zeros(BigFloat, 4)
	a[1] = real(z[1])
	a[2] = imag(z[1])
	a[3] = real(z[2])
	a[4] = imag(z[2])
	M = zeros(BigFloat, 4, 4)
	for i in 0:1
		for j in 1:4
			M[2*i+1, j] = real(P[i+1 , j])
			M[2*i+2, j] = imag(P[i+1, j])
		end
	end
	c = inv(M) * a
	c = c - round.(c)
	return c[1] * P[:, 1] + c[2] * P[:, 2] + c[3] * P[:, 3] + c[4] * P[:, 4]
end

function PeriodReductionAlt(z, P)
	a = zeros(Float64, 4)
	a[1] = real(z[1])
	a[2] = imag(z[1])
	a[3] = real(z[2])
	a[4] = imag(z[2])
	M = zeros(Float64, 4, 4)
	for i in 0:1
		for j in 1:4
			M[2*i+1, j] = real(P[i+1 , j])
			M[2*i+2, j] = imag(P[i+1, j])
		end
	end
	c = inv(M) * a
	c = c - round.(c)
	for j in 1:4
		if(abs(c[j]) < 1e-6)
			c[j] = 0
		end
		if(c[j] < 0)
			c[j] = c[j] + 1
		end
	end
	println(c)
	return c[1] * P[:, 1] + c[2] * P[:, 2] + c[3] * P[:, 3] + c[4] * P[:, 4]
end
