using WeierstrassElliptic
using EllipticFunctions
using Random

function compareWithEF(m, n)
    g2 = rand(Complex{Float64}, m)
    g3 = rand(Complex{Float64}, m)
	t_constr = time_ns()
    C = [WCurve((g2[j],g3[j]), source = "Invariants", e = 1e-16) for j in 1:m]
	t_constr = time_ns() - t_constr
    w1 = [C[j].Periods[1] for j in 1:m]
    w2 = [C[j].Periods[2] for j in 1:m]
	
	#to exclude compile time
    wsigma(0.1, g = (g2[1],g3[1]))
    wzeta(0.1, g = (g2[1],g3[1]))
    wp(0.1, g = (g2[1],g3[1]))
	wp(0.1, g = (g2[1],g3[1]), derivative = 1)
    time_ns()
	for j in 1:m
		C[j].WeierstrassFunc(0.1)
	end
	#################
	
    z = [(rand(Float64, n) .- 0.5).*w1[j] + (rand(Float64, n) .- 0.5).*w2[j] for j in 1:m]
    t_my = time_ns()
    w = [C[j].WeierstrassFunc.(z[j]) for j in 1:m]
    t_my = time_ns() - t_my
    t_sigma_EF = time_ns()
    sig = [wsigma.(z[j], g = (g2[j],g3[j])) for j in 1:m]
    t_sigma_EF = time_ns() - t_sigma_EF
    t_zeta_EF = time_ns()
    zet = [wzeta.(z[j], g = (g2[j],g3[j])) for j in 1:m]
    t_zeta_EF = time_ns() - t_zeta_EF
    t_p_EF = time_ns()
    p = [wp.(z[j], g = (g2[j],g3[j])) for j in 1:m]
    t_p_EF = time_ns() - t_p_EF
    t_pz_EF = time_ns()
    pz = [wp.(z[j], g = (g2[j],g3[j]), derivative = 1) for j in 1:m]
    t_pz_EF = time_ns() - t_pz_EF
    my_sigma = [[w[j][k][1] for k in 1:n] for j in 1:m]
    my_zeta = [[w[j][k][2] for k in 1:n] for j in 1:m]
    my_p = [[w[j][k][3] for k in 1:n] for j in 1:m]
    my_pz = [[w[j][k][4] for k in 1:n] for j in 1:m]
    diff_sigma = sum([sum(abs.(my_sigma[j] .- sig[j]))/n for j in 1:m])/m
    diff_zeta = sum([sum(abs.(my_zeta[j] .- zet[j]))/n for j in 1:m])/m
    diff_p = sum([sum(abs.(my_p[j] .- p[j]))/n for j in 1:m])/m
    diff_pz = sum([sum(abs.(my_pz[j] .- pz[j]))/n for j in 1:m])/m
    println("My time construction: ", t_constr * 1e-6 / m)
	println("My time func: ", t_my * 1e-6 / (m * n))
    println("EF time sigma: ", t_sigma_EF * 1e-6 / (m * n))
    println("EF time zeta: ", t_zeta_EF * 1e-6 / (m * n))
    println("EF time p: ", t_p_EF * 1e-6 / (m * n))
    println("EF time pz: ", t_pz_EF * 1e-6 / (m * n))
    println("avg errors:")
    println("sigma: ", diff_sigma)
    println("zeta: ", diff_zeta)
    println("p: ", diff_p)
    println("pz: ", diff_pz)
end

function testWithMultiprecision(m,n)
	setprecision(1024)
	g2 = rand(Complex{BigFloat}, m)
    g3 = rand(Complex{BigFloat}, m)
	C = [WCurve((g2[j],g3[j]), source = "Invariants", e = 1e-280, n = 50) for j in 1:m]
	w1 = [C[j].Periods[1] for j in 1:m]
    w2 = [C[j].Periods[2] for j in 1:m]
	er = 1/BigFloat(100000)^2 #1e-10
	z = [(rand(BigFloat, n) .- 0.5).*w1[j] + (rand(BigFloat, n) .- 0.5).*w2[j] for j in 1:m]
	w_max_prec = [C[j].WeierstrassFunc.(z[j]) for j in 1:m]
	s_max_prec = [[w_max_prec[j][k][1] for k in 1:n] for j in 1:m]
	zeta_max_prec = [[w_max_prec[j][k][2] for k in 1:n] for j in 1:m]
	p_max_prec = [[w_max_prec[j][k][3] for k in 1:n] for j in 1:m]
	pz_max_prec = [[w_max_prec[j][k][4] for k in 1:n] for j in 1:m]
	
	t_constr = zeros(UInt64, 5)
	t_func = zeros(UInt64, 5)
	t_ab = zeros(UInt64, 5)
	
	err_s = zeros(BigFloat, 5)
	err_zeta = zeros(BigFloat, 5)
	err_p = zeros(BigFloat, 5)
	err_pz = zeros(BigFloat, 5)
	err_ab = zeros(BigFloat, 5)
	
	for k = 1:5
		t_constr[k] = time_ns()
		Ctmp = [WCurve((g2[j],g3[j]), source = "Invariants", e = er, n = 50) for j in 1:m]
		t_constr[k] = time_ns() - t_constr[k]
		t_func[k] = time_ns()
		wtmp = [Ctmp[j].WeierstrassFunc.(z[j]) for j in 1:m]
		t_func[k] = time_ns() - t_func[k]
		t_ab[k] = time_ns()
		abtmp = [[Ctmp[j].AbelMap(p_max_prec[j][k], pz_max_prec[j][k]) for k in 1:n] for j in 1:m]
		t_ab[k] = time_ns() - t_ab[k]
		s_tmp = [[wtmp[j][k][1] for k in 1:n] for j in 1:m]
		zeta_tmp = [[wtmp[j][k][2] for k in 1:n] for j in 1:m]
		p_tmp = [[wtmp[j][k][3] for k in 1:n] for j in 1:m]
		pz_tmp = [[wtmp[j][k][4] for k in 1:n] for j in 1:m]
		
		err_s[k] = sum([sum(abs.(s_max_prec[j] .- s_tmp[j]))/n for j in 1:m])/m
		err_zeta[k] = sum([sum(abs.(zeta_max_prec[j] .- zeta_tmp[j]))/n for j in 1:m])/m
		err_p[k] = sum([sum(abs.(p_max_prec[j] .- p_tmp[j]))/n for j in 1:m])/m
		err_pz[k] = sum([sum(abs.(pz_max_prec[j] .- pz_tmp[j]))/n for j in 1:m])/m
		
		diff_ab = [[PeriodReduction(z[j][k] - abtmp[j][k], w1[j], w2[j]) for k in 1:n] for j in 1:m]
		err_ab[k] = sum([sum(abs.(diff_ab[j]))/n for j in 1:m])/m
		##### twice to exclude compile time.	
		t_constr[k] = time_ns()
		Ctmp = [WCurve((g2[j],g3[j]), source = "Invariants", e = er, n = 50) for j in 1:m]
		t_constr[k] = time_ns() - t_constr[k]
		t_func[k] = time_ns()
		wtmp = [Ctmp[j].WeierstrassFunc.(z[j]) for j in 1:m]
		t_func[k] = time_ns() - t_func[k]
		t_ab[k] = time_ns()
		abtmp = [[Ctmp[j].AbelMap(p_max_prec[j][k], pz_max_prec[j][k]) for k in 1:n] for j in 1:m]
		t_ab[k] = time_ns() - t_ab[k]
		s_tmp = [[wtmp[j][k][1] for k in 1:n] for j in 1:m]
		zeta_tmp = [[wtmp[j][k][2] for k in 1:n] for j in 1:m]
		p_tmp = [[wtmp[j][k][3] for k in 1:n] for j in 1:m]
		pz_tmp = [[wtmp[j][k][4] for k in 1:n] for j in 1:m]
		
		err_s[k] = sum([sum(abs.(s_max_prec[j] .- s_tmp[j]))/n for j in 1:m])/m
		err_zeta[k] = sum([sum(abs.(zeta_max_prec[j] .- zeta_tmp[j]))/n for j in 1:m])/m
		err_p[k] = sum([sum(abs.(p_max_prec[j] .- p_tmp[j]))/n for j in 1:m])/m
		err_pz[k] = sum([sum(abs.(pz_max_prec[j] .- pz_tmp[j]))/n for j in 1:m])/m
		
		diff_ab = [[PeriodReduction(z[j][k] - abtmp[j][k], w1[j], w2[j]) for k in 1:n] for j in 1:m]
		err_ab[k] = sum([sum(abs.(diff_ab[j]))/n for j in 1:m])/m
	end
	return (t_constr .* 1e-6 ./ m, t_func .* 1e-6 ./ (m * n), t_ab .* 1e-6 ./ (m * n), err_s, err_zeta, err_p, err_pz, err_ab)
end

function PeriodReduction(z, w1, w2)
	M = [w1 w2; conj(w1) conj(w2)]
	c = inv(M)*[z, conj(z)]
	c = real.(c)
	c = c - round.(c)
	return c[1]*w1 + c[2]*w2
end
