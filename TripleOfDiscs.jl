struct Circle
	center
	radius
end

function is_intersecting(C1, C2)
	return ((sign(C1.radius) * sign(C2.radius) * abs(C1.center - C2.center) - C1.radius - C2.radius) <= 0)
end

function belongs_to_circle(z, C)
	return ((C.radius - sign(C.radius)*abs(z - C.center)) >= 0)
end

function common_strip_positive_exists(z, w)
	ang = z/abs(z)
	if(real(w/ang) <= 0)
		return false
	end
	t = (z + w)/abs(z + w)
	if(imag(z/t) > 1 || imag(w/t) > 1)
		return false
	end
	return true
end

function belongs_to_strip(z, a)
	w = z/a
	if(abs(imag(w)) <= 1)
		return true
	end
	return false
end

function belongs_to_strip_positive(z, a)
	w = z/a
	if(abs(imag(w)) <= 1 && real(w) >= 0)
		return true
	end
	return false
end

function circle_on_given_diameter(z, w)
	return Circle((z + w)/2, abs(z-w)/2)
end

function strip_separates(e, a)
	arr = e./a
	if(minimum(abs.(imag.(arr))) <= 1)
		return false
	end
	if(minimum(imag.(arr)) > 0 || maximum(imag.(arr)) < 0)
		return false
	end
	return true
end

function generate_separating_strip_if_exists(e)
	f = z->z/abs(z)
	arr = vcat(f.(e .* im), f.([e[mod1(j+1, 3)] + e[mod1(j+2, 3)] for j in 1:3]), f.([e[mod1(j+1, 3)] - e[mod1(j+2, 3)] for j in 1:3]))
	d = [minimum(abs.(imag.(e ./ arr[j])))*Int64(strip_separates(e, arr[j])) for j in 1:length(arr)]
	if(maximum(d) < 1)
		return nothing
	end
	j = findmax(d)[2]
	return arr[j]
end

function find_separating_circles_given_a_strip(e, dir)
	if(sum(Int64.(imag.(e./dir) .> 0)) == 1)
		dir = -dir
	end
	e = sort(e./dir, by = imag)
	t = max(((real(e[2] - e[3])/2)^2 + imag(e[2])^2)/(2*imag(e[2]) - 2), imag(e[2] + e[3])/2) + 1
	center1 = (real(e[2] + e[3])/2) + t*im
	radius1 = (sqrt((real(e[2] - e[3])/2)^2 + (t - imag(e[2]))^2) + t - 1)/2
	center2 = 0
	eps = min(abs(imag(e[1])) - 1, imag(center1) - radius1 - 1)/4
	radius2 = 1 + eps
	(center3, radius3) = separate_point_from_two_circles_given_direction(e[1], Circle(center1, radius1), Circle(center2, radius2), 1im)
	return (Circle(center1*dir, radius1),Circle(center2*dir, radius2), Circle(center3*dir, -radius3))
end

function find_direction_from_point_to_convex_hull_of_two_circles(z, C1, C2)
	c1 = C1.center - z
	c2 = C2.center - z
	r1 = C1.radius
	r2 = C2.radius
	arr = [c1/abs(c1), c2/abs(c2)]
	if(abs(r1 - r2) <= abs(c1 - c2))
		arr = vcat(arr, [(c2 - c1)/(r2 - r1 + im*sqrt(abs(c2 - c1)^2 - (r2 - r1)^2)), (c2 - c1)/(r2 - r1 - im*sqrt(abs(c2 - c1)^2 - (r2 - r1)^2))])
	end
	d = [min(real(c1/arr[j]) - r1, real(c2/arr[j]) - r2) for j in 1:length(arr)]
	(eps, j) = findmax(d)
	return (arr[j], eps)
end

function separate_point_from_two_circles_given_direction(z, C1, C2, dir)
	c1 = C1.center/dir - z/dir
	r1 = C1.radius
	c2 = C2.center/dir - z/dir
	r2 = C2.radius
	eps = min(real(c1) - r1, real(c2) - r2)/2
	t = max((imag(c1)^2 + real(c1)^2 - (eps + r1)^2)/(2*real(c1) - 2*r1 - 2*eps), 
			(imag(c2)^2 + real(c2)^2 - (eps + r2)^2)/(2*real(c2) - 2*r2 - 2*eps)) +  1
	r = t - eps
	c = dir*(t + z/dir)
	return (c,r)
end

function separate_point_from_two_circles(z, C1, C2)
	(dir,) = find_direction_from_point_to_convex_hull_of_two_circles(z, C1, C2)
	return separate_point_from_two_circles_given_direction(z, C1, C2, dir)
end

function small_inflation_and_separation(z, C1, C2)
	(dir, eps) = find_direction_from_point_to_convex_hull_of_two_circles(z, C1, C2)
	t = min(eps, abs(C1.center - C2.center) - C1.radius - C2.radius)/4
	D1 = Circle(C1.center, C1.radius + t)
	D2 = Circle(C2.center, C2.radius + t)
	(c, r) = separate_point_from_two_circles(z, D1, D2)
	return (D1, D2, Circle(c, -r))
end

function find_strips_with_given_point_on_boundary(z)
	return (z/(sqrt(abs(z)^2 - 1) + 1im), z/(sqrt(abs(z)^2 - 1) - 1im))
end

function check_availability_of_tangent_circle(e, dir, j, k)
	a = e[j]/dir
	b = e[k]/dir
	l = setdiff([1,2,3], [j,k])[1]
	c = e[l]/dir
	if(real(a) >= real(b) && imag(b)*imag(c) < 0)
		return true
	end
	return false
end

function construct_tangent_circle(e, dir, j, k)
	a = e[j]/dir
	b = e[k]/dir
	z = (a + b)/2
	v = (a - b)*im*sign(imag(b))
	rts = roots([abs(a - b)^2/4 - (imag(a+b)/2 + sign(imag(b)))^2, -2*imag(v)*(imag(a + b)/2 + sign(imag(b))), (real(v))^2])
	arr = sort(real.(rts))
	arr = arr[arr .> 0]
	t = minimum(arr)
	center = (a + b)/2 + t*v
	radius = sqrt(t^2 * abs(v)^2 + abs(a - b)^2/4)
	return Circle(dir*center, radius)
end

function separate(e)
	dir = generate_separating_strip_if_exists(e)
	if(!isnothing(dir))
		return find_separating_circles_given_a_strip(e, dir)
	end
	CD = [Circle((e[mod1(j+1, 3)] + e[mod1(j+2, 3)])/2, abs(e[mod1(j+1, 3)] - e[mod1(j+2, 3)])/2) for j in 1:3]
	refCirc = Circle(0, 1)
	ib = [Int64(is_intersecting(CD[j], refCirc)) for j in 1:3]
	b = [Int64(common_strip_positive_exists(e[mod1(j+1, 3)], e[mod1(j+2, 3)])) for j in 1:3]
	if(sum(b .* ib) > 0)
		j = findmax(b .* ib)[2]
		C1 = Circle((e[mod1(j+1, 3)] + e[mod1(j+2, 3)])/2, abs(e[mod1(j+1, 3)] - e[mod1(j+2, 3)])/2)
		C2 = Circle(-im*sign(imag(C1.center))/sqrt(3), 2/sqrt(3))
		return small_inflation_and_separation(e[j], C1, C2)
	end
	if(sum(b) == 3)
		j = findmax(abs.(e))[2]
		C1 = Circle((e[mod1(j+1, 3)] + e[mod1(j+2, 3)])/2, abs(e[mod1(j+1, 3)] - e[mod1(j+2, 3)])/2)
		C2 = Circle(0, 1)
		return small_inflation_and_separation(e[j], C1, C2)
	end
	j = findmin(b)[2]
	for i in 1:2
		d = find_strips_with_given_point_on_boundary(e[mod1(j+i, 3)])
		for k in 1:2
			if(check_availability_of_tangent_circle(e, d[k], j, mod1(j+i, 3)))
				C1 = construct_tangent_circle(e, d[k], j, mod1(j+i, 3))
				l = setdiff([1,2,3], [j,mod1(j+i, 3)])[1]
				C2 = Circle(0, 1)
				return small_inflation_and_separation(e[l], C1, C2)
			end
		end
	end
	k = findmin(abs.([e[mod1(j+1, 3)], e[mod1(j+2, 3)]]))[2]
	k = mod1(j+k, 3)
	C1 = Circle((e[j] + e[k])/2, abs(e[j] - e[k])/2)
	C2 = Circle(0, 1)
	l = setdiff([1,2,3], [j,k])[1]
	return small_inflation_and_separation(e[l], C1, C2)
end

#=function check_separation(C,e)
	arr = vcat(e, [-1,1])
	if(is_intersecting(C[1], C[2]) || is_intersecting(C[1], C[3]) || is_intersecting(C[2], C[3]))
		return false
	end
	ind = [Int64(belongs_to_circle(arr[j], C[i])) for i in 1:3, j in 1:5]
	s = [sum(ind[i, :]) for i in 1:3]
	if(sum(s) != 5 || maximum(s) != 2)
		return false
	end
	j = findmin(s)[2]
	if(C[j].radius > 0)
		return false
	end
	return true
end=#

function check_separation(C, arr)
	if(is_intersecting(C[1], C[2]) || is_intersecting(C[1], C[3]) || is_intersecting(C[2], C[3]))
		return 1
	end
	ind = [Int64(belongs_to_circle(arr[j], C[i])) for i in 1:3, j in 1:5]
	s = [sum(ind[i, :]) for i in 1:3]
	if(sum(s) != 5 || maximum(s) != 2)
		return 2
	end
	j = findmin(s)[2]
	if(C[j].radius > 0)
		return 2
	end
	return 0
end

function generateTripleOfDiscs(R)
	M = [abs(R[j] - R[k]) for j in 1:5, k in 1:5]
	m = maximum(M)
	for i in 1:5
		M[i,i] = m + 1
	end
	ind = findmin(M)[2]
	j = ind[1]
	k = ind[2]
	a = 2/(R[k] - R[j])
	b = -(R[k] + R[j])/(R[k] - R[j])
	f = z->(a*z + b)
	remaining_ind = setdiff(collect(1:5), [j,k])
	e = f.(R[remaining_ind])
	C = separate(e)
	return [((C[j].center - b)/a, C[j].radius/abs(a)) for j in 1:3]
end
