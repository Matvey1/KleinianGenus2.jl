struct Genus2WCurveDivisor{T<:AbstractFloat}
	P1::Array{Complex{T}}
	P2::Array{Complex{T}}
	
	function Genus2WCurveDivisor()
		return new{Float64}([], [])
	end
	function Genus2WCurveDivisor(p1)
		return new{typeof(real(p1[1]))}([p1[1], p1[2]], [])
	end
	function Genus2WCurveDivisor(p1,p2)
		return new{typeof(real(p1[1]))}([p1[1], p1[2]], [p2[1], p2[2]])
	end
end

Base.show(io::IO,D::Genus2WCurveDivisor) = print(io, Genus2WCurveDivisorString(D))

function Genus2WCurveDivisorString(D::Genus2WCurveDivisor)
	s = "Genus2WCurveDivisor{" * string(typeof(real(D.P1[1]))) * "}\n"
	if(length(D.P2) > 0)
		s = s * "P1 + P2 - 2*inf, where\n"
		s = s * "P1 = " * string(D.P1) * "\n"
		s = s * "P2 = " * string(D.P2) * "\n"
		return s
	end
	if(length(D.P1) > 0)
		s = s * "P1 - inf, where\n"
		s = s * "P1 = " * string(D.P1) * "\n"
		return s
	end
	s = s*"0\n"
	return s
end

@inline function Rexpr(p,q,r,pp,qq,rr,a,b,c,d,e)
	return 0.5*(pp[a]*p[b]*p[c]*(q[d]*r[e] + q[e]*r[d]) +
				qq[a]*q[b]*q[c]*(p[d]*r[e] + p[e]*r[d]) + 
				rr[a]*r[b]*r[c]*(p[d]*q[e] + p[e]*q[d]))
end

function RichelotTriple(p,q,r)
	PP = chop(derivative(q)*r - derivative(r)*q)
	QQ = chop(derivative(r)*p - derivative(p)*r)
	RR = chop(derivative(p)*q - derivative(q)*p)
	return (PP,QQ,RR)
end

@inline function Dis(p)
	return p[1]^2 - 4*p[0]*p[2]
end

@inline function DDD(p,q,r)
	return det([p[0] q[0] r[0]; p[1] q[1] r[1]; p[2] q[2] r[2]])
end

@inline function Bilin(A,x,y)
	return sum(y.*(A[:,:]*x))
end

function DiskTest(disc)
	t = z->(abs(z - disc[1])*sign(disc[2]) < disc[2])
	return t
end

function GenerateSortingMethod(triple_of_discs)
	cr = [DiskTest(triple_of_discs[i]) for i in 1:3]
	sm = e->[[e[j] for j in 1:length(e) if cr[i](e[j])] for i in 1:3]
end

function SplitToLimitRoots(R)
	return [(R[1][1] + R[1][2])/2, (R[2][1] + R[2][2])/2, (R[3][1] + R[3][2])/2]
end

function LeftRealSortMethod(e)
	E = sort(real.(e))
	s = length(E)
	return [E[1:s-4], E[s-3:s-2], E[s-1:s]]
end

function RightRealSortMethod(e)
	E = sort(real.(e))
	s = length(E)
	return [E[1:2], E[3:4], E[5:s]]
end

function HalfPeriodCoords(c,u,v)
	return [u[2], -u[1], -u[0], -c*(v[0]*u[2]^2 + u[2]*u[0]*v[2] + u[0]^2*v[4])/4]
end

function AuxDiff(T)
	if(length(T) == 1)
		return typeof(real(T[1]))(1.0)
	end
	return T[1] - T[2]
end

function KummerCoord(D::Genus2WCurveDivisor, f::Polynomial)
	if(length(D.P2) > 0)
		x1 = D.P1[1]
		x2 = D.P2[1]
		y1 = D.P1[2]
		y2 = D.P2[2]
		F = 2*f[0] + f[1]*(x1 + x2) + 2*f[2]*x1*x2 + f[3]*x1*x2*(x1 + x2) + 2*f[4]*x1^2*x2^2 + f[5]*x1^2*x2^2*(x1 + x2) + 2*f[6]*x1^3*x2^3
		return [1, x1 + x2, -x1*x2, (F - 2 * y1 * y2)/(4*(x1 - x2)^2)]
	end
	if(length(D.P1) > 0)
		return [0, 1, -D.P1[1], (D.P1[1])^2]
	end
	return typeof(f[0]).([0, 0, 0, 1.0])
end

function truedist(u,v)
	return norm(u)^2 + norm(v)^2 - 2*abs(dot(u,v))
end

function KummerEquationCheck(v, f)
	return (16*(v[2]^2 + 4*v[3]*v[1])*v[4]^2 - 8*v[4]*(2f[0]*v[1]^3 + f[1]*v[2]*v[1]^2 - 2*f[2]*v[3]*v[1]^2 - f[3]*v[3]*v[2]*v[1]
	+ 2*f[4]*v[3]^2*v[1] + f[5]*v[2]*v[3]^2 - 2*f[6]*v[3]^3) + (f[1]^2 - 4*f[0]*f[2])*v[1]^4 - 4*f[0]*f[3]*v[2]*v[1]^3 + 2*f[1]*f[3]*v[3]*v[1]^3
		- 4*f[0]*f[4]*v[2]^2*v[1]^2 - 4*(f[0]*f[5] - f[1]*f[4])*v[2]*v[3]*v[1]^2 + (f[3]^2 + 2*f[1]*f[5] - 4*f[2]*f[4] - 4*f[0]*f[6])*v[3]^2*v[1]^2
		- 4*f[0]*f[5]*v[2]^3*v[1] - 4*(2*f[0]*f[6] - f[1]*f[5])*v[2]^2*v[3]*v[1] + 4*(f[1]*f[6] - f[2]*f[5])*v[2]*v[3]^2*v[1] 
		+ 2*f[3]*f[5]*v[3]^3*v[1] - 4*f[0]*f[6]*v[2]^4 + 4*f[1]*f[6]*v[2]^3*v[3] - 4*f[2]*f[6]*v[2]^2*v[3]^2 + 4*f[3]*f[6]*v[2]*v[3]^3
		+ (f[5]^2 - 4*f[4]*f[6])*v[3]^4)
end

function GenerateRandomPoint(Roots)
	T = typeof(real(Roots[1]))
	while(true)
		a = 10*(rand(Complex{T}) - 0.5 - 0.5im)
		if( minimum([abs(Roots[j] - a) for j in 1:length(Roots)]) > 1)
			return a
		end
	end
end

function normalizeToPeriods(Per, z)
  a = zeros(typeof(real(z[1])), 4)
  a[1] = real(z[1])
  a[2] = imag(z[1])
  a[3] = real(z[2])
  a[4] = imag(z[2])
  M = zeros(typeof(real(z[1])), 4, 4)
  for i in 0:1
    for j in 1:4
      M[2*i+1, j] = real(Per[i+1 , j])
      M[2*i+2, j] = imag(Per[i+1, j])
    end
  end
  
  w = inv(M) * a
  for i in 1:4
    w[i] = round(w[i])
  end

  return z - (w[1] * Per[:, 1] + w[2] * Per[:, 2] + w[3] * Per[:, 3] + w[4] * Per[:, 4])
end

function verifySortingMethod(Roots, triple_of_discs)
	C = [Circle(triple_of_discs[i][1], triple_of_discs[i][2]) for i in 1:3]
	return check_separation(C, Roots)
end
