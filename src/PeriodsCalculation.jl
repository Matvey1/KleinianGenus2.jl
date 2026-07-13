function DegeneratePeriods(lcoeff, p, mode)
	if(mode == "right")
		q = sqrt(lcoeff + 0*im)
		w1 = pi*[1/((p[1] - p[2])*(p[1] - p[3])), p[1]/((p[1] - p[2])*(p[1] - p[3]))]*(2*im)/q
		w2 = pi*[1/((p[3] - p[1])*(p[3] - p[2])), p[3]/((p[3] - p[1])*(p[3] - p[2]))]*((-2)*im)/q
		M = lcoeff*[p[1]*p[2]*p[3]*(p[1] + p[2] + p[3]) -p[1]*p[2]*p[3];
					-p[1]*p[2]*p[3] (p[1]*p[2] + p[1]*p[3] + p[2]*p[3])]/2
		eta1 = M*w1
		eta2 = M*w2
		return hcat(w1, w2), hcat(eta1, eta2)
	end
	q = sqrt(lcoeff + 0*im)
	w1 = pi*[1/((p[2] - p[1])*(p[2] - p[3])), p[2]/((p[2] - p[1])*(p[2] - p[3]))]*2*im/q
	w2 = pi*[1/((p[3] - p[1])*(p[3] - p[2])), p[3]/((p[3] - p[1])*(p[3] - p[2]))]*2*im/q
	M = lcoeff*[p[1]*p[2]*p[3]*(p[1] + p[2] + p[3]) -p[1]*p[2]*p[3];
			    -p[1]*p[2]*p[3] (p[1]*p[2] + p[1]*p[3] + p[2]*p[3])]/2
	eta1 = M*w1
	eta2 = M*w2
	return hcat(w1, w2), hcat(eta1, eta2)
end

function ComputeEta(W, E, ArrChi)
	CurrE = copy(E)
	n = length(ArrChi)
	for j in n:-1:1
		CurrE = 2 * CurrE + ArrChi[j] * W
	end
	return CurrE
end

function CalculatePeriods(ArrChi,lcoeff,p;mode = "left")
	W,E = DegeneratePeriods(lcoeff,p,mode)
	E = ComputeEta(W,E,ArrChi)
	return W,E
end
