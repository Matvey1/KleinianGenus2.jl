function DegenerateS(w,lcoeff,p)
    q = sqrt(lcoeff + 0*im)
    M = lcoeff*[p[1]*p[2]*p[3]*(p[1] + p[2] + p[3]) -p[1]*p[2]*p[3];
               -p[1]*p[2]*p[3] (p[1]*p[2] + p[1]*p[3] + p[2]*p[3])]/2
    E = exp(Bilin(M,w,w))/((p[1] - p[2])^2*(p[1] - p[3])^2*(p[2] - p[3])^2)
    s3 = sin(q*im*(p[2] - p[1])*(w[2] - p[3]*w[1])/2)^2
    s2 = sin(q*im*(p[3] - p[1])*(w[2] - p[2]*w[1])/2)^2
    s1 = sin(q*im*(p[2] - p[3])*(w[2] - p[1]*w[1])/2)^2
	S31 = q*im*(p[2] - p[1])*(- p[3]/2)*sin(q*im*(p[2] - p[1])*(w[2] - p[3]*w[1]))
	S21 = q*im*(p[3] - p[1])*(- p[2]/2)*sin(q*im*(p[3] - p[1])*(w[2] - p[2]*w[1]))
	S11 = q*im*(p[2] - p[3])*(- p[1]/2)*sin(q*im*(p[2] - p[3])*(w[2] - p[1]*w[1]))
	S32 = q*im*(p[2] - p[1])*sin(q*im*(p[2] - p[1])*(w[2] - p[3]*w[1]))/2
	S22 = q*im*(p[3] - p[1])*sin(q*im*(p[3] - p[1])*(w[2] - p[2]*w[1]))/2
	S12 = q*im*(p[2] - p[3])*sin(q*im*(p[2] - p[3])*(w[2] - p[1]*w[1]))/2
	
	S = [(-4*E/lcoeff)*((p[3] - p[1])*(p[3] - p[2])*s3 + (p[2] - p[1])*(p[2] - p[3])*s2 + (p[1] - p[2])*(p[1] - p[3])*s1),
			(-4*E/lcoeff)*((p[3] - p[1])*(p[3] - p[2])*(p[1] + p[2])*s3 + (p[2] - p[1])*(p[2] - p[3])*(p[1] + p[3])*s2 + 
																						(p[1] - p[2])*(p[1] - p[3])*(p[2] + p[3])*s1),
			(4*E/lcoeff)*((p[3] - p[1])*(p[3] - p[2])*p[1]*p[2]*s3 + (p[2] - p[1])*(p[2] - p[3])*p[1]*p[3]*s2 + 
																							(p[1] - p[2])*(p[1] - p[3])*p[2]*p[3]*s1),
			2*E*(0.5 * (p[1] - p[2])^2*(p[1] - p[3])^2*(p[2] - p[3])^2 + 
			(p[3] - p[1])*(p[3] - p[2])*(p[1]*p[2]*p[3]*(p[1] + p[2] + p[3]) + p[1]^2*p[2]^2)*s3 +
			(p[2] - p[1])*(p[2] - p[3])*(p[1]*p[2]*p[3]*(p[1] + p[2] + p[3]) + p[1]^2*p[3]^2)*s2 + 
			(p[1] - p[2])*(p[1] - p[3])*(p[1]*p[2]*p[3]*(p[1] + p[2] + p[3]) + p[2]^2*p[3]^2)*s1)]
			
	dS1 = Bilin(M, [2,0], w) .* S .+ 
			[(-4*E/lcoeff)*((p[3] - p[1])*(p[3] - p[2])*S31 + (p[2] - p[1])*(p[2] - p[3])*S21 + (p[1] - p[2])*(p[1] - p[3])*S11),
			(-4*E/lcoeff)*((p[3] - p[1])*(p[3] - p[2])*(p[1] + p[2])*S31 + (p[2] - p[1])*(p[2] - p[3])*(p[1] + p[3])*S21 + 
																						(p[1] - p[2])*(p[1] - p[3])*(p[2] + p[3])*S11),
			(4*E/lcoeff)*((p[3] - p[1])*(p[3] - p[2])*p[1]*p[2]*S31 + (p[2] - p[1])*(p[2] - p[3])*p[1]*p[3]*S21 + 
																							(p[1] - p[2])*(p[1] - p[3])*p[2]*p[3]*S11),
			2*E*((p[3] - p[1])*(p[3] - p[2])*(p[1]*p[2]*p[3]*(p[1] + p[2] + p[3]) + p[1]^2*p[2]^2)*S31 +
			(p[2] - p[1])*(p[2] - p[3])*(p[1]*p[2]*p[3]*(p[1] + p[2] + p[3]) + p[1]^2*p[3]^2)*S21 + 
			(p[1] - p[2])*(p[1] - p[3])*(p[1]*p[2]*p[3]*(p[1] + p[2] + p[3]) + p[2]^2*p[3]^2)*S11)]
			
	dS2 = Bilin(M, [0,2], w) .* S .+ 
			[(-4*E/lcoeff)*((p[3] - p[1])*(p[3] - p[2])*S32 + (p[2] - p[1])*(p[2] - p[3])*S22 + (p[1] - p[2])*(p[1] - p[3])*S12),
			(-4*E/lcoeff)*((p[3] - p[1])*(p[3] - p[2])*(p[1] + p[2])*S32 + (p[2] - p[1])*(p[2] - p[3])*(p[1] + p[3])*S22 + 
																						(p[1] - p[2])*(p[1] - p[3])*(p[2] + p[3])*S12),
			(4*E/lcoeff)*((p[3] - p[1])*(p[3] - p[2])*p[1]*p[2]*S32 + (p[2] - p[1])*(p[2] - p[3])*p[1]*p[3]*S22 + 
																							(p[1] - p[2])*(p[1] - p[3])*p[2]*p[3]*S12),
			2*E*((p[3] - p[1])*(p[3] - p[2])*(p[1]*p[2]*p[3]*(p[1] + p[2] + p[3]) + p[1]^2*p[2]^2)*S32 +
			(p[2] - p[1])*(p[2] - p[3])*(p[1]*p[2]*p[3]*(p[1] + p[2] + p[3]) + p[1]^2*p[3]^2)*S22 + 
			(p[1] - p[2])*(p[1] - p[3])*(p[1]*p[2]*p[3]*(p[1] + p[2] + p[3]) + p[2]^2*p[3]^2)*S12)]
	return (S,dS1,dS2)
end

# Symmetric bilinear contractions; one A*S per component serves all three outputs.
function TransformWeight2(S, dS1, dS2, A, Chi, w)
    T = promote_type(eltype(S), eltype(A), eltype(w))
    R, R1, R2 = zeros(T,4), zeros(T,4), zeros(T,4)
    @inbounds for k in 1:4, i in 1:4
        v = zero(T)
        for j in 1:4
            v += A[k,i,j]*S[j]
        end
        R[k] += S[i]*v
        R1[k] += dS1[i]*v
        R2[k] += dS2[i]*v
    end
    factor = exp(Bilin(Chi,w,w))/A[4,4,4]
    h1 = 2*(Chi[1,1]*w[1]+Chi[1,2]*w[2])
    h2 = 2*(Chi[2,1]*w[1]+Chi[2,2]*w[2])
    return (factor .* R, factor .* (h1 .* R .+ 2 .* R1),
            factor .* (h2 .* R .+ 2 .* R2))
end

function SFuncGeneric(ArrA, ArrChi, lcoeff, p, w)
    length(w) == 2 || throw(ArgumentError("The Abel coordinate must have length two"))
    all(isfinite, w) || throw(ArgumentError("The Abel coordinate must be finite"))
    S = DegenerateS(w, lcoeff, p)
    for j in length(ArrA):-1:1
        S = TransformWeight2(S...,ArrA[j],ArrChi[j],w)
    end
    return S
end

function SigmaDuplication(S)
	return S[1][3]*S[2][2] - S[1][2]*S[2][3] + S[1][4]*S[2][1] - S[1][1]*S[2][4]
end

function SFuncGenericDuplication(SS,w,A2,Chi2,A1,Chi1)
    S = (SS[1] .* [4,4,4,1], SS[2] .* [2,2,2,1/2], SS[3] .* [2,2,2,1/2])
    S = TransformWeight2(S...,A2,Chi2,w)
    return TransformWeight2(S...,A1,Chi1,w)
end

function GenericKleinianDuplication(ArrA, ArrChi, lcoeff, p, AA, ChiChi, w)
	S = SFuncGeneric(ArrA,ArrChi,lcoeff,p,w/2)
	sigma = SigmaDuplication(S)
	S = SFuncGenericDuplication(S,w,AA,ChiChi,ArrA[1],ArrChi[1])
	return (sigma, S[2][1]/(2*S[1][1]), S[3][1]/(2*S[1][1]), S[1][2]/S[1][1], S[1][3]/S[1][1], S[1][4]/S[1][1])
end

function GenericKleinianSigmaSquared(S)
	return (S[1][1], S[2][1]/(2*S[1][1]), S[3][1]/(2*S[1][1]), S[1][2]/S[1][1], S[1][3]/S[1][1], S[1][4]/S[1][1])
end
