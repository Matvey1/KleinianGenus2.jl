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

function SFuncGeneric(ArrA, ArrChi, lcoeff, p, w)
	(S, dS1, dS2) = DegenerateS(w, lcoeff, p)
	for j in length(ArrA):-1:1
		R = [Bilin(ArrA[j][k,:,:], S, S) for k in 1:4]
		R1 = [Bilin(ArrA[j][k,:,:], S, dS1) for k in 1:4]
		R2 = [Bilin(ArrA[j][k,:,:], S, dS2) for k in 1:4]
		S = (exp(Bilin(ArrChi[j], w, w))/ArrA[j][4,4,4]) .* R
		dS1 = (exp(Bilin(ArrChi[j], w, w))/ArrA[j][4,4,4]) .* (Bilin(ArrChi[j], [2,0], w) .* R + 2 .* R1)
		dS2 = (exp(Bilin(ArrChi[j], w, w))/ArrA[j][4,4,4]) .* (Bilin(ArrChi[j], [0,2], w) .* R + 2 .* R2)
	end
	return (S, dS1, dS2)
end

function SigmaDuplication(S)
	return S[1][3]*S[2][2] - S[1][2]*S[2][3] + S[1][4]*S[2][1] - S[1][1]*S[2][4]
end

function SFuncGenericDuplication(SS,w,A2,Chi2,A1,Chi1)
	S, dS1, dS2 = SS[1] .* [4, 4, 4, 1], SS[2] .* [2, 2, 2, 1/2], SS[3] .* [2, 2, 2, 1/2]
	R = [Bilin(A2[k,:,:], S, S) for k in 1:4]
	R1 = [Bilin(A2[k,:,:], S, dS1) for k in 1:4]
	R2 = [Bilin(A2[k,:,:], S, dS2) for k in 1:4]
	S = (exp(Bilin(Chi2, w, w))/A2[4,4,4]) .* R
	dS1 = (exp(Bilin(Chi2, w, w))/A2[4,4,4]) .* (Bilin(Chi2, [2,0], w) .* R + 2 .* R1)
	dS2 = (exp(Bilin(Chi2, w, w))/A2[4,4,4]) .* (Bilin(Chi2, [0,2], w) .* R + 2 .* R2)
	R = [Bilin(A1[k,:,:], S, S) for k in 1:4]
	R1 = [Bilin(A1[k,:,:], S, dS1) for k in 1:4]
	R2 = [Bilin(A1[k,:,:], S, dS2) for k in 1:4]
	S = (exp(Bilin(Chi1, w, w))/A1[4,4,4]) .* R
	dS1 = (exp(Bilin(Chi1, w, w))/A1[4,4,4]) .* (Bilin(Chi1, [2,0], w) .* R + 2 .* R1)
	dS2 = (exp(Bilin(Chi1, w, w))/A1[4,4,4]) .* (Bilin(Chi1, [0,2], w) .* R + 2 .* R2)
	return (S, dS1, dS2)
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
