module CalcS

using ..Definitions: PGTF, Basis
using SpecialFunctions: binomial, gamma

export Sij


function S1D(l1::Int, l2::Int, A::Float64, B::Float64, alpha1::Float64, alpha2::Float64)
	p=alpha1+alpha2
	u=(alpha1 * alpha2) / p
	P = (alpha1 * A + alpha2 * B) / p
	PA = P - A
	PB = P - B
	PAB = A - B
	S00=(pi/p)^(0.5)*exp(-u * PAB^2)
	if l1<0 || l2<0
		return 0.0
	elseif l1==0 && l2==0
		return S00
	elseif l1==0 && l2>=1
		return PB*S1D(0, l2-1, A, B, alpha1, alpha2)+(l2-1)*S1D(0, l2-2, A, B, alpha1, alpha2)/(2*p)
	elseif l1>=1
		return PA*S1D(l1-1, l2, A, B, alpha1, alpha2)+(l1-1)*S1D(l1-2, l2, A, B, alpha1, alpha2)/(2*p)+l2*S1D(l1-1, l2-1, A, B, alpha1, alpha2)/(2*p)
	end
end

function Sij(basis1::Basis, basis2::Basis)
	S_total = 0.0

	R1 = basis1.position
	R2 = basis2.position
	l1 = basis1.Type
	l2 = basis2.Type
	for pgtf1 in basis1.GTFs
		alpha1 = pgtf1.alpha
		coeff1 = pgtf1.coeff
		norm1 = pgtf1.norms
		for pgtf2 in basis2.GTFs
			alpha2 = pgtf2.alpha
			coeff2 = pgtf2.coeff
			norm2 = pgtf2.norms


			Sx = S1D(l1[1], l2[1], R1[1], R2[1], alpha1, alpha2)
			Sy = S1D(l1[2], l2[2], R1[2], R2[2], alpha1, alpha2)
			Sz = S1D(l1[3], l2[3], R1[3], R2[3], alpha1, alpha2)

			SPrim = Sx * Sy * Sz

			S_total += coeff1 * coeff2 * norm1 * norm2 * SPrim
		end
	end

	return S_total
end

end
