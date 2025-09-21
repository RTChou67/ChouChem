using LinearAlgebra
using Printf

struct PGTF
	alpha::Float64
	coeff::Float64
	norms::Float64
end

struct Basis
	GTFs::Vector{PGTF}
	Center::NTuple{3, Float64}
	Z::Int64
end

function NormalizePGTF(alpha::Float64)
	return (2alpha/π)^(3/4)
end

PGTF(alpha::Float64, coeff::Float64)=PGTF(alpha, coeff, NormalizePGTF(alpha))

STO6G1=PGTF(0.3552322122e02, 0.9163596281e-02)
STO6G2=PGTF(0.6513143725e01, 0.4936149294e-01)
STO6G3=PGTF(0.1822142904e01, 0.1685383049e+00)
STO6G4=PGTF(0.6259552659e00, 0.3705627997e+00)
STO6G5=PGTF(0.2430767471e00, 0.4164915298e+00)
STO6G6=PGTF(0.1001124280e00, 0.1303340841e+00)

H1Basis1s=Basis([STO6G1, STO6G2, STO6G3, STO6G4, STO6G5, STO6G6], (0.0, 0.0, -0.529), 1)
H2Basis1s=Basis([STO6G1, STO6G2, STO6G3, STO6G4, STO6G5, STO6G6], (0.0, 0.0, +0.529), 1)

BasisSet=[H1Basis1s, H2Basis1s]


function STij(GTF1::Basis, GTF2::Basis)
	Sij=0.0
	Tij=0.0
	R1=GTF1.Center
	R2=GTF2.Center
	R12=sqrt((R1[1]-R2[1])^2+(R1[2]-R2[2])^2+(R1[3]-R2[3])^2)
	for pgtf1 in GTF1.GTFs
		for pgtf2 in GTF2.GTFs
			alpha1=pgtf1.alpha
			alpha2=pgtf2.alpha
			p=alpha1+alpha2
			m=alpha1*alpha2/p
			S00=(π/p)^(3/2)*exp(-m*R12^2)
			Sij+=pgtf1.coeff*pgtf2.coeff*pgtf1.norms*pgtf2.norms*S00
			Tij+=pgtf1.coeff*pgtf2.coeff*pgtf1.norms*pgtf2.norms*m*(3-2*m*R12^2)*S00
		end
	end
	return Sij, Tij
end








Num=length(BasisSet)

S=[STij(BasisSet[i], BasisSet[j])[1] for i in 1:Num, j in 1:Num]
T=[STij(BasisSet[i], BasisSet[j])[2] for i in 1:Num, j in 1:Num]


#=
Calc S
Calc T
Calc V
Calc 2e
Calc P
for i in Iter:
	Calc F
	Calc X
	Orthogonalize F=F'
	F'C'=C'E
	C=XC'
	Calc P, E
	if DeltaE<10e-8
		break
	end
end
=#