using LinearAlgebra
using Printf
using SpecialFunctions

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

struct Atom
	Z::Int
	center::NTuple{3, Float64}
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

H1=Atom(1, (0.0, 0.0, -0.529))
H2=Atom(1, (0.0, 0.0, +0.529))


H1Basis1s=Basis([STO6G1, STO6G2, STO6G3, STO6G4, STO6G5, STO6G6], (0.0, 0.0, -0.529), 1)
H2Basis1s=Basis([STO6G1, STO6G2, STO6G3, STO6G4, STO6G5, STO6G6], (0.0, 0.0, +0.529), 1)

BasisSet=[H1Basis1s, H2Basis1s]


function STVij(GTF1::Basis, GTF2::Basis, AllAtom::Vector{Atom})
	Sij=0.0
	Tij=0.0
	Vij=0.0
	R1=GTF1.Center
	R2=GTF2.Center
	R12=sqrt((R1[1]-R2[1])^2+(R1[2]-R2[2])^2+(R1[3]-R2[3])^2)
	for pgtf1 in GTF1.GTFs
		for pgtf2 in GTF2.GTFs
			alpha1, coeff1, norms1=pgtf1.alpha, pgtf1.coeff, pgtf1.norms
			alpha2, coeff2, norms2=pgtf2.alpha, pgtf2.coeff, pgtf2.norms
			p=alpha1+alpha2
			m=alpha1*alpha2/p
			K12=exp(-m*R12^2)
			S00=(π/p)^(3/2)*K12
			Sij+=coeff1*coeff2*norms1*norms2*S00
			Tij+=coeff1*coeff2*norms1*norms2*m*(3-2*m*R12^2)*S00
			Rxc=(alpha1*R1[1]+alpha2*R2[1])/p
			Ryc=(alpha1*R1[2]+alpha2*R2[2])/p
			Rzc=(alpha1*R1[3]+alpha2*R2[3])/p
			for atom in AllAtom
				Ra=atom.center
				Rac=sqrt((Ra[1]-Rxc)^2+(Ra[2]-Ryc)^2+(Ra[3]-Rzc)^2)
				Vij+=-atom.Z*coeff1*coeff2*norms1*norms2*K12*(2π/p)*BoysF0(p*Rac^2)
			end
		end
	end
	return Sij, Tij, Vij
end

function BoysF0(t::Float64)
	if t==0
		return 1.0
	elseif t<10e-8
		return 1.0 - t/3 + t^2/10 - t^3/42 + t^4/216
	else
		return 0.5*sqrt(π/t)*erf(sqrt(t))
	end
end

function CalcERI(Index1::Int, Index2::Int, Index3::Int, Index4::Int)
	eri = 0.0
	GTF1 = BasisSet[Index1]
	GTF2 = BasisSet[Index2]
	GTF3 = BasisSet[Index3]
	GTF4 = BasisSet[Index4]
	Ri, Rj, Rk, Rl = GTF1.Center, GTF2.Center, GTF3.Center, GTF4.Center
	Rijsq = sum((Ri .- Rj) .^ 2)
	Rklsq = sum((Rk .- Rl) .^ 2)
	for i in GTF1.GTFs
		for j in GTF2.GTFs
			alpha_i, coeff_i, norm_i = i.alpha, i.coeff, i.norms
			alpha_j, coeff_j, norm_j = j.alpha, j.coeff, j.norms
			p1 = alpha_i + alpha_j
			R1 = (alpha_i .* Ri .+ alpha_j .* Rj) ./ p1
			Kij = exp(-alpha_i * alpha_j / p1 * Rijsq)
			for k in GTF3.GTFs
				for l in GTF4.GTFs
					alpha_k, coeff_k, norm_k = k.alpha, k.coeff, k.norms
					alpha_l, coeff_l, norm_l = l.alpha, l.coeff, l.norms

					p2 = alpha_k + alpha_l
					R2 = (alpha_k .* Rk .+ alpha_l .* Rl) ./ p2
					Kkl = exp(-alpha_k * alpha_l / p2 * Rklsq)

					R12sq = sum((R1 .- R2) .^ 2)
					t = (p1 * p2) / (p1 + p2) * R12sq

					prefactor = (2 * π^2.5) / (p1 * p2 * sqrt(p1 + p2))
					prim_eri = prefactor * Kij * Kkl * BoysF0(t)

					norms_prod = norm_i * norm_j * norm_k * norm_l
					coeffs_prod = coeff_i * coeff_j * coeff_k * coeff_l

					eri += norms_prod * coeffs_prod * prim_eri
				end
			end
		end
	end

	return eri
end



Num=length(BasisSet)

S=[STVij(BasisSet[i], BasisSet[j], [H1, H2])[1] for i in 1:Num, j in 1:Num]
T=[STVij(BasisSet[i], BasisSet[j], [H1, H2])[2] for i in 1:Num, j in 1:Num]
V=[STVij(BasisSet[i], BasisSet[j], [H1, H2])[3] for i in 1:Num, j in 1:Num]
ERI=[CalcERI(i, j, k, l) for i in 1:Num, j in 1:Num, k in 1:Num, l in 1:Num]
Hcore=T+V





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