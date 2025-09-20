struct PGTF
	alpha::Float64
	coeff::Float64
	norms::Float64
end

struct Atom
	Z::Int
	center::NTuple{3, Float64}
end

struct CGTF
	GTFs::Vector{PGTF}
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

STO6GH=CGTF([STO6G1, STO6G2, STO6G3, STO6G4, STO6G5, STO6G6])

struct AtomicBasis
	atom::Atom
	basis::Vector{CGTF}
end


H1=Atom(1, (0.0, 0.0, -0.529))
H2=Atom(1, (0.0, 0.0, +0.529))

H1AtomicBasis=AtomicBasis(H1, [STO6GH])
H2AtomicBasis=AtomicBasis(H2, [STO6GH])


function Sij(Atom1::Atom, Atom2::Atom, GTF1::CGTF, GTF2::CGTF)
	println(length(GTF1.GTFs))
	println(length(GTF2.GTFs))
end

Sij(H1, H2, STO6GH, STO6GH)


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