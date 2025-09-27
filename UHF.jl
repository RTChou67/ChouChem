using LinearAlgebra
using Printf
using SpecialFunctions

include("Definitions.jl")
include("GetBasisList.jl")
include("CalcS.jl")
include("CalcT.jl")
include("CalcV.jl")
include("CalcG.jl")
using .Definitions: PGTF, CGTF, Basis, Atom
using .GetBasisList: generate_basis_list, get_basis_set
using .CalcS: Sij
using .CalcT: Tij
using .CalcV: Vij
using .CalcG: Gijkl



MolInAng = [
	Atom("H", 1, "STO-3G", (0.0, 0.0, +1.0)),
	Atom("F", 9, "STO-3G", (0.0, 0.0, -1.0)),
]
Charge=1
Multiplicity=2

Bohr2Ang=0.52917721092
MolInBohr = [Atom(atom.symbol, atom.Z, atom.basis_set, atom.position ./ Bohr2Ang) for atom in MolInAng]
Molecule = MolInBohr

BasisSet = generate_basis_list(Molecule)



Num=length(BasisSet)
ENum=sum(atom.Z for atom in Molecule)-Charge
NNum=length(Molecule)

ENumAlpha=(ENum+Multiplicity-1)÷2
ENumBeta=(ENum-Multiplicity+1)÷2





S=[Sij(BasisSet[i], BasisSet[j]) for i in 1:Num, j in 1:Num]
T=[Tij(BasisSet[i], BasisSet[j]) for i in 1:Num, j in 1:Num]
V=[Vij(BasisSet[i], BasisSet[j], Molecule) for i in 1:Num, j in 1:Num]
ERI=[Gijkl(BasisSet[i], BasisSet[j], BasisSet[k], BasisSet[l]) for i in 1:Num, j in 1:Num, k in 1:Num, l in 1:Num]




function format_to_custom_eng(x::Float64)
	if x == 0.0
		return "0.000000E+00"
	end
	std_eng_str = @sprintf("%1.6e", x)
	parts = split(std_eng_str, 'e')
	mantissa_std = parse(Float64, parts[1])
	exponent_std = parse(Int, parts[2])
	mantissa_custom = mantissa_std / 10.0
	exponent_custom = exponent_std + 1
	sign_char = mantissa_std < 0 ? "-" : ""
	mantissa_str = @sprintf("%.6f", abs(mantissa_custom))
	exponent_str = @sprintf("%+03d", exponent_custom)
	return string(sign_char, mantissa_str, "E", exponent_str)
end

function print_formatted_matrix(matrix::Matrix{Float64})
	n = size(matrix, 1)
	labels = [string(i) for i in 1:n]
	@printf("%-5s", "")
	for j in 1:n
		@printf("%15s", labels[j])
	end
	println()
	println(repeat("-", 5 + 15*n))
	for i in 1:n
		@printf("%-5s", labels[i] * "|")
		for j in 1:n
			formatted_num = format_to_custom_eng(matrix[i, j])
			@printf("%15s", formatted_num)
		end
		println()
	end
end


println("\nOverlap Matrix S:")
print_formatted_matrix(S)
println("\nKinetic Energy Matrix T:")
print_formatted_matrix(T)
println("\nNuclear Attraction Matrix V:")
print_formatted_matrix(V)




Hcore=T+V

println("\nCore Hamiltonian Hcore:")
print_formatted_matrix(Hcore)


F=Hcore

MaxIter=100
Threshold=1e-10


X=S^(-0.5)
Fprime=X'*F*X
E, Cprime=eigen(Fprime)
C=X*Cprime
p = sortperm(E)
E=E[p]
Cprime=Cprime[:, p]
PAlpha=[sum(C[i, m]*C[j, m] for m in 1:ENumAlpha) for i in 1:Num, j in 1:Num]
PBeta=[sum(C[i, m]*C[j, m] for m in 1:ENumBeta) for i in 1:Num, j in 1:Num]



for i in 1:MaxIter
	global PAlpha, PBeta, FAlpha, FBeta, FprimeAlpha, FprimeBeta, EAlpha, EBeta, CprimeAlpha, CprimeBeta, CAlpha, CBeta, GAlpha, GBeta, PnewAlpha, PnewBeta
	Ptotal=PAlpha+PBeta
	#global PAlpha, PBeta, FAlpha, FBeta, Fprime, E, Cprime, C, G, Pnew, p
	GAlpha = [sum(Ptotal[k, l]*ERI[i, j, k, l] - PAlpha[k, l]*ERI[i, l, k, j] for k in 1:Num, l in 1:Num) for i in 1:Num, j in 1:Num]
	GBeta = [sum(Ptotal[k, l]*ERI[i, j, k, l] - PBeta[k, l]*ERI[i, l, k, j] for k in 1:Num, l in 1:Num) for i in 1:Num, j in 1:Num]
	FAlpha=Hcore+GAlpha
	FBeta=Hcore+GBeta
	FprimeAlpha=X'*FAlpha*X
	FprimeBeta=X'*FBeta*X
	EAlpha, CprimeAlpha=eigen(FprimeAlpha)
	EBeta, CprimeBeta=eigen(FprimeBeta)
	pAlpha=sortperm(EAlpha)
	pBeta=sortperm(EBeta)
	EAlpha=EAlpha[pAlpha]
	EBeta=EBeta[pBeta]
	CprimeAlpha=CprimeAlpha[:, pAlpha]
	CprimeBeta=CprimeBeta[:, pBeta]
	CAlpha=X*CprimeAlpha
	CBeta=X*CprimeBeta
	PnewAlpha=[sum(CAlpha[i, m]*CAlpha[j, m] for m in 1:ENumAlpha) for i in 1:Num, j in 1:Num]
	PnewBeta=[sum(CBeta[i, m]*CBeta[j, m] for m in 1:ENumBeta) for i in 1:Num, j in 1:Num]
	delta_PAlpha = sqrt(sum((PnewAlpha - PAlpha) .^ 2))
	delta_PBeta = sqrt(sum((PnewBeta - PBeta) .^ 2))
	delta_P = max(delta_PAlpha, delta_PBeta)
	if delta_P < Threshold
		println("SCF")
		PAlpha=PnewAlpha
		PBeta=PnewBeta
		break
	end
	PAlpha=PnewAlpha
	PBeta=PnewBeta
end

Ptotal=PAlpha+PBeta

VNN = sum(Molecule[i].Z * Molecule[j].Z / sqrt(sum((Molecule[i].position .- Molecule[j].position) .^ 2)) for i in 1:NNum for j in (i+1):NNum)
print_formatted_matrix(Ptotal)

Ee = 0.5 * sum(PAlpha .* (Hcore + FAlpha)) + 0.5 * sum(PBeta .* (Hcore + FBeta))
Etot = Ee + VNN


Te = sum(Ptotal .* T)
VNe = sum(Ptotal .* V)
Vee = 0.5 * (sum(PAlpha .* GAlpha) + sum(PBeta .* GBeta))
Etot_components = Te + VNe + Vee + VNN


@printf("\n--- Molecular Structure ---\n")
for i in 1:NNum
	atom = Molecule[i]
	@printf("Atom %d: %s (Z=%d) at (%+.6f, %+.6f, %+.6f) Å\n", i, atom.symbol, atom.Z, atom.position[1]*0.52918, atom.position[2]*0.52918, atom.position[3]*0.52918)
end



@printf("\n--- Final Energy Results ---\n")
@printf("Total Energy      = %.10f Hartree\n", Etot)
@printf("Electronic Energy = %.10f Hartree\n", Ee)
@printf("Nuclear Repulsion =  %.10f Hartree\n", VNN)
#