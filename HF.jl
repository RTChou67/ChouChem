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
	Atom("H", 1, "STO-3G", (0.0, 0.0, +1)),
	Atom("O", 8, "STO-3G", (0.0, +1, 0.0)),
	Atom("H", 1, "STO-3G", (0.0, 0.0, -1)),
]

Bohr2Ang=0.52917721092
MolInBohr = [Atom(atom.symbol, atom.Z, atom.basis_set, atom.position ./ Bohr2Ang) for atom in MolInAng]
Molecule = MolInBohr

BasisSet = generate_basis_list(Molecule)



Num=length(BasisSet)
ENum=sum(atom.Z for atom in Molecule)
NNum=length(Molecule)

S=[Sij(BasisSet[i], BasisSet[j]) for i in 1:Num, j in 1:Num]
T=[Tij(BasisSet[i], BasisSet[j]) for i in 1:Num, j in 1:Num]
V=[Vij(BasisSet[i], BasisSet[j], Molecule) for i in 1:Num, j in 1:Num]
ERI=[Gijkl(BasisSet[i], BasisSet[j], BasisSet[k], BasisSet[l]) for i in 1:Num, j in 1:Num, k in 1:Num, l in 1:Num]



#=
function print_formatted_matrix(matrix::Matrix{Float64})
	n = size(matrix, 1)
	labels = [string(i) for i in 1:n] # 动态生成标签 "1", "2", ...

	# 打印列标题
	@printf("%-5s", "") # 左对齐的标签列
	for j in 1:n
		@printf("%10s", labels[j])
	end
	println()
	println(repeat("-", 5 + 10*n)) # 打印分隔线

	# 打印矩阵内容
	for i in 1:n
		@printf("%-5s", labels[i] * "|") # 打印行标题和分隔符
		for j in 1:n
			@printf("%10.6f", matrix[i, j]) # 格式化浮点数
		end
		println()
	end
end



Hcore=T+V
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
P=[2*sum(C[i, m]*C[j, m] for m in 1:(ENum÷2)) for i in 1:Num, j in 1:Num]

for i in 1:MaxIter
	global P, F, Fprime, E, Cprime, C, G, Pnew, p
	G=[sum(P[k, l]*(ERI[i, j, k, l]-0.5*ERI[i, l, k, j]) for k in 1:Num, l in 1:Num) for i in 1:Num, j in 1:Num]
	F=Hcore+G
	Fprime=X'*F*X
	E, Cprime=eigen(Fprime)
	p = sortperm(E)
	E=E[p]
	Cprime=Cprime[:, p]
	C=X*Cprime
	Pnew=[2*sum(C[i, m]*C[j, m] for m in 1:(ENum÷2)) for i in 1:Num, j in 1:Num]


	if sum(P-Pnew)<Threshold
		println("22222")
		P=Pnew
		break
	end
	P=Pnew
end

VNN = sum(Molecule[i].Z * Molecule[j].Z / sqrt(sum((Molecule[i].position .- Molecule[j].position) .^ 2)) for i in 1:NNum for j in (i+1):NNum)


Ee = 0.5 * sum(P .* (Hcore + F))
Etot = Ee + VNN


Te = sum(P .* T)
VNe = sum(P .* V)
Vee = 0.5 * sum(P .* G)
Etot_components = Te + VNe + Vee + VNN


@printf("\n--- Molecular Structure ---\n")
for i in 1:NNum
	atom = Molecule[i]
	@printf("Atom %d: %s (Z=%d) at (%.3f, %.3f, %.3f) Å\n", i, atom.symbol, atom.Z, atom.position[1]*0.52918, atom.position[2]*0.52918, atom.position[3]*0.52918)
end



@printf("\n--- Final Energy Results ---\n")
@printf("Total Energy      = %.10f Hartree\n", Etot)
@printf("Electronic Energy = %.10f Hartree\n", Ee)
@printf("Nuclear Repulsion =  %.10f Hartree\n", VNN)
=#