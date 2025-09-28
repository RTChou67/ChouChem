module HartreeFock

using LinearAlgebra
using Printf
using SpecialFunctions

# It is assumed that the following files are in the same directory.
include("Definitions.jl")
include("GetBasisList.jl")
include("CalcS.jl")
include("CalcT.jl")
include("CalcV.jl")
include("CalcG.jl")

using .Definitions: Atom
using .GetBasisList: generate_basis_list
using .CalcS: Sij
using .CalcT: Tij
using .CalcV: Vij
using .CalcG: Gijkl

export RHF, UHF

function format_to_custom_eng(x::Float64)
	iszero(x) && return "0.000000E+00"
	exponent = floor(Int, log10(abs(x))) + 1
	mantissa = x / 10.0^exponent
	return @sprintf("%1.6fE%+03d", mantissa, exponent)
end


function print_formatted_matrix(matrix::Matrix{Float64})
	n = size(matrix, 1)
	col_width = 15
	header = "     " * join(lpad.(1:n, col_width))
	println(header)
	println(repeat("-", 5 + col_width * n))
	for i in 1:n
		print(rpad("$i|", 5))
		row_str = join(lpad.(format_to_custom_eng.(matrix[i, :]), col_width))
		println(row_str)
	end
end






function RHF(MolInAng, Charge)
	println("--- Starting RHF Calculation ---")
	Bohr2Ang=0.52917721092
	MolInBohr = [Atom(atom.symbol, atom.Z, atom.basis_set, atom.position ./ Bohr2Ang) for atom in MolInAng]
	Molecule = MolInBohr
	BasisSet = generate_basis_list(Molecule)
	Num=length(BasisSet)
	ENum=sum(atom.Z for atom in Molecule)-Charge
	NNum=length(Molecule)
	println("Number of Atoms: $NNum")
	println("Number of Electrons: $ENum")
	println("Number of Basis Functions: $Num")

	S=[Sij(BasisSet[i], BasisSet[j]) for i in 1:Num, j in 1:Num]
	println("\nOverlap Matrix S:")
	print_formatted_matrix(S)

	T=[Tij(BasisSet[i], BasisSet[j]) for i in 1:Num, j in 1:Num]
	println("\nKinetic Energy Matrix T:")
	print_formatted_matrix(T)

	V=[Vij(BasisSet[i], BasisSet[j], Molecule) for i in 1:Num, j in 1:Num]
	println("\nNuclear Attraction Matrix V:")
	print_formatted_matrix(V)

	ERI=[Gijkl(BasisSet[i], BasisSet[j], BasisSet[k], BasisSet[l]) for i in 1:Num, j in 1:Num, k in 1:Num, l in 1:Num]

	HCore=T+V

	MaxIter=1000
	EConv=1e-10
	PConv=1e-10

	F=HCore
	X=S^(-0.5)
	Fprime=X'*F*X
	E, Cprime=eigen(Fprime)
	C=X*Cprime
	p=sortperm(E)
	E=E[p]

	Cprime=Cprime[:, p]
	P=2*C[:, 1:(ENum÷2)]*C[:, 1:(ENum÷2)]'
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
		delta_P = sqrt(sum((Pnew - P) .^ 2))
        
		if delta_P < PConv
			println("SCF")
			P=Pnew
			break
		end
		P=Pnew
	end
end







end
