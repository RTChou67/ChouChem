using LinearAlgebra
using Printf
using Dates

include("Definitions.jl")
include("GetBasisList.jl")
include("CalcS.jl")
include("CalcT.jl")
include("CalcV.jl")
include("CalcG.jl")
using .Definitions: PGTF, CGTF, Basis, Atom
using .GetBasisList: generate_basis_list
using .CalcS: Sij
using .CalcT: Tij
using .CalcV: Vij
using .CalcG: Gijkl

module UHF_DIIS
using LinearAlgebra
export DIISController, insert!, extrapolate

struct DIISController
	DIISMax::Int
	FockListAlpha::Vector{Matrix{Float64}}
	FockListBeta::Vector{Matrix{Float64}}
	ErrList::Vector{Matrix{Float64}}

	function DIISController(MaxSize::Int = 8)
		new(MaxSize, [], [], [])
	end
end

function insert!(DC::DIISController, FockAlpha::Matrix{Float64}, FockBeta::Matrix{Float64}, ErrMat::Matrix{Float64})
	push!(DC.FockListAlpha, FockAlpha)
	push!(DC.FockListBeta, FockBeta)
	push!(DC.ErrList, ErrMat)
	if length(DC.ErrList) > DC.DIISMax
		popfirst!(DC.FockListAlpha)
		popfirst!(DC.FockListBeta)
		popfirst!(DC.ErrList)
	end
end

function extrapolate(DC::DIISController)
	n = length(DC.ErrList)
	if n < 2
		return DC.FockListAlpha[end], DC.FockListBeta[end]
	end

	B = [dot(DC.ErrList[i], DC.ErrList[j]) for i in 1:n, j in 1:n]

	A = ones(Float64, n + 1, n + 1) * -1.0
	A[1:n, 1:n] = B
	A[n+1, n+1] = 0.0

	b = zeros(Float64, n + 1)
	b[n+1] = -1.0

	coeffs = pinv(A) * b
	c = coeffs[1:n]

	FockAlphaNext = sum(c[i] * DC.FockListAlpha[i] for i in 1:n)
	FockBetaNext = sum(c[i] * DC.FockListBeta[i] for i in 1:n)

	return FockAlphaNext, FockBetaNext
end

end


module UHF
using LinearAlgebra
using Printf
using ..Definitions: Atom
using ..UHF_DIIS
using ..GetBasisList: generate_basis_list
using ..CalcS: Sij
using ..CalcT: Tij
using ..CalcV: Vij
using ..CalcG: Gijkl

export SCF

function CalcMatrices(BasisSet, Molecule)
	BNum = length(BasisSet)
	TimeS1=time_ns()
	S = [Sij(BasisSet[i], BasisSet[j]) for i in 1:BNum, j in 1:BNum]
	TimeS2=time_ns()
	println("Calculation for S Matrix took $( (TimeS2-TimeS1)/1e6 ) ms")
	TimeT1=time_ns()
	T = [Tij(BasisSet[i], BasisSet[j]) for i in 1:BNum, j in 1:BNum]
	TimeT2=time_ns()
	println("Calculation for T Matrix took $( (TimeT2-TimeT1)/1e6 ) ms")
	TimeV1=time_ns()
	V = [Vij(BasisSet[i], BasisSet[j], Molecule) for i in 1:BNum, j in 1:BNum]
	TimeV2=time_ns()
	println("Calculation for V Matrix took $( (TimeV2-TimeV1)/1e6 ) ms")
	TimeERI1=time_ns()
	ERI = [Gijkl(BasisSet[i], BasisSet[j], BasisSet[k], BasisSet[l]) for i in 1:BNum, j in 1:BNum, k in 1:BNum, l in 1:BNum]
	TimeERI2=time_ns()
	println("Calculation for ERI Tensor took $( (TimeERI2-TimeERI1)/1e6 ) ms")
	return S, T, V, ERI
end

function SCF(Molecule::Vector{Atom}, charge::Int, multiplicity::Int; MaxIter = 100, Threshold = 1e-10)
	BasisSet = generate_basis_list(Molecule)
	BNum = length(BasisSet)
	ENum = sum(atom.Z for atom in Molecule) - charge
	ENumAlpha = (ENum + multiplicity - 1) ÷ 2
	ENumBeta = (ENum - multiplicity + 1) ÷ 2

	println("--- System Information ---")
	@printf("Basis functions: %d\n", BNum)
	@printf("Alpha electrons: %d\n", ENumAlpha)
	@printf("Beta electrons:  %d\n", ENumBeta)
	println("--------------------------\n")

	S, T, V, ERI = CalcMatrices(BasisSet, Molecule)
	Hcore = T + V
	X = S^(-0.5)
	EGuess, CGuess = eigen(X' * Hcore * X)
	C = X * CGuess[:, sortperm(EGuess)]
	PAlpha = C[:, 1:ENumAlpha] * C[:, 1:ENumAlpha]'
	PBeta = C[:, 1:ENumBeta] * C[:, 1:ENumBeta]'

	DIIS = UHF_DIIS.DIISController()

	println("--- Starting SCF Iterations ---")
	Etot_old = 0.0

	for i in 1:MaxIter
		Ptotal = PAlpha + PBeta
		GAlpha = [sum(Ptotal[k, l]*ERI[i, j, k, l] - PAlpha[k, l]*ERI[i, l, k, j] for k in 1:BNum, l in 1:BNum) for i in 1:BNum, j in 1:BNum]
		GBeta = [sum(Ptotal[k, l]*ERI[i, j, k, l] - PBeta[k, l]*ERI[i, l, k, j] for k in 1:BNum, l in 1:BNum) for i in 1:BNum, j in 1:BNum]
		FAlpha_current = Hcore + GAlpha
		FBeta_current = Hcore + GBeta

		ErrMat = X' * (FAlpha_current * PAlpha * S - S * PAlpha * FAlpha_current) * X +
				 X' * (FBeta_current * PBeta * S - S * PBeta * FBeta_current) * X

		UHF_DIIS.insert!(DIIS, FAlpha_current, FBeta_current, ErrMat)

		FAlpha, FBeta = UHF_DIIS.extrapolate(DIIS)

		EAlphaVec, CprimeAlpha = eigen(X' * FAlpha * X)
		EBetaVec, CprimeBeta = eigen(X' * FBeta * X)

		CAlpha = X * CprimeAlpha[:, sortperm(EAlphaVec)]
		CBeta = X * CprimeBeta[:, sortperm(EBetaVec)]

		PnewAlpha = CAlpha[:, 1:ENumAlpha] * CAlpha[:, 1:ENumAlpha]'
		PnewBeta = CBeta[:, 1:ENumBeta] * CBeta[:, 1:ENumBeta]'

		VNN = sum(Molecule[i].Z * Molecule[j].Z / norm(Molecule[i].position .- Molecule[j].position) for i in 1:length(Molecule) for j in (i+1):length(Molecule))
		Ee = 0.5 * sum(PAlpha .* (Hcore + FAlpha)) + 0.5 * sum(PBeta .* (Hcore + FBeta))
		Etot = Ee + VNN

		delta_E = abs(Etot - Etot_old)
		delta_P = max(sqrt(sum((PnewAlpha - PAlpha) .^ 2)), sqrt(sum((PnewBeta - PBeta) .^ 2)))
		@printf("Iteration %3d: E = %-16.10f  ΔE = %-12.2e  ΔP = %.2e\n", i, Etot, delta_E, delta_P)

		PAlpha = PnewAlpha
		PBeta = PnewBeta
		Etot_old = Etot

		if delta_E < Threshold && delta_P < Threshold
			println("\nSCF converged in $i iterations.")
			@printf("\n--- Final Energy Results ---\n")
			@printf("Electronic Energy = %.10f Hartree\n", Ee)
			@printf("Nuclear Repulsion = %.10f Hartree\n", VNN)
			@printf("Total Energy      = %.10f Hartree\n", Etot)
			println("----------------------------\n")
			return Etot
		end
	end
	println("\nSCF failed to converge after $MaxIter iterations.")
	return nothing
end

end

function main()
	TStart=time_ns()
	MolInAng = [
		Atom("H", 1, "6-31G", (0.0, 0.0, -1.0)),
		Atom("F", 9, "6-31G", (0.0, 0.0, 1.0)),
	]
	Charge = 0
	Multiplicity = 1

	Bohr2Ang = 0.52917721092
	Molecule = [Atom(atom.symbol, atom.Z, atom.basis_set, atom.position ./ Bohr2Ang) for atom in MolInAng]

	@printf("\n--- Molecular Structure ---\n")
	for atom in MolInAng
		@printf("Atom: %-2s at (%8.4f, %8.4f, %8.4f) Å\n", atom.symbol, atom.position...)
	end
	println("---------------------------\n")

	UHF.SCF(Molecule, Charge, Multiplicity, MaxIter = 128, Threshold = 1e-8)
	TEnd=time_ns()
	TSeconds = (TEnd - TStart) / 1e9
	days = floor(Int, TSeconds / 86400)
	hours = floor(Int, (TSeconds % 86400) / 3600)
	minutes = floor(Int, (TSeconds % 3600) / 60)
	seconds = TSeconds % 60
	DateTime = Dates.format(now(), "e u dd HH:MM:SS yyyy")
	@printf(" Job cpu time:       %d days %2d hours %2d minutes %5.1f seconds.\n", days, hours, minutes, seconds)
	@printf(" Elapsed time:       %d days %2d hours %2d minutes %5.1f seconds.\n", days, hours, minutes, seconds)
	println(" Normal termination of Julia HF at $(DateTime).")
end

main()
