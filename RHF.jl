# RHF.jl

module RHF

using LinearAlgebra
using Printf
using ..Definitions: Atom, CGTF, Basis
using ..GetBasisList: generate_basis_list
using ..CalcS: Sij
using ..CalcT: Tij
using ..CalcV: Vij
using ..CalcG: Gijkl

export SCF, RHFResults

module RHF_DIIS
	using LinearAlgebra
	export DIISController, insert!, extrapolate

	struct DIISController
		DIISMax::Int
		FockList::Vector{Matrix{Float64}}
		ErrList::Vector{Matrix{Float64}}

		function DIISController(MaxSize::Int = 8)
			new(MaxSize, [], [])
		end
	end

	function insert!(DC::DIISController, Fock::Matrix{Float64}, ErrMat::Matrix{Float64})
		push!(DC.FockList, Fock)
		push!(DC.ErrList, ErrMat)
		if length(DC.ErrList) > DC.DIISMax
			popfirst!(DC.FockList)
			popfirst!(DC.ErrList)
		end
	end

	function extrapolate(DC::DIISController)
		n = length(DC.ErrList)
		if n < 2
			return DC.FockList[end]
		end
		B = [dot(DC.ErrList[i], DC.ErrList[j]) for i in 1:n, j in 1:n]
		A = -ones(Float64, n + 1, n + 1)
		A[1:n, 1:n] = B
		A[n+1, n+1] = 0.0

		b = zeros(Float64, n + 1)
		b[n+1] = -1.0

		coeffs = pinv(A) * b
		c = coeffs[1:n]

		FockNext = sum(c[i] * DC.FockList[i] for i in 1:n)

		return FockNext
	end

end

struct RHFResults
	Molecule::Vector{Atom}
	BasisSet::Vector{Basis}
	ENum::Int
	S::Matrix{Float64}
	Hcore::Matrix{Float64}
	ERI::Array{Float64, 4}
	C::Matrix{Float64}
	E::Vector{Float64}
	P::Matrix{Float64}
	Ee::Float64
	VNN::Float64
	Etot::Float64
end

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

function SCF(Molecule::Vector{Atom}, charge::Int; MaxIter = 100, Threshold = 1e-10)
	BasisSet = generate_basis_list(Molecule)
	BNum = length(BasisSet)
	ENum = sum(atom.Z for atom in Molecule) - charge
	Nocc = ENum ÷ 2

	println("--- System Information ---")
	@printf("Basis functions: %d\n", BNum)
	@printf("Electrons:       %d\n", ENum)
	println("--------------------------\n")

	S, T, V, ERI = CalcMatrices(BasisSet, Molecule)
	Hcore = T + V
	X = S^(-0.5)

	E_guess, C_guess = eigen(X' * Hcore * X)
	p = sortperm(E_guess)
	C = X * C_guess[:, p]
	P = 2 * C[:, 1:Nocc] * C[:, 1:Nocc]'

	DIIS = RHF_DIIS.DIISController()

	println("--- Starting SCF Iterations (with DIIS) ---")
	Etot_old = 0.0

	for i in 1:MaxIter
		G = [sum(P[k, l] * (ERI[i, j, k, l] - 0.5 * ERI[i, l, k, j]) for k in 1:BNum, l in 1:BNum) for i in 1:BNum, j in 1:BNum]
		F_current = Hcore + G

		ErrMat = X' * (F_current * P * S - S * P * F_current) * X

		RHF_DIIS.insert!(DIIS, F_current, ErrMat)

		F = RHF_DIIS.extrapolate(DIIS)

		Fprime = X' * F * X
		E, Cprime = eigen(Fprime)

		p = sortperm(E)
		E = E[p]
		C = X * Cprime[:, p]

		Pnew = 2 * C[:, 1:Nocc] * C[:, 1:Nocc]'

		VNN = sum(Molecule[i].Z * Molecule[j].Z / norm(Molecule[i].position .- Molecule[j].position) for i in 1:length(Molecule) for j in (i+1):length(Molecule))
		Ee = 0.5 * sum(P .* (Hcore + F))
		Etot = Ee + VNN

		delta_E = abs(Etot - Etot_old)
		delta_P = sqrt(sum((Pnew - P) .^ 2))
		@printf("Iteration %3d: E = %-16.10f  ΔE = %-12.2e  ΔP = %.2e\n", i, Etot, delta_E, delta_P)

		P = Pnew
		Etot_old = Etot

		if delta_E < Threshold && delta_P < Threshold
			println("\nSCF converged in $i iterations.")
			@printf("\n--- Final Energy Results ---\n")
			@printf("Electronic Energy = %.10f Hartree\n", Ee)
			@printf("Nuclear Repulsion = %.10f Hartree\n", VNN)
			@printf("Total Energy      = %.10f Hartree\n", Etot)
			println("----------------------------\n")
			return RHFResults(Molecule, BasisSet, ENum, S, Hcore, ERI, C, E, P, Ee, VNN, Etot)
		end
	end

	println("\nSCF failed to converge after $MaxIter iterations.")
	return nothing
end



end
