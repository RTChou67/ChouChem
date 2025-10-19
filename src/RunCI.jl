module RunCI

using ..Definitions
include("GetBasisList.jl")
include("CalcS.jl")
include("CalcT.jl")
include("CalcV.jl")
include("CalcG.jl")
include("UHF.jl")
include("CI.jl")

export main

function main(MolInAng::Vector{Atom}, Charge::Int, Multiplicity::Int, MaxExcitation::Int)
	TStart=time_ns()


	Bohr2Ang = 0.52917721092
	Molecule = [Atom(atom.symbol, atom.Z, atom.basis_set, atom.position ./ Bohr2Ang) for atom in MolInAng]

	@printf("\n--- Molecular Structure ---\n")
	for atom in MolInAng
		@printf("Atom: %-2s at (%8.4f, %8.4f, %8.4f) Å\n", atom.symbol, atom.position...)
	end
	println("---------------------------\n")

	SCF_Results=UHF.SCF(Molecule, Charge, Multiplicity, MaxIter = 128, Threshold = 1e-8)
	if isnothing(SCF_Results)
		error("UHF calculation did not converge. Aborting.")
		return
	else
		CI_Results=CI.RunCI(SCF_Results, MaxExcitation)
		if isnothing(CI_Results)
			error("CI calculation did not converge. Aborting.")
			return
		end
	end

	println("\n--- Post-CI Analysis ---")


	println("\nGround State Wavefunction Analysis (Top 5 components):")
	GroundStateVector = CI_Results.CIVectors[:, 1]

	sorted_indices = sortperm(abs.(GroundStateVector), rev = true)

	RefDeterminant = CI_Results.DeterminantSpace[1]
	@printf("  Coeff    | Contribution  | Determinant Configuration\n")
	println("  -----------------------------------------------------")
	for i in 1:min(5, CI_Results.NumberOfDeterminants)
		idx = sorted_indices[i]
		coeff = GroundStateVector[idx]
		determinant = CI_Results.DeterminantSpace[idx]
		is_ref = (determinant == RefDeterminant) ? "(Ref)" : ""
		@printf("  %-8.4f | %-12.2f%% | %s %s\n", coeff, coeff^2 * 100, determinant, is_ref)
	end
	println("  -----------------------------------------------------\n")



	TEnd=time_ns()
	TSeconds = (TEnd - TStart) / 1e9
	days = floor(Int, TSeconds / 86400)
	hours = floor(Int, (TSeconds % 86400) / 3600)
	minutes = floor(Int, (TSeconds % 3600) / 60)
	seconds = TSeconds % 60
	DateTime = Dates.format(now(), "e u dd HH:MM:SS yyyy")
	@printf(" Job cpu time:       %d days %2d hours %2d minutes %5.1f seconds.\n", days, hours, minutes, seconds)
	@printf(" Elapsed time:       %d days %2d hours %2d minutes %5.1f seconds.\n", days, hours, minutes, seconds)
	println(" Normal termination of Julia UHF-CI at $(DateTime).")



end



end