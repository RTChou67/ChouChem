module RunRMPn
using ..Definitions

include("GetBasisList.jl")
include("CalcS.jl")
include("CalcT.jl")
include("CalcV.jl")
include("CalcG.jl")
include("RHF.jl")
include("RMPn.jl")

using .RMPn

using Dates

export main


function main(MolInAng::Vector{Atom}, Charge::Int, order::Int)
	TStart=time_ns()


	Bohr2Ang = 0.52917721092
	Molecule = [Atom(atom.symbol, atom.Z, atom.basis_set, atom.position ./ Bohr2Ang) for atom in MolInAng]

	RHF_Results = RHF.SCF(Molecule, Charge; MaxIter = 100, Threshold = 1e-8)

	if RHF_Results !== nothing
		println("RHF Calculation Successful. Proceeding to MPn.")

		MPn_Results = RMPn.RunRMPn(RHF_Results, order)
		println("\n--- Calculation Summary ---")
		println("RHF Total Energy: ", RHF_Results.Etot, " Hartree")
		if MPn_Results !== nothing
			println("MP2 Total Energy: ", MPn_Results.E_RMP2_total, " Hartree")
			TEnd=time_ns()
			TSeconds = (TEnd - TStart) / 1e9
			days = floor(Int, TSeconds / 86400)
			hours = floor(Int, (TSeconds % 86400) / 3600)
			minutes = floor(Int, (TSeconds % 3600) / 60)
			seconds = TSeconds % 60
			DateTime = Dates.format(now(), "e u dd HH:MM:SS yyyy")
			@printf(" Job cpu time:       %d days %2d hours %2d minutes %5.1f seconds.\n", days, hours, minutes, seconds)
			@printf(" Elapsed time:       %d days %2d hours %2d minutes %5.1f seconds.\n", days, hours, minutes, seconds)
			println(" Normal termination of Julia RMPn at $(DateTime).")
		end
		println("-------------------------\n")

	else
		println("RHF calculation failed to converge. Cannot perform MPn calculation.")
	end
end

end