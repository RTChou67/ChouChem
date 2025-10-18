include("Definitions.jl")
include("GetBasisList.jl")
include("CalcS.jl")
include("CalcT.jl")
include("CalcV.jl")
include("CalcG.jl")
include("RHF.jl")
include("RMPn.jl")

using .RMPn
using .Definitions


function main()
	MolInAng = [
		Atom("H", 1, "STO-3G", (0.0, 0.0, -1.0)),
		Atom("F", 9, "STO-3G", (0.0, 0.0, +1.0)),
	]
	Charge = 0
	Multiplicity = 1
	order=2

	Bohr2Ang = 0.52917721092
	Molecule = [Atom(atom.symbol, atom.Z, atom.basis_set, atom.position ./ Bohr2Ang) for atom in MolInAng]

	RHF_Results = RHF.SCF(Molecule, Charge; MaxIter = 100, Threshold = 1e-10)

	if RHF_Results !== nothing
		println("RHF Calculation Successful. Proceeding to MPn.")

		MPn_Results = RMPn.RunRMPn(RHF_Results, order)
		println("\n--- Calculation Summary ---")
		println("RHF Total Energy: ", RHF_Results.Etot, " Hartree")
		if MPn_Results !== nothing
			println("MP2 Total Energy: ", MPn_Results.E_RMP2_total, " Hartree")
		end
		println("-------------------------\n")

	else
		println("RHF calculation failed to converge. Cannot perform MPn calculation.")
	end
end

main()