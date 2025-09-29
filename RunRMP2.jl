include("Definitions.jl")
include("GetBasisList.jl")
include("CalcS.jl")
include("CalcT.jl")
include("CalcV.jl")
include("CalcG.jl")
include("RHF.jl")
include("RMP2.jl")

using .RMP2
using .Definitions


function main()

	MolInAng = [
		Atom("H", 1, "STO-3G", (0.0, 0.0, 0.0)),
		Atom("F", 9, "STO-3G", (0.0, 0.0, 1.0)),
	]
	Charge = 0
	Multiplicity = 1
	Bohr2Ang = 0.52917721092
	Molecule = [Atom(atom.symbol, atom.Z, atom.basis_set, atom.position ./ Bohr2Ang) for atom in MolInAng]

	RHF_Results = RHF.SCF(Molecule, Charge; MaxIter = 100, Threshold = 1e-10)

	if RHF_Results !== nothing
		println("RHF Calculation Successful. Proceeding to MP2.")

		MP2_Results = RMP2.RunRMP2(RHF_Results)

		println("--- Final Summary ---")
		@printf("RHF Total Energy:   %.10f Hartree\n", RHF_Results.Etot)
		@printf("MP2 Total Energy:   %.10f Hartree\n", MP2_Results.MP2_total)
		println("---------------------\n")

	else
		println("RHF calculation failed to converge. Cannot perform MPn calculation.")
	end
end

main()