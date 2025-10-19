using Pkg
Pkg.activate("..")

using ChouChem

MolInAng = [
	Atom("O", 8, "STO-3G", (0.0, 0.0, 0.0)),
	Atom("O", 8, "STO-3G", (0.0, 0.0, 0.96)),
]

Charge = 0
Multiplicity = 1
order = 2
MaxExcitation = 2



ci_results = run_ci(MolInAng, Charge, Multiplicity, MaxExcitation)