include("Definitions.jl")

include("RunRHF.jl")
include("RunUHF.jl")
include("RunRMPn.jl")
include("RunCI.jl")

using Dates
using .Definitions
using .RunRMPn


MolInAng = [
	Atom("H", 1, "6-31G", (0.0, 0.0, -1.0)),
	Atom("F", 9, "6-31G", (0.0, 0.0, 1.0)),
]

Charge = 0
Multiplicity = 1
order = 2
MaxExcitation = 2

RunRHF.main(MolInAng, Charge)
RunUHF.main(MolInAng, Charge, Multiplicity)
RunRMPn.main(MolInAng, Charge, order)
RunCI.main(MolInAng, Charge, Multiplicity, MaxExcitation)