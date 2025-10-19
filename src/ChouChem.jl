module ChouChem

using LinearAlgebra
using Printf
using Combinatorics
using Dates

include("Definitions.jl")
include("GetBasisList.jl")
include("CalcS.jl")
include("CalcT.jl")
include("CalcV.jl")
include("CalcG.jl")

include("RHF.jl")
include("RunRHF.jl")

include("UHF.jl")
include("RunUHF.jl")

include("RMPn.jl")
include("RunRMPn.jl")

include("CI.jl")
include("RunCI.jl")

run_rhf(MolInAng::Vector{Atom}, Charge::Int) = RunRHF.main(MolInAng, Charge)
run_uhf(MolInAng::Vector{Atom}, Charge::Int, Multiplicity::Int) = RunUHF.main(MolInAng, Charge, Multiplicity)
run_rmpn(MolInAng::Vector{Atom}, Charge::Int, order::Int) = RunRMPn.main(MolInAng, Charge, order)
run_ci(MolInAng::Vector{Atom}, Charge::Int, Multiplicity::Int, MaxExcitation::Int) = RunCI.main(MolInAng, Charge, Multiplicity, MaxExcitation)

export Atom
export run_rhf
export run_uhf
export run_rmpn
export run_ci


end