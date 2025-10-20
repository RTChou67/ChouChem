module ChouChem

using LinearAlgebra
using Printf
using Combinatorics
using Dates
using Quadmath
using SpecialFunctions

include("Definitions.jl")
include("MathFunctions.jl")

include("CalcS.jl")
include("CalcT.jl")
include("CalcV.jl")
include("CalcG.jl")

include("GetBasisList.jl")

include("RHF.jl")
include("UHF.jl")
include("RMPn.jl")
#include("RunRMPn.jl")

#include("CI.jl")
#include("RunCI.jl")

run_rhf(MolInAng::Vector{Atom}, Charge::Int) = RunRHF(MolInAng, Charge)
run_uhf(MolInAng::Vector{Atom}, Charge::Int, Multiplicity::Int) = RunUHF(MolInAng, Charge, Multiplicity)
run_rmpn(MolInAng::Vector{Atom}, Charge::Int, order::Int) = RunRMPn(MolInAng, Charge, order)
#run_ci(MolInAng::Vector{Atom}, Charge::Int, Multiplicity::Int, MaxExcitation::Int) = RunCI.main(MolInAng, Charge, Multiplicity, MaxExcitation)

export RunRHF
export RunUHF
export RunRMPn

export Atom, Basis, CGTF, PGTF
export generate_basis_list
export Sij, Tij, Vij, Gijkl

end