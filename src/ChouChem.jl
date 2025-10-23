module ChouChem

using LinearAlgebra
using Printf
using Combinatorics
using BenchmarkTools
using Dates
using SpecialFunctions
using Arpack
using SparseArrays

include("Definitions.jl")


include("CalcS.jl")
include("CalcT.jl")
include("CalcV.jl")
include("CalcG.jl")

include("GetBasisList.jl")

include("RHF.jl")
include("UHF.jl")
include("Functions.jl")
include("RMPn.jl")
include("CI.jl")


export RunRHF
export RunUHF
export RunRMPn
export RunUCI
export RunRCI

export Atom, Basis, CGTF, PGTF
export generate_basis_list
export Sij, Tij, Vij, Gijkl
export RHF_SCF, RHFResults
export UHF_SCF, UHFResults

end