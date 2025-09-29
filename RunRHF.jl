
using Printf
using Dates

include("Definitions.jl")
include("GetBasisList.jl")
include("CalcS.jl")
include("CalcT.jl")
include("CalcV.jl")
include("CalcG.jl")
include("RHF.jl")

using .RHF
using .Definitions

function main()
    TStart = time_ns()

    MolInAng = [
        Atom("H", 1, "STO-3G", (0.0, 0.0, +1.0)),
        Atom("F", 9, "STO-3G", (0.0, 0.0, -1.0)),
    ]
    Charge = 0


    Bohr2Ang = 0.52917721092
    Molecule = [Atom(atom.symbol, atom.Z, atom.basis_set, atom.position ./ Bohr2Ang) for atom in MolInAng]

    @printf("\n--- Molecular Structure ---\n")
    for atom in MolInAng
        @printf("Atom: %-2s at (%8.4f, %8.4f, %8.4f) Å\n", atom.symbol, atom.position...)
    end
    println("---------------------------\n")

    SCF_Results = RHF.SCF(Molecule, Charge, MaxIter=100, Threshold=1e-10)
    
    if isnothing(SCF_Results)
        error("RHF calculation did not converge. Aborting.")
        return
    end

    TEnd = time_ns()
    TSeconds = (TEnd - TStart) / 1e9
    days = floor(Int, TSeconds / 86400)
    hours = floor(Int, (TSeconds % 86400) / 3600)
    minutes = floor(Int, (TSeconds % 3600) / 60)
    seconds = TSeconds % 60
    DateTime = Dates.format(now(), "e u dd HH:MM:SS yyyy")
    
    @printf(" Job cpu time:       %d days %2d hours %2d minutes %5.1f seconds.\n", days, hours, minutes, seconds)
    @printf(" Elapsed time:       %d days %2d hours %2d minutes %5.1f seconds.\n", days, hours, minutes, seconds)
    println(" Normal termination of Julia RHF at $(DateTime).")
end

main()