using Pkg
Pkg.activate("..")

using ChouChem

function mkMol(Coord, BasisSet::String)
	Mol = Vector{Atom}()
	for atom in Coord
		push!(Mol, Atom(atom[1], atom[2], BasisSet, atom[3]))
	end
	return Mol
end


O2Coord=[
	("O", 8, (0.000, 0.000, 0.000)),
	("O", 8, (0.000, 0.000, 1.208)),
]


O2STO3G=mkMol(O2Coord, "STO-3G")
run_rmpn(O2STO3G, 0, 2)


#=
HFCoord=[
	("H", 1, (0.000, 0.000, 0.000)),
	("F", 9, (0.000, 0.000, 0.917)),
]

H2OCoord=[
	("O", 8, (0.0000, 0.0000, 0.0000)),
	("H", 1, (0.0000, 0.7572, 0.5865)),
	("H", 1, (0.0000, -0.7572, 0.5865)),
]

O2STO3G=mkMol(O2Coord, "STO-3G")
O2631G=mkMol(O2Coord, "6-31G")
HFSTO3G=mkMol(HFCoord, "STO-3G")
HF631G=mkMol(HFCoord, "6-31G")
H2OSTO3G=mkMol(H2OCoord, "STO-3G")
H2O631G=mkMol(H2OCoord, "6-31G")




O2STO3G_RHF=run_rhf(O2STO3G, 0)
O2631G_RHF=run_rhf(O2631G, 0)
HFSTO3G_RHF=run_rhf(HFSTO3G, 0)
HF631G_RHF=run_rhf(HF631G, 0)
H2OSTO3G_RHF=run_rhf(H2OSTO3G, 0)
H2O631G_RHF=run_rhf(H2O631G, 0)

O2STO3G_UHF=run_uhf(O2STO3G, 0, 1)
O2631G_UHF=run_uhf(O2631G, 0, 1)
HFSTO3G_UHF=run_uhf(HFSTO3G, 0, 1)
HF631G_UHF=run_uhf(HF631G, 0, 1)
H2OSTO3G_UHF=run_uhf(H2OSTO3G, 0, 1)
H2O631G_UHF=run_uhf(H2O631G, 0, 1)

O2STO3G_RMP2=run_rmpn(O2STO3G, 0, 2)
O2631G_RMP2=run_rmpn(O2631G, 0, 2)
HFSTO3G_RMP2=run_rmpn(HFSTO3G, 0, 2)
HF631G_RMP2=run_rmpn(HF631G, 0, 2)
H2OSTO3G_RMP2=run_rmpn(H2OSTO3G, 0, 2)
H2O631G_RMP2=run_rmpn(H2O631G, 0, 2)

O2STO3G_CISD=run_ci(O2STO3G, 0, 1, 2)
O2631G_CISD=run_ci(O2631G, 0, 1, 2)
HFSTO3G_CISD=run_ci(HFSTO3G, 0, 1, 2)
HF631G_CISD=run_ci(HF631G, 0, 1, 2)
H2OSTO3G_CISD=run_ci(H2OSTO3G, 0, 1, 2)
H2O631G_CISD=run_ci(H2O631G, 0, 1, 2)
=#