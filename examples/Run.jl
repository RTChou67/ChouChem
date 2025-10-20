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


HFCoord=[
	("H", 1, (0.000, 0.000, 0.000)),
	("F", 9, (0.000, 0.000, 0.917)),
]

H2OCoord=[
	("O", 8, (0.0000, 0.0000, 0.0000)),
	("H", 1, (0.0000, 0.7572, 0.5865)),
	("H", 1, (0.0000, -0.7572, 0.5865)),
]

HFSTO3G=mkMol(HFCoord, "STO-3G")
HF631G=mkMol(HFCoord, "6-31G")
H2OSTO3G=mkMol(H2OCoord, "STO-3G")
H2O631G=mkMol(H2OCoord, "6-31G")


HFSTO3G_RHF=RunRHF(HFSTO3G, 0)
HFSTO3G_UHF=RunUHF(HFSTO3G, 0, 1)
HFSTO3G_RMP2=RunRMPn(HFSTO3G, 0, 2)
HFSTO3G_CI=RunCI(HFSTO3G, 0, 1, 2)

HF631G_RHF=RunRHF(HF631G, 0)
HF631G_UHF=RunUHF(HF631G, 0, 1)
HF631G_RMP2=RunRMPn(HF631G, 0, 2)
HF631G_CI=RunCI(HF631G, 0, 1, 2)

H2OSTO3G_RHF=RunRHF(H2OSTO3G, 0)
H2OSTO3G_UHF=RunUHF(H2OSTO3G, 0, 1)
H2OSTO3G_RMP2=RunRMPn(H2OSTO3G, 0, 2)
H2OSTO3G_CI=RunCI(H2OSTO3G, 0, 1, 2)

H2O631G_RHF=RunRHF(H2O631G, 0)
H2O631G_UHF=RunUHF(H2O631G, 0, 1)
H2O631G_RMP2=RunRMPn(H2O631G, 0, 2)
H2O631G_CI=RunCI(H2O631G, 0, 1, 2)