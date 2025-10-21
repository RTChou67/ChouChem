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
	("F", 9, (0.000, 0.000, 1.500)),
]

H2OCoord=[
	("O", 8, (0.0000, 0.0000, 0.0000)),
	("H", 1, (0.0000, 1.0000, 0.8000)),
	("H", 1, (0.0000, -1.0000, 0.8000)),
]

HFSTO3G=mkMol(HFCoord, "STO-3G")
HF631G=mkMol(HFCoord, "6-31G")
H2OSTO3G=mkMol(H2OCoord, "STO-3G")
H2O631G=mkMol(H2OCoord, "6-31G")


HFSTO3G_RHF=RunRHF(HFSTO3G, 0, 1)
if HFSTO3G_RHF==nothing
	error("RHF calculation failed for HF with STO-3G basis set")
end
HFSTO3G_UHF=RunUHF(HFSTO3G, 0, 1)
if HFSTO3G_UHF==nothing
	error("UHF calculation failed for HF with STO-3G basis set")
end
HFSTO3G_RMP2=RunRMPn(HFSTO3G, 0, 1, 2)
if HFSTO3G_RMP2==nothing
	error("RMP2 calculation failed for HF with STO-3G basis set")
end
HFSTO3G_UCI=RunUCI(HFSTO3G, 0, 1, 2)
if HFSTO3G_UCI==nothing
	error("UCI calculation failed for HF with STO-3G basis set")
end
HFSTO3G_RCI=RunRCI(HFSTO3G, 0, 1, 2)
if HFSTO3G_RCI==nothing
	error("RCI calculation failed for HF with STO-3G basis set")
end

HF631G_RHF=RunRHF(HF631G, 0, 1)
if HF631G_RHF==nothing
	error("RHF calculation failed for HF with 6-31G basis set")
end
HF631G_UHF=RunUHF(HF631G, 0, 1)
if HF631G_UHF==nothing
	error("UHF calculation failed for HF with 6-31G basis set")
end
HF631G_RMP2=RunRMPn(HF631G, 0, 1, 2)
if HF631G_RMP2==nothing
	error("RMP2 calculation failed for HF with 6-31G basis set")
end
HF631G_UCI=RunUCI(HF631G, 0, 1, 2)
if HF631G_UCI==nothing
	error("UCI calculation failed for HF with 6-31G basis set")
end
HF631G_RCI=RunRCI(HF631G, 0, 1, 2)
if HF631G_RCI==nothing
	error("RCI calculation failed for HF with 6-31G basis set")
end

H2OSTO3G_RHF=RunRHF(H2OSTO3G, 0, 1)
if H2OSTO3G_RHF==nothing
	error("RHF calculation failed for H2O with STO-3G basis set")
end
H2OSTO3G_UHF=RunUHF(H2OSTO3G, 0, 1)
if H2OSTO3G_UHF==nothing
	error("UHF calculation failed for H2O with STO-3G basis set")
end
H2OSTO3G_RMP2=RunRMPn(H2OSTO3G, 0, 1, 2)
if H2OSTO3G_RMP2==nothing
	error("RMP2 calculation failed for H2O with STO-3G basis set")
end
H2OSTO3G_UCI=RunUCI(H2OSTO3G, 0, 1, 2)
if H2OSTO3G_UCI==nothing
	error("UCI calculation failed for H2O with STO-3G basis set")
end
H2OSTO3G_RCI=RunRCI(H2OSTO3G, 0, 1, 2)
if H2OSTO3G_RCI==nothing
	error("RCI calculation failed for H2O with STO-3G basis set")
end

H2O631G_RHF=RunRHF(H2O631G, 0, 1)
if H2O631G_RHF==nothing
	error("RHF calculation failed for H2O with 6-31G basis set")
end
H2O631G_UHF=RunUHF(H2O631G, 0, 1)
if H2O631G_UHF==nothing
	error("UHF calculation failed for H2O with 6-31G basis set")
end
H2O631G_RMP2=RunRMPn(H2O631G, 0, 1, 2)
if H2O631G_RMP2==nothing
	error("RMP2 calculation failed for H2O with 6-31G basis set")
end
H2O631G_UCI=RunUCI(H2O631G, 0, 1, 2)
if H2O631G_UCI==nothing
	error("UCI calculation failed for H2O with 6-31G basis set")
end
H2O631G_RCI=RunRCI(H2O631G, 0, 1, 2)
if H2O631G_RCI==nothing
	error("RCI calculation failed for H2O with 6-31G basis set")
end
