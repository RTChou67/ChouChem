using LinearAlgebra
using Printf
using SpecialFunctions
using Combinatorics



include("Definitions.jl")
include("GetBasisList.jl")
include("CalcS.jl")
include("CalcT.jl")
include("CalcV.jl")
include("CalcG.jl")
using .Definitions: PGTF, CGTF, Basis, Atom
using .GetBasisList: generate_basis_list, get_basis_set
using .CalcS: Sij
using .CalcT: Tij
using .CalcV: Vij
using .CalcG: Gijkl



MolInAng = [
	Atom("H", 1, "STO-3G", (0.0, 0.0, +1.0)),
	Atom("F", 9, "STO-3G", (0.0, 0.0, -1.0)),
]
Charge=0
Multiplicity=1

Bohr2Ang=0.52917721092
MolInBohr = [Atom(atom.symbol, atom.Z, atom.basis_set, atom.position ./ Bohr2Ang) for atom in MolInAng]
Molecule = MolInBohr

BasisSet = generate_basis_list(Molecule)



ONum=length(BasisSet)
ENum=sum(atom.Z for atom in Molecule)-Charge
NNum=length(Molecule)

ENumAlpha=(ENum+Multiplicity-1)÷2
ENumBeta=(ENum-Multiplicity+1)÷2





S=[Sij(BasisSet[i], BasisSet[j]) for i in 1:ONum, j in 1:ONum]

T=[Tij(BasisSet[i], BasisSet[j]) for i in 1:ONum, j in 1:ONum]

V=[Vij(BasisSet[i], BasisSet[j], Molecule) for i in 1:ONum, j in 1:ONum]

ERI_AO=[Gijkl(BasisSet[i], BasisSet[j], BasisSet[k], BasisSet[l]) for i in 1:ONum, j in 1:ONum, k in 1:ONum, l in 1:ONum]



function format_to_custom_eng(x::Float64)
	if x == 0.0
		return "0.000000E+00"
	end
	std_eng_str = @sprintf("%1.6e", x)
	parts = split(std_eng_str, 'e')
	mantissa_std = parse(Float64, parts[1])
	exponent_std = parse(Int, parts[2])
	mantissa_custom = mantissa_std / 10.0
	exponent_custom = exponent_std + 1
	sign_char = mantissa_std < 0 ? "-" : ""
	mantissa_str = @sprintf("%.6f", abs(mantissa_custom))
	exponent_str = @sprintf("%+03d", exponent_custom)
	return string(sign_char, mantissa_str, "E", exponent_str)
end

function print_formatted_matrix(matrix::Matrix{Float64})
	n = size(matrix, 1)
	labels = [string(i) for i in 1:n]
	@printf("%-5s", "")
	for j in 1:n
		@printf("%15s", labels[j])
	end
	println()
	println(repeat("-", 5 + 15*n))
	for i in 1:n
		@printf("%-5s", labels[i] * "|")
		for j in 1:n
			formatted_num = format_to_custom_eng(matrix[i, j])
			@printf("%15s", formatted_num)
		end
		println()
	end
end


println("\nOverlap Matrix S:")
print_formatted_matrix(S)
println("\nKinetic Energy Matrix T:")
print_formatted_matrix(T)
println("\nNuclear Attraction Matrix V:")
print_formatted_matrix(V)




Hcore=T+V

println("\nCore Hamiltonian Hcore:")
print_formatted_matrix(Hcore)


F=Hcore

MaxIter=100
Threshold=1e-10



X=S^(-0.5)
Fprime=X'*F*X
E, Cprime=eigen(Fprime)
C=X*Cprime
p = sortperm(E)
E=E[p]
Cprime=Cprime[:, p]
PAlpha=[sum(C[i, m]*C[j, m] for m in 1:ENumAlpha) for i in 1:ONum, j in 1:ONum]
PBeta=[sum(C[i, m]*C[j, m] for m in 1:ENumBeta) for i in 1:ONum, j in 1:ONum]



for i in 1:MaxIter
	global PAlpha, PBeta, FAlpha, FBeta, FprimeAlpha, FprimeBeta, EAlpha, EBeta, CprimeAlpha, CprimeBeta, CAlpha, CBeta, GAlpha, GBeta, PnewAlpha, PnewBeta
	Ptotal=PAlpha+PBeta
	GAlpha = [sum(Ptotal[k, l]*ERI_AO[i, j, k, l] - PAlpha[k, l]*ERI_AO[i, l, k, j] for k in 1:ONum, l in 1:ONum) for i in 1:ONum, j in 1:ONum]
	GBeta = [sum(Ptotal[k, l]*ERI_AO[i, j, k, l] - PBeta[k, l]*ERI_AO[i, l, k, j] for k in 1:ONum, l in 1:ONum) for i in 1:ONum, j in 1:ONum]
	FAlpha=Hcore+GAlpha
	FBeta=Hcore+GBeta
	FprimeAlpha=X'*FAlpha*X
	FprimeBeta=X'*FBeta*X
	EAlpha, CprimeAlpha=eigen(FprimeAlpha)
	EBeta, CprimeBeta=eigen(FprimeBeta)
	pAlpha=sortperm(EAlpha)
	pBeta=sortperm(EBeta)
	EAlpha=EAlpha[pAlpha]
	EBeta=EBeta[pBeta]
	CprimeAlpha=CprimeAlpha[:, pAlpha]
	CprimeBeta=CprimeBeta[:, pBeta]
	CAlpha=X*CprimeAlpha
	CBeta=X*CprimeBeta
	PnewAlpha=[sum(CAlpha[i, m]*CAlpha[j, m] for m in 1:ENumAlpha) for i in 1:ONum, j in 1:ONum]
	PnewBeta=[sum(CBeta[i, m]*CBeta[j, m] for m in 1:ENumBeta) for i in 1:ONum, j in 1:ONum]
	delta_PAlpha = sqrt(sum((PnewAlpha - PAlpha) .^ 2))
	delta_PBeta = sqrt(sum((PnewBeta - PBeta) .^ 2))
	delta_P = max(delta_PAlpha, delta_PBeta)
	println("Iteration $i: ΔP = $delta_P")
	if delta_P < Threshold
		println("UHF converged in $i iterations with ΔP = $delta_P")
		PAlpha=PnewAlpha
		PBeta=PnewBeta
		break
	end
	PAlpha=PnewAlpha
	PBeta=PnewBeta
end

Ee_UHF = 0.5 * (sum(PAlpha .* (Hcore + FAlpha)) + sum(PBeta .* (Hcore + FBeta)))
VNN_UHF = sum(Molecule[i].Z * Molecule[j].Z / sqrt(sum((Molecule[i].position .- Molecule[j].position) .^ 2)) for i in 1:NNum for j in (i+1):NNum)
Etot_UHF = Ee_UHF + VNN_UHF



@printf("\n--- Molecular Structure ---\n")
for i in 1:NNum
	atom = Molecule[i]
	@printf("Atom %d: %s (Z=%d) at (%+.6f, %+.6f, %+.6f) Å\n", i, atom.symbol, atom.Z, atom.position[1]*0.52918, atom.position[2]*0.52918, atom.position[3]*0.52918)
end



@printf("\n--- Final Energy Results (UHF) ---\n")
@printf("Total Energy      = %.10f Hartree\n", Etot_UHF)
@printf("Electronic Energy = %.10f Hartree\n", Ee_UHF)
@printf("Nuclear Repulsion =  %.10f Hartree\n", VNN_UHF)



SONum=2*ONum

HcoreAlphaMO=CAlpha'*Hcore*CAlpha
HcoreBetaMO=CBeta'*Hcore*CBeta

hMO=[HcoreAlphaMO zeros(ONum, ONum); zeros(ONum, ONum) HcoreBetaMO]

function TransERI(ERI_AO::Array{Float64, 4}, c1::Matrix{Float64}, c2::Matrix{Float64}, c3::Matrix{Float64}, c4::Matrix{Float64})
	return [sum(c1[i, p] * c2[j, q] * c3[k, r] * c4[l, s] * ERI_AO[i, j, k, l] for i in 1:ONum, j in 1:ONum, k in 1:ONum, l in 1:ONum) for p in 1:ONum, q in 1:ONum, r in 1:ONum, s in 1:ONum]
end
println("\nTransforming ERIs to MO basis...")
ERI_aaaa=TransERI(ERI_AO, CAlpha, CAlpha, CAlpha, CAlpha)
ERI_bbbb=TransERI(ERI_AO, CBeta, CBeta, CBeta, CBeta)
ERI_aabb=TransERI(ERI_AO, CAlpha, CAlpha, CBeta, CBeta)
ERI_abab=TransERI(ERI_AO, CAlpha, CBeta, CAlpha, CBeta)
println("ERI transformation complete.")

function GetERI(p::Int, r::Int, q::Int, s::Int, ONum::Int, aaaa, bbbb, aabb)
	IsPA = p <= ONum
	IsRA = r <= ONum
	IsQA = q <= ONum
	IsSA = s <= ONum

	PSpa = IsPA ? p : p - ONum
	RSpa = IsRA ? r : r - ONum
	QSpa = IsQA ? q : q - ONum
	SSpa = IsSA ? s : s - ONum

	if IsPA && IsRA
		if IsQA && IsSA
			return aaaa[PSpa, RSpa, QSpa, SSpa]
		elseif !IsQA && !IsSA
			return aabb[PSpa, RSpa, QSpa, SSpa]
		else
			return 0.0
		end
	elseif !IsPA && !IsRA
		if !IsQA && !IsSA
			return bbbb[PSpa, RSpa, QSpa, SSpa]
		elseif IsQA && IsSA
			return aabb[QSpa, SSpa, PSpa, RSpa]
		else
			return 0.0
		end
	else
		return 0.0
	end
end




Determinants=collect(combinations(1:SONum, ENum))

function GetPhase(p::Int, q::Int, CommonOrb::Vector{Int})
	NInterOrbs = count(i -> (p < i < q) || (q < i < p), CommonOrb)
	return (-1)^NInterOrbs
end

function GetPhase(S1::Vector{Int}, p::Int, q::Int, r::Int, s::Int)
	PositionP = findfirst(==(p), S1)
	PositionQ = findfirst(==(q), S1)

	PhaseAnnihilate = (-1)^(PositionP - 1 + PositionQ - 2)

	S2 = sort(vcat(collect(setdiff(S1, [p, q])), [r, s]))

	PositionR = findfirst(==(r), S2)
	PositionS = findfirst(==(s), S2)

	PhaseCreate = (-1)^(PositionR - 1 + PositionS - 2)

	return PhaseAnnihilate * PhaseCreate
end




function CalcCI(S1::Vector{Int64}, S2::Vector{Int64}, hMO::Matrix{Float64}, ONum::Int, ERI_aaaa, ERI_bbbb, ERI_aabb)
	S1Set=Set(S1)
	S2Set=Set(S2)
	Diff1=sort(collect(setdiff(S1Set, S2Set)))
	Diff2=sort(collect(setdiff(S2Set, S1Set)))
	NDiff=length(Diff1)

	if NDiff > 2
		return 0.0
	elseif NDiff == 0
		E1e = sum(hMO[i, i] for i in S1)
		E2e = sum(GetERI(i, i, j, j, ONum, ERI_aaaa, ERI_bbbb, ERI_aabb) - GetERI(i, j, i, j, ONum, ERI_aaaa, ERI_bbbb, ERI_aabb) for (idx, i) in enumerate(S1) for j in S1[(idx+1):end])
		return E1e + E2e
	elseif NDiff == 1
		p=Diff1[1]
		q=Diff2[1]
		OrbCommon=sort(collect(intersect(S1Set, S2Set)))
		Phase=GetPhase(p, q, OrbCommon)
		E1e=hMO[p, q]
		E2e = sum(GetERI(p, q, i, i, ONum, ERI_aaaa, ERI_bbbb, ERI_aabb) - GetERI(p, i, q, i, ONum, ERI_aaaa, ERI_bbbb, ERI_aabb) for i in OrbCommon)
		return (E1e + E2e)*Phase
	elseif NDiff == 2
		p, q=Diff1[1], Diff1[2]
		r, s=Diff2[1], Diff2[2]
		OrbCommon=sort(collect(intersect(S1Set, S2Set)))
		Phase = GetPhase(S1, p, q, r, s)
		return (GetERI(p, r, q, s, ONum, ERI_aaaa, ERI_bbbb, ERI_aabb) - GetERI(p, s, q, r, ONum, ERI_aaaa, ERI_bbbb, ERI_aabb)) * Phase
	end
end


CI=[CalcCI(Determinants[i], Determinants[j], hMO, ONum, ERI_aaaa, ERI_bbbb, ERI_aabb) for i in eachindex(Determinants), j in eachindex(Determinants)]


E_CI, C_CI=eigen(CI)

Ee_FCI=E_CI[1]
VNN_FCI = sum(Molecule[i].Z * Molecule[j].Z / sqrt(sum((Molecule[i].position .- Molecule[j].position) .^ 2)) for i in 1:NNum for j in (i+1):NNum)
Etot_FCI = Ee_FCI + VNN_FCI

@printf("\n--- Final Energy Results (FCI) ---\n")
@printf("Total Energy      = %.10f Hartree\n", Etot_FCI)
@printf("Electronic Energy = %.10f Hartree\n", Ee_FCI)
@printf("Nuclear Repulsion =  %.10f Hartree\n", VNN_FCI)


