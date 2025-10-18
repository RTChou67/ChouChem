module CI

using LinearAlgebra
using Printf
using Combinatorics


using ..UHF: UHFResults
using ..Definitions: Atom


export RunCI, CIResults

struct CIResults
	EtotUHF::Float64
	EtotCI::Float64
	Ecorr::Float64
	StateEnergies::Vector{Float64}

	CIVectors::Matrix{Float64}
	DeterminantSpace::Vector{Vector{Int}}

	MaxExcitation::Int
	NumberOfDeterminants::Int
end

function TransERI(ERI_AO::Array{Float64, 4}, c1::Matrix{Float64}, c2::Matrix{Float64}, c3::Matrix{Float64}, c4::Matrix{Float64}, ONum::Int)
	return [sum(c1[i, p] * c2[j, q] * c3[k, r] * c4[l, s] * ERI_AO[i, j, k, l] for i in 1:ONum, j in 1:ONum, k in 1:ONum, l in 1:ONum) for p in 1:ONum, q in 1:ONum, r in 1:ONum, s in 1:ONum]
end


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



function GenCISpace(RefDet::Vector{Int}, OrbRange::UnitRange{Int}, MaxExcit::Int)
	OccOrbs=RefDet
	VirOrbs=collect(setdiff(Set(collect(OrbRange)), Set(OccOrbs)))

	println("\n--- Excitation Analysis ---")
	DetAllLvl=Set([RefDet])
	@printf("Number of 0-fold excitations: %d\n", 1)
	for Lvl in 1:MaxExcit
		if Lvl>length(OccOrbs) || Lvl>length(VirOrbs)
			break
		end
		DetThisLvl=[]

		OccComb=combinations(OccOrbs, Lvl)
		VirComb=combinations(VirOrbs, Lvl)
		for Occ in OccComb
			BaseDetSet=setdiff(Set(OccOrbs), Set(Occ))
			for Vir in VirComb
				NewDetSet=union(BaseDetSet, Set(Vir))
				push!(DetThisLvl, sort(collect(NewDetSet)))
			end
		end
		UniqueDetListThisLvl=unique(DetThisLvl)
		CountThisLvl=length(UniqueDetListThisLvl)

		@printf("Number of %-2d-fold excitations: %d\n", Lvl, CountThisLvl)
		union!(DetAllLvl, UniqueDetListThisLvl)
	end
	println("--------------------------------------------\n")
	return collect(DetAllLvl)
end







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
	elseif NDiff == 1 # Single excitation difference
		p = Diff1[1]
		q = Diff2[1]
		PPos = findfirst(==(p), S1)
		QPos = findfirst(==(q), S2)
		Phase = (-1)^(PPos + QPos)

		OrbCommon=sort(collect(intersect(S1Set, S2Set)))
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

function RunCI(SCF_Results::UHFResults, MaxExcitation::Int)
	Molecule = SCF_Results.Molecule
	ONum = length(SCF_Results.BasisSet)
	ENumAlpha = SCF_Results.ENumAlpha
	ENumBeta = SCF_Results.ENumBeta
	Hcore = SCF_Results.Hcore
	CAlpha = SCF_Results.CAlpha
	CBeta = SCF_Results.CBeta
	ERI_AO = SCF_Results.ERI
	Ee_UHF=SCF_Results.Ee
	VNN_UHF=SCF_Results.VNN
	Etot_UHF=SCF_Results.Etot


	ENum = ENumAlpha + ENumBeta
	NNum = length(Molecule)

	SONum=2*ONum
	SONum=2*ONum

	HcoreAlphaMO=CAlpha'*Hcore*CAlpha
	HcoreBetaMO=CBeta'*Hcore*CBeta

	hMO=[HcoreAlphaMO zeros(ONum, ONum); zeros(ONum, ONum) HcoreBetaMO]
	println("\nTransforming ERIs to MO basis...")
	ERI_aaaa=TransERI(ERI_AO, CAlpha, CAlpha, CAlpha, CAlpha, ONum)
	ERI_bbbb=TransERI(ERI_AO, CBeta, CBeta, CBeta, CBeta, ONum)
	ERI_aabb=TransERI(ERI_AO, CAlpha, CAlpha, CBeta, CBeta, ONum)
	println("ERI transformation complete.")

	OccAlpha=1:ENumAlpha
	OccBeta=(ONum+1):(ONum+ENumBeta)
	RefDeterminants=sort(vcat(collect(OccAlpha), collect(OccBeta)))

	if MaxExcitation==-1
		NVir=SONum-ENum
		FCILvl=min(ENum, NVir)
		@printf("\nMaxExcit = -1 detected. Setting excitation level to %d for Full CI.\n", FCILvl)
		MaxExcitation=FCILvl
	end

	@printf("Generating CI space for excitation level <= %d...\n", MaxExcitation)
	Determinants=GenCISpace(RefDeterminants, 1:SONum, MaxExcitation)
	@printf("Number of determinants in the final CI space: %d\n", length(Determinants))

	println("\nBuilding and diagonalizing CI Matrix...")
	CI_Matrix = [CalcCI(Determinants[i], Determinants[j], hMO, ONum, ERI_aaaa, ERI_bbbb, ERI_aabb) for i in eachindex(Determinants), j in eachindex(Determinants)]
	E_CI, C_CI = eigen(CI_Matrix)
	println("CI Matrix diagonalization complete.")

	Ee_CI = E_CI[1]
	VNN_CI = VNN_UHF
	Etot_CI = Ee_CI + VNN_CI
	Ecorr = Etot_CI - Etot_UHF
	EStates = E_CI .+ VNN_CI

	@printf("\n--- Final Energy Results (CI) ---\n")
	@printf("Total Energy      = %.10f Hartree\n", Etot_CI)
	@printf("Electronic Energy = %.10f Hartree\n", Ee_CI)
	@printf("Nuclear Repulsion =  %.10f Hartree\n", VNN_CI)
	return CIResults(
		Etot_UHF,
		Etot_CI,
		Ecorr,
		EStates,
		C_CI,
		Determinants,
		MaxExcitation,
		length(Determinants),
	)
end



end




