struct CIResults
	EtotHF::Float64
	EtotCI::Float64
	Ecorr::Float64
	StateEnergies::Vector{Float64}
	CIVectors::Matrix{Float64}
	DeterminantSpace::Vector{Vector{Int}}
	MaxExcitation::Int
	NumberOfDeterminants::Int
end



function TransERI_blas(ERI_AO::Array{Float64, 4}, c1::Matrix{Float64}, c2::Matrix{Float64}, c3::Matrix{Float64}, c4::Matrix{Float64}, ONum::Int)
	N = ONum
	M_AO = reshape(ERI_AO, (N, N*N*N))

	M_tmp1 = zeros(N, N*N*N)
	mul!(M_tmp1, c1', M_AO)
	tmp_pjkl = reshape(M_tmp1, (N, N, N, N))
	tmp_pjkl_perm = permutedims(tmp_pjkl, (1, 3, 4, 2))
	M_tmp1_perm = reshape(tmp_pjkl_perm, (N*N*N, N))

	M_tmp2 = zeros(N*N*N, N)
	mul!(M_tmp2, M_tmp1_perm, c2)
	tmp_pqkl_perm = reshape(M_tmp2, (N, N, N, N))
	tmp_pqkl = permutedims(tmp_pqkl_perm, (1, 4, 2, 3))
	tmp_pqkl_perm = permutedims(tmp_pqkl, (1, 2, 4, 3))
	M_tmp2_perm = reshape(tmp_pqkl_perm, (N*N*N, N))

	M_tmp3 = zeros(N*N*N, N)
	mul!(M_tmp3, M_tmp2_perm, c3)
	tmp_pqrl_perm = reshape(M_tmp3, (N, N, N, N))
	tmp_pqrl = permutedims(tmp_pqrl_perm, (1, 2, 4, 3))
	M_tmp3_perm = reshape(tmp_pqrl, (N*N*N, N))

	ERI_MO_M = zeros(N*N*N, N)
	mul!(ERI_MO_M, M_tmp3_perm, c4)
	ERI_MO = reshape(ERI_MO_M, (N, N, N, N))
	return ERI_MO
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
	SortedRefDet = sort(RefDet)
	OccOrbs = SortedRefDet
	VirOrbs = collect(setdiff(Set(collect(OrbRange)), Set(OccOrbs)))
	println("\n--- Excitation Analysis ---")
	DetAllLvl = Set([SortedRefDet])
	@printf("Number of 0 -fold excitations: %d\n", 1)
	for Lvl in 1:MaxExcit
		if Lvl > length(OccOrbs) || Lvl > length(VirOrbs)
			@printf("Max excitation level (%d) reached capacity at Lvl=%d. Stopping.\n", MaxExcit, Lvl-1)
			break
		end
		DetThisLvl = Vector{Vector{Int}}()
		OccComb = combinations(OccOrbs, Lvl)
		VirComb = combinations(VirOrbs, Lvl)
		for Occ in OccComb
			BaseDetSet = setdiff(Set(OccOrbs), Set(Occ))
			for Vir in VirComb
				NewDetSet = union(BaseDetSet, Set(Vir))
				push!(DetThisLvl, sort(collect(NewDetSet)))
			end
		end
		CountThisLvl = length(DetThisLvl)
		@printf("Number of %-2d-fold excitations: %d\n", Lvl, CountThisLvl)
		union!(DetAllLvl, DetThisLvl)
	end
	println("--------------------------------------------")
	@printf("Total number of configurations: %d\n", length(DetAllLvl))
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
	S2 = sort(vcat(collect(setdiff(S1, [p, q])), [r, s]))
	PositionR = findfirst(==(r), S2)
	PositionS = findfirst(==(s), S2)
	Phase = (-1)^(PositionP + PositionQ + PositionR + PositionS)
	return Phase
end

function CalcCIHij(S1::Vector{Int64}, S2::Vector{Int64}, NDiff::Int, Diff1::Vector{Int}, Diff2::Vector{Int}, hMO::Matrix{Float64}, ONum::Int, ERI_aaaa, ERI_bbbb, ERI_aabb)
	if NDiff > 2
		return 0.0
	elseif NDiff == 0
		E1e = sum(hMO[i, i] for i in S1)
		E2e = sum(GetERI(i, i, j, j, ONum, ERI_aaaa, ERI_bbbb, ERI_aabb) - GetERI(i, j, i, j, ONum, ERI_aaaa, ERI_bbbb, ERI_aabb) for (idx, i) in enumerate(S1) for j in S1[(idx+1):end])
		return E1e + E2e
	elseif NDiff == 1
		p = Diff1[1]
		q = Diff2[1]
		PPos = findfirst(==(p), S1)
		QPos = findfirst(==(q), S2)
		Phase = (-1)^(PPos + QPos)

		S1Set = Set(S1)
		OrbCommon = sort(collect(setdiff(S1Set, Set(Diff1))))

		E1e = hMO[p, q]
		E2e = sum(GetERI(p, q, i, i, ONum, ERI_aaaa, ERI_bbbb, ERI_aabb) - GetERI(p, i, q, i, ONum, ERI_aaaa, ERI_bbbb, ERI_aabb) for i in OrbCommon)
		return (E1e + E2e)*Phase

	elseif NDiff == 2
		p, q = Diff1[1], Diff1[2]
		r, s = Diff2[1], Diff2[2]

		Phase = GetPhase(S1, p, q, r, s) # [cite: 8]
		return (GetERI(p, r, q, s, ONum, ERI_aaaa, ERI_bbbb, ERI_aabb) - GetERI(p, s, q, r, ONum, ERI_aaaa, ERI_bbbb, ERI_aabb)) * Phase
	end
end

function CalcUCI(SCF_Results::UHFResults, MaxExcitation::Int)
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
	HcoreAlphaMO=CAlpha'*Hcore*CAlpha
	HcoreBetaMO=CBeta'*Hcore*CBeta
	hMO=[HcoreAlphaMO zeros(ONum, ONum); zeros(ONum, ONum) HcoreBetaMO]
	println("\nTransforming ERIs to MO basis...")
	@time ERI_aaaa=TransERI_blas(ERI_AO, CAlpha, CAlpha, CAlpha, CAlpha, ONum)
	@time ERI_bbbb=TransERI_blas(ERI_AO, CBeta, CBeta, CBeta, CBeta, ONum)
	@time ERI_aabb=TransERI_blas(ERI_AO, CAlpha, CAlpha, CBeta, CBeta, ONum)
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
	Ndet = length(Determinants)
	I = Int[]
	J = Int[]
	V = Float64[]
	S_Sets = [Set(det) for det in Determinants]
	for i in 1:Ndet
		S1Set = S_Sets[i]
		S1 = Determinants[i]
		for j in i:Ndet
			S2Set = S_Sets[j]
			S2 = Determinants[j]
			Diff1 = sort(collect(setdiff(S1Set, S2Set)))
			NDiff = length(Diff1)
			if NDiff<=2
				Diff2 = sort(collect(setdiff(S2Set, S1Set)))
				if NDiff != length(Diff2)
					continue
				end
				hij = CalcCIHij(S1, S2, NDiff, Diff1, Diff2, hMO, ONum, ERI_aaaa, ERI_bbbb, ERI_aabb)
				if abs(hij) > 1e-12
					push!(I, i)
					push!(J, j)
					push!(V, hij)
					if i != j
						push!(I, j)
						push!(J, i)
						push!(V, hij)
					end
				end
			end
		end
	end


	println("Constructing sparse matrix from $(length(V)) triplets...")
	CI_Matrix = sparse(I, J, V, Ndet, Ndet)
	E_CI, C_CI = eigs(CI_Matrix, nev = 30, which = :SR)
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


function RunUCI(MolInAng::Vector{Atom}, Charge::Int, Multiplicity::Int, MaxExcitation::Int)
	TStart=time_ns()
	Bohr2Ang = 0.52917721092
	Molecule = [Atom(atom.symbol, atom.Z, atom.basis_set, atom.position ./ Bohr2Ang) for atom in MolInAng]
	@printf("\n--- Molecular Structure ---\n")
	for atom in MolInAng
		@printf("Atom: %-2s at (%8.4f, %8.4f, %8.4f) Å\n", atom.symbol, atom.position...)
	end
	println("---------------------------\n")
	SCF_Results=UHF_SCF(Molecule, Charge, Multiplicity, MaxIter = 128, Threshold = 1e-10)
	if isnothing(SCF_Results)
		error("UHF calculation did not converge. Aborting.")
		return
	else
		CI_Results=CalcUCI(SCF_Results, MaxExcitation)
		if isnothing(CI_Results)
			error("CI calculation did not converge. Aborting.")
			return
		end
	end
	println("\n--- Post-CI Analysis ---")
	println("\nGround State Wavefunction Analysis (Top 5 components):")
	GroundStateVector = CI_Results.CIVectors[:, 1]
	sorted_indices = sortperm(abs.(GroundStateVector), rev = true)
	RefDeterminant = CI_Results.DeterminantSpace[1]
	@printf("  Coeff    | Contribution  | Determinant Configuration\n")
	println("  -----------------------------------------------------")
	for i in 1:min(5, CI_Results.NumberOfDeterminants)
		idx = sorted_indices[i]
		coeff = GroundStateVector[idx]
		determinant = CI_Results.DeterminantSpace[idx]
		is_ref = (determinant == RefDeterminant) ? "(Ref)" : ""
		@printf("  %-8.4f | %-12.2f%% | %s %s\n", coeff, coeff^2 * 100, determinant, is_ref)
	end
	println("  -----------------------------------------------------\n")
	TEnd=time_ns()
	TSeconds = (TEnd - TStart) / 1e9
	days = floor(Int, TSeconds / 86400)
	hours = floor(Int, (TSeconds % 86400) / 3600)
	minutes = floor(Int, (TSeconds % 3600) / 60)
	seconds = TSeconds % 60
	DateTime = Dates.format(now(), "e u dd HH:MM:SS yyyy")
	@printf(" Job cpu time:       %d days %2d hours %2d minutes %5.1f seconds.\n", days, hours, minutes, seconds)
	@printf(" Elapsed time:       %d days %2d hours %2d minutes %5.1f seconds.\n", days, hours, minutes, seconds)
	println(" Normal termination of Julia UHF-CI at $(DateTime).")
	return (Etot = CI_Results.EtotCI,)
end

function RunRCI(MolInAng::Vector{Atom}, Charge::Int, Multiplicity::Int, MaxExcitation::Int)
	TStart=time_ns()
	Bohr2Ang = 0.52917721092
	Molecule = [Atom(atom.symbol, atom.Z, atom.basis_set, atom.position ./ Bohr2Ang) for atom in MolInAng]
	@printf("\n--- Molecular Structure ---\n")
	for atom in MolInAng
		@printf("Atom: %-2s at (%8.4f, %8.4f, %8.4f) Å\n", atom.symbol, atom.position...)
	end
	println("---------------------------\n")
	SCF_Results=RHF2UHF(RHF_SCF(Molecule, Charge, Multiplicity; MaxIter = 128, Threshold = 1e-10))
	if isnothing(SCF_Results)
		error("RHF calculation did not converge. Aborting.")
		return
	else
		CI_Results=CalcUCI(SCF_Results, MaxExcitation)
		if isnothing(CI_Results)
			error("CI calculation did not converge. Aborting.")
			return
		end
	end
	println("\n--- Post-CI Analysis ---")
	println("\nGround State Wavefunction Analysis (Top 5 components):")
	GroundStateVector = CI_Results.CIVectors[:, 1]
	sorted_indices = sortperm(abs.(GroundStateVector), rev = true)
	RefDeterminant = CI_Results.DeterminantSpace[1]
	@printf("  Coeff    | Contribution  | Determinant Configuration\n")
	println("  -----------------------------------------------------")
	for i in 1:min(5, CI_Results.NumberOfDeterminants)
		idx = sorted_indices[i]
		coeff = GroundStateVector[idx]
		determinant = CI_Results.DeterminantSpace[idx]
		is_ref = (determinant == RefDeterminant) ? "(Ref)" : ""
		@printf("  %-8.4f | %-12.2f%% | %s %s\n", coeff, coeff^2 * 100, determinant, is_ref)
	end
	println("  -----------------------------------------------------\n")
	TEnd=time_ns()
	TSeconds = (TEnd - TStart) / 1e9
	days = floor(Int, TSeconds / 86400)
	hours = floor(Int, (TSeconds % 86400) / 3600)
	minutes = floor(Int, (TSeconds % 3600) / 60)
	seconds = TSeconds % 60
	DateTime = Dates.format(now(), "e u dd HH:MM:SS yyyy")
	@printf(" Job cpu time:       %d days %2d hours %2d minutes %5.1f seconds.\n", days, hours, minutes, seconds)
	@printf(" Elapsed time:       %d days %2d hours %2d minutes %5.1f seconds.\n", days, hours, minutes, seconds)
	println(" Normal termination of Julia RHF-CI at $(DateTime).")
	return (Etot = CI_Results.EtotCI,)
end