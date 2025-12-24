const DetType = Int128 

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

function pack_det(indices::Vector{Int})::DetType
    d = zero(Int128)
    for i in indices
        d |= (one(Int128) << (i - 1))
    end
    return d
end

function unpack_det(d::DetType)::Vector{Int}
    res = Int[]
    idx = 1
    while d > 0
        if (d & 1) == 1
            push!(res, idx)
        end
        d >>= 1
        idx += 1
    end
    return res
end

# --- [修复核心] 绝对稳健的相位计算函数 ---

# 计算轨道 k 在行列式 det 中的正则位置 (即 k 前面有多少个电子)
@inline function get_position(det::DetType, k::Int)
    # 掩码：覆盖 k 之前的所有位
    mask = (one(DetType) << (k - 1)) - 1
    # 计算 k 之前的 1 的个数
    return count_ones(det & mask)
end

# 单激发相位 (-1)^(pos_p + pos_q)
function get_phase_single(det_i::DetType, det_j::DetType, p::Int, q::Int)
    pos_p = get_position(det_i, p)
    pos_q = get_position(det_j, q)
    return iseven(pos_p + pos_q) ? 1.0 : -1.0
end

# 双激发相位 (-1)^(pos_p + pos_q + pos_r + pos_s)
function get_phase_double(det_i::DetType, det_j::DetType, p::Int, q::Int, r::Int, s::Int)
    # p, q 是 det_i 中的空穴
    pos_p = get_position(det_i, p)
    pos_q = get_position(det_i, q)
    
    # r, s 是 det_j 中的粒子
    pos_r = get_position(det_j, r)
    pos_s = get_position(det_j, s)
    
    return iseven(pos_p + pos_q + pos_r + pos_s) ? 1.0 : -1.0
end

# ==============================================================================
# 2. 积分变换与预处理
# ==============================================================================

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
    return reshape(ERI_MO_M, (N, N, N, N))
end

function build_spin_eri(aaaa, bbbb, aabb, ONum::Int)
    N = ONum
    N2 = 2 * N
    ERI = zeros(Float64, N2, N2, N2, N2)
    ERI[1:N, 1:N, 1:N, 1:N] .= aaaa
    ERI[N+1:2N, N+1:2N, N+1:2N, N+1:2N] .= bbbb
    ERI[1:N, 1:N, N+1:2N, N+1:2N] .= aabb
    ERI[N+1:2N, N+1:2N, 1:N, 1:N] .= permutedims(aabb, (3,4,1,2))
    return ERI
end

# ==============================================================================
# 3. CI Space Generation
# ==============================================================================
function GenCISpace_Optimized(RefDet::DetType, NSpinOrb::Int, MaxExcit::Int)
    println("\n--- Generating CI Space (Optimized) ---")
    
    # 1. 解析参考态
    Occ = unpack_det(RefDet) # 占据轨道列表
    Vir = setdiff(1:NSpinOrb, Occ) # 空轨道列表
    nocc = length(Occ)
    nvir = length(Vir)
    
    # 2. 预计算掩码表 (避免在循环中重复做位移运算)
    # Masks[i] = 1 << (i-1)
    # 注意：Orbital index 从 1 开始
    Masks = [one(DetType) << (i-1) for i in 1:NSpinOrb]
    
    # 3. 精确计算总行列式数量 (Binomial Coefficients)
    # N_det = 1 (Ref) + N_S + N_D ...
    n_singles = nocc * nvir
    n_doubles = (nocc * (nocc - 1) ÷ 2) * (nvir * (nvir - 1) ÷ 2)
    
    total_dets = 1 # Reference
    if MaxExcit >= 1; total_dets += n_singles; end
    if MaxExcit >= 2; total_dets += n_doubles; end
    
    # 4. 一次性分配内存 (Zero Allocation during loop)
    Dets = Vector{DetType}(undef, total_dets)
    Dets[1] = RefDet
    
    # 指针，指向当前填入位置
    idx = 2 
    
    # 5. 生成单激发 (Singles): i -> a
    if MaxExcit >= 1
        @inbounds for i_idx in 1:nocc
            i = Occ[i_idx]
            # 移除 i: Ref & ~mask[i]
            # 这里利用异或 (XOR) 更快，因为 Ref 中肯定有 i
            det_minus_i = RefDet ⊻ Masks[i] 
            
            for a_idx in 1:nvir
                a = Vir[a_idx]
                # 加入 a: det | mask[a] -> 同样可用 XOR 或 OR
                Dets[idx] = det_minus_i | Masks[a]
                idx += 1
            end
        end
    end
    
    # 6. 生成双激发 (Doubles): ij -> ab
    if MaxExcit >= 2
        @inbounds for i_idx in 1:nocc
            i = Occ[i_idx]
            det_minus_i = RefDet ⊻ Masks[i]
            
            for j_idx in (i_idx+1):nocc
                j = Occ[j_idx]
                det_minus_ij = det_minus_i ⊻ Masks[j]
                
                for a_idx in 1:nvir
                    a = Vir[a_idx]
                    det_plus_a = det_minus_ij | Masks[a]
                    
                    for b_idx in (a_idx+1):nvir
                        b = Vir[b_idx]
                        Dets[idx] = det_plus_a | Masks[b]
                        idx += 1
                    end
                end
            end
        end
    end
    
    # 7. 排序
    # 虽然生成是有序的，但按整数值排序有助于后续二分查找或块处理
    # Julia 的 sort! 对 Int128 非常快
    sort!(Dets)
    
    @printf("CI Space Size: %d (Ref: 1, S: %d, D: %d)\n", length(Dets), n_singles, n_doubles)
    return Dets
end
# ==============================================================================
# 4. Matrix Element Calculation
# ==============================================================================

function CalcUCI(SCF_Results::UHFResults, MaxExcitation::Int)
    ONum = length(SCF_Results.BasisSet)
    ENum = SCF_Results.ENumAlpha + SCF_Results.ENumBeta
    SONum = 2 * ONum
    
    C_a, C_b = SCF_Results.CAlpha, SCF_Results.CBeta
    H_a = C_a' * SCF_Results.Hcore * C_a
    H_b = C_b' * SCF_Results.Hcore * C_b
    
    h_so = zeros(SONum, SONum)
    h_so[1:ONum, 1:ONum] .= H_a
    h_so[ONum+1:end, ONum+1:end] .= H_b
    
    println("\nTransforming ERIs to MO basis...")
    @time begin
        ERI_aa = TransERI_blas(SCF_Results.ERI, C_a, C_a, C_a, C_a, ONum)
        ERI_bb = TransERI_blas(SCF_Results.ERI, C_b, C_b, C_b, C_b, ONum)
        ERI_ab = TransERI_blas(SCF_Results.ERI, C_a, C_a, C_b, C_b, ONum)
        V_so = build_spin_eri(ERI_aa, ERI_bb, ERI_ab, ONum)
    end
    println("ERI transformation complete.")
    
    RefOcc = vcat(collect(1:SCF_Results.ENumAlpha), 
                  collect(ONum+1 : ONum+SCF_Results.ENumBeta))
    RefDet = pack_det(sort(RefOcc))
    
    if MaxExcitation == -1
        NVir = SONum - ENum
        FCILvl = min(ENum, NVir)
        MaxExcitation = FCILvl
    end
    
    @printf("Generating CI space for excitation level <= %d...\n", MaxExcitation)
    Dets = GenCISpace_Optimized(RefDet, SONum, MaxExcitation)
    Ndet = length(Dets)
    
    println("\nBuilding and diagonalizing CI Matrix...")
    
    I_idx = Int[]
    J_idx = Int[]
    V_val = Float64[]
    sizehint!(I_idx, Ndet * 50)
    sizehint!(J_idx, Ndet * 50)
    sizehint!(V_val, Ndet * 50)
    
    for i in 1:Ndet
        Di = Dets[i]
        
        for j in i:Ndet
            Dj = Dets[j]
            diff = Di ⊻ Dj
            n_diff_bits = count_ones(diff)
            
            if n_diff_bits > 4; continue; end
            
            val = 0.0
            
            # === Case 0: Diagonal (i == j) ===
            if n_diff_bits == 0
                occ_i = unpack_det(Di)
                for orb in occ_i
                    val += h_so[orb, orb]
                end
                for p_idx in 1:length(occ_i)
                    p = occ_i[p_idx]
                    for q_idx in (p_idx+1):length(occ_i)
                        q = occ_i[q_idx]
                        val += V_so[p, p, q, q] - V_so[p, q, q, p]
                    end
                end
                
            # === Case 1: Single Excitation (2 bits diff) ===
            elseif n_diff_bits == 2
                # p in Di, q in Dj
                p = trailing_zeros(Di & diff) + 1
                q = trailing_zeros(Dj & diff) + 1
                
                phase = get_phase_single(Di, Dj, p, q)
                
                val = h_so[p, q]
                
                # Sum over common orbitals k
                common = Di & Dj
                temp_comm = common
                while temp_comm > 0
                    k = trailing_zeros(temp_comm) + 1
                    val += V_so[p, q, k, k] - V_so[p, k, k, q]
                    temp_comm &= ~(one(DetType) << (k-1))
                end
                val *= phase
                
            # === Case 2: Double Excitation (4 bits diff) ===
            elseif n_diff_bits == 4
                # p, q in Di (Holes)
                holes = Di & diff
                p = trailing_zeros(holes) + 1
                q = trailing_zeros(holes & ~(one(DetType) << (p-1))) + 1
                
                # r, s in Dj (Particles)
                parts = Dj & diff
                r = trailing_zeros(parts) + 1
                s = trailing_zeros(parts & ~(one(DetType) << (r-1))) + 1
                
                phase = get_phase_double(Di, Dj, p, q, r, s)
                
                # Formula: <pq||rs> = (pr|qs) - (ps|qr)
                val = (V_so[p, r, q, s] - V_so[p, s, q, r]) * phase
            end
            
            if abs(val) > 1e-12
                push!(I_idx, i)
                push!(J_idx, j)
                push!(V_val, val)
                if i != j
                    push!(I_idx, j)
                    push!(J_idx, i)
                    push!(V_val, val)
                end
            end
        end
    end
    
    println("Constructing sparse matrix from $(length(V_val)) triplets...")
    @time CI_Matrix = sparse(I_idx, J_idx, V_val, Ndet, Ndet)
    println("Sparse matrix construction complete.")
    println("Diagonalizing CI Matrix ...")
    n_roots = min(Ndet, 30)
    @time E_CI_raw, C_CI = eigs(CI_Matrix, nev = n_roots, which = :SR)
    println("CI Matrix diagonalization complete.")
    
    Ee_CI = E_CI_raw[1]
    VNN_CI = SCF_Results.VNN
    Etot_CI = Ee_CI + VNN_CI
    Ecorr = Etot_CI - SCF_Results.Etot
    EStates = E_CI_raw .+ VNN_CI
    
    @printf("\n--- Final Energy Results (CI) ---\n")
    @printf("Total Energy      = %.10f Hartree\n", Etot_CI)
    @printf("Electronic Energy = %.10f Hartree\n", Ee_CI)
    @printf("Nuclear Repulsion =  %.10f Hartree\n", VNN_CI)
    
    Dets_Vec = [unpack_det(d) for d in Dets]
    
    return CIResults(
        SCF_Results.Etot,
        Etot_CI,
        Ecorr,
        EStates,
        C_CI,
        Dets_Vec,
        MaxExcitation,
        Ndet
    )
end

function RunUCI(MolInAng::Vector{Atom}, Charge::Int, Multiplicity::Int, MaxExcitation::Int)
    TStart = time_ns()
    Bohr2Ang = 0.52917721092
    Molecule = [Atom(atom.symbol, atom.Z, atom.basis_set, atom.position ./ Bohr2Ang) for atom in MolInAng]
    
    @printf("\n--- Molecular Structure ---\n")
    for atom in MolInAng
        @printf("Atom: %-2s at (%8.4f, %8.4f, %8.4f) Å\n", atom.symbol, atom.position...)
    end
    println("---------------------------\n")
    
    SCF_Results = UHF_SCF(Molecule, Charge, Multiplicity, MaxIter = 128, Threshold = 1e-8)
    if isnothing(SCF_Results)
        error("UHF calculation did not converge. Aborting.")
        return
    end
    
    CI_Results = CalcUCI(SCF_Results, MaxExcitation)
    if isnothing(CI_Results)
        error("CI calculation did not converge. Aborting.")
        return
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
        det_str = string(determinant)
        is_ref = (determinant == RefDeterminant) ? "(Ref)" : ""
        @printf("  %-8.4f | %-12.2f%% | %s %s\n", coeff, coeff^2 * 100, det_str, is_ref)
    end
    println("  -----------------------------------------------------\n")
    
    TEnd = time_ns()
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
    TStart = time_ns()
    Bohr2Ang = 0.52917721092
    Molecule = [Atom(atom.symbol, atom.Z, atom.basis_set, atom.position ./ Bohr2Ang) for atom in MolInAng]
    
    @printf("\n--- Molecular Structure ---\n")
    for atom in MolInAng
        @printf("Atom: %-2s at (%8.4f, %8.4f, %8.4f) Å\n", atom.symbol, atom.position...)
    end
    println("---------------------------\n")
    
    RHF_Res = RHF_SCF(Molecule, Charge, Multiplicity; MaxIter = 128, Threshold = 1e-8)
    if isnothing(RHF_Res)
        error("RHF calculation did not converge. Aborting.")
        return
    end
    SCF_Results = RHF2UHF(RHF_Res)
    
    CI_Results = CalcUCI(SCF_Results, MaxExcitation)
    if isnothing(CI_Results)
        error("CI calculation did not converge. Aborting.")
        return
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
        det_str = string(determinant)
        is_ref = (determinant == RefDeterminant) ? "(Ref)" : ""
        @printf("  %-8.4f | %-12.2f%% | %s %s\n", coeff, coeff^2 * 100, det_str, is_ref)
    end
    println("  -----------------------------------------------------\n")
    
    TEnd = time_ns()
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