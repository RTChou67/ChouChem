@inline function Hij_1D(Idx::Int, l1::Int, l2::Int, R1D1::Float64, R1D2::Float64, alpha1::Float64, alpha2::Float64, exp_pre::Float64)
    p = alpha1 + alpha2
    inv_2p = 0.5 / p
    P1D = (alpha1 * R1D1 + alpha2 * R1D2) / p
    PA1D = P1D - R1D1
    PB1D = P1D - R1D2
    if Idx == 0 && l1 == 0 && l2 == 0
        return exp_pre 
    elseif Idx < 0 || Idx > (l1 + l2)
        return 0.0
    end
    val = 0.0
    if l1 == 0
        val = PB1D * Hij_1D(Idx, 0, l2-1, R1D1, R1D2, alpha1, alpha2, exp_pre) + 
              (l2 > 1 ? (l2-1) * inv_2p * Hij_1D(Idx, 0, l2-2, R1D1, R1D2, alpha1, alpha2, exp_pre) : 0.0) + 
              inv_2p * Hij_1D(Idx-1, 0, l2-1, R1D1, R1D2, alpha1, alpha2, exp_pre)
    else
        val = PA1D * Hij_1D(Idx, l1-1, l2, R1D1, R1D2, alpha1, alpha2, exp_pre) + 
              (l1 > 1 ? (l1-1) * inv_2p * Hij_1D(Idx, l1-2, l2, R1D1, R1D2, alpha1, alpha2, exp_pre) : 0.0) + 
              l2 * inv_2p * Hij_1D(Idx, l1-1, l2-1, R1D1, R1D2, alpha1, alpha2, exp_pre) + 
              inv_2p * Hij_1D(Idx-1, l1-1, l2, R1D1, R1D2, alpha1, alpha2, exp_pre)
    end
    return val
end

@inline function RtuvG(t::Int, u::Int, v::Int, n::Int, p::Float64, X::Float64, Y::Float64, Z::Float64, F_vals::Vector{Float64})
    if t == 0 && u == 0 && v == 0
        return (-2.0 * p)^n * @inbounds F_vals[n+1]
    end
    
    val = 0.0
    if v > 0
        val = (v-1) * RtuvG(t, u, v-2, n+1, p, X, Y, Z, F_vals) + 
               Z * RtuvG(t, u, v-1, n+1, p, X, Y, Z, F_vals)
    elseif u > 0
        val = (u-1) * RtuvG(t, u-2, v, n+1, p, X, Y, Z, F_vals) + 
               Y * RtuvG(t, u-1, v, n+1, p, X, Y, Z, F_vals)
    elseif t > 0
        val = (t-1) * RtuvG(t-2, u, v, n+1, p, X, Y, Z, F_vals) + 
               X * RtuvG(t-1, u, v, n+1, p, X, Y, Z, F_vals)
    end
    return val
end

struct ERIWorkspace
    Ex1::Vector{Float64}
    Ey1::Vector{Float64}
    Ez1::Vector{Float64}
    Ex2::Vector{Float64}
    Ey2::Vector{Float64}
    Ez2::Vector{Float64}
    Cx::Vector{Float64}
    Cy::Vector{Float64}
    Cz::Vector{Float64}
    F_vals::Vector{Float64}
    
    function ERIWorkspace(max_L::Int=20)
        new(zeros(max_L), zeros(max_L), zeros(max_L),
            zeros(max_L), zeros(max_L), zeros(max_L),
            zeros(2*max_L), zeros(2*max_L), zeros(2*max_L),
            zeros(4*max_L))
    end
end

function Gijkl(basis1::Basis, basis2::Basis, basis3::Basis, basis4::Basis)
    eri = 0.0
    PGTFs1, PGTFs2 = basis1.GTFs, basis2.GTFs
    PGTFs3, PGTFs4 = basis3.GTFs, basis4.GTFs
    R1, R2, R3, R4 = basis1.position, basis2.position, basis3.position, basis4.position
    l1, l2, l3, l4 = basis1.Type, basis2.Type, basis3.Type, basis4.Type
    ws = ERIWorkspace(20)

    lim1_x = l1[1] + l2[1]
    lim1_y = l1[2] + l2[2]
    lim1_z = l1[3] + l2[3]
    
    lim2_x = l3[1] + l4[1]
    lim2_y = l3[2] + l4[2]
    lim2_z = l3[3] + l4[3]
    
    Max_N = lim1_x + lim2_x
    Max_M = lim1_y + lim2_y
    Max_K = lim1_z + lim2_z
    L_sum_max = Max_N + Max_M + Max_K

    @inbounds for pgtf1 in PGTFs1
        alpha1 = pgtf1.alpha
        c1 = pgtf1.coeff * pgtf1.norms
        
        for pgtf2 in PGTFs2
            alpha2 = pgtf2.alpha
            c2 = pgtf2.coeff * pgtf2.norms

            p1_exp = alpha1 + alpha2
            P1 = (alpha1 .* R1 .+ alpha2 .* R2) ./ p1_exp

            u1 = (alpha1 * alpha2) / p1_exp
            exp_pre1_x = exp(-u1 * (R1[1]-R2[1])^2)
            exp_pre1_y = exp(-u1 * (R1[2]-R2[2])^2)
            exp_pre1_z = exp(-u1 * (R1[3]-R2[3])^2)
            
            for t in 0:lim1_x; ws.Ex1[t+1] = Hij_1D(t, l1[1], l2[1], R1[1], R2[1], alpha1, alpha2, exp_pre1_x); end
            for u_idx in 0:lim1_y; ws.Ey1[u_idx+1] = Hij_1D(u_idx, l1[2], l2[2], R1[2], R2[2], alpha1, alpha2, exp_pre1_y); end
            for v in 0:lim1_z; ws.Ez1[v+1] = Hij_1D(v, l1[3], l2[3], R1[3], R2[3], alpha1, alpha2, exp_pre1_z); end
            
            coeff12 = c1 * c2

            for pgtf3 in PGTFs3
                alpha3 = pgtf3.alpha
                c3 = pgtf3.coeff * pgtf3.norms
                
                for pgtf4 in PGTFs4
                    alpha4 = pgtf4.alpha
                    c4 = pgtf4.coeff * pgtf4.norms

                    p2_exp = alpha3 + alpha4
                    P2 = (alpha3 .* R3 .+ alpha4 .* R4) ./ p2_exp
                    
                    u2 = (alpha3 * alpha4) / p2_exp
                    exp_pre2_x = exp(-u2 * (R3[1]-R4[1])^2)
                    exp_pre2_y = exp(-u2 * (R3[2]-R4[2])^2)
                    exp_pre2_z = exp(-u2 * (R3[3]-R4[3])^2)
                    
                    for tau in 0:lim2_x; ws.Ex2[tau+1] = Hij_1D(tau, l3[1], l4[1], R3[1], R4[1], alpha3, alpha4, exp_pre2_x); end
                    for niu in 0:lim2_y; ws.Ey2[niu+1] = Hij_1D(niu, l3[2], l4[2], R3[2], R4[2], alpha3, alpha4, exp_pre2_y); end
                    for phi in 0:lim2_z; ws.Ez2[phi+1] = Hij_1D(phi, l3[3], l4[3], R3[3], R4[3], alpha3, alpha4, exp_pre2_z); end

                    p = p1_exp * p2_exp / (p1_exp + p2_exp)
                    RP1P2 = P1 .- P2
                    X_P1P2, Y_P1P2, Z_P1P2 = RP1P2
                    T_val = p * sum(RP1P2 .^ 2)

                    for n in 0:L_sum_max
                        ws.F_vals[n+1] = boys(n, T_val)
                    end

                    fill!(view(ws.Cx, 1:Max_N+1), 0.0)
                    fill!(view(ws.Cy, 1:Max_M+1), 0.0)
                    fill!(view(ws.Cz, 1:Max_K+1), 0.0)
                    
                    for t in 0:lim1_x
                        val1 = ws.Ex1[t+1]
                        @simd for tau in 0:lim2_x
                            sign = 1 - 2*(tau&1)
                            ws.Cx[t+tau+1] += val1 * ws.Ex2[tau+1] * sign
                        end
                    end
                    
                    for u_idx in 0:lim1_y
                        val1 = ws.Ey1[u_idx+1]
                        @simd for niu in 0:lim2_y
                            sign = 1 - 2*(niu&1)
                            ws.Cy[u_idx+niu+1] += val1 * ws.Ey2[niu+1] * sign
                        end
                    end
                    
                    for v in 0:lim1_z
                        val1 = ws.Ez1[v+1]
                        @simd for phi in 0:lim2_z
                            sign = 1 - 2*(phi&1)
                            ws.Cz[v+phi+1] += val1 * ws.Ez2[phi+1] * sign
                        end
                    end
                    
                    I_val = 0.0
                    for N in 0:Max_N
                        val_x = ws.Cx[N+1]
                        if abs(val_x) > 1e-15 
                            for M in 0:Max_M
                                val_y = ws.Cy[M+1]
                                if abs(val_y) > 1e-15
                                    for K in 0:Max_K
                                        val_z = ws.Cz[K+1]
                                        if abs(val_z) > 1e-15
                                            R_val = RtuvG(N, M, K, 0, p, X_P1P2, Y_P1P2, Z_P1P2, ws.F_vals)
                                            I_val += val_x * val_y * val_z * R_val
                                        end
                                    end
                                end
                            end
                        end
                    end
                    
                    prefactor = (2 * π^2.5) / (p1_exp * p2_exp * sqrt(p1_exp + p2_exp))
                    eri += coeff12 * c3 * c4 * I_val * prefactor
                end
            end
        end
    end
    return eri
end