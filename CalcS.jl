module CalcS
using ..Definitions: PGTF, CGTF, Basis, Atom
using Printf

export Sij

function Sij(basis1::Basis, basis2::Basis)
	Stot = 0.0
	for PGTF1 in basis1.GTFs
		for PGTF2 in basis2.GTFs
			alpha1, coeff1, norm1 = PGTF1.alpha, PGTF1.coeff, PGTF1.norms
			alpha2, coeff2, norm2 = PGTF2.alpha, PGTF2.coeff, PGTF2.norms
			R1 = basis1.position
			R2 = basis2.position
			L1 = basis1.Type
			L2 = basis2.Type
			MemoX = Dict{Tuple{Int, Int}, Float64}()
			MemoY = Dict{Tuple{Int, Int}, Float64}()
			MemoZ = Dict{Tuple{Int, Int}, Float64}()
			Sx = calculate_s_1d(L1[1], L2[1], alpha1, alpha2, R1[1], R2[1], MemoX)
			Sy = calculate_s_1d(L1[2], L2[2], alpha1, alpha2, R1[2], R2[2], MemoY)
			Sz = calculate_s_1d(L1[3], L2[3], alpha1, alpha2, R1[3], R2[3], MemoZ)
			SPrim = Sx * Sy * Sz
			Stot += coeff1 * coeff2 * norm1 * norm2 * SPrim
		end
	end
	return Stot
end



function calculate_s_1d(a::Int, b::Int, alpha1::Float64, alpha2::Float64, R1::Float64, R2::Float64, memo::Dict)
    if haskey(memo, (a, b))
        return memo[(a, b)]
    end
    if a < 0 || b < 0
        return 0.0
    end
    p = alpha1 + alpha2
    u = (alpha1 * alpha2) / p

    R12 = R1 - R2
    Rp = (alpha1 * R1 + alpha2 * R2)

    if a == 0 && b == 0
        result = sqrt(π / p) * exp(-u * R12^2)
        memo[(a, b)] = result
        return result
    end
    Rpa = Rp - R1
    Rpb = Rp - R2
    result = 0.0
    if a > b
        term1 = Rpa * calculate_s_1d(a - 1, b, alpha1, alpha2, R1, R2, memo)
        term2 = (a - 1) / (2p) * calculate_s_1d(a - 2, b, alpha1, alpha2, R1, R2, memo)
        term3 = b / (2p) * calculate_s_1d(a - 1, b - 1, alpha1, alpha2, R1, R2, memo)
        result = term1 + term2 + term3
    else
        term1 = Rpb * calculate_s_1d(a, b - 1, alpha1, alpha2, R1, R2, memo)
        term2 = (b - 1) / (2p) * calculate_s_1d(a, b - 2, alpha1, alpha2, R1, R2, memo)
        term3 = a / (2p) * calculate_s_1d(a - 1, b, alpha1, alpha2, R1, R2, memo)
        result = term1 + term2 + term3
    end
    memo[(a, b)] = result
    return result
end

end