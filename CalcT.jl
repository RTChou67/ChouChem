module CalcT

using ..Definitions: PGTF, Basis
using SpecialFunctions: binomial, gamma

export Tij


function double_factorial(m::Int)
	if m <= 1
		return 1.0
	end
	if isodd(m)
		n = (m + 1) ÷ 2
		return 2.0^n * gamma(n + 0.5) / sqrt(pi)
	else
		k = m ÷ 2
		return 2.0^k * factorial(k)
	end
end


function overlap_1D(l1::Int, l2::Int, A::Float64, B::Float64, alpha1::Float64, alpha2::Float64, p::Float64)
	# l < 0 的情况对应于求导后角动量为负，其积分为0
	if l1 < 0 || l2 < 0
		return 0.0
	end

	# 根据高斯乘积定理计算新中心 P 和相对位置 PA, PB
	P = (alpha1 * A + alpha2 * B) / p
	PA = P - A
	PB = P - B

	sum_val = 0.0
	for i in 0:l1
		for j in 0:l2
			if (i + j) % 2 == 0
				term = binomial(l1, i) * binomial(l2, j) *
					   PA^(l1 - i) * PB^(l2 - j) *
					   double_factorial(i + j - 1) / (2.0*p)^((i + j) / 2.0)
				sum_val += term
			end
		end
	end
	return sum_val
end


function kinetic_1D(l1::Int, l2::Int, A::Float64, B::Float64, alpha1::Float64, alpha2::Float64, p::Float64)
	# <∂g(l1)/∂x | ∂g(l2)/∂x> = <l1*g(l1-1) - 2*α1*g(l1+1) | l2*g(l2-1) - 2*α2*g(l2+1)>
	# 展开后得到四项重叠积分

	term1 = l1 * l2 * overlap_1D(l1 - 1, l2 - 1, A, B, alpha1, alpha2, p)
	term2 = -2.0 * alpha2 * l1 * overlap_1D(l1 - 1, l2 + 1, A, B, alpha1, alpha2, p)
	term3 = -2.0 * alpha1 * l2 * overlap_1D(l1 + 1, l2 - 1, A, B, alpha1, alpha2, p)
	term4 = 4.0 * alpha1 * alpha2 * overlap_1D(l1 + 1, l2 + 1, A, B, alpha1, alpha2, p)

	return term1 + term2 + term3 + term4
end


"""
	Tij(basis1::Basis, basis2::Basis)

计算两个收缩高斯基函数 (CGTF) 之间的动能积分 Tij。
"""
function Tij(basis1::Basis, basis2::Basis)
	T_total = 0.0

	R1 = basis1.position
	R2 = basis2.position
	l_vec1 = basis1.Type # (lx, ly, lz) for basis1
	l_vec2 = basis2.Type # (lx, ly, lz) for basis2

	# 遍历第一个基函数的所有原初高斯函数 (PGTF)
	for pgtf1 in basis1.GTFs
		alpha1 = pgtf1.alpha
		coeff1 = pgtf1.coeff
		norm1 = pgtf1.norms

		# 遍历第二个基函数的所有原初高斯函数 (PGTF)
		for pgtf2 in basis2.GTFs
			alpha2 = pgtf2.alpha
			coeff2 = pgtf2.coeff
			norm2 = pgtf2.norms

			# --- 根据高斯乘积定理计算参数 ---
			p = alpha1 + alpha2
			R12_sq = sum((R1 .- R2) .^ 2)
			K_ab = exp(-alpha1 * alpha2 / p * R12_sq)
			prefactor = K_ab * sqrt((π/p)^3)

			# --- 计算各维度上的 1D 重叠积分与 1D 动能积分项 ---
			Sx = overlap_1D(l_vec1[1], l_vec2[1], R1[1], R2[1], alpha1, alpha2, p)
			Sy = overlap_1D(l_vec1[2], l_vec2[2], R1[2], R2[2], alpha1, alpha2, p)
			Sz = overlap_1D(l_vec1[3], l_vec2[3], R1[3], R2[3], alpha1, alpha2, p)

			Tx_term = kinetic_1D(l_vec1[1], l_vec2[1], R1[1], R2[1], alpha1, alpha2, p)
			Ty_term = kinetic_1D(l_vec1[2], l_vec2[2], R1[2], R2[2], alpha1, alpha2, p)
			Tz_term = kinetic_1D(l_vec1[3], l_vec2[3], R1[3], R2[3], alpha1, alpha2, p)

			Tx_primitive = 0.5 * Tx_term * Sy * Sz
			Ty_primitive = 0.5 * Ty_term * Sx * Sz
			Tz_primitive = 0.5 * Tz_term * Sx * Sy

			T_primitive = prefactor * (Tx_primitive + Ty_primitive + Tz_primitive)

			T_total += coeff1 * coeff2 * norm1 * norm2 * T_primitive
		end
	end

	return T_total
end

end