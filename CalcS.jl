module CalcS

using ..Definitions: PGTF, Basis
using SpecialFunctions: binomial

export Sij

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
	# 根据高斯乘积定理计算新中心 P 和相对位置 PA, PB
	P = (alpha1 * A + alpha2 * B) / p
	PA = P - A
	PB = P - B

	sum_val = 0.0
	for i in 0:l1
		for j in 0:l2
			# 仅当 i + j 为偶数时，积分结果不为零
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

"""
	Sij(basis1::Basis, basis2::Basis)

计算两个收缩高斯基函数 (CGTF) 之间的重叠积分 Sij。
这两个基函数可以具有任意的角动量。
"""
function Sij(basis1::Basis, basis2::Basis)
	S_total = 0.0

	R1 = basis1.position
	R2 = basis2.position
	l1 = basis1.Type # (l_x, l_y, l_z) for basis1
	l2 = basis2.Type # (l_x, l_y, l_z) for basis2

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

			# 预指数因子
			R12_sq = sum((R1 .- R2) .^ 2)
			K_ab = exp(-alpha1 * alpha2 / p * R12_sq)

			# --- 计算三个一维方向上的重叠积分 (使用修正后的函数) ---
			# X方向
			Sx = overlap_1D(l1[1], l2[1], R1[1], R2[1], alpha1, alpha2, p)
			# Y方向
			Sy = overlap_1D(l1[2], l2[2], R1[2], R2[2], alpha1, alpha2, p)
			# Z方向
			Sz = overlap_1D(l1[3], l2[3], R1[3], R2[3], alpha1, alpha2, p)

			# --- 计算原初高斯函数对的重叠积分 ---
			S_primitive = K_ab * sqrt((π/p)^3) * Sx * Sy * Sz

			# --- 累加到总的重叠积分中 ---
			# 乘以收缩系数和归一化系数
			S_total += coeff1 * coeff2 * norm1 * norm2 * S_primitive
		end
	end

	return S_total
end

end
