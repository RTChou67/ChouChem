module CalcV

using ..Definitions: PGTF, Basis, Atom
using SpecialFunctions: binomial, gamma, erf
using LinearAlgebra

export Vij

function double_factorial(m::Int)
	if m <= 1
		return 1.0
	end
	if isodd(m)
		n = (m + 1) / 2
		return 2.0^n * gamma(n + 0.5) / sqrt(pi)
	else
		k = m / 2
		return 2.0^k * factorial(k)
	end
end

function boys(m::Int, x::Float64)
    if x == 0.0
        return 1.0 / (2m + 1)
    end


    if x < 1e-8
        return 1.0 / (2m + 1) - x / (2m + 3)
    end

    s = m + 0.5

    return 0.5 * x^(-s) * (gamma(s) - gamma(s, x))
end

# ----------------- 1D Hermite-like 系数 -----------------
function Hij_1D(Idx, l1::Int64, l2::Int64, R1D1::Float64, R1D2::Float64, alpha1, alpha2)
	p=alpha1+alpha2
	u = (alpha1 * alpha2) / p
	P1D = (alpha1 * R1D1 + alpha2 * R1D2) / p
	PA1D = P1D - R1D1
	PB1D = P1D - R1D2
	PAB = R1D1 - R1D2
	E0_00=exp(-u * PAB^2)
	if Idx<0 || l1<0 || l2<0 || Idx>(l1+l2)
		return 0.0
	elseif Idx==0 && l1==0 && l2==0
		return E0_00
	elseif l1==0 && l2>=1
		return PB1D*Hij_1D(Idx, 0, l2-1, R1D1, R1D2, alpha1, alpha2)+((l2-1)*Hij_1D(Idx, 0, l2-2, R1D1, R1D2, alpha1, alpha2)+Hij_1D(Idx-1, 0, l2-1, R1D1, R1D2, alpha1, alpha2))/(2*p)
	elseif l1>=1
		return PA1D*Hij_1D(Idx, l1-1, l2, R1D1, R1D2, alpha1, alpha2)+((l1-1)*Hij_1D(Idx, l1-2, l2, R1D1, R1D2, alpha1, alpha2)+l2*Hij_1D(Idx-1, l1-1, l2, R1D1, R1D2, alpha1, alpha2)+Hij_1D(Idx-1, l1-1, l2, R1D1, R1D2, alpha1, alpha2))/(2*p)
	end
end

function Rtuv(t::Int, u::Int, v::Int, n::Int, pgtf1::PGTF, pgtf2::PGTF, R1::NTuple{3, Float64}, R2::NTuple{3, Float64}, atom::Atom)
	p=pgtf1.alpha+pgtf2.alpha
	PX=(pgtf1.alpha*R1[1]+pgtf2.alpha*R2[1])/p
	PY=(pgtf1.alpha*R1[2]+pgtf2.alpha*R2[2])/p
	PZ=(pgtf1.alpha*R1[3]+pgtf2.alpha*R2[3])/p
	XPC=PX-atom.position[1]
	YPC=PY-atom.position[2]
	ZPC=PZ-atom.position[3]
	RPC2=XPC^2+YPC^2+ZPC^2
	if t<0 || u<0 || v<0 || n<0
		return 0.0
	elseif t==0 && u==0 && v==0 && n>=0
		return (-2*p)^n*boys(n, RPC2*p)
	elseif t==0 && u==0 && v>=1
		return (v-1)*Rtuv(0, 0, v-2, n+1, pgtf1, pgtf2, R1, R2, atom)+ZPC*Rtuv(0, 0, v-1, n+1, pgtf1, pgtf2, R1, R2, atom)
	elseif t==0 && u>=1
		return (u-1)*Rtuv(0, u-2, v, n+1, pgtf1, pgtf2, R1, R2, atom)+YPC*Rtuv(0, u-1, v, n+1, pgtf1, pgtf2, R1, R2, atom)
	elseif t>=1
		return (t-1)*Rtuv(t-2, u, v, n+1, pgtf1, pgtf2, R1, R2, atom)+XPC*Rtuv(t-1, u, v, n+1, pgtf1, pgtf2, R1, R2, atom)
	end
end

function Vij(basis1::Basis, basis2::Basis, atoms::Vector{Atom})
	V_total = 0.0
	XIdx1, YIdx1, ZIdx1 = basis1.Type
	XIdx2, YIdx2, ZIdx2 = basis2.Type

	R1 = basis1.position
	R2 = basis2.position

	for pgtf1 in basis1.GTFs
		alpha1 = pgtf1.alpha
		coeff1 = pgtf1.coeff
		norm1 = pgtf1.norms

		for pgtf2 in basis2.GTFs
			alpha2 = pgtf2.alpha
			coeff2 = pgtf2.coeff
			norm2 = pgtf2.norms

			p = alpha1 + alpha2
			for t in 0:(XIdx1+XIdx2)
				for u in 0:(YIdx1+YIdx2)
					for v in 0:(ZIdx1+ZIdx2)
						Hx=Hij_1D(t, XIdx1, XIdx2, R1[1], R2[1], alpha1, alpha2)
						Hy=Hij_1D(u, YIdx1, YIdx2, R1[2], R2[2], alpha1, alpha2)
						Hz=Hij_1D(v, ZIdx1, ZIdx2, R1[3], R2[3], alpha1, alpha2)
						for atom in atoms
							R=Rtuv(t, u, v, 0, pgtf1, pgtf2, R1, R2, atom)
							V_total+=coeff1*coeff2*norm1*norm2*Hx*Hy*Hz*R*(-atom.Z*(2π/p))
						end
					end
				end
			end
		end

	end
	return V_total
end

end
