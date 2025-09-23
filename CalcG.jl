module CalcG

using ..Definitions: PGTF, Basis
using SpecialFunctions: binomial, gamma

export CalcERI

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

function RtuvG(t::Int, u::Int, v::Int, n::Int, pgtf1::PGTF, pgtf2::PGTF, pgtf3::PGTF, pgtf4::PGTF, R1::NTuple{3, Float64}, R2::NTuple{3, Float64}, R3::NTuple{3, Float64}, R4::NTuple{3, Float64})
	p1=pgtf1.alpha+pgtf2.alpha
	p2=pgtf3.alpha+pgtf4.alpha
	p=p1*p2/(p1+p2)
	P1X=(pgtf1.alpha*R1[1]+pgtf2.alpha*R2[1])/p1
	P1Y=(pgtf1.alpha*R1[2]+pgtf2.alpha*R2[2])/p1
	P1Z=(pgtf1.alpha*R1[3]+pgtf2.alpha*R2[3])/p1
	P2X=(pgtf3.alpha*R3[1]+pgtf4.alpha*R4[1])/p2
	P2Y=(pgtf3.alpha*R3[2]+pgtf4.alpha*R4[2])/p2
	P2Z=(pgtf3.alpha*R3[3]+pgtf4.alpha*R4[3])/p2
	XP1P2=P1X-P2X
	YP1P2=P1Y-P2Y
	ZP1P2=P1Z-P2Z
	RP1P22=XP1P2^2+YP1P2^2+ZP1P2^2
	if t<0 || u<0 || v<0 || n<0
		return 0.0
	elseif t==0 && u==0 && v==0 && n>=0
		return (-2*p)^n*boys(n, p*RP1P22)
	elseif t==0 && u==0 && v>=1
		return (v-1)*RtuvG(0, 0, v-2, n+1, pgtf1, pgtf2, pgtf3, pgtf4, R1, R2, R3, R4)+ZP1P2*RtuvG(0, 0, v-1, n+1, pgtf1, pgtf2, pgtf3, pgtf4, R1, R2, R3, R4)
	elseif t==0 && u>=1
		return (u-1)*RtuvG(0, u-2, v, n+1, pgtf1, pgtf2, pgtf3, pgtf4, R1, R2, R3, R4)+YP1P2*RtuvG(0, u-1, v, n+1, pgtf1, pgtf2, pgtf3, pgtf4, R1, R2, R3, R4)
	elseif t>=1
		return (t-1)*RtuvG(t-2, u, v, n+1, pgtf1, pgtf2, pgtf3, pgtf4, R1, R2, R3, R4)+XP1P2*RtuvG(t-1, u, v, n+1, pgtf1, pgtf2, pgtf3, pgtf4, R1, R2, R3, R4)
	end
end



function Gijkl(basis1::Basis, basis2::Basis, basis3::Basis, basis4::Basis)
	eri=0
	l1=basis1.Type
	l2=basis2.Type
	l3=basis3.Type
	l4=basis4.Type
	R1=basis1.position
	R2=basis2.position
	R3=basis3.position
	R4=basis4.position
	for pgtf1 in basis1.GTFs
		alpha1 = pgtf1.alpha
		coeff1 = pgtf1.coeff
		norm1 = pgtf1.norms

		for pgtf2 in basis2.GTFs
			alpha2 = pgtf2.alpha
			coeff2 = pgtf2.coeff
			norm2 = pgtf2.norms

			for pgtf3 in basis3.GTFs
				alpha3 = pgtf3.alpha
				coeff3 = pgtf3.coeff
				norm3 = pgtf3.norms

				for pgtf4 in basis4.GTFs
					alpha4 = pgtf4.alpha
					coeff4 = pgtf4.coeff
					norm4 = pgtf4.norms

					p1=alpha1 + alpha2
					q2=alpha3 + alpha4
					I1=0.0
					for t in 0:(l1[1]+l2[1])
						for u in 0:(l1[2]+l2[2])
							for v in 0:(l1[3]+l2[3])
								E1x=Hij_1D(t, l1[1], l2[1], R1[1], R2[1], alpha1, alpha2)
								E1y=Hij_1D(u, l1[2], l2[2], R1[2], R2[2], alpha1, alpha2)
								E1z=Hij_1D(v, l1[3], l2[3], R1[3], R2[3], alpha1, alpha2)
								E1=E1x*E1y*E1z
								for tau in 0:(l3[1]+l4[1])
									for niu in 0:(l3[2]+l4[2])
										for phi in 0:(l3[3]+l4[3])
											E2x=Hij_1D(tau, l3[1], l4[1], R3[1], R4[1], alpha3, alpha4)
											E2y=Hij_1D(niu, l3[2], l4[2], R3[2], R4[2], alpha3, alpha4)
											E2z=Hij_1D(phi, l3[3], l4[3], R3[3], R4[3], alpha3, alpha4)
											E2=E2x*E2y*E2z
											if (tau+niu+phi) % 2 == 0
												I1+=E1*E2*RtuvG(t+tau, u+niu, v+phi, 0, pgtf1, pgtf2, pgtf3, pgtf4, R1, R2, R3, R4)
											else
												I1-=E1*E2*RtuvG(t+tau, u+niu, v+phi, 0, pgtf1, pgtf2, pgtf3, pgtf4, R1, R2, R3, R4)
											end
										end
									end
								end
							end
						end
					end
					prefactor=(2*π^2.5)/(p1*q2*sqrt(p1+q2))
					eri+=norm1*norm2*norm3*norm4*coeff1*coeff2*coeff3*coeff4*prefactor*(I1)
				end
			end
		end
	end
	return eri
end




end
