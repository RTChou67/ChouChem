module RMP2

using LinearAlgebra
using Printf

using ..RHF: RHFResults
using ..Definitions: Atom

export RunMP2


function TransG(G_AO::Array{Float64, 4}, C::Matrix{Float64})
	NMO = size(C, 2)

	println("\n--- Starting Integral Transformation ---")
	println("[1/4] Transforming First Index...")
	Tmp1 = zeros(NMO, NMO, NMO, NMO)
	for p in 1:NMO, j in 1:NMO, k in 1:NMO, l in 1:NMO
		sum_val = 0.0
		for i in 1:NMO
			sum_val += C[i, p] * G_AO[i, j, k, l]
		end
		Tmp1[p, j, k, l] = sum_val
	end

	println("[2/4] Transforming Second Index...")
	Tmp2 = zeros(NMO, NMO, NMO, NMO)
	for p in 1:NMO, q in 1:NMO, k in 1:NMO, l in 1:NMO
		sum_val = 0.0
		for j in 1:NMO
			sum_val += C[j, q] * Tmp1[p, j, k, l]
		end
		Tmp2[p, q, k, l] = sum_val
	end

	println("[3/4] Transforming Third Index...")
	Tmp3 = zeros(NMO, NMO, NMO, NMO)
	for p in 1:NMO, q in 1:NMO, r in 1:NMO, l in 1:NMO
		sum_val = 0.0
		for k in 1:NMO
			sum_val += C[k, r] * Tmp2[p, q, k, l]
		end
		Tmp3[p, q, r, l] = sum_val
	end

	println("[4/4] Transforming Fourth Index...")
	G_MO = zeros(NMO, NMO, NMO, NMO)
	for p in 1:NMO, q in 1:NMO, r in 1:NMO, s in 1:NMO
		sum_val = 0.0
		for l in 1:NMO
			sum_val += C[l, s] * Tmp3[p, q, r, l]
		end
		G_MO[p, q, r, s] = sum_val
	end

	println("--- Integral Transformation Finished ---\n")
	return G_MO
end


function CalcRMP2(G_MO::Array{Float64, 4}, E_MO::Vector{Float64}, NOcc::Int)
	NMO = length(E_MO)
	NVir = NMO - NOcc

	EOcc = E_MO[1:NOcc]
	EVir = E_MO[(NOcc+1):NMO]

	EMP2 = 0.0

	for i in 1:NOcc
		for j in 1:NOcc
			for a in 1:NVir
				for b in 1:NVir
					aIdx = a + NOcc
					bIdx = b + NOcc

					Denom = EOcc[i] + EOcc[j] - EVir[a] - EVir[b]

					Int_iajb = G_MO[i, aIdx, j, bIdx]
					Int_ibja = G_MO[i, bIdx, j, aIdx]

					EMP2 += Int_iajb * (2 * Int_iajb - Int_ibja) / Denom
				end
			end
		end
	end
	return EMP2
end


function RunRMP2(RHF_Results::RHFResults)
	println("--- Starting Restricted MP2 Calculation ---")

	C = RHF_Results.C
	e = RHF_Results.E
	G_AO = RHF_Results.ERI
	ENum = RHF_Results.ENum
	EtotRHF = RHF_Results.Etot

	NOcc = ENum ÷ 2
	G_MO = TransG(G_AO, C)

	E_RMP2_Corr = CalcRMP2(G_MO, e, NOcc)
	E_RMP2_total = EtotRHF + E_RMP2_Corr

	@printf("\n--- MP2 Energy Results ---\n")
	@printf("MP2 Correlation Energy = %.10f Hartree\n", E_RMP2_Corr)
	@printf("Total MP2 Energy       = %.10f Hartree\n", E_RMP2_total)
	println("----------------------------\n")

	return (MP2_corr = E_RMP2_Corr, MP2_total = E_RMP2_total)
end

end
