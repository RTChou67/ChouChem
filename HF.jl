using LinearAlgebra
using Printf
using SpecialFunctions

include("Definitions.jl")
include("GetBasisList.jl")
include("CalcS.jl")
using .Definitions: PGTF, CGTF, Basis, Atom
using .GetBasisList: generate_basis_list, get_basis_set
using .CalcS: Sij




molecule = [
	Atom("H", 1, "STO-3G", (0.0, 0.0, +1.0)),
	Atom("O", 8, "STO-3G", (0.0, +1.0, 0.0)),
	Atom("H", 1, "STO-3G", (0.0, 0.0, -1.0)),
]

BasisSet = generate_basis_list(molecule)


Num=length(BasisSet)
ENum=sum(atom.Z for atom in molecule)
NNum=length(molecule)

S=[Sij(BasisSet[i], BasisSet[j]) for i in 1:Num, j in 1:Num]

function print_formatted_matrix(matrix::Matrix{Float64})
    n = size(matrix, 1)
    labels = [string(i) for i in 1:n] # 动态生成标签 "1", "2", ...
    
    # 打印列标题
    @printf("%-5s", "") # 左对齐的标签列
    for j in 1:n
        @printf("%10s", labels[j])
    end
    println()
    println(repeat("-", 5 + 10*n)) # 打印分隔线
    
    # 打印矩阵内容
    for i in 1:n
        @printf("%-5s", labels[i] * "|") # 打印行标题和分隔符
        for j in 1:n
            @printf("%10.6f", matrix[i, j]) # 格式化浮点数
        end
        println()
    end
end

print_formatted_matrix(S)

#=
function STVij(GTF1::Basis, GTF2::Basis, AllAtom::Vector{Atom})
	Sij=0.0
	Tij=0.0
	Vij=0.0
	R1=GTF1.position
	R2=GTF2.position
	R12=sqrt((R1[1]-R2[1])^2+(R1[2]-R2[2])^2+(R1[3]-R2[3])^2)
	for pgtf1 in GTF1.GTFs
		for pgtf2 in GTF2.GTFs
			alpha1, coeff1, norms1=pgtf1.alpha, pgtf1.coeff, pgtf1.norms
			alpha2, coeff2, norms2=pgtf2.alpha, pgtf2.coeff, pgtf2.norms
			p=alpha1+alpha2
			m=alpha1*alpha2/p
			K12=exp(-m*R12^2)
			S00=(π/p)^(3/2)*K12
			Sij+=coeff1*coeff2*norms1*norms2*S00
			Tij+=coeff1*coeff2*norms1*norms2*m*(3-2*m*R12^2)*S00
			Rxc=(alpha1*R1[1]+alpha2*R2[1])/p
			Ryc=(alpha1*R1[2]+alpha2*R2[2])/p
			Rzc=(alpha1*R1[3]+alpha2*R2[3])/p
			for atom in AllAtom
				Ra=atom.position
				Rac=sqrt((Ra[1]-Rxc)^2+(Ra[2]-Ryc)^2+(Ra[3]-Rzc)^2)
				Vij+=-atom.Z*coeff1*coeff2*norms1*norms2*K12*(2π/p)*BoysF0(p*Rac^2)
			end
		end
	end
	return Sij, Tij, Vij
end

function BoysF0(t::Float64)
	if t==0
		return 1.0
	elseif t<10e-8
		return 1.0 - t/3 + t^2/10 - t^3/42 + t^4/216
	else
		return 0.5*sqrt(π/t)*erf(sqrt(t))
	end
end

function CalcERI(Index1::Int, Index2::Int, Index3::Int, Index4::Int)
	eri = 0.0
	GTF1 = BasisSet[Index1]
	GTF2 = BasisSet[Index2]
	GTF3 = BasisSet[Index3]
	GTF4 = BasisSet[Index4]
	Ri, Rj, Rk, Rl = GTF1.position, GTF2.position, GTF3.position, GTF4.position
	Rijsq = sum((Ri .- Rj) .^ 2)
	Rklsq = sum((Rk .- Rl) .^ 2)
	for i in GTF1.GTFs
		for j in GTF2.GTFs
			alpha_i, coeff_i, norm_i = i.alpha, i.coeff, i.norms
			alpha_j, coeff_j, norm_j = j.alpha, j.coeff, j.norms
			p1 = alpha_i + alpha_j
			R1 = (alpha_i .* Ri .+ alpha_j .* Rj) ./ p1
			Kij = exp(-alpha_i * alpha_j / p1 * Rijsq)
			for k in GTF3.GTFs
				for l in GTF4.GTFs
					alpha_k, coeff_k, norm_k = k.alpha, k.coeff, k.norms
					alpha_l, coeff_l, norm_l = l.alpha, l.coeff, l.norms

					p2 = alpha_k + alpha_l
					R2 = (alpha_k .* Rk .+ alpha_l .* Rl) ./ p2
					Kkl = exp(-alpha_k * alpha_l / p2 * Rklsq)

					R12sq = sum((R1 .- R2) .^ 2)
					t = (p1 * p2) / (p1 + p2) * R12sq

					prefactor = (2 * π^2.5) / (p1 * p2 * sqrt(p1 + p2))
					prim_eri = prefactor * Kij * Kkl * BoysF0(t)

					norms_prod = norm_i * norm_j * norm_k * norm_l
					coeffs_prod = coeff_i * coeff_j * coeff_k * coeff_l

					eri += norms_prod * coeffs_prod * prim_eri
				end
			end
		end
	end

	return eri
end


println(BasisSet[1])


S=[STVij(BasisSet[i], BasisSet[j], molecule)[1] for i in 1:Num, j in 1:Num]
T=[STVij(BasisSet[i], BasisSet[j], molecule)[2] for i in 1:Num, j in 1:Num]
V=[STVij(BasisSet[i], BasisSet[j], molecule)[3] for i in 1:Num, j in 1:Num]
ERI=[CalcERI(i, j, k, l) for i in 1:Num, j in 1:Num, k in 1:Num, l in 1:Num]
Hcore=T+V
F=Hcore

MaxIter=100
Threshold=1e-6


X=S^(-0.5)
Fprime=X'*F*X
E, Cprime=eigen(Fprime)
C=X*Cprime
P=[2*sum(C[i, m]*C[j, m] for m in 1:1) for i in 1:Num, j in 1:Num]
for i in 1:MaxIter
	global P, F, Fprime, E, Cprime, C, G
	G=[sum(P[k, l]*(ERI[i, j, k, l]-0.5*ERI[i, k, j, l]) for k in 1:Num, l in 1:Num) for i in 1:Num, j in 1:Num]
	F=Hcore+G
	Fprime=X'*F*X
	E, Cprime=eigen(Fprime)
	C=X*Cprime
	Pnew=[2*sum(C[i, m]*C[j, m] for m in 1:1) for i in 1:Num, j in 1:Num]
	if norm(P-Pnew)<Threshold
		break
	end
	P=Pnew
end


VNN = sum(molecule[i].Z * molecule[j].Z / sqrt(sum((molecule[i].position .- molecule[j].position) .^ 2)) for i in 1:NNum for j in (i+1):NNum)


Ee = 0.5 * sum(P .* (Hcore + F))
Etot = Ee + VNN


Te = sum(P .* T)
VNe = sum(P .* V)
Vee = 0.5 * sum(P .* G)
Etot_components = Te + VNe + Vee + VNN


@printf("\n--- Molecular Structure ---\n")
for i in 1:NNum
	atom = molecule[i]
	@printf("Atom %d: %s (Z=%d) at (%.3f, %.3f, %.3f) Å\n", i, atom.symbol, atom.Z, atom.position[1]*0.52918, atom.position[2]*0.52918, atom.position[3]*0.52918)
end



@printf("\n--- Final Energy Results ---\n")
@printf("Total Energy      = %.10f Hartree\n", Etot)
@printf("Electronic Energy = %.10f Hartree\n", Ee)
@printf("Nuclear Repulsion =  %.10f Hartree\n", VNN)
=#
