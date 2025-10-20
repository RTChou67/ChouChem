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

function CalcMatrices(BasisSet, Molecule)
	BNum = length(BasisSet)
	TimeS1=time_ns()
	S = [Sij(BasisSet[i], BasisSet[j]) for i in 1:BNum, j in 1:BNum]
	TimeS2=time_ns()
	println("Calculation for S Matrix took $( (TimeS2-TimeS1)/1e6 ) ms")
	TimeT1=time_ns()
	T = [Tij(BasisSet[i], BasisSet[j]) for i in 1:BNum, j in 1:BNum]
	TimeT2=time_ns()
	println("Calculation for T Matrix took $( (TimeT2-TimeT1)/1e6 ) ms")
	TimeV1=time_ns()
	V = [Vij(BasisSet[i], BasisSet[j], Molecule) for i in 1:BNum, j in 1:BNum]
	TimeV2=time_ns()
	println("Calculation for V Matrix took $( (TimeV2-TimeV1)/1e6 ) ms")
	TimeERI1=time_ns()
	ERI = [Gijkl(BasisSet[i], BasisSet[j], BasisSet[k], BasisSet[l]) for i in 1:BNum, j in 1:BNum, k in 1:BNum, l in 1:BNum]
	TimeERI2=time_ns()
	println("Calculation for ERI Tensor took $( (TimeERI2-TimeERI1)/1e6 ) ms")
	return S, T, V, ERI
end