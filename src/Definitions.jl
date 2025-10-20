struct PGTF
	alpha::Float64
	coeff::Float64
	norms::Float64
end

struct CGTF
	Type::NTuple{3, Int64}
	GTFs::Vector{PGTF}
end

struct Basis
	Type::NTuple{3, Int64}
	GTFs::Vector{PGTF}
	position::NTuple{3, Float64}
end

struct Atom
	symbol::String
	Z::Int
	basis_set::String
	position::NTuple{3, Float64}
end
