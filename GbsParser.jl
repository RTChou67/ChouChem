using Printf

include("Definitions.jl")
using ..Definitions: PGTF, CGTF


const ELEMENTS = Dict(
	"H" => 1, "He" => 2, "Li" => 3, "Be" => 4, "B" => 5, "C" => 6, "N" => 7, "O" => 8, "F" => 9, "Ne" => 10,
	"Na" => 11, "Mg" => 12, "Al" => 13, "Si" => 14, "P" => 15, "S" => 16, "Cl" => 17, "Ar" => 18, "K" => 19, "Ca" => 20,
	"Sc" => 21, "Ti" => 22, "V" => 23, "Cr" => 24, "Mn" => 25, "Fe" => 26, "Co" => 27, "Ni" => 28, "Cu" => 29, "Zn" => 30,
	"Ga" => 31, "Ge" => 32, "As" => 33, "Se" => 34, "Br" => 35, "Kr" => 36, "Rb" => 37, "Sr" => 38, "Y" => 39, "Zr" => 40,
	"Nb" => 41, "Mo" => 42, "Tc" => 43, "Ru" => 44, "Rh" => 45, "Pd" => 46, "Ag" => 47, "Cd" => 48, "In" => 49, "Sn" => 50,
	"Sb" => 51, "Te" => 52, "I" => 53, "Xe" => 54, "Cs" => 55, "Ba" => 56, "La" => 57, "Ce" => 58, "Pr" => 59, "Nd" => 60,
	"Pm" => 61, "Sm" => 62, "Eu" => 63, "Gd" => 64, "Tb" => 65, "Dy" => 66, "Ho" => 67, "Er" => 68, "Tm" => 69, "Yb" => 70,
	"Lu" => 71, "Hf" => 72, "Ta" => 73, "W" => 74, "Re" => 75, "Os" => 76, "Ir" => 77, "Pt" => 78, "Au" => 79, "Hg" => 80,
	"Tl" => 81, "Pb" => 82, "Bi" => 83, "Po" => 84, "At" => 85, "Rn" => 86,
)

const SHELL_TYPE = Dict(
	"S" => [(0, 0, 0)],
	"P" => [(1, 0, 0), (0, 1, 0), (0, 0, 1)], # 顺序通常是 Px, Py, Pz
	"D" => [(2, 0, 0), (0, 2, 0), (0, 0, 2), (1, 1, 0), (1, 0, 1), (0, 1, 1)], # Dxx, Dyy, Dzz, Dxy, Dxz, Dyz
	"F" => [(3, 0, 0), (0, 3, 0), (0, 0, 3), (1, 2, 0), (2, 1, 0), (1, 0, 2), (2, 0, 1), (0, 1, 2), (0, 2, 1), (1, 1, 1)], # Fxxx, Fyyy, Fzzz, Fxyy, Fxxy, ...
)


using SpecialFunctions
function doublefactorial(m::Int)
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


function normalize_primitive(alpha::Float64, l::NTuple{3, Int}; pure::Bool = false)
	lx, ly, lz = l
	L = lx + ly + lz
	df_product = doublefactorial(2*lx - 1) * doublefactorial(2*ly - 1) * doublefactorial(2*lz - 1)
	norm_const = (2 * alpha / pi)^(0.75) * sqrt((4 * alpha)^L / df_product)

	return norm_const
end

function GetBasisName(filepath::String)
	try
		lines = readlines(filepath)
		if length(lines) >= 6
			line6 = lines[6]
			parts = split(line6, ':', limit = 2)
			if length(parts) == 2
				return String(strip(parts[2]))
			end
		end
	catch e
		println("读取基组名称时出错: $e")
	end
	@error "无法从文件第六行提取基组名称。请确保文件格式正确。"
	return nothing
end

function FindEcpElems(filepath::String)
	EcpElems = Set{String}()
	try
		for line in eachline(filepath)
			if contains(line, "-ECP")
				parts = split(strip(line))
				if !isempty(parts)
					SymbolParts = split(parts[1], '-')
					if !isempty(SymbolParts)
						ElemSymbol = uppercasefirst(lowercase(SymbolParts[1]))
						if haskey(ELEMENTS, ElemSymbol)
							push!(EcpElems, ElemSymbol)
						end
					end
				end
			end
		end
	catch e
		println("扫描ECP元素时出错: $e")
	end
	return EcpElems
end

function GBS2CGTF(filepath::String)
	BasisData = Dict{Int, Vector{CGTF}}()
	lines = readlines(filepath)
	LineIdx = 1
	while LineIdx <= length(lines)
		line = strip(lines[LineIdx])
		if isempty(line) || startswith(line, "!")
			LineIdx += 1
			continue
		end
		parts = split(line)
		if length(parts) == 2 && haskey(ELEMENTS, parts[1])
			ElemSymbol = parts[1]
			z = ELEMENTS[ElemSymbol]
			if !haskey(BasisData, z)
				BasisData[z] = CGTF[]
				LineIdx += 1
				while LineIdx <= length(lines)
					shell_line = strip(lines[LineIdx])
					if startswith(shell_line, "****")
						break
					end
					if isempty(shell_line)
						LineIdx += 1
						continue
					end

					shell_parts = split(shell_line)
					shell_type_str = uppercase(shell_parts[1])
					num_primitives = parse(Int, shell_parts[2])

					primitives_data = []
					for _ in 1:num_primitives
						LineIdx += 1
						prim_line = replace(strip(lines[LineIdx]), "D" => "E")
						prim_parts = [parse(Float64, p) for p in split(prim_line)]
						push!(primitives_data, prim_parts)
					end

					if haskey(SHELL_TYPE, shell_type_str)
						am_representative = SHELL_TYPE[shell_type_str][1]
						pgtf_vec = [PGTF(prim[1], prim[2], normalize_primitive(prim[1], am_representative)) for prim in primitives_data]

						for am_tuple in SHELL_TYPE[shell_type_str]
							push!(BasisData[z], CGTF(am_tuple, pgtf_vec))
						end

					elseif shell_type_str in ["SP", "L"]
						s_am = (0, 0, 0)
						s_pgtf_vec = [PGTF(prim[1], prim[2], normalize_primitive(prim[1], s_am)) for prim in primitives_data]
						push!(BasisData[z], CGTF(s_am, s_pgtf_vec))
						p_am_representative = (1, 0, 0)
						p_pgtf_vec = [PGTF(prim[1], prim[3], normalize_primitive(prim[1], p_am_representative)) for prim in primitives_data]

						for am_tuple in SHELL_TYPE["P"]
							push!(BasisData[z], CGTF(am_tuple, p_pgtf_vec))
						end
					end
					LineIdx += 1
				end
			end
		end
		LineIdx += 1
	end
	return BasisData
end

function WriteBasisLib(data::Dict, BasisName::String, OutputFile::String)
	open(OutputFile, "w") do f
		println(f, "# Basis Set Data for: $BasisName")
		println(f, "# This file returns a Dict{Int, Vector{CGTF}} when included.")
		println(f, "# Generated by a corrected version of GbsParser.jl")
		println(f, "Dict(")

		sorted_zs = sort(collect(keys(data)))
		for (z_idx, z) in enumerate(sorted_zs)
			ElemSymbol = findfirst(isequal(z), ELEMENTS)
			println(f, "\t$z => [ # $ElemSymbol")

			CGTFs = data[z]
			for (cgtf_idx, cgtf) in enumerate(CGTFs)
				print(f, "\t\tCGTF($(cgtf.Type), [")
				if !isempty(cgtf.GTFs)
					println(f)
					for (pgtf_idx, pgtf) in enumerate(cgtf.GTFs)
						line = @sprintf("\t\t\tPGTF(%.10e, %.10e, %.10e)", pgtf.alpha, pgtf.coeff, pgtf.norms)
						print(f, line)
						println(f, pgtf_idx < length(cgtf.GTFs) ? "," : "")
					end
					print(f, "\t\t])")
				else
					print(f, "])")
				end
				println(f, cgtf_idx < length(CGTFs) ? "," : "")
			end
			println(f, z_idx < length(sorted_zs) ? "\t]," : "\t]")
		end
		println(f, ")")
	end
end


function main()
	if length(ARGS) != 1
		println("用法: julia GbsParser.jl <input_gbs_file>")
		return
	end
	InputFile = ARGS[1]

	BasisName = GetBasisName(InputFile)
	if isnothing(BasisName)
		return
	end

	SafeBasisName = replace(BasisName, r"[\s\(\)/]" => "-")
	OutputFileName = SafeBasisName * ".jl"

	EcpElems = FindEcpElems(InputFile)

	println("正在为基组 '$BasisName' 解析文件 '$InputFile'...")
	ParsedData = GBS2CGTF(InputFile)

	if !isempty(EcpElems)
		println("检测到赝势(ECP)定义，将从库中移除以下元素: ", join(sort(collect(EcpElems)), ", "))
		for symbol in EcpElems
			if haskey(ELEMENTS, symbol)
				z = ELEMENTS[symbol]
				delete!(ParsedData, z)
			end
		end
	end

	println("正在将全电子基组库写入到 '$OutputFileName'...")
	WriteBasisLib(ParsedData, BasisName, OutputFileName)

	println("完成. 文件 '$OutputFileName' 已成功创建。")
end

main()
