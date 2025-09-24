# GetBasisList.jl
using Printf

module GetBasisList
using ..Definitions: PGTF, CGTF, Basis, Atom
export generate_basis_list, get_basis_set

const LoadedBasisSets = Dict{String, Dict{Int, Vector{CGTF}}}()
const BASIS_SET_DIR = "Basis"
const BasisNameToFile = include(joinpath(BASIS_SET_DIR, "BasisList.jl"))


function get_basis_set(name::String)
	if haskey(LoadedBasisSets, name)
		return LoadedBasisSets[name]
	end
	if !haskey(BasisNameToFile, name)
		error("在 BasisList.jl 中未找到基组名 '$name' 的路径。")
	end
	filename = BasisNameToFile[name]
	filepath = joinpath(BASIS_SET_DIR, filename)


	if !isfile(filepath)
		error("基组文件未找到: '$filepath'")
	end


	try
		basis_data = include(filepath)
		LoadedBasisSets[name] = basis_data
		return basis_data
	catch e

		error("加载或解析基组文件 '$filepath' 失败: $e")
	end
end



function generate_basis_list(molecule::Vector{Atom})

	BasisList = Basis[]

	println("--- 开始为分子生成基函数列表 ---")

	for atom in molecule
		println("处理原子 '$(atom.symbol)'，基组为 '$(atom.basis_set)'...")

		basis_set_data = get_basis_set(atom.basis_set)
		if isnothing(basis_set_data)

			error("无法为原子 $(atom.symbol) 加载基组 $(atom.basis_set)，计算终止。")
		end


		if !haskey(basis_set_data, atom.Z)

			error("在基组 '$(atom.basis_set)' 中未找到原子序数为 $(atom.Z) ('$(atom.symbol)') 的数据。计算终止。")
		end

		element_data = basis_set_data[atom.Z]

		for cgtf in element_data
			new_basis_function = Basis(
				cgtf.Type,
				cgtf.GTFs,
				atom.position,
			)
			push!(BasisList, new_basis_function)
		end
		println("  -> 为 '$(atom.symbol)' 添加了 $(length(element_data)) 个基函数。")
	end

	println("--- 基函数列表生成完毕 ---\n")
	return BasisList
end

end