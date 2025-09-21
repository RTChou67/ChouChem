using Printf



module GetBasisList
using ..Definitions: PGTF, CGTF, Basis, Atom, CGTF
export generate_basis_list, get_basis_set
# --- 2. 基组库的按需加载系统 ---

# 用于缓存已加载的基组
const LoadedBasisSets = Dict{String, Dict{Int, Vector{CGTF}}}()
const BASIS_SET_DIR = "Basis"

# 基组管理器函数
function get_basis_set(name::String)
	if haskey(LoadedBasisSets, name)
		return LoadedBasisSets[name]
	end

	safe_filename = replace(name, r"[\s\(\)/]" => "-") * ".jl"
	filepath = joinpath(BASIS_SET_DIR, safe_filename)

	if !isfile(filepath)
		@error "基组文件未找到: $filepath"
		return nothing
	end

	try
		# 这里的 include 会返回一个 Dict{Int, Vector{CGTF}}
		basis_data = include(filepath)
		LoadedBasisSets[name] = basis_data
		return basis_data
	catch e
		@error "加载或解析基组文件 '$filepath' 失败: $e"
		return nothing
	end
end

# --- 3. 新的核心函数：生成基函数列表 ---

"""
	generate_basis_list(molecule::Vector{Atom})

接收一个分子（原子列表），返回一个包含该分子所有基函数的平面列表（Vector{Basis}）。
"""
function generate_basis_list(molecule::Vector{Atom})
	# 初始化一个空的 Basis 列表
	BasisList = Basis[]

	println("--- 开始为分子生成基函数列表 ---")

	# 遍历分子中的每一个原子
	for atom in molecule
		println("处理原子 '$(atom.symbol)'，基组为 '$(atom.basis_set)'...")

		# 1. 获取整个基组的数据（例如 "STO-3G" 的所有元素数据）
		basis_set_data = get_basis_set(atom.basis_set)
		if isnothing(basis_set_data)
			@warn "跳过原子 $(atom.symbol)，因为无法加载基组 $(atom.basis_set)。"
			continue
		end

		# 2. 从基组中提取该原子对应的数据
		element_data = basis_set_data[atom.Z]

		# 3. 遍历该原子的所有CGTF，并结合原子坐标创建新的 Basis 对象
		for CGTF in element_data
			# 核心步骤：将解析器中的数据与原子位置结合
			new_basis_function = Basis(
				CGTF.Type,
				CGTF.GTFs,
				atom.position,  # 使用当前原子的坐标
			)
			# 将新创建的基函数添加到总列表中
			push!(BasisList, new_basis_function)
		end
		println("  -> 为 '$(atom.symbol)' 添加了 $(length(element_data)) 个基函数。")
	end

	println("--- 基函数列表生成完毕 ---\n")
	# 返回最终的列表
	return BasisList
end

end