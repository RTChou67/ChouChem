using Printf

# --- 1. 定义您的分子数据 ---
# (格式: "符号", (X, Y, Z))
molecules = Dict(
	"HF" => (
		coords =  HFCoord=[
			("H", (0.000, 0.000, 0.000)),
			("F", (0.000, 0.000, 1.500)),
		] ,
		charge = 0,
		mult = 1,
	),
	"H2O" => (
		coords = [
			("O", (0.0000, 0.0000, 0.0000)),
			("H", (0.0000, 1.0000, 0.8000)),
			("H", (0.0000, -1.0000, 0.8000)),
		],
		charge = 0,
		mult = 1,
	),
)

# --- 2. 定义您的计算设置 ---
basis_sets = ["STO-3G", "6-31G"]

# 将您的运行函数名映射到Gaussian的"Route Section"关键字
# 这是最关键的一步
calc_types = Dict(
	"RHF" => "RHF",
	"UHF" => "UHF",
	"RMP2" => "MP2",      # Gaussian中MP2默认使用RHF参考态
	"RCI" => "RCISD", # 强制使用RHF参考态进行CISD 
	"UCI" => "UCISD",  # 强制使用UHF参考态进行CISD 
)

"""
	generate_gaussian_input(mol_name, mol_data, basis, calc_key, calc_route)

一个辅助函数，用于生成单个Gaussian输入文件。
"""
function generate_gaussian_input(mol_name, mol_data, basis, calc_key, calc_route)
	charge = mol_data.charge
	mult = mol_data.mult

	# 格式化文件名，例如: "H2O_6-31G_RCI.gjf"
	filename = "$(mol_name)_$(replace(basis, "-"=>"_"))_$(calc_key).gjf"

	try
		open(filename, "w") do f
			# 1. Checkpoint 和资源
			write(f, "%Chk=$(mol_name)_$(replace(basis, "-"=>"_"))_$(calc_key).chk\n")
			write(f, "%Mem=4GB\n")
			write(f, "%NProcShared=8\n")

			# 2. Route Section (方法/基组)
			write(f, "# $(calc_route)/$(basis)\n\n")

			# 3. 标题
			write(f, "Title: $(mol_name) $(basis) $(calc_route) calculation\n\n")

			# 4. 电荷和多重度
			write(f, "$(charge) $(mult)\n")

			# 5. 坐标
			for (atom_symbol, atom_coords) in mol_data.coords
				line = @sprintf(" %-2s %12.8f %12.8f %12.8f\n",
					atom_symbol, atom_coords[1], atom_coords[2], atom_coords[3])
				write(f, line)
			end

			# 6. 末尾空行
			write(f, "\n")
		end
		println("Generated: $(filename)")
	catch e
		println("Error generating $(filename): $e")
	end
end

# --- 3. 主循环：遍历所有组合 ---
function generate_all_files()
	println("Starting Gaussian input file generation...")
	total_files = 0

	# 遍历每种分子
	for (mol_name, mol_data) in molecules
		# 遍历每种基组
		for basis in basis_sets
			# 遍历每种计算方法
			for (calc_key, calc_route) in calc_types
				generate_gaussian_input(mol_name, mol_data, basis, calc_key, calc_route)
				total_files += 1
			end
		end
	end

	println("\nDone. Generated $total_files input files.")
end

# --- 运行脚本 ---
generate_all_files()
