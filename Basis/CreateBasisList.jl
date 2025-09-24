const SCRIPT_NAME = basename(@__FILE__)
const OUTPUT_FILE = "BasisList.jl"

function main()
	println("--- 开始生成 $OUTPUT_FILE ---")
	basis_dict = Dict{String, String}()
	name_regex = r"^\s*BasisSetName\s*=\s*\"(.*?)\"\s*$"
	try
		for filename in readdir(".")
			if filename == SCRIPT_NAME || filename == OUTPUT_FILE || !endswith(filename, ".jl") || !isfile(filename)
				continue
			end
			println("正在处理文件: $filename")
			try
				first_line = readline(filename)
				mat = match(name_regex, first_line)
				if mat !== nothing
					basis_name = mat.captures[1]
					basis_dict[basis_name] = filename
					println("  -> 找到基组名: '$basis_name'")
				else
					@warn "在 '$filename' 的第一行中未能找到有效的 'BasisSetName'。已跳过此文件。"
				end
			catch e
				if isa(e, EOFError)
					@warn "'$filename' 是一个空文件。已跳过。"
				else
					@error "读取文件 '$filename' 时出错: $e"
				end
			end
		end
		if isempty(basis_dict)
			@warn "未找到任何有效的基组文件。将不会创建 '$OUTPUT_FILE'。"
			return
		end
		try
			open(OUTPUT_FILE, "w") do f
				write(f, "# 此文件由 $SCRIPT_NAME 自动生成。\n")
				write(f, "# 它包含了基组名称到其源文件名的映射。\n\n")
				write(f, "BasisList = Dict(\n")
				sorted_keys = sort(collect(keys(basis_dict)))
				for key in sorted_keys
					value = basis_dict[key]
					write(f, "    \"$key\" => \"$value\",\n")
				end
				write(f, ")\n")
			end
			println("\n--- 成功创建 '$OUTPUT_FILE'，共包含 $(length(basis_dict)) 个条目。 ---")
		catch e
			@error "写入文件 '$OUTPUT_FILE' 失败: $e"
		end
	catch e
		@error "读取目录时发生错误: $e"
	end
end

main()
