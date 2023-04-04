using Clang.Generators

for pkg in ["GNC"]

    library_name = "lib$pkg"
    module_name = "Lib$pkg"

    options = Dict{String,Any}("general" => Dict{String,Any}())
    options["general"]["library_name"] = library_name
    options["general"]["module_name"] = module_name
    options["general"]["use_julia_native_enum_type"] = true


    # Compile into C shared library
    #
    # Optionally generate code from Matlab

    src_dir = normpath(joinpath(@__DIR__, "..", "src", "C"))
    C_src = joinpath(src_dir, "$pkg.c")

    run(`gcc -c -fPIC $C_src -I$src_dir -lm -o $pkg.o`)
    run(`gcc -shared -o $library_name.so $pkg.o`)

    mv("$library_name.so", joinpath(@__DIR__, "..", "lib", "$library_name.so"); force=true)

    lib_path = normpath(joinpath(@__DIR__, "..", "lib", "$library_name.so"))


    # Generate prologue.jl file
    #
    # Inserts arbitrary code at beginning of new Julia library file

    prologue_file_path = normpath(joinpath(@__DIR__, "prologue.jl"))
    if isfile(prologue_file_path)
        rm(prologue_file_path)
    end

    open(prologue_file_path, "w") do io
        write(io, "const lib$pkg = joinpath(@__DIR__, \"..\", \"lib\", \"$library_name.so\")")
    end

    options["general"]["prologue_file_path"] = prologue_file_path

    options["general"]["output_file_path"] = normpath(joinpath(@__DIR__, "..", "src", "$module_name.jl"))


    args = get_default_args()
    push!(args, "-I$src_dir")

    headers = detect_headers(src_dir, args)

    ctx = create_context(headers, args, options)

    build!(ctx)

end