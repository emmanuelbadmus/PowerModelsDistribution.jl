using PackageCompiler

println("Starting sysimage build...")
println("This will take a few minutes (time to grab a coffee ☕).")

create_sysimage(
    [:PowerModelsDistribution, :SCS, :JuMP],
    sysimage_path="pmd_sysimage.so",
    precompile_execution_file="run_simple_case.jl"
)

println("\nSysimage built successfully! Saved as 'pmd_sysimage.so'")
println("You can now run fast using:")
println("  julia --project=. --sysimage pmd_sysimage.so run_simple_case.jl")
