using Revise
using Pkg
Pkg.activate(".")

println("Loading PowerModelsDistribution...")
using PowerModelsDistribution

# Include your runner script (Revise will track changes to it)
includet("run_simple_case.jl")

println("\n--- Environment Ready ---")
println("Run `run_case3_balanced()` to execute the test case.")
println("Edit `run_simple_case.jl` and run `run_case3_balanced()` again to see changes instantly.")
