using PowerModelsDistribution
# Use the native fixed-point power flow solver (much faster for PF tasks)
function run_case3_balanced()
    dss_file = joinpath(dirname(pathof(PowerModelsDistribution)), "../test/data/opendss/case3_balanced.dss")
    println("Running Power Flow on: $dss_file")

    START_TIME = time()
    data = parse_file(dss_file)
    # Prefer the native (compute_mc_pf) solver for power flow — it's much faster
    # than building and solving a JuMP optimization model. Tune `max_iter` and
    # `stat_tol` to trade off speed vs accuracy.
    result = compute_mc_pf(data; max_iter=50, stat_tol=1e-6)

    println("Solver Status: ", get(result, "termination_status", "unknown"))
    println("Iterations: ", get(result, "iterations", "N/A"))
    println("Total PF Time: ", get(result, "time_total", "N/A"))
    println("\nCase run successfully!")
    println("Total Execution Time: ", time() - START_TIME, " seconds")
    
    return result
end

# Execute the function so clicking the editor "Run/Play" button executes the case.
# This file is intended as a small runnable script, so call it unconditionally.
try
    run_case3_balanced()
catch err
    # Print error and rethrow to make failures visible in editors/REPL
    println("Error running run_case3_balanced():", err)
    rethrow(err)
end
