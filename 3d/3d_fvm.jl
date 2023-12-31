module therml

using DifferentialEquations, LinearAlgebra, SparseArrays
using PlotlyJS
using Symbolics
using JSON
using Dates
using LoopVectorization
using CSV, DataFrames
using Interpolations
using Logging: global_logger
using TerminalLoggers: TerminalLogger
using JLD2
using HDF5
using Logging
using LoggingExtras

function julia_main()::Cint
    try
        real_main()
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end

function real_main()
    working_dir = ARGS[4]
    scenario_name = ARGS[6]
    run_name = ARGS[8]
    installation_dir = ENV["THERML"]

    if run_name == ""
        run_name = now()
    end

    run_wd = joinpath(working_dir, run_name)
    log_wd = joinpath(run_wd,"Logs")
    sol_wd = joinpath(run_wd,"Solution")
    temp_wd = joinpath(run_wd,"Temp")

    progress_file = open(joinpath(temp_wd, scenario_name*"__progress.txt"), "w")
    global_logger(TerminalLogger(progress_file))

    include("./functions_fvm_3d.jl")

    f = open(joinpath(working_dir,"settings.json"), "r")
    settings = JSON.parse(f)
    temperature_base_units = settings["units"]["temperature"]
    println(Dates.format(now(), "HH:MM:SS"))
    sol = solve_(working_dir, scenario_name);
    do_plotting(sol, false);
    save_fields(sol);
    println(Dates.format(now(), "HH:MM:SS"))
end

end