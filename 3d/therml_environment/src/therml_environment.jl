module therml_environment

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
    file = open("./log.txt","w")
    println(file, ARGS)
    println(file,"working_dir:  ",ARGS[4])
    working_dir = ARGS[4]
    println(file,"scenario_name:  ",ARGS[6])
    scenario_name = ARGS[6]
    println(file,"run_name:  ", ARGS[8])
    run_name = ARGS[8]

    if run_name == ""
        run_name = now()
    end

    run_wd = joinpath(working_dir, run_name)
    log_wd = joinpath(run_wd,"Logs")
    sol_wd = joinpath(run_wd,"Solution")
    temp_wd = joinpath(run_wd,"Temp")
    println(file, "temp file:  ",temp_wd)

    progress_file_name = joinpath(temp_wd, scenario_name*"__progress.txt")
    println(file,progress_file_name)
    progress_file = open(progress_file_name, "w")
    global_logger(TerminalLogger(progress_file))

    include("/Users/aniket/Documents/MarlinSim/03_code/therml/3d/therml_environment/src/functions_fvm_3d.jl")

    f = open(joinpath(working_dir,"settings.json"), "r")
    settings = JSON.parse(f)
    println(file,Dates.format(now(), "HH:MM:SS"))
    # sol = solve_(working_dir, scenario_name);
    sol = Base.invokelatest(solve_, working_dir, scenario_name, settings);
    # do_plotting(sol, sol_wd, false);
    Base.invokelatest(do_plotting, sol, sol_wd, settings, scenario_name, false)
    # save_fields(sol,sol_wd);
    Base.invokelatest(save_fields, sol, sol_wd, scenario_name)
    println(file,Dates.format(now(), "HH:MM:SS"))
end

end