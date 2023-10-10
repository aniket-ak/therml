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

working_dir = ARGS[4]
scenario_name = ARGS[6]
run_name = ARGS[8]

installation_dir = ENV["THERML"]

if run_name == ""
    run_name = now()
end

# julia -t 4 -working_dir "./" -scenario_name "name" -run_name "run1"
run_wd = joinpath(working_dir, run_name)
progress_file = open(joinpath(run_wd, "progress.txt"), "w")
global_logger(TerminalLogger(progress_file))

include("./functions_fvm_3d.jl")

f = open(joinpath(working_dir,"settings.json"), "r")
settings = JSON.parse(f)

temperature_base_units = settings["units"]["temperature"]

println(Dates.format(now(), "HH:MM:SS"))

with_logger(logger) do
    @info "\n-----\n"
    @info "Start time: "*Dates.format(now(), "HH:MM:SS")*"\n"
end

function solve_(working_dir, power_file)

    k = settings["model"]["bodies"]["die"]["material"]["k"]
    rho = settings["model"]["bodies"]["die"]["material"]["rho"]
    cp = settings["model"]["bodies"]["die"]["material"]["cp"]
    
    (delta_x, delta_y, delta_z), (Nx,Ny,Nz), (x_mesh, y_mesh, z_mesh) = create_mesh(settings)

    u0 = build_domain(Nx,Ny,Nz)
    
    source = define_volume_sources(working_dir, power_file,settings,Nx,Ny,Nz)

    p = (source, k, rho, cp, (delta_x, delta_y, delta_z))
    u0 = initialize_domain!(u0)

    # problem, algorithm = configure_problem_klu!(u0,p)
    # sol = solve(problem, algorithm, saveat=1.0,progress=true, maxiters=100, abstol=1e-3, reltol=1e-3)

    problem,_ = configure_problem_ode!(u0,p)
    sol = solve(problem, saveat=1.0,progress=true, progress_steps = 1,
                maxiters=1000, abstol=1e-4, reltol=1e-4)

    return sol
end

sol = solve_(working_dir, scenario_name);

println(Dates.format(now(), "HH:MM:SS"))
with_logger(logger) do
    @info "\n-----\n"
    @info Dates.format(now(), "HH:MM:SS")*"\n"
end

do_plotting(sol, false);
save_fields(sol);