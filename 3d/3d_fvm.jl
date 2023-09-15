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

global_logger(TerminalLogger())

include("./functions_fvm_3d.jl")

fw = open("./run.log","w")

f = open("/Users/aniket/Documents/MarlinSim/03_code/therml/3d/settings.json", "r")
settings = JSON.parse(f)

temperature_base_units = settings["units"]["temperature"]

println(Dates.format(now(), "HH:MM:SS"))

write(fw, "\n-----\n")
write(fw, "Start time:", Dates.format(now(), "HH:MM:SS"),"\n")

function solve_()

    k = settings["model"]["bodies"]["die"]["material"]["k"]
    rho = settings["model"]["bodies"]["die"]["material"]["rho"]
    cp = settings["model"]["bodies"]["die"]["material"]["cp"]
    
    (delta_x, delta_y, delta_z), (Nx,Ny,Nz), (x_mesh, y_mesh, z_mesh) = create_mesh(settings)

    u0 = build_domain(Nx,Ny,Nz)
    
    source = define_volume_sources(settings,Nx,Ny,Nz)

    p = (source, k, rho, cp, (delta_x, delta_y, delta_z))
    u0 = initialize_domain!(u0)

    # problem, algorithm = configure_problem_klu!(u0,p)
    # sol = solve(problem, algorithm, saveat=1.0,progress=true, maxiters=100, abstol=1e-3, reltol=1e-3)

    problem,_ = configure_problem_ode!(u0,p)
    sol = solve(problem, saveat=1.0,progress=true, progress_steps = 1,
                maxiters=1000, abstol=1e-4, reltol=1e-4)

    return sol
end

sol = solve_();

println(Dates.format(now(), "HH:MM:SS"))

do_plotting(sol, false);
save_fields(sol);