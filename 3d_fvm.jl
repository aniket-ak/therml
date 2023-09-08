using DifferentialEquations, LinearAlgebra, SparseArrays
using PlotlyJS
using Symbolics
using JSON
using Dates
using LoopVectorization
include("./functions_fvm_3d.jl")

fw = open("./run.log","w")

f = open("settings.json", "r")
settings = JSON.parse(f)

temperature_base_units = settings["units"]["temperature"]

function create_mesh(settings)
    x_length = settings["model"]["bodies"]["die"]["size"]["X"]
    y_length = settings["model"]["bodies"]["die"]["size"]["Y"]
    z_length = settings["model"]["bodies"]["die"]["size"]["Z"]

    x_normalized = x_length/z_length
    y_normalized = y_length/z_length
    z_normalized = z_length/z_length

    Nx = 2*settings["model"]["num_sources"]["X"] + 2 # including ghost cells
    Ny = 2*settings["model"]["num_sources"]["Y"] + 2 # including ghost cells
    Nz = max(2, 2*int(settings["model"]["bodies"]["die"]["size"]["Z"]/settings["model"]["smallest_thickness"]))
    
    x_mesh = range(0, stop = x_normalized, length = Nx+1)
    y_mesh = range(0, stop = y_normalized, length = Ny+1)
    z_mesh = range(0, stop = z_normalized, length = Nz+1)

    delta_x, delta_y, delta_z = step(x_mesh), step(y_mesh), step(z_mesh)
    return (delta_x, delta_y, delta_z), (Nx,Ny,Nz), (x_mesh, y_mesh, z_mesh)
end

function build_domain(Nx,Ny,Nz)
    u = zeros(Nx, Ny, Nz)
    println("total mesh count: ", (Nx-2)*(Ny-2)*(Nz-2)/1e3, " k")
    return u
end

function initialize_domain!(u0)
    T_initial = convert_units("temperature", settings["IC"], temperature_base_units, "K")
    u0 = u0 .+ T_initial
    return u0
end

function define_volume_sources(Nx,Ny,Nz)
    source = zeros(Nx,Ny,Nz)
    source[10:20, 10:20, :] .= 1
    return source
end

function do_plotting(settings, sol)
    result_ = convert_units("temperature", sol[end][2:end-1,2:end-1,2], "K", "C")'
    # result_ = sol[end][2:end-1,2:end-1,2]
    (delta_x, delta_y, delta_z), (Nx,Ny,Nz), (x_mesh, y_mesh, z_mesh) = create_mesh(settings)
    z_length = settings["model"]["bodies"]["die"]["size"]["Z"]
    plot(contour(x=x_mesh*z_length,y=y_mesh*z_length,z=result_, colorscale="Jet",colorbar=attr(
        title="Temperature", titleside="right",titlefont=attr(size=14,family="Arial, sans-serif")
    )),Layout(autosize=true))
end

println(Dates.format(now(), "HH:MM:SS"))

write(fw, "\n-----\n")
write(fw, "Start time:", Dates.format(now(), "HH:MM:SS"),"\n")

function conduction_3d_loop!(du, u, p, t)
    source, k, rho, cp, (delta_x, delta_y, delta_z) = p
    alpha_x = (k/(rho*cp)) / delta_x^2
    alpha_y = (k/(rho*cp)) / delta_y^2
    alpha_z = (k/(rho*cp)) / delta_z^2

    Nx, Ny, Nz = size(u)

    Threads.@threads for k_ in range(2, Nz-1)
        Threads.@threads for j_ in range(2, Ny-1)
            Threads.@threads for i_ in range(2, Nx-1)
                ip1, im1, jp1, jm1, kp1, km1 = i_ + 1, i_ - 1, j_ + 1, j_ - 1, k_+1, k_-1
                du[i_, j_, k_] =  alpha_x * (u[im1, j_, k_] + u[ip1, j_, k_] - 2u[i_,j_, k_])+
                                    alpha_y * (u[i_, jp1, k_] + u[i_, jm1, k_] - 2u[i_, j_, k_]) + 
                                    alpha_z* (u[i_, j_, kp1] + u[i_, j_, km1] - 2u[i_, j_, k_]) +  
                                    source[i_, j_, k_]
            end
        end
    end

    apply_bc(u, settings, delta_x, delta_y, delta_z)

end

function solve_!()

    k = settings["model"]["bodies"]["die"]["material"]["k"]
    rho = settings["model"]["bodies"]["die"]["material"]["rho"]
    cp = settings["model"]["bodies"]["die"]["material"]["cp"]
    
    (delta_x, delta_y, delta_z), (Nx,Ny,Nz), (x_mesh, y_mesh, z_mesh) = create_mesh(settings)

    u0 = build_domain(Nx,Ny,Nz)
    
    source = define_volume_sources(Nx,Ny,Nz)

    p = (source, k, rho, cp, (delta_x, delta_y, delta_z))
    u0 = initialize_domain!(u0)

    # du0 = zeros(size(u0))
    # jac_sparsity = Symbolics.jacobian_sparsity((du, u) -> conduction_3d_loop!(du, u, p, 0.0), du0, u0)

    # f = ODEFunction(conduction_3d_loop!; jac_prototype = float.(jac_sparsity))
    # prob_conduction_3d_sparse = ODEProblem(f, u0, (0.0, 10.0), p)

    # sol = solve(prob_conduction_3d_sparse, 
    #             KenCarp47(linsolve = KLUFactorization()),
    #             save_everystep = false,
    #             progress=true,
    #             maxiters=100,
    #             abstol=1e-3, 
    #             reltol=1e-3);

    problem = ODEProblem(conduction_3d_loop!, u0, (0.0, 10.0), p)
    sol = solve(problem, saveat=1.0)

    # return sol
end

# sol = solve_!();

println(Dates.format(now(), "HH:MM:SS"))


# do_plotting(settings, sol)
using BenchmarkTools

@btime solve_!()