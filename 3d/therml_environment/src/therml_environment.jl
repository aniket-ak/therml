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
int(x) = floor(Int, x)

function nodes_from_vertices(vertices)
    # find a large multiplier to convert the float vertices to int
    large_multiplier = scaling_factor(vertices)

    # scale the vertices
    scaled_up_vertices = [round(Int, i*large_multiplier) for i in vertices]

    # compute the differences in the vertices. In some sense these differences are the edges of the model
    differences = compute_differences(scaled_up_vertices)

    # Find the GCD of the these differences. We cannot enforce a dx in our case. Enforcing dx may not capture the edges (differences) correctly
    gcd_ = find_gcd_of_list(differences)

    # Scale the dx back to original dimension
    discretization = gcd_/large_multiplier

    # Include ghost cells on extreme ends and discretize the domain and finally collect nodes
    nodes = range(vertices[1]-discretization, vertices[end]+discretization, step=discretization) |> collect
    return nodes, discretization
end

function get_vertices_from_bodies(settings)
    bodies_ = settings[:model]["bodies"]

    x_start_solder = -bodies_["solder"]["size"]["X"]/2
    x_end_solder = bodies_["solder"]["size"]["X"]/2

    x_start_substrate = -bodies_["substrate"]["size"]["X"]/2
    x_end_substrate = bodies_["substrate"]["size"]["X"]/2
    
    x_start_bumps = -bodies_["bumps"]["size"]["X"]/2
    x_end_bumps = bodies_["bumps"]["size"]["X"]/2

    x_start_underfill = -bodies_["underfill"]["size"]["X"]/2
    x_end_underfill = bodies_["underfill"]["size"]["X"]/2

    x_start_die = -bodies_["die"]["size"]["X"]/2
    x_end_die = bodies_["die"]["size"]["X"]/2
    
    x_start_mold = -bodies_["mold"]["size"]["X"]/2
    x_end_mold = bodies_["mold"]["size"]["X"]/2

    y_start_solder = -bodies_["solder"]["size"]["Y"]/2
    y_end_solder = bodies_["solder"]["size"]["Y"]/2

    y_start_substrate = -bodies_["substrate"]["size"]["Y"]/2
    y_end_substrate = bodies_["substrate"]["size"]["Y"]/2

    y_start_bumps = -bodies_["bumps"]["size"]["Y"]/2
    y_end_bumps = bodies_["bumps"]["size"]["Y"]/2

    y_start_underfill = -bodies_["underfill"]["size"]["Y"]/2
    y_end_underfill = bodies_["underfill"]["size"]["Y"]/2

    y_start_die = -bodies_["die"]["size"]["Y"]/2
    y_end_die = bodies_["die"]["size"]["Y"]/2

    y_start_mold = -bodies_["mold"]["size"]["Y"]/2
    y_end_mold = bodies_["mold"]["size"]["Y"]/2

    z_start_solder = 0
    z_end_solder = bodies_["solder"]["size"]["Z"]

    z_start_substrate = z_end_solder
    z_end_substrate = z_end_solder + bodies_["substrate"]["size"]["Z"]

    z_start_mold = z_end_substrate
    z_end_mold = z_start_mold + bodies_["mold"]["size"]["Z"]

    z_start_bumps = z_start_mold
    z_end_bumps = z_start_mold + bodies_["bumps"]["size"]["Z"]

    z_start_underfill = z_start_mold
    z_end_underfill = z_start_underfill + bodies_["underfill"]["size"]["Z"]

    z_start_die = z_end_bumps
    z_end_die = z_start_die + bodies_["die"]["size"]["Z"]

    X = [x_start_solder, x_end_solder, x_start_substrate, x_end_substrate, x_start_bumps, x_end_bumps, x_start_underfill, x_end_underfill, x_start_die, x_end_die, x_start_mold, x_end_mold]
    Y = [y_start_solder, y_end_solder, y_start_substrate, y_end_substrate, y_start_bumps, y_end_bumps, y_start_underfill, y_end_underfill, y_start_die, y_end_die, y_start_mold, y_end_mold]
    Z = [z_start_solder, z_end_solder, z_start_substrate, z_end_substrate, z_start_mold, z_end_mold, z_start_bumps, z_end_bumps, z_start_underfill, z_end_underfill, z_start_die, z_end_die]

    return X, Y, Z
end

function generate_mesh(settings)
    X,Y,Z = get_vertices_from_bodies(settings)
    X, Y, Z = sort(unique(X)), sort(unique(Y)), sort(unique(Z))

    X_nodes, dx = nodes_from_vertices(X)
    Y_nodes, dy = nodes_from_vertices(Y)
    Z_nodes, dz = nodes_from_vertices(Z)

    num_cells_x, num_cells_y, num_cells_z = (size(X_nodes)[1]-1, size(Y_nodes)[1]-1, size(Z_nodes)[1]-1)

    return (num_cells_x, num_cells_y, num_cells_z), (dx, dy, dz)
end

function cell_belongs_to_bbox(cell, bbox)
    x_,y_,z_ = cell
    ((x_s, x_e), (y_s, y_e), (z_s, z_e)) = bbox

    if x_ < x_e && x_ > x_s && y_ < y_e && y_ > y_s && z_ < z_e && z_ > z_s
        return true
    else
        return false
    end
end

function get_k_by_rho_cp(num_cells, discretization, settings)
    Nx,Ny,Nz = num_cells
    dx,dy,dz = discretization
    k_by_rho_cp = zeros(Nx, Ny, Nz)
    X,Y,Z = get_vertices_from_bodies(settings)

    x_start_solder, x_end_solder, x_start_substrate, x_end_substrate, x_start_bumps, x_end_bumps, x_start_underfill, x_end_underfill, x_start_die, x_end_die, x_start_mold, x_end_mold = X
    y_start_solder, y_end_solder, y_start_substrate, y_end_substrate, y_start_bumps, y_end_bumps, y_start_underfill, y_end_underfill, y_start_die, y_end_die, y_start_mold, y_end_mold = Y
    z_start_solder, z_end_solder, z_start_substrate, z_end_substrate, z_start_mold, z_end_mold, z_start_bumps, z_end_bumps, z_start_underfill, z_end_underfill, z_start_die, z_end_die = Z

    min_x_domain = minimum(X)
    min_y_domain = minimum(Y)
    min_z_domain = minimum(Z)

    for k in range(1, Nz)
        for j in range(1, Ny)
            for i in range(1, Nx)
                x_ = min_x_domain + i * dx/2
                y_ = min_y_domain + j * dy/2
                z_ = min_z_domain + k * dz/2
                if cell_belongs_to_bbox((x_,y_,z_), ((x_start_solder, x_end_solder), (y_start_solder, y_end_solder), (z_start_solder, z_end_solder)))
                    k_by_rho_cp[i,j,k] = settings[:model]["bodies"]["solder"]["material"]["k"] / (settings[:model]["bodies"]["solder"]["material"]["rho"] * settings[:model]["bodies"]["solder"]["material"]["cp"])
                elseif cell_belongs_to_bbox((x_,y_,z_), ((x_start_substrate, x_end_substrate), (y_start_substrate, y_end_substrate), (z_start_substrate, z_end_substrate)))
                    k_by_rho_cp[i,j,k] = settings[:model]["bodies"]["substrate"]["material"]["k"] / (settings[:model]["bodies"]["substrate"]["material"]["rho"] * settings[:model]["bodies"]["substrate"]["material"]["cp"])
                elseif cell_belongs_to_bbox((x_,y_,z_), ((x_start_bumps, x_end_bumps), (y_start_bumps, y_end_bumps), (z_start_bumps, z_end_bumps)))
                    k_by_rho_cp[i,j,k] = settings[:model]["bodies"]["bumps"]["material"]["k"] / (settings[:model]["bodies"]["bumps"]["material"]["rho"] * settings[:model]["bodies"]["bumps"]["material"]["cp"])
                elseif cell_belongs_to_bbox((x_,y_,z_), ((x_start_mold, x_end_mold), (y_start_mold, y_end_mold), (z_start_mold, z_end_mold)))
                    k_by_rho_cp[i,j,k] = settings[:model]["bodies"]["mold"]["material"]["k"] / (settings[:model]["bodies"]["mold"]["material"]["rho"] * settings[:model]["bodies"]["mold"]["material"]["cp"])
                elseif cell_belongs_to_bbox((x_,y_,z_), ((x_start_underfill, x_end_underfill), (y_start_underfill, y_end_underfill), (z_start_underfill, z_end_underfill)))
                    k_by_rho_cp[i,j,k] = settings[:model]["bodies"]["underfill"]["material"]["k"] / (settings[:model]["bodies"]["underfill"]["material"]["rho"] * settings[:model]["bodies"]["underfill"]["material"]["cp"])
                elseif cell_belongs_to_bbox((x_,y_,z_), ((x_start_die, x_end_die), (y_start_die, y_end_die), (z_start_die, z_end_die)))
                    k_by_rho_cp[i,j,k] = settings[:model]["bodies"]["die"]["material"]["k"] / (settings[:model]["bodies"]["die"]["material"]["rho"] * settings[:model]["bodies"]["die"]["material"]["cp"])
                end
            end
        end
    end
    
    return k_by_rho_cp
end

function initialize_domain!(u0, settings)
    temperature_base_units = settings[:units]["temperature"]
    T_initial = convert_units("temperature", settings[:IC], temperature_base_units, "K")
    u0 = u0 .+ T_initial
    return u0
end

function define_volume_sources(working_dir, power_file, settings, Nx, Ny, Nz)
    source = zeros(Nx,Ny,Nz)
    # source[10:20, 50:70, 1:5] .= 1

    x_length = settings[:model]["bodies"]["die"]["size"]["X"]
    y_length = settings[:model]["bodies"]["die"]["size"]["Y"]
    z_length = settings[:model]["bodies"]["die"]["size"]["Z"]

    x_normalized = x_length/z_length
    y_normalized = y_length/z_length
    z_normalized = z_length/z_length

    x_mesh = range(0, stop = x_normalized, length = Nx+1)
    y_mesh = range(0, stop = y_normalized, length = Ny+1)
    z_mesh = range(0, stop = z_normalized, length = Nz+1)

    cell_centers_x = range(step(x_mesh)/2, stop = x_normalized-step(x_mesh)/2, length = Nx)
    cell_centers_y = range(step(y_mesh)/2, stop = x_normalized-step(y_mesh)/2, length = Ny)

    heat_values = read_csv(joinpath(working_dir,power_file))
    source_Nx, source_Ny = size(heat_values)
    source_x = range(0, x_normalized, source_Nx)
    source_y = range(0, y_normalized, source_Ny)

    interpolated_heat = interpolate_(source_x,source_y,heat_values,cell_centers_x,cell_centers_y)
    source[:,:,2] = interpolated_heat
    return source
end

function do_plotting(sol, sol_wd, settings, scenario_name, interpolation)
    result_ = convert_units("temperature", sol[end, 2:end-1,2:end-1,1], "K", "C")'
    # result_ = sol[end][2:end-1,2:end-1,2]
    (delta_x, delta_y, delta_z), (Nx,Ny,Nz), (x_mesh, y_mesh, z_mesh) = create_mesh(settings)
    z_length = settings[:model]["bodies"]["die"]["size"]["Z"]
    
    if interpolation == true
        p = plot(
            contour(
                x=x_mesh*z_length,
                y=y_mesh*z_length,
                z=result_, 
                colorscale="Jet",
                colorbar=attr(width=80, height=80, automargin=true,title="Temperature", 
                            titleside="right",
                            titlefont=attr(size=14,family="Arial, sans-serif"))
                ),
            )
    else
        p = plot(
            heatmap(
                x=x_mesh*z_length,
                y=y_mesh*z_length,
                z=result_, 
                colorscale="Jet",
                colorbar=attr(width=80, height=80, automargin=true,title="Temperature", 
                            titleside="right",
                            titlefont=attr(size=14,family="Arial, sans-serif"))
                ),
            )
    end

    open(joinpath(sol_wd,scenario_name*"__plot.html"), "w") do io
        PlotlyBase.to_html(io, p.plot)
    end
end

function convert_units(quantity, value, from_units, to_units)
    if quantity == "temperature"
        if from_units == "C" && to_units == "K"
            return_value = value .+ 273.15
            return return_value
        elseif from_units == "K" && to_units == "C"
            return_value = value .- 273.15
            return return_value
        else
            return value
        end
    else
        return value
    end
end

function apply_bc(u, settings, delta_x, delta_y, delta_z)
    k = settings[:model]["bodies"]["die"]["material"]["k"]
    rho = settings[:model]["bodies"]["die"]["material"]["rho"]
    cp = settings[:model]["bodies"]["die"]["material"]["cp"]

    # At X-
    side = "X-"
    if settings[:BC]["X-"]["type"] == "constant_T"
        const_temp = convert_units("temperature", settings[:BC]["X-"]["value"]["value"], settings[:units]["temperature"], "K") 
        u[1,:,:] .= 2*const_temp .- u[2,:,:]
    elseif settings[:BC]["X-"]["type"] == "insulated"
        u[1,:,:] .= u[2,:,:]
    elseif settings[:BC]["X-"]["type"] == "const_flux"
        u[1,:,:] .= u[2,:,:] .+ settings[:BC]["X-"]["value"]["value"] * (delta_x/k)
    elseif settings[:BC]["X-"]["type"] == "HTC"
        h = settings[:BC]["X-"]["value"]["value"]
        t_amb = convert_units("temperature", settings[:BC]["X-"]["value"]["t_amb"], settings[:units]["temperature"], "K")
        u[1,:,:] .= u[2,:,:]*(1 - delta_x * h/k) .+ (h*delta_x/k) * t_amb
    end

    # At X+
    side = "X+"
    if settings[:BC]["X+"]["type"] == "constant_T"
        const_temp = convert_units("temperature", settings[:BC]["X+"]["value"]["value"], settings[:units]["temperature"], "K") 
        u[end,:,:] .= 2*settings[:BC]["X+"]["value"] .- u[end-1,:,:]
    elseif settings[:BC]["X+"]["type"] == "insulated"
        u[end,:,:] .= u[end-1,:,:]
    elseif settings[:BC]["X+"]["type"] == "const_flux"
        u[end,:,:] .= u[end-1,:,:] .+ settings[:BC]["X+"]["value"]["value"] * (delta_x/k)
    elseif settings[:BC]["X+"]["type"] == "HTC"
        h = settings[:BC]["X+"]["value"]["value"]
        t_amb = convert_units("temperature", settings[:BC]["X+"]["value"]["t_amb"], settings[:units]["temperature"], "K")
        u[end,:,:] .= u[end-1,:,:]*(1 - delta_x * h/k) .+ (h*delta_x/k) * t_amb
    end

    # Y-
    side = "Y-"
    if settings[:BC]["Y-"]["type"] == "constant_T"
        const_temp = convert_units("temperature", settings[:BC]["Y-"]["value"]["value"], settings[:units]["temperature"], "K") 
        u[:,1,:] .= 2*const_temp .- u[:,2,:]
    elseif settings[:BC]["Y-"]["type"] == "insulated"
        u[:,1,:] .= u[:,2,:]
    elseif settings[:BC]["Y-"]["type"] == "const_flux"
        u[:,1,:] .= u[:,2,:] .+ settings[:BC]["Y-"]["value"]["value"] * (delta_y/k)
    elseif settings[:BC]["Y-"]["type"] == "HTC"
        h = settings[:BC]["Y-"]["value"]["value"]
        t_amb = convert_units("temperature", settings[:BC]["Y-"]["value"]["t_amb"], settings[:units]["temperature"], "K")
        u[:,1,:] .= u[:,2,:]*(1 - delta_y * h/k) .+ (h*delta_y/k) * t_amb
    end
    
    # Y+
    side = "Y+"
    if settings[:BC]["Y+"]["type"] == "constant_T"
        const_temp = convert_units("temperature", settings[:BC]["Y+"]["value"]["value"], settings[:units]["temperature"], "K") 
        u[:,end,:] .= 2*const_temp .- u[:,end-1,:]
    elseif settings[:BC]["Y+"]["type"] == "insulated"
        u[:,end,:] .= u[:,end-1,:]
    elseif settings[:BC]["Y+"]["type"] == "const_flux"
        u[:,end,:] .= u[:,end-1,:] .+ settings[:BC]["Y+"]["value"]["value"] * (delta_y/k)
    elseif settings[:BC]["Y+"]["type"] == "HTC"
        h = settings[:BC]["Y+"]["value"]["value"]
        t_amb = convert_units("temperature", settings[:BC]["Y+"]["value"]["t_amb"], settings[:units]["temperature"], "K")
        u[:,end,:] .= u[:,end-1,:]*(1 - delta_y * h/k) .+ (h*delta_y/k) * t_amb
    end

    # Z-
    side = "Z-"
    if settings[:BC]["Z-"]["type"] == "constant_T"
        const_temp = convert_units("temperature", settings[:BC]["Z-"]["value"]["value"], settings[:units]["temperature"], "K") 
        u[:,:,1] .= 2*const_temp .- u[:,:,2]
    elseif settings[:BC]["Z-"]["type"] == "insulated"
        u[:,:,1] .= u[:,:,2]
    elseif settings[:BC]["Z-"]["type"] == "const_flux"
        u[:,:,1] .= u[:,:,2] .+ settings[:BC]["Z-"]["value"]["value"] * (delta_z/k)
    elseif settings[:BC]["Z-"]["type"] == "HTC"
        h = settings[:BC]["Z-"]["value"]["value"]
        t_amb = convert_units("temperature", settings[:BC]["Z-"]["value"]["t_amb"], settings[:units]["temperature"], "K")
        u[:,:,1] .= u[:,:,2]*(1 - delta_z * h/k) .+ (h*delta_z/k) * t_amb
    end

    # Z+
    side = "Z+"
    if settings[:BC]["Z+"]["type"] == "constant_T"
        const_temp = convert_units("temperature", settings[:BC]["Z+"]["value"]["value"], settings[:units]["temperature"], "K") 
        u[:,:,end] .= 2*const_temp .- u[:,:,end-1]
    elseif settings[:BC]["Z+"]["type"] == "insulated"
        u[:,:,end] .= u[:,:,end-1]
    elseif settings[:BC]["Z+"]["type"] == "const_flux"
        u[:,:,end] .= u[:,:,end-1] .+ settings[:BC]["Z+"]["value"]["value"] * (delta_z/k)
    elseif settings[:BC]["Z+"]["type"] == "HTC"
        h = settings[:BC]["Z+"]["value"]["value"]
        t_amb = convert_units("temperature", settings[:BC]["Z+"]["value"]["t_amb"], settings[:units]["temperature"], "K")
        u[:,:,end] .= u[:,:,end-1]*(1 - delta_z * h/k) .+ (h*delta_z/k) * t_amb
    end
end

function conduction_3d_loop!(du, u, p, t)
    k_by_rho_cp, (dx, dy, dz), settings = p

    Nx, Ny, Nz = size(u)

    source = define_volume_sources(working_dir, power_file,settings,Nx,Ny,Nz)

    Threads.@threads for k_ in range(2, Nz-1)
        Threads.@threads for j_ in range(2, Ny-1)
            Threads.@threads for i_ in range(2, Nx-1)
                ip1, im1, jp1, jm1, kp1, km1 = i_ + 1, i_ - 1, j_ + 1, j_ - 1, k_+1, k_-1
                du[i_, j_, k_] =  k_by_rho_cp[i_, j_, k_]/(dx^2) * (u[im1, j_, k_] + u[ip1, j_, k_] - 2u[i_,j_, k_])+
                                    k_by_rho_cp[i_, j_, k_]/(dy^2) * (u[i_, jp1, k_] + u[i_, jm1, k_] - 2u[i_, j_, k_]) + 
                                    k_by_rho_cp[i_, j_, k_]/(dz^2) * (u[i_, j_, kp1] + u[i_, j_, km1] - 2u[i_, j_, k_]) +  
                                    source[i_, j_, k_]
            end
        end
    end

    apply_bc(u, settings, 0.1, 0.1, 0.1)

end

function configure_problem_klu!(u0,p)
    du0 = zeros(size(u0))
    jac_sparsity = Symbolics.jacobian_sparsity((du, u) -> conduction_3d_loop!(du, u, p, 0.0), du0, u0)

    f = ODEFunction(conduction_3d_loop!; jac_prototype = float.(jac_sparsity))

    start_time = settings[:start_time]
    end_time = settings[:end_time]

    prob_conduction_3d_sparse = ODEProblem(f, u0, (start_time, end_time), p)

    algorithm = KenCarp47(linsolve = KLUFactorization())

    return prob_conduction_3d_sparse, algorithm
end

function configure_problem_ode!(u0,settings)
    start_time = settings[:start_time]
    end_time = settings[:end_time]

    # problem = ODEProblem(conduction_3d_loop!, u0, (start_time, end_time), p)
    problem = ODEProblem(conduction_3d_loop!, u0, (start_time, end_time))
    algorithm = ""

    return problem, algorithm
end

function read_csv(file_name)
    data = CSV.read(file_name, DataFrame)
    mat = Matrix(data)
    return mat
end

function read_excel(file_name)
    xf = XLSX.readxlsx(file_name)
    XLSX.sheetnames(xf)
    sh = xf["Sheet1"]
    return sh[:]

end

function interpolate_(x,y,z,new_x,new_y)
    itp = LinearInterpolation((x, y), z)
    new_z = [itp(x,y) for y in new_y, x in new_x]
    return new_z
end

function save_fields(sol, t_, sol_wd, scenario_name)
    h5open(joinpath(sol_wd,scenario_name*"__solution.sol"), "w") do file
        g = HDF5.create_group(file, "solution")

        for (i,t) in enumerate(t_)
            g[string(t)] = sol[i,:,:,:]
        end
    end
end

function solve_(working_dir, power_file, settings, progress_file_name)

    progress_file = open(progress_file_name, "w")
    close(progress_file)

    (Nx,Ny,Nz), (dx, dy, dz) = generate_mesh(settings)
    k_by_rho_cp = get_k_by_rho_cp((Nx,Ny,Nz), (dx, dy, dz), settings)

    println("Mesh details - Total cells :"(Nx-2)*(Ny-2)*(Nz-2)/1e6, " k", "with ", Nx-2, ",", Ny-2, ", and ", Nz-2, " in X,Y and Z")

    u0 = zeros(Nx,Ny,Nz)
    
    # source = define_volume_sources(working_dir, power_file,settings,Nx,Ny,Nz)

    p = (k_by_rho_cp, (dx,dy,dz), settings)

    # p = NamedTuple([pair for pair in settings])
    u0 = initialize_domain!(u0, settings)

    start_time = settings[:start_time]
    end_time = settings[:end_time]

    problem = ODEProblem(conduction_3d_loop!, u0, (start_time, end_time), p)

    dt = settings[:dt]
    n_steps = round(Int, settings[:end_time]/dt)

    sol = zeros((n_steps, size(u0)...))
    t_ = zeros(n_steps)

    integrator = init(problem; reltol=1e-4, abstol=1e-4, maxiters=1000)
    for i in range(1,n_steps)
        step!(integrator, dt, true)
        t,u = integrator.t, integrator.u
        sol[i,:,:,:] = u
        t_[i] = t

        progress_file = open(progress_file_name, "a")
        write(progress_file, string(i/n_steps*100))
        write(progress_file,"\n")
        close(progress_file)
    end

    return sol, t_
end

function real_main()
    working_dir = ARGS[1]
    scenario_name = ARGS[2]
    run_name = ARGS[3]

    println("working dir : ", working_dir)
    println("scenario name : ", scenario_name)
    println("run name : ", run_name)

    if run_name == ""
        run_name = now()
    end

    run_wd = joinpath(working_dir, run_name)
    log_wd = joinpath(run_wd,"Logs")
    sol_wd = joinpath(run_wd,"Solution")
    temp_wd = joinpath(run_wd,"Temp")

    progress_file_name = joinpath(temp_wd, scenario_name*"__progress.txt")
    # progress_file = open(progress_file_name, "w")
    # global_logger(TerminalLogger(progress_file))

    #include("/root/therml/therml/3d/therml_environment/src/unctions_fvm_3d.jl")

    f = open(joinpath(working_dir,"settings.json"), "r")
    s1 = JSON.parse(f)

    settings = Dict{Symbol, Any}()
    for (k,v) in s1
        settings[Symbol(k)] = v
    end

    date_start = now()
    println("Start solution at : ", Dates.format(date_start, "HH:MM:SS"))
    sol,t_ = solve_(working_dir, scenario_name, settings, progress_file_name);
    
    # do_plotting(sol, sol_wd, settings, scenario_name, false);
    
    save_fields(sol, t_, sol_wd, scenario_name);
    
    date_end = now()
    println("End solution at : ", Dates.format(date_end, "HH:MM:SS"))
    println("Time taken for solution : ", (date_end-date_start)/Millisecond(1000), " [s]")

    flush(stdout)
    exit()
end

end
