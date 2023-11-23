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
using Sundials
using Nettle

function get_serial_number()
    if Sys.iswindows()
        command = `wmic bios get serialnumber`
    elseif Sys.isapple()
        command = pipeline(`ioreg -l`, `grep IOPlatformSerialNumber`)
    elseif Sys.islinux()
        command = `cat /var/lib/dbus/machine-id`
    else
        return "Unsupported OS"
    end

    try
        if Sys.isapple()
            # For macOS, the output processing is a bit different
            serial = read(pipeline(command, `awk '/IOPlatformSerialNumber/ {print $4}'`, `tr -d '"'`), String)
        else
            serial = read(pipeline(command, `tr -d '\r'`, `awk NR==2`), String)
        end
        return strip(serial)
    catch e
        return string(e)
    end
end

function read_license_file(license_file_location)
    return "34d7d82e6b4d3870a9ccc3a18111b6442bae42cd612c07bea075cb4e60b30685"
end

function check_license()
    serial_num = get_serial_number()
    key_ = "MarlinSim_" * string(serial_num)
    license_key = hexdigest("sha256", key_)

    license_file_location = ""
    if haskey(ENV, "THERML_ENVIRONMENT")
        license_file_location = ENV["THERML_ENVIRONMENT"]
    end
    given_license = read_license_file(license_file_location)

    if license_key == given_license
        return true
    else
        return false
    end
end

function julia_main()::Cint
    try
        if check_license()
            real_main()
        else
            println("-------License file not found-------")
        end
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end
int(x) = floor(Int, x)

function write_mesh_and_conductivity_to_file(k, X_nodes, Y_nodes, Z_nodes, path)
    k_file = joinpath(path, "conductivity.jld")
    save(k_file, "k", k)

    f = open(joinpath(path, "x_mesh.txt"), "w")
    for i in X_nodes
        write(f, string(i))
        write(f,"\n")
    end
    close(f)

    f = open(joinpath(path, "y_mesh.txt"), "w")
    for i in Y_nodes
        write(f, string(i))
        write(f,"\n")
    end
    close(f)

    f = open(joinpath(path, "z_mesh.txt"), "w")
    for i in Z_nodes
        write(f, string(i))
        write(f,"\n")
    end
    close(f)
end

function round_off(list_)
    list_out = [round(i, digits=6) for i in list_]
    return list_out
end

function get_vertices_from_bodies(settings)
    bodies_ = settings["model"]["bodies"]

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

    X = round_off(X)
    Y = round_off(Y)
    Z = round_off(Z)

    return X, Y, Z
end

function generate_mesh(settings)
    X,Y,Z = get_vertices_from_bodies(settings)
    bodies_ = settings["model"]["bodies"]
    
    x_start_solder, x_end_solder, x_start_substrate, x_end_substrate, x_start_bumps, x_end_bumps, x_start_underfill, x_end_underfill, x_start_die, x_end_die, x_start_mold, x_end_mold = X
    y_start_solder, y_end_solder, y_start_substrate, y_end_substrate, y_start_bumps, y_end_bumps, y_start_underfill, y_end_underfill, y_start_die, y_end_die, y_start_mold, y_end_mold = Y
    z_start_solder, z_end_solder, z_start_substrate, z_end_substrate, z_start_mold, z_end_mold, z_start_bumps, z_end_bumps, z_start_underfill, z_end_underfill, z_start_die, z_end_die = Z

    X, Y, Z = sort(unique(X)), sort(unique(Y)), sort(unique(Z))

    X_mesh, Y_mesh, Z_mesh = [], [], []

    for (i, mesh_point) in enumerate(X)
        if i > 1
            mid_point = (mesh_point+X[i-1])/2

            largest_precedence, sizing = 0, 0
            if mid_point > x_start_solder && mid_point < x_end_solder
                if largest_precedence < bodies_["solder"]["mesh"]["precedence"]
                    largest_precedence = bodies_["solder"]["mesh"]["precedence"]
                    sizing = bodies_["solder"]["mesh"]["size"]["dx"]
                    # println("Solder  - ",i, "\t", mid_point, "\t", largest_precedence, "\t", sizing)
                end
            end

            if mid_point > x_start_substrate && mid_point < x_end_substrate
                if largest_precedence < bodies_["substrate"]["mesh"]["precedence"]
                    largest_precedence = bodies_["substrate"]["mesh"]["precedence"]
                    sizing = bodies_["substrate"]["mesh"]["size"]["dx"]
                    # println("substrate  - ",i, "\t", mid_point, "\t", largest_precedence, "\t", sizing)
                end
            end

            if mid_point > x_start_bumps && mid_point < x_end_bumps
                if largest_precedence < bodies_["bumps"]["mesh"]["precedence"]
                    largest_precedence = bodies_["bumps"]["mesh"]["precedence"]
                    sizing = bodies_["bumps"]["mesh"]["size"]["dx"]
                    # println("bumps  - ",i, "\t", mid_point, "\t", largest_precedence, "\t", sizing)
                end
            end

            if mid_point > x_start_underfill && mid_point < x_end_underfill
                if largest_precedence < bodies_["underfill"]["mesh"]["precedence"]
                    largest_precedence = bodies_["underfill"]["mesh"]["precedence"]
                    sizing = bodies_["underfill"]["mesh"]["size"]["dx"]
                    # println("underfill  - ",i, "\t", mid_point, "\t", largest_precedence, "\t", sizing)
                end
            end

            if mid_point > x_start_die && mid_point < x_end_die
                if largest_precedence < bodies_["die"]["mesh"]["precedence"]
                    largest_precedence = bodies_["die"]["mesh"]["precedence"]
                    sizing = bodies_["die"]["mesh"]["size"]["dx"]
                    # println("die  - ",i, "\t", mid_point, "\t", largest_precedence, "\t", sizing)
                end
            end

            if mid_point > x_start_mold && mid_point < x_end_mold
                if largest_precedence < bodies_["mold"]["mesh"]["precedence"]
                    largest_precedence = bodies_["mold"]["mesh"]["precedence"]
                    sizing = bodies_["mold"]["mesh"]["size"]["dx"]
                    # println("mold  - ",i, "\t", mid_point, "\t", largest_precedence, "\t", sizing)
                end
            end

            min_num_nodes = max(2, round(Int, abs(mesh_point-X[i-1])/sizing)) + 1
            # println(X[i-1],"\t", mesh_point,"\t", min_num_nodes, "\t", sizing)
            mesh_ = range(X[i-1], mesh_point, min_num_nodes)|>collect
            # println(X[i-1], "\t", mesh_point, "\t", sizing, "\t", min_num_nodes, "\t", mesh_[end])
            append!(X_mesh, mesh_)
        end

    end

    # println("Mesh for X done")

    for (i, mesh_point) in enumerate(Y)
        if i > 1
            mid_point = (mesh_point+Y[i-1])/2

            largest_precedence, sizing = 0, 0
            if mid_point > y_start_solder && mid_point < y_end_solder
                if largest_precedence < bodies_["solder"]["mesh"]["precedence"]
                    largest_precedence = bodies_["solder"]["mesh"]["precedence"]
                    sizing = bodies_["solder"]["mesh"]["size"]["dy"]
                    # println("Solder  - ",i, "\t", mid_point, "\t", largest_precedence, "\t", sizing)
                end
            end

            if mid_point > y_start_substrate && mid_point < y_end_substrate
                if largest_precedence < bodies_["substrate"]["mesh"]["precedence"]
                    largest_precedence = bodies_["substrate"]["mesh"]["precedence"]
                    sizing = bodies_["substrate"]["mesh"]["size"]["dy"]
                    # println("substrate  - ",i, "\t", mid_point, "\t", largest_precedence, "\t", sizing)
                end
            end

            if mid_point > y_start_bumps && mid_point < y_end_bumps
                if largest_precedence < bodies_["bumps"]["mesh"]["precedence"]
                    largest_precedence = bodies_["bumps"]["mesh"]["precedence"]
                    sizing = bodies_["bumps"]["mesh"]["size"]["dy"]
                    # println("bumps  - ",i, "\t", mid_point, "\t", largest_precedence, "\t", sizing)
                end
            end

            if mid_point > y_start_underfill && mid_point < y_end_underfill
                if largest_precedence < bodies_["underfill"]["mesh"]["precedence"]
                    largest_precedence = bodies_["underfill"]["mesh"]["precedence"]
                    sizing = bodies_["underfill"]["mesh"]["size"]["dy"]
                    # println("underfill  - ",i, "\t", mid_point, "\t", largest_precedence, "\t", sizing)
                end
            end

            if mid_point > y_start_die && mid_point < y_end_die
                if largest_precedence < bodies_["die"]["mesh"]["precedence"]
                    largest_precedence = bodies_["die"]["mesh"]["precedence"]
                    sizing = bodies_["die"]["mesh"]["size"]["dy"]
                    # println("die  - ",i, "\t", mid_point, "\t", largest_precedence, "\t", sizing)
                end
            end

            if mid_point > y_start_mold && mid_point < y_end_mold
                if largest_precedence < bodies_["mold"]["mesh"]["precedence"]
                    largest_precedence = bodies_["mold"]["mesh"]["precedence"]
                    sizing = bodies_["mold"]["mesh"]["size"]["dy"]
                    # println("mold  - ",i, "\t", mid_point, "\t", largest_precedence, "\t", sizing)
                end
            end

            min_num_nodes = max(2, round(Int, abs(mesh_point-Y[i-1])/sizing)) + 1
            mesh_ = range(Y[i-1], mesh_point, min_num_nodes)|>collect
            # println(X[i-1], "\t", mesh_point, "\t", sizing, "\t", min_num_nodes, "\t", mesh_[end])
            append!(Y_mesh, mesh_)
        end

    end

    # println("Mesh for Y done")

    for (i, mesh_point) in enumerate(Z)
        if i > 1
            mid_point = (mesh_point+Z[i-1])/2

            largest_precedence, sizing = 0, 0
            if mid_point > z_start_solder && mid_point < z_end_solder
                if largest_precedence < bodies_["solder"]["mesh"]["precedence"]
                    largest_precedence = bodies_["solder"]["mesh"]["precedence"]
                    sizing = bodies_["solder"]["mesh"]["size"]["dz"]
                    # println("Solder  - ",i, "\t", mid_point, "\t", largest_precedence, "\t", sizing)
                end
            end

            if mid_point > z_start_substrate && mid_point < z_end_substrate
                if largest_precedence < bodies_["substrate"]["mesh"]["precedence"]
                    largest_precedence = bodies_["substrate"]["mesh"]["precedence"]
                    sizing = bodies_["substrate"]["mesh"]["size"]["dz"]
                    # println("substrate  - ",i, "\t", mid_point, "\t", largest_precedence, "\t", sizing)
                end
            end

            if mid_point > z_start_bumps && mid_point < z_end_bumps
                if largest_precedence < bodies_["bumps"]["mesh"]["precedence"]
                    largest_precedence = bodies_["bumps"]["mesh"]["precedence"]
                    sizing = bodies_["bumps"]["mesh"]["size"]["dz"]
                    # println("bumps  - ",i, "\t", mid_point, "\t", largest_precedence, "\t", sizing)
                end
            end

            if mid_point > z_start_underfill && mid_point < z_end_underfill
                if largest_precedence < bodies_["underfill"]["mesh"]["precedence"]
                    largest_precedence = bodies_["underfill"]["mesh"]["precedence"]
                    sizing = bodies_["underfill"]["mesh"]["size"]["dz"]
                    # println("underfill  - ",i, "\t", mid_point, "\t", largest_precedence, "\t", sizing)
                end
            end

            if mid_point > z_start_die && mid_point < z_end_die
                if largest_precedence < bodies_["die"]["mesh"]["precedence"]
                    largest_precedence = bodies_["die"]["mesh"]["precedence"]
                    sizing = bodies_["die"]["mesh"]["size"]["dz"]
                    # println("die  - ",i, "\t", mid_point, "\t", largest_precedence, "\t", sizing)
                end
            end

            if mid_point > z_start_mold && mid_point < z_end_mold
                if largest_precedence < bodies_["mold"]["mesh"]["precedence"]
                    largest_precedence = bodies_["mold"]["mesh"]["precedence"]
                    sizing = bodies_["mold"]["mesh"]["size"]["dz"]
                    # println("mold  - ",i, "\t", mid_point, "\t", largest_precedence, "\t", sizing)
                end
            end

            min_num_nodes = max(2, round(Int, abs(mesh_point-Z[i-1])/sizing)) + 1
            # println(Z[i-1], "\t", mesh_point, "\t", sizing, "\t", min_num_nodes)
            mesh_ = range(Z[i-1], mesh_point, min_num_nodes)|>collect
            append!(Z_mesh, mesh_)
        end

    end

    X_mesh, Y_mesh, Z_mesh = sort(unique(X_mesh)), sort(unique(Y_mesh)), sort(unique(Z_mesh))
    ghost_x_min = X_mesh[1] - abs(X_mesh[2] - X_mesh[1])
    ghost_x_max = X_mesh[end] + abs(X_mesh[end] - X_mesh[end-1])
    ghost_y_min = Y_mesh[1] - abs(Y_mesh[2] - Y_mesh[1])
    ghost_y_max = Y_mesh[end] + abs(Y_mesh[end] - Y_mesh[end-1])
    ghost_z_min = Z_mesh[1] - abs(Z_mesh[2] - Z_mesh[1])
    ghost_z_max = Z_mesh[end] + abs(Z_mesh[end] - Z_mesh[end-1])

    append!(X_mesh, ghost_x_min)
    append!(X_mesh, ghost_x_max)
    append!(Y_mesh, ghost_y_min)
    append!(Y_mesh, ghost_y_max)
    append!(Z_mesh, ghost_z_min)
    append!(Z_mesh, ghost_z_max)
    
    X_mesh, Y_mesh, Z_mesh = sort(unique(X_mesh)), sort(unique(Y_mesh)), sort(unique(Z_mesh))

    num_cells_x, num_cells_y, num_cells_z = size(X_mesh)[1]-1, size(Y_mesh)[1]-1, size(Z_mesh)[1]-1

    return (X_mesh, Y_mesh, Z_mesh), (num_cells_x, num_cells_y, num_cells_z)
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

function get_k_by_rho_cp(nodes, (Nx, Ny, Nz), settings)
    nodes_x, nodes_y, nodes_z = nodes
    k_by_rho_cp = ones(Nx, Ny, Nz)
    k_ = ones(Nx, Ny, Nz)
    
    X,Y,Z = get_vertices_from_bodies(settings)

    x_start_solder, x_end_solder, x_start_substrate, x_end_substrate, x_start_bumps, x_end_bumps, x_start_underfill, x_end_underfill, x_start_die, x_end_die, x_start_mold, x_end_mold = X
    y_start_solder, y_end_solder, y_start_substrate, y_end_substrate, y_start_bumps, y_end_bumps, y_start_underfill, y_end_underfill, y_start_die, y_end_die, y_start_mold, y_end_mold = Y
    z_start_solder, z_end_solder, z_start_substrate, z_end_substrate, z_start_mold, z_end_mold, z_start_bumps, z_end_bumps, z_start_underfill, z_end_underfill, z_start_die, z_end_die = Z

    solder_precedence = settings["model"]["bodies"]["solder"]["mesh"]["precedence"]
    substrate_precedence = settings["model"]["bodies"]["substrate"]["mesh"]["precedence"]
    bumps_precedence = settings["model"]["bodies"]["bumps"]["mesh"]["precedence"]
    mold_precedence = settings["model"]["bodies"]["mold"]["mesh"]["precedence"]
    underfill_precedence = settings["model"]["bodies"]["underfill"]["mesh"]["precedence"]
    die_precedence = settings["model"]["bodies"]["die"]["mesh"]["precedence"]

    for k in range(2, Nz-1)
        for j in range(2, Ny-1)
            for i in range(2, Nx-1)
                x_ = nodes_x[i] + (nodes_x[i+1] - nodes_x[i])/2
                y_ = nodes_y[j] + (nodes_y[j+1] - nodes_y[j])/2
                z_ = nodes_z[k] + (nodes_z[k+1] - nodes_z[k])/2

                precedence = -1

                if cell_belongs_to_bbox((x_,y_,z_), ((x_start_solder, x_end_solder), (y_start_solder, y_end_solder), (z_start_solder, z_end_solder)))
                    if precedence < solder_precedence
                        k_by_rho_cp[i,j,k] = settings["model"]["bodies"]["solder"]["material"]["k"] / (settings["model"]["bodies"]["solder"]["material"]["rho"] * settings["model"]["bodies"]["solder"]["material"]["cp"])
                        k_[i,j,k] = settings["model"]["bodies"]["solder"]["material"]["k"]
                        precedence = solder_precedence
                    end
                    # println("Point falls in solder ", x_*1e3, ",", y_*1e3, ",", z_*1e3, " Solder bounds x", x_start_solder*1e3, ",", x_end_solder*1e3)
                end

                if cell_belongs_to_bbox((x_,y_,z_), ((x_start_substrate, x_end_substrate), (y_start_substrate, y_end_substrate), (z_start_substrate, z_end_substrate)))
                    if precedence < substrate_precedence
                        k_by_rho_cp[i,j,k] = settings["model"]["bodies"]["substrate"]["material"]["k"] / (settings["model"]["bodies"]["substrate"]["material"]["rho"] * settings["model"]["bodies"]["substrate"]["material"]["cp"])
                        k_[i,j,k] = settings["model"]["bodies"]["substrate"]["material"]["k"]
                        precedence = substrate_precedence
                    end
                end

                if cell_belongs_to_bbox((x_,y_,z_), ((x_start_bumps, x_end_bumps), (y_start_bumps, y_end_bumps), (z_start_bumps, z_end_bumps)))
                    if precedence < bumps_precedence
                        k_by_rho_cp[i,j,k] = settings["model"]["bodies"]["bumps"]["material"]["k"] / (settings["model"]["bodies"]["bumps"]["material"]["rho"] * settings["model"]["bodies"]["bumps"]["material"]["cp"])
                        k_[i,j,k] = settings["model"]["bodies"]["bumps"]["material"]["k"]
                        precedence = bumps_precedence
                    end
                end

                if cell_belongs_to_bbox((x_,y_,z_), ((x_start_mold, x_end_mold), (y_start_mold, y_end_mold), (z_start_mold, z_end_mold)))
                    if precedence < mold_precedence
                        k_by_rho_cp[i,j,k] = settings["model"]["bodies"]["mold"]["material"]["k"] / (settings["model"]["bodies"]["mold"]["material"]["rho"] * settings["model"]["bodies"]["mold"]["material"]["cp"])
                        k_[i,j,k] = settings["model"]["bodies"]["mold"]["material"]["k"]
                        precedence = mold_precedence
                    end
                end

                if cell_belongs_to_bbox((x_,y_,z_), ((x_start_underfill, x_end_underfill), (y_start_underfill, y_end_underfill), (z_start_underfill, z_end_underfill)))
                    if precedence < underfill_precedence
                        k_by_rho_cp[i,j,k] = settings["model"]["bodies"]["underfill"]["material"]["k"] / (settings["model"]["bodies"]["underfill"]["material"]["rho"] * settings["model"]["bodies"]["underfill"]["material"]["cp"])
                        k_[i,j,k] = settings["model"]["bodies"]["underfill"]["material"]["k"]
                        precedence = underfill_precedence
                    end
                end

                if cell_belongs_to_bbox((x_,y_,z_), ((x_start_die, x_end_die), (y_start_die, y_end_die), (z_start_die, z_end_die)))
                    if precedence < die_precedence
                        k_by_rho_cp[i,j,k] = settings["model"]["bodies"]["die"]["material"]["k"] / (settings["model"]["bodies"]["die"]["material"]["rho"] * settings["model"]["bodies"]["die"]["material"]["cp"])
                        k_[i,j,k] = settings["model"]["bodies"]["die"]["material"]["k"]
                        precedence = die_precedence
                    end
                end

            end
        end
    end
    
    return k_by_rho_cp, k_
end

function initialize_domain!(u0, settings)
    temperature_base_units = settings["units"]["temperature"]
    T_initial = convert_units("temperature", settings["IC"], temperature_base_units, "K")
    u0 = u0 .+ T_initial
    return u0
end

function define_volume_sources(power_file, settings, Nx, Ny, Nz, X_nodes, Y_nodes, Z_nodes)
    source = zeros(Nx,Ny,Nz)
    # source[10:20, 50:70, 1:5] .= 1

    x_length = settings["model"]["bodies"]["die"]["size"]["X"]
    y_length = settings["model"]["bodies"]["die"]["size"]["Y"]
    source_z_val = settings["model"]["bodies"]["solder"]["size"]["Z"] + 
                    settings["model"]["bodies"]["substrate"]["size"]["Z"] + 
                    settings["model"]["bodies"]["bumps"]["size"]["Z"] + 
                    settings["model"]["bodies"]["die"]["size"]["Z"]

    die_nodes_x = (findall(x -> x == -x_length/2, X_nodes)[1], findall(x -> x == x_length/2, X_nodes)[1])
    die_nodes_y = (findall(x -> x == -y_length/2, Y_nodes)[1], findall(x -> x == y_length/2, Y_nodes)[1])

    cell_centers_x = [(X_nodes[i]+X_nodes[i+1])/2 for i in die_nodes_x[1]:die_nodes_x[2]-1]
    cell_centers_y = [(Y_nodes[i]+Y_nodes[i+1])/2 for i in die_nodes_y[1]:die_nodes_y[2]-1]

    heat_values = read_csv(power_file)
    source_Nx, source_Ny = size(heat_values)
    source_x = range(-x_length/2, x_length/2, source_Nx) |> collect
    source_y = range(-y_length/2, y_length/2, source_Ny) |> collect

    interpolated_heat = transpose(interpolate_(source_x,source_y,heat_values,cell_centers_x,cell_centers_y))
    diff_ = abs.(Z_nodes.-source_z_val)
    source_index = findall(x -> x < 1e-6, diff_)[1]-1
    source[die_nodes_x[1]:die_nodes_x[2]-1,die_nodes_y[1]:die_nodes_y[2]-1,source_index] = interpolated_heat
    return source
end

function do_plotting(sol, sol_wd, settings, scenario_name, interpolation)
    result_ = convert_units("temperature", sol[end, 2:end-1,2:end-1,1], "K", "C")'
    # result_ = sol[end][2:end-1,2:end-1,2]
    (delta_x, delta_y, delta_z), (Nx,Ny,Nz), (x_mesh, y_mesh, z_mesh) = create_mesh(settings)
    z_length = settings["model"]["bodies"]["die"]["size"]["Z"]
    
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

function apply_bc(settings, u, k, delta_x_s, delta_x_e, delta_y_s, delta_y_e, delta_z_s, delta_z_e)
    
    # At X-
    side = "X-"
    if settings["BC"]["X-"]["type"] == "constant_T"
        const_temp = convert_units("temperature", settings["BC"]["X-"]["value"]["value"], settings["units"]["temperature"], "K") 
        u[1,:,:] .= 2*const_temp .- u[2,:,:]
    elseif settings["BC"]["X-"]["type"] == "insulated"
        u[1,:,:] .= u[2,:,:]
    elseif settings["BC"]["X-"]["type"] == "const_flux"
        @. u[1,:,:] = u[2,:,:] + settings["BC"]["X-"]["value"]["value"] * (delta_x_s / k[1, :, :])
    elseif settings["BC"]["X-"]["type"] == "HTC"
        h = settings["BC"]["X-"]["value"]["value"]
        t_amb = convert_units("temperature", settings["BC"]["X-"]["value"]["t_amb"], settings["units"]["temperature"], "K")
        @. u[1,:,:] .= u[2,:,:] * (1 - delta_x_s * h / k[1,:,:]) + (h*delta_x_s / k[1, :, :]) * t_amb
    end

    # At X+
    side = "X+"
    if settings["BC"]["X+"]["type"] == "constant_T"
        const_temp = convert_units("temperature", settings["BC"]["X+"]["value"]["value"], settings["units"]["temperature"], "K") 
        u[end,:,:] .= 2*settings["BC"]["X+"]["value"] .- u[end-1,:,:]
    elseif settings["BC"]["X+"]["type"] == "insulated"
        u[end,:,:] .= u[end-1,:,:]
    elseif settings["BC"]["X+"]["type"] == "const_flux"
        @. u[end,:,:] = u[end-1,:,:] + settings["BC"]["X+"]["value"]["value"] * (delta_x_e / k[end, :, :])
    elseif settings["BC"]["X+"]["type"] == "HTC"
        h = settings["BC"]["X+"]["value"]["value"]
        t_amb = convert_units("temperature", settings["BC"]["X+"]["value"]["t_amb"], settings["units"]["temperature"], "K")
        @. u[end,:,:] .= u[end-1,:,:]*(1 - delta_x_e * h ./ k[end, : ,:]) .+ (h*delta_x_e ./ k[end, :, :]) * t_amb
    end

    # Y-
    side = "Y-"
    if settings["BC"]["Y-"]["type"] == "constant_T"
        const_temp = convert_units("temperature", settings["BC"]["Y-"]["value"]["value"], settings["units"]["temperature"], "K") 
        u[:,1,:] .= 2*const_temp .- u[:,2,:]
    elseif settings["BC"]["Y-"]["type"] == "insulated"
        u[:,1,:] .= u[:,2,:]
    elseif settings["BC"]["Y-"]["type"] == "const_flux"
        @. u[:,1,:] = u[:,2,:] + settings["BC"]["Y-"]["value"]["value"] * (delta_y_s / k[:, 1, :])
    elseif settings["BC"]["Y-"]["type"] == "HTC"
        h = settings["BC"]["Y-"]["value"]["value"]
        t_amb = convert_units("temperature", settings["BC"]["Y-"]["value"]["t_amb"], settings["units"]["temperature"], "K")
        @. u[:,1,:] .= u[:,2,:]*(1 - delta_y_s * h ./ k[:, 1, :]) .+ (h*delta_y_s ./ k[:, 1, :]) * t_amb
    end
    
    # Y+
    side = "Y+"
    if settings["BC"]["Y+"]["type"] == "constant_T"
        const_temp = convert_units("temperature", settings["BC"]["Y+"]["value"]["value"], settings["units"]["temperature"], "K") 
        u[:,end,:] .= 2*const_temp .- u[:,end-1,:]
    elseif settings["BC"]["Y+"]["type"] == "insulated"
        u[:,end,:] .= u[:,end-1,:]
    elseif settings["BC"]["Y+"]["type"] == "const_flux"
        @. u[:,end,:] = u[:,end-1,:] + settings["BC"]["Y+"]["value"]["value"] * (delta_y_e / k[:, end, :])
    elseif settings["BC"]["Y+"]["type"] == "HTC"
        h = settings["BC"]["Y+"]["value"]["value"]
        t_amb = convert_units("temperature", settings["BC"]["Y+"]["value"]["t_amb"], settings["units"]["temperature"], "K")
        @. u[:,end,:] .= u[:,end-1,:]*(1 - delta_y_e * h ./ k[:, end, :]) .+ (h*delta_y_e ./ k[:, end, :]) * t_amb
    end

    # Z-
    side = "Z-"
    if settings["BC"]["Z-"]["type"] == "constant_T"
        const_temp = convert_units("temperature", settings["BC"]["Z-"]["value"]["value"], settings["units"]["temperature"], "K") 
        u[:,:,1] .= 2*const_temp .- u[:,:,2]
    elseif settings["BC"]["Z-"]["type"] == "insulated"
        u[:,:,1] .= u[:,:,2]
    elseif settings["BC"]["Z-"]["type"] == "const_flux"
        @. u[:,:,1] = u[:,:,2] + settings["BC"]["Z-"]["value"]["value"] * (delta_z_s / k[:, :, 1])
    elseif settings["BC"]["Z-"]["type"] == "HTC"
        h = settings["BC"]["Z-"]["value"]["value"]
        t_amb = convert_units("temperature", settings["BC"]["Z-"]["value"]["t_amb"], settings["units"]["temperature"], "K")
        @. u[:,:,1] .= u[:,:,2]*(1 - delta_z_s * h ./ k[:, :, 1]) .+ (h*delta_z_s ./ k[:, :, 1]) * t_amb
    end

    # Z+
    side = "Z+"
    if settings["BC"]["Z+"]["type"] == "constant_T"
        const_temp = convert_units("temperature", settings["BC"]["Z+"]["value"]["value"], settings["units"]["temperature"], "K") 
        u[:,:,end] .= 2*const_temp .- u[:,:,end-1]
    elseif settings["BC"]["Z+"]["type"] == "insulated"
        u[:,:,end] .= u[:,:,end-1]
    elseif settings["BC"]["Z+"]["type"] == "const_flux"
        @. u[:,:,end] = u[:,:,end-1] + settings["BC"]["Z+"]["value"]["value"] * (delta_z_e / k[:, :, end])
    elseif settings["BC"]["Z+"]["type"] == "HTC"
        h = settings["BC"]["Z+"]["value"]["value"]
        t_amb = convert_units("temperature", settings["BC"]["Z+"]["value"]["t_amb"], settings["units"]["temperature"], "K")
        @. u[:,:,end] = u[:,:,end-1]*(1 - delta_z_e * h / k[:, :, end]) + (h*delta_z_e / k[:, :, end]) * t_amb
    end
end

function get_alpha((Nx,Ny,Nz), nodes, k_by_rho_cp)
    X_nodes, Y_nodes, Z_nodes = nodes

    alpha_x_i_minus_half, alpha_x_i_plus_half, alpha_y_i_minus_half, alpha_y_i_plus_half, alpha_z_i_minus_half, alpha_z_i_plus_half = ones(Nx,Ny,Nz), 
    ones(Nx,Ny,Nz), ones(Nx,Ny,Nz), ones(Nx,Ny,Nz), ones(Nx,Ny,Nz), ones(Nx,Ny,Nz)

    for k_ in range(2, Nz-1)
        for j_ in range(2, Ny-1)
            for i_ in range(2, Nx-1)
                dx_i_minus_half = (X_nodes[i_+1] - X_nodes[i_-1])/2
                dy_i_minus_half = (Y_nodes[j_+1] - Y_nodes[j_-1])/2
                dz_i_minus_half = (Z_nodes[k_+1] - Z_nodes[k_-1])/2

                dx_i_plus_half = (X_nodes[i_+2] - X_nodes[i_])/2
                dy_i_plus_half = (Y_nodes[j_+2] - Y_nodes[j_])/2
                dz_i_plus_half = (Z_nodes[k_+2] - Z_nodes[k_])/2

                dx = X_nodes[i_+1] - X_nodes[i_]
                dy = Y_nodes[j_+1] - Y_nodes[j_]
                dz = Z_nodes[k_+1] - Z_nodes[k_]

                alpha_x_i_minus_half[i_, j_, k_] = k_by_rho_cp[i_, j_, k_]/(dx * dx_i_minus_half)
                alpha_x_i_plus_half[i_, j_, k_] = k_by_rho_cp[i_, j_, k_]/(dx * dx_i_plus_half)

                alpha_y_i_minus_half[i_, j_, k_] = k_by_rho_cp[i_, j_, k_]/(dy * dy_i_minus_half)
                alpha_y_i_plus_half[i_, j_, k_] = k_by_rho_cp[i_, j_, k_]/(dy * dy_i_plus_half)

                alpha_z_i_minus_half[i_, j_, k_] = k_by_rho_cp[i_, j_, k_]/(dz * dz_i_minus_half)
                alpha_z_i_plus_half[i_, j_, k_] = k_by_rho_cp[i_, j_, k_]/(dz * dz_i_plus_half)
            end
        end
    end

    return alpha_x_i_minus_half, alpha_x_i_plus_half, alpha_y_i_minus_half, alpha_y_i_plus_half, alpha_z_i_minus_half, alpha_z_i_plus_half
end

function conduction_3d_loop!(du, u, p, t)
    (X_nodes, Y_nodes, Z_nodes), power_file, alpha_x_i_minus_half, alpha_x_i_plus_half, alpha_y_i_minus_half, alpha_y_i_plus_half, alpha_z_i_minus_half, alpha_z_i_plus_half, k, settings = p

    Nx, Ny, Nz = size(u)

    source = define_volume_sources(power_file, settings, Nx,Ny,Nz, X_nodes, Y_nodes, Z_nodes)

    apply_bc(settings, u, k, abs(X_nodes[2]-X_nodes[1]), abs(X_nodes[end] - X_nodes[end-1]), 
                abs(Y_nodes[2] - Y_nodes[1]), abs(Y_nodes[end] - Y_nodes[end-1]), 
                abs(Z_nodes[2]-Z_nodes[1]), abs(Z_nodes[end] - Z_nodes[end-1]))

    Threads.@threads for k_ in range(2, Nz-1)
        for j_ in range(2, Ny-1)
            for i_ in range(2, Nx-1)
                ip1, im1, jp1, jm1, kp1, km1 = i_ + 1, i_ - 1, j_ + 1, j_ - 1, k_+1, k_-1
                du[i_, j_, k_] =  alpha_x_i_plus_half[i_, j_, k_] * (u[ip1, j_, k_] - u[i_, j_, k_]) - alpha_x_i_minus_half[i_, j_, k_] * (u[i_, j_, k_] - u[im1, j_, k_]) + 
                                  alpha_y_i_plus_half[i_, j_, k_] * (u[i_, jp1, k_] - u[i_, j_, k_]) - alpha_y_i_minus_half[i_, j_, k_] * (u[i_, j_, k_] - u[i_, jm1, k_]) + 
                                  alpha_z_i_plus_half[i_, j_, k_] * (u[i_, j_, kp1] - u[i_, j_, k_]) - alpha_z_i_minus_half[i_, j_, k_] * (u[i_, j_, k_] - u[i_, j_, km1]) + 
                                  source[i_, j_, k_]
            end
        end
    end
    
    nothing
end

function read_csv(file_name)
    data = CSV.read(file_name, DataFrame, header=false)
    mat = Matrix(data)
    return mat
end

function interpolate_(x,y,z,new_x,new_y)
    itp = LinearInterpolation((x, y), z)
    new_z = [itp(x,y) for y in new_y, x in new_x]
    return new_z
end

function save_fields(sol, t_, sol_wd, scenario_name)
    sol = sol .- 273.15
    h5open(joinpath(sol_wd,scenario_name*"__solution.sol"), "w") do file
        g = HDF5.create_group(file, "solution")

        for (i,t) in enumerate(t_)
            g[string(t)] = sol[i,:,:,:]
        end
    end
end

function solve_(working_dir, sol_wd, power_file, settings, progress_file_name)

    progress_file = open(progress_file_name, "w")
    close(progress_file)

    (X_nodes, Y_nodes, Z_nodes) , (Nx,Ny,Nz) = generate_mesh(settings)

    min_dx = minimum([X_nodes[i]-X_nodes[i-1] for i in range(2,size(X_nodes)[1])])
    min_dy = minimum([Y_nodes[i]-Y_nodes[i-1] for i in range(2,size(Y_nodes)[1])])
    min_dz = minimum([Z_nodes[i]-Z_nodes[i-1] for i in range(2,size(Z_nodes)[1])])

    println("Mesh generation done")
    println("Mesh details - Total cells :", (Nx-2)*(Ny-2)*(Nz-2)/1e6, " M ", "with ", Nx-2, ",", Ny-2, ", and ", Nz-2, " in X,Y and Z")
    println("Smallest cells in X, Y and Z are ", min_dx, " , ", min_dy, " , ", min_dz)

    k_by_rho_cp, k = get_k_by_rho_cp((X_nodes, Y_nodes, Z_nodes), (Nx, Ny, Nz), settings)

    write_mesh_and_conductivity_to_file(k, X_nodes, Y_nodes, Z_nodes, sol_wd)

    u0 = zeros(Nx,Ny,Nz)

    alpha_x_i_minus_half, alpha_x_i_plus_half, alpha_y_i_minus_half, alpha_y_i_plus_half, alpha_z_i_minus_half, alpha_z_i_plus_half = get_alpha((Nx,Ny,Nz), (X_nodes, Y_nodes, Z_nodes), k_by_rho_cp)

    p = ((X_nodes, Y_nodes, Z_nodes), joinpath(working_dir,power_file), alpha_x_i_minus_half, alpha_x_i_plus_half, alpha_y_i_minus_half, alpha_y_i_plus_half, alpha_z_i_minus_half, alpha_z_i_plus_half, k, settings)

    u0 = initialize_domain!(u0, settings)

    start_time = settings["start_time"]
    end_time = settings["end_time"]

    problem = ODEProblem(conduction_3d_loop!, u0, (start_time, end_time), p)

    dt = settings["dt"]
    n_steps = round(Int, settings["end_time"]/dt)

    sol = zeros((n_steps+1, size(u0)...))
    t_ = zeros(n_steps+1)

    sol[1, :,:,:] = u0
    t_ = 0

    integrator = init(problem, CVODE_BDF(linear_solver=:GMRES) ; reltol=1e-6, abstol=1e-3, maxiters=1000, progress=true, save_everystep=false)
    # integrator = init(problem, reltol=1e-, abstol=1e-3, maxiters=10000, progress=true)
    for i in range(1,n_steps)
        step!(integrator, dt, true)
        # step!(integrator)
        t,u = integrator.t, integrator.u
        println("Solution time : ",t, " , at ", Dates.format(now(), "HH:MM:SS"))
        sol[i+1,:,:,:] = u
        t_[i+1] = t

        if t > end_time
            terminate!(integrator)
            break
        end

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
    settings = JSON.parse(f)
    # s1 = JSON.parse(f)

    # settings = Dict{Symbol, Any}()
    # for (k,v) in s1
    #     settings[Symbol(k)] = v
    # end

    date_start = now()
    println("Start solution at : ", Dates.format(date_start, "HH:MM:SS"))
    sol,t_ = solve_(working_dir, sol_wd, scenario_name, settings, progress_file_name);

    #do_plotting(sol, sol_wd, settings, scenario_name, false);
    
    save_fields(sol, t_, sol_wd, scenario_name);
    
    date_end = now()
    println("End solution at : ", Dates.format(date_end, "HH:MM:SS"))
    println("Time taken for solution : ", (date_end-date_start)/Millisecond(1000), " [s]")

    flush(stdout)
    exit()
end

end
