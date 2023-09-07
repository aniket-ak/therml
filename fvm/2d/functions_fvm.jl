function build_domain(settings, num_sources)

    N = 2*50 + 2 # including ghost cells

    xy_mesh = range(0, stop = 1, length = N-1)

    u = zeros(N, N)

    return xy_mesh, u, N
end

function convert_units(quantity, value, from_units, to_units)
    if quantity == "temperature"
        if from_units == "C" && to_units == "K"
            return_value = value .+ 273.15
            return return_value
        elseif from_units == "K" && to_units == "C"
            return_value = value .- 273.15
            return return_value
        end
    else
        return value
    end
end

function apply_bc(u, settings, mesh_size)
    k = settings["model"]["bodies"]["die"]["material"]["k"]
    rho = settings["model"]["bodies"]["die"]["material"]["rho"]
    cp = settings["model"]["bodies"]["die"]["material"]["cp"]
    
    # bc_config = Dict(
    #     "X-" => Dict(
    #         "ghost_cells" => u[1,:],
    #         "b_cells" => u[2,:]
    #     ),
    #     "X+" => Dict(
    #         "ghost_cells" => u[end,:],
    #         "b_cells" => u[end-1,:]
    #     ),
    #     "Y-" => Dict(
    #         "ghost_cells" => u[:,1],
    #         "b_cells" => u[:,2]
    #     ),
    #     "Y+" => Dict(
    #         "ghost_cells" => u[:,end],
    #         "b_cells" => u[:,end-1]
    #     ),
    #     # "Z-" => Dict(
    #     #     "ghost_cells" => u[1,:],
    #     #     "b_cells" => u[2,:]
    #     # ),
    #     # "Z+" => Dict(
    #     #     "ghost_cells" => u[1,:],
    #     #     "b_cells" => u[2,:]
    #     # ),
    # )

    # for side in ["X-", "X+", "Y-", "Y+"]
        
    #     ghost_cells = bc_config[side]["ghost_cells"]
    #     b_cells = bc_config[side]["b_cells"]

    #     if settings["BC"][side]["type"] == "constant_T"
    #         println(side, "T")
    #         const_temp = convert_units("temperature", settings["BC"][side]["value"]["value"], settings["units"]["temperature"], "K") 
    #         ghost_cells .= 2*const_temp .- b_cells
    #     elseif settings["BC"][side]["type"] == "insulated"
    #         println(side, "q=0")
    #         ghost_cells .= b_cells
    #     elseif settings["BC"][side]["type"] == "const_flux"
    #         println(side, "q=c")
    #         ghost_cells .= b_cells .+ settings["BC"][side]["value"]["value"] * (mesh_size/k)
    #     elseif settings["BC"][side]["type"] == "HTC"
    #         println(side, "HTC")
    #         h = settings["BC"][side]["value"]["value"]
    #         t_amb = convert_units("temperature", settings["BC"][side]["value"]["t_amb"], settings["units"]["temperature"], "K")
    #         ghost_cells .= b_cells*(mesh_size * h/k - 1) .- (h*mesh_size/k) * t_amb
    #     end
    # end

    # At X-
    side = "X-"
    if settings["BC"]["X-"]["type"] == "constant_T"
        const_temp = convert_units("temperature", settings["BC"]["X-"]["value"]["value"], settings["units"]["temperature"], "K") 
        u[1,:] .= 2*const_temp .- u[2,:]
    elseif settings["BC"]["X-"]["type"] == "insulated"
        u[1,:] .= u[2,:]
    elseif settings["BC"]["X-"]["type"] == "const_flux"
        u[1,:] .= u[2,:] .+ settings["BC"]["X-"]["value"]["value"] * (mesh_size/k)
    elseif settings["BC"]["X-"]["type"] == "HTC"
        h = settings["BC"]["X-"]["value"]["value"]
        t_amb = convert_units("temperature", settings["BC"]["X-"]["value"]["t_amb"], settings["units"]["temperature"], "K")
        u[1,:] .= u[2,:]*(1 - mesh_size * h/k) .+ (h*mesh_size/k) * t_amb
    end

    # At X+
    side = "X+"
    if settings["BC"]["X+"]["type"] == "constant_T"
        const_temp = convert_units("temperature", settings["BC"]["X+"]["value"]["value"], settings["units"]["temperature"], "K") 
        u[end,:] .= 2*settings["BC"]["X+"]["value"] .- u[end-1,:]
    elseif settings["BC"]["X+"]["type"] == "insulated"
        u[end,:] .= u[end-1,:]
    elseif settings["BC"]["X+"]["type"] == "const_flux"
        u[end,:] .= u[end-1,:] .+ settings["BC"]["X+"]["value"]["value"] * (mesh_size/k)
    elseif settings["BC"]["X+"]["type"] == "HTC"
        h = settings["BC"]["X+"]["value"]["value"]
        t_amb = convert_units("temperature", settings["BC"]["X+"]["value"]["t_amb"], settings["units"]["temperature"], "K")
        u[end,:] .= u[end-1,:]*(1 - mesh_size * h/k) .+ (h*mesh_size/k) * t_amb
    end

    # Y-
    side = "Y-"
    if settings["BC"]["Y-"]["type"] == "constant_T"
        const_temp = convert_units("temperature", settings["BC"]["Y-"]["value"]["value"], settings["units"]["temperature"], "K") 
        u[:,1] .= 2*const_temp .- u[:,2]
    elseif settings["BC"]["Y-"]["type"] == "insulated"
        u[:,1] .= u[:,2]
    elseif settings["BC"]["Y-"]["type"] == "const_flux"
        u[:,1] .= u[:,2] .+ settings["BC"]["Y-"]["value"]["value"] * (mesh_size/k)
    elseif settings["BC"]["Y-"]["type"] == "HTC"
        h = settings["BC"]["Y-"]["value"]["value"]
        t_amb = convert_units("temperature", settings["BC"]["Y-"]["value"]["t_amb"], settings["units"]["temperature"], "K")
        u[:,1] .= u[:,1]*(1 - mesh_size * h/k) .+ (h*mesh_size/k) * t_amb
    end
    
    # Y+
    side = "Y+"
    if settings["BC"]["Y+"]["type"] == "constant_T"
        const_temp = convert_units("temperature", settings["BC"]["Y+"]["value"]["value"], settings["units"]["temperature"], "K") 
        u[:,end] .= 2*const_temp .- u[:,end-1]
    elseif settings["BC"]["Y+"]["type"] == "insulated"
        u[:,end] .= u[:,end-1]
    elseif settings["BC"]["Y+"]["type"] == "const_flux"
        u[:,end] .= u[:,end-1] .+ settings["BC"]["Y+"]["value"]["value"] * (mesh_size/k)
    elseif settings["BC"]["Y+"]["type"] == "HTC"
        h = settings["BC"]["Y+"]["value"]["value"]
        t_amb = convert_units("temperature", settings["BC"]["Y+"]["value"]["t_amb"], settings["units"]["temperature"], "K")
        u[:,end] .= u[:,end-1]*(1 - mesh_size * h/k) .+ (h*mesh_size/k) * t_amb
    end
end

function read_csv(file_name)

    data = CSV.read(file_name, DataFrame,)

    mat = Matrix(data)
    return mat
end

function read_excel(file_name)
    xf = XLSX.readxlsx(file_name)
    XLSX.sheetnames(xf)
    sh = xf["Sheet1"]
    return sh[:]

end

function interpolate()
    x = range(-2, 3, length=20)
    y = range(3, 4, length=10)
    z = @. cos(x) + sin(y')
    # Interpolatant object
    itp = LinearInterpolation((x, y), z)
    # Fine grid
    x2 = range(extrema(x)..., length=300)
    y2 = range(extrema(y)..., length=200)
    # Interpolate
    z2 = [itp(x,y) for y in y2, x in x2]
    # Plot
    p = heatmap(x2, y2, z2, clim=(-2,2), title="Interpolated heatmap")
    scatter!(p, [x for _ in y for x in x], [y for y in y for _ in x], zcolor=z[:]; lab="original data", clim=(-2,2))

    # plot(z)
    plot(z2)

    z2
end
