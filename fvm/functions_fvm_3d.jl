int(x) = floor(Int, x)

function build_domain(settings)

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

    u = zeros(Nx, Ny, Nz)

    return (x_mesh, y_mesh, z_mesh), u, (Nx, Ny, Nz), z_length
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

function apply_bc(u, settings, delta_x, delta_y, delta_z)
    k = settings["model"]["bodies"]["die"]["material"]["k"]
    rho = settings["model"]["bodies"]["die"]["material"]["rho"]
    cp = settings["model"]["bodies"]["die"]["material"]["cp"]

    ## TODO convert the bc application to 3d
    # At X-
    side = "X-"
    if settings["BC"]["X-"]["type"] == "constant_T"
        const_temp = convert_units("temperature", settings["BC"]["X-"]["value"]["value"], settings["units"]["temperature"], "K") 
        u[1,:,:] .= 2*const_temp .- u[2,:,:]
    elseif settings["BC"]["X-"]["type"] == "insulated"
        u[1,:,:] .= u[2,:,:]
    elseif settings["BC"]["X-"]["type"] == "const_flux"
        u[1,:,:] .= u[2,:,:] .+ settings["BC"]["X-"]["value"]["value"] * (delta_x/k)
    elseif settings["BC"]["X-"]["type"] == "HTC"
        h = settings["BC"]["X-"]["value"]["value"]
        t_amb = convert_units("temperature", settings["BC"]["X-"]["value"]["t_amb"], settings["units"]["temperature"], "K")
        u[1,:,:] .= u[2,:,:]*(1 - delta_x * h/k) .+ (h*delta_x/k) * t_amb
    end

    # At X+
    side = "X+"
    if settings["BC"]["X+"]["type"] == "constant_T"
        const_temp = convert_units("temperature", settings["BC"]["X+"]["value"]["value"], settings["units"]["temperature"], "K") 
        u[end,:,:] .= 2*settings["BC"]["X+"]["value"] .- u[end-1,:,:]
    elseif settings["BC"]["X+"]["type"] == "insulated"
        u[end,:,:] .= u[end-1,:,:]
    elseif settings["BC"]["X+"]["type"] == "const_flux"
        u[end,:,:] .= u[end-1,:,:] .+ settings["BC"]["X+"]["value"]["value"] * (delta_x/k)
    elseif settings["BC"]["X+"]["type"] == "HTC"
        h = settings["BC"]["X+"]["value"]["value"]
        t_amb = convert_units("temperature", settings["BC"]["X+"]["value"]["t_amb"], settings["units"]["temperature"], "K")
        u[end,:,:] .= u[end-1,:,:]*(1 - delta_x * h/k) .+ (h*delta_x/k) * t_amb
    end

    # Y-
    side = "Y-"
    if settings["BC"]["Y-"]["type"] == "constant_T"
        const_temp = convert_units("temperature", settings["BC"]["Y-"]["value"]["value"], settings["units"]["temperature"], "K") 
        u[:,1,:] .= 2*const_temp .- u[:,2,:]
    elseif settings["BC"]["Y-"]["type"] == "insulated"
        u[:,1,:] .= u[:,2,:]
    elseif settings["BC"]["Y-"]["type"] == "const_flux"
        u[:,1,:] .= u[:,2,:] .+ settings["BC"]["Y-"]["value"]["value"] * (delta_y/k)
    elseif settings["BC"]["Y-"]["type"] == "HTC"
        h = settings["BC"]["Y-"]["value"]["value"]
        t_amb = convert_units("temperature", settings["BC"]["Y-"]["value"]["t_amb"], settings["units"]["temperature"], "K")
        u[:,1,:] .= u[:,2,:]*(1 - delta_y * h/k) .+ (h*delta_y/k) * t_amb
    end
    
    # Y+
    side = "Y+"
    if settings["BC"]["Y+"]["type"] == "constant_T"
        const_temp = convert_units("temperature", settings["BC"]["Y+"]["value"]["value"], settings["units"]["temperature"], "K") 
        u[:,end,:] .= 2*const_temp .- u[:,end-1,:]
    elseif settings["BC"]["Y+"]["type"] == "insulated"
        u[:,end,:] .= u[:,end-1,:]
    elseif settings["BC"]["Y+"]["type"] == "const_flux"
        u[:,end,:] .= u[:,end-1,:] .+ settings["BC"]["Y+"]["value"]["value"] * (delta_y/k)
    elseif settings["BC"]["Y+"]["type"] == "HTC"
        h = settings["BC"]["Y+"]["value"]["value"]
        t_amb = convert_units("temperature", settings["BC"]["Y+"]["value"]["t_amb"], settings["units"]["temperature"], "K")
        u[:,end,:] .= u[:,end-1,:]*(1 - delta_y * h/k) .+ (h*delta_y/k) * t_amb
    end

    # Z-
    side = "Z-"
    if settings["BC"]["Z-"]["type"] == "constant_T"
        const_temp = convert_units("temperature", settings["BC"]["Z-"]["value"]["value"], settings["units"]["temperature"], "K") 
        u[:,:,1] .= 2*const_temp .- u[:,:,2]
    elseif settings["BC"]["Z-"]["type"] == "insulated"
        u[:,:,1] .= u[:,:,2]
    elseif settings["BC"]["Z-"]["type"] == "const_flux"
        u[:,:,1] .= u[:,:,2] .+ settings["BC"]["Z-"]["value"]["value"] * (delta_z/k)
    elseif settings["BC"]["Z-"]["type"] == "HTC"
        h = settings["BC"]["Z-"]["value"]["value"]
        t_amb = convert_units("temperature", settings["BC"]["Z-"]["value"]["t_amb"], settings["units"]["temperature"], "K")
        u[:,:,1] .= u[:,:,2]*(1 - delta_z * h/k) .+ (h*delta_z/k) * t_amb
    end

    # Z+
    side = "Z+"
    if settings["BC"]["Z+"]["type"] == "constant_T"
        const_temp = convert_units("temperature", settings["BC"]["Z+"]["value"]["value"], settings["units"]["temperature"], "K") 
        u[:,:,end] .= 2*const_temp .- u[:,:,end-1]
    elseif settings["BC"]["Z+"]["type"] == "insulated"
        u[:,:,end] .= u[:,:,end-1]
    elseif settings["BC"]["Z+"]["type"] == "const_flux"
        u[:,:,end] .= u[:,:,end-1] .+ settings["BC"]["Z+"]["value"]["value"] * (delta_z/k)
    elseif settings["BC"]["Z+"]["type"] == "HTC"
        h = settings["BC"]["Z+"]["value"]["value"]
        t_amb = convert_units("temperature", settings["BC"]["Z+"]["value"]["t_amb"], settings["units"]["temperature"], "K")
        u[:,:,end] .= u[:,:,end-1]*(1 - delta_z * h/k) .+ (h*delta_z/k) * t_amb
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