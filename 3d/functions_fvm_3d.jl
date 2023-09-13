int(x) = floor(Int, x)

function create_mesh(settings)
    x_length = settings["model"]["bodies"]["die"]["size"]["X"]
    y_length = settings["model"]["bodies"]["die"]["size"]["Y"]
    z_length = settings["model"]["bodies"]["die"]["size"]["Z"]

    x_normalized = x_length/z_length
    y_normalized = y_length/z_length
    z_normalized = z_length/z_length

    Nx = 2*settings["model"]["num_sources"]["X"] + 2 # including ghost cells
    Ny = 2*settings["model"]["num_sources"]["Y"] + 2 # including ghost cells
    Nz = max(2, 2*int(settings["model"]["bodies"]["die"]["size"]["Z"]/settings["model"]["smallest_thickness"])) + 2
    
    x_mesh = range(0, stop = x_normalized, length = Nx+1)
    y_mesh = range(0, stop = y_normalized, length = Ny+1)
    z_mesh = range(0, stop = z_normalized, length = Nz+1)

    delta_x, delta_y, delta_z = step(x_mesh), step(y_mesh), step(z_mesh)
    return (delta_x, delta_y, delta_z), (Nx,Ny,Nz), (x_mesh, y_mesh, z_mesh)
end

function build_domain(Nx,Ny,Nz)
    u = zeros(Nx, Ny, Nz)
    println("total mesh count: ", (Nx-2)*(Ny-2)*(Nz-2)/1e3, " k")
    println("Mesh elements across X, Y and Z ", Nx-2," , ",Ny-2," , ",Nz-2)
    return u
end

function initialize_domain!(u0)
    T_initial = convert_units("temperature", settings["IC"], temperature_base_units, "K")
    u0 = u0 .+ T_initial
    return u0
end

function define_volume_sources(Nx,Ny,Nz)
    source = zeros(Nx,Ny,Nz)
    # source[10:20, 50:70, 1:5] .= 1
    heat_values = read_csv("./heat_sources.csv")
    source[:,:,trunc(Int, Nz/2)] = heat_values
    return source
end

function do_plotting(sol)
    result_ = convert_units("temperature", sol[end][2:end-1,2:end-1,2], "K", "C")'
    # result_ = sol[end][2:end-1,2:end-1,2]
    (delta_x, delta_y, delta_z), (Nx,Ny,Nz), (x_mesh, y_mesh, z_mesh) = create_mesh(settings)
    z_length = settings["model"]["bodies"]["die"]["size"]["Z"]
    
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
        # Layout(autosize=true)
    )

    open("./plot.html", "w") do io
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

function configure_problem_klu!(u0,p)
    du0 = zeros(size(u0))
    jac_sparsity = Symbolics.jacobian_sparsity((du, u) -> conduction_3d_loop!(du, u, p, 0.0), du0, u0)

    f = ODEFunction(conduction_3d_loop!; jac_prototype = float.(jac_sparsity))

    start_time = settings["start_time"]
    end_time = settings["end_time"]

    prob_conduction_3d_sparse = ODEProblem(f, u0, (start_time, end_time), p)

    algorithm = KenCarp47(linsolve = KLUFactorization())

    return prob_conduction_3d_sparse, algorithm
end

function configure_problem_ode!(u0,p)
    start_time = settings["start_time"]
    end_time = settings["end_time"]

    problem = ODEProblem(conduction_3d_loop!, u0, (start_time, end_time), p)
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