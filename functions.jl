function get_indices(settings, N)
    if settings["BC"]["X-"]["type"] == "constant_T"
        x_start = 2
    elseif settings["BC"]["X-"]["type"] == "adiabatic"
        x_start = 1
    elseif settings["BC"]["X-"]["type"] == "HTC"
        x_start = 1
    end
    
    if settings["BC"]["Y-"]["type"] == "constant_T"
        y_start = 2
    elseif settings["BC"]["Y-"]["type"] == "adiabatic"
        y_start = 1
    elseif settings["BC"]["Y-"]["type"] == "HTC"
        y_start = 1
    end

    if settings["BC"]["Z-"]["type"] == "constant_T"
        z_start = 2
    elseif settings["BC"]["Z-"]["type"] == "adiabatic"
        z_start = 1
    elseif settings["BC"]["Z-"]["type"] == "HTC"
        z_start = 1
    end

    if settings["BC"]["X+"]["type"] == "constant_T"
        x_end = N-1
    elseif settings["BC"]["X+"]["type"] == "adiabatic"
        x_end = N
    elseif settings["BC"]["X+"]["type"] == "HTC"
        x_end = N
    end

    if settings["BC"]["Y+"]["type"] == "constant_T"
        y_end = N-1
    elseif settings["BC"]["Y+"]["type"] == "adiabatic"
        y_end = N
    elseif settings["BC"]["Y+"]["type"] == "HTC"
        y_end = N
    end

    if settings["BC"]["Z+"]["type"] == "constant_T"
        z_end = N-1
    elseif settings["BC"]["Z+"]["type"] == "adiabatic"
        z_end = N
    elseif settings["BC"]["Z+"]["type"] == "HTC"
        z_end = N
    end

    return x_start, y_start, z_start, x_end, y_end, z_end
end