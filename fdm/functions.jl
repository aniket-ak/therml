function build_domain(settings, num_sources)

    N = 2*num_sources+1 # excluding ghost cells

    if settings["BC"]["X-"]["type"] == "constant_T"
        x_start = 2
        x_start_ng = 0
    else
        x_start = 2
        x_start_ng = 1
    end
    
    if settings["BC"]["Y-"]["type"] == "constant_T"
        y_start = 2
        y_start_ng = 0
    else
        y_start = 1
        y_start_ng = 1
    end

    if settings["BC"]["Z-"]["type"] == "constant_T"
        z_start = 2
        z_start_ng = 0
    else
        z_start = 1
        z_start_ng = 1
    end

    if settings["BC"]["X+"]["type"] == "constant_T"
        x_end = N-1
        x_end_ng = 0
    else
        x_end = N
        x_end_ng = 1
    end

    if settings["BC"]["Y+"]["type"] == "constant_T"
        y_end = N-1
        y_end_ng = 0
    else
        y_end = N
        y_end_ng = 1
    end

    if settings["BC"]["Z+"]["type"] == "constant_T"
        z_end = N-1
        z_end_ng = 0
    else
        z_end = N
        z_end_ng = 1
    end

    xy_mesh = range(0, stop = 1, length = N)

    u = zeros(N+x_start_ng+x_end_ng, N+y_start_ng+y_end_ng)
    ghost_cells = (x_start_ng, x_end_ng, y_start_ng, y_end_ng, z_start_ng, z_end_ng)

    return (x_start,x_end), (y_start, y_end), (z_start, z_end), xy_mesh, u, N, ghost_cells
end