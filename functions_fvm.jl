function build_domain(settings, num_sources)

    N = 2*num_sources + 2 # including ghost cells

    xy_mesh = range(0, stop = 1, length = N-1)

    u = zeros(N, N)

    return xy_mesh, u, N
end