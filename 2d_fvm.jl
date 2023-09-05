using DifferentialEquations, LinearAlgebra, SparseArrays
using PlotlyJS
using Symbolics
using JSON
using Dates
include("./functions_fvm.jl")

f = open("settings.json", "r")
settings = JSON.parse(f)

num_sources = settings["model"]["num_sources"]
xy_mesh, u0, N = build_domain(settings, num_sources)

source = zeros(N,N)
source[10:20, 10:20] .= 10.0

println(Dates.format(now(), "HH:MM:SS"))

function conduction_2d_loop!(du, u, p, t)
    source, alpha, dx = p
    alpha = alpha / dx^2
    
    if settings["BC"]["X-"]["type"] == "constant_T"
        u[1,:] .= 2*settings["BC"]["X-"]["value"] .- u[2,:]
    elseif settings["BC"]["X-"]["type"] == "insulated"
        u[1,:] .= u[2,:]
    end
    if settings["BC"]["X+"]["type"] == "constant_T"
        u[end,:] .= 2*settings["BC"]["X+"]["value"] .- u[end-1,:]
    elseif settings["BC"]["X+"]["type"] == "insulated"
        u[end,:] .= u[end-1,:]
    end

    if settings["BC"]["Y-"]["type"] == "constant_T"
        u[:,1] .= 2*settings["BC"]["Y-"]["value"] .- u[:,2]
    elseif settings["BC"]["Y-"]["type"] == "insulated"
        u[:,1] .= u[:,2]
    end
    if settings["BC"]["Y+"]["type"] == "constant_T"
        u[:,end] .= 2*settings["BC"]["Y+"]["value"] .- u[:,end-1]
    elseif settings["BC"]["Y+"]["type"] == "insulated"
        u[:,end] .= u[:,end-1]
    end

    @inbounds for i in range(2, N-1)
        @inbounds for j in range(2, N-1)
            # ip1, im1, jp1, jm1 = limit(i + 1, N), limit(i - 1, N), limit(j + 1, N), limit(j - 1, N)
            ip1, im1, jp1, jm1 = i + 1, i - 1, j + 1, j - 1
            # if im1 == 0
            #     du[i, j] = alpha * (u[im1, j] + u[ip1, j] + u[i, jp1] + u[i, jm1] - 4u[i, j]) + source[i, j]
            # elseif ip1 == x_end
            #     du[i, j] = alpha * (u[im1, j] + u[ip1, j] + u[i, jp1] + u[i, jm1] - 4u[i, j]) + source[i, j]
            # elseif 
            du[i, j] = alpha * (u[im1, j] + u[ip1, j] + u[i, jp1] + u[i, jm1] - 4u[i, j]) + source[i, j]
        end
    end
    # du[:,1] = du[:, end] .= 1

end
p = (source, 1, step(xy_mesh))

# function init_conduction_2d!(xyd)
#     N = length(xyd)
#     u = zeros(N, N)
#     return u
# end
# u0 = init_conduction_2d!(xy_mesh)
# prob_conduction_2d = ODEProblem(conduction_2d_loop!, u0, (0.0, 10), p)


du0 = copy(u0)
jac_sparsity = Symbolics.jacobian_sparsity((du, u) -> conduction_2d_loop!(du, u, p, 0.0),
                                           du0, u0)

f = ODEFunction(conduction_2d_loop!; jac_prototype = float.(jac_sparsity))
prob_conduction_2d_sparse = ODEProblem(f, u0, (0.0, 10.0), p, saveat=1.0)
# using BenchmarkTools # for @btime
# solve(prob_conduction_2d, TRBDF2(), saveat=1.0);
# sol = solve(prob_conduction_2d_sparse, TRBDF2(), save_everystep = false, saveat=1.0)
sol = solve(prob_conduction_2d_sparse, 
            KenCarp47(linsolve = KLUFactorization()),
            save_everystep = false, 
            saveat=1.0);

println(Dates.format(now(), "HH:MM:SS"))


plot(contour(x=xy_mesh,y=xy_mesh,z=sol[end][2:end-1,2:end-1]', colorscale="Jet",colorbar=attr(
    title="Temperature", # title here
    titleside="right",
    titlefont=attr(
        size=14,
        family="Arial, sans-serif"
    )
)),Layout(autosize=true))