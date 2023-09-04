using DifferentialEquations, LinearAlgebra, SparseArrays
using PlotlyJS
using Symbolics
using JSON
using Dates
include("./functions.jl")

f = open("settings.json", "r")
settings = JSON.parse(f)

N = settings["model"]["num_sources"]
xy_mesh = range(0, stop = 1, length = N)
source = zeros(N,N)
source[2:3, 2:3] .= 1.0

println(Dates.format(now(), "HH:MM:SS"))

limit(a, N) = a == N + 1 ? 1 : a == 0 ? N : a
# if a == N+1:
    # return 1
# elif a == 0:
    # return N
# else:
# return a

x_start, y_start, z_start, x_end, y_end, z_end = get_indices(settings, N)

function conduction_2d_loop!(du, u, p, t)
    source, alpha, dx = p
    alpha = alpha / dx^2
    @inbounds for i in range(x_start, x_end)
        @inbounds for j in range(y_start, y_end)
            ip1, im1, jp1, jm1 = limit(i + 1, N), limit(i - 1, N), limit(j + 1, N), limit(j - 1, N)
            du[i, j] = alpha * (u[im1, j] + u[ip1, j] + u[i, jp1] + u[i, jm1] - 4u[i, j]) + source[i, j]
        end
    end

    du[1,:] = du[end,:] .= 0
    du[:,1] = du[:, end] .= 1

    # if settings["BC"]["X-"]["type"] == "adiabatic":

end
p = (source, 0.001, step(xy_mesh))

function init_conduction_2d!(xyd)
    N = length(xyd)
    u = zeros(N, N)
    # for I in CartesianIndices((N, N))
    #     i,j = Tuple(I)
    #     if i > 45 && i<50 && j<50 && j>45
    #         u[i, j] = 10
    #     else
    #         u[i,j] = 0
    #     end
    # end
    u
end
u0 = init_conduction_2d!(xy_mesh)
prob_conduction_2d = ODEProblem(conduction_2d_loop!, u0, (0.0, 10), p)


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
            save_everystep = true, 
            saveat=1.0);

println(Dates.format(now(), "HH:MM:SS"))


plot(contour(x=xy_mesh,y=xy_mesh,z=sol[end]', colorscale="Jet",colorbar=attr(
    title="Temperature", # title here
    titleside="right",
    titlefont=attr(
        size=14,
        family="Arial, sans-serif"
    )
)),Layout(autosize=true))