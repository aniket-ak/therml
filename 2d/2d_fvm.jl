using DifferentialEquations, LinearAlgebra, SparseArrays
using PlotlyJS
using Symbolics
using JSON
using Dates
using CUDA
include("./functions_fvm.jl")

f = open("settings.json", "r")
settings = JSON.parse(f)

temperature_base_units = settings["units"]["temperature"]
num_sources = settings["model"]["num_sources"]
xy_mesh, u0, N = build_domain(settings, num_sources)
mesh_size = step(xy_mesh)

u0 .= cu(convert_units("temperature", settings["IC"], temperature_base_units, "K"))

source = zeros(N,N)
source[10:20, 10:20] .= 1

println(Dates.format(now(), "HH:MM:SS"))

function conduction_2d_loop!(du, u, p, t)
    source, alpha, k, rho, cp, mesh_size = p
    alpha = alpha / mesh_size^2

    @inbounds for i in range(2, N-1)
        @inbounds for j in range(2, N-1)
            ip1, im1, jp1, jm1 = i + 1, i - 1, j + 1, j - 1
            du[i, j] = alpha * (u[im1, j] + u[ip1, j] + u[i, jp1] + u[i, jm1] - 4u[i, j]) + source[i, j]
        end
    end

    apply_bc(u, settings, mesh_size)

end

k = settings["model"]["bodies"]["die"]["material"]["k"]
rho = settings["model"]["bodies"]["die"]["material"]["rho"]
cp = settings["model"]["bodies"]["die"]["material"]["cp"]
alpha = k/(rho*cp)

p = (source, alpha, k, rho, cp, mesh_size)

du0 = cu(zeros(size(u0)))

#jac_sparsity = Symbolics.jacobian_sparsity((du, u) -> conduction_2d_loop!(du, u, p, 0.0), du0, u0)

#f = ODEFunction(conduction_2d_loop!; jac_prototype = float.(jac_sparsity))
#prob_conduction_2d_sparse = ODEProblem(f, u0, (0.0, 1000.0), p)

#sol = solve(prob_conduction_2d_sparse, 
#            KenCarp47(linsolve = KLUFactorization()),
#            save_everystep = false, 
#            saveat=100.0);
problem = ODEProblem(conduction_2d_loop!, u0, (0.0f0,10.0f0), p)
sol = solve(problem)

println(Dates.format(now(), "HH:MM:SS"))
#result_ = convert_units("temperature", sol[end][2:end-1,2:end-1], "K", "C")'

#plot(contour(x=xy_mesh,y=xy_mesh,z=result_, colorscale="Jet",colorbar=attr(
#    title="Temperature", # title here
#    titleside="right",
#    titlefont=attr(
#       size=14,
#        family="Arial, sans-serif"
#    )
#)),Layout(autosize=true))
