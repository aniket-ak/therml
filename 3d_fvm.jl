using DifferentialEquations, LinearAlgebra, SparseArrays
using PlotlyJS
using Symbolics
using JSON
using Dates
include("./functions_fvm_3d.jl")

f = open("settings.json", "r")
settings = JSON.parse(f)

temperature_base_units = settings["units"]["temperature"]

mesh, u0, N = build_domain(settings)

x_mesh, y_mesh, z_mesh = mesh
Nx,Ny,Nz = N

delta_x, delta_y, delta_z = step(x_mesh), step(y_mesh), step(z_mesh)

u0 .= convert_units("temperature", settings["IC"], temperature_base_units, "K")

source = zeros(Nx,Ny,Nz)
source[10:20, 10:20, :] .= 1

println(Dates.format(now(), "HH:MM:SS"))

function conduction_3d_loop!(du, u, p, t)
    source, k, rho, cp, (delta_x, delta_y, delta_z) = p
    alpha_x = (k/(rho*cp)) / delta_x^2
    alpha_y = (k/(rho*cp)) / delta_y^2
    alpha_z = (k/(rho*cp)) / delta_z^2

    @inbounds for i_ in range(2, Nx-1)
        @inbounds for j_ in range(2, Ny-1)
            @inbounds for k_ in range(2, Nz-1)
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

k = settings["model"]["bodies"]["die"]["material"]["k"]
rho = settings["model"]["bodies"]["die"]["material"]["rho"]
cp = settings["model"]["bodies"]["die"]["material"]["cp"]

p = (source, k, rho, cp, (delta_x, delta_y, delta_z))

du0 = zeros(size(u0))
jac_sparsity = Symbolics.jacobian_sparsity((du, u) -> conduction_3d_loop!(du, u, p, 0.0), du0, u0)

f = ODEFunction(conduction_3d_loop!; jac_prototype = float.(jac_sparsity))
prob_conduction_3d_sparse = ODEProblem(f, u0, (0.0, 10.0), p)

sol = solve(prob_conduction_3d_sparse, 
            KenCarp47(linsolve = KLUFactorization()),
            save_everystep = false, 
            saveat=100.0, 
            progress=true, 
            maxiters=1e3, 
            abstol=1e-4, 
            reltol=1e-4);

# prob_conduction_3d = ODEProblem(conduction_3d_loop!, u0, (0,10), p);
# sol = solve(prob_conduction_3d, TRBDF2(), saveat=1.0);

println(Dates.format(now(), "HH:MM:SS"))
result_ = convert_units("temperature", sol[end][2:end-1,2:end-1,1], "K", "C")'

plot(contour(x=xy_mesh,y=xy_mesh,z=result_, colorscale="Jet",colorbar=attr(
    title="Temperature", # title here
    titleside="right",
    titlefont=attr(
        size=14,
        family="Arial, sans-serif"
    )
)),Layout(autosize=true))