using Test
using JSON
using therml_environment

f = open("./3d/therml_environment/test/settings_helper_functions.json", "r")
settings = JSON.parse(f)

X,Y,Z = therml_environment.get_vertices_from_bodies(settings)

x_start_solder, x_end_solder, x_start_substrate, x_end_substrate, x_start_bumps, x_end_bumps, x_start_underfill, x_end_underfill, x_start_die, x_end_die, x_start_mold, x_end_mold = X
y_start_solder, y_end_solder, y_start_substrate, y_end_substrate, y_start_bumps, y_end_bumps, y_start_underfill, y_end_underfill, y_start_die, y_end_die, y_start_mold, y_end_mold = Y
z_start_solder, z_end_solder, z_start_substrate, z_end_substrate, z_start_mold, z_end_mold, z_start_bumps, z_end_bumps, z_start_underfill, z_end_underfill, z_start_die, z_end_die = Z

@test x_start_die == -0.005
@test x_end_die == 0.005 
@test x_start_mold == -0.0075 
@test x_end_mold == 0.0075 
@test x_start_underfill == -0.006 
@test x_end_underfill == 0.006 
@test x_start_bumps == -0.005 
@test x_end_bumps == 0.005 
@test x_start_solder == -0.0075 
@test x_end_solder == 0.0075 
@test x_start_substrate == -0.0075 
@test x_end_substrate == 0.0075

@test y_start_die == -0.005
@test y_end_die == 0.005 
@test y_start_mold == -0.0075 
@test y_end_mold == 0.0075 
@test y_start_underfill == -0.006 
@test y_end_underfill == 0.006 
@test y_start_bumps == -0.005 
@test y_end_bumps == 0.005 
@test y_start_solder == -0.0075 
@test y_end_solder == 0.0075 
@test y_start_substrate == -0.0075 
@test y_end_substrate == 0.0075

@test z_start_die == 0.009
@test z_end_die == 0.012 
@test z_start_mold == 0.007
@test z_end_mold == 0.012 
@test z_start_underfill == 0.007 
@test z_end_underfill == 0.011
@test z_start_bumps == 0.007
@test z_end_bumps == 0.009
@test z_start_solder == 0.0 
@test z_end_solder == 0.002 
@test z_start_substrate == 0.002
@test z_end_substrate == 0.007

k_by_rho_cp_actual = load("./3d/therml_environment/test/k_by_rho_cp.jld")["k_by_rho_cp"]

(X_nodes, Y_nodes, Z_nodes) , (Nx,Ny,Nz) = therml_environment.generate_mesh(settings)
k_by_rho_cp = therml_environment.get_k_by_rho_cp((X_nodes, Y_nodes, Z_nodes), (Nx, Ny, Nz), settings)

@test k_by_rho_cp == k_by_rho_cp_actual