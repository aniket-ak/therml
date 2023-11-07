using Test
using JSON
using therml_environment

f = open("./3d/therml_environment/test/settings_mesher_test.json", "r")
settings = JSON.parse(f)

(X_nodes, Y_nodes, Z_nodes) , (Nx,Ny,Nz) = therml_environment.generate_mesh(settings)

@test X_nodes == [-0.006, -0.005, -0.004, -0.003, -0.002, -0.001, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006]
@test Y_nodes == [-0.006, -0.005, -0.004, -0.003, -0.002, -0.001, 0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006]
@test Z_nodes == -0.001:0.001:0.041 |> collect