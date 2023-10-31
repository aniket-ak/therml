using JSON

function generate_mesh(working_dir)
    f = open(joinpath(working_dir,"settings.json"), "r")
    settings = JSON.parse(f)
    
    bodies_ = settings["model"]["bodies"]

    x_start_solder = -bodies_["solder"]["size"]["X"]/2
    x_end_solder = bodies_["solder"]["size"]["X"]/2

    x_start_substrate = -bodies_["substrate"]["size"]["X"]/2
    x_end_substrate = bodies_["substrate"]["size"]["X"]/2
    
    x_start_bumps = -bodies_["bumps"]["size"]["X"]/2
    x_end_bumps = bodies_["bumps"]["size"]["X"]/2

    x_start_underfill = -bodies_["underfill"]["size"]["X"]/2
    x_end_underfill = bodies_["underfill"]["size"]["X"]/2

    x_start_die = -bodies_["die"]["size"]["X"]/2
    x_end_die = bodies_["die"]["size"]["X"]/2
    
    x_start_mold = -bodies_["mold"]["size"]["X"]/2
    x_end_mold = bodies_["mold"]["size"]["X"]/2

    y_start_solder = -bodies_["solder"]["size"]["Y"]/2
    y_end_solder = bodies_["solder"]["size"]["Y"]/2

    y_start_substrate = -bodies_["substrate"]["size"]["Y"]/2
    y_end_substrate = bodies_["substrate"]["size"]["Y"]/2

    y_start_bumps = -bodies_["bumps"]["size"]["Y"]/2
    y_end_bumps = bodies_["bumps"]["size"]["Y"]/2

    y_start_underfill = -bodies_["underfill"]["size"]["Y"]/2
    y_end_underfill = bodies_["underfill"]["size"]["Y"]/2

    y_start_die = -bodies_["die"]["size"]["Y"]/2
    y_end_die = bodies_["die"]["size"]["Y"]/2

    y_start_mold = -bodies_["mold"]["size"]["Y"]/2
    y_end_mold = bodies_["mold"]["size"]["Y"]/2

    z_start_solder = 0
    z_end_solder = bodies_["solder"]["size"]["Z"]

    z_start_substrate = z_end_solder
    z_end_substrate = z_end_solder + bodies_["substrate"]["size"]["Z"]

    z_start_mold = z_end_substrate
    z_end_mold = z_start_mold + bodies_["mold"]["size"]["Z"]

    z_start_bumps = z_start_mold
    z_end_bumps = z_start_mold + bodies_["bumps"]["size"]["Z"]

    z_start_underfill = z_start_mold
    z_end_underfill = z_start_underfill + bodies_["underfill"]["size"]["Z"]

    z_start_die = z_end_bumps
    z_end_die = z_start_die + bodies_["die"]["size"]["Z"]

    X = [x_start_solder, x_end_solder, x_start_substrate, x_end_substrate, x_start_bumps, x_end_bumps, x_start_underfill, x_end_underfill, x_start_die, x_end_die, x_start_mold, x_end_mold]
    Y = [y_start_solder, y_end_solder, y_start_substrate, y_end_substrate, y_start_bumps, y_end_bumps, y_start_underfill, y_end_underfill, y_start_die, y_end_die, y_start_mold, y_end_mold]
    Z = [z_start_solder, z_end_solder, z_start_substrate, z_end_substrate, z_start_mold, z_end_mold, z_start_bumps, z_end_bumps, z_start_underfill, z_end_underfill, z_start_die, z_end_die]

    # println("\n-----------------------------------------")
    # println("Domain extents : ")
    # println("Direction\t","Min           \t", "Max")
    # println("X        \t", minimum(X) ," \t", maximum(X))
    # println("Y        \t", minimum(Y) ," \t", maximum(Y))
    # println("Z        \t", minimum(Z) ," \t", maximum(Z))
    # println("-----------------------------------------\n")

    X, Y, Z = sort(unique(X)), sort(unique(Y)), sort(unique(Z))

    ghost_x_min = X[1] - bodies_["die"]["mesh"]["size"]["dx"]
    ghost_x_max = X[end] + bodies_["die"]["mesh"]["size"]["dx"]
    ghost_y_min = Y[1] - bodies_["die"]["mesh"]["size"]["dy"]
    ghost_y_max = Y[end] - bodies_["die"]["mesh"]["size"]["dy"]
    ghost_z_min = Z[1] - bodies_["die"]["mesh"]["size"]["dz"]
    ghost_z_max =  Z[end] - bodies_["die"]["mesh"]["size"]["dz"]
    

    X_mesh = [ghost_x_min, X[1],X[end], ghost_x_max]
    Y_mesh = [ghost_y_min, Y[1],Y[end], ghost_y_max]
    Z_mesh = [ghost_z_min, Z[1],Z[end], ghost_z_max]

    for (i, mesh_point) in enumerate(X)
        if i > 1
            mid_point = (mesh_point+X[i-1])/2

            largest_precedence, sizing = 0, 0
            if mid_point > x_start_solder && mid_point < x_end_solder
                if largest_precedence < bodies_["solder"]["mesh"]["precedence"]
                    largest_precedence = bodies_["solder"]["mesh"]["precedence"]
                    sizing = bodies_["solder"]["mesh"]["size"]["dx"]
                    # println("Solder  - ",i, "\t", mid_point, "\t", largest_precedence, "\t", sizing)
                end
            end

            if mid_point > x_start_substrate && mid_point < x_end_substrate
                if largest_precedence < bodies_["substrate"]["mesh"]["precedence"]
                    largest_precedence = bodies_["substrate"]["mesh"]["precedence"]
                    sizing = bodies_["substrate"]["mesh"]["size"]["dx"]
                    # println("substrate  - ",i, "\t", mid_point, "\t", largest_precedence, "\t", sizing)
                end
            end

            if mid_point > x_start_bumps && mid_point < x_end_bumps
                if largest_precedence < bodies_["bumps"]["mesh"]["precedence"]
                    largest_precedence = bodies_["bumps"]["mesh"]["precedence"]
                    sizing = bodies_["bumps"]["mesh"]["size"]["dx"]
                    # println("bumps  - ",i, "\t", mid_point, "\t", largest_precedence, "\t", sizing)
                end
            end

            if mid_point > x_start_underfill && mid_point < x_end_underfill
                if largest_precedence < bodies_["underfill"]["mesh"]["precedence"]
                    largest_precedence = bodies_["underfill"]["mesh"]["precedence"]
                    sizing = bodies_["underfill"]["mesh"]["size"]["dx"]
                    # println("underfill  - ",i, "\t", mid_point, "\t", largest_precedence, "\t", sizing)
                end
            end

            if mid_point > x_start_die && mid_point < x_end_die
                if largest_precedence < bodies_["die"]["mesh"]["precedence"]
                    largest_precedence = bodies_["die"]["mesh"]["precedence"]
                    sizing = bodies_["die"]["mesh"]["size"]["dx"]
                    # println("die  - ",i, "\t", mid_point, "\t", largest_precedence, "\t", sizing)
                end
            end

            if mid_point > x_start_mold && mid_point < x_end_mold
                if largest_precedence < bodies_["mold"]["mesh"]["precedence"]
                    largest_precedence = bodies_["mold"]["mesh"]["precedence"]
                    sizing = bodies_["mold"]["mesh"]["size"]["dx"]
                    # println("mold  - ",i, "\t", mid_point, "\t", largest_precedence, "\t", sizing)
                end
            end

            min_num_nodes = max(2, round(Int, abs(mesh_point-X[i-1])/sizing))
            # println(X[i-1],"\t", mesh_point,"\t", min_num_nodes, "\t", sizing)
            mesh_ = range(X[i-1], mesh_point, min_num_nodes)|>collect
            # println(X[i-1], "\t", mesh_point, "\t", sizing, "\t", min_num_nodes, "\t", mesh_[end])
            append!(X_mesh, mesh_)
        end

    end

    # println("Mesh for X done")

    for (i, mesh_point) in enumerate(Y)
        if i > 1
            mid_point = (mesh_point+Y[i-1])/2

            largest_precedence, sizing = 0, 0
            if mid_point > y_start_solder && mid_point < y_end_solder
                if largest_precedence < bodies_["solder"]["mesh"]["precedence"]
                    largest_precedence = bodies_["solder"]["mesh"]["precedence"]
                    sizing = bodies_["solder"]["mesh"]["size"]["dy"]
                    # println("Solder  - ",i, "\t", mid_point, "\t", largest_precedence, "\t", sizing)
                end
            end

            if mid_point > y_start_substrate && mid_point < y_end_substrate
                if largest_precedence < bodies_["substrate"]["mesh"]["precedence"]
                    largest_precedence = bodies_["substrate"]["mesh"]["precedence"]
                    sizing = bodies_["substrate"]["mesh"]["size"]["dy"]
                    # println("substrate  - ",i, "\t", mid_point, "\t", largest_precedence, "\t", sizing)
                end
            end

            if mid_point > y_start_bumps && mid_point < y_end_bumps
                if largest_precedence < bodies_["bumps"]["mesh"]["precedence"]
                    largest_precedence = bodies_["bumps"]["mesh"]["precedence"]
                    sizing = bodies_["bumps"]["mesh"]["size"]["dy"]
                    # println("bumps  - ",i, "\t", mid_point, "\t", largest_precedence, "\t", sizing)
                end
            end

            if mid_point > y_start_underfill && mid_point < y_end_underfill
                if largest_precedence < bodies_["underfill"]["mesh"]["precedence"]
                    largest_precedence = bodies_["underfill"]["mesh"]["precedence"]
                    sizing = bodies_["underfill"]["mesh"]["size"]["dy"]
                    # println("underfill  - ",i, "\t", mid_point, "\t", largest_precedence, "\t", sizing)
                end
            end

            if mid_point > y_start_die && mid_point < y_end_die
                if largest_precedence < bodies_["die"]["mesh"]["precedence"]
                    largest_precedence = bodies_["die"]["mesh"]["precedence"]
                    sizing = bodies_["die"]["mesh"]["size"]["dy"]
                    # println("die  - ",i, "\t", mid_point, "\t", largest_precedence, "\t", sizing)
                end
            end

            if mid_point > y_start_mold && mid_point < y_end_mold
                if largest_precedence < bodies_["mold"]["mesh"]["precedence"]
                    largest_precedence = bodies_["mold"]["mesh"]["precedence"]
                    sizing = bodies_["mold"]["mesh"]["size"]["dy"]
                    # println("mold  - ",i, "\t", mid_point, "\t", largest_precedence, "\t", sizing)
                end
            end

            min_num_nodes = max(2, round(Int, abs(mesh_point-Y[i-1])/sizing))
            mesh_ = range(Y[i-1], mesh_point, min_num_nodes)|>collect
            # println(X[i-1], "\t", mesh_point, "\t", sizing, "\t", min_num_nodes, "\t", mesh_[end])
            append!(Y_mesh, mesh_)
        end

    end

    # println("Mesh for Y done")

    for (i, mesh_point) in enumerate(Z)
        if i > 1
            mid_point = (mesh_point+Z[i-1])/2

            largest_precedence, sizing = 0, 0
            if mid_point > z_start_solder && mid_point < z_end_solder
                if largest_precedence < bodies_["solder"]["mesh"]["precedence"]
                    largest_precedence = bodies_["solder"]["mesh"]["precedence"]
                    sizing = bodies_["solder"]["mesh"]["size"]["dz"]
                    # println("Solder  - ",i, "\t", mid_point, "\t", largest_precedence, "\t", sizing)
                end
            end

            if mid_point > z_start_substrate && mid_point < z_end_substrate
                if largest_precedence < bodies_["substrate"]["mesh"]["precedence"]
                    largest_precedence = bodies_["substrate"]["mesh"]["precedence"]
                    sizing = bodies_["substrate"]["mesh"]["size"]["dz"]
                    # println("substrate  - ",i, "\t", mid_point, "\t", largest_precedence, "\t", sizing)
                end
            end

            if mid_point > z_start_bumps && mid_point < z_end_bumps
                if largest_precedence < bodies_["bumps"]["mesh"]["precedence"]
                    largest_precedence = bodies_["bumps"]["mesh"]["precedence"]
                    sizing = bodies_["bumps"]["mesh"]["size"]["dz"]
                    # println("bumps  - ",i, "\t", mid_point, "\t", largest_precedence, "\t", sizing)
                end
            end

            if mid_point > z_start_underfill && mid_point < z_end_underfill
                if largest_precedence < bodies_["underfill"]["mesh"]["precedence"]
                    largest_precedence = bodies_["underfill"]["mesh"]["precedence"]
                    sizing = bodies_["underfill"]["mesh"]["size"]["dz"]
                    # println("underfill  - ",i, "\t", mid_point, "\t", largest_precedence, "\t", sizing)
                end
            end

            if mid_point > z_start_die && mid_point < z_end_die
                if largest_precedence < bodies_["die"]["mesh"]["precedence"]
                    largest_precedence = bodies_["die"]["mesh"]["precedence"]
                    sizing = bodies_["die"]["mesh"]["size"]["dz"]
                    # println("die  - ",i, "\t", mid_point, "\t", largest_precedence, "\t", sizing)
                end
            end

            if mid_point > z_start_mold && mid_point < z_end_mold
                if largest_precedence < bodies_["mold"]["mesh"]["precedence"]
                    largest_precedence = bodies_["mold"]["mesh"]["precedence"]
                    sizing = bodies_["mold"]["mesh"]["size"]["dz"]
                    # println("mold  - ",i, "\t", mid_point, "\t", largest_precedence, "\t", sizing)
                end
            end

            min_num_nodes = max(2, round(Int, abs(mesh_point-Z[i-1])/sizing))
            # println(Z[i-1], "\t", mesh_point, "\t", sizing, "\t", min_num_nodes)
            mesh_ = range(Z[i-1], mesh_point, min_num_nodes)|>collect
            append!(Z_mesh, mesh_)
        end

    end


    nodes = X_mesh, Y_mesh, Z_mesh = sort(unique(X_mesh)), sort(unique(Y_mesh)), sort(unique(Z_mesh))

    num_elements_x, num_elements_y, num_elements_z = (size(X_mesh)[1]-1, size(Y_mesh)[1]-1, size(Z_mesh)[1]-1)
    elements = zeros(num_elements_x, num_elements_y, num_elements_z)

    return nodes, elements
end

nodes, elements = generate_mesh("/Users/aniket/Documents/MarlinSim/04_testing/scenarios/")
