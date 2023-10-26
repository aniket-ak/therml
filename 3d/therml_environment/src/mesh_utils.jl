struct Vertex
    x::Float64
    y::Float64
    z::Float64
end

struct Edge
    vertex1::Vertex
    vertex2::Vertex
end

struct Surface
    vertices::Vector{Vertex}
    function Surface(vertices::Vector{Vertex})
        if length(vertices) != 4
            throw(ArgumentError("A surface requires 4 vertices."))
        end
        return new(vertices)
    end
end

struct MaterialProperties
    k::Float64
    rho::Float64
    cp::Float64
end

mutable struct MeshDetails
    elements_per_edge_x::Int
    elements_per_edge_y::Int
    elements_per_edge_z::Int
    precedence::Int
end

mutable struct Box
    vertices::Vector{Vertex}
    edges::Vector{Edge}
    faces::Vector{Surface}
    material_properties::MaterialProperties
    mesh_details::MeshDetails

    function Box(vertex1::Vertex, vertex2::Vertex, material_properties::MaterialProperties, mesh_details::MeshDetails)
        vertices = [
            vertex1,
            Vertex(vertex1.x, vertex1.y, vertex2.z),
            Vertex(vertex1.x, vertex2.y, vertex1.z),
            Vertex(vertex1.x, vertex2.y, vertex2.z),
            Vertex(vertex2.x, vertex1.y, vertex1.z),
            Vertex(vertex2.x, vertex1.y, vertex2.z),
            Vertex(vertex2.x, vertex2.y, vertex1.z),
            vertex2
        ]

        edges = [
            Edge(vertices[1], vertices[2]),
            Edge(vertices[1], vertices[3]),
            Edge(vertices[1], vertices[5]),
            Edge(vertices[4], vertices[2]),
            Edge(vertices[4], vertices[3]),
            Edge(vertices[4], vertices[8]),
            Edge(vertices[7], vertices[3]),
            Edge(vertices[7], vertices[5]),
            Edge(vertices[7], vertices[8]),
            Edge(vertices[6], vertices[2]),
            Edge(vertices[6], vertices[5]),
            Edge(vertices[6], vertices[8]),
        ]

        faces = [
            Surface([vertices[1], vertices[2], vertices[4], vertices[3]]),
            Surface([vertices[1], vertices[2], vertices[6], vertices[5]]),
            Surface([vertices[1], vertices[3], vertices[7], vertices[5]]),
            Surface([vertices[8], vertices[4], vertices[2], vertices[6]]),
            Surface([vertices[8], vertices[4], vertices[3], vertices[7]]),
            Surface([vertices[8], vertices[6], vertices[5], vertices[7]])
        ]

        return new(vertices, edges, faces, material_properties, mesh_details)
    end
end

function get_bounds(b::Box)
    return Dict(
        :xmin => minimum(v.x for v in b.vertices),
        :xmax => maximum(v.x for v in b.vertices),
        :ymin => minimum(v.y for v in b.vertices),
        :ymax => maximum(v.y for v in b.vertices),
        :zmin => minimum(v.z for v in b.vertices),
        :zmax => maximum(v.z for v in b.vertices)
    )
end

function intersects(b1::Box, b2::Box)
    bounds1 = get_bounds(b1)
    bounds2 = get_bounds(b2)

    if bounds1[:xmax] < bounds2[:xmin] || bounds1[:xmin] > bounds2[:xmax]
        return false
    end
    if bounds1[:ymax] < bounds2[:ymin] || bounds1[:ymin] > bounds2[:ymax]
        return false
    end
    if bounds1[:zmax] < bounds2[:zmin] || bounds1[:zmin] > bounds2[:zmax]
        return false
    end

    return true
end

struct Node
    x::Float64
    y::Float64
    z::Float64
end

mutable struct HexElement
    nodes::Vector{Node}
end

mutable struct Mesh
    nodes::Vector{Node}
    elements::Vector{HexElement}
end

function generate_box_mesh(box::Box)
    bounds = get_bounds(box)
    local_mesh = Mesh([], [])
    divisions = (box.mesh_details.elements_per_edge_x, box.mesh_details.elements_per_edge_y, box.mesh_details.elements_per_edge_z)

    dx = (bounds[:xmax] - bounds[:xmin]) / divisions[1]
    dy = (bounds[:ymax] - bounds[:ymin]) / divisions[2]
    dz = (bounds[:zmax] - bounds[:zmin]) / divisions[3]

    node_grid = []
    for i in 1:divisions[1] + 1
        push!(node_grid, [])
        for j in 1:divisions[2] + 1
            push!(node_grid[end], [])
            for k in 1:divisions[3] + 1
                x = bounds[:xmin] + (i-1) * dx
                y = bounds[:ymin] + (j-1) * dy
                z = bounds[:zmin] + (k-1) * dz

                node = Node(x, y, z)
                push!(local_mesh.nodes, node)
                push!(node_grid[end][end], node)
            end
        end
    end

    for i in 1:divisions[1]
        for j in 1:divisions[2]
            for k in 1:divisions[3]
                nodes = [
                    node_grid[i][j][k],
                    node_grid[i+1][j][k],
                    node_grid[i+1][j+1][k],
                    node_grid[i][j+1][k],
                    node_grid[i][j][k+1],
                    node_grid[i+1][j][k+1],
                    node_grid[i+1][j+1][k+1],
                    node_grid[i][j+1][k+1]
                ]
                hex_element = HexElement(nodes)
                push!(local_mesh.elements, hex_element)
            end
        end
    end

    return local_mesh
end

import Base: ==

function ==(n1::Node, n2::Node)
    return isapprox(n1.x, n2.x) && isapprox(n1.y, n2.y) && isapprox(n1.z, n2.z)
end

function ==(e1::HexElement, e2::HexElement)
    return all(e1.nodes .== e2.nodes)
end


# Define hash for Node based on its coordinates
function Base.hash(n::Node, h::UInt)
    return hash(n.x, hash(n.y, hash(n.z, h)))
end

# Define hash for HexElement based on the hash of its nodes
function Base.hash(e::HexElement, h::UInt)
    node_hashes = [hash(node) for node in e.nodes]
    return foldl(hash, node_hashes, init=h)
end


function assemble_meshes(boxes::Vector{Box})
    global_mesh = Mesh([], [])

    for box in boxes
        box_mesh = generate_box_mesh(box)

        # Check node duplicates and add to global mesh
        for node in box_mesh.nodes
            if !any(existing_node -> existing_node == node, global_mesh.nodes)
                push!(global_mesh.nodes, node)
            end
        end

        # Just add the elements directly to the global mesh
        for element in box_mesh.elements
            if !any(existing_element -> existing_element == element, global_mesh.elements)
                push!(global_mesh.elements, element)
            end
        end
    end

    return global_mesh
end

function centroid(element::HexElement)
    avg_x = sum(node.x for node in element.nodes) / 8.0
    avg_y = sum(node.y for node in element.nodes) / 8.0
    avg_z = sum(node.z for node in element.nodes) / 8.0
    return avg_x, avg_y, avg_z
end


# Example usage:
vertex_a = Vertex(0.0, 0.0, 0.0);
vertex_b = Vertex(1.0, 1.0, 1.0);
mat_ = MaterialProperties(0.3,1000,10);
mesh_ = MeshDetails(2,2,2,100)
box_1 = Box(vertex_a, vertex_b, mat_, mesh_);

box_2 = Box(Vertex(1, 0.0, 0.0), Vertex(2.0, 2.0, 2.0), mat_, MeshDetails(4,2,2,100))
# println(box_1)

mesh_ = assemble_meshes([box_1, box_2])

cell_grid = []
centroids_x, centroids_y, centroids_z = [], [], []
for element in mesh_.elements
    e_ = centroid_x, centroid_y, centroid_z = centroid(element)
    push!(cell_grid, e_)
    push!(centroids_x, centroid_x)
    push!(centroids_y, centroid_y)
    push!(centroids_z, centroid_z)
end

global_matrix_size = (size(unique(centroids_x))[1], size(unique(centroids_y))[1], size(unique(centroids_z))[1])

function split_cube(box_1::Box, box_2::Box)
    box_1_bounds = get_bounds(box_1)
    box_2_bounds = get_bounds(box_2)

    if !intersects(box_1, box_2)
        return box_1,box_2
    else
        "a"
    end
end

function generate_mesh(boxes::Vector{Box})
    println(boxes)
    
end