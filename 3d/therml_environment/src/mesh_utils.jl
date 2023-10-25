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

mutable struct Node
    x::Float64
    y::Float64
    z::Float64
    id::Int

    Node(x, y, z) = new(x, y, z, 1)
end

Base.isequal(n1::Node, n2::Node) = n1.x == n2.x && n1.y == n2.y && n1.z == n2.z
Base.hash(n::Node, h::UInt) = hash((n.x, n.y, n.z), h)

mutable struct HexElement
    nodes::Vector{Node}
end

mutable struct Mesh
    nodes::Vector{Node}
    elements::Vector{HexElement}
    
    Mesh() = new(Node[], HexElement[])
end

function generate_box_mesh(box::Box)
    mesh = Mesh()
    bounds = get_bounds(box)
    divisions = (box.mesh_details.elements_per_edge_x, box.mesh_details.elements_per_edge_y, box.mesh_details.elements_per_edge_z)
    
    dx = (bounds[:xmax] - bounds[:xmin]) / divisions[1]
    dy = (bounds[:ymax] - bounds[:ymin]) / divisions[2]
    dz = (bounds[:zmax] - bounds[:zmin]) / divisions[3]
    
    node_grid = []

    for i in 1:(divisions[1] + 1)
        push!(node_grid, [])
        for j in 1:(divisions[2] + 1)
            push!(node_grid[end], [])
            for k in 1:(divisions[3] + 1)
                x = bounds[:xmin] + (i - 1) * dx
                y = bounds[:ymin] + (j - 1) * dy
                z = bounds[:zmin] + (k - 1) * dz
                node = Node(x, y, z)
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
                push!(mesh.elements, hex_element)
            end
        end
    end

    return mesh
end

function assemble_meshes(boxes::Vector{Box})
    global_nodes = Node[]
    global_mesh = Mesh()
    
    for box in boxes
        box_mesh = generate_box_mesh(box)
        
        # Check nodes and avoid duplicates
        for node in box_mesh.nodes
            println("node --> ", node)
            idx = findfirst(isequal(node), global_nodes)
            if isnothing(idx)
                push!(global_nodes, node)
            else
                node = global_nodes[idx]
            end
            push!(global_mesh.nodes, node)
        end
        
        # Add elements
        for elem in box_mesh.elements
            push!(global_mesh.elements, elem)
        end
    end

    return global_mesh
end

# Example usage:
vertex_a = Vertex(0.0, 0.0, 0.0)
vertex_b = Vertex(1.0, 1.0, 1.0)
mat_ = MaterialProperties(0.3,1000,10)
mesh_ = MeshDetails(2,2,2)
box_1 = Box(vertex_a, vertex_b, mat_, mesh_)

box_2 = Box(Vertex(0.5, 0.0, 0.0), Vertex(1.0, 1.0, 1.0), mat_, mesh_)
# println(box_1)

mesh_ = assemble_meshes([box_1, box_2])