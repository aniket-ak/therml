import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import plotly.graph_objects as go

class Vertex:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def __repr__(self):
        return f"Vertex({self.x}, {self.y}, {self.z})"


class Edge:
    def __init__(self, vertex1, vertex2):
        self.vertex1 = vertex1
        self.vertex2 = vertex2

    def __repr__(self):
        return f"Edge({self.vertex1}, {self.vertex2})"


class Surface:
    def __init__(self, vertices):
        if len(vertices) != 4:
            raise ValueError("A surface requires 4 vertices.")
        self.vertices = vertices

    def __repr__(self):
        return f"Surface({', '.join(map(str, self.vertices))})"


class Box:
    def __init__(self, vertex1, vertex2):
        # Define the 8 vertices of the box based on the provided diagonally opposite corners
        self.vertices = [
            vertex1,
            Vertex(vertex1.x, vertex1.y, vertex2.z),
            Vertex(vertex1.x, vertex2.y, vertex1.z),
            Vertex(vertex1.x, vertex2.y, vertex2.z),
            Vertex(vertex2.x, vertex1.y, vertex1.z),
            Vertex(vertex2.x, vertex1.y, vertex2.z),
            Vertex(vertex2.x, vertex2.y, vertex1.z),
            vertex2
        ]

        # Define the 12 edges of the box
        self.edges = [
            Edge(self.vertices[0], self.vertices[1]),
            Edge(self.vertices[0], self.vertices[2]),
            Edge(self.vertices[0], self.vertices[4]),
            Edge(self.vertices[3], self.vertices[1]),
            Edge(self.vertices[3], self.vertices[2]),
            Edge(self.vertices[3], self.vertices[7]),
            Edge(self.vertices[6], self.vertices[2]),
            Edge(self.vertices[6], self.vertices[4]),
            Edge(self.vertices[6], self.vertices[7]),
            Edge(self.vertices[5], self.vertices[1]),
            Edge(self.vertices[5], self.vertices[4]),
            Edge(self.vertices[5], self.vertices[7]),
        ]

        # Define the 6 faces of the box
        self.faces = [
            Surface([self.vertices[0], self.vertices[1], self.vertices[3], self.vertices[2]]),
            Surface([self.vertices[0], self.vertices[1], self.vertices[5], self.vertices[4]]),
            Surface([self.vertices[0], self.vertices[2], self.vertices[6], self.vertices[4]]),
            Surface([self.vertices[7], self.vertices[3], self.vertices[1], self.vertices[5]]),
            Surface([self.vertices[7], self.vertices[3], self.vertices[2], self.vertices[6]]),
            Surface([self.vertices[7], self.vertices[5], self.vertices[4], self.vertices[6]])
        ]

    def __repr__(self):
        return f"Box(Vertices: {self.vertices}, Edges: {self.edges}, Faces: {self.faces})"
    
    def get_bounds(self):
        """Return the min and max bounds of the box."""
        return {
            'xmin': min(v.x for v in self.vertices),
            'xmax': max(v.x for v in self.vertices),
            'ymin': min(v.y for v in self.vertices),
            'ymax': max(v.y for v in self.vertices),
            'zmin': min(v.z for v in self.vertices),
            'zmax': max(v.z for v in self.vertices),
        }

    def intersects(self, other_box):
        """Check if this box intersects with another box."""
        bounds1 = self.get_bounds()
        bounds2 = other_box.get_bounds()

        # Check if boxes are separated along x, y, or z axis
        if bounds1['xmax'] < bounds2['xmin'] or bounds1['xmin'] > bounds2['xmax']:
            return False
        if bounds1['ymax'] < bounds2['ymin'] or bounds1['ymin'] > bounds2['ymax']:
            return False
        if bounds1['zmax'] < bounds2['zmin'] or bounds1['zmin'] > bounds2['zmax']:
            return False

        return True


# Example usage:
vertex_a = Vertex(0, 0, 0)
vertex_b = Vertex(1, 1, 1)
box_1 = Box(vertex_a, vertex_b)
print(box_1)

# Refactored code for mesh generation

class Node:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def __repr__(self):
        return f"Node({self.x}, {self.y}, {self.z})"
    
    def __eq__(self, other):
        return self.x == other.x and self.y == other.y and self.z == other.z

    def __hash__(self):
        return hash((self.x, self.y, self.z))
    
class HexElement:
    def __init__(self, nodes):
        if len(nodes) != 8:
            raise ValueError("A hexahedral element requires 8 nodes.")
        self.nodes = nodes

    def __repr__(self):
        node_indices = [node.id for node in self.nodes]
        return f"HexElement({node_indices})"

class Mesh:
    def __init__(self):
        self.nodes = []
        self.elements = []

    def add_node(self, node):
        node.id = len(self.nodes) + 1  # Assign a unique ID to the node
        self.nodes.append(node)

    def add_element(self, element):
        self.elements.append(element)

def generate_box_mesh(box, divisions=(4, 4, 4), global_nodes=None):
    global_nodes = set()
    bounds = box.get_bounds()
    mesh = Mesh()

    # If global_nodes isn't provided, initialize it
    if global_nodes is None:
        global_nodes = set()

    dx = (bounds['xmax'] - bounds['xmin']) / divisions[0]
    dy = (bounds['ymax'] - bounds['ymin']) / divisions[1]
    dz = (bounds['zmax'] - bounds['zmin']) / divisions[2]

    node_grid = []
    for i in range(divisions[0] + 1):
        node_grid.append([])
        for j in range(divisions[1] + 1):
            node_grid[i].append([])
            for k in range(divisions[2] + 1):
                x = bounds['xmin'] + i * dx
                y = bounds['ymin'] + j * dy
                z = bounds['zmin'] + k * dz

                # Check if a node already exists at this location
                node = Node(x, y, z)
                if node in global_nodes:
                    node = [n for n in global_nodes if n == node][0]
                else:
                    global_nodes.add(node)
                    mesh.add_node(node)

                node_grid[i][j].append(node)

    for i in range(divisions[0]):
        for j in range(divisions[1]):
            for k in range(divisions[2]):
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
                mesh.add_element(hex_element)

    return mesh

def assemble_meshes(boxes, divisions=(1, 1, 1)):
    global_nodes = set()
    global_mesh = Mesh()

    for box in boxes:
        box_mesh = generate_box_mesh(box, divisions, global_nodes)
        global_mesh.nodes.extend(box_mesh.nodes)
        global_mesh.elements.extend(box_mesh.elements)

    return global_mesh

# Example usage:
# box = Box(Vertex(0, 0, 0), Vertex(2, 2, 2))
# mesh = generate_box_mesh(box, divisions=(2, 2, 2))

# print("Nodes:", mesh.nodes)
# print("Hex elements:", mesh.elements)

box1 = Box(Vertex(0, 0, 0), Vertex(1, 1, 1))
box2 = Box(Vertex(1, 0, 0), Vertex(2, 1, 1))
# boxes = [box1, box2]
# mesh = assemble_meshes(boxes, divisions=(2, 2, 2))

mesh1 = generate_box_mesh(box1, divisions=(2,3,4))
mesh2 = generate_box_mesh(box2, divisions=(4,2,2))

global_nodes = set()
global_mesh = Mesh()

global_mesh.nodes.extend(mesh1.nodes)
global_mesh.elements.extend(mesh1.elements)

global_mesh.nodes.extend(mesh2.nodes)
global_mesh.elements.extend(mesh2.elements)


print(len(global_mesh.elements))

def visualize_mesh_plotly(mesh):
    # For the mesh
    x, y, z = [], [], []
    i, j, k = [], [], []

    # For the edges
    lines = []

    # Loop through each hex element to extract vertices and faces
    for element in mesh.elements:
        # Append node coordinates to lists
        for node in element.nodes:
            x.append(node.x)
            y.append(node.y)
            z.append(node.z)
        
        base_index = len(x) - 8
        face_indices = [
            [0, 1, 5, 4],
            [7, 6, 2, 3],
            [0, 3, 2, 1],
            [7, 4, 5, 6],
            [0, 4, 7, 3],
            [1, 2, 6, 5]
        ]

        # Append face indices to lists
        for face in face_indices:
            i.append(base_index + face[0])
            j.append(base_index + face[1])
            k.append(base_index + face[2])
            i.append(base_index + face[0])
            j.append(base_index + face[2])
            k.append(base_index + face[3])
        
        # Define edges of the hex element
        edges = [
            [0, 1], [1, 2], [2, 3], [3, 0],
            [4, 5], [5, 6], [6, 7], [7, 4],
            [0, 4], [1, 5], [2, 6], [3, 7]
        ]

        # Append the edges
        for edge in edges:
            lines.append(go.Scatter3d(x=[x[base_index + edge[0]], x[base_index + edge[1]]],
                                      y=[y[base_index + edge[0]], y[base_index + edge[1]]],
                                      z=[z[base_index + edge[0]], z[base_index + edge[1]]],
                                      mode='lines',
                                      line=dict(color='red', width=2)))

    # Create mesh3d plot
    fig = go.Figure(
        data=[
            go.Mesh3d(
                x=x, y=y, z=z,
                i=i, j=j, k=k,
                opacity=0.6,
                color='cyan'
            )
        ]
    )

    # Add the lines (edges)
    for line in lines:
        fig.add_trace(line)

    fig.show()

def visualize_mesh(mesh):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # Draw each hex element
    for element in mesh.elements:
        # Extract node coordinates for easier referencing
        coords = [(node.x, node.y, node.z) for node in element.nodes]
        
        # Faces of the hexahedron
        faces = [
            [coords[0], coords[1], coords[5], coords[4]],
            [coords[7], coords[6], coords[2], coords[3]],
            [coords[0], coords[3], coords[2], coords[1]],
            [coords[7], coords[4], coords[5], coords[6]],
            [coords[0], coords[4], coords[7], coords[3]],
            [coords[1], coords[2], coords[6], coords[5]]
        ]
        
        ax.add_collection3d(Poly3DCollection(faces, facecolors='cyan', linewidths=1, edgecolors='r', alpha=0.5))

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()

# Example usage:
# box1 = Box(Vertex(0, 0, 0), Vertex(0.5, 0.5, 0.5))
# box2 = Box(Vertex(0.6, 0, 0), Vertex(1, 1, 1))
# boxes = [box1, box2]
# mesh = assemble_meshes(boxes, divisions=(2, 2, 2))

visualize_mesh_plotly(global_mesh)
