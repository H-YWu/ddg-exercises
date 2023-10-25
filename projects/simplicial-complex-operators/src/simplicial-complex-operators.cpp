// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

/*
 * Assign a unique index to each vertex, edge, and face of a mesh.
 * All elements are 0-indexed.
 *
 * Input: None. Access geometry via the member variable <geometry>, and pointer to the mesh via <mesh>.
 * Returns: None.
 */
void SimplicialComplexOperators::assignElementIndices() {

    // Needed to access geometry->vertexIndices, etc. as cached quantities.
    // Not needed if you're just using v->getIndex(), etc.
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    // You can set the index field of a vertex via geometry->vertexIndices[v], where v is a Vertex object (or an
    // integer). Similarly you can do edges and faces via geometry->edgeIndices, geometry->faceIndices, like so:
    size_t idx = 0;
    for (Vertex v : mesh->vertices()) {
        idx = geometry->vertexIndices[v];
    }

    for (Edge e : mesh->edges()) {
        idx = geometry->edgeIndices[e];
    }

    for (Face f : mesh->faces()) {
        idx = geometry->faceIndices[f];
    }

    // You can more easily get the indices of mesh elements using the function getIndex(), albeit less efficiently and
    // technically less safe (although you don't need to worry about it), like so:
    //
    //      v.getIndex()
    //
    // where v can be a Vertex, Edge, Face, Halfedge, etc. For example:

    for (Vertex v : mesh->vertices()) {
        idx = v.getIndex(); // == geometry->vertexIndices[v])
    }

    // Geometry Central already sets the indices for us, though, so this function is just here for demonstration.
    // You don't have to do anything :)
}

/*
 * Construct the unsigned vertex-edge adjacency matrix A0.
 *
 * Input:
 * Returns: The sparse vertex-edge adjacency matrix which gets stored in the global variable A0.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildVertexEdgeAdjacencyMatrix() const {

    // Note: You can build an Eigen sparse matrix from triplets, then return it as a Geometry Central SparseMatrix.
    // See <https://eigen.tuxfamily.org/dox/group__TutorialSparse.html> for documentation.

    // construct the triplets of (row, column, 1) where row is edge index and column is vertex index 
    std::vector<Eigen::Triplet<size_t>> coefficients;
    for (Edge e : mesh->edges()) {
        std::array<Vertex, 2> vertices = e.adjacentVertices();
        size_t e_idx = e.getIndex();
        for (Vertex v : vertices) {
            size_t v_idx = v.getIndex();
            Eigen::Triplet<size_t> tri(e_idx, v_idx, 1);
            coefficients.push_back(tri);
        }
    }

    // build the sparse vertex-edge adjacency matrix from triplets
    Eigen::SparseMatrix<size_t> A(mesh->nEdges(), mesh->nVertices());
    A.setFromTriplets(coefficients.begin(), coefficients.end());

    return A;
}

/*
 * Construct the unsigned face-edge adjacency matrix A1.
 *
 * Input:
 * Returns: The sparse face-edge adjacency matrix which gets stored in the global variable A1.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildFaceEdgeAdjacencyMatrix() const {

    // construct the triplets of (row, column, 1) where row is face index and column is edge index 
    std::vector<Eigen::Triplet<size_t>> coefficients;
    for (Face f : mesh->faces()) {
        auto edges = f.adjacentEdges();
        size_t f_idx = f.getIndex();
        for (Edge e : edges) {
            size_t e_idx = e.getIndex();
            Eigen::Triplet<size_t> tri(f_idx, e_idx, 1);
            coefficients.push_back(tri);
        }
    }

    // build the sparse face-edge adjacency matrix from triplets
    Eigen::SparseMatrix<size_t> A(mesh->nFaces(), mesh->nEdges());
    A.setFromTriplets(coefficients.begin(), coefficients.end());

    return A;
}

/*
 * Construct a vector encoding the vertices in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |V|, where |V| = # of vertices in the mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildVertexVector(const MeshSubset& subset) const {

    Vector<size_t> vertices;
    vertices.resize(mesh->nVertices());
    vertices.setZero();
    for (size_t v_idx : subset.vertices) {
        vertices(v_idx) = 1;
    }
    return vertices;
}

/*
 * Construct a vector encoding the edges in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |E|, where |E| = # of edges in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildEdgeVector(const MeshSubset& subset) const {

    Vector<size_t> edges;
    edges.resize(mesh->nEdges());
    edges.setZero();
    for (size_t e_idx : subset.edges) {
        edges(e_idx) = 1;
    }
    return edges;
}

/*
 * Construct a vector encoding the faces in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |F|, where |F| = # of faces in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildFaceVector(const MeshSubset& subset) const {

    Vector<size_t> faces;
    faces.resize(mesh->nFaces());
    faces.setZero();
    for (size_t f_idx : subset.faces) {
        faces(f_idx) = 1;
    }
    return faces;
}

/*
 * Compute the simplicial star St(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The star of the given subset.
 */
MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {

    auto vertexVecter = buildVertexVector(subset);
    auto edgeVecter = buildEdgeVector(subset);
    // edges containing all vertices in subset
    std::set<size_t> stEdges(subset.edges.begin(), subset.edges.end());
    for (auto e : mesh->edges()) {
        size_t e_idx = e.getIndex();
        auto adj_vertices = e.adjacentVertices();
        for (auto v : adj_vertices) {
            size_t v_idx = v.getIndex();
            if (vertexVecter(v_idx) == 1) {
                stEdges.insert(e_idx);
            }
        }
    }
    // faces containing all edges and vertices in subset
    std::set<size_t> stFaces(subset.faces.begin(), subset.faces.end());
    for (auto f : mesh->faces()) {
        size_t f_idx = f.getIndex();
        auto adj_edges = f.adjacentEdges();
        for (auto e : adj_edges) {
            size_t e_idx = e.getIndex();
            if (edgeVecter(e_idx) == 1) {
                stFaces.insert(f_idx);
            }
        }
        auto adj_vertices = f.adjacentVertices();
        for (auto v : adj_vertices) {
            size_t v_idx = v.getIndex();
            if (vertexVecter(v_idx) == 1) {
                stFaces.insert(f_idx);
            }
        }
    }
    return MeshSubset(subset.vertices, stEdges, stFaces);
}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {

    auto Cl = subset;
    // check faces first
    for (size_t f_idx : Cl.faces) {
        Face f = mesh->face(f_idx); 
        auto adj_edges = f.adjacentEdges();
        for (auto e : adj_edges) {
            size_t e_idx = e.getIndex();
            Cl.addEdge(e_idx);
        }
        auto adj_vertices = f.adjacentVertices();
        for (auto v : adj_vertices) {
            size_t v_idx = v.getIndex();
            Cl.addVertex(v_idx);
        }
    }
    // then edges
    for (size_t e_idx : Cl.edges) {
        Edge e = mesh->edge(e_idx);
        auto adj_vertices = e.adjacentVertices();
        for (auto v : adj_vertices) {
            size_t v_idx = v.getIndex();
            Cl.addVertex(v_idx);
        }
    }
    return Cl;
}

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {

    auto StCl = star(closure(subset));
    auto ClSt = closure(star(subset));
    std::set<size_t> LkVertices;
    std::set<size_t> LkEdges;
    std::set<size_t> LkFaces;
    std::set_difference(ClSt.vertices.begin(), ClSt.vertices.end(), StCl.vertices.begin(), StCl.vertices.end(),
        std::inserter(LkVertices, LkVertices.end()));
    std::set_difference(ClSt.edges.begin(), ClSt.edges.end(), StCl.edges.begin(), StCl.edges.end(),
        std::inserter(LkEdges, LkEdges.end()));
    std::set_difference(ClSt.faces.begin(), ClSt.faces.end(), StCl.faces.begin(), StCl.faces.end(),
        std::inserter(LkFaces, LkFaces.end()));
    return MeshSubset(LkVertices, LkEdges, LkFaces);
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */
bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {

    for (size_t f_idx : subset.faces) {
        Face f = mesh->face(f_idx); 
        auto adj_edges = f.adjacentEdges();
        for (auto e : adj_edges) {
            size_t e_idx = e.getIndex();
            if (subset.edges.find(e_idx) != subset.edges.end()) {
                continue;
            }
            else {
                return false;
            }
        }
        auto adj_vertices = f.adjacentVertices();
        for (auto v : adj_vertices) {
            size_t v_idx = v.getIndex();
             if (subset.vertices.find(v_idx) != subset.vertices.end()) {
                continue;
            }
            else {
                return false;
            } 
        }
    }
    // then edges
    for (size_t e_idx : subset.edges) {
        Edge e = mesh->edge(e_idx);
        auto adj_vertices = e.adjacentVertices();
        for (auto v : adj_vertices) {
            size_t v_idx = v.getIndex();
             if (subset.vertices.find(v_idx) != subset.vertices.end()) {
                continue;
            }
            else {
                return false;
            } 
        } 
    }
    return true;
}

/*
 * Check if the given subset S is a pure simplicial complex. If so, return the degree of the complex. Otherwise, return
 * -1.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: int representing the degree of the given complex (-1 if not pure)
 */
int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {

    // TODO
    if (subset.faces.size() > 0) {
        size_t face_edges_count = 0;
        size_t face_vertices_count = 0;
        for (size_t f_idx : subset.faces) {
            Face f = mesh->face(f_idx); 
            auto adj_edges = f.adjacentEdges();
            for (auto e : adj_edges) {
                size_t e_idx = e.getIndex();
                if (subset.edges.find(e_idx) != subset.edges.end()) {
                    face_edges_count ++;
                }
                else {
                    return -1;
                }
            }
            auto adj_vertices = f.adjacentVertices();
            for (auto v : adj_vertices) {
                size_t v_idx = v.getIndex();
                if (subset.vertices.find(v_idx) != subset.vertices.end()) {
                    face_vertices_count ++;
                }
                else {
                    return -1;
                } 
            }
        }
        if (face_edges_count != subset.edges.size()) return -1;
        if (face_vertices_count != subset.vertices.size()) return -1;
        return 3;
    }
    else if (subset.edges.size() > 0) {
        size_t edge_vertices_count = 0;
        for (size_t e_idx : subset.edges) {
            Edge e = mesh->edge(e_idx);
            auto adj_vertices = e.adjacentVertices();
            for (auto v : adj_vertices) {
                size_t v_idx = v.getIndex();
                if (subset.vertices.find(v_idx) != subset.vertices.end()) {
                    edge_vertices_count ++;
                }
                else {
                    return -1;
                } 
            }
        }
        if (edge_vertices_count != subset.vertices.size()) return -1;
        return 2;
    }
    else if (subset.vertices.size() > 0) return 1;
    return 0;
}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The boundary of the given subset.
 */
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {

    // TODO
    return subset; // placeholder
}