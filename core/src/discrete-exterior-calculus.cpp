// PLEASE READ:
//
// This file additional geometry routines for the VertexPositionGeometry class in Geometry Central. Because we are
// "inside" the class, we no longer have to call
//
//          geometry->inputVertexPositions[v], etc.
//
// We can just call
//
//          this->inputVertexPositions[v], etc.
//
// or simply
//
//          inputVertexPositions[v], etc.
//
// In addition, we no longer access the corresponding surface mesh via
//
//          mesh->vertices(), etc.
//
// but instead <mesh> is not a pointer anymore, so we use
//
//          mesh.vertices(), etc.
//
// Functions in this file can be called from other projects simply by using geometry->buildHodgeStar0Form(), etc. where
// "geometry" is a pointer to a VertexPositionGeometry. This avoids having to declare a GeometryRoutines object in every
// project, and also mimics the way that geometry routines are normally called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.

#include "geometrycentral/surface/vertex_position_geometry.h"

namespace geometrycentral {
namespace surface {


/*
 * Build Hodge operator on 0-forms.
 * By convention, the area of a vertex is 1.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar0Form() const {
    int n = mesh.nVertices();
    SparseMatrix<double> HodgeStar0Form(n, n);
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(n);
    for (Vertex v : mesh.vertices()) {
        double area = barycentricDualArea(v);
        int i = v.getIndex();
        tripletList.emplace_back(i, i, area);
    }
    HodgeStar0Form.setFromTriplets(tripletList.begin(), tripletList.end());
    return HodgeStar0Form;
}

/*
 * Build Hodge operator on 1-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar1Form() const {
    int n = mesh.nEdges();
    SparseMatrix<double> HodgeStar1Form(n, n);
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(n);
    for (Edge e : mesh.edges()) {
        Halfedge he = e.halfedge();
        double cota = cotan(he);
        double cotb = cotan(he.twin());
        int i = e.getIndex();
        tripletList.emplace_back(i, i, 0.5 * (cota + cotb));
    }
    HodgeStar1Form.setFromTriplets(tripletList.begin(), tripletList.end());
    return HodgeStar1Form;
}

/*
 * Build Hodge operator on 2-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 2-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar2Form() const {
    int n = mesh.nFaces();
    SparseMatrix<double> HodgeStar2Form(n, n);
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(n);
    for (Face f : mesh.faces()) {
        double area = faceArea(f);
        int i = f.getIndex();
        tripletList.emplace_back(i, i, 1.0/area);
    }
    HodgeStar2Form.setFromTriplets(tripletList.begin(), tripletList.end());
    return HodgeStar2Form;
}

/*
 * Build exterior derivative on 0-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative0Form() const {
    size_t nv = mesh.nVertices();
    size_t ne = mesh.nEdges(); 
    SparseMatrix<double> Derivative0Form(ne, nv);
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(2*ne);

    for (Edge e : mesh.edges()) {
        Halfedge he = e.halfedge();
        Vertex vs = he.vertex();
        Vertex vt = he.twin().vertex();
        tripletList.emplace_back(e.getIndex(), vs.getIndex(), -1.0);
        tripletList.emplace_back(e.getIndex(), vt.getIndex(), 1.0);
    }

    Derivative0Form.setFromTriplets(tripletList.begin(), tripletList.end());

    return Derivative0Form;
}

/*
 * Build exterior derivative on 1-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative1Form() const {
    size_t ne = mesh.nEdges();
    size_t nf = mesh.nFaces();
    SparseMatrix<double> Derivative1Form(nf, ne);
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(3*nf);

    for (Face f : mesh.faces()) {
        Halfedge he = f.halfedge();
        for (int i = 0; i < 3; ++i) {
            Edge e = he.edge();
            bool isSameOrientation = (he == e.halfedge());
            tripletList.emplace_back(f.getIndex(), e.getIndex(), isSameOrientation ? 1.0 : -1.0);
            he = he.next();
        }
    }

    Derivative1Form.setFromTriplets(tripletList.begin(), tripletList.end());

    return Derivative1Form;
}

} // namespace surface
} // namespace geometrycentral