#include <igl/triangle/triangulate.h>

#include "CRLHelper/HeatMethod.h"
#include "CRLHelper/GeometryHelper.h"

#include <set>

#include <Eigen/SparseLU>
typedef Eigen::SparseLU<SparseMatrixF> SparseLU;

#include <iostream>

bool HeatMethod::solveHeat(const MatrixXF &V, const MatrixXI &T, const VectorXI &B, VectorXF &sol) {
    TripletListF triplets_A;

    /// Non-boundary vertices are inverse-distance-weighted average of neighbors.
    for (int i = 0; i < T.rows(); i++) {
        for (int j = 0; j < 3; j++) {
            int iv0 = T(i, j);
            int iv1 = T(i, (j + 1) % 3);
            int iv2 = T(i, (j + 2) % 3);

            Vector2F v0 = V.row(iv0);
            Vector2F v1 = V.row(iv1);
            Vector2F v2 = V.row(iv2);

            Vector2F v20 = v0 - v2;
            Vector2F v21 = v1 - v2;

            F dot = v20.dot(v21);
            F cross = v20(0) * v21(1) - v20(1) * v21(0);
            F cot = dot / cross;

            if (B(iv0) == 0) {
                triplets_A.emplace_back(iv0, iv0, cot);
                triplets_A.emplace_back(iv0, iv1, -cot);
            }
            if (B(iv1) == 0) {
                triplets_A.emplace_back(iv1, iv1, cot);
                triplets_A.emplace_back(iv1, iv0, -cot);
            }
        }
    }

    /// Boundary vertices are fixed.
    for (int iv = 0; iv < V.rows(); iv++) {
        if (B(iv) != 0) {
            triplets_A.emplace_back(iv, iv, 1);
        }
    }

    SparseMatrixF A(V.rows(), V.rows());
    A.setFromTriplets(triplets_A.begin(), triplets_A.end());

    /// RHS is zero except for boundary vertex rows (i.e. boundary conditions).
    VectorXF b = VectorXF::Zero(V.rows());
    Vector4F b_vals = {0, 0, 1, 10};  // Domain, Target, Farfield, Inside.
    for (int iv = 0; iv < V.rows(); iv++) {
        b(iv) = b_vals(B(iv));
    }

    /// Solve.
    SparseLU solver;
    solver.compute(A);
    sol = solver.solve(b);

    return true;
}

struct edgeCompare {
    bool operator()(const Vector2I &a, const Vector2I &b) const {
        for (int i = 0; i < std::min(a.size(), b.size()); ++i) {
            if (a[i] < b[i]) return true;
            if (a[i] > b[i]) return false;
        }
        return a.size() < b.size();
    }
};

bool HeatMethod::mesh2D(const MatrixXF &V_in, const MatrixXI &E_in, MatrixXF &V_out, MatrixXI &T_out, VectorXI &B_out,
                        F far_field_radius, F max_area) {
    // No H points.
    MatrixXF H(0, 2);

    // Add far field boundary.
    int N_ff = 20;
    MatrixXF V_in_with_ff(V_in.rows() + N_ff, 2);
    V_in_with_ff.topRows(V_in.rows()) = V_in;
    MatrixXI E_in_with_ff(E_in.rows() + N_ff, 2);
    E_in_with_ff.topRows(E_in.rows()) = E_in;
    for (int i = 0; i < N_ff; i++) {
        F angle = i * 2 * M_PI / N_ff;
        V_in_with_ff.row(V_in.rows() + i) = Vector2F(cos(angle), sin(angle)) * far_field_radius;
        E_in_with_ff.row(E_in.rows() + i) = Vector2I(i, (i + 1) % N_ff).array() + V_in.rows();
    }

    // Label shape and far field edges.
    VectorXI VM, EM, VM_out, EM_out;
    MatrixXI E_out;
    EM.resize(E_in_with_ff.rows());
    EM.head(E_in.rows()).setConstant(1);
    EM.tail(N_ff).setConstant(2);

    // Triangulate with Triangle using libigl
    std::string flags = "pqa" + std::to_string(max_area);  // 'a' followed by the max area constraint
    // std::string flags = "pa0.2";
    igl::triangle::triangulate(V_in_with_ff, E_in_with_ff, H, VM, EM, flags, V_out, T_out, VM_out, E_out, EM_out);

    // Mark vertices
    B_out = VectorXI::Zero(V_out.rows());
    for (int i = 0; i < E_out.rows(); i++) {
        int i0 = E_out(i, 0);
        int i1 = E_out(i, 1);
        int label = EM_out(i);
        if (label != 0) {
            B_out(i0) = label;
            B_out(i1) = label;
        }
    }
    for (int i = 0; i < V_out.rows(); i++) {
        if (B_out(i) != 0) continue;

        Vector2F v = V_out.row(i);
        F winding = 0;
        for (int j = 0; j < E_in.rows(); j++) {
            winding += GeometryHelper::windingNumber2D(v, V_in.row(E_in(j, 0)), V_in.row(E_in(j, 1)));
        }
        if (fabs(winding) > M_PI) {
            B_out(i) = 3;
        }
    }

    return true;
}
