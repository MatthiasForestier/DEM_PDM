#include "Algebra.h"
#include <Eigen/Eigenvalues>
#include <limits>
#include <cmath>

float computeConditionNumber(const MatrixXF &denseMatrix) {
    Eigen::SelfAdjointEigenSolver<MatrixXF> eigensolver(denseMatrix);
    if (eigensolver.info() != Eigen::Success) {
        // In case of failure, return a large number.
        return std::numeric_limits<float>::max();
    }
    float minEigen = eigensolver.eigenvalues().minCoeff();
    float maxEigen = eigensolver.eigenvalues().maxCoeff();

    // Avoid division by zero.
    if (std::abs(minEigen) < 1e-8f)
        return std::numeric_limits<float>::max();
    return maxEigen / minEigen;
}
