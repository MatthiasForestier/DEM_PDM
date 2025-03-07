#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <memory>

template <typename T, int dim>
using Vector = Eigen::Matrix<T, dim, 1, 0, dim, 1>;

template <typename T, int n, int m>
using Matrix = Eigen::Matrix<T, n, m, 0, n, m>;

using F = double;
using I = int;

using Vector2F = Vector<F, 2>;
using Vector3F = Vector<F, 3>;
using Vector4F = Vector<F, 4>;
using Vector6F = Vector<F, 6>;
using Vector2I = Vector<I, 2>;
using Vector3I = Vector<I, 3>;
using Vector4I = Vector<I, 4>;
using Matrix2F = Matrix<F, 2, 2>;
using Matrix3F = Matrix<F, 3, 3>;

using VectorXF = Matrix<F, Eigen::Dynamic, 1>;
using MatrixXF = Matrix<F, Eigen::Dynamic, Eigen::Dynamic>;
using VectorXI = Vector<I, Eigen::Dynamic>;
using MatrixXI = Matrix<I, Eigen::Dynamic, Eigen::Dynamic>;

using TripletF = Eigen::Triplet<F>;
using TripletListF = std::vector<TripletF>;
using SparseVectorF = Eigen::SparseVector<F>;
using SparseMatrixF = Eigen::SparseMatrix<F>;
using TripletI = Eigen::Triplet<I>;
using TripletListI = std::vector<TripletI>;
using SparseVectorI = Eigen::SparseVector<I>;
using SparseMatrixI = Eigen::SparseMatrix<I>;
