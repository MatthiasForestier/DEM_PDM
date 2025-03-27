#include <Eigen/Dense>
#include <Eigen/CholmodSupport>

#include "CRLHelper/Optimization.h"
#include "CRLHelper/CRLTimer.h"

#include <iostream>
#include <cstdio>
#include <unistd.h>
#include <fstream>
#include <fcntl.h>

#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

typedef Eigen::CholmodSupernodalLLT<SparseMatrixF, Eigen::Upper> CholmodSolver;

int suppress_stdout() {
    fflush(stdout);

    int ret = dup(1);
    int nullfd = open("/dev/null", O_WRONLY);
    // check nullfd for error omitted
    dup2(nullfd, 1);
    close(nullfd);

    return ret;
}

void resume_stdout(int fd) {
    fflush(stdout);
    dup2(fd, 1);
    close(fd);
}

Optimization::OptimizationStatus Optimization::step(VectorXF &y) {
    Optimization::OptimizationStatus status;
    switch (optimizer) {
        case Optimization::GRADIENT_DESCENT:
            status = stepGradientDescent(y);
            break;
        case Optimization::NEWTON:
            status = stepNewton(y, false);
            break;
        case Optimization::NEWTON_WOODBURY:
            status = stepNewton(y, true);
            break;
        case Optimization::BFGS:
            status = stepBFGS(y);
            break;
        case Optimization::LBFGS:
            status = stepLBFGS(y);
            break;
        case Optimization::NLOPT:
            status = stepNLopt(y);
            break;
        default:
            assert(0);
    }
    return status;
}

Optimization::OptimizationStatus Optimization::stepGradientDescent(VectorXF &y) {
    VectorXF gradient;
    F initial_objective_value;
    if (!gradient_function(y, initial_objective_value, gradient)) {
        assert(0);
        return FAILURE;
    }

    if (gradient.squaredNorm() / y.rows() < pow(10.0, tolerance_exponent)) return CONVERGED;

    bool line_search_success = lineSearch(y, -gradient, initial_objective_value);
    return line_search_success ? SUCCESS : FAILURE;
}

Optimization::OptimizationStatus Optimization::stepBFGS(VectorXF &y) {
    if (reset_optimization) {
        int ny = y.rows();
        bfgs_inverse_hessian = MatrixXF::Zero(ny, ny);
        reset_optimization = false;
    }

    VectorXF y0 = y;

    F initial_objective_value;
    VectorXF g0;
    if (!gradient_function(y, initial_objective_value, g0)) {
        assert(0);
        return FAILURE;
    }

    if (g0.squaredNorm() / y.rows() < pow(10.0, tolerance_exponent)) return CONVERGED;

    /// Compute search direction using current inverse Hessian approximation and update y using line search.
    VectorXF search_direction = -bfgs_inverse_hessian * g0;
    bool line_search_success = lineSearch(y, search_direction, initial_objective_value);
    if (!line_search_success) {
        return FAILURE;
    }

    /// Update inverse Hessian approximation for next step. This is somewhat inefficient due to the extra gradient
    /// function call, but avoids storing previous y and g like the alternative (updating inverse hessian pre-step)
    /// would require.
    VectorXF g;
    if (!gradient_function(y, initial_objective_value, g)) {  /// initial_objective_value not used here.
        assert(0);
        return FAILURE;
    }

    VectorXF deltaY = y - y0;
    VectorXF deltaG = g - g0;

    if (deltaG.norm() > 1e-10 && deltaY.norm() > 1e-10) {
        MatrixXF I = MatrixXF::Identity(y.rows(), y.rows());
        MatrixXF A = (I - (deltaY * deltaG.transpose()) / (deltaG.transpose() * deltaY));
        MatrixXF B = (I - (deltaG * deltaY.transpose()) / (deltaG.transpose() * deltaY));
        MatrixXF C = (deltaY * deltaY.transpose()) / (deltaG.transpose() * deltaY);
        bfgs_inverse_hessian = A * bfgs_inverse_hessian * B + C;
    }
    return SUCCESS;
}

Optimization::OptimizationStatus Optimization::stepLBFGS(VectorXF &y) {
    if (reset_optimization) {
        int ny = y.rows();
        lbfgs_mat.reset(ny, lbfgs_m);
        reset_optimization = false;
    }

    VectorXF y0 = y;

    F initial_objective_value;
    VectorXF g0;
    if (!gradient_function(y, initial_objective_value, g0)) {
        assert(0);
        return FAILURE;
    }

    if (g0.squaredNorm() / y.rows() < pow(10.0, tolerance_exponent)) return CONVERGED;

    /// Compute search direction using current inverse Hessian approximation and update y using line search.
    VectorXF search_direction;
    lbfgs_mat.apply_Hv(g0, -1.0, search_direction);
    bool line_search_success = lineSearch(y, search_direction, initial_objective_value);
    if (!line_search_success) {
        return FAILURE;
    }

    /// Update inverse Hessian approximation for next step. This is somewhat inefficient due to the extra gradient
    /// function call, but avoids storing previous y and g like the alternative (updating inverse hessian pre-step)
    /// would require.
    VectorXF g;
    if (!gradient_function(y, initial_objective_value, g)) {  /// initial_objective_value not used here.
        assert(0);
        return FAILURE;
    }

    VectorXF deltaY = y - y0;
    VectorXF deltaG = g - g0;

    if (deltaG.norm() > 1e-10 && deltaY.norm() > 1e-10) {
        lbfgs_mat.add_correction(deltaY, deltaG);
    }
    return SUCCESS;
}

Optimization::OptimizationStatus Optimization::stepNewton(VectorXF &y, bool use_woodbury) {
    F initial_objective_value;
    VectorXF gradient;
    HessianF hessian;

    if (!hessian_function(y, initial_objective_value, gradient, hessian)) {
        assert(0);
        return FAILURE;
    }

    // std::cout << "Newton " << gradient.squaredNorm() / y.rows() << " " << initial_objective_value << std::endl;

    if (gradient.squaredNorm() / y.rows() < pow(10.0, tolerance_exponent)) return CONVERGED;

    /// Linear solve to obtain search direction.
    bool linear_solve_success;
    VectorXF dy;
    if (use_woodbury) {
        linear_solve_success = linearSolveWoodbury(hessian, -gradient, dy);
    } else {
        linear_solve_success = linearSolve(hessian, -gradient, dy);
    }
    if (!linear_solve_success) return FAILURE;

    bool line_search_success = lineSearch(y, dy, initial_objective_value);
    return line_search_success ? SUCCESS : FAILURE;
}

double nlopt_objective_function(const std::vector<double> &y_std, std::vector<double> &grad_std, void *f_data) {
    Optimization *optimization = reinterpret_cast<Optimization *>(f_data);
    VectorXF y = Eigen::Map<const VectorXF>(y_std.data(), y_std.size());

    F objective;
    VectorXF grad;
    optimization->gradient_function(y, objective, grad);

    grad_std.resize(grad.rows());
    for (int i = 0; i < grad.rows(); i++) {
        grad_std[i] = grad(i);
    }

    return objective;
}

Optimization::OptimizationStatus Optimization::stepNLopt(VectorXF &y) {
    if (reset_optimization) {
        int ny = y.rows();
        nlopt_opt = nlopt::opt(nlopt::algorithm(nlopt_algorithm), ny);
        nlopt_opt.set_min_objective(nlopt_objective_function, this);
        nlopt_opt.set_maxeval(100);
        reset_optimization = false;
    }

    /// Just for convergence check.
    F initial_objective_value;
    VectorXF gradient;
    if (!gradient_function(y, initial_objective_value, gradient)) {
        assert(0);
        return FAILURE;
    }
    if (gradient.squaredNorm() / y.rows() < pow(10.0, tolerance_exponent)) return CONVERGED;

    std::vector<double> y_std(y.rows());
    for (int i = 0; i < y.rows(); i++) {
        y_std[i] = y(i);
    }
    y_std = nlopt_opt.optimize(y_std);

    y = Eigen::Map<const VectorXF>(y_std.data(), y_std.size());
    return SUCCESS;
}

bool Optimization::lineSearch(VectorXF &y, const VectorXF &dy, F initial_objective_value) {
    F alpha = 1.0;
    for (int i = 0; i < 50; i++) {
        VectorXF y_line_search = y + alpha * dy;

        F new_objective_value;
        if (objective_function(y_line_search, new_objective_value)) {
            if (new_objective_value < initial_objective_value) {
                y = y_line_search;
                return true;
            }
        } else {
            /// Model generation failed. Go to next loop iteration.
        }
        alpha *= 0.5;
    }
    return false;
}

bool Optimization::lineSearch(VectorXF &y, const VectorXF &dy) {
    F initial_objective_value;
    objective_function(y, initial_objective_value);
    return lineSearch(y, dy, initial_objective_value);
}

bool Optimization::linearSolve(const HessianF &hessian, const VectorXF &b, VectorXF &x) {
    CholmodSolver solver;

    SparseMatrixF A = hessian.A;
    if (hessian.U.cols() > 0) {
        A += (hessian.U * hessian.V.transpose()).sparseView();
    }

    /// Start with very small diagonal regularizer to ensure diagonal is included in sparsity pattern.
    SparseMatrixF H(A.rows(), A.cols());
    H.setIdentity();
    H.diagonal().array() = 1e-10;
    solver.analyzePattern(A + H);

    int indefinite_count_reg_cnt = 0, invalid_search_dir_cnt = 0, invalid_residual_cnt = 0;
    F alpha = 1e-6;
    for (int i = 0; i < 50; i++) {
        int fd = suppress_stdout();
        solver.factorize(A + H);
        resume_stdout(fd);

        /// Should have NumericalIssue if K is not positive definite.
        if (solver.info() == Eigen::NumericalIssue) {
            H.diagonal().array() = alpha;
            alpha *= 10;
            indefinite_count_reg_cnt++;
            continue;
        }

        x = solver.solve(b);

        /// Check that step has positive dot product with residual.
        F dot_x_b = x.normalized().dot(b.normalized());
        // F dot_x_b = x.dot(b);
        bool search_dir_correct_sign = (dot_x_b > 0); //> 1e-6;
        if (!search_dir_correct_sign) {
            invalid_search_dir_cnt++;
        }

        /// Check for reasonable step size.
        bool solve_success = x.norm() < 1e3;
        if (!solve_success) {
            invalid_residual_cnt++;
        }

        /// In case of problem, increase regularizer by an order of magnitude and try again.
        if (search_dir_correct_sign && solve_success) {
            // std::cout << i << " regularization steps" << std::endl;
            return true;
        } else {
            H.diagonal().array() = alpha;
            alpha *= 10;
        }
    }
    return false;
}

bool Optimization::linearSolveWoodbury(const HessianF &hessian, const VectorXF &b, VectorXF &x) {
    CholmodSolver solver;

    SparseMatrixF A = hessian.A;
    MatrixXF U = hessian.U;
    MatrixXF VT = hessian.V.transpose();

    /// Add an initial very small diagonal regularizer?
    SparseMatrixF H(hessian.A.rows(), hessian.A.cols());
    H.setIdentity();
    H.diagonal().array() = 1e-10;
    solver.analyzePattern(hessian.A + H);

    int indefinite_count_reg_cnt = 0, invalid_search_dir_cnt = 0, invalid_residual_cnt = 0;
    F alpha = 1e-6;
    for (int i = 0; i < 50; i++) {
        int fd = suppress_stdout();
        solver.factorize(A + H);
        resume_stdout(fd);

        if (solver.info() == Eigen::NumericalIssue) {
            H.diagonal().array() = alpha;
            alpha *= 10;
            indefinite_count_reg_cnt++;
            continue;
        }

        if (U.cols() == 1) {
            assert((U - hessian.V).squaredNorm() < 1e-10);

            /// Sherman-Morrison Formula
            VectorXF u = U.col(0);
            MatrixXF rhs(A.rows(), 2);
            rhs.col(0) = b;
            rhs.col(1) = u;
            MatrixXF A_inv_bu = solver.solve(rhs);

            F dem = 1.0 + u.dot(A_inv_bu.col(1));

            x.noalias() = A_inv_bu.col(0) - (A_inv_bu.col(0).dot(u)) * A_inv_bu.col(1) / dem;
        } else {
            /// Woodbury Matrix Identity
            VectorXF A_inv_b = solver.solve(b);

            MatrixXF A_inv_U(U.rows(), U.cols());
            A_inv_U = solver.solve(U);

            MatrixXF C(U.cols(), U.cols());
            C.setIdentity();
            C += VT * A_inv_U;

            x.noalias() = A_inv_b - A_inv_U * C.inverse() * VT * A_inv_b;
        }

        /// Check that step has positive dot product with residual.
        F dot_x_b = x.normalized().dot(b.normalized());
        bool search_dir_correct_sign = dot_x_b > 0;
        if (!search_dir_correct_sign) {
            invalid_search_dir_cnt++;
        }

        /// Check for reasonable step size.
        bool solve_success = x.norm() < 1e3;
        if (!solve_success) {
            invalid_residual_cnt++;
        }

        /// In case of problem, increase regularizer by an order of magnitude and try again.
        if (search_dir_correct_sign && solve_success) {
            return true;
        } else {
            H.diagonal().array() = alpha;
            alpha *= 10;
        }
    }
    return false;
}

bool Optimization::checkGradient(const VectorXF &y, VectorXF &error,
                                 const std::function<bool(const VectorXF &, F &)> &func,
                                 const std::function<bool(const VectorXF &, VectorXF &)> &grad_func, F epsilon,
                                 int print_level) {
    int num_variables = (int)y.rows();

    VectorXF grad;
    if (!(grad_func(y, grad))) {
        /// Gradient evaluation failed.
        error = VectorXF::Constant(num_variables, NAN);
        if (print_level > 0)
            std::cout << "Optimization::checkGradient: Gradient evaluation failed at given state!" << std::endl;
        return false;
    }

    VectorXF grad_fd(num_variables);
    for (int i = 0; i < num_variables; ++i) {
        VectorXF y_plus = y, y_minus = y, y_plus2 = y, y_minus2 = y;
        y_plus(i) += epsilon;
        y_minus(i) -= epsilon;
        y_plus2(i) += 2.0 * epsilon;
        y_minus2(i) -= 2.0 * epsilon;

        F f_plus, f_minus, f_plus2, f_minus2;
        if (func(y_plus, f_plus) && func(y_minus, f_minus) && func(y_plus2, f_plus2) && func(y_minus2, f_minus2)) {
            F grad_fd_i = (f_plus - f_minus) / (2.0 * epsilon);
            F error_ratio =
                ((f_plus2 - f_minus2) / (4.0 * epsilon) - grad(i)) / ((f_plus - f_minus) / (2.0 * epsilon) - grad(i));

            if ((print_level == 1 && fabs(grad_fd_i - grad(i)) > epsilon) || print_level > 1)
                std::cout << "Optimization::checkGradient: grad_fd[" << i << "] = " << grad_fd_i << ", grad[" << i
                          << "] = " << grad(i) << ", error ratio = " << error_ratio << std::endl;

            grad_fd(i) = grad_fd_i;
        } else {
            if (print_level > 0) std::cout << "Optimization::checkGradient: Function evaluation failed!" << std::endl;
            grad_fd(i) = NAN;
        }
    }

    error = grad - grad_fd;
    return error.allFinite() && error.cwiseAbs().maxCoeff() < epsilon;
}

bool Optimization::checkHessian(const VectorXF &y, MatrixXF &error,
                                const std::function<bool(const VectorXF &, VectorXF &)> &grad_func,
                                const std::function<bool(const VectorXF &, MatrixXF &)> &hess_func, F epsilon,
                                int print_level) {
    int num_variables = (int)y.rows();

    MatrixXF hess;
    if (!(hess_func(y, hess))) {
        /// Hessian function failed.
        if (print_level > 0)
            std::cout << "Optimization::checkHessian: Hessian evaluation failed at given state!" << std::endl;
        error = MatrixXF::Constant(num_variables, num_variables, NAN);
        return false;
    }

    MatrixXF hess_fd(num_variables, num_variables);
    for (int i = 0; i < num_variables; ++i) {
        VectorXF y_plus = y, y_minus = y, y_plus2 = y, y_minus2 = y;
        y_plus(i) += epsilon;
        y_minus(i) -= epsilon;
        y_plus2(i) += 2.0 * epsilon;
        y_minus2(i) -= 2.0 * epsilon;

        VectorXF grad_plus, grad_minus, grad_plus2, grad_minus2;
        if (grad_func(y_plus, grad_plus) && grad_func(y_minus, grad_minus) && grad_func(y_plus2, grad_plus2) &&
            grad_func(y_minus2, grad_minus2)) {
            for (int j = i; j < num_variables; ++j) {  /// Upper triangular hessian
                F hess_fd_i_j = (grad_plus(j) - grad_minus(j)) / (2.0 * epsilon);
                F error_ratio = ((grad_plus2(j) - grad_minus2(j)) / (4.0 * epsilon) - hess(i, j)) /
                                ((grad_plus(j) - grad_minus(j)) / (2.0 * epsilon) - hess(i, j));

                if ((print_level == 1 && fabs(hess_fd_i_j - hess(i, j)) > epsilon) || print_level > 1)
                    std::cout << "Optimization::checkHessian: hess_fd[" << i << ", " << j << "] = " << hess_fd_i_j
                              << ", hess[" << i << ", " << j << "] = " << hess(i, j)
                              << ", error ratio = " << error_ratio << std::endl;

                hess_fd(i, j) = hess_fd_i_j;
            }
        } else {
            if (print_level > 0) std::cout << "Optimization::checkHessian: Gradient evaluation failed!" << std::endl;
            hess_fd.row(i).setConstant(NAN);
        }
    }

    error = hess - hess_fd;
    return error.allFinite() && error.cwiseAbs().maxCoeff() < epsilon;
}

void Optimization::makeConfigMenu() {
    std::vector<std::string> optimizer_names = getOptimizerNames();
    if (ImGui::Combo("Optimizer", reinterpret_cast<int *>(&optimizer), optimizer_names)) {
        reset_optimization = true;
    }
    if (optimizer == NLOPT) {
        std::vector<std::string> nlopt_names;
        for (int i = 0; i < nlopt::NUM_ALGORITHMS; i++) {
            nlopt_names.emplace_back(nlopt::algorithm_name((nlopt::algorithm)i));
        }
        if (ImGui::Combo("Algorithm##NLOptOptimizer", reinterpret_cast<int *>(&nlopt_algorithm), nlopt_names)) {
            reset_optimization = true;
        }
    }
    ImGui::Text("Tolerance: 10^");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(100);
    ImGui::InputInt("##ToleranceForwardSimConvergence", &tolerance_exponent);
}
