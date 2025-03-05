#include "CRLHelper/SensitivityAnalysis.h"

#include <Eigen/SparseLU>
#include <Eigen/LU>
#include <Eigen/PardisoSupport>
#include <Eigen/SparseQR>
#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <Eigen/CholmodSupport>
#include <Eigen/UmfPackSupport>
#include <Eigen/KLUSupport>

#include <iostream>

#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

#include "ThirdParty/mma/src/mma/MMASolver.h"
#include "ThirdParty/mma/src/gcmma/GCMMASolver.h"

#include "CRLHelper/CRLTimer.h"

typedef Eigen::SparseLU<SparseMatrixF> SparseLU;
typedef Eigen::SparseQR<SparseMatrixF, Eigen::COLAMDOrdering<int>> SparseQR;
typedef Eigen::PardisoLU<SparseMatrixF> PardisoLU;
typedef Eigen::UmfPackLU<SparseMatrixF> UmfPackLU;
typedef Eigen::KLU<SparseMatrixF> KLU;
typedef Eigen::LLT<MatrixXF> DenseLLT;
typedef Eigen::CholmodSupernodalLLT<SparseMatrixF, Eigen::Upper> CholmodSolver;

Optimization::OptimizationStatus SensitivityAnalysis::step() {
    if (partial_steps) {
        return partialStep();
    } else {
        return fullStep();
    }
}

Optimization::OptimizationStatus SensitivityAnalysis::fullStep() {
    VectorXF u;
    getUFunction(u);

    Optimization::OptimizationStatus status = runEnergyMinimization();
    if (status == Optimization::FAILURE) {
        return Optimization::FAILURE;
    }

    VectorXF du;
    F initial_objective_value;
    bool use_line_search = true;
    switch (optimizer) {
        case SensitivityAnalysis::SA_GRADIENT_DESCENT:
            status = getDirectionGradientDescent(u, du, initial_objective_value);
            break;
        case SensitivityAnalysis::SA_BFGS:
            status = getDirectionBFGS(u, du, initial_objective_value);
            break;
        case SensitivityAnalysis::SA_LBFGS:
            status = getDirectionLBFGS(u, du, initial_objective_value);
            break;
        case SensitivityAnalysis::SA_GAUSS_NEWTON:
            status = getDirectionGaussNewton(u, du, initial_objective_value);
            break;
        case SensitivityAnalysis::SA_SGN:
            status = getDirectionSGN(u, du, initial_objective_value);
            break;
        case SensitivityAnalysis::SA_GCMMA:
            status = stepGCMMA(u);
            use_line_search = false;  /// No line search.
            break;
        case SensitivityAnalysis::SA_NLOPT:
            status = stepNLopt(u);
            use_line_search = false;  /// No line search.
            break;
        default:
            assert(0);
    }

    if (status == Optimization::SUCCESS && use_line_search) {
        bool line_search_success = lineSearch(u, du, initial_objective_value);
        status = line_search_success ? Optimization::SUCCESS : Optimization::FAILURE;
    }

    if (status == Optimization::CONVERGED) {
        std::cout << "line search converged" << std::endl;
    }

    return status;
}

Optimization::OptimizationStatus SensitivityAnalysis::partialStep() {
    if (partial_step_state == PARTIAL_STEP_START) {
        partial_step_state = PARTIAL_STEP_INITIAL_EQUILIBRIUM;
    }

    if (partial_step_state == PARTIAL_STEP_INITIAL_EQUILIBRIUM) {
        Optimization::OptimizationStatus status = energyMinimizationStepFunction();
        switch (status) {
            case Optimization::FAILURE:
                return Optimization::FAILURE;
                break;
            case Optimization::CONVERGED:
                getYFunction(partial_step_y0);
                partial_step_state = PARTIAL_STEP_SEARCH_DIRECTION;
                break;
            case Optimization::SUCCESS:
                break;
            default:
                assert(0);
        }
    }

    if (partial_step_state == PARTIAL_STEP_SEARCH_DIRECTION) {
        Optimization::OptimizationStatus status;
        partial_step_state = PARTIAL_STEP_LINE_SEARCH_OUTER;
        partial_step_line_search_step = -1;

        VectorXF u;
        getUFunction(u);
        partial_step_u0 = u;

        switch (optimizer) {
            case SensitivityAnalysis::SA_GRADIENT_DESCENT:
                status = getDirectionGradientDescent(u, partial_step_search_direction, partial_step_initial_objective);
                break;
            case SensitivityAnalysis::SA_BFGS:
                status = getDirectionBFGS(u, partial_step_search_direction, partial_step_initial_objective);
                break;
            case SensitivityAnalysis::SA_LBFGS:
                status = getDirectionLBFGS(u, partial_step_search_direction, partial_step_initial_objective);
                break;
            case SensitivityAnalysis::SA_GAUSS_NEWTON:
                status = getDirectionGaussNewton(u, partial_step_search_direction, partial_step_initial_objective);
                break;
            case SensitivityAnalysis::SA_SGN:
                status = getDirectionSGN(u, partial_step_search_direction, partial_step_initial_objective);
                break;
            case SensitivityAnalysis::SA_GCMMA:
            case SensitivityAnalysis::SA_NLOPT:
            default:
                partial_step_state = PARTIAL_STEP_START;  /// No partial steps implemented for these optimizers.
                assert(0);
        }

        switch (status) {
            case Optimization::FAILURE:
                return Optimization::FAILURE;
                break;
            case Optimization::CONVERGED:
                return Optimization::CONVERGED;
                break;
            default:
                break;
        }
    }

    if (partial_step_state == PARTIAL_STEP_LINE_SEARCH_OUTER) {
        partial_step_line_search_step++;
        if (partial_step_line_search_step >= 50) {
            partial_step_state = PARTIAL_STEP_START;
        } else {
            F alpha = pow(0.5, partial_step_line_search_step);
            VectorXF u_line_search = partial_step_u0 + alpha * partial_step_search_direction;
            if (applyUFunction(u_line_search)) {
                applyYFunction(partial_step_y0);
                partial_step_state = PARTIAL_STEP_LINE_SEARCH_EQUILIBRIUM;
            }  /// Else take another line search step.
        }
    }

    if (partial_step_state == PARTIAL_STEP_LINE_SEARCH_EQUILIBRIUM) {
        Optimization::OptimizationStatus status = energyMinimizationStepFunction();
        VectorXF y;
        F alpha = pow(0.5, partial_step_line_search_step);
        VectorXF u = partial_step_u0 + alpha * partial_step_search_direction;
        F new_objective_value;
        switch (status) {
            case Optimization::FAILURE:
                return Optimization::FAILURE;
                break;
            case Optimization::CONVERGED:
                partial_step_state = PARTIAL_STEP_LINE_SEARCH_OUTER;
                getYFunction(y);
                if (objectiveFunction(u, y, new_objective_value)) {
                    if (new_objective_value < partial_step_initial_objective) {
                        partial_step_state = PARTIAL_STEP_START;
                        return Optimization::SUCCESS;
                    }
                } else {
                    assert(0);
                    return Optimization::FAILURE;
                }
                break;
            case Optimization::SUCCESS:
                break;
            default:
                assert(0);
        }
    }

    return Optimization::PARTIAL;
}

Optimization::OptimizationStatus SensitivityAnalysis::getDirectionGradientDescent(const VectorXF &u, VectorXF &du,
                                                                                  F &initial_objective_value) {
    VectorXF y;
    getYFunction(y);
    VectorXF gradient;
    bool success = computeSAGradient(u, y, gradient);
    if (!success) {
        assert(0);
        return Optimization::FAILURE;
    }
    if (gradient.squaredNorm() / u.rows() < pow(10.0, tolerance_exponent)) return Optimization::CONVERGED;

    success = success && objectiveFunction(u, y, initial_objective_value);
    if (!success) {
        assert(0);
        return Optimization::FAILURE;
    }

    std::cout << gradient.squaredNorm() / u.rows() << " " << initial_objective_value << std::endl;

    du = -gradient.normalized();
    return Optimization::SUCCESS;
}

Optimization::OptimizationStatus SensitivityAnalysis::getDirectionBFGS(const VectorXF &u, VectorXF &du,
                                                                       F &initial_objective_value) {
    VectorXF y;
    getYFunction(y);
    VectorXF gradient;
    bool success = computeSAGradient(u, y, gradient);
    if (!success) {
        assert(0);
        return Optimization::FAILURE;
    }
    if (gradient.squaredNorm() / u.rows() < pow(10.0, tolerance_exponent)) return Optimization::CONVERGED;

    if (reset_optimization) {
        int nu = u.rows();
        bfgs_inverse_hessian = MatrixXF::Zero(nu, nu);
        reset_optimization = false;
    } else {
        VectorXF deltaU = u - prev_u;
        VectorXF deltaG = gradient - prev_g;
        if (deltaG.norm() > 1e-10 && deltaU.norm() > 1e-10) {
            MatrixXF I = MatrixXF::Identity(u.rows(), u.rows());
            MatrixXF A = (I - (deltaU * deltaG.transpose()) / (deltaG.transpose() * deltaU));
            MatrixXF B = (I - (deltaG * deltaU.transpose()) / (deltaG.transpose() * deltaU));
            MatrixXF C = (deltaU * deltaU.transpose()) / (deltaG.transpose() * deltaU);
            bfgs_inverse_hessian = A * bfgs_inverse_hessian * B + C;
        }
    }
    prev_u = u;
    prev_g = gradient;

    success = success && objectiveFunction(u, y, initial_objective_value);
    if (!success) {
        assert(0);
        return Optimization::FAILURE;
    }

    std::cout << gradient.squaredNorm() / u.rows() << " " << initial_objective_value << std::endl;

    du = -bfgs_inverse_hessian * gradient;
    if (du.dot(gradient) >= 0) {
        std::cout << "BADDDD" << std::endl;
        return Optimization::FAILURE;
    }

    return Optimization::SUCCESS;
}

Optimization::OptimizationStatus SensitivityAnalysis::getDirectionLBFGS(const VectorXF &u, VectorXF &du,
                                                                        F &initial_objective_value) {
    VectorXF y;
    getYFunction(y);
    VectorXF gradient;
    bool success = computeSAGradient(u, y, gradient);
    if (!success) {
        assert(0);
        return Optimization::FAILURE;
    }
    if (gradient.squaredNorm() / u.rows() < pow(10.0, tolerance_exponent)) return Optimization::CONVERGED;

    if (reset_optimization) {
        int nu = u.rows();
        lbfgs_mat.reset(nu, lbfgs_m);
        reset_optimization = false;
    } else {
        VectorXF deltaU = u - prev_u;
        VectorXF deltaG = gradient - prev_g;
        if (deltaG.norm() > 1e-10 && deltaU.norm() > 1e-10) {
            lbfgs_mat.add_correction(deltaU, deltaG);
        }
    }
    prev_u = u;
    prev_g = gradient;

    success = success && objectiveFunction(u, y, initial_objective_value);
    if (!success) {
        assert(0);
        return Optimization::FAILURE;
    }

    std::cout << gradient.squaredNorm() / u.rows() << " " << initial_objective_value << std::endl;

    lbfgs_mat.apply_Hv(gradient, -1.0, du);
    if (du.dot(gradient) >= 0) {
        std::cout << "BADDDD" << std::endl;
        lbfgs_mat.reset(u.rows(), lbfgs_m);
        return Optimization::FAILURE;
    }

    return Optimization::SUCCESS;
}

Optimization::OptimizationStatus SensitivityAnalysis::getDirectionGaussNewton(const VectorXF &u, VectorXF &du,
                                                                              F &initial_objective_value) {
    VectorXF y;
    getYFunction(y);
    MatrixXF dydu;
    computeSensitivityMatrix(u, y, dydu);
    VectorXF dLdu;
    bool success = computeSAGradient(u, y, dydu, dLdu);
    if (!success) {
        assert(0);
        return Optimization::FAILURE;
    }
    if (dLdu.squaredNorm() / u.rows() < pow(10.0, tolerance_exponent)) return Optimization::CONVERGED;

    success = success && objectiveFunction(u, y, initial_objective_value);
    if (!success) {
        assert(0);
        return Optimization::FAILURE;
    }

    std::cout << dLdu.squaredNorm() / u.rows() << " " << initial_objective_value << std::endl;

    SparseMatrixF p2Lpu2, p2Lpypu, p2Lpy2;
    hessianFunction(u, y, p2Lpu2, p2Lpypu, p2Lpy2);
    MatrixXF gauss_newton_hessian =
        dydu.transpose() * p2Lpy2 * dydu + p2Lpypu.transpose() * dydu + dydu.transpose() * p2Lpypu + p2Lpu2;
    MatrixXF identity = MatrixXF::Identity(u.rows(), u.rows());
    F alpha = 1e-10;

    DenseLLT solver;
    solver.compute(gauss_newton_hessian + alpha * identity);
    while (solver.info() == Eigen::NumericalIssue) {
        alpha *= 10;
        solver.compute(gauss_newton_hessian + alpha * identity);
    }

    du = solver.solve(-dLdu).normalized() * 1e-1;
    if (du.dot(dLdu) >= 0) {
        std::cout << "BADDDD" << std::endl;
        return Optimization::FAILURE;
    }

    return Optimization::SUCCESS;
}

Optimization::OptimizationStatus SensitivityAnalysis::getDirectionSGN(const VectorXF &u, VectorXF &du,
                                                                      F &initial_objective_value) {
    VectorXF y;
    getYFunction(y);
    VectorXF dLdu;
    bool success = computeSAGradient(u, y, dLdu);
    if (!success) {
        assert(0);
        return Optimization::FAILURE;
    }
    if (dLdu.squaredNorm() / u.rows() < pow(10.0, tolerance_exponent)) return Optimization::CONVERGED;

    success = success && objectiveFunction(u, y, initial_objective_value);
    if (!success) {
        assert(0);
        return Optimization::FAILURE;
    }

    std::cout << dLdu.squaredNorm() / u.rows() << " " << initial_objective_value << std::endl;

    success = success && computeSGNSearchDirection(u, y, dLdu, du);
    if (!success) {
        assert(0);
        return Optimization::FAILURE;
    }
    if (du.dot(dLdu) >= 0) {
        std::cout << "BADDDD" << std::endl;
        return Optimization::FAILURE;
    }

    return Optimization::SUCCESS;
}

double sa_nlopt_objective_function(const std::vector<double> &u_std, std::vector<double> &grad_std, void *f_data) {
    SensitivityAnalysis *sensitivity_analysis = reinterpret_cast<SensitivityAnalysis *>(f_data);
    VectorXF u = Eigen::Map<const VectorXF>(u_std.data(), u_std.size());

    bool valid = sensitivity_analysis->applyUFunction(u);
    if (!valid) {
        std::cout << "SA: invalid u" << std::endl;
        return 1e10;
    }
    Optimization::OptimizationStatus status = sensitivity_analysis->runEnergyMinimization();
    if (status == Optimization::FAILURE) {
        std::cout << "SA: energy minimization failed" << std::endl;
        return 1e10;
    }

    VectorXF y;
    sensitivity_analysis->getYFunction(y);
    F objective;
    sensitivity_analysis->objectiveFunction(u, y, objective);
    VectorXF grad;
    sensitivity_analysis->computeSAGradient(u, y, grad);

    grad_std.resize(grad.rows());
    for (int i = 0; i < grad.rows(); i++) {
        grad_std[i] = grad(i);
    }

    return objective;
}

Optimization::OptimizationStatus SensitivityAnalysis::stepGCMMA(VectorXF &u) {
    VectorXF y0;
    getYFunction(y0);
    VectorXF dLdu;
    bool success = computeSAGradient(u, y0, dLdu);
    if (!success) {
        assert(0);
        return Optimization::FAILURE;
    }
    if (dLdu.squaredNorm() / u.rows() < pow(10.0, tolerance_exponent)) return Optimization::CONVERGED;

    F initial_objective;
    success = success && objectiveFunction(u, y0, initial_objective);
    if (!success) {
        assert(0);
        return Optimization::FAILURE;
    }

    std::cout << dLdu.squaredNorm() / u.rows() << " " << initial_objective << std::endl;

    int nu = u.rows();
    if (reset_optimization) {
        gcmma_solver = std::make_shared<GCMMASolver>(nu, 0);
        reset_optimization = false;
    }

    VectorXF placeholder(1);
    VectorXF umin = VectorXF::Constant(nu, -0.01);
    VectorXF umax = VectorXF::Constant(nu, 50);

    VectorXF unew(nu);
    gcmma_solver->OuterUpdate(unew.data(), u.data(), initial_objective, dLdu.data(), placeholder.data(),
                              placeholder.data(), umin.data(), umax.data());

    F new_objective;
    success = false;
    for (int i = 0; i < 10; i++) {
        applyYFunction(y0);
        if (!applyUFunction(unew)) {
            new_objective = 1e10;
        } else {
            Optimization::OptimizationStatus status = runEnergyMinimization();
            if (status == Optimization::FAILURE) {
                new_objective = 1e10;
            } else {
                VectorXF ynew;
                getYFunction(ynew);
                objectiveFunction(unew, ynew, new_objective);
            }
        }
        if (new_objective < initial_objective) {
            u = unew;
            success = true;
            break;
        }

        gcmma_solver->InnerUpdate(unew.data(), new_objective, placeholder.data(), u.data(), initial_objective,
                                  dLdu.data(), placeholder.data(), placeholder.data(), umin.data(), umax.data());
    }

    if (success) {
        return Optimization::SUCCESS;
    } else {
        applyYFunction(y0);
        applyUFunction(u);
        return Optimization::FAILURE;
    }
}

Optimization::OptimizationStatus SensitivityAnalysis::stepNLopt(VectorXF &u) {
    if (reset_optimization) {
        int nu = u.rows();
        nlopt_opt = nlopt::opt(nlopt::algorithm(nlopt_algorithm), nu);
        nlopt_opt.set_min_objective(sa_nlopt_objective_function, this);
        nlopt_opt.set_maxeval(100);
        reset_optimization = false;
    }

    /// Just for convergence check.
    VectorXF y;
    getYFunction(y);
    VectorXF dLdu;
    bool success = computeSAGradient(u, y, dLdu);
    if (!success) {
        assert(0);
        return Optimization::FAILURE;
    }

    if (dLdu.squaredNorm() / u.rows() < pow(10.0, tolerance_exponent)) return Optimization::CONVERGED;

    F initial_objective_value;
    success = success && objectiveFunction(u, y, initial_objective_value);
    if (!success) {
        assert(0);
        return Optimization::FAILURE;
    }

    std::cout << dLdu.squaredNorm() / u.rows() << " " << initial_objective_value << std::endl;

    std::vector<double> u_std(u.rows());
    for (int i = 0; i < u.rows(); i++) {
        u_std[i] = u(i);
    }
    u_std = nlopt_opt.optimize(u_std);

    u = Eigen::Map<const VectorXF>(u_std.data(), u_std.size());
    return Optimization::SUCCESS;
}

bool SensitivityAnalysis::computeSensitivityMatrix(const VectorXF &u, const VectorXF &y, MatrixXF &dydu) {
    if (computeSensitivityMatrixOverride(u, y, dydu)) {
        return true;
    }

    SparseMatrixF d2Edy2, d2Edydu;
    energyYHessianFunction(u, y, d2Edy2);
    energyMixedHessianFunction(u, y, d2Edydu);
    MatrixXF rhs = -d2Edydu.toDense();

    CRLTimer timer;

    // SparseLU solver;
    // PardisoLU solver;
    // SparseQR solver;
    // UmfPackLU solver;
    KLU solver;

    solver.compute(d2Edy2);
    // timer.printTime("SA solve compute");
    dydu = solver.solve(rhs);
    // timer.printTime("SA solve solve");

    return true;
}

bool SensitivityAnalysis::computeSAGradient(const VectorXF &u, const VectorXF &y, const MatrixXF &dydu,
                                            VectorXF &gradient) {
    VectorXF pLpu, pLpy;
    gradientFunction(u, y, pLpu, pLpy);

    gradient = pLpy.transpose() * dydu + pLpu.transpose();
    return true;
}
bool SensitivityAnalysis::computeSAGradient(const VectorXF &u, const VectorXF &y, VectorXF &gradient) {
    MatrixXF dydu;
    computeSensitivityMatrix(u, y, dydu);

    return computeSAGradient(u, y, dydu, gradient);
}

void addTripletsBlock(TripletListF &v, const SparseMatrixF &M, int row_offset, int col_offset) {
    TripletListF temp;
    for (int i = 0; i < M.outerSize(); i++) {
        for (typename SparseMatrixF::InnerIterator it(M, i); it; ++it) {
            temp.emplace_back(it.row() + row_offset, it.col() + col_offset, it.value());
        }
    }
    v.insert(v.end(), temp.begin(), temp.end());
}
void addTripletsBlockTranspose(TripletListF &v, const SparseMatrixF &M, int row_offset, int col_offset) {
    TripletListF temp;
    for (int i = 0; i < M.outerSize(); i++) {
        for (typename SparseMatrixF::InnerIterator it(M, i); it; ++it) {
            temp.emplace_back(it.col() + row_offset, it.row() + col_offset, it.value());
        }
    }
    v.insert(v.end(), temp.begin(), temp.end());
}
void addIdentityBlock(TripletListF &v, int row_offset, int col_offset, int size, F val) {
    TripletListF temp;
    for (int i = 0; i < size; i++) {
        temp.emplace_back(i + row_offset, i + col_offset, val);
    }
    v.insert(v.end(), temp.begin(), temp.end());
}

bool SensitivityAnalysis::computeSGNSearchDirection(const VectorXF &u, const VectorXF &y, const VectorXF &dLdu,
                                                    VectorXF &search_direction) {
    VectorXF du_gauss_newton;
    {
        VectorXF y;
        getYFunction(y);
        MatrixXF dydu;
        computeSensitivityMatrix(u, y, dydu);
        VectorXF dLdu;
        bool success = computeSAGradient(u, y, dydu, dLdu);

        SparseMatrixF p2Lpu2, p2Lpypu, p2Lpy2;
        hessianFunction(u, y, p2Lpu2, p2Lpypu, p2Lpy2);
        MatrixXF gauss_newton_hessian =
            dydu.transpose() * p2Lpy2 * dydu + p2Lpypu.transpose() * dydu + dydu.transpose() * p2Lpypu + p2Lpu2;
        MatrixXF identity = MatrixXF::Identity(u.rows(), u.rows());
        F alpha = 1e-10;

        DenseLLT solver;
        solver.compute(gauss_newton_hessian + alpha * identity);
        while (solver.info() == Eigen::NumericalIssue) {
            std::cout << "Numbor" << std::endl;
            alpha *= 10;
            solver.compute(gauss_newton_hessian + alpha * identity);
        }

        du_gauss_newton = solver.solve(-dLdu);
        std::cout << "norm GN " << (gauss_newton_hessian * du_gauss_newton + dLdu).norm() << std::endl;
    }

    SparseMatrixF d2Edy2, d2Edydu;
    energyYHessianFunction(u, y, d2Edy2);
    energyMixedHessianFunction(u, y, d2Edydu);
    VectorXF pLpu, pLpy;
    gradientFunction(u, y, pLpu, pLpy);
    SparseMatrixF p2Lpu2, p2Lpypu, p2Lpy2;
    hessianFunction(u, y, p2Lpu2, p2Lpypu, p2Lpy2);

    int nu = u.rows(), ny = y.rows();

    TripletListF sgn_triplets;
    addTripletsBlock(sgn_triplets, p2Lpy2, 0, 0);
    addTripletsBlock(sgn_triplets, p2Lpypu, 0, ny);
    addTripletsBlock(sgn_triplets, d2Edy2, 0, ny + nu);

    addTripletsBlockTranspose(sgn_triplets, p2Lpypu, ny, 0);
    addTripletsBlock(sgn_triplets, p2Lpu2, ny, ny);
    addTripletsBlockTranspose(sgn_triplets, d2Edydu, ny, ny + nu);

    addTripletsBlock(sgn_triplets, d2Edy2, ny + nu, 0);
    addTripletsBlock(sgn_triplets, d2Edydu, ny + nu, ny);

    int size = 2 * ny + nu;
    SparseMatrixF sgn_matrix(size, size);
    sgn_matrix.setFromTriplets(sgn_triplets.begin(), sgn_triplets.end());

    SparseMatrixF H(sgn_matrix.rows(), sgn_matrix.cols());
    TripletListF H_triplets;
    addIdentityBlock(H_triplets, 0, 0, ny, 1e-8);
    addIdentityBlock(H_triplets, ny, ny, nu, 0);
    addIdentityBlock(H_triplets, ny + nu, ny + nu, ny, -1e-8);
    // addIdentityBlock(H_triplets, 0, 0, 2 * ny + nu, 0.0001);
    H.setFromTriplets(H_triplets.begin(), H_triplets.end());

    VectorXF sgn_rhs = VectorXF::Zero(size);
    sgn_rhs.segment(ny, nu) = -dLdu;

    // CholmodSolver solver;
    // solver.analyzePattern(sgn_matrix + H);
    // F alpha = 1e-10;
    // for (int i = 0; i < 50; i++) {
    //     solver.factorize(sgn_matrix + alpha * H);
    //     /// Should have NumericalIssue if K is not positive definite.
    //     if (solver.info() == Eigen::NumericalIssue) {
    //         alpha *= 10;
    //         continue;
    //     }
    // }
    // search_direction = solver.solve(sgn_rhs).segment(ny, nu);

    // PardisoLU solver;
    // solver.compute(sgn_matrix);
    // auto solution = solver.solve(sgn_rhs);
    // search_direction = solution.segment(ny, nu);
    // auto BB = sgn_matrix * solution - sgn_rhs;
    // std::cout << "normo " << BB.norm() << std::endl;

    PardisoLU solver;
    solver.compute(sgn_matrix + H);
    auto solution = solver.solve(sgn_rhs);
    search_direction = solution.segment(ny, nu);
    auto BB = (sgn_matrix + H) * solution - sgn_rhs;
    std::cout << "normo " << BB.norm() << std::endl;
    std::cout << "normo2 " << (du_gauss_newton - search_direction).norm() << std::endl;

    return true;
}

bool SensitivityAnalysis::lineSearch(VectorXF &u, const VectorXF &du, F initial_objective_value) {
    VectorXF y0;
    getYFunction(y0);

    bool success = false;
    F alpha = 2.0;
    for (int i = 0; i < 50; i++) {
        applyYFunction(y0);
        alpha *= 0.5;

        VectorXF u_line_search = u + alpha * du;
        if (!applyUFunction(u_line_search)) {
            continue;
        }

        Optimization::OptimizationStatus status = runEnergyMinimization();
        if (status == Optimization::FAILURE) {
            continue;
        }
        VectorXF y;
        getYFunction(y);

        F new_objective_value;
        if (objectiveFunction(u_line_search, y, new_objective_value)) {
            if (new_objective_value < initial_objective_value) {
                u = u_line_search;
                success = true;
                break;
            }
        } else {
            /// Model generation failed. Go to next loop iteration.
        }
    }

    if (!success) {
        applyYFunction(y0);
        applyUFunction(u);
    }
    return success;
}

bool SensitivityAnalysis::lineSearch(VectorXF &u, const VectorXF &du) {
    F initial_objective_value;
    applyUFunction(u);
    Optimization::OptimizationStatus status = runEnergyMinimization();
    if (status == Optimization::FAILURE) {
        return false;
    }
    VectorXF y;
    getYFunction(y);
    objectiveFunction(u, y, initial_objective_value);

    return lineSearch(u, du, initial_objective_value);
}

Optimization::OptimizationStatus SensitivityAnalysis::runEnergyMinimization() {
    Optimization::OptimizationStatus status;
    do {
        status = energyMinimizationStepFunction();
    } while (status == Optimization::SUCCESS);
    return status;
}

bool SensitivityAnalysis::checkGradient(const VectorXF &u, const VectorXF &du, const VectorXF &y0, F epsilon,
                                        int print_level) {
    VectorXF y = y0;
    applyUFunction(u);
    applyYFunction(y);

    bool success = false;

    do {
        Optimization::OptimizationStatus status;
        status = energyMinimizationStepFunction();

        if (status != Optimization::CONVERGED) {
            std::cout << "SA gradient check: y0 not an equilibrium state for given u. Running energy minimization."
                      << std::endl;
            status = runEnergyMinimization();
            getYFunction(y);
        }

        if (status == Optimization::FAILURE) {
            std::cout << "SA gradient check: failed to find equilibrium state." << std::endl;
            break;  // return false
        }
        assert(status == Optimization::CONVERGED);

        F L_0;
        objectiveFunction(u, y, L_0);
        VectorXF grad;
        computeSAGradient(u, y, grad);

        VectorXF u_plus = u + du.normalized() * epsilon;
        VectorXF u_minus = u - du.normalized() * epsilon;
        VectorXF y_plus, y_minus;
        F L_plus, L_minus;

        applyUFunction(u_plus);
        applyYFunction(y);
        status = runEnergyMinimization();
        if (status == Optimization::FAILURE) {
            std::cout << "SA gradient check: failed to find equilibrium state for u_plus." << std::endl;
            break;  // return false
        }
        getYFunction(y_plus);
        objectiveFunction(u_plus, y_plus, L_plus);

        applyUFunction(u_minus);
        applyYFunction(y);
        status = runEnergyMinimization();
        if (status == Optimization::FAILURE) {
            std::cout << "SA gradient check: failed to find equilibrium state for u_minus." << std::endl;
            break;  // return false
        }
        getYFunction(y_minus);
        objectiveFunction(u_minus, y_minus, L_minus);

        F val_computed = grad.transpose() * (u_plus - u_minus);
        F val_fd = L_plus - L_minus;
        std::cout << "Gradient check: Computed = " << val_computed << ", FD = " << val_fd << std::endl;

        success = true;
    } while (0);

    // Restore previous state.
    applyUFunction(u);
    applyYFunction(y0);

    return success;
}

void SensitivityAnalysis::makeConfigMenu() {
    std::vector<std::string> optimizer_names = getOptimizerNames();
    if (ImGui::Combo("Optimizer##SA", reinterpret_cast<int *>(&optimizer), optimizer_names)) {
        reset_optimization = true;
        partial_step_state = PARTIAL_STEP_START;
    }
    if (optimizer == SA_NLOPT) {
        std::vector<std::string> nlopt_names;
        for (int i = 0; i < nlopt::NUM_ALGORITHMS; i++) {
            nlopt_names.emplace_back(nlopt::algorithm_name((nlopt::algorithm)i));
        }
        if (ImGui::Combo("Algorithm##SA_NLOptOptimizer", reinterpret_cast<int *>(&nlopt_algorithm), nlopt_names)) {
            reset_optimization = true;
            partial_step_state = PARTIAL_STEP_START;
        }
    }
    ImGui::Text("Tolerance: 10^");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(100);
    ImGui::InputInt("##ToleranceSAConvergence", &tolerance_exponent);
    ImGui::Checkbox("Show partial steps", &partial_steps);
}
