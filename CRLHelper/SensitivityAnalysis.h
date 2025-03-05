#pragma once

#include "CRLHelper/Optimization.h"

#include "ThirdParty/mma/src/mma/MMASolver.h"
#include "ThirdParty/mma/src/gcmma/GCMMASolver.h"

class MMASolver;
class GCMMASolver;

class SensitivityAnalysis {
   public:
    enum SAOptimizer {
        SA_GRADIENT_DESCENT,  //
        SA_BFGS,
        SA_LBFGS,
        SA_GAUSS_NEWTON,
        SA_SGN,
        SA_GCMMA,
        SA_NLOPT,
        NUM_OPTIMIZERS
    };

    enum PartialStepState {
        PARTIAL_STEP_START,  //
        PARTIAL_STEP_INITIAL_EQUILIBRIUM,
        PARTIAL_STEP_SEARCH_DIRECTION,
        PARTIAL_STEP_LINE_SEARCH_OUTER,
        PARTIAL_STEP_LINE_SEARCH_EQUILIBRIUM
    };

   public:
    static std::vector<std::string> getOptimizerNames() {
        return {"Gradient Descent",          //
                "BFGS",                      //
                "LBFGS",                     //
                "Gauss-Newton",              //
                "SGN (Not working)",         //
                "GCMMA (No partial steps)",  //
                "NLopt (No partial steps)"};
    }

   public:
    SAOptimizer optimizer = SA_GRADIENT_DESCENT;
    int tolerance_exponent = -10;
    MatrixXF bfgs_inverse_hessian;
    VectorXF prev_u;  /// For (L-)BFGS
    VectorXF prev_g;  /// For (L-)BFGS
    LBFGSpp::BFGSMat<F> lbfgs_mat;
    int lbfgs_m = 6;
    bool reset_optimization = true;  /// Reset e.g. BFGS Hessian on next step.

    std::shared_ptr<GCMMASolver> gcmma_solver;

    nlopt::opt nlopt_opt;
    int nlopt_algorithm = 24;

    bool partial_steps = false;
    PartialStepState partial_step_state = PARTIAL_STEP_START;
    VectorXF partial_step_search_direction;
    F partial_step_initial_objective;
    int partial_step_line_search_step;
    VectorXF partial_step_u0;
    VectorXF partial_step_y0;

   public:
    std::function<bool(VectorXF &)> getUFunction;          // Args u
    std::function<bool(const VectorXF &)> applyUFunction;  // Args u
    std::function<bool(VectorXF &)> getYFunction;          // Args y
    std::function<bool(const VectorXF &)> applyYFunction;  // Args y
    std::function<Optimization::OptimizationStatus()> energyMinimizationStepFunction;
    std::function<bool(const VectorXF &, const VectorXF &, F &)> objectiveFunction;  // Args u, y, L
    std::function<bool(const VectorXF &, const VectorXF &, VectorXF &, VectorXF &)>
        gradientFunction;  // Args u, y, pLpu, pLpy
    std::function<bool(const VectorXF &, const VectorXF &, SparseMatrixF &, SparseMatrixF &, SparseMatrixF &)>
        hessianFunction;  // Args u, y, p2Lpu2, p2Lpypu, p2Lpy2
    std::function<bool(const VectorXF &, const VectorXF &, SparseMatrixF &)>
        energyYHessianFunction;  // Args u, y, p2Epy2
    std::function<bool(const VectorXF &, const VectorXF &, SparseMatrixF &)>
        energyMixedHessianFunction;  // Args u, y, p2Epypu
    std::function<bool(const VectorXF &, const VectorXF &, MatrixXF &)>
        computeSensitivityMatrixOverride;  // Args u, y, dydu

    SensitivityAnalysis() {
        getUFunction = [&](VectorXF &) { return false; };
        applyUFunction = [&](const VectorXF &) { return false; };
        getYFunction = [&](VectorXF &) { return false; };
        applyYFunction = [&](const VectorXF &) { return false; };
        energyMinimizationStepFunction = [&]() { return Optimization::FAILURE; };
        objectiveFunction = [&](const VectorXF &, const VectorXF &, F &) { return false; };
        gradientFunction = [&](const VectorXF &, const VectorXF &, VectorXF &, VectorXF &) { return false; };
        hessianFunction = [&](const VectorXF &, const VectorXF &, SparseMatrixF &, SparseMatrixF &, SparseMatrixF &) {
            return false;
        };
        energyYHessianFunction = [&](const VectorXF &, const VectorXF &, SparseMatrixF &) { return false; };
        energyMixedHessianFunction = [&](const VectorXF &, const VectorXF &, SparseMatrixF &) { return false; };
        computeSensitivityMatrixOverride = [&](const VectorXF &, const VectorXF &, MatrixXF &) { return false; };
    }

   public:
    /// These optimization routines take lambda functions for computing the objective value and its gradients.
    /// Gradient and hessian functions also compute objective value and lower-order derivatives.
    Optimization::OptimizationStatus getDirectionGradientDescent(const VectorXF &u, VectorXF &du,
                                                                 F &initial_objective_value);

    Optimization::OptimizationStatus getDirectionBFGS(const VectorXF &u, VectorXF &du, F &initial_objective_value);

    Optimization::OptimizationStatus getDirectionLBFGS(const VectorXF &u, VectorXF &du, F &initial_objective_value);

    Optimization::OptimizationStatus getDirectionGaussNewton(const VectorXF &u, VectorXF &du,
                                                             F &initial_objective_value);

    Optimization::OptimizationStatus getDirectionSGN(const VectorXF &u, VectorXF &du, F &initial_objective_value);

    Optimization::OptimizationStatus stepGCMMA(VectorXF &u);

    Optimization::OptimizationStatus stepNLopt(VectorXF &u);

    Optimization::OptimizationStatus step();
    Optimization::OptimizationStatus fullStep();
    Optimization::OptimizationStatus partialStep();

    bool computeSensitivityMatrix(const VectorXF &u, const VectorXF &y, MatrixXF &dydu);

    bool computeSAGradient(const VectorXF &u, const VectorXF &y, VectorXF &gradient);
    bool computeSAGradient(const VectorXF &u, const VectorXF &y, const MatrixXF &dydu, VectorXF &gradient);

    bool computeSGNSearchDirection(const VectorXF &u, const VectorXF &y, const VectorXF &dLdu,
                                   VectorXF &search_direction);

    /// Find step along search direction du that decreases the objective value. Update u accordingly.
    bool lineSearch(VectorXF &u, const VectorXF &du, F initial_objective_value);

    bool lineSearch(VectorXF &u, const VectorXF &du);

    bool checkGradient(const VectorXF &u, const VectorXF &du, const VectorXF &y0, F epsilon, int print_level = 0);

    /// Run energy minimization to convergence.
    Optimization::OptimizationStatus runEnergyMinimization();

   public:
    void makeConfigMenu();
};
