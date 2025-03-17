#pragma once

#include "CRLHelper/VecMatDef.h"
#include "CRLHelper/Hessian.h"

#include "ThirdParty/LBFGSpp/include/LBFGSpp/BFGSMat.h"

#include "nlopt.hpp"

class Optimization {
   public:
    enum OptimizationStatus {
        CONVERGED,  //
        SUCCESS,
        FAILURE,
        PARTIAL
    };

    enum Optimizer {
        GRADIENT_DESCENT,  //
        NEWTON,
        NEWTON_WOODBURY,
        BFGS,
        LBFGS,
        NLOPT,
        NUM_OPTIMIZERS
    };

    
   public:
    static std::vector<std::string> getOptimizerNames() {
        return {"Gradient Descent",      //
                "Newton",                //
                "Newton with Woodbury",  //
                "BFGS",                  //
                "LBFGS",                 //
                "NLOpt"};
    }

   public:
    Optimizer optimizer = NEWTON;
    int tolerance_exponent = -16;
    MatrixXF bfgs_inverse_hessian;
    LBFGSpp::BFGSMat<F> lbfgs_mat;
    int lbfgs_m = 6;
    bool reset_optimization = true;  /// Reset e.g. BFGS Hessian on next step.

    nlopt::opt nlopt_opt;
    int nlopt_algorithm = 24;

   public:
    /// Gradient and hessian functions also compute objective value and lower-order derivatives.
    std::function<bool(const VectorXF &, F &)> objective_function;
    std::function<bool(const VectorXF &, F &, VectorXF &)> gradient_function;
    std::function<bool(const VectorXF &, F &, VectorXF &, HessianF &)> hessian_function;

    Optimization() {
        objective_function = [&](const VectorXF &y, F &energy) { return false; };
        gradient_function = [&](const VectorXF &y, F &energy, VectorXF &gradient) { return false; };
        hessian_function = [&](const VectorXF &y, F &energy, VectorXF &gradient, HessianF &hessian) { return false; };
    }

   public:
    /// Solve for x in (A + U * V^T)x = b. Return true on success.
    bool linearSolve(const HessianF &hessian, const VectorXF &b, VectorXF &x);

    /// Solve for x in (A + U * V^T)x = b, using Woodbury formula. Return true on success.
    bool linearSolveWoodbury(const HessianF &hessian, const VectorXF &b, VectorXF &x);

   private:
    OptimizationStatus stepGradientDescent(VectorXF &y);

    OptimizationStatus stepBFGS(VectorXF &y);

    OptimizationStatus stepLBFGS(VectorXF &y);

    OptimizationStatus stepNewton(VectorXF &y, bool use_woodbury);

    OptimizationStatus stepNLopt(VectorXF &y);

   public:
    OptimizationStatus step(VectorXF &y);

    /// Find step along search direction dy that decreases the objective value. Update y accordingly.
    bool lineSearch(VectorXF &y, const VectorXF &dy, F initial_objective_value);

    bool lineSearch(VectorXF &y, const VectorXF &dy);

    /// Check gradient for scalar function of vector argument.
    static bool checkGradient(const VectorXF &y, VectorXF &error,
                              const std::function<bool(const VectorXF &, F &)> &func,
                              const std::function<bool(const VectorXF &, VectorXF &)> &grad_func, F epsilon,
                              int print_level = 0);

    /// Check hessian for scalar function of vector argument.
    static bool checkHessian(const VectorXF &y, MatrixXF &error,
                             const std::function<bool(const VectorXF &, VectorXF &)> &grad_func,
                             const std::function<bool(const VectorXF &, MatrixXF &)> &hess_func, F epsilon,
                             int print_level = 0);

   public:
    void makeConfigMenu();
};
