#pragma once

#include "CRLHelper/CRLSubApp.h"
#include "CRLHelper/Optimization.h"
#include "CRLHelper/SensitivityAnalysis.h"
#include "CRLHelper/Solvers.h"

#include "Projects/Template/include/Model/MassSpring.h"
#include "Projects/Template/include/App/SplashScreen.h" 

class TemplateApp : public CRLSubApp {
   public:
    /// Model
    // MassSpring model;

    ///
    Simulation sim;

    /// Optimizer
    bool optimize = false;
    Optimization optimization;
    int convergence_tolerance_exponent = -16;

    /// Solver
    bool runSimulation = false;     // A toggle to decide if we should step the simulation
    Solver* solver = nullptr;       // Pointer to the chosen solver
    std::string solverType = "Forward Euler"; // or "Backward Euler", etc.

    /// Analysis
    int check_gradient_epsilon_exponent = -2;
    bool check_gradient_print_all = false;

    /// Camera state
    CRLCamera app_camera;

   public:
   Optimization::OptimizationStatus energyMinimizationStep();
   Optimization::OptimizationStatus energyMinimizationStepDyn();

   public:
    /// Template Builder functions
    TemplateApp() = default;
    
    TemplateApp(const SplashScreenResult &initParam);

    /// Template related functions
    void initializeSubApp() override;

    void mainLoop() override;

    void makeConfigWindow() override;

    void makeRunCheckWindow() override;

    void checkGradient(int order = 1);

    bool callbackKeyPressed(const CRLControlState &control_state, int key) override;

    [[nodiscard]] std::string getName() const override { return "Template"; }

    void createOrUpdateSolver();

   public:
    void getViewerData(std::vector<CRLViewerData> &viewer_data, CRLCamera &viewer_camera) override;
    
};
