#pragma once

#include "CRLHelper/CRLSubApp.h"
#include "CRLHelper/Optimization.h"
#include "CRLHelper/SensitivityAnalysis.h"
#include "CRLHelper/Logger.h"

#include "Projects/Template/include/Model/MassSpring.h"
#include "Projects/Template/include/App/SplashScreen.h" 

class TemplateApp : public CRLSubApp {
   public:
    /// Model
    // MassSpring model;

    /// Main simulation object
    Simulation sim;

    /// Logger
    Logger logger;

    /// Optimizer
    bool optimize = false;
    Optimization optimization;
    I convergence_tolerance_exponent = -16;

    /// Dynamic Convergence
    I exponent_convergence_threshold = -2;
    F dynamic_convergence_threshold;
    I maxIter = 100; 
    I breach = 0;
    I calls = 0;

    VectorXF globalState_0;
    VectorXF globalState_1;
    VectorXF globalState_2;

    /// Analysis
    I check_gradient_epsilon_exponent = -4;
    bool check_gradient_print_all = false;

    /// Camera state
    CRLCamera app_camera;

   public:
   Optimization::OptimizationStatus energyMinimizationStep();
   void reinitializeGlobalState();
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

    void getViewerData(std::vector<CRLViewerData> &viewer_data, CRLCamera &viewer_camera) override;

    void showLoggerWindow() override;
    
};
