#pragma once

#include "CRLHelper/CRLSubApp.h"
#include "CRLHelper/Optimization.h"
#include "CRLHelper/SensitivityAnalysis.h"
#include "CRLHelper/Logger.h"
#include "CRLHelper/Algebra.h"
#include "CRLHelper/VecMatDef.h"

#include "Projects/Template/include/Model/MassSpring.h"

struct SplashScreenResult {
    // For interactive mode:
    bool use3D;         // Select 2D (false) or 3D (true)
    bool exportMode;    // If true, run export (batch) mode instead of launching the interactive app.
    double timeStep;          // Simulation time step.
    bool dynamic;             // Whether to use the dynamic formulation.
    int exponent_convergence_threshold;
    bool viscosity;            // Whether friction is enabled.
    int numSimulations;       // Number of simulation runs.
    int numParticles;         // Number of particles in the simulation.
    double simulationTime;    // Duration (in seconds) for each simulation run.
    double animationStartTime;
    int scenarioShapeIndex = 0;
    bool periodicX;
    double V0;
    double L;
};

class TemplateApp : public CRLSubApp {
   public:
    /// Model
    // MassSpring model;

    /// Main simulation object
    Simulation sim;

    /// Global state
    VectorXF globalState_0;

    /// Logger
    Logger logger;

    /// Optimizer
    bool optimize = false;
    Optimization optimization;
    I convergence_tolerance_exponent = -16;

    /// Dynamic Convergence
    I exponent_convergence_threshold = -3;
    F dynamic_convergence_threshold;
    I maxIter = 1000; 
    I breach = 0;
    I calls = 0;

    /// Analysis
    I check_gradient_epsilon_exponent = -6;
    bool check_gradient_print_all = false;

    /// Camera state
    CRLCamera app_camera;

   public:
   Optimization::OptimizationStatus energyMinimizationStep();
   void reinitializeGlobalState();
   Optimization::OptimizationStatus energyMinimizationStepDyn();

   public:
    /// Template Builder functions
    TemplateApp(const SplashScreenResult &initParam);
    TemplateApp() = default;

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
