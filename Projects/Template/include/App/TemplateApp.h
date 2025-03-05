#pragma once

#include "CRLHelper/CRLSubApp.h"
#include "CRLHelper/Optimization.h"
#include "CRLHelper/SensitivityAnalysis.h"

#include "Projects/Template/include/Model/MassSpring.h"

class TemplateApp : public CRLSubApp {
   public:
    /// Model
    MassSpring model;

    /// Optimizer
    bool optimize = false;
    Optimization optimization;
    int convergence_tolerance_exponent = -16;

    /// Analysis
    int check_gradient_epsilon_exponent = -4;
    bool check_gradient_print_all = false;

    /// Camera state
    CRLCamera app_camera;

   public:
    Optimization::OptimizationStatus energyMinimizationStep();

   public:
    void initializeSubApp() override;

    void mainLoop() override;

    void makeConfigWindow() override;

    void makeAnalysisWindow() override;
    void checkGradient(int order = 1);

    bool callbackKeyPressed(const CRLControlState &control_state, int key) override;

    [[nodiscard]] std::string getName() const override { return "Template"; }

   public:
    void getViewerData(std::vector<CRLViewerData> &viewer_data, CRLCamera &viewer_camera) override;
};
