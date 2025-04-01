#ifndef SPLASH_SCREEN_RESULT_H
#define SPLASH_SCREEN_RESULT_H

struct SplashScreenResult {
    // For interactive mode:
    bool use3D;         // Select 2D (false) or 3D (true)
    bool exportMode;    // If true, run export (batch) mode instead of launching the interactive app.
    double timeStep;          // Simulation time step.
    bool dynamic;             // Whether to use the dynamic formulation.
    int exponent_convergence_threshold;
    bool viscosity;            // Whether friction is enabled.
    int numSimulations;       // Number of simulation runs.
    double simulationTime;    // Duration (in seconds) for each simulation run.
    double animationStartTime;
    int scenarioShapeIndex = 0;
};

#endif // SPLASH_SCREEN_RESULT_H