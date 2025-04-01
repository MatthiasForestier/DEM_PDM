#include <chrono>
#include <ctime>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <filesystem>
namespace fs = std::filesystem;

#include "Projects/Template/include/Model/MassSpring.h"       // Your Simulation class definitions
#include "Projects/Template/include/App/ExportSimulation.h"
#include "Projects/Template/include/App/TemplateApp.h"

// Helper function to get the current date/time as a string in the format YYYY-MM-DD_HH-MM-SS.
std::string getCurrentDateTimeString() {
    auto now = std::chrono::system_clock::now();
    std::time_t now_time = std::chrono::system_clock::to_time_t(now);
    
    std::tm tm;
#ifdef _WIN32
    localtime_s(&tm, &now_time);
#else
    localtime_r(&now_time, &tm);
#endif

    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y-%m-%d_%H-%M-%S");
    return oss.str();
}

void exportSimulationBehavior(const SplashScreenResult &params) {
    // Create a TemplateApp instance using the splash screen parameters.
    TemplateApp app(params);
    
    // Override simulation parameters with export-mode values.
    app.sim.timeStep = params.timeStep;
    app.sim.dynamic = params.dynamic;
    app.dynamic_convergence_threshold = pow(10.0, params.exponent_convergence_threshold);
    app.sim.viscosity = params.viscosity;
    // IMPORTANT: If you're using 2D particles, make sure use3D is false.
    app.sim.use3D = false;  
    app.sim.endTime = params.simulationTime;

    app.sim.scenarioObjects.clear();
    if (params.scenarioShapeIndex == 0) {
        // Create a Circle scenario object with 64 segments, radius 1.1, centered at the origin.
        app.sim.scenarioObjects.push_back(std::make_unique<Circle>(1.1f, 64, Vector3F(0.0f, 0.0f, 0.0f)));
    } else {
        // Create a Square scenario object with half-dimensions 1.1 and 0.9, centered at the origin.
        app.sim.scenarioObjects.push_back(std::make_unique<Square>(1.1f, 0.9f, Vector3F(0.0f, 0.0f, 0.0f)));
    }
    

    // Ensure output directory exists.
    const std::string outputDir = "output";
    if (!fs::exists(outputDir)) {
        if (!fs::create_directories(outputDir)) {
            std::cerr << "Error creating output directory: " << outputDir << std::endl;
            return;
        }
    }
    std::string dateTime = getCurrentDateTimeString();
    std::string outputFile = outputDir + "/simulation_export_" + dateTime + ".txt";
    std::ofstream outfile(outputFile);
    if (!outfile) {
        std::cerr << "Error opening output file: " << outputFile << std::endl;
        return;
    }
    
    // Write header information.
    outfile << "Exporting Simulation Behavior\n";
    outfile << "Parameters:\n";
    outfile << "  Number of Particles: " << app.sim.numParticles << "\n";
    outfile << "  Gravity: " << app.sim.gravity(1) << "\n";
    outfile << "  Time Step: " << app.sim.timeStep << "\n";
    outfile << "  Dynamic: " << (app.sim.dynamic ? "true" : "false") << "\n";
    outfile << "  Viscosity: " << (app.sim.viscosity ? "true" : "false") << "\n";
    outfile << "  Simulation Duration: " << app.sim.endTime << "\n";
    outfile << "  Animation Start Time: " << params.animationStartTime << "\n";
    outfile << "  Number of Simulation Runs: " << params.numSimulations << "\n";
    outfile << "\n";
    
    // Outer loop: run the simulation multiple times.
    for (int run = 0; run < params.numSimulations; run++) {
        outfile << "=== Run " << (run + 1) << " ===\n";
        
        // Reinitialize simulation state:
        // (Recreate particles, reinsert into grid, and update the global state.)
        app.sim.buildGridDataStructure();
        app.sim.particles2D = app.sim.createRandomParticles2D();
        app.sim.insertParticlesIntoGrid();
        app.reinitializeGlobalState();
        
        // Reset time and animation flag for each run.
        double t = 0.0;
        bool animationTriggered = false;
        
        // Simulation loop for one run.
        while (t < app.sim.endTime) {
            if (!animationTriggered && t >= params.animationStartTime) {
                app.sim.startScenarioAnimation();
                animationTriggered = true;
                std::cout << "Scenario animation started at t = " << t << " seconds.\n";
            }
            
            // Perform optimization step(s).
            if (app.sim.dynamic) {
                Optimization::OptimizationStatus status = app.energyMinimizationStepDyn();
            } else {
                Optimization::OptimizationStatus status = app.energyMinimizationStep();
            }
            
            // Update scenario animation and auxiliary structures.
            app.sim.updateScenarioAnimation(app.sim.timeStep);
            double energy;
            app.sim.compute_energy(energy);
            
            // Log the current state.
            auto state = app.sim.getGlobalState();
            outfile << "Time: " << t << "  Energy: " << energy << "\n";
            outfile << "State: " << state.transpose() << "\n";
            
            // Advance simulation time.
            t += app.sim.timeStep;
        }
        
        outfile << "\n";
    }
    
    outfile.close();
    std::cout << "Simulation behavior exported to " << outputFile << std::endl;
}
