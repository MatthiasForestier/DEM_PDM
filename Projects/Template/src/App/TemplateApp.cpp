#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

#include "Projects/Template/include/App/TemplateApp.h"

#include "CRLHelper/CameraHelper.h"
#include "CRLHelper/Optimization.h"
#include "CRLHelper/CRLTimer.h"
#include "CRLHelper/Solvers.h"

#include <iostream>
#include <random>


TemplateApp::TemplateApp(const SplashScreenResult &initParam) {
    // Simply set use3D in sim based on the splash screen's flag
    sim.use3D = initParam.use3D;
}

void TemplateApp::initializeSubApp() {
    /// Initialize camera.
    app_camera.eye = Vector3F(0, 0, 5);
    app_camera.center = Vector3F(0, 0, 0);
    app_camera.set_up_direction(Vector3F(0, 1, 0));
    app_camera.height = 1.5; //Change the distance from sim
    //Initialize particles
    sim.particles2D = sim.createRandomParticles2D();
    createOrUpdateSolver();
}

void TemplateApp::makeConfigWindow() {
    if (ImGui::CollapsingHeader("Optimization", ImGuiTreeNodeFlags_DefaultOpen)) {
        optimization.makeConfigMenu();
        ImGui::Checkbox("Optimize", &optimize);
    }

    ImGui::Spacing();
    ImGui::Spacing();

    if (ImGui::CollapsingHeader("Simulation", ImGuiTreeNodeFlags_DefaultOpen)) {
        sim.makeConfigMenu();
    }
}

void TemplateApp::makeRunCheckWindow() {
    // --- Existing gradient-check controls ---
    if (ImGui::Button("Check Gradient")) {
        checkGradient(1);
    }
    ImGui::SameLine();
    if (ImGui::Button("Check Hessian")) {
        checkGradient(2);
    }

    ImGui::Text("Epsilon: 10^");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(100);
    ImGui::InputInt("##check_gradient_epsilon", &check_gradient_epsilon_exponent);
    ImGui::SameLine();
    ImGui::Checkbox("Print All##0", &check_gradient_print_all);

    ImGui::Separator(); // a horizontal separator

    // --- New UI for Coloring a Particle Red ---
    static int particleID = 0;
    ImGui::Separator();
    ImGui::Text("Color Particle Red");
    ImGui::InputInt("Particle ID", &particleID);

    if (ImGui::Button("Set Particle Red")) {
        if (particleID >= 0 && particleID < static_cast<int>(sim.particles2D.size())) {
            sim.colorParticleRed(particleID);
        } else {
            std::cerr << "Invalid Particle ID: " << particleID << "\n";
        }
    }

    ImGui::Separator(); // a horizontal separator

    if (ImGui::CollapsingHeader("Dynamic Formulation", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::Checkbox("Dynamic", &sim.dynamic);
    }

    // --- New solver selection controls ---
    // static const char* solverNames[] = {"Forward Euler", "Backward Euler"};
    // static int currentSolverIndex = 0;
    // if (ImGui::Combo("Solver Type", &currentSolverIndex, solverNames, IM_ARRAYSIZE(solverNames))) {
    //     // Update solverType based on user choice
    //     solverType = solverNames[currentSolverIndex];
    //     // Recreate or update the solver if needed
    //     createOrUpdateSolver();
    // }

    // // --- Run / Pause controls ---
    // if (ImGui::Button("Run Simulation")) {
    //     runSimulation = true;  // or any other boolean you use to indicate "simulation running"
    // }
    // ImGui::SameLine();
    // if (ImGui::Button("Pause Simulation")) {
    //     runSimulation = false;
    // }
}

void TemplateApp::createOrUpdateSolver() {
    // Clean up if we have an existing solver
    if (solver) {
        delete solver;
        solver = nullptr;
    }

    // Use the static factory method we defined in Solver.h
    solver = Solver::createSolver(solverType, &sim);
}


void TemplateApp::mainLoop() {
    if (optimize) {
        Optimization::OptimizationStatus status = energyMinimizationStep();
    }
    if (runSimulation && solver) {
        solver->stepMultipleTimes();
    }
    if (optimize && sim.dynamic){
        Optimization::OptimizationStatus status = energyMinimizationStepDyn();
    }
}

Optimization::OptimizationStatus TemplateApp::energyMinimizationStep() {
    // Get the global state vector representing all particles.
    VectorXF globalState = sim.getGlobalState();

    // Set up the optimization functions to work on the global state:
    optimization.objective_function = [&](const VectorXF &y, F &value) {
        sim.setGlobalState(y);
        sim.updateAuxiliaryStructures();
        sim.compute_energy(value);
        return true;
    };
    optimization.gradient_function = [&](const VectorXF &y, F &value, VectorXF &gradient) {
        sim.setGlobalState(y);
        sim.updateAuxiliaryStructures();
        sim.compute_energy(value);
        sim.compute_gradient(gradient);
        return true;
    };
    optimization.hessian_function = [&](const VectorXF &y, F &value, VectorXF &gradient, HessianF &hessian) {
        sim.setGlobalState(y);
        sim.updateAuxiliaryStructures();
        sim.compute_energy(value);
        sim.compute_gradient(gradient);
        sim.compute_hessian(hessian.A);
        return true;
    };

    // Run one optimization step on the global state.
    auto status = optimization.step(globalState);

    // Update the simulation with the optimized state.
    sim.setGlobalState(globalState);
    sim.updateAuxiliaryStructures();

    return status;
}

Optimization::OptimizationStatus TemplateApp::energyMinimizationStepDyn() {
    // Get the global state vector representing all particles.
    VectorXF globalState = sim.getGlobalState();

    // Set up the optimization functions to work on the global state:
    optimization.objective_function = [&](const VectorXF &y, F &value) {
        sim.setGlobalState(y);
        sim.updateAuxiliaryStructures();
        sim.compute_energy(value);
        return true;
    };
    optimization.gradient_function = [&](const VectorXF &y, F &value, VectorXF &gradient) {
        sim.setGlobalState(y);
        sim.updateAuxiliaryStructures();
        sim.compute_energy(value);
        sim.compute_gradient(gradient);
        return true;
    };
    optimization.hessian_function = [&](const VectorXF &y, F &value, VectorXF &gradient, HessianF &hessian) {
        sim.setGlobalState(y);
        sim.updateAuxiliaryStructures();
        sim.compute_energy(value);
        sim.compute_gradient(gradient);
        sim.compute_hessian(hessian.A);
        return true;
    };
    
    // Run one optimization step on the global state.
    auto status = optimization.step(globalState);

    // Update the simulation with the optimized state.
    sim.setGlobalState(globalState);
    sim.updateAuxiliaryStructures();

    return status;
}

void TemplateApp::checkGradient(int order) {
    // Define the objective function to work on the global state.
    auto objective_function = [&](const VectorXF &y, F &value) {
        sim.setGlobalState(y);        // Update the simulation's particles
        sim.compute_energy(value);    // Compute global energy for all particles
        return true;
    };

    // Define the gradient function for the global state.
    auto gradient_function = [&](const VectorXF &y, VectorXF &gradient) {
        sim.setGlobalState(y);
        sim.compute_gradient(gradient);
        return true;
    };

    // Define the Hessian function.
    auto hessian_function = [&](const VectorXF &y, MatrixXF &hessian) {
        sim.setGlobalState(y);
        SparseMatrixF hessian_sparse;
        sim.compute_hessian(hessian_sparse);
        hessian = hessian_sparse.toDense();
        return true;
    };

    // Get the current global state vector (which contains all particles' states).
    VectorXF curr_y = sim.getGlobalState();

    // Epsilon for finite differences is set based on a given exponent.
    F epsilon = std::pow(10.0, check_gradient_epsilon_exponent);
    int print_level = (check_gradient_print_all ? 2 : 1);

    if (order == 1) {
        VectorXF error;
        Optimization::checkGradient(curr_y, error, objective_function, gradient_function, epsilon, print_level);
    }
    else if (order == 2) {
        MatrixXF error;
        Optimization::checkHessian(curr_y, error, gradient_function, hessian_function, epsilon, print_level);
    }
}

void TemplateApp::getViewerData(std::vector<CRLViewerData>& viewer_data, CRLCamera& viewer_camera) {
    // Copy camera state.
    viewer_camera = app_camera;

    /////////////
    // PARTICLE DISKS
    /////////////
    viewer_data.emplace_back();
    CRLViewerData& particleViewerData = viewer_data.back();

    // Use the actual number of particles.
    int numParticles = static_cast<int>(sim.particles2D.size());
    int numSegments = 32;  // number of segments per disk
    // Each disk: (numSegments + 1) vertices (center + perimeter) and numSegments faces.
    int totalParticleVertices = numParticles * (numSegments + 1);
    int totalParticleFaces = numParticles * numSegments;

    // Preallocate the matrices.
    particleViewerData.mesh_v.resize(totalParticleVertices, 3);
    particleViewerData.mesh_f.resize(totalParticleFaces, 3);

    int vertexOffset = 0;
    int faceOffset = 0;
    for (int i = 0; i < numParticles; i++) {
        const Particle2D& p = sim.particles2D[i];
        Vector3F center = p.pos;
        F radius = p.radius;

        // Write the center vertex.
        particleViewerData.mesh_v.row(vertexOffset) = center.transpose();
        int diskStart = vertexOffset; // starting index for this disk
        vertexOffset++;

        // Compute and write the perimeter vertices.
        for (int j = 0; j < numSegments; j++) {
            F theta = 2.0 * M_PI * j / numSegments;
            F x = radius * std::cos(theta);
            F y = radius * std::sin(theta);
            Vector3F vertex = center + Vector3F(x, y, 0.0f);
            particleViewerData.mesh_v.row(vertexOffset) = vertex.transpose();
            vertexOffset++;
        }

        // Create faces using a triangle fan.
        for (int j = 1; j < numSegments; j++) {
            particleViewerData.mesh_f.row(faceOffset) =
                Vector3I(diskStart, diskStart + j, diskStart + j + 1).transpose();
            faceOffset++;
        }
        // Close the fan.
        particleViewerData.mesh_f.row(faceOffset) =
            Vector3I(diskStart, diskStart + numSegments, diskStart + 1).transpose();
        faceOffset++;
    }

    // Optionally color the particle disks using each particle's own color.
    particleViewerData.mesh_c = MatrixXF::Constant(totalParticleFaces, 3, 1.0);
    int currentFace = 0;
    for (int i = 0; i < numParticles; i++) {
        const Particle2D& p = sim.particles2D[i];
        for (int j = 0; j < numSegments; j++) {
            particleViewerData.mesh_c.row(currentFace) << p.color.r, p.color.g, p.color.b;
            currentFace++;
        }
    }

    /////////////
    // SCENARIO OBJECTS
    /////////////
    viewer_data.emplace_back();
    CRLViewerData& scenarioViewerData = viewer_data.back();

    int totalScenarioVertices = 0;
    int totalScenarioFaces = 0;
    for (size_t soIdx = 0; soIdx < sim.scenarioObjects.size(); soIdx++) {
        const ScenarioObject& so = *(sim.scenarioObjects[soIdx]);
        if (so.vertices.empty()) continue;
        totalScenarioVertices += (1 + static_cast<int>(so.vertices.size()));
        totalScenarioFaces += static_cast<int>(so.vertices.size());
    }

    scenarioViewerData.mesh_v.resize(totalScenarioVertices, 3);
    scenarioViewerData.mesh_f.resize(totalScenarioFaces, 3);

    int scenarioVertexOffset = 0;
    int scenarioFaceOffset = 0;
    for (size_t soIdx = 0; soIdx < sim.scenarioObjects.size(); soIdx++) {
        const ScenarioObject& so = *(sim.scenarioObjects[soIdx]);
        if (so.vertices.empty()) continue;

        int currentStart = scenarioVertexOffset;
        scenarioViewerData.mesh_v.row(scenarioVertexOffset) = so.position.transpose();
        scenarioVertexOffset++;

        for (size_t i = 0; i < so.vertices.size(); i++) {
            scenarioViewerData.mesh_v.row(scenarioVertexOffset) = so.vertices[i].transpose();
            scenarioVertexOffset++;
        }

        int N = static_cast<int>(so.vertices.size());
        for (int i = 1; i < N; i++) {
            scenarioViewerData.mesh_f.row(scenarioFaceOffset) =
                Vector3I(currentStart, currentStart + i, currentStart + i + 1).transpose();
            scenarioFaceOffset++;
        }
        scenarioViewerData.mesh_f.row(scenarioFaceOffset) =
            Vector3I(currentStart, currentStart + N, currentStart + 1).transpose();
        scenarioFaceOffset++;
    }

    scenarioViewerData.mesh_c = MatrixXF::Constant(totalScenarioFaces, 3, 0.0);
}




bool TemplateApp::callbackKeyPressed(const CRLControlState& control_state, int key) {
    switch (key) {
        case GLFW_KEY_SPACE:
            optimize = !optimize;
            break;
        default:
            break;
    }
    return false;
}



