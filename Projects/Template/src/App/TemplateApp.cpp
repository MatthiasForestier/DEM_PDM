#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

#include "Projects/Template/include/App/TemplateApp.h"

#include "CRLHelper/CameraHelper.h"
#include "CRLHelper/Optimization.h"
#include "CRLHelper/CRLTimer.h"

#include <iostream>
#include <random>


TemplateApp::TemplateApp(const SplashScreenResult &initParam) {
    // Initialize your simulation parameters using initParam.
    // For instance:
    sim.use3D = initParam.use3D;
    sim.timeStep = initParam.timeStep;
    sim.dynamic = initParam.dynamic;
    sim.friction = initParam.viscosity;
    // ... initialize any other simulation parameters ...
}

void TemplateApp::initializeSubApp() {
    /// Initialize camera.
    app_camera.eye = Vector3F(0, 0, 5);
    app_camera.center = Vector3F(0, 0, 0);
    app_camera.set_up_direction(Vector3F(0, 1, 0));
    app_camera.height = 1.5; //Change the distance from sim
    //Initialize particles
    sim.particles2D = sim.createRandomParticles2D();
    for (auto &p : sim.particles2D) p.ix = 0;
    reinitializeGlobalState();
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

    if (ImGui::Button("Reload Particles")) {
        sim.particles2D = sim.createRandomParticles2D();
        for (auto &p : sim.particles2D) p.ix = 0;
        sim.insertParticlesIntoGrid();
        reinitializeGlobalState();
        breach = 0;
        calls = 0;
        logger.clear();
    }

    // --- New UI for Scenario Objects ---

    // Toggle flag for activation.
    static bool activateScenario = false;
    static bool lastFrameToggleState = false;

    // Provide a few shape choices.
    static int shapeIndex = 0;
    const char* shapes[] = { "Circle", "Square", "Tunnel" };

    // ImGui Checkbox for toggling on/off.
    if (ImGui::Checkbox("Activate Scenario Object", &activateScenario)) {
        // If the user has just toggled OFF, remove the scenario object(s) from the simulation.
        if (!activateScenario && lastFrameToggleState) {
            sim.scenarioObjects.clear();
            std::cout << "Scenario object(s) removed from simulation.\n";
        }
        lastFrameToggleState = activateScenario;
    }

    // If activated, show a combo for shape selection + create button.
    if (activateScenario) {
        ImGui::Combo("Shape Choice", &shapeIndex, shapes, IM_ARRAYSIZE(shapes));

        // Use a lambda function for creating the scenario object.
        auto scenarioObjectCreation = [&]() {
            sim.scenarioObjects.clear();
            switch (shapeIndex) {
                case 0:  sim.scenarioObjects.push_back(
                            std::make_unique<Circle>(1.1f, 64,
                                                    Vector3F(0,0,0)));        break;
                case 1:  sim.scenarioObjects.push_back(
                            std::make_unique<Square>(1.1f, 0.9f,
                                                    Vector3F(0,0,0)));        break;
                case 2:  sim.scenarioObjects.push_back(
                            std::make_unique<Tunnel2D>( /*halfLen =*/ 2.0f,
                                                        /*halfWid =*/ sim.L/2,
                                                        Vector3F(0,0,0)));     break;
            }
        };

        if (ImGui::Button("Create Scenario")) {
            scenarioObjectCreation();
            sim.buildGridDataStructure();
            sim.insertParticlesIntoGrid();
            std::cout << "New scenario object created and stored!\n";
        }

        if (ImGui::CollapsingHeader("Scenario Animation", ImGuiTreeNodeFlags_DefaultOpen)) {
            // Toggle animation on/off (this flag is also updated by startScenarioAnimation, so the checkbox is optional)
            // ImGui::Checkbox("Animate Scenario Object", &sim.animateScenario);
            // Button to start/reset the animation.
            if (ImGui::Button("Start Animation")) {
                scenarioObjectCreation();
                sim.startScenarioAnimation();
            }
            // Allow user to change the duration and distance.
            ImGui::InputDouble("Animation Duration (s)", &sim.animationDuration, 1.0, 5.0, "%.1f");
            ImGui::InputDouble("Animation Distance", &sim.animationDistance, 0.1, 1.0, "%.1f");
        }
    }


}

void TemplateApp::reinitializeGlobalState() {
    globalState_0 = sim.getGlobalState();
    sim.globalState_1 = globalState_0; // start with zero velocity/dynamic change
    sim.globalState_2 = globalState_0;
    // Optionally, update viewer data immediately here or in the next frame.
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
        ImGui::Text("Convergence Threshold: 10^");
        ImGui::SameLine();
        ImGui::SetNextItemWidth(100);
        ImGui::InputInt("## Convergence exponent:", &exponent_convergence_threshold);
        dynamic_convergence_threshold = pow(10.0, exponent_convergence_threshold);
        ImGui::InputInt("max Iteration for Time Stepping", &maxIter, 10.0, 100.0);
        ImGui::Checkbox("Viscosity", &sim.viscosity);
        ImGui::Checkbox("Friction", &sim.friction);
    }

    // --- New button to print breach and calls to the terminal ---
    if (ImGui::Button("Print breach and calls")) {
        std::cout << "Breach: " << breach << ", Calls: " << calls << std::endl;
    }
}

void TemplateApp::showLoggerWindow() {
    // Use a child region to embed the logger content inside the existing window.
    ImGui::BeginChild("LoggerChild", ImVec2(0, 450), true); // Adjust size as needed

    // Objective Function Plot
    if (!logger.objectiveHistory.empty()) {
        ImGui::Text("Objective Function");
        ImGui::PlotLines("##Objective", logger.objectiveHistory.data(),
                         static_cast<int>(logger.objectiveHistory.size()),
                         0, nullptr, FLT_MAX, FLT_MAX, ImVec2(0,80));
    }
    
    // Gradient Norm Plot
    if (!logger.gradientNormHistory.empty()) {
        ImGui::Text("Gradient Norm");
        ImGui::PlotLines("##Gradient", logger.gradientNormHistory.data(),
                         static_cast<int>(logger.gradientNormHistory.size()),
                         0, nullptr, FLT_MAX, FLT_MAX, ImVec2(0,80));
    }
    
    // Condition Number Plot
    // if (!logger.conditionHistory.empty()) {
    //     ImGui::Text("Hessian Condition Number");
    //     ImGui::PlotLines("##Condition", logger.conditionHistory.data(),
    //                      static_cast<int>(logger.conditionHistory.size()),
    //                      0, nullptr, FLT_MAX, FLT_MAX, ImVec2(0,80));
    // }
    
    // Latest values (optional)
    if (!logger.objectiveHistory.empty() && !logger.gradientNormHistory.empty()){ //&& !logger.conditionHistory.empty()) {
        ImGui::Text("Latest Objective: %f", logger.objectiveHistory.back());
        ImGui::Text("Latest Gradient Norm: %f", logger.gradientNormHistory.back());
        //ImGui::Text("Latest Condition Number: %f", logger.conditionHistory.back());
    }
    
    ImGui::EndChild();
}


void TemplateApp::mainLoop() {
    if (optimize) {
        if (sim.dynamic) {
            Optimization::OptimizationStatus status = energyMinimizationStepDyn();
            calls += 1;
        } else {
            Optimization::OptimizationStatus status = energyMinimizationStep();
        }
    }
    // Update the scenario object animation each frame using the simulation time step.
    sim.updateScenarioAnimation(sim.timeStep);
}


Optimization::OptimizationStatus TemplateApp::energyMinimizationStep() {
    // Get the global state vector representing all particles.
    Simulation simCopy = sim;
    globalState_0 = simCopy.getGlobalState();
    // Set up the optimization functions to work on the global state:
    optimization.objective_function = [&](const VectorXF &y, F &value) {
        simCopy.setGlobalState(y);
        simCopy.updateAuxiliaryStructures();
        simCopy.compute_energy(value);
        return true;
    };
    optimization.gradient_function = [&](const VectorXF &y, F &value, VectorXF &gradient) {
        simCopy.setGlobalState(y);
        simCopy.updateAuxiliaryStructures();
        simCopy.compute_energy(value);
        simCopy.compute_gradient(gradient);
        return true;
    };
    optimization.hessian_function = [&](const VectorXF &y, F &value, VectorXF &gradient, HessianF &hessian) {
        simCopy.setGlobalState(y);
        simCopy.updateAuxiliaryStructures();
        simCopy.compute_energy(value);
        simCopy.compute_gradient(gradient);
        simCopy.compute_hessian(hessian.A);
        return true;
    };

    // Run one optimization step on the global state.
    auto status = optimization.step(globalState_0);
    // Update the simulation with the optimized state.
    sim.setGlobalState(simCopy.getGlobalState());
    sim.renormalise();      
    sim.updateAuxiliaryStructures();
    // Log the objective, gradient norm, and Hessian condition number.
    F currentObjective;
    optimization.objective_function(globalState_0, currentObjective);
    VectorXF currentGradient;
    optimization.gradient_function(globalState_0, currentObjective, currentGradient);
    float gradNorm = currentGradient.norm();
    HessianF hessian(globalState_0.size());
    optimization.hessian_function(globalState_0, currentObjective, currentGradient, hessian);
    // MatrixXF hessianDense = hessian.evalDense();
    // float condition = computeConditionNumber(hessianDense);

    logger.logStep(static_cast<float>(currentObjective), gradNorm); //, condition);

    return status;
}

Optimization::OptimizationStatus TemplateApp::energyMinimizationStepDyn() {
    Simulation simCopy = sim;
    globalState_0 = simCopy.getGlobalState();
    Optimization::OptimizationStatus status;
    I iter = 0;
    
    F lambda = 1e0;
    F dynamicWeight = lambda / (sim.timeStep * sim.timeStep);

    optimization.objective_function = [&](const VectorXF &y, F &value) {
        simCopy.setGlobalState(y);
        simCopy.updateAuxiliaryStructures();
        // simCopy.updateEffectiveNeighborCounts();
        simCopy.compute_energy_dyn(value);
        return true;
    };

    optimization.gradient_function = [&](const VectorXF &y, F &value, VectorXF &gradient) {
        simCopy.setGlobalState(y);
        simCopy.updateAuxiliaryStructures();
        // simCopy.updateEffectiveNeighborCounts();
        simCopy.compute_energy_dyn(value);
        simCopy.compute_gradient_dyn(gradient);
        return true;
    };

    optimization.hessian_function = [&](const VectorXF &y, F &value, VectorXF &gradient, HessianF &hessian) {
        simCopy.setGlobalState(y);
        simCopy.updateAuxiliaryStructures();
        // simCopy.updateEffectiveNeighborCounts();
        simCopy.compute_energy_dyn(value);
        simCopy.compute_gradient_dyn(gradient);
        simCopy.compute_hessian_dyn(hessian.A);
        return true;
    };

    // Compute initial gradient.
    F dummyValue;
    // F diff;
    // VectorXF grad_old;
    VectorXF grad;
    optimization.gradient_function(globalState_0, dummyValue, grad);

    // Run optimization until convergence or maximum iterations are reached.
    do{
        // grad_old = grad;
        status = optimization.step(globalState_0);
        optimization.gradient_function(globalState_0, dummyValue, grad);
        // diff = (grad - grad_old).norm();
        //std::cout<<grad.norm()<<std::endl;
        iter++;
        if (iter >= maxIter) {
            std::cout << "Maximum iterations reached. Stopping simulation." << std::endl;
            // std::cout << iter << std::endl;
            optimize = false;
            break;
        }
    }while (grad.norm() > dynamic_convergence_threshold && iter < maxIter);
    
    sim.setGlobalState(simCopy.getGlobalState());
    sim.renormalise();      
    sim.updateAuxiliaryStructures();
    sim.updateEffectiveNeighborCountsFinal();
    sim.globalState_2 = sim.globalState_1;
    sim.globalState_1 = globalState_0;

    F currentObjective;
    optimization.objective_function(globalState_0, currentObjective);
    float gradNorm = grad.norm();
    HessianF hessian(globalState_0.size());
    optimization.hessian_function(globalState_0, currentObjective, grad, hessian);
    MatrixXF hessianDense = hessian.evalDense();
    //float condition = computeConditionNumber(hessianDense);

    logger.logStep(static_cast<float>(currentObjective), gradNorm);  //, condition);

    return status;
}



void TemplateApp::checkGradient(int order) {
    F lambda = 1e0;
    F dynamicWeight = lambda / (sim.timeStep * sim.timeStep);
    // Define the objective function to work on the global state.
    auto objective_function = [&](const VectorXF &y, F &value) {
        sim.setGlobalState(y);
        sim.compute_energy(value);
        if (sim.dynamic) {
            sim.compute_energy_dyn(value);
        }
        return true;
    };

    // Define the gradient function for the global state.
    auto gradient_function = [&](const VectorXF &y, VectorXF &gradient) {
        sim.setGlobalState(y);
        sim.compute_gradient(gradient);
        if (sim.dynamic) {
            sim.compute_gradient_dyn(gradient);
        }
        return true;
    };

    // Define the Hessian function.
    auto hessian_function = [&](const VectorXF &y, MatrixXF &hessian) {
        sim.setGlobalState(y);
        SparseMatrixF hessian_sparse;
        sim.compute_hessian(hessian_sparse);
        if (sim.dynamic) {
            sim.compute_hessian_dyn(hessian_sparse);
        }
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
        center(0) = sim.wrapX(center(0));   // visual copy only
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



