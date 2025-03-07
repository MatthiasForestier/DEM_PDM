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
    app_camera.height = 1;
    //Initialize particles
    sim.particles2D = sim.createRandomParticles2D();
    createOrUpdateSolver();
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

Optimization::OptimizationStatus TemplateApp::energyMinimizationStep() {
    optimization.objective_function = [&](const VectorXF& y, F& value) {
        model.y = y;
        model.compute_energy(value);
        return true;
    };
    optimization.gradient_function = [&](const VectorXF& y, F& value, VectorXF& gradient) {
        model.y = y;
        model.compute_energy(value);
        model.compute_gradient(gradient);
        return true;
    };
    optimization.hessian_function = [&](const VectorXF& y, F& value, VectorXF& gradient, HessianF& hessian) {
        model.y = y;
        model.compute_energy(value);
        model.compute_gradient(gradient);
        model.compute_hessian(hessian.A);
        return true;
    };

    VectorXF opt_y = model.y;
    Optimization::OptimizationStatus status = optimization.step(opt_y);
    model.y = opt_y;

    return status;
}

void TemplateApp::mainLoop() {
    if (optimize) {
        Optimization::OptimizationStatus status = energyMinimizationStep();
    }
    if (runSimulation && solver) {
        solver->step(sim.timeStep);
    }
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

    // --- New solver selection controls ---
    static const char* solverNames[] = {"Forward Euler", "Backward Euler"};
    static int currentSolverIndex = 0;
    if (ImGui::Combo("Solver Type", &currentSolverIndex, solverNames, IM_ARRAYSIZE(solverNames))) {
        // Update solverType based on user choice
        solverType = solverNames[currentSolverIndex];
        // Recreate or update the solver if needed
        createOrUpdateSolver();
    }

    // --- Run / Pause controls ---
    if (ImGui::Button("Run Simulation")) {
        runSimulation = true;  // or any other boolean you use to indicate "simulation running"
    }
    ImGui::SameLine();
    if (ImGui::Button("Pause Simulation")) {
        runSimulation = false;
    }
}

void TemplateApp::checkGradient(int order) {
    auto objective_function = [&](const VectorXF& y, F& value) {
        model.y = y;
        model.compute_energy(value);
        return true;
    };
    auto gradient_function = [&](const VectorXF& y, VectorXF& gradient) {
        model.y = y;
        model.compute_gradient(gradient);
        return true;
    };
    auto hessian_function = [&](const VectorXF& y, MatrixXF& hessian) {
        model.y = y;
        SparseMatrixF hessian_sparse;
        model.compute_hessian(hessian_sparse);
        hessian = hessian_sparse.toDense();
        return true;
    };

    VectorXF curr_y = model.y;
    F epsilon = pow(10.0, check_gradient_epsilon_exponent);
    int print_level = (check_gradient_print_all ? 2 : 1);
    if (order == 1) {
        VectorXF error;
        Optimization::checkGradient(curr_y, error, objective_function, gradient_function, epsilon, print_level);
    }
    if (order == 2) {
        MatrixXF error;
        Optimization::checkHessian(curr_y, error, gradient_function, hessian_function, epsilon, print_level);
    }
}

void TemplateApp::getViewerData(std::vector<CRLViewerData>& viewer_data, CRLCamera& viewer_camera) {
    // Copy camera state.
    viewer_camera = app_camera;

    // Create a new CRLViewerData for particle disks.
    viewer_data.emplace_back();
    CRLViewerData& particleViewerData = viewer_data.back();
    
    // Clear any existing mesh data.
    particleViewerData.mesh_v.resize(0, 3);
    particleViewerData.mesh_f.resize(0, 3);
    
    // Variables to keep track of global indices.
    int vertexOffset = 0;
    
    int numSegments = 32;
    
    for (size_t i = 0; i < sim.particles2D.size(); i++) {
        const Particle2D& p = sim.particles2D[i];
        
        // Create disk mesh for this particle.
        std::vector<Vector3F> diskVertices;
        std::vector<Vector3I> diskFaces;
        Vector3F center = p.pos;
        F radius = p.radius;
        
        // Center vertex.
        diskVertices.push_back(center);
        
        // Perimeter vertices.
        for (int j = 0; j < numSegments; j++) {
            F theta = 2.0 * M_PI * j / numSegments;
            diskVertices.push_back(center + Vector3F(radius * cos(theta), radius * sin(theta), 0.0));
        }
        
        // Create faces using triangle fan.
        for (int j = 1; j < numSegments; j++) {
            diskFaces.push_back(Vector3I(0, j, j + 1));
        }
        diskFaces.push_back(Vector3I(0, numSegments, 1));
        
        // Append diskVertices into the global mesh.
        for (const auto& v : diskVertices) {
            particleViewerData.mesh_v.conservativeResize(particleViewerData.mesh_v.rows() + 1, 3);
            particleViewerData.mesh_v.row(particleViewerData.mesh_v.rows() - 1) = v.transpose();
        }
        
        // Append diskFaces into the global face array.
        for (const auto& f : diskFaces) {
            particleViewerData.mesh_f.conservativeResize(particleViewerData.mesh_f.rows() + 1, 3);
            // Adjust face indices by vertexOffset.
            particleViewerData.mesh_f.row(particleViewerData.mesh_f.rows() - 1) = (f.array() + vertexOffset).matrix().transpose();
        }
        
        // Update the vertex offset.
        vertexOffset += diskVertices.size();
    }
    
    // Optionally set a color for the particle disks.
    // For example, you could fill particleViewerData.mesh_c with a uniform color.
    particleViewerData.mesh_c = MatrixXF::Ones(particleViewerData.mesh_f.rows(), 3); // white color
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



