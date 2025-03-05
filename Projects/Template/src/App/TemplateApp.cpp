#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

#include "Projects/Template/include/App/TemplateApp.h"

#include "CRLHelper/CameraHelper.h"
#include "CRLHelper/Optimization.h"
#include "CRLHelper/CRLTimer.h"

#include <iostream>

void TemplateApp::initializeSubApp() {
    /// Initialize camera.
    app_camera.eye = Vector3F(0, 0, 5);
    app_camera.center = Vector3F(0, 0, 0);
    app_camera.set_up_direction(Vector3F(0, 1, 0));
    app_camera.height = 1;
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
}

void TemplateApp::makeConfigWindow() {
    if (ImGui::CollapsingHeader("Optimization", ImGuiTreeNodeFlags_DefaultOpen)) {
        optimization.makeConfigMenu();
        ImGui::Checkbox("Optimize", &optimize);
    }

    ImGui::Spacing();
    ImGui::Spacing();

    if (ImGui::CollapsingHeader("Model", ImGuiTreeNodeFlags_DefaultOpen)) {
        model.makeConfigMenu();
    }
}

void TemplateApp::makeAnalysisWindow() {
    if (ImGui::Button("Check Gradient")) {
        checkGradient(1);
    }
    if (ImGui::Button("Check Hessian")) {
        checkGradient(2);
    }
    ImGui::Text("Epsilon: 10^");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(100);
    ImGui::InputInt("##check_gradient_epsilon", &check_gradient_epsilon_exponent);
    ImGui::SameLine();
    ImGui::Checkbox("Print All##0", &check_gradient_print_all);
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
    /// Copy camera state.
    viewer_camera = app_camera;

    /// Initialize a single CRLViewerData struct.
    viewer_data.emplace_back();
    CRLViewerData& viewer_data_mass_spring = viewer_data.back();

    /// Render fixed spring endpoints as triangle mesh.
    viewer_data_mass_spring.mesh_v = MatrixXF::Zero(8, 3);
    viewer_data_mass_spring.mesh_f.resize(4, 3);
    viewer_data_mass_spring.mesh_c = MatrixXF::Zero(4, 3);
    F a = 0.02;
    viewer_data_mass_spring.mesh_v.row(0).head(2) = model.endpoint0 + Vector2F(-a, -a);
    viewer_data_mass_spring.mesh_v.row(1).head(2) = model.endpoint0 + Vector2F(a, -a);
    viewer_data_mass_spring.mesh_v.row(2).head(2) = model.endpoint0 + Vector2F(a, a);
    viewer_data_mass_spring.mesh_v.row(3).head(2) = model.endpoint0 + Vector2F(-a, a);
    viewer_data_mass_spring.mesh_v.row(4).head(2) = model.endpoint1 + Vector2F(-a, -a);
    viewer_data_mass_spring.mesh_v.row(5).head(2) = model.endpoint1 + Vector2F(a, -a);
    viewer_data_mass_spring.mesh_v.row(6).head(2) = model.endpoint1 + Vector2F(a, a);
    viewer_data_mass_spring.mesh_v.row(7).head(2) = model.endpoint1 + Vector2F(-a, a);
    viewer_data_mass_spring.mesh_f.row(0) << 0, 1, 2;
    viewer_data_mass_spring.mesh_f.row(1) << 0, 2, 3;
    viewer_data_mass_spring.mesh_f.row(2) << 4, 5, 6;
    viewer_data_mass_spring.mesh_f.row(3) << 4, 6, 7;

    /// Render free-floating midpoint as point cloud.
    viewer_data_mass_spring.points = Vector3F(model.y(0), model.y(1), 0).transpose();
    viewer_data_mass_spring.points_c = Vector3F(1, 0, 0).transpose();

    /// Render springs as lines.
    viewer_data_mass_spring.lines_v = MatrixXF::Zero(3, 3);
    viewer_data_mass_spring.lines_e.resize(2, 2);
    viewer_data_mass_spring.lines_c = MatrixXF::Zero(2, 3);
    viewer_data_mass_spring.lines_v.row(0).head(2) = model.endpoint0;
    viewer_data_mass_spring.lines_v.row(1).head(2) = model.y;
    viewer_data_mass_spring.lines_v.row(2).head(2) = model.endpoint1;
    viewer_data_mass_spring.lines_e.row(0) << 0, 1;
    viewer_data_mass_spring.lines_e.row(1) << 1, 2;

    /// Copy camera state.
    viewer_camera = app_camera;
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
