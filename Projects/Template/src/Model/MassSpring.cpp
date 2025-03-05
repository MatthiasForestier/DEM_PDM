#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

#include "Projects/Template/include/Model/MassSpring.h"

#include "CRLHelper/MapleHelper.h"

void MassSpring::makeConfigMenu() {
    ImGui::InputDouble("Rest Length", &rest_length, 0.05, 0.25, "%.4f");

    ImGui::Text("Endpoints:");

    ImGui::SetNextItemWidth(150);
    ImGui::InputDouble("##x0", &endpoint0(0), 0.05, 0.25, "%.4f");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(150);
    ImGui::InputDouble("##y0", &endpoint0(1), 0.05, 0.25, "%.4f");

    ImGui::SetNextItemWidth(150);
    ImGui::InputDouble("##x1", &endpoint1(0), 0.05, 0.25, "%.4f");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(150);
    ImGui::InputDouble("##y1", &endpoint1(1), 0.05, 0.25, "%.4f");
}

void MassSpring::compute_energy(F& value) const {
    VectorXF inputs(6);
    inputs << y, endpoint0, endpoint1;

    // clang-format off
    F y0 = inputs[0];
    F y1 = inputs[1];
    F e00 = inputs[2];
    F e01 = inputs[3];
    F e10 = inputs[4];
    F e11 = inputs[5];

    F t2 = pow(y0 - e00, 0.2e1);
    F t4 = pow(y1 - e01, 0.2e1);
    F t6 = sqrt(t2 + t4);
    F t8 = pow(t6 - rest_length, 0.2e1);
    F t10 = pow(y0 - e10, 0.2e1);
    F t12 = pow(y1 - e11, 0.2e1);
    F t14 = sqrt(t10 + t12);
    F t16 = pow(t14 - rest_length, 0.2e1);

    value = t8 / 0.2e1 + t16 / 0.2e1;
    // clang-format on
}

void MassSpring::compute_gradient(VectorXF& gradient) const {
    VectorXF inputs(6);
    inputs << y, endpoint0, endpoint1;

    // clang-format off
    F y0 = inputs[0];
    F y1 = inputs[1];
    F e00 = inputs[2];
    F e01 = inputs[3];
    F e10 = inputs[4];
    F e11 = inputs[5];

    F t1 = y0 - e00;
    F t2 = t1 * t1;
    F t3 = y1 - e01;
    F t4 = t3 * t3;
    F t6 = sqrt(t2 + t4);
    F t9 = 0.1e1 / t6 * (t6 - rest_length);
    F t11 = y0 - e10;
    F t12 = t11 * t11;
    F t13 = y1 - e11;
    F t14 = t13 * t13;
    F t16 = sqrt(t12 + t14);
    F t19 = 0.1e1 / t16 * (t16 - rest_length);

    F unknown[6];

    unknown[0] = t1 * t9 + t11 * t19;
    unknown[1] = t13 * t19 + t3 * t9;
    unknown[2] = -t1 * t9;
    unknown[3] = -t3 * t9;
    unknown[4] = -t11 * t19;
    unknown[5] = -t13 * t19;

    VectorXF gradient_full;
    processMapleOutput(reinterpret_cast<F *>(unknown), gradient_full, 6, 1);
    gradient = gradient_full.segment(0, 2);
    // clang-format on
}

void MassSpring::compute_hessian(SparseMatrixF& hessian) const {
    VectorXF inputs(6);
    inputs << y, endpoint0, endpoint1;

    // clang-format off
    F y0 = inputs[0];
    F y1 = inputs[1];
    F e00 = inputs[2];
    F e01 = inputs[3];
    F e10 = inputs[4];
    F e11 = inputs[5];

    F t1 = y0 - e00;
    F t2 = t1 * t1;
    F t3 = y1 - e01;
    F t4 = t3 * t3;
    F t5 = t2 + t4;
    F t6 = 0.1e1 / t5;
    F t7 = 0.4e1 * t1 * t1;
    F t9 = t7 * t6 / 0.4e1;
    F t10 = sqrt(t5);
    F t11 = t10 - rest_length;
    F t14 = 0.1e1 / t10 / t5 * t11;
    F t16 = t7 * t14 / 0.4e1;
    F t18 = 0.1e1 / t10 * t11;
    F t19 = y0 - e10;
    F t20 = t19 * t19;
    F t21 = y1 - e11;
    F t22 = t21 * t21;
    F t23 = t20 + t22;
    F t24 = 0.1e1 / t23;
    F t25 = 0.4e1 * t19 * t19;
    F t27 = t25 * t24 / 0.4e1;
    F t28 = sqrt(t23);
    F t29 = t28 - rest_length;
    F t32 = 0.1e1 / t28 / t23 * t29;
    F t34 = t25 * t32 / 0.4e1;
    F t36 = 0.1e1 / t28 * t29;
    F t41 = 0.4e1 * t3 * t1 * t14;
    F t45 = 0.4e1 * t21 * t19 * t32;
    F t46 = 0.4e1 * t1 * t3 * t6 + 0.4e1 * t19 * t21 * t24 - t41 - t45;
    F t47 = -0.2e1 * t1 * t6;
    F t53 = t1 * t47 / 0.2e1 + t1 * t1 * t14 - t18;
    F t54 = -0.2e1 * t3 * t6;
    F t57 = -0.4e1 * t3 * t1 * t14;
    F t58 = 0.2e1 * t1 * t54 - t57;
    F t59 = -0.2e1 * t19 * t24;
    F t65 = t19 * t59 / 0.2e1 + t19 * t19 * t32 - t36;
    F t66 = -0.2e1 * t21 * t24;
    F t69 = -0.4e1 * t21 * t19 * t32;
    F t70 = 0.2e1 * t19 * t66 - t69;
    F t71 = 0.4e1 * t3 * t3;
    F t73 = t71 * t6 / 0.4e1;
    F t75 = t71 * t14 / 0.4e1;
    F t76 = 0.4e1 * t21 * t21;
    F t78 = t76 * t24 / 0.4e1;
    F t80 = t76 * t32 / 0.4e1;
    F t83 = 0.2e1 * t3 * t47 - t57;
    F t89 = t3 * t54 / 0.2e1 + t3 * t3 * t14 - t18;
    F t91 = 0.2e1 * t21 * t59 - t69;
    F t97 = t21 * t66 / 0.2e1 + t21 * t21 * t32 - t36;
    F t100 = -0.2e1 * t1 * t54 - t41;
    F t104 = -0.2e1 * t19 * t66 - t45;

    F unknown[6][6];

    unknown[0][0] = t9 - t16 + t18 + t27 - t34 + t36;
    unknown[0][1] = t46 / 0.4e1;
    unknown[0][2] = t53;
    unknown[0][3] = t58 / 0.4e1;
    unknown[0][4] = t65;
    unknown[0][5] = t70 / 0.4e1;
    unknown[1][0] = t46 / 0.4e1;
    unknown[1][1] = t73 - t75 + t18 + t78 - t80 + t36;
    unknown[1][2] = t83 / 0.4e1;
    unknown[1][3] = t89;
    unknown[1][4] = t91 / 0.4e1;
    unknown[1][5] = t97;
    unknown[2][0] = t53;
    unknown[2][1] = t83 / 0.4e1;
    unknown[2][2] = t9 - t16 + t18;
    unknown[2][3] = t100 / 0.4e1;
    unknown[2][4] = 0.0e0;
    unknown[2][5] = 0.0e0;
    unknown[3][0] = t58 / 0.4e1;
    unknown[3][1] = t89;
    unknown[3][2] = t100 / 0.4e1;
    unknown[3][3] = t73 - t75 + t18;
    unknown[3][4] = 0.0e0;
    unknown[3][5] = 0.0e0;
    unknown[4][0] = t65;
    unknown[4][1] = t91 / 0.4e1;
    unknown[4][2] = 0.0e0;
    unknown[4][3] = 0.0e0;
    unknown[4][4] = t27 - t34 + t36;
    unknown[4][5] = t104 / 0.4e1;
    unknown[5][0] = t70 / 0.4e1;
    unknown[5][1] = t97;
    unknown[5][2] = 0.0e0;
    unknown[5][3] = 0.0e0;
    unknown[5][4] = t104 / 0.4e1;
    unknown[5][5] = t78 - t80 + t36;

    MatrixXF hessian_dense;
    processMapleOutput(reinterpret_cast<F *>(unknown), hessian_dense, 6, 6);
    hessian = hessian_dense.block(0, 0, 2, 2).sparseView();
    // clang-format on
}
