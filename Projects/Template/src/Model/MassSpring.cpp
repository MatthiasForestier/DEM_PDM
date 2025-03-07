#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

#include "Projects/Template/include/Model/MassSpring.h"

#include "CRLHelper/MapleHelper.h"

#include <iostream>
#include <random>   

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

void Simulation::makeConfigMenu() {
    // World parameters
    ImGui::InputDouble("Time Step", &timeStep, 1e-5, 1e-4, "%.6f");

    // Gravity (Vector3F)
    ImGui::InputDouble("Gravity", &gravity(1), 0.1, 0.5, "%.2f");

    // Number of particles
    ImGui::InputInt("Number of Particles", &numParticles);

    // Particle-related parameters
    ImGui::InputDouble("Mean", &mean, 0.01, 0.05, "%.2f");
    ImGui::InputDouble("Standard Deviation", &std, 0.01, 0.05, "%.2f");
    ImGui::InputDouble("Density", &density, 10.0, 100.0, "%.2f");
    ImGui::InputDouble("Contact Bond Normal Stiffness", &contactBondNormalStiffness, 1.0, 5.0, "%.2f");
    ImGui::InputDouble("Contact Stiffness Ratio", &contactStiffnessRatio, 0.01, 0.05, "%.2f");
    ImGui::InputDouble("Interparticle Friction", &interparticleFriction, 0.01, 0.05, "%.2f");
    ImGui::InputDouble("Contact Bond Normal Strength", &contactBondNormalStrength, 1.0, 5.0, "%.2f");
    ImGui::InputDouble("Contact Bond Shear Strength", &contactBondShearStrength, 1.0, 5.0, "%.2f");
    ImGui::InputDouble("Particle Wall Contact Normal Stiffness", &particleWallContactNormalStiffness, 1.0, 5.0, "%.2f");
    ImGui::InputDouble("Particle Wall Contact Tangential Stiffness", &particleWallContactTangentialStiffness, 1.0, 5.0, "%.2f");
    ImGui::InputDouble("Particle Wall Friction", &particleWallFriction, 0.01, 0.05, "%.2f");
    ImGui::InputDouble("Translational Damping", &translationalDamping, 0.01, 0.05, "%.2f");
    ImGui::InputDouble("Rotational Damping", &rotationalDamping, 0.01, 0.05, "%.2f");
    ImGui::InputDouble("Young's Modulus Min", &youngsModulusMin, 1e9, 1e10, "%.2e");
    ImGui::InputDouble("Young's Modulus Max", &youngsModulusMax, 1e9, 1e10, "%.2e");
    ImGui::InputDouble("Loading Velocity", &loadingVelocity, 0.01, 0.05, "%.2f");


    // Static labels for display.
    static const std::vector<const char*> poissonChoiceLabels = {"0.05", "0.15", "0.25", "0.35", "0.45"};

    // Determine the current index based on sim.poissonRatio and sim.poissonChoices.
    int currentIndex = 0;
    for (int i = 0; i < poissonChoices.size(); i++) {
        if (std::abs(poissonChoices(i) - poissonRatio) < 1e-6) {
            currentIndex = i;
            break;
        }
    }

    // Create a single combo box for Poisson Ratio.
    if (ImGui::Combo("Poisson Ratio", &currentIndex, poissonChoiceLabels.data(), static_cast<int>(poissonChoices.size()))) {
        // Update the simulation's poissonRatio based on the selected index.
        poissonRatio = poissonChoices(currentIndex);
    }

    if (ImGui::Button("↻ Reload Particles")) {
        // Here you would call your particle creation function.
        // For example, if you have a function createRandomParticles2D() accessible from this context:
        particles2D = createRandomParticles2D();
    }

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

Particle2D::Particle2D(F radius, const Simulation& simParams)
    : pos(Vector3F::Zero()), vel(Vector3F::Zero()), acc(Vector3F::Zero()), radius(radius)
{
    // For a 2D disc, mass = area * density.
    mass = M_PI * radius * radius * simParams.density;
    // Moment of inertia for a uniform disc about its center: I = 1/2 * m * r^2.
    inertia = 0.5 * mass * radius * radius;
}

Particle3D::Particle3D(F radius, const Simulation& sim)
    : pos(Vector6F::Zero()), vel(Vector6F::Zero()), acc(Vector6F::Zero()), radius(radius)
{
    // For a sphere, mass = volume * density.
    mass = (4.0 / 3.0) * M_PI * std::pow(radius, 3) * sim.density;
    // Moment of inertia for a solid sphere: I = 2/5 * m * r^2.
    inertia = (2.0 / 5.0) * mass * radius * radius;
}

// Function to create a specified number of Particle2D objects with random positions and radii.
std::vector<Particle2D> Simulation::createRandomParticles2D() {
    std::vector<Particle2D> particles;
    particles.reserve(numParticles);

    // Set up random number generators.
    std::random_device rd;
    std::mt19937 gen(rd());

    // Uniform distributions for x and y coordinates within the display.
    std::uniform_real_distribution<F> distX(-1.0, 1.0);
    std::uniform_real_distribution<F> distY(-0.5, 0.5);

    // Normal distribution for disk radii.
    std::normal_distribution<F> radiusDist(mean, std); // mean 0.05, std 0.01

    for (int i = 0; i < numParticles; i++) {
        bool validCandidate = false;
        int attempts = 0;
        Particle2D candidate(0.0, *this); // Use *this instead of 'sim'

        // Try until we find a candidate that doesn't overlap or reach the maximum attempts.
        while (!validCandidate && attempts < maxAttemptsPerParticle) {
            attempts++;

            // Generate a random radius and ensure it's positive.
            F radius = radiusDist(gen);
            if (radius <= 0)
                radius = 0.01;

            // Create a candidate particle with the new radius.
            candidate = Particle2D(radius, *this);
            candidate.pos = Vector3F(distX(gen), distY(gen), 0.0);

            // Check for overlap with all previously accepted particles.
            validCandidate = true;
            for (const auto &existing : particles) {
                // Only the x and y coordinates are considered.
                F distance = (candidate.pos.head(2) - existing.pos.head(2)).norm();
                if (distance < (candidate.radius + existing.radius)) {
                    validCandidate = false;
                    break;
                }
            }
        }

        if (validCandidate) {
            particles.push_back(candidate);
        } else {
            std::cerr << "Warning: Could not place particle " << i + 1 << " without overlapping after "
                      << maxAttemptsPerParticle << " attempts.\n";
            // Optionally, break out of the loop if placement becomes too difficult.
        }
    }

    return particles;
}