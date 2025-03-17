#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

#include "Projects/Template/include/Model/MassSpring.h"

#include "CRLHelper/MapleHelper.h"

#include <iostream>
#include <random>   
#include <cmath>
#include <algorithm>
#include <unordered_set>


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
    ImGui::InputDouble("Boundary Overlap", &overlapParam, 10.0, 100.0, "%.2f");
    ImGui::InputDouble("Interaction Overlap", &interactionParam, 10.0, 100.0, "%.2f");
    // ImGui::InputDouble("Contact Bond Normal Stiffness", &contactBondNormalStiffness, 1.0, 5.0, "%.2f");
    // ImGui::InputDouble("Contact Stiffness Ratio", &contactStiffnessRatio, 0.01, 0.05, "%.2f");
    // ImGui::InputDouble("Interparticle Friction", &interparticleFriction, 0.01, 0.05, "%.2f");
    // ImGui::InputDouble("Contact Bond Normal Strength", &contactBondNormalStrength, 1.0, 5.0, "%.2f");
    // ImGui::InputDouble("Contact Bond Shear Strength", &contactBondShearStrength, 1.0, 5.0, "%.2f");
    // ImGui::InputDouble("Particle Wall Contact Normal Stiffness", &particleWallContactNormalStiffness, 1.0, 5.0, "%.2f");
    // ImGui::InputDouble("Particle Wall Contact Tangential Stiffness", &particleWallContactTangentialStiffness, 1.0, 5.0, "%.2f");
    // ImGui::InputDouble("Particle Wall Friction", &particleWallFriction, 0.01, 0.05, "%.2f");
    // ImGui::InputDouble("Translational Damping", &translationalDamping, 0.01, 0.05, "%.2f");
    // ImGui::InputDouble("Rotational Damping", &rotationalDamping, 0.01, 0.05, "%.2f");
    // ImGui::InputDouble("Young's Modulus Min", &youngsModulusMin, 1e9, 1e10, "%.2e");
    // ImGui::InputDouble("Young's Modulus Max", &youngsModulusMax, 1e9, 1e10, "%.2e");
    // ImGui::InputDouble("Loading Velocity", &loadingVelocity, 0.01, 0.05, "%.2f");

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
    if (ImGui::Combo("Poisson Ratio", &currentIndex,
                     poissonChoiceLabels.data(),
                     static_cast<int>(poissonChoices.size()))) {
        // Update the simulation's poissonRatio based on the selected index.
        poissonRatio = poissonChoices(currentIndex);
    }

    if (ImGui::Button("↻ Reload Particles")) {
        // Here you would call your particle creation function.
        // For example, if you have a function createRandomParticles2D() accessible from this context:
        particles2D = createRandomParticles2D();
        insertParticlesIntoGrid();
    }

    // --- New UI for Scenario Objects ---

    // Toggle flag for activation.
    static bool activateScenario = false;
    static bool lastFrameToggleState = false;

    // Provide a few shape choices.
    static int shapeIndex = 0;
    const char* shapes[] = {"Circle", "Square"};

    // ImGui Checkbox for toggling on/off.
    if (ImGui::Checkbox("Activate Scenario Object", &activateScenario)) {
        // If the user has just toggled OFF, remove the scenario object(s) from the simulation.
        if (!activateScenario && lastFrameToggleState) {
            scenarioObjects.clear();
            std::cout << "Scenario object(s) removed from simulation.\n";
        }
        lastFrameToggleState = activateScenario;
    }

    // If activated, show a combo for shape selection + create button.
    if (activateScenario) {
        ImGui::Combo("Shape Choice", &shapeIndex, shapes, IM_ARRAYSIZE(shapes));
        
        if (ImGui::Button("Create Scenario")) {
            // Clear any previously stored scenario objects.
            scenarioObjects.clear();
            
            if (shapeIndex == 0) {
                // Create a Circle scenario object with 64 segments, radius 1.1, centered at the origin.
                scenarioObjects.push_back(std::make_unique<Circle>(1.1f, 64, Vector3F(0.0f, 0.0f, 0.0f)));
            } else {
                // Create a Square scenario object with half-dimensions 1.1 and 0.9, centered at the origin.
                scenarioObjects.push_back(std::make_unique<Square>(1.1f, 0.9f, Vector3F(0.0f, 0.0f, 0.0f)));
            }
            buildGridDataStructure();
            insertParticlesIntoGrid();
            std::cout << "New scenario object created and stored!\n";
        }
    }
}

void Simulation::buildGridDataStructure() 
{
    if (scenarioObjects.empty()) {
    std::cerr << "No scenario objects defined in simulation.\n";
    return;
    }
    // Use the cached bounding box (BB) from the first scenario object.
    const BoundingBox &bbox = scenarioObjects[0]->getBoundingBox();
    ;
    minX = bbox.min_x;
    maxX = bbox.max_x;
    minY = bbox.min_y;
    maxY = bbox.max_y;

    // Compute the number of cells along each axis.
    numCellsX = static_cast<int>(std::ceil((maxX - minX) / cellSize));
    numCellsY = static_cast<int>(std::ceil((maxY - minY) / cellSize));

    // Resize the grid vector to hold all cells.
    grid.clear();
    grid.resize(numCellsX * numCellsY);
}

// Updated: Insert Particle2D objects into the grid.
// Now each particle is inserted into all cells that its circle overlaps.
void Simulation::insertParticlesIntoGrid()
{
    // Clear any previous data in the grid.
    for (auto &cell : grid) {
        cell.clear();
    }

    // Loop over all Particle2D objects.
    for (size_t i = 0; i < particles2D.size(); i++) {
        const Particle2D &p = particles2D[i];

        // Compute the extent of the particle (its bounding box).
        F xMinParticle = p.pos(0) - p.radius;
        F xMaxParticle = p.pos(0) + p.radius;
        F yMinParticle = p.pos(1) - p.radius;
        F yMaxParticle = p.pos(1) + p.radius;

        // Compute the grid cell range that covers this bounding box.
        int minCellX = static_cast<int>(std::floor((xMinParticle - minX) / cellSize));
        int maxCellX = static_cast<int>(std::floor((xMaxParticle - minX) / cellSize));
        int minCellY = static_cast<int>(std::floor((yMinParticle - minY) / cellSize));
        int maxCellY = static_cast<int>(std::floor((yMaxParticle - minY) / cellSize));

        // Clamp indices to grid bounds.
        minCellX = std::max(minCellX, 0);
        maxCellX = std::min(maxCellX, numCellsX - 1);
        minCellY = std::max(minCellY, 0);
        maxCellY = std::min(maxCellY, numCellsY - 1);

        // Insert the particle into each overlapping cell.
        for (int cx = minCellX; cx <= maxCellX; cx++) {
            for (int cy = minCellY; cy <= maxCellY; cy++) {
                int index = cy * numCellsX + cx;
                grid[index].push_back(static_cast<int>(i));
            }
        }
    }
}

void Simulation::updateNeighborLists()
{
    // Clear neighbor lists for all particles.
    for (auto &p : particles2D) {
        p.neighborIndices.clear();
    }

    // Loop over each particle.
    for (size_t i = 0; i < particles2D.size(); i++) {
        Particle2D &p = particles2D[i];

        // Compute the bounding cell range for particle p (using its full circle).
        int minCellX = static_cast<int>(std::floor((p.pos(0) - p.radius - minX) / cellSize));
        int maxCellX = static_cast<int>(std::floor((p.pos(0) + p.radius - minX) / cellSize));
        int minCellY = static_cast<int>(std::floor((p.pos(1) - p.radius - minY) / cellSize));
        int maxCellY = static_cast<int>(std::floor((p.pos(1) + p.radius - minY) / cellSize));

        // Clamp indices.
        minCellX = std::max(minCellX, 0);
        maxCellX = std::min(maxCellX, numCellsX - 1);
        minCellY = std::max(minCellY, 0);
        maxCellY = std::min(maxCellY, numCellsY - 1);

        // Use a set to avoid duplicate candidate neighbors.
        std::unordered_set<int> neighborCandidates;
        for (int cx = minCellX; cx <= maxCellX; cx++) {
            for (int cy = minCellY; cy <= maxCellY; cy++) {
                int gridIndex = cy * numCellsX + cx;
                for (int j : grid[gridIndex]) {
                    if (j != static_cast<int>(i)) { // avoid self
                        neighborCandidates.insert(j);
                    }
                }
            }
        }

        // Check for actual overlap with each candidate.
        for (int j : neighborCandidates) {
            // To avoid duplicate symmetric entries, you can choose to only add if j > i.
            if (j <= static_cast<int>(i)) continue;
            const Particle2D &q = particles2D[j];
            F dx = p.pos(0) - q.pos(0);
            F dy = p.pos(1) - q.pos(1);
            F distance = std::sqrt(dx * dx + dy * dy);
            if (distance < (p.radius + q.radius)) {
                p.neighborIndices.push_back(j);
            }
        }
    }
}

void Simulation::updateGlobalPositions() {
    int n = 0;
    if (use3D) {
        n = static_cast<int>(particles3D.size());
        globalPositions.resize(6 * n);
        for (int i = 0; i < n; i++) {
            globalPositions.segment<6>(6*i) = particles3D[i].pos;
        }
    }else{
        n = static_cast<int>(particles2D.size());
        globalPositions.resize(3 * n);
        for (int i = 0; i < n; i++) {
            globalPositions.segment<3>(3*i) = particles2D[i].pos;
        }
    }
}

// Similarly, when updating particles from the global state:
void Simulation::applyGlobalPositions(const VectorXF &positions) {
    if (use3D) {
        int n = static_cast<int>(particles3D.size());
        for (int i = 0; i < n; i++) {
            particles3D[i].pos = positions.segment<6>(6*i);
        }
    } else {
        int n = static_cast<int>(particles2D.size());
        for (int i = 0; i < n; i++) {
            particles2D[i].pos = positions.segment<3>(3*i);
        }
    }
}

VectorXF Simulation::getGlobalState() {
    updateGlobalPositions();
    return globalPositions;
}

void Simulation::setGlobalState(const VectorXF &state) {
    globalPositions = state;
    applyGlobalPositions(globalPositions);
}

void Simulation::updateAuxiliaryStructures() {
    // Build the grid and update neighbor lists based on the current particle positions.
    insertParticlesIntoGrid();
    updateNeighborLists();
}

void Simulation::compute_energy(F &value) const {
    value = 0;
    // First, add gravitational and boundary collision energy contributions.
    for (size_t i = 0; i < particles2D.size(); i++) {
        // Copy particle so that detectBoundaryCollision2D can update its state.
        Particle2D p = particles2D[i];
        detectBoundaryCollision2D(p);

        // Gravitational potential energy: m * g * y.
        value += gravity(1) * p.pos(1); //p.mass *

        // Energy contributions due to boundary collisions.
        if (p.BoundaryCollision == 1) {
            for (const auto &so : scenarioObjects) {
                if (Circle* circle = dynamic_cast<Circle*>(so.get())) {
                    F px = p.pos(0);
                    F py = p.pos(1);
                    F radius = p.radius;
                    F Rb = circle->radius;
                    F so_x = so->position(0);
                    F so_y = so->position(1);
                    // Distance between particle and circle center.
                    F dx = px - so_x;
                    F dy = py - so_y;
                    F dist = std::sqrt(dx * dx + dy * dy);
                    // Penetration offset (here simplified as the particle radius minus the circle radius)
                    F delta = radius - Rb;
                    // Penalty energy contribution.
                    value += 0.5 * overlapParam * (dist + delta) * (dist + delta);
                }
            }
        } else if (p.BoundaryCollision == 2) {
            for (const auto &so : scenarioObjects) {
                if (Square* square = dynamic_cast<Square*>(so.get())) {
                    F px = p.pos(0);
                    F py = p.pos(1);
                    F radius = p.radius;

                    F min_x = square->min_x;
                    F max_x = square->max_x;
                    F min_y = square->min_y;
                    F max_y = square->max_y;

                    F overlap_x_r = 0, overlap_x_l = 0;
                    if (px - radius < min_x) {
                        overlap_x_l = min_x - (px - radius);
                    } else if (px + radius > max_x) {
                        overlap_x_r = (px + radius) - max_x;
                    }
                    F overlap_y_b = 0, overlap_y_t = 0;
                    if (py - radius < min_y) {
                        overlap_y_b = min_y - (py - radius);
                    } else if (py + radius > max_y) {
                        overlap_y_t = (py + radius) - max_y;
                    }
                    // Total overlap computed as the Euclidean norm of the directional overlaps.
                    F total_overlap = std::sqrt(overlap_x_l * overlap_x_l +
                                                overlap_y_t * overlap_y_t +
                                                overlap_x_r * overlap_x_r +
                                                overlap_y_b * overlap_y_b);
                    if (total_overlap > 0) {
                        value += 0.5 * overlapParam * total_overlap * total_overlap;
                    }
                }
            }
        }
    }

    // --- Interparticle overlap energy ---
    // Loop over each particle and its neighbor list.
    for (size_t i = 0; i < particles2D.size(); i++) {
        const Particle2D &p = particles2D[i];
        for (int j : p.neighborIndices) {
            // Each neighbor index j is guaranteed to be > i (avoid duplicate work).
            const Particle2D &q = particles2D[j];
            F dx = p.pos(0) - q.pos(0);
            F dy = p.pos(1) - q.pos(1);
            F d = std::sqrt(dx * dx + dy * dy);
            // Compute overlap only if particles are close.
            F overlap = (p.radius + q.radius) - d;
            if (overlap > 0) {
                value += 0.5 * interactionParam * overlap * overlap;
            }
        }
    }
}

void Simulation::compute_gradient(VectorXF &gradient) const {
    // Global state has 3 entries per particle: [x, y, theta]
    gradient.resize(3 * particles2D.size());
    gradient.setZero();

    // First, add gravitational and boundary collision gradient contributions.
    for (size_t i = 0; i < particles2D.size(); i++) {
        Particle2D p = particles2D[i];
        detectBoundaryCollision2D(p);

        // Gravitational gradient: only affects the y-component.
        gradient(3 * i)     = 0;
        gradient(3 * i + 1) = gravity(1); //p.mass *
        gradient(3 * i + 2) = 0;

        if (p.BoundaryCollision == 1) {
            for (const auto &so : scenarioObjects) {
                if (Circle* circle = dynamic_cast<Circle*>(so.get())) {
                    F px = p.pos(0);
                    F py = p.pos(1);
                    F boundaryRadius = circle->radius;
                    F so_x = so->position(0);
                    F so_y = so->position(1);
                    
                    // Penetration depth.
                    F penetration = (p.radius - boundaryRadius);
                    F stiffness = overlapParam;
                    
                    F dx = px - so_x;
                    F dy = py - so_y;
                    F dist = std::sqrt(dx * dx + dy * dy);
                    if (dist == 0) continue;
                    // The gradient contribution (simplified).
                    F factor = stiffness * penetration / dist;
                    gradient(3 * i)     += dx * factor + dx * stiffness;
                    gradient(3 * i + 1) += dy * factor + dy * stiffness;
                }
            }
        } else if (p.BoundaryCollision == 2) {
            for (const auto &so : scenarioObjects) {
                if (Square* square = dynamic_cast<Square*>(so.get())) {
                    F px = p.pos(0);
                    F py = p.pos(1);
                    F radius = p.radius;

                    F min_x = square->min_x;
                    F max_x = square->max_x;
                    F min_y = square->min_y;
                    F max_y = square->max_y;

                    F overlap_x_r = 0, overlap_x_l = 0;
                    if (px - radius < min_x) {
                        overlap_x_l = min_x - (px - radius);
                    } else if (px + radius > max_x) {
                        overlap_x_r = (px + radius) - max_x;
                    }
                    F overlap_y_b = 0, overlap_y_t = 0;
                    if (py - radius < min_y) {
                        overlap_y_b = min_y - (py - radius);
                    } else if (py + radius > max_y) {
                        overlap_y_t = (py + radius) - max_y;
                    }
                    F total_overlap = std::sqrt(overlap_x_l * overlap_x_l +
                                                overlap_y_t * overlap_y_t +
                                                overlap_x_r * overlap_x_r +
                                                overlap_y_b * overlap_y_b);
                    if (total_overlap > 0) {
                        if (overlap_x_l > 0) {
                            gradient(3 * i) += overlapParam * (-min_x + px - radius);
                        } else if (overlap_x_r > 0) {
                            gradient(3 * i) += overlapParam * (px - max_x + radius);
                        }
                        if (overlap_y_b > 0) {
                            gradient(3 * i + 1) += overlapParam * (-min_y + py - radius);
                        } else if (overlap_y_t > 0) {
                            gradient(3 * i + 1) +=  overlapParam * (py - max_y + radius);
                        }
                    }
                }
            }
        }
    }

    // --- Interparticle gradient contributions ---
    // Loop over each particle and its neighbor list.
    for (size_t i = 0; i < particles2D.size(); i++) {
        const Particle2D &p = particles2D[i];
        for (int j : p.neighborIndices) {
            const Particle2D &q = particles2D[j];
            F dx = p.pos(0) - q.pos(0);
            F dy = p.pos(1) - q.pos(1);
            F d = std::sqrt(dx * dx + dy * dy);
            const F eps = 1e-6; // or another appropriate threshold
            if (d < eps) continue; // avoid division by zero
            F overlap = (p.radius + q.radius) - d;
            if (overlap > 0) {
                // Gradient contribution: ∇E = -k_overlap * Δ * ( (dx, dy)/d ).
                F factor = interactionParam * overlap / d;
                // Update gradient for particle i.
                gradient(3 * i)     += -factor * dx;
                gradient(3 * i + 1) += -factor * dy;
                // Update gradient for particle j (opposite sign).
                gradient(3 * j)     += factor * dx;
                gradient(3 * j + 1) += factor * dy;
            }
        }
    }
}

void Simulation::compute_hessian(SparseMatrixF &hessian) const {
    int n = 3 * particles2D.size();
    hessian.resize(n, n);
    hessian.setZero();

    // Boundary collision Hessian contributions (existing code remains unchanged).
    for (size_t i = 0; i < particles2D.size(); i++) {
        Particle2D p = particles2D[i];
        detectBoundaryCollision2D(p);
        if (p.BoundaryCollision == 1) {
            for (const auto &so : scenarioObjects) {
                if (Circle* circle = dynamic_cast<Circle*>(so.get())) {
                    F px = p.pos(0), py = p.pos(1);
                    F boundaryRadius = circle->radius;
                    F so_x = so->position(0), so_y = so->position(1);
                    F t3 = (p.radius - boundaryRadius) * overlapParam;
                    F t4 = px * px;
                    F t7 = py * py;
                    F t10 = so_x * so_x;
                    F t11 = so_y * so_y;
                    F t12 = t4 - 2.0 * px * so_x + t7 - 2.0 * py * so_y + t10 + t11;
                    F t13 = std::sqrt(t12);
                    F t15 = 1.0 / (t13 * t12);
                    F t18 = 2.0 * px - 2.0 * so_x;
                    F t19 = t18 * t18;
                    F t24 = (1.0 / t13) * t3;
                    F t25 = overlapParam;
                    F t30 = 2.0 * py - 2.0 * so_y;
                    F t33 = (t30 * t18 * t15 * t3) / 4.0;
                    F t34 = t30 * t30;

                    hessian.coeffRef(3 * i, 3 * i)       = -t19 * t15 * t3 / 4.0 + t24 + t25; // 4.0 
                    hessian.coeffRef(3 * i, 3 * i + 1)   = -t33;
                    hessian.coeffRef(3 * i, 3 * i + 2)   = 0.0;
                    hessian.coeffRef(3 * i + 1, 3 * i)   = -t33;
                    hessian.coeffRef(3 * i + 1, 3 * i + 1) = -t34 * t15 * t3 / 4.0  + t24 + t25; // 4.0 
                    hessian.coeffRef(3 * i + 1, 3 * i + 2) = 0.0;
                    hessian.coeffRef(3 * i + 2, 3 * i)   = 0.0;
                    hessian.coeffRef(3 * i + 2, 3 * i + 1) = 0.0;
                    hessian.coeffRef(3 * i + 2, 3 * i + 2) = 0.0;
                }
            }
        } else if (p.BoundaryCollision == 2) {
            for (const auto &so : scenarioObjects) {
                if (Square* square = dynamic_cast<Square*>(so.get())) {
                    F px = p.pos(0), py = p.pos(1);
                    F min_x = square->min_x, max_x = square->max_x;
                    F min_y = square->min_y, max_y = square->max_y;
                    F overlap_x_r = 0, overlap_x_l = 0;
                    if (px - p.radius < min_x) {
                        overlap_x_l = min_x - (px - p.radius);
                    } else if (px + p.radius > max_x) {
                        overlap_x_r = (px + p.radius) - max_x;
                    }
                    F overlap_y_b = 0, overlap_y_t = 0;
                    if (py - p.radius < min_y) {
                        overlap_y_b = min_y - (py - p.radius);
                    } else if (py + p.radius > max_y) {
                        overlap_y_t = (py + p.radius) - max_y;
                    }
                    F total_overlap = std::sqrt(overlap_x_l * overlap_x_l +
                                               overlap_y_t * overlap_y_t +
                                               overlap_x_r * overlap_x_r +
                                               overlap_y_b * overlap_y_b);
                    F t1 = overlapParam;
                    if (total_overlap > 0) {
                        if (overlap_x_l > 0 || overlap_x_r > 0){
                            hessian.coeffRef(3 * i, 3 * i)       = t1;
                        }
                        hessian.coeffRef(3 * i, 3 * i + 1)   = 0.0;
                        hessian.coeffRef(3 * i, 3 * i + 2)   = 0.0;
                        hessian.coeffRef(3 * i + 1, 3 * i)   = 0.0;
                        if (overlap_y_t > 0 || overlap_y_b > 0){
                            hessian.coeffRef(3 * i + 1, 3 * i + 1) = t1;
                        }
                        hessian.coeffRef(3 * i + 1, 3 * i + 2) = 0.0;
                        hessian.coeffRef(3 * i + 2, 3 * i)   = 0.0;
                        hessian.coeffRef(3 * i + 2, 3 * i + 1) = 0.0;
                        hessian.coeffRef(3 * i + 2, 3 * i + 2) = 0.0;
                    }
                }
            }
        }
    }

    //--- Interparticle Hessian contributions ---
    //Loop over each particle and its neighbor list.
    for (size_t i = 0; i < particles2D.size(); i++) {
        const Particle2D &p = particles2D[i];
        for (int j : p.neighborIndices) {
            const Particle2D &q = particles2D[j];
            F dx = p.pos(0) - q.pos(0);
            F dy = p.pos(1) - q.pos(1);
            F d = std::sqrt(dx * dx + dy * dy);
            const F eps = 1e-6; // or another appropriate threshold
            if (d < eps) continue; // avoid division by zero
            F overlap = (p.radius + q.radius) - d;
            if (overlap > 0) {
                // Unit vector from q to p.
                Eigen::Vector2d u(dx / d, dy / d);
                // Compute the 2x2 Hessian block.
                Eigen::Matrix2d H_block = interactionParam * (u * u.transpose()) -
                    (interactionParam * overlap / d) * (Eigen::Matrix2d::Identity() - u * u.transpose());

                // Update the Hessian blocks for particles i and j (only for x and y components).
                for (int a = 0; a < 2; a++) {
                    for (int b = 0; b < 2; b++) {
                        hessian.coeffRef(3 * i + a, 3 * i + b) += H_block(a, b);
                        hessian.coeffRef(3 * j + a, 3 * j + b) += H_block(a, b);
                        hessian.coeffRef(3 * i + a, 3 * j + b) -= H_block(a, b);
                        hessian.coeffRef(3 * j + a, 3 * i + b) -= H_block(a, b);
                    }
                }
            }
        }
    }
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

std::vector<Particle2D> Simulation::createRandomParticles2D() {
    std::vector<Particle2D> particles;
    particles.reserve(numParticles);

    // Set up random number generators.
    std::random_device rd;
    std::mt19937 gen(rd());

    // Uniform distributions for x and y coordinates within the display.
    std::uniform_real_distribution<F> distX(-0.8, 0.8);
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

void Simulation::colorParticleRed(int particleID) {
    if (particleID < 0 || particleID >= static_cast<int>(particles2D.size())) {
        std::cerr << "Error: Invalid particle ID: " << particleID << std::endl;
        return;
    }
    // Set the particle's color to red (RGB: 1, 0, 0)
    particles2D[particleID].color = Color(1.0f, 0.0f, 0.0f);
}


void Simulation::detectBoundaryCollision2D(Particle2D &p) const {
    p.BoundaryCollision = 0;
    // Loop over each scenario object and use the virtual function.
    for (const auto &so : scenarioObjects) {
        int collision = so->detectCollision(p);
        if (collision != 0) {
            p.BoundaryCollision = collision;
            // Optionally break out if you only need to detect the first collision.
            break;
        }
    }
}

Square::Square(F halfLength, F halfWidth, const Vector3F& pos)
    : halfLength(halfLength), halfWidth(halfWidth)
{
    position = pos;
    generateVertices();
}

void Square::generateVertices() {
    vertices.clear();
    // Compute bounds.
    min_x = position(0) - halfLength;
    max_x = position(0) + halfLength;
    min_y = position(1) - halfWidth;
    max_y = position(1) + halfWidth;
    
    // Update the cached bounding box.
    BB.min_x = min_x;
    BB.max_x = max_x;
    BB.min_y = min_y;
    BB.max_y = max_y;
    
    // Compute the vertices.
    vertices.push_back(Vector3F(min_x, min_y, position(2)));
    vertices.push_back(Vector3F(min_x, max_y, position(2)));
    vertices.push_back(Vector3F(max_x, max_y, position(2)));
    vertices.push_back(Vector3F(max_x, min_y, position(2)));
}

int Square::detectCollision(const Particle2D &p) const {
    // Use the cached bounding box (BB) for collision detection.
    if ((p.pos(0) - p.radius) < BB.min_x || (p.pos(0) + p.radius) > BB.max_x ||
        (p.pos(1) - p.radius) < BB.min_y || (p.pos(1) + p.radius) > BB.max_y)
    {
        return 2; // Collision with square boundary.
    }
    return 0;
}

// Implementation for Circle.

Circle::Circle(F radius, I numSegments, const Vector3F& pos)
    : radius(radius), numSegments(numSegments)
{
    position = pos;
    generateVertices();
}

void Circle::generateVertices() {
    vertices.clear();
    // Generate vertices approximating the circle.
    for (I i = 0; i < numSegments; i++) {
        F theta = 2.0f * F(M_PI) * F(i) / F(numSegments);
        F x = radius * std::cos(theta);
        F y = radius * std::sin(theta);
        vertices.push_back(Vector3F(position(0) + x, position(1) + y, position(2)));
    }
    
    // Update the bounding box for the circle.
    BB.min_x = position(0) - radius;
    BB.max_x = position(0) + radius;
    BB.min_y = position(1) - radius;
    BB.max_y = position(1) + radius;
}

int Circle::detectCollision(const Particle2D &p) const {
    // Use the circle's center and radius to detect collision.
    F dx = p.pos(0) - position(0);
    F dy = p.pos(1) - position(1);
    F distSq = dx * dx + dy * dy;
    F limit = radius - p.radius;
    if (distSq > (limit * limit))
        return 1; // Collision with circle boundary.
    return 0;
}