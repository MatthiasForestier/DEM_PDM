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
    //Circular particles
    
    ImGui::InputDouble("Mean", &mean, 0.01, 0.05, "%.2f");
    ImGui::InputDouble("Standard Deviation", &std, 0.01, 0.05, "%.4f");

    //Polygonal particles

    //Ellipsoidal particles


    // Physics-related parameters
    ImGui::InputDouble("Density", &density, 10.0, 100.0, "%.2f");
    ImGui::InputDouble("Boundary Overlap", &overlapParam, 10.0, 100.0, "%.2f");
    ImGui::InputDouble("Interaction Overlap", &interactionParam, 10.0, 100.0, "%.2f");
    ImGui::InputDouble("Viscosity Coefficient", &viscosity_coeff, 0.01, 0.05, "%.3f");
    ImGui::InputDouble("Sigma", &sigma, 0.01, 0.05, "%.3f");

    // experiment selector
    const char* expNames[] = { "Default", "Shear Flow" };
    int expIdx = int(experiment);
    if (ImGui::Combo("Experiment", &expIdx, expNames, IM_ARRAYSIZE(expNames))) {
        experiment = Experiment(expIdx);
        gravity    = (experiment==Experiment::ShearFlow)
                    ? Vector3F::Zero()
                    : Vector3F(0,1,0);
    }
    if (experiment==Experiment::ShearFlow) {
        ImGui::Checkbox("Periodic X",   &periodicX);
        ImGui::InputDouble("ν (drag)",    &fluidViscosity, 0.01f,0.1f,"%.3f");
        ImGui::InputDouble("V₀ (max speed)", &V0,         0.1f,1.0f,"%.2f");
        ImGui::InputDouble("L (half‑height)", &L,         0.1f,1.0f,"%.2f");
    }
    static double  prevL         = L;
    static bool    prevPeriodic  = periodicX;

    bool needRebuild = false;

    if (experiment == Experiment::ShearFlow)
    {
        if (L != prevL)            { prevL = L;           needRebuild = true; }
        if (periodicX != prevPeriodic)
                                   { prevPeriodic = periodicX; needRebuild = true; }
    }

    if (needRebuild && !scenarioObjects.empty())
    {
        /* regenerate tunnel geometry (if one is present) */
        if (auto* tun = dynamic_cast<Tunnel2D*>(scenarioObjects[0].get()))
        {
            tun->halfWidth = L * 0.5;     // our convention
            tun->generateVertices();
        }           

        buildGridDataStructure();         // BB changed → new grid size
        insertParticlesIntoGrid();        // refill the grid
        renormalise();
        updateNeighborLists();            // update neighbour caches
    }

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
    // static const std::vector<const char*> poissonChoiceLabels = {"0.05", "0.15", "0.25", "0.35", "0.45"};

    // Determine the current index based on sim.poissonRatio and sim.poissonChoices.
    // int currentIndex = 0;
    // for (int i = 0; i < poissonChoices.size(); i++) {
    //     if (std::abs(poissonChoices(i) - poissonRatio) < 1e-6) {
    //         currentIndex = i;
    //         break;
    //     }
    // }

    // Create a single combo box for Poisson Ratio.
    // if (ImGui::Combo("Poisson Ratio", &currentIndex,
    //                  poissonChoiceLabels.data(),
    //                  static_cast<int>(poissonChoices.size()))) {
    //     // Update the simulation's poissonRatio based on the selected index.
    //     poissonRatio = poissonChoices(currentIndex);
    // }
}

void Simulation::updateCellSizeFromParticles() {
    maxParticleRadius = 0.0;
    for (auto &p : particles2D)
        maxParticleRadius = std::max(maxParticleRadius, p.radius);
    // pick whatever factor keeps each disc safely within one neighbour cell
    cellSize = 2.0 * maxParticleRadius;
}

void Simulation::buildGridDataStructure() 
{
    if (scenarioObjects.empty()) {
    std::cerr << "No scenario objects defined in simulation.\n";
    return;
    }
    updateCellSizeFromParticles();
    // Use the cached bounding box (BB) from the first scenario object.
    const BoundingBox &bbox = scenarioObjects[0]->getBoundingBox();
    F factor = 1.1; // Add a small buffer around the BB.
    minX = bbox.min_x; // * factor;
    maxX = bbox.max_x; // * factor;
    minY = bbox.min_y; // * factor;
    maxY = bbox.max_y; // * factor;

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
    // 1) Clear all cells
    for (auto &cell : grid)
        cell.clear();

    // 2) Standard insertion (possibly wrapping in X)
    for (I i = 0; i < (I)particles2D.size(); ++i) {
        const Particle2D &p = particles2D[i];
        // compute world‐space AABB of this particle
        F x0 = p.pos(0) - p.radius;
        F x1 = p.pos(0) + p.radius;
        F y0 = p.pos(1) - p.radius;
        F y1 = p.pos(1) + p.radius;

        int minCX = (int)std::floor((x0 - minX) / cellSize);
        int maxCX = (int)std::floor((x1 - minX) / cellSize);
        int minCY = (int)std::floor((y0 - minY) / cellSize);
        int maxCY = (int)std::floor((y1 - minY) / cellSize);

        for (int cy = minCY; cy <= maxCY; ++cy) {
            if (cy < 0 || cy >= numCellsY) continue;
            for (int cxRaw = minCX; cxRaw <= maxCX; ++cxRaw) {
                int cxWrapped;
                if (periodicX) {
                    // wrap X
                    cxWrapped = ((cxRaw % numCellsX) + numCellsX) % numCellsX;
                } else {
                    if (cxRaw < 0 || cxRaw >= numCellsX) continue;
                    cxWrapped = cxRaw;
                }
                grid[cy * numCellsX + cxWrapped].push_back(i);
            }
        }
    }

    // 3) **Only if** we're in periodic‑X mode, merge the two edge columns
    if (periodicX && numCellsX > 1) {
        for (int cy = 0; cy < numCellsY; ++cy) {
            int idxL = cy * numCellsX + 0;
            int idxR = cy * numCellsX + (numCellsX - 1);

            auto &leftCell  = grid[idxL];
            auto &rightCell = grid[idxR];

            // append right‐cell contents into left‐cell
            leftCell.insert(leftCell.end(),
                            rightCell.begin(),
                            rightCell.end());

            // (optional) keep symmetry so neighbor queries on the right edge see the same list
            rightCell = leftCell;
        }
    }
}



void Simulation::updateNeighborLists()
{
    /* clear old lists */
    for (auto& p : particles2D)
        p.neighborIndices.clear();

    /* constants */
    const int W = numCellsX;                // horizontal cell count

    /* -------- one pass per particle -------- */
    for (size_t i = 0; i < particles2D.size(); ++i)
    {
        Particle2D& p = particles2D[i];

        /* bounding box in *cell* coordinates – NOT clamped */
        F cxWrapped  = periodicX ? wrapX(p.pos(0)) : p.pos(0);
        int cxMinRaw = (int)std::floor((cxWrapped - p.radius - minX) / cellSize);
        int cxMaxRaw = (int)std::floor((cxWrapped + p.radius - minX) / cellSize);
        int cyMin    = std::max(0,  (int)std::floor((p.pos(1)-p.radius - minY)/cellSize));
        int cyMax    = std::min(numCellsY-1,
                                (int)std::floor((p.pos(1)+p.radius - minY)/cellSize));

        /* gather candidate indices (avoid duplicates with a set) */
        std::unordered_set<int> cand;
        for (int cy = cyMin; cy <= cyMax; ++cy)
        {
            for (int cxRaw = cxMinRaw; cxRaw <= cxMaxRaw; ++cxRaw)
            {
                int cx = cxRaw;

                /* wrap or reject in X */
                if (periodicX)
                    cx = ((cx % W) + W) % W;          // modulo wrap
                else if (cx < 0 || cx >= W)
                    continue;                         // outside → skip

                int g = cy * W + cx;                  // 1‑D cell index
                for (int j : grid[g])
                    if (j != (int)i) cand.insert(j);
            }
        }

        /* real‑overlap filtering (and symmetry) */
        for (int j : cand)
        {
            if (j <= (int)i) continue;               // keep ≤ once
            const Particle2D& q = particles2D[j];

            F dx = periodicDx(p.pos(0), p.ix, q.pos(0), q.ix);
            F dy = p.pos(1) - q.pos(1);
            if (std::sqrt(dx*dx + dy*dy) < p.radius + q.radius)
                p.neighborIndices.push_back(j);
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
        int n = (int)particles2D.size();
        for (int i = 0; i < n; ++i) {
            particles2D[i].pos = positions.segment<3>(3*i);

            if (periodicX)   // keep it canonical
                particles2D[i].pos(0) = wrapX(particles2D[i].pos(0));
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
        if (experiment == Experiment::Default){
            // Check if at least one scenario object is available.
            if (!scenarioObjects.empty()) {
                // Use the bounding box from the first scenario object.
                const BoundingBox &bbox = scenarioObjects[0]->getBoundingBox();
                F factor = 1.1; // Optional: apply a small buffer if needed.
                // You can adjust the baseline using the bounding box. For instance:
                F minYAdjusted = bbox.min_y; // Optionally, multiply by factor if required.
                value += gravity(1) * (p.pos(1) - minYAdjusted); // Compute energy relative to minY.
            } else {
                // Fallback: compute gravitational energy without a scenario object.
                value += gravity(1) * p.pos(1);
            }
        }
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
        } else if (p.BoundaryCollision == 3) {                    // ← NEW
            for (const auto& so : scenarioObjects) {
                if (auto* tun = dynamic_cast<Tunnel2D*>(so.get())) {
                    F py     = p.pos(1);
                    F r      = p.radius;
                    F min_y  = tun->BB.min_y;
                    F max_y  = tun->BB.max_y;
    
                    /* vertical penetration */
                    F overlap = 0.0;
                    if (py - r < min_y)       overlap = min_y - (py - r);
                    else if (py + r > max_y)  overlap = (py + r) - max_y;
    
                    if (overlap > 0.0)
                        value += 0.5 * overlapParam * overlap * overlap;
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
            F dx = periodicDx(p.pos(0), p.ix, q.pos(0), q.ix);   // <── single‑line change
            F dy = p.pos(1) - q.pos(1);
            F d = std::sqrt(dx * dx + dy * dy);
            // Compute overlap only if particles are close.
            F overlap = (p.radius + q.radius) - d;
            if (overlap > 0) {
                value += 0.5 * interactionParam * overlap * overlap / (1 + overlap * overlap);
            }
        }
    }

    //Prevent rigid body motion
    if (periodicX && !particles2D.empty())
    {
        F dx = particles2D[0].pos(0) - pinXref;
        value += 0.5 * pinK * dx * dx;
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
        if (experiment == Experiment::Default){
            // Gravitational gradient: only affects the y-component.
            gradient(3 * i)     = 0;
            gradient(3 * i + 1) = gravity(1); //p.mass *
            gradient(3 * i + 2) = 0;
        }
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
        } else if (p.BoundaryCollision == 3) {                    // ← NEW
            for (const auto& so : scenarioObjects) {
                if (auto* tun = dynamic_cast<Tunnel2D*>(so.get())) {
                    F py     = p.pos(1);
                    F r      = p.radius;
                    F min_y  = tun->BB.min_y;
                    F max_y  = tun->BB.max_y;
    
                    if (py - r < min_y)
                        gradient(3*i + 1) += overlapParam * (-min_y + py - r);
                    else if (py + r > max_y)
                        gradient(3*i + 1) += overlapParam * (py - max_y + r);
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
            F dx = periodicDx(p.pos(0), p.ix, q.pos(0), q.ix);   // <── single‑line change
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

    //Prevent rigid body motion
    if (periodicX && !particles2D.empty())
    {
        gradient(0) += pinK * (particles2D[0].pos(0) - pinXref);
        /*  (only x‑dof of particle 0; no effect on y or θ)  */
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
        } else if (p.BoundaryCollision == 3) {                    // ← NEW
            for (const auto& so : scenarioObjects) {
                if (dynamic_cast<Tunnel2D*>(so.get())) {
                    /* same constant stiffness as square walls, but only in y */
                    hessian.coeffRef(3*i + 1, 3*i + 1) += overlapParam;
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
            F dx = periodicDx(p.pos(0), p.ix, q.pos(0), q.ix);   // <── single‑line change
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

    //Prevent rigid body motion
    if (periodicX && !particles2D.empty())
        hessian.coeffRef(0,0) += pinK;
}

void Simulation::compute_energy_dyn(F &value) {

    F dynamicWeight = lambda / (timeStep * timeStep);
    compute_energy(value);
    
    if (globalState_1.size() == globalPositions.size()) {
        value += dynamicWeight * (0.5 * (globalPositions - globalState_1).squaredNorm() 
                                   - globalPositions.dot(globalState_1 - globalState_2));
    } else {
        globalState_1 = globalPositions;
        globalState_2 = globalPositions;
        value += dynamicWeight * (0.5 * (globalPositions - globalState_1).squaredNorm() 
                                   - globalPositions.dot(globalState_1 - globalState_2));
    }
    
    if (viscosity) {
        for (size_t i = 0; i < particles2D.size(); i++) {
            // Use the filtered effective neighbor count stored in the particle.
            F effectiveCountFiltered = particles2D[i].prevEffectiveCount;
            VectorXF diff = globalPositions.segment(3 * i, 3) - 
                            globalState_1.segment(3 * i, 3);
            value += 0.5 * viscosity_coeff * effectiveCountFiltered / (timeStep * timeStep) * diff.squaredNorm();
        }
    }
    if (experiment == Experiment::ShearFlow)
    {
        for (int i = 0; i < (int)particles2D.size(); ++i)
        {
            const int base = 3 * i;

            /* current & previous positions */
            F x1 = globalPositions(base    );
            F y1 = globalPositions(base + 1);
            F t1 = globalPositions(base + 2);
            F x0 = globalState_1  (base    );
            F y0 = globalState_1  (base + 1);
            F t0 = globalState_1  (base + 2);

            /* fluid velocity */
            F vfx , dvf_dy , dummy;
            shearFlowProfile(y1, vfx, dvf_dy, dummy);

            /* particle velocity components */
            F dvx = (x1 - x0)/timeStep - vfx;
            F dvy = (y1 - y0)/timeStep;
            F dvt = (t1 - t0)/timeStep;

            value += 0.5 * fluidViscosity *
                     (dvx*dvx + dvy*dvy + dvt*dvt);
        }
    }
    
}

void Simulation::compute_gradient_dyn(VectorXF &gradient) {
    
    F dynamicWeight = lambda / (timeStep * timeStep);
    compute_gradient(gradient);
    
    if (globalState_1.size() == globalPositions.size()) {
        gradient += dynamicWeight * (globalPositions - 2 * globalState_1 + globalState_2);
    } else {
        globalState_1 = globalPositions;
        globalState_2 = globalPositions;
        gradient += dynamicWeight * (globalPositions - 2 * globalState_1 + globalState_2);
    }
    
    if (viscosity) {
        for (size_t i = 0; i < particles2D.size(); i++) {
            F effectiveCountFiltered = particles2D[i].prevEffectiveCount;
            gradient.segment(3 * i, 3) += viscosity_coeff * effectiveCountFiltered / (timeStep * timeStep) *
                                          (globalPositions.segment(3 * i, 3) - 
                                           globalState_1.segment(3 * i, 3));
        }
    }
    if (experiment == Experiment::ShearFlow)
    {
        const F nu      = fluidViscosity;
        const F inv_dt  = 1.0 / timeStep;
        const F nu_over_dt  = nu * inv_dt;
        const F nu_over_dt2 = nu * inv_dt * inv_dt;

        for (int i = 0; i < (int)particles2D.size(); ++i)
        {
            const int base = 3 * i;

            /* positions */
            F x1 = globalPositions(base    );
            F y1 = globalPositions(base + 1);
            F t1 = globalPositions(base + 2);
            F x0 = globalState_1  (base    );
            F y0 = globalState_1  (base + 1);
            F t0 = globalState_1  (base + 2);

            /* fluid profile & derivative */
            F vfx , dvf_dy , dummy;
            shearFlowProfile(y1, vfx, dvf_dy, dummy);

            /* velocity differences */
            F dvx = (x1 - x0)*inv_dt - vfx;   // (v_p − v_f)_x
            F dvy = (y1 - y0)*inv_dt;         // v_p,y
            F dvt = (t1 - t0)*inv_dt;         // ω

            // -------- grad_x --------
            gradient(base) += nu_over_dt * dvx;

            // -------- grad_y --------
            gradient(base + 1) += nu_over_dt * dvy  // from v_p,y
                                - nu * dvx * dvf_dy; // from v_f(y)

            // -------- grad_theta ----
            gradient(base + 2) += nu_over_dt * dvt;
        }
    }
    
    
}

void Simulation::compute_hessian_dyn(SparseMatrixF &hessian) {
    
    F dynamicWeight = lambda / (timeStep * timeStep);
    compute_hessian(hessian);
    
    for (int i = 0; i < (3 * static_cast<int>(particles2D.size())); i++) {
        hessian.coeffRef(i, i) += dynamicWeight;
    }
    
    if (viscosity) {
        for (size_t i = 0; i < particles2D.size(); i++) {
            F effectiveCountFiltered = particles2D[i].prevEffectiveCount;
            for (int j = 0; j < 3; j++) {
                int index = static_cast<int>(3 * i + j);
                hessian.coeffRef(index, index) += viscosity_coeff * effectiveCountFiltered / (timeStep * timeStep);
            }
        }
    }
    if (experiment == Experiment::ShearFlow)
    {
        const F nu       = fluidViscosity;
        const F inv_dt   = 1.0 / timeStep;
        const F nu_dt2   = nu * inv_dt * inv_dt;
        const F k        = M_PI / L;           // π/L

        for (int i = 0; i < (int)particles2D.size(); ++i)
        {
            const int base = 3 * i;

            /* positions & velocities */
            F x1 = globalPositions(base    );
            F y1 = globalPositions(base + 1);
            F x0 = globalState_1  (base    );

            F sin_k_y , cos_k_y;
            {
                F vfx , dvf_dy , d2vf_dy2;
                shearFlowProfile(y1, vfx, dvf_dy, d2vf_dy2);
                sin_k_y = std::sin(k*y1);
                cos_k_y = std::cos(k*y1);
            }

            /* pre‑compute   dvx   */
            F dvx = (x1 - x0)*inv_dt - V0*sin_k_y;

            /* ===== diagonal blocks ===== */
            hessian.coeffRef(base    , base    ) += nu_dt2;   // d²E/dx²
            hessian.coeffRef(base + 2, base + 2) += nu_dt2;   // d²E/dθ²

            /* d²E/dy² :
               nu_dt2  from v_p,y term
             + nu * (V0 k)^2 * cos²(·)   from ∂dvx/∂y
             + nu * dvx * V0 * k^2 * sin(·)   from ∂²v_f/∂y²
            */
            F term1 = nu_dt2;
            F term2 = nu * (V0*V0) * k*k * cos_k_y*cos_k_y;
            F term3 = nu * dvx * V0 * k*k * sin_k_y;
            hessian.coeffRef(base + 1, base + 1) += term1 + term2 + term3;

            /* ===== off‑diagonal  d²E/dxdy  (symmetric) =====
               ∂grad_x/∂y = nu * (-V0 k cos) / dt
            */
            F off = -nu * V0 * k * cos_k_y * inv_dt;
            hessian.coeffRef(base    , base + 1) += off;
            hessian.coeffRef(base + 1, base    ) += off;
        }
    }
    
    
}

inline void Simulation::shearFlowProfile(F y,
    F& v_fx,      //  sin(π y/L)
    F& dvf_dy,    //  (π/L) cos(π y/L)
    F& d2vf_dy2)  // −(π/L)^2 sin(π y/L)
{
    const F k = M_PI / L;            // π / L
    v_fx     = V0 * std::sin(k * y);
    F c      = std::cos(k * y);
    F s      = std::sin(k * y);
    dvf_dy   = V0 * k * c;
    d2vf_dy2 = -V0 * k * k * s;
}

F Simulation::wrapX(F x) const {
    const F W = maxX - minX;

    // If the domain has not been initialised yet (W == 0)
    // simply return the original coordinate to avoid Inf / NaN.
    if (W == 0 || !periodicX)
        return x;

    return x - std::floor((x - minX) / W) * W;
}

inline F Simulation::periodicDx(F x1,int ix1, F x2,int ix2) const
{
    if (!periodicX) return x1 - x2;

    const F W = maxX - minX;
    if (W <= std::numeric_limits<F>::epsilon())
        return x1 - x2;

    /* true signed separation in the covering space */
    F dx = (x1 - x2) + (F)(ix1 - ix2) * W;

    /* wrap onto (‑½W , ½W]  –  round() gives the nearest integer */
    dx -= std::round(dx / W) * W;
    return dx;
}

void Simulation::renormalise()
{
    const F W = maxX - minX;

    for (int i = 0; i < (int)particles2D.size(); ++i) {
        auto &p = particles2D[i];

        /* old wrap counter */
        int old_ix = p.ix;

        /* keep pos.x in [minX,maxX) and update p.ix */
        while (p.pos(0) <  minX) { p.pos(0) += W; --p.ix; }
        while (p.pos(0) >= maxX) { p.pos(0) -= W; ++p.ix; }

        /* Δix since last frame */
        int dix = p.ix - old_ix;
        if (dix != 0) {
            F shift = (F)dix * W;
            globalState_1(3*i    ) += shift;
            globalState_2(3*i    ) += shift;
        }
    }
}

void Simulation::updateEffectiveNeighborCounts() {
    // Loop over each particle.
    for (size_t i = 0; i < particles2D.size(); i++) {
        F rawEffectiveCount = 0;

        // Fetch mod‑position + wrap counter for convenience
        const F xi   = globalPositions(3*i    );  // already ∈ [minX,maxX)
        const int ixi = particles2D[i].ix;

        // Compute the raw effective count using a soft kernel w(d) = exp(−d²/σ²)
        for (int j : particles2D[i].neighborIndices) {
            // x‑difference with correct periodic image
            const F xj   = globalPositions(3*j    );
            const int ixj = particles2D[j].ix;
            F dx = periodicDx(xi, ixi, xj, ixj);

            // y is non‑periodic
            F dy = globalPositions(3*i + 1)
                 - globalPositions(3*j + 1);

            F d2 = dx*dx + dy*dy;  // ignore θ
            rawEffectiveCount += std::exp(-d2 / (sigma * sigma));
        }

        // Store the raw count for this frame
        particles2D[i].effectiveCountCurrent = rawEffectiveCount;
        
        // Temporal filtering
        F filteredCount = alpha * rawEffectiveCount
                        + (1 - alpha) * particles2D[i].prevEffectiveCount;
        // Clamp to avoid runaway stiffness
        filteredCount = std::min(filteredCount, F(1.0));

        // **Remember to write it back** so next frame can filter against it
        particles2D[i].prevEffectiveCount = filteredCount;
    }
}

void Simulation::updateEffectiveNeighborCountsFinal() {
    // Loop over each particle.
    for (size_t i = 0; i < particles2D.size(); i++) {
        F rawEffectiveCount = 0;

        // Cached mod‐position + wrap counter
        const F  xi   = globalPositions(3*i    );  // ∈ [minX,maxX)
        const int ixi = particles2D[i].ix;

        // Soft‐kernel sum over neighbours
        for (int j : particles2D[i].neighborIndices) {
            const F  xj   = globalPositions(3*j    );
            const int ixj = particles2D[j].ix;

            F dx = periodicDx(xi, ixi, xj, ixj);
            F dy = globalPositions(3*i + 1)
                 - globalPositions(3*j + 1);

            F d2 = dx*dx + dy*dy;
            rawEffectiveCount += std::exp(-d2 / (sigma * sigma));
        }

        // Store the raw count
        particles2D[i].effectiveCountCurrent = rawEffectiveCount;
        
        // Temporal filtering
        F filteredCount = alpha * rawEffectiveCount
                        + (1 - alpha) * particles2D[i].prevEffectiveCount;
        filteredCount = std::min(filteredCount, F(1.0));

        // Write it back for the next frame
        particles2D[i].prevEffectiveCount = filteredCount;
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

    // ------------------------------------------------------------
    // ❶  Find the area where we may drop particles
    // ------------------------------------------------------------
    F xMin = -0.8 , xMax = 0.8;          // ← fall‑back values
    F yMin = -0.5 , yMax = 0.5;

    if (!scenarioObjects.empty()) {
        if (auto* tun = dynamic_cast<Tunnel2D*>(scenarioObjects[0].get())) {
            const BoundingBox& BB = tun->getBoundingBox();
            xMin = BB.min_x;
            xMax = BB.max_x;
            yMin = BB.min_y;
            yMax = BB.max_y;
        }
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<F> radiusDist(mean, std);   
    std::uniform_real_distribution<F> distX(xMin, xMax);

    // ------------------------------------------------------------
    // ❷  Rejection sampling with exact periodic distance
    // ------------------------------------------------------------
    for (int i = 0; i < numParticles; ++i) {
        bool ok = false;
        int attempts = 0;
        Particle2D cand(0.0, *this);
        cand.ix = 0;  // start in the base cell

        while (!ok && attempts < maxAttemptsPerParticle) {
            ++attempts;

            // sample radius
            F r = std::max<F>(F(0.01), radiusDist(gen));
            cand = Particle2D(r, *this);
            cand.ix = 0;

            // sample position in the fundamental domain
            std::uniform_real_distribution<F> distY(yMin + r, yMax - r);
            cand.pos(0) = distX(gen);
            cand.pos(1) = distY(gen);
            cand.pos(2) = 0.0f;

            // enforce pos.x ∈ [xMin, xMax)
            // (distX already does this, but just to be safe)
            cand.pos(0) = std::clamp(cand.pos(0), xMin, std::nextafter(xMax, xMin));

            // overlap test using the new 4‑arg periodicDx
            ok = true;
            for (auto& ex : particles) {
                F dx = periodicDx(
                    cand.pos(0), cand.ix,
                    ex.pos(0),   ex.ix
                );
                F dy = cand.pos(1) - ex.pos(1);
                F dist = std::sqrt(dx*dx + dy*dy);
                if (dist < cand.radius + ex.radius) {
                    ok = false;
                    break;
                }
            }
        }

        if (ok) {
            particles.push_back(std::move(cand));
        } else {
            std::cerr << "createRandomParticles2D: could not place particle "
                      << (i+1) << " after " << maxAttemptsPerParticle << " tries.\n";
        }
    }

    // All new particles have ix == 0 and pos.x in [xMin, xMax).
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

void Simulation::startScenarioAnimation() {
    if (!scenarioObjects.empty()) {
        animationTimer = 0;
        // Store the current position of the first scenario object as the starting point.
        initialScenarioPosition = scenarioObjects[0]->position;
        animateScenario = true;
    }
}

void Simulation::updateScenarioAnimation(F dt) {
    if (!animateScenario || scenarioObjects.empty())
        return;

    // Increment the timer.
    animationTimer += dt;
    if (animationTimer > animationDuration) {
        animationTimer = animationDuration;
        animateScenario = false; // Stop the animation once complete.
    }

    // Compute progress as a fraction between 0 and 1.
    F t = animationTimer / animationDuration;

    // Start from the initial position.
    Vector3F newPos = initialScenarioPosition;
    
    // Apply linear motion along the x-axis.
    // newPos(0) += t * animationDistance;
    
    // Superimpose a shaking motion along the y-axis.
    // Here, the y-offset is given by 0.5*sin(4PI*t) (t goes from 0 to 1).
    newPos(0) += animationDistance * sin(8.0 * M_PI * t);

    // Update each scenario object (regenerate its vertices, and update auxiliary structures).
    for (auto &so : scenarioObjects) {
         so->position = newPos;
         so->generateVertices();
         buildGridDataStructure();
         updateAuxiliaryStructures();
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

Tunnel2D::Tunnel2D(F halfLen, F halfWid, const Vector3F& pos)
    : halfWidth(halfWid), halfLength(halfLen)
{
    position = pos;
    generateVertices();
}

void Tunnel2D::generateVertices()
{
    vertices.clear();

    /* true physical Y–walls */
    BB.min_y = position(1) - halfWidth;
    BB.max_y = position(1) + halfWidth;

    /* we expose an *artificially long* X span so the viewer can draw it.
       It *doesn’t* constrain particles – X is handled by periodic wrapping
       inside Simulation. */
    BB.min_x = position(0) - halfLength;
    BB.max_x = position(0) + halfLength;

    /* rectangle mesh (purely visual) */
    vertices.emplace_back(BB.min_x, BB.min_y, position(2));
    vertices.emplace_back(BB.min_x, BB.max_y, position(2));
    vertices.emplace_back(BB.max_x, BB.max_y, position(2));
    vertices.emplace_back(BB.max_x, BB.min_y, position(2));
}

int Tunnel2D::detectCollision(const Particle2D& p) const
{
    /* only top / bottom walls act as barriers */
    if (p.pos(1) - p.radius < BB.min_y || p.pos(1) + p.radius > BB.max_y)
        return 3;                 // same collision code used by Square
    return 0;
}