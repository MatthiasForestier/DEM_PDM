#include "Solvers.h"
#include "CRLHelper/VecMatDef.h"

#include <string>
#include <iostream>
#include <chrono>

// ---------------------------
// Newton's method (template)
// ---------------------------
// template <typename VectorType,
//           typename ResidualFunc,
//           typename JacobianFunc>

// VectorType Solver::newtonMethod(const VectorType &initialGuess,
//                                 ResidualFunc residual,
//                                 JacobianFunc jacobian,
//                                 int maxIterations,
//                                 double tolerance)
// {
//     VectorType x = initialGuess;
//     for (int iter = 0; iter < maxIterations; ++iter) {
//         // Compute the residual
//         VectorType r = residual(x);

//         // Check convergence
//         if (r.norm() < tolerance) {
//             break;
//         }

//         // Compute the Jacobian
//         auto J = jacobian(x);

//         // Solve J * delta = -r  for delta
//         VectorType delta = J.fullPivLu().solve(-r);

//         // Update the iterate
//         x += delta;

//         // Check the update size
//         if (delta.norm() < tolerance) {
//             break;
//         }
//     }
//     return x;
// }


// Factory: create the correct solver based on string key
Solver* Solver::createSolver(const std::string& solverType, Simulation* simulation) {
    if (solverType == "Forward Euler") {
        return new ForwardEulerSolver(simulation);
    } else if (solverType == "Backward Euler") {
        return new BackwardEulerSolver(simulation);
    }
    // Default or fallback if unknown
    return new ForwardEulerSolver(simulation);
}

// Default collision detection (can be overridden per solver if needed)
void Solver::collisionDetection() {
    // Example: check collisions for 2D or 3D
    // ...
}

// Default force calculation (can be overridden if needed)
void Solver::forceCalculation() {
    // Example: apply gravity, call collision detection, etc.
    if (sim->use3D) {
        for (auto &particle : sim->particles3D) {
            particle.force = Vector3F::Zero();
            particle.force += sim->gravity * particle.mass;
        }
    } else {
        for (auto &particle : sim->particles2D) {
            particle.force = Vector3F::Zero();
            particle.force += sim->gravity * particle.mass;
        }
    }
}

void Solver::stepMultipleTimes() {
    // Desired display frequency (frame time in seconds)
    F freqDisplay = 1.0 / 100.0;
    int numSteps = static_cast<int>(freqDisplay / sim->timeStep);
    if (numSteps < 1) {
        numSteps = 1;
    }
    // Record the starting time
    auto startTime = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < numSteps; ++i) {
        // Perform one simulation step
        step(sim->timeStep);

        // Measure elapsed time
        auto currentTime = std::chrono::high_resolution_clock::now();
        std::chrono::duration<F> elapsed = currentTime - startTime;
        if (elapsed.count() >= freqDisplay) {
            // If the computation time has exceeded the display interval, break early
            break;
        }
    }
}

// -----------------------------------------------------------------------------
// Forward Euler Implementation
// -----------------------------------------------------------------------------
void ForwardEulerSolver::step(F dt) {
    // If it's a 2D simulation, update each particle in 2D
    if (!sim->use3D) {
        forceCalculation();
        for (auto &particle : sim->particles2D) {
            particle.acc = particle.force / particle.mass;
            particle.vel += dt * particle.acc;
            particle.pos += dt * particle.vel;
        }
    }
    else {
        forceCalculation();
        for (auto &particle : sim->particles3D) {
            particle.acc.head<3>() = particle.force / particle.mass;
            particle.vel.head<3>() += dt * particle.acc.head<3>();
            particle.pos.head<3>() += dt * particle.vel.head<3>();
        }
    }
}

// -----------------------------------------------------------------------------
// Backward Euler Implementation
// (For demonstration, this is just a placeholder)
// -----------------------------------------------------------------------------
void BackwardEulerSolver::step(F dt) {

    //if (!sim->use3D) {
        // 2D case
    //     forceCalculation();
    //     for (auto &particle : sim->particles2D) {
    //         // Old state
    //         Vector3F v_old = particle.vel;
    //         Vector3F x_old = particle.pos;
    //         F mass = particle.mass;
            
    //         // Define the residual R(v_new): WRONG!
    //         auto residual = [&](const Vector3F &v_new) -> Vector3F {
    //             return v_new - v_old - dt * particle.acc;
    //         };

    //         // Define the Jacobian dR/dv_new; for constant force => Identity: WRONG!!
    //         auto jacobian = [&](const Vector3F &) -> Matrix3F {
    //             return Matrix3F::Identity();
    //         };

    //         // Choose an initial guess for v_new (e.g. forward Euler guess)
    //         Vector3F v_guess = v_old + dt * sim->gravity;

    //         // Solve via Newton's method
    //         Vector3F v_new = newtonMethod(v_guess, residual, jacobian, /*maxIter*/ 10, /*tol*/ 1e-6);

    //         // Update velocity and position
    //         particle.vel = v_new;
    //         particle.pos = x_old + dt * v_new;
    //     }
    // } else {
    //     // 3D case
    //     for (auto &particle : sim->particles3D) {
    //         Vector3F v_old = particle.vel;
    //         Vector3F x_old = particle.pos;
    //         F mass = particle.mass;

    //         auto residual = [&](const Vector3F &v_new) -> Vector3F {
    //             // Same logic in 3D
    //             return v_new - v_old - dt * sim->gravity;
    //         };

    //         auto jacobian = [&](const Vector3F &) -> Matrix3F {
    //             return Matrix3F::Identity();
    //         };

    //         Vector3F v_guess = v_old + dt * sim->gravity;
    //         Vector3F v_new   = newtonMethod(v_guess, residual, jacobian);

    //         particle.vel = v_new;
    //         particle.pos = x_old + dt * v_new;
    //     }
    //}
}

// auto jacobian = [&](const Vector3F &v_new) -> Matrix3F {
//     // 1) Compute x_new = x_old + dt*v_new
//     // 2) Compute df/dx at x_new
//     // 3) Multiply by chain rule factor dt/m (in a velocity-based scheme)
//     // 
//     // So the Jacobian might look like:
//     // I - (dt/m) * dF/dX(x_new)*dt, etc.

//     // For example, if f is a spring force: f(x) = -k (x - x_anchor),
//     // then dF/dX = -k I.
//     // So the Jacobian becomes: I - (dt^2 k / m) I.

//     Matrix3F df_dx = computeDfDx(x_new);  // your function for partial derivatives
//     Matrix3F I = Matrix3F::Identity();

//     return I - (dt * dt / mass) * df_dx;
// };
