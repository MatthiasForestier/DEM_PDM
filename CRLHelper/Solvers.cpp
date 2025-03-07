#include "Solvers.h"
#include <string>

Solver* Solver::createSolver(const std::string& solverType, Simulation* simulation) {
    if (solverType == "Forward Euler") {
        return new ForwardEulerSolver(simulation);
    } else if (solverType == "Backward Euler") {
        return new BackwardEulerSolver(simulation);
    }
    // Add more else-if for other solvers, or default:
    return new ForwardEulerSolver(simulation);
}

void ForwardEulerSolver::step(double dt) {
    if(sim->use3D){
        // Implement the forward Euler step here.
        // e.g., update simulation state using forward Euler method.
    }else{
        // Implement the forward Euler step here.
        // e.g., update simulation state using forward Euler method.
    }
    // Implement the forward Euler step here.
    // e.g., update simulation state using forward Euler method.
}

void BackwardEulerSolver::step(double dt) {
    // Implement the backward Euler step here.
    // e.g., update simulation state using backward Euler method.
}

void Solver::collisionDetection(){
    // Implement the collision detection here.
    // e.g., check for collisions between particles.
}

// Force Calculation per particle
void Solver::forceCalculation(){
    // Implement the force calculation here.
    // e.g., calculate forces between particles.
}