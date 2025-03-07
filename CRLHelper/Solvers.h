#pragma once

#include "CRLHelper/VecMatDef.h"

#include "Projects/Template/include/Model/MassSpring.h"  
#include <string>

// Abstract base class for solvers
class Solver {
protected:
    Simulation* sim;  // Pointer to the simulation instance
public:
    // Constructor takes a pointer to a Simulation object
    Solver(Simulation* simulation) : sim(simulation) {}

    // Virtual destructor for proper cleanup of derived classes
    virtual ~Solver() {}

    // Pure virtual function to perform one simulation step
    virtual void step(F dt) = 0;

    // Optional: Factory method to create a solver based on a string identifier
    static Solver* createSolver(const std::string& solverType, Simulation* simulation);

    // Collision Detection per particle
    virtual void collisionDetection();
    // Force Calculation per particle
    virtual void forceCalculation();
};

// Forward Euler Solver: updates state using the explicit method
class ForwardEulerSolver : public Solver {
public:
    ForwardEulerSolver(Simulation* simulation) : Solver(simulation) {}

    // Override the step function with Forward Euler implementation
    virtual void step(F dt) override;
};

// Backward Euler Solver: updates state using the implicit method
class BackwardEulerSolver : public Solver {
public:
    BackwardEulerSolver(Simulation* simulation) : Solver(simulation) {}

    // Override the step function with Backward Euler implementation
    virtual void step(F dt) override;
};

