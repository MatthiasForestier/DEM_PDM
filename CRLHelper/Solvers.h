#pragma once

#include "CRLHelper/VecMatDef.h"
#include "Projects/Template/include/Model/MassSpring.h"
#include <string>

// Forward declaration
class Simulation;

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
    // -> each solver implements its own version
    virtual void step(F dt) = 0;

    // A virtual or non-virtual function that performs multiple steps in a row
    // Default behavior: just call step(dt) repeatedly
    virtual void stepMultipleTimes();

    // Factory method to create a solver based on a string identifier
    static Solver* createSolver(const std::string& solverType, Simulation* simulation);

    // Collision Detection: optionally override in derived classes if needed
    virtual void collisionDetection();

    // Force Calculation: optionally override in derived classes if needed
    virtual void forceCalculation();

    template <typename VectorType,
    typename ResidualFunc,
    typename JacobianFunc>
    VectorType newtonMethod(const VectorType &initialGuess,
                    ResidualFunc residual,
                    JacobianFunc jacobian,
                    int maxIterations = 10,
                    double tolerance = 1e-6);
    };

// Forward Euler Solver: updates state using the explicit (Forward Euler) method
class ForwardEulerSolver : public Solver {
public:
    ForwardEulerSolver(Simulation* simulation) : Solver(simulation) {}

    // Override the step function with Forward Euler implementation
    void step(F dt) override;
};

// Backward Euler Solver: updates state using an implicit (Backward Euler) method
class BackwardEulerSolver : public Solver {
public:
    BackwardEulerSolver(Simulation* simulation) : Solver(simulation) {}

    // Override the step function with Backward Euler implementation
    void step(F dt) override;
};

