#pragma once

#include "CRLHelper/VecMatDef.h"
#include <vector>
#include <string>
#include <cmath>
#include <memory>

//------------------------------------------------------------------------------
// Basic types
//------------------------------------------------------------------------------
struct BoundingBox {
    F min_x, max_x, min_y, max_y;
};

struct Color {
    float r, g, b;
    Color(float r_ = 1.0f, float g_ = 1.0f, float b_ = 1.0f)
      : r(r_), g(g_), b(b_) {}
};

//------------------------------------------------------------------------------
// Forward declaration for Simulation (used in Particle constructors)
//------------------------------------------------------------------------------
class Simulation;

//------------------------------------------------------------------------------
// Scenario Objects
//------------------------------------------------------------------------------

class ScenarioObject {
public:
    std::vector<Vector3F> vertices;
    Vector3F position;
    
    virtual ~ScenarioObject() = default;
    virtual void generateVertices() = 0;
    virtual int detectCollision(const class Particle2D &p) const = 0;
    virtual BoundingBox getBoundingBox() const = 0;
    // Virtual clone function for deep copying.
    virtual std::unique_ptr<ScenarioObject> clone() const = 0;
};

class Square : public ScenarioObject {
public:
    F halfLength;
    F halfWidth;
    // Cached bounds.
    F min_x, max_x, min_y, max_y;
    BoundingBox BB;

    // Constructor.
    Square(F halfLength, F halfWidth, const Vector3F& pos);
    
    // Overridden functions.
    virtual void generateVertices() override;
    virtual int detectCollision(const class Particle2D &p) const override;
    virtual BoundingBox getBoundingBox() const override { return BB; }
    // Clone method.
    virtual std::unique_ptr<ScenarioObject> clone() const override {
        return std::make_unique<Square>(*this);
    }
};

class Circle : public ScenarioObject {
public:
    F radius;      // Radius of the circle.
    I numSegments; // Number of segments to approximate the circle.
    BoundingBox BB;

    // Constructor.
    Circle(F radius, I numSegments, const Vector3F& pos);
    
    // Overridden functions.
    virtual void generateVertices() override;
    virtual int detectCollision(const class Particle2D &p) const override;
    virtual BoundingBox getBoundingBox() const override { return BB; }
    // Clone method.
    virtual std::unique_ptr<ScenarioObject> clone() const override {
        return std::make_unique<Circle>(*this);
    }
};

//------------------------------------------------------------------------------
// Particles
//------------------------------------------------------------------------------

class Particle2D {
public:
    // Position: [x, y, theta]
    Vector3F pos;
    // Velocity: [x_dot, y_dot, theta_dot]
    Vector3F vel;
    // Acceleration: [x_ddot, y_ddot, theta_ddot]
    Vector3F acc;
    // Forces.
    Vector3F force;
    // Moments.
    Vector3F moment;

    int BoundaryCollision = 0;
    int SquareCollision = 0;
    std::vector<I> neighborIndices;
    F effectiveCountCurrent = 0;  
    F prevEffectiveCount = 0; 
    F radius;
    F mass;
    F inertia; // Moment of inertia for a disc: (1/2)*mass*radius^2
    Color color = Color(0.678f, 0.847f, 0.902f);

    // Constructor (implementation in cpp file).
    Particle2D(F radius, const Simulation& simParams);
};

class Particle3D {
public:
    // Position: [x, y, z, theta, phi, psi]
    Vector6F pos;
    // Velocity: [x_dot, y_dot, z_dot, theta_dot, phi_dot, psi_dot]
    Vector6F vel;
    // Acceleration: [x_ddot, y_ddot, z_ddot, theta_ddot, phi_ddot, psi_ddot]
    Vector6F acc;
    // Forces.
    Vector3F force;
    // Moments.
    Vector3F moment;

    int BoundaryCollision = 0;
    int SquareCollision = 0;
    std::vector<I> neighborIndices;
    F effectiveCountCurrent = 0;  
    F prevEffectiveCount = 0; 
    F radius;
    F mass;
    F inertia; // Moment of inertia for a sphere: (2/5)*mass*radius^2
    Color color = Color(0.678f, 0.847f, 0.902f);

    // Constructor (implementation in cpp file).
    Particle3D(F radius, const Simulation& simParams);
};

//------------------------------------------------------------------------------
// Simulation
//------------------------------------------------------------------------------

class Simulation {
public:
    // World parameters.
    F timeStep = 1e-2;
    F lambda = 1e0;
    Vector3F gravity = Vector3F(0.0, 1, 0.0);
    F endTime = 1.0;
    F outputInterval = 1e-2;
    std::string outputDirectory = "output";
    I numParticles = 100;
    const int maxAttemptsPerParticle = 100;
    std::vector<std::unique_ptr<ScenarioObject>> scenarioObjects;
    bool use3D = false;

    // Particle-related parameters.
    F density = 2780.0; // [kg/m^3]
    F overlapParam = 1000000;
    F interactionParam = 100000;
    F mean = 0.03;  // [m]
    F std = 0.01;  // [m]

    // Grid-related parameters.
    F cellSize = mean * 3;
    int numCellsX;
    int numCellsY;
    F minX;
    F maxX;
    F minY;
    F maxY;
    std::vector<std::vector<int>> grid;

    /// Dynamic formulation.
    bool dynamic = false;
    bool friction = false;
    bool viscosity = false;
    F viscosity_coeff = 0.1;
    F sigma = 3 * mean;  
    F alpha = 0.5;

    /// Animation parameters.
    bool animateScenario = false;
    F animationTimer = 0;         // Elapsed time.
    F animationDuration = 2.0;      // Duration (seconds).
    F animationDistance = 0.2;      // Distance along x.
    Vector3F initialScenarioPosition = Vector3F::Zero();

    // Particles.
    VectorXF globalPositions;
    VectorXF globalState_1;
    VectorXF globalState_2;
    std::vector<Particle2D> particles2D;
    std::vector<Particle3D> particles3D;

public:
    Simulation() = default;

    // Custom copy constructor with members initialized in the order they are declared.
    Simulation(const Simulation &other)
        : timeStep(other.timeStep),
          lambda(other.lambda),
          gravity(other.gravity),
          endTime(other.endTime),
          outputInterval(other.outputInterval),
          outputDirectory(other.outputDirectory),
          numParticles(other.numParticles),
          maxAttemptsPerParticle(other.maxAttemptsPerParticle),
          scenarioObjects(),  // We'll deep-copy these below.
          use3D(other.use3D),
          density(other.density),
          overlapParam(other.overlapParam),
          interactionParam(other.interactionParam),
          mean(other.mean),
          std(other.std),
          cellSize(other.cellSize),
          numCellsX(other.numCellsX),
          numCellsY(other.numCellsY),
          minX(other.minX),
          maxX(other.maxX),
          minY(other.minY),
          maxY(other.maxY),
          grid(other.grid),
          dynamic(other.dynamic),
          friction(other.friction),
          viscosity(other.viscosity),
          viscosity_coeff(other.viscosity_coeff),
          sigma(other.sigma),
          alpha(other.alpha),
          animateScenario(other.animateScenario),
          animationTimer(other.animationTimer),
          animationDuration(other.animationDuration),
          animationDistance(other.animationDistance),
          initialScenarioPosition(other.initialScenarioPosition),
          globalPositions(other.globalPositions),
          globalState_1(other.globalState_1),
          globalState_2(other.globalState_2),
          particles2D(other.particles2D),
          particles3D(other.particles3D)
    {
        // Deep copy scenario objects using clone().
        for (const auto &obj : other.scenarioObjects) {
            if (obj) {
                scenarioObjects.push_back(obj->clone());
            }
        }
    }

    // Member function declarations.
    void makeConfigMenu();
    std::vector<Particle2D> createRandomParticles2D();
    void updateEffectiveNeighborCounts();
    void updateEffectiveNeighborCountsFinal();
    ScenarioObject createScenarioCircle2D(int numSegments, F radius);
    ScenarioObject createScenarioSquare2D();
    void detectBoundaryCollision2D(Particle2D &p) const;
    void detectParticleCollision();
    VectorXF getGlobalState();
    void setGlobalState(const VectorXF &state);

    void compute_energy(F &value) const;
    void compute_gradient(VectorXF &gradient) const;
    void compute_hessian(SparseMatrixF &hessian) const;

    void compute_energy_dyn(F &value);
    void compute_gradient_dyn(VectorXF &gradient);
    void compute_hessian_dyn(SparseMatrixF &hessian);

    void updateGlobalPositions();
    void applyGlobalPositions(const VectorXF &positions);
    void buildGridDataStructure();
    void insertParticlesIntoGrid();
    void updateNeighborLists();
    void updateAuxiliaryStructures();
    void colorParticleRed(int particleID);

    void startScenarioAnimation();
    void updateScenarioAnimation(F dt);
};
