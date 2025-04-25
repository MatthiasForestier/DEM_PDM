#pragma once

#include "CRLHelper/VecMatDef.h"
#include <vector>
#include <string>
#include <cmath>
#include <memory>

constexpr int DOF_FULL   = 3;   // x, y, θ  (internal)
constexpr int DOF    = 2;   // x, y     (optimizer)

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

class Tunnel2D : public ScenarioObject {
    public:
        /* halfWidth  = ±Y  extent  (the tunnel “height”)
           halfLength = only for drawing ‑– pick something long (e.g. 3 × halfWidth) */
        F halfWidth;
        F halfLength;
    
        BoundingBox BB;
    
        Tunnel2D(F halfLen, F halfWid, const Vector3F& pos);
    
        /* overrides */
        void generateVertices()               override;
        int  detectCollision(const Particle2D& p) const override;
        BoundingBox getBoundingBox()   const   override { return BB; }
    
        std::unique_ptr<ScenarioObject> clone() const override {
            return std::make_unique<Tunnel2D>(*this);
        }
    };

//------------------------------------------------------------------------------
// Particles
//------------------------------------------------------------------------------

class Particle2D {
public:
    // Position: [x, y]
    Vector2F pos;
    // Velocity: [x_dot, y_dot]
    Vector2F vel;
    // Acceleration: [x_ddot, y_ddot]
    Vector2F acc;

    I BoundaryCollision = 0;
    I SquareCollision = 0;
    std::vector<I> neighborIndices;
    F effectiveCountCurrent = 0;  
    F prevEffectiveCount = 0; 
    F radius;
    I ix = 0; 
    F mass;
    F inertia; // Moment of inertia for a disc: (1/2)*mass*radius^2
    Color color = Color(0.678f, 0.847f, 0.902f);

    // Constructor (implementation in cpp file).
    Particle2D(F radius, const Simulation& simParams);
};

class Particle3D {
public:
    // Position: [x, y, z]
    Vector3F pos;
    // Velocity: [x_dot, y_dot, z_dot]
    Vector3F vel;
    // Acceleration: [x_ddot, y_ddot, z_ddot]
    Vector3F acc;

    I BoundaryCollision = 0;
    I SquareCollision = 0;
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
    const int maxAttemptsPerParticle = 10000;
    std::vector<std::unique_ptr<ScenarioObject>> scenarioObjects;
    bool use3D = false;

    // Particle-related parameters.
    F density = 2780.0; // [kg/m^3]
    F overlapParam = 100000000;
    F interactionParam = 10000000;

    //Particles shape related parameters.
    bool is_circular = true;
    bool is_square = false;
    bool is_elliptical = false;
    F radiusMean = 0.05;  // [m]
    F radiusStd = 0.005;  // [m]

    // Grid-related parameters.
    F maxParticleRadius = 0.0;
    F cellSize = 0.0;
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
    F viscosityCoeff = 0.1;
    F kernelSigma = 3 * radiusMean;  
    F alpha = 0.5;

    /// Animation parameters.
    bool animateScenario = false;
    F animationTimer = 0;         // Elapsed time.
    F animationDuration = 5.0;      // Duration (seconds).
    F animationDistance = 0.2;      // Distance along x.
    Vector3F initialScenarioPosition = Vector3F::Zero();

    // Particles.
    VectorXF globalPositions;
    VectorXF globalState_1;
    VectorXF globalState_2;
    std::vector<Particle2D> particles2D;
    std::vector<Particle3D> particles3D;
    F minParticleDiam = 100;
    SparseMatrixF M; //mass matrix

    //Experiment parameters.
    enum class Experiment { Default = 0, ShearFlow = 1 };
    Experiment experiment = Experiment::Default;

    // Shear‐flow parameters
    bool periodicX       = false;    // wrap X
    F  fluidViscosity    = 1.0f;     // ν in F_drag = ν ||v_p−v_f||²
    F  V0                = 0.1f;     // max fluid speed
    F  L                 = 1.0f;     // half‐height of tunnel
    F     pinK      = 1e-2;       // << spring stiffness
    F     pinXref   = 0.0;        // << reference position


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
          is_circular(other.is_circular),
          is_square(other.is_square),
          is_elliptical(other.is_elliptical),
          radiusMean(other.radiusMean),
          radiusStd(other.radiusStd),
          maxParticleRadius(other.maxParticleRadius),
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
          viscosityCoeff(other.viscosityCoeff),
          kernelSigma(other.kernelSigma),
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
          particles3D(other.particles3D),
          minParticleDiam(other.minParticleDiam),
          M(other.M),
          experiment(other.experiment),
          periodicX(other.periodicX),
          fluidViscosity(other.fluidViscosity),
          V0(other.V0),
          L(other.L),
          pinK(other.pinK),
          pinXref(other.pinXref)
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

    inline void shearFlowProfile(F y,F& v_fx,F& dvf_dy,F& d2vf_dy2);
    F wrapX(F x) const;           // defined in .cpp
    inline F periodicDx(F x1, int ix1, F x2, int ix2) const;
    void renormalise();
    void updateCellSizeFromParticles();
    void minParticles();
    void computeCFL ();
    void buildMassMatrix(SparseMatrixF& M) const;

};
