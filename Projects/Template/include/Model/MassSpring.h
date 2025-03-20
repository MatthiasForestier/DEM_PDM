#pragma once

#include "CRLHelper/VecMatDef.h"
#include <vector>
#include <string>
#include <cmath>

// Particle2D is defined first.
class Simulation; 
class ScenarioObject;

struct BoundingBox {
    F min_x, max_x, min_y, max_y;
};

struct Color {
    float r, g, b;
    Color(float r_ = 1.0f, float g_ = 1.0f, float b_ = 1.0f)
      : r(r_), g(g_), b(b_) {}
};

class Particle2D {
public:
    // Position: [x, y, theta]
    Vector3F pos;
    // Velocity: [x_dot, y_dot, theta_dot]
    Vector3F vel;
    // Acceleration: [x_ddot, y_ddot, theta_ddot]
    Vector3F acc;
    // Forces
    Vector3F force;
    // Moments
    Vector3F moment;

    int BoundaryCollision = 0;
    int SquareCollision = 0;
    std::vector<I> neighborIndices;
    F radius;
    F mass;
    F inertia; // Moment of inertia for a disc: (1/2)*mass*radius^2
    Color color = Color(0.678f, 0.847f, 0.902f);

    // Constructor that computes mass and inertia using the density from simParams.
    Particle2D(F radius, const Simulation& simParams);
};

// Particle3D is defined similarly.
class Particle3D {
public:
    // Position: [x, y, z, theta, phi, psi]
    Vector6F pos;
    // Velocity: [x_dot, y_dot, z_dot, theta_dot, phi_dot, psi_dot]
    Vector6F vel;
    // Acceleration: [x_ddot, y_ddot, z_ddot, theta_ddot, phi_ddot, psi_ddot]
    Vector6F acc;
    // Forces   
    Vector3F force;
    // Moments
    Vector3F moment;

    int BoundaryCollision = 0;
    int SquareCollision = 0;
    std::vector<I> neighborIndices;
    F radius;
    F mass;
    F inertia; // Moment of inertia for a sphere: (2/5)*mass*radius^2
    Color color = Color(0.678f, 0.847f, 0.902f);

    // Constructor that computes mass and inertia using the density from simParams.
    Particle3D(F radius, const Simulation& simParams);
};

// Now define Simulation. Note that Particle2D is already fully defined.
class Simulation {
public:
    // World parameters
    F timeStep = 1e-3;
    Vector3F gravity = Vector3F(0.0, 100, 0.0);
    F endTime = 1.0;
    F outputInterval = 1e-2;
    std::string outputDirectory = "output";
    I numParticles = 10;
    const int maxAttemptsPerParticle = 100;
    // Display dimensions commented out for now.
    // int displayWidth = 800;
    // int displayHeight = 800;
    std::vector<std::unique_ptr<ScenarioObject>> scenarioObjects; 
    bool use3D = false;

    // Particle-related parameters
    F density = 2780.0; // [kg/m^3] used for mass computation
    F overlapParam = 1000000;
    F interactionParam = 100000;
    F contactBondNormalStiffness = 35.0;   // [MN/m]
    F contactStiffnessRatio = 0.5;         // (range: 0 to 1)
    F interparticleFriction = 0.83;
    F contactBondNormalStrength = 180.0;   // [N]
    F contactBondShearStrength = 180.0;    // [N]
    F particleWallContactNormalStiffness = 950.0;   // [MN/m]
    F particleWallContactTangentialStiffness = 950.0; // [MN/m]
    F particleWallFriction = 0.0;
    F translationalDamping = 0.2;
    F rotationalDamping = 0.2;
    F youngsModulusMin = 7e10;     // [Pa]
    F youngsModulusMax = 1e11;     // [Pa]
    VectorXF poissonChoices = (VectorXF(5) << 0.05, 0.15, 0.25, 0.35, 0.45).finished();
    F poissonRatio = 0.25;
    F loadingVelocity = 0.5;       // [m/s]
    F mean = 0.04;                 // [m]
    F std = 0.001;                  // [m]

    // Particles 
    VectorXF globalPositions;
    std::vector<Particle2D> particles2D;
    std::vector<Particle3D> particles3D;

    // Grid related parameters
    F cellSize = mean * 3;
    std::vector<std::vector<int>> grid;
    int numCellsX;
    int numCellsY;
    F minX;
    F maxX;
    F minY;
    F maxY;

    ///Dynamic formulation
    bool dynamic = false;

public:
    void makeConfigMenu();
    std::vector<Particle2D> createRandomParticles2D();
    ScenarioObject createScenarioCircle2D(int numSegments,F radius);
    ScenarioObject createScenarioSquare2D();
    void detectBoundaryCollision2D(Particle2D &p) const;
    void detectParticleCollision();
    VectorXF getGlobalState();
    void setGlobalState(const VectorXF &state);
    void compute_energy(F &value) const;
    void compute_gradient(VectorXF &gradient) const;
    void compute_hessian(SparseMatrixF &hessian) const;
    void updateGlobalPositions();
    void applyGlobalPositions(const VectorXF &positions);
    void buildGridDataStructure();
    void insertParticlesIntoGrid();
    void updateNeighborLists();
    void updateAuxiliaryStructures();
    void colorParticleRed(int particleID);
};

class ScenarioObject {
    public:
        std::vector<Vector3F> vertices;
        Vector3F position;
        
        virtual ~ScenarioObject() = default;
        virtual void generateVertices() = 0;
        
        // Virtual function for collision detection.
        // Return 0 for no collision, 1 for circle, 2 for square, etc.
        virtual int detectCollision(const Particle2D &p) const = 0;
        virtual BoundingBox getBoundingBox() const = 0;
    };
    
    class Square : public ScenarioObject {
    public:
        F halfLength;
        F halfWidth;
        // Cache bounds to avoid recomputation.
        F min_x, max_x, min_y, max_y;
        BoundingBox BB;

        // Constructor.
        Square(F halfLength, F halfWidth, const Vector3F& pos);
        
        // Overridden functions.
        virtual void generateVertices() override;
        virtual int detectCollision(const Particle2D &p) const override;
        virtual BoundingBox getBoundingBox() const override { return BB; }
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
        virtual int detectCollision(const Particle2D &p) const override;
        virtual BoundingBox getBoundingBox() const override { return BB; }
    };
