#pragma once

#include <vector>
#include <string>
#include <cmath>
#include <memory>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <algorithm>
#include <iostream>
#include <random>
#include <unordered_set>

#include "CRLHelper/MapleHelper.h"
#include "CRLHelper/VecMatDef.h"

constexpr int DOF_FULL   = 3;   // x, y, θ  (internal)
constexpr int DOF       =  2;

/// number of active degrees of freedom per 2-D particle

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

    // ---------- deformation of particles -----------------
    F epsV = 0.0;
    F epsVdot = 0.0;

    // ---------- collision bookkeeping -----------------
    I BoundaryCollision = 0;
    I SquareCollision = 0;

    // ---------- neighbor bookkeeping -----------------
    std::vector<I> neighborIndices;
    F effectiveCountCurrent = 0;  
    F prevEffectiveCount = 0;
    
    // ---------- particle properties -----------------
    F radius;
    I ix = 0; 
    F mass;
    F inertia; // Moment of inertia for a disc: (1/2)*mass*radius^2
    Color color = Color(0.678f, 0.847f, 0.902f);

    // ---------- deformation‑viscosity bookkeeping -----------------
    Matrix2F F_prev = Matrix2F::Identity();   ///< F at previous step
    Matrix2F D      = Matrix2F::Zero();       ///< symmetric rate‑of‑def.
    Vector2F X0     = Vector2F::Zero();       ///< reference position

    // Constructor (implementation in cpp file).
    Particle2D(F radius, const Simulation& simParams);
    F effectiveRadius() const;
};

class Particle3D {
public:
    // Position: [x, y, z]
    Vector3F pos;

    F epsV = 0.0;
    F epsVdot = 0.0;
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
    Vector3F gravity = Vector3F(0.0, 0.0, 0.0);
    F endTime = 1.0;
    F outputInterval = 1e-2;
    std::string outputDirectory = "output";
    I numParticles = 42;
    const int maxAttemptsPerParticle = 10000;
    std::vector<std::unique_ptr<ScenarioObject>> scenarioObjects;
    bool use3D = false;

    // Particle-related parameters.
    F density = 100.0; // [kg/m^3]
    F overlapParam = 10000000;
    F interactionParam = 10000000;
    F Young  = 215000;     // material E (could be global)
    F Poisson = 0.3;    // ν  (or exactly 0.5 for incompressible trick)
    F K;
    bool boolSoftDEM = false; // Soft-DEM flag
    F a = 0.58f; 
    F b = 2.0f;

    //Particles shape related parameters.
    bool is_circular = true;
    bool is_square = false;
    bool is_elliptical = false;
    F radiusMean = 0.04;  // [m]
    F radiusStd = 0.005;  // [m]
    bool densify = false;
    bool bidispersed = false;  ///< NEW – toggle two‑size mode
    F  radiusBig   = 0.05;     ///< used only if bidispersed==true
    F  radiusSmall = 0.02;
    int numBig     = 0;        ///< ditto
    int numSmall   = 0;

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
    bool deformation_viscosity = false;
    bool viscosity = false;
    F viscosityCoeff = 0.1;
    F kernelSigma = 3 * radiusMean;  
    F alpha = 0.5;
    F   eta_def   = 10.0;     ///< dynamic viscosity [Pa·s]
    F   F_kernel  = 2.5;      ///< h = F_kernel * maxParticleRadius

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
    bool periodicY       = false;    // wrap Y
    F  fluidViscosity    = 10.0f;     // ν in F_drag = ν ||v_p−v_f||²
    F  V0                = 0.02f;     // max fluid speed
    F  L                 = 0.5f;     // half‐height of tunnel
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
          Young(other.Young),
          Poisson(other.Poisson),
          K(other.K),
          boolSoftDEM(other.boolSoftDEM),
          a(other.a),
          b(other.b),
          is_circular(other.is_circular),
          is_square(other.is_square),
          is_elliptical(other.is_elliptical),
          radiusMean(other.radiusMean),
          radiusStd(other.radiusStd),
          densify(other.densify),
          bidispersed(other.bidispersed),  
          radiusBig(other.radiusBig),
          radiusSmall(other.radiusSmall), 
          numBig(other.numBig),
          numSmall(other.numSmall),    
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
          deformation_viscosity(other.deformation_viscosity),
          viscosity(other.viscosity),
          viscosityCoeff(other.viscosityCoeff),
          kernelSigma(other.kernelSigma),
          alpha(other.alpha),
          eta_def(other.eta_def),
          F_kernel(other.F_kernel),
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
          periodicY(other.periodicY),
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
    F eps_n(F delta, F Lc) const;
    F depsn_dDelta(F delta, F Rbar) const;
    F depsn_dEps(F delta,F Rbar,F dRbar_dEps) const;
    F min_image(F dx, F W);
    void makeConfigMenu();
    std::vector<Particle2D> createRandomParticles2D();
    void updateEffectiveNeighborCounts();
    void updateEffectiveNeighborCountsFinal();

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
    void updateDeformationGradients();
    void updateAuxiliaryStructures();
    void colorParticleRed(int particleID);

    void startScenarioAnimation();
    void updateScenarioAnimation(F dt);

    void shearFlowProfile(F y,F& v_fx,F& dvf_dy,F& d2vf_dy2);
    F wrapX(F x) const;           // defined in .cpp
    F periodicDx(F x1, int ix1, F x2, int ix2) const;
    void renormalise();
    void updateCellSizeFromParticles();
    void minParticles();
    void computeCFL ();
    void buildMassMatrix(SparseMatrixF& M) const;
    int stride () const;
    void rebuildGlobalVectors();

};
