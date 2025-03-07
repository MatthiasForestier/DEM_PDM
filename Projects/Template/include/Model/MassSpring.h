#pragma once

#include "CRLHelper/VecMatDef.h"
#include <vector>
#include <string>
#include <cmath>

// Particle2D is defined first.
class Simulation; // Forward declaration needed for the constructor below.

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

    F radius;
    F mass;
    F inertia; // Moment of inertia for a disc: (1/2)*mass*radius^2

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

    F radius;
    F mass;
    F inertia; // Moment of inertia for a sphere: (2/5)*mass*radius^2

    // Constructor that computes mass and inertia using the density from simParams.
    Particle3D(F radius, const Simulation& simParams);
};

// Now define Simulation. Note that Particle2D is already fully defined.
class Simulation {
public:
    // World parameters
    F timeStep = 1e-4;
    Vector3F gravity = Vector3F(0.0, -9.81, 0.0);
    F endTime = 1.0;
    F outputInterval = 1e-2;
    std::string outputDirectory = "output";
    I numParticles = 100;
    const int maxAttemptsPerParticle = 100;
    // Display dimensions commented out for now.
    // int displayWidth = 800;
    // int displayHeight = 800;
    bool use3D = false;

    // Particle-related parameters
    F density = 2780.0; // [kg/m^3] used for mass computation
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
    F mean = 0.05;                 // [m]
    F std = 0.01;                  // [m]
    // Particles
    std::vector<Particle2D> particles2D;
    std::vector<Particle3D> particles3D;

public:
    void makeConfigMenu();
    std::vector<Particle2D> createRandomParticles2D();
};

// Finally, define MassSpring.
class MassSpring {
public:
    /// Fixed parameters
    F rest_length = 0.9;
    Vector2F endpoint0 = Vector2F(-1, 0);
    Vector2F endpoint1 = Vector2F(1, 0);

    /// Degrees of freedom
    Vector2F y = Vector2F(0, 0);

public:
    void makeConfigMenu();

public:
    void compute_energy(F &value) const;
    void compute_gradient(VectorXF &gradient) const;
    void compute_hessian(SparseMatrixF &hessian) const;
};
