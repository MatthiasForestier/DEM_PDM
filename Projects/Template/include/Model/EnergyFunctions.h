#pragma once

#include "Projects/Template/include/Model/MassSpring.h"  // for Simulation, Particle2D, ScenarioObject
#include "CRLHelper/MapleHelper.h"

using AddFunc = std::function<void(int /*row*/, int /*col*/, F /*value*/)>;

namespace EnergyFunctions {

// Gravitational potential for one particle
F gravity(const Particle2D& p, const Simulation& sim, F y0);

// Volumetric strain energy for one particle
F volumetricStrain(const Particle2D& p, const Simulation& sim);

// Boundary energies (circle, square, tunnel)
F boundaryEnergy(const Particle2D& p, const Simulation& sim);

// Inter-particle contact (old DEM)
F interParticleOld(const Particle2D& p, const Particle2D& q, const Simulation& sim);

// Inter-particle contact (Soft-DEM)
F interParticleSoft(const Particle2D& p, const Particle2D& q, const Simulation& sim);

// Pin spring energy
F pinSpring(const Simulation& sim);

} 

namespace GradientFunctions
{
    // d/dy of gravitational potential (only y-component)
F gravityGradient(const Particle2D& p, const Simulation& sim);

// derivative wrt epsV of volumetric strain energy (only epsV component)
F volumetricStrainGradient(const Particle2D& p, const Simulation& sim);

// boundary force contributions: returns (fx, fy, feps)
void boundaryGradient(const Particle2D& p,
                      const Simulation& sim,
                      F& fx, F& fy, F& feps);

// inter-particle gradient (old DEM): returns (fx, fy)
std::pair<F,F> interParticleOldGradient(const Particle2D& p,
                                         const Particle2D& q,
                                         const Simulation& sim);

// inter-particle gradient (Soft-DEM): returns (fx, fy, feps_p, feps_q)
std::tuple<F,F,F,F> interParticleSoftGradient(const Particle2D& p,
                                              const Particle2D& q,
                                              const Simulation& sim);

// pin spring in x only: returns dU/dx at particle 0's x DOF
F pinSpringGradient(const Simulation& sim);

}

namespace HessianFunctions {
  
    /// boundary Hessian contributions at particle i
    void boundaryHessian(
      const Particle2D& p,
      const Simulation& sim,
      int             i,   // particle index
      AddFunc         add  // to accumulate into H
    );
  
    void interParticleOldHessian(
        const Particle2D& p,
        const Particle2D& q,
        const Simulation& sim,
        int               i,  // index of p
        int               j,  // index of q
        AddFunc           add
      );

    /// inter‐particle Hessian for pair (i<j)
    void interParticleSoftHessian(
      const Particle2D& p,
      const Particle2D& q,
      const Simulation& sim,
      int               i,  // index of p
      int               j,  // index of q
      AddFunc           add
    );
  
    /// pin spring on x for periodicX
    void pinSpringHessian(
      const Simulation& sim,
      AddFunc           add
    );
  
}

namespace EnergyDynFunctions {

    F deformationViscous(const Particle2D& p, const Simulation& sim);

}

namespace GradientDynFunctions
{
    void deformationViscous(const Particle2D& p,
                            const Simulation& sim,
                            F &fx, F &fy, F &fe);

}

namespace HessianDynFunctions {

    void deformationViscous(const Particle2D& p,
        const Simulation& sim,
        int i,
        AddFunc add);

}