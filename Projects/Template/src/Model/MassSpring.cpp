#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

#include "Projects/Template/include/Model/MassSpring.h"
#include "Projects/Template/include/Model/EnergyFunctions.h"

/* ------------------------------------------------------------------ */
/*  Local helpers / constants                                          */
/* ------------------------------------------------------------------ */
namespace {

using Vec2  = Eigen::Matrix<F,2,1>;
using Vec3  = Eigen::Matrix<F,3,1>;
using Trip  = Eigen::Triplet<F>;

inline F cube(F x)                   { return x * x * x; }
inline F sqr(F x)                   { return x * x; }
constexpr F oneThird   = static_cast<F>(1.0 / 3.0);
constexpr F half       = static_cast<F>(0.5);
constexpr F two        = static_cast<F>(2.0);
constexpr F tiny       = static_cast<F>(1e-12);

/* uniform wrapper for ImGui scalar inputs -------------------------------- */
void scalarInput(const char* label, F& var, F step, const char* fmt = "%.3f")
{ ImGui::InputDouble(label, &var, step, step*5, fmt); }

} // anonymous namespace

// SDEM helpers --------------------------------------------------------------
inline F segmentLength(F Rbar, F delta)                // eq. (55) Mollon
{   return 2.0f * std::sqrt(Rbar * delta); }

F Simulation::eps_n(F delta, F Lc) const                // eq. (54)
{   return delta / (a*Lc + b*delta); }

F Simulation::depsn_dDelta(const F delta, const F Rbar) const
{
    /* ε_n(δ,R̄) = δ / ( a·L_c + b·δ ) ,  L_c = 2√(R̄ δ)  */
    const F Lc  = 2.0f * std::sqrt(Rbar * delta);
    const F den = a * Lc + b * delta;
    // ∂ε_n/∂δ = ( a √(R̄/δ) + b ) / (a L_c + b δ)^2  – simplified form:
    return a * std::sqrt(Rbar * delta) / (den * den);
}

F Simulation::depsn_dEps(const F delta,
                           const F Rbar,
                           const F dRbar_dEps) const
{
    /* chain rule: ∂ε_n/∂ε = (∂ε_n/∂δ) ∂δ/∂ε + (∂ε_n/∂R̄) ∂R̄/∂ε         */
    const F Lc   = 2.0f * std::sqrt(Rbar * delta);
    const F den  = a * Lc + b * delta;

    /* ∂ε_n/∂R̄ */
    const F dLc_dRbar  = delta / Lc;          // since L_c = 2 √(R̄ δ)
    const F dden_dRbar = a * dLc_dRbar;
    const F depsn_dRbar = -delta * dden_dRbar / (den * den);

    /* ∂δ/∂ε = 2 R₀/2 = R₀  (undeformed radius) * 1  */
    const F ddelta_dEps = dRbar_dEps * 2.0f;

    /* total derivative */
    return depsn_dDelta(delta, Rbar) * ddelta_dEps + depsn_dRbar * dRbar_dEps;
}
/* ---------------------------------------------------------------------- */
/*  CONFIG MENU                                                           */
/* ---------------------------------------------------------------------- */
void Simulation::makeConfigMenu()
{
    /* world ------------------------------------------------------------- */
    if (ImGui::CollapsingHeader("Solver / World"))
    {
        scalarInput("dt (time step)", timeStep, 1e-5, "%.6f");
        scalarInput("Gravity y"      , gravity(1), 0.1, "%.2f");
    }

    /* particles --------------------------------------------------------- */
    if (ImGui::CollapsingHeader("Particles"))
    {
        ImGui::InputInt("Count", &numParticles);
        scalarInput("Radius mean", radiusMean, 0.01, "%.3f");
        scalarInput("Radius std", radiusStd , 0.01, "%.4f");
    }

    /* contact / material ----------------------------------------------- */
    if (ImGui::CollapsingHeader("Material & Contact"))
    {
        scalarInput("Density"          , density          , 10.0);
        scalarInput("Wall stiff k"       , overlapParam     , 10.0);
        scalarInput("Particle stiff k"   , interactionParam , 10.0);
        scalarInput("Viscosity c"        , viscosityCoeff   , 0.01, "%.4f");
        scalarInput("Kernel sigma"           , kernelSigma      , 0.01, "%.4f");
        scalarInput("Young modulus" , Young   , 0.01, "%.4f");
        scalarInput("Poisson ratio" , Poisson , 0.01, "%.4f");
        K = Young / (3*(1-2*Poisson));

        bool old = boolSoftDEM;
        ImGui::Checkbox("Soft-DEM (1 DoF/particle)", &boolSoftDEM);

        /* if the user toggled, rebuild the global state */
        if (old != boolSoftDEM)
            rebuildGlobalVectors();   
            // std::cout << K << std::endl;       
    }
}



/* ---------------------------------------------------------------------- */
/*  GRID / NEIGHBOUR STRUCTURES                                           */
/* ---------------------------------------------------------------------- */
void Simulation::updateCellSizeFromParticles()
{
    maxParticleRadius = F(0);
    for (auto& p : particles2D) maxParticleRadius = std::max(maxParticleRadius, p.radius);
    cellSize = two * maxParticleRadius;
}

void Simulation::buildGridDataStructure()
{
    if (scenarioObjects.empty())
    {
        std::cerr << "[Grid] No scenario objects – grid not built.\n";
        return;
    }

    updateCellSizeFromParticles();

    const BoundingBox& bb = scenarioObjects.front()->getBoundingBox();
    minX = bb.min_x; maxX = bb.max_x;
    minY = bb.min_y; maxY = bb.max_y;

    numCellsX = int(std::ceil((maxX-minX)/cellSize));
    numCellsY = int(std::ceil((maxY-minY)/cellSize));

    grid.assign(numCellsX * numCellsY, {});
}

void Simulation::minParticles() {
    for (auto &p : particles2D)
        minParticleDiam = std::min(minParticleDiam, p.radius);
    // pick whatever factor keeps each disc safely within one neighbour cell
    minParticleDiam = 2*minParticleDiam;
}

// Updated: Insert Particle2D objects into the grid.
// Now each particle is inserted into all cells that its circle overlaps.
void Simulation::insertParticlesIntoGrid()
{
    // 1) Clear all cells
    for (auto &cell : grid)
        cell.clear();

    // 2) Standard insertion (possibly wrapping in X)
    for (I i = 0; i < (I)particles2D.size(); ++i) {
        const Particle2D &p = particles2D[i];
        // compute world‐space AABB of this particle
        F x0 = p.pos(0) - p.effectiveRadius();
        F x1 = p.pos(0) + p.effectiveRadius();
        F y0 = p.pos(1) - p.effectiveRadius();
        F y1 = p.pos(1) + p.effectiveRadius();

        int minCX = (int)std::floor((x0 - minX) / cellSize);
        int maxCX = (int)std::floor((x1 - minX) / cellSize);
        int minCY = (int)std::floor((y0 - minY) / cellSize);
        int maxCY = (int)std::floor((y1 - minY) / cellSize);

        for (int cy = minCY; cy <= maxCY; ++cy) {
            if (cy < 0 || cy >= numCellsY) continue;
            for (int cxRaw = minCX; cxRaw <= maxCX; ++cxRaw) {
                int cxWrapped;
                if (periodicX) {
                    // wrap X
                    cxWrapped = ((cxRaw % numCellsX) + numCellsX) % numCellsX;
                } else {
                    if (cxRaw < 0 || cxRaw >= numCellsX) continue;
                    cxWrapped = cxRaw;
                }
                grid[cy * numCellsX + cxWrapped].push_back(i);
            }
        }
    }

    // 3) **Only if** we're in periodic‑X mode, merge the two edge columns
    if (periodicX && numCellsX > 1) {
        for (int cy = 0; cy < numCellsY; ++cy) {
            int idxL = cy * numCellsX + 0;
            int idxR = cy * numCellsX + (numCellsX - 1);

            auto &leftCell  = grid[idxL];
            auto &rightCell = grid[idxR];

            // append right‐cell contents into left‐cell
            leftCell.insert(leftCell.end(),
                            rightCell.begin(),
                            rightCell.end());

            // (optional) keep symmetry so neighbor queries on the right edge see the same list
            rightCell = leftCell;
        }
    }
}



void Simulation::updateNeighborLists()
{
    /* clear old lists */
    for (auto& p : particles2D)
        p.neighborIndices.clear();

    /* constants */
    const int W = numCellsX;                // horizontal cell count

    /* -------- one pass per particle -------- */
    for (size_t i = 0; i < particles2D.size(); ++i)
    {
        Particle2D& p = particles2D[i];

        /* bounding box in *cell* coordinates – NOT clamped */
        F cxWrapped  = periodicX ? wrapX(p.pos(0)) : p.pos(0);
        int cxMinRaw = (int)std::floor((cxWrapped - p.effectiveRadius() - minX) / cellSize);
        int cxMaxRaw = (int)std::floor((cxWrapped + p.effectiveRadius() - minX) / cellSize);
        int cyMin    = std::max(0,  (int)std::floor((p.pos(1)-p.effectiveRadius() - minY)/cellSize));
        int cyMax    = std::min(numCellsY-1,
                                (int)std::floor((p.pos(1)+p.effectiveRadius() - minY)/cellSize));

        /* gather candidate indices (avoid duplicates with a set) */
        std::unordered_set<int> cand;
        for (int cy = cyMin; cy <= cyMax; ++cy)
        {
            for (int cxRaw = cxMinRaw; cxRaw <= cxMaxRaw; ++cxRaw)
            {
                int cx = cxRaw;

                /* wrap or reject in X */
                if (periodicX)
                    cx = ((cx % W) + W) % W;          // modulo wrap
                else if (cx < 0 || cx >= W)
                    continue;                         // outside → skip

                int gradient = cy * W + cx;                  // 1‑D cell index
                for (int j : grid[gradient])
                    if (j != (int)i) cand.insert(j);
            }
        }

        /* real‑overlap filtering (and symmetry) */
        for (int j : cand)
        {
            if (j <= (int)i) continue;               // keep ≤ once
            const Particle2D& q = particles2D[j];

            F dx = periodicDx(p.pos(0), p.ix, q.pos(0), q.ix);
            F dy = p.pos(1) - q.pos(1);
            if (std::sqrt(dx*dx + dy*dy) < p.effectiveRadius() + q.effectiveRadius())
                p.neighborIndices.push_back(j);
        }
    }
}

void Simulation::updateDeformationGradients()
{
    if (!deformation_viscosity) return;          // fast exit if disabled

    const F h2 = std::pow(F_kernel * maxParticleRadius, 2);
    const F inv_dt = 1.0 / timeStep;

    for (std::size_t i = 0; i < particles2D.size(); ++i)
    {
        auto &p = particles2D[i];

        Matrix2F Apq = Matrix2F::Zero();
        Matrix2F Aqq = Matrix2F::Zero();

        /* accumulate weighted moments over neighbours of i */
        for (I j : p.neighborIndices)
        {
            const auto &q = particles2D[j];

            Vector2F Xij = q.X0 - p.X0;         // rest‑space offset
            Vector2F xij = q.pos - p.pos;       // current offset

            F w = std::exp( -xij.squaredNorm() / h2 );

            Apq += w * (xij * Xij.transpose());
            Aqq += w * (Xij * Xij.transpose());
        }

        /* small regulariser for robustness */
        const F eps = 1e-6 * Aqq.trace();
        Aqq(0,0) += eps;  Aqq(1,1) += eps;

        Matrix2F F_now = Apq * Aqq.inverse();

        /* rate‑of‑deformation tensor */
        Matrix2F Fdot = (F_now - p.F_prev) * inv_dt;
        p.D  = 0.5 * (Fdot + Fdot.transpose());
        p.F_prev = F_now;                        // roll to next step
    }
}

/* ------------------------------------------------------------------ */
/*  write particle data → globalPositions                             */
/* ------------------------------------------------------------------ */
void Simulation::updateGlobalPositions()
{
    int n = use3D ? int(particles3D.size())
                  : int(particles2D.size());

    const int S = stride();                     // 2 or 3
    globalPositions.resize(S * n);

    if (use3D) {
        const int n = static_cast<int>(particles3D.size());
        globalPositions.resize(DOF_FULL * n);

        for (int i = 0; i < n; ++i) {
            globalPositions.segment<DOF_FULL>(DOF_FULL * i) = particles3D[i].pos;
        }
        return;
    }
    else {
        for (int i = 0; i < n; ++i)
        {
            const Particle2D &p = particles2D[i];

            /* x , y */
            globalPositions.segment( S*i, DOF) = p.pos;

            /* optional ε_V */
            if (boolSoftDEM)
                globalPositions[S*i + 2] = p.epsV;
        }
    }
}

/* ------------------------------------------------------------------ */
/*  read → particles                                                  */
/* ------------------------------------------------------------------ */
void Simulation::applyGlobalPositions(const VectorXF &P)
{
    if (use3D) {
        const int n = static_cast<int>(particles3D.size());
        for (int i = 0; i < n; ++i)
            particles3D[i].pos = P.segment<DOF_FULL>(DOF_FULL * i);
        return;
    }
    else {
        int n = int(particles2D.size());
        const int S = stride();

        for (int i = 0; i < n; ++i)
        {
            Particle2D &p = particles2D[i];

            /* x , y (always) */
            p.pos = P.segment( S*i, DOF);

            if (periodicX)
                p.pos(0) = wrapX(p.pos(0));
                //globalPositions[S*i] = p.pos(0);

            /* ε_V if present */
            if (boolSoftDEM)
                p.epsV = P[S*i + 2];
            else
                p.epsV = 0.0;                     // muted
        }
    }
}


VectorXF Simulation::getGlobalState() {
    updateGlobalPositions();
    return globalPositions;
}

void Simulation::setGlobalState(const VectorXF &state) {
    globalPositions = state;
    applyGlobalPositions(globalPositions); 
    //renormalise();
}

void Simulation::updateAuxiliaryStructures() {
    // Build the grid and update neighbor lists based on the current particle positions.
    insertParticlesIntoGrid();
    updateNeighborLists();
    // updateDeformationGradients();
}

void Simulation::rebuildGlobalVectors()
{
    /* ------------------------------------------------------------------ */
    /* 0) cache the old global state (whatever layout it had)             */
    /* ------------------------------------------------------------------ */
    VectorXF oldState = getGlobalState();          // uses *old* stride
    const int Sold    = boolSoftDEM ? 2 : 3;       // old stride (before toggle)

    /* ------------------------------------------------------------------ */
    /* 1) update stride after the flag has been flipped outside           */
    /* ------------------------------------------------------------------ */
    const int Snew = stride();                     // 2 or 3

    /* ------------------------------------------------------------------ */
    /* 2) resize globalPositions and repack x , y ( + ε_V if active )     */
    /* ------------------------------------------------------------------ */
    const int n = static_cast<int>(particles2D.size());
    globalPositions.setZero(Snew * n);

    for (int i = 0; i < n; ++i)
    {
        /* copy back x , y ------------------------------------------------ */
        globalPositions.segment<2>(Snew * i) =            // always 2 entries
            oldState.segment<2>(Sold * i);

        /* set / reset ε_V ------------------------------------------------- */
        if (boolSoftDEM) {
            //  old version may or may not have ε_V – treat missing as 0
            F epsOld = (Sold > 2) ? oldState(Sold * i + 2) : F(0);
            particles2D[i].epsV = epsOld;                 // keep continuity
            globalPositions(Snew * i + 2) = epsOld;
        } else {
            particles2D[i].epsV = F(0);                   // mute deformation
        }
    }

    /* ------------------------------------------------------------------ */
    /* 3) history buffers, mass-matrix, neighbour structures, CFL …       */
    /* ------------------------------------------------------------------ */
    globalState_1 = globalPositions;
    globalState_2 = globalPositions;

    /* mass matrix (diagonal) */
    buildMassMatrix(M);

    /* neighbour lists depend only on geometry (unchanged) but            */
    /* grid cell width uses the new *physical* radii → recompute          */
    updateCellSizeFromParticles();
    buildGridDataStructure();
    insertParticlesIntoGrid();
    updateNeighborLists();

    /*  keep CFL diagnostic up-to-date                                    */
    minParticles();
}


/* -------------------------------------------------------------------------- */
/*  ENERGY                                                                    */
/* -------------------------------------------------------------------------- */
void Simulation::compute_energy(F &value) const
{
    value = F(0);

    for (std::size_t ii = 0; ii < particles2D.size(); ++ii)
    {
        /* ------------------------------------------------------------------ */
        /*  particle → local copy (for collision tagging)                     */
        /* ------------------------------------------------------------------ */
        Particle2D p = particles2D[ii];
        detectBoundaryCollision2D(p);

        /* ---------------- gravitational potential ------------------------- */
        const F y0 = scenarioObjects.empty()
                    ? F(0)                       /* fallback datum       */
                    : scenarioObjects[0]->getBoundingBox().min_y;

        value += EnergyFunctions::gravity(p, *this, y0);
        
        /* volumetric strain energy ----------------- */
        if (boolSoftDEM)
        { 
            value += EnergyFunctions::volumetricStrain(p, *this);
        }
        /* ----------------- boundary energy ----------------- */
        value += EnergyFunctions::boundaryEnergy(p, *this);
    }

    /*  -------------inter-particle contacts--------------------------------*/
    for (std::size_t i = 0; i < particles2D.size(); ++i)
    {
        const Particle2D &p = particles2D[i];
        for (int j : p.neighborIndices)
        {
            if (j <= (int)i) continue;
            const Particle2D &q = particles2D[j];

            if (!boolSoftDEM){
                value += EnergyFunctions::interParticleOld(p, q, *this);
            }else{
                value += EnergyFunctions::interParticleSoft(p, q, *this);
            }
        }
    }

    /* ------------- pin spring (prevent rigid motion) -------------------- */
    if (periodicX && !particles2D.empty())
    {
        value += EnergyFunctions::pinSpring(*this);
    }
}

/* -------------------------------------------------------------------------- */
/*  GRADIENT                                                                  */
/* -------------------------------------------------------------------------- */

void Simulation::compute_gradient(VectorXF &g) const {
    const int S = stride();
    g.resize(S * particles2D.size());
    g.setZero();

    // per-particle terms
    for (size_t ii = 0; ii < particles2D.size(); ++ii) {
        Particle2D p = particles2D[ii];
        detectBoundaryCollision2D(p);

        int xIdx = S*ii, yIdx = xIdx+1, epsIdx = xIdx+2;
        // gravity
        g[yIdx] += GradientFunctions::gravityGradient(p, *this);

        // volumetric
        if (boolSoftDEM)
            g[epsIdx] += GradientFunctions::volumetricStrainGradient(p, *this);

        // boundary
        F fx, fy, feps;
        GradientFunctions::boundaryGradient(p, *this, fx, fy, feps);
        g[xIdx]   += fx;
        g[yIdx]   += fy;
        if (boolSoftDEM) g[epsIdx] += feps;
    }

    // inter-particle
    for (size_t i = 0; i < particles2D.size(); ++i) {
        auto &p = particles2D[i];
        int ix = S*i, iy = ix+1, ie = ix+2;

        for (int j : p.neighborIndices) {
            if (j <= (int)i) continue;
            auto &q = particles2D[j];
            int jx = S*j, jy = jx+1, je = jx+2;

            if (!boolSoftDEM) {
                auto [fx,fy] = GradientFunctions::interParticleOldGradient(p, q, *this);
                g[ix] += fx;  g[iy] += fy;
                g[jx] -= fx;  g[jy] -= fy;
            } else {
                auto [fx, fy, feps_i, feps_j]
                    = GradientFunctions::interParticleSoftGradient(p, q, *this);
                g[ix] += fx;    g[iy] += fy;
                g[jx] -= fx;    g[jy] -= fy;
                g[ie] += feps_i;
                g[je] += feps_j;
            }
        }
    }

    // pin spring
    if (periodicX && !particles2D.empty()) {
        g[0] += GradientFunctions::pinSpringGradient(*this);
    }
}

/* -------------------------------------------------------------------------- */
/*  HESSIAN                                                                   */
/*  (sparse COO additions – same pattern everywhere)                          */
/* -------------------------------------------------------------------------- */
void Simulation::compute_hessian(SparseMatrixF &H) const
{
  const int S = stride();
  const int N = S * (int)particles2D.size();
  H.resize(N,N);
  H.setZero();

  auto add = [&](int r,int c,F v){
    H.coeffRef(r,c) += v;
  };

  // 1) boundary walls
  for (int i = 0; i < (int)particles2D.size(); ++i) {
    Particle2D p = particles2D[i];
    detectBoundaryCollision2D(p);
    HessianFunctions::boundaryHessian(p, *this, i, add);
  }

  // 2) inter‐particle
  for (int i = 0; i < (int)particles2D.size(); ++i) {
    for (int j : particles2D[i].neighborIndices) {
      if (j <= i) continue;
      if (!boolSoftDEM){
        HessianFunctions::interParticleOldHessian(
            particles2D[i], particles2D[j], *this, i, j, add);
      }else{
        HessianFunctions::interParticleSoftHessian(
            particles2D[i], particles2D[j], *this, i, j, add);
      }
    }
  }

  // 3) pin spring
  HessianFunctions::pinSpringHessian(*this, add);
}

/* ======================================================================= */
/*  DYNAMIC  ENERGY                                                        */
/* ======================================================================= */
void Simulation::compute_energy_dyn(F& value)
{
    /* --- static part ---------------------------------------------------- */
    compute_energy(value);

    /* -------------------------------------------------------------------- */
    /*  make sure history vectors have the right size                       */
    /* -------------------------------------------------------------------- */
    if (globalState_1.size() != globalPositions.size()) {
        globalState_1 = globalPositions;
        globalState_2 = globalPositions;
    }

    /* Δ²-like difference  a = xⁿ⁺¹ − 2xⁿ + xⁿ⁻¹ ------------------------- */
    const VectorXF a = globalPositions
                     - 2 * globalState_1
                     + globalState_2;

    const F inv_h2 = 1.0 / (timeStep * timeStep);

    /* ½/h² · aᵀ M a  (M already includes the ε_V row when Soft-DEM on)   */
    value += half * inv_h2 * a.dot(M * a);

    /* ---------------------------------------------------------------- */
    /*  neighbour-dependent dash-pot (x,y only, same as before)         */
    /* ---------------------------------------------------------------- */
    if (viscosity) {
        const int S = stride();                 // 2 or 3
        for (size_t i = 0; i < particles2D.size(); ++i) {
            F nEff = particles2D[i].prevEffectiveCount;

            VectorXF diff = globalPositions.segment(S * i, DOF)      // x,y
                          - globalState_1.segment(S * i, DOF);

            value += half * viscosityCoeff * nEff * inv_h2 * diff.squaredNorm();
        }
    }

    if (deformation_viscosity) {
        for (auto &p : particles2D)
            value += EnergyDynFunctions::deformationViscous(p, *this);
    }
    

    /* ----------------------------- shear-flow drag -------------------- */
    if (experiment == Experiment::ShearFlow) {
        for (int i = 0; i < (int)particles2D.size(); ++i) {
            const int base = stride() * i;          // works for 2 or 3 DOF

            F x1 = globalPositions(base    );
            F y1 = globalPositions(base + 1);
            F x0 = globalState_1  (base    );
            F y0 = globalState_1  (base + 1);

            F vfx , dvf_dy , dummy;
            shearFlowProfile(y1, vfx, dvf_dy, dummy);

            F dvx = (x1 - x0) / timeStep - vfx;
            F dvy = (y1 - y0) / timeStep;

            value += half * fluidViscosity * (dvx*dvx + dvy*dvy);
        }
    }
}


/* ======================================================================= */
/*  DYNAMIC  GRADIENT                                                      */
/* ======================================================================= */
void Simulation::compute_gradient_dyn(VectorXF& grad)
{
    /* static part -------------------------------------------------------- */
    compute_gradient(grad);

    if (globalState_1.size() != globalPositions.size()) {
        globalState_1 = globalPositions;
        globalState_2 = globalPositions;
    }

    const VectorXF a      = globalPositions
                          - 2 * globalState_1
                          + globalState_2;
    const F inv_h2 = 1.0 / (timeStep * timeStep);

    /* 1/h² · M a --------------------------------------------------------- */
    grad += inv_h2 * (M * a);

    /* neighbour-dependent dash-pot (x,y only) --------------------------- */
    if (viscosity) {
        const int S = stride();
        for (size_t i = 0; i < particles2D.size(); ++i) {
            F nEff = particles2D[i].prevEffectiveCount;
            grad.segment(S * i, DOF) += viscosityCoeff * nEff * inv_h2 *
                                         (globalPositions.segment(S * i, DOF) -
                                          globalState_1  .segment(S * i, DOF));
        }
    }

    if (deformation_viscosity) {
        if (deformation_viscosity) {
            const int S = stride();
            for (int i = 0; i < (int)particles2D.size(); ++i) {
              F fx, fy, fe;
              GradientDynFunctions::deformationViscous(particles2D[i], *this, fx, fy, fe);
              grad[S*i+0] += fx;
              grad[S*i+1] += fy;
              if (S==3)    grad[S*i+2] += fe;
            }
        }
    }

    /* shear-flow drag ---------------------------------------------------- */
    if (experiment == Experiment::ShearFlow) {
        const F nu      = fluidViscosity;
        const F inv_dt  = 1.0 / timeStep;
        const int S     = stride();

        for (int i = 0; i < (int)particles2D.size(); ++i) {
            const int base = S * i;

            F x1 = globalPositions(base    );
            F y1 = globalPositions(base + 1);
            F x0 = globalState_1  (base    );
            F y0 = globalState_1  (base + 1);

            F vfx , dvf_dy , dummy;
            shearFlowProfile(y1, vfx, dvf_dy, dummy);

            F dvx = (x1 - x0)*inv_dt - vfx;
            F dvy = (y1 - y0)*inv_dt;

            grad(base    ) += nu * dvx * inv_dt;             // ∂E/∂x
            grad(base + 1) += nu * dvy * inv_dt
                             - nu * dvx * dvf_dy;            // ∂E/∂y
        }
    }
}


/* ======================================================================= */
/*  DYNAMIC  HESSIAN                                                       */
/* ======================================================================= */
void Simulation::compute_hessian_dyn(SparseMatrixF& H)
{
    /* static part -------------------------------------------------------- */
    compute_hessian(H);

    const F inv_h2 = 1.0 / (timeStep * timeStep);

    /* add 1/h² · M (diagonal) ------------------------------------------- */
    H += inv_h2 * M;

    /* neighbour-dependent dash-pot (diag on x,y) ------------------------ */
    if (viscosity) {
        const int S = stride();
        for (size_t i = 0; i < particles2D.size(); ++i) {
            F nEff = particles2D[i].prevEffectiveCount;
            const F coeff = viscosityCoeff * nEff * inv_h2;

            for (int d = 0; d < DOF; ++d) {          // x,y only
                int id = int(S * i + d);
                H.coeffRef(id, id) += coeff;
            }
        }
    }


    if (deformation_viscosity) {
        for (int i = 0; i < (int)particles2D.size(); ++i)
            HessianDynFunctions::deformationViscous(particles2D[i], *this, i,
                                                    [&](int r,int c,F v){ H.coeffRef(r,c) += v; });
    }

    /* shear-flow contribution ------------------------------------------- */
    if (experiment == Experiment::ShearFlow) {
        const F nu     = fluidViscosity;
        const F inv_dt = 1.0 / timeStep;
        const F k      = M_PI / L;
        const int S    = stride();

        for (int i = 0; i < (int)particles2D.size(); ++i) {
            const int base = S * i;

            F x1 = globalPositions(base    );
            F y1 = globalPositions(base + 1);
            F x0 = globalState_1  (base    );

            /* v_f and derivatives */
            F vfx , dvf_dy , d2vf_dy2;
            shearFlowProfile(y1, vfx, dvf_dy, d2vf_dy2);

            F sin_k_y = std::sin(k * y1);
            F cos_k_y = std::cos(k * y1);
            F dvx     = (x1 - x0)*inv_dt - V0*sin_k_y;

            /* diagonal -------------------------------------------------- */
            H.coeffRef(base, base) += nu * inv_dt * inv_dt;          // d²E/dx²

            F term1 = nu * inv_dt * inv_dt;                          // from v_p,y
            F term2 = nu * (V0*V0) * k*k * cos_k_y * cos_k_y;        // ∂dvx/∂y
            F term3 = nu * dvx * V0 * k*k * sin_k_y;                 // 2nd deriv
            H.coeffRef(base + 1, base + 1) += term1 + term2 + term3;

            /* off-diag (symmetric) ------------------------------------ */
            F off = -nu * V0 * k * cos_k_y * inv_dt;
            H.coeffRef(base,     base + 1) += off;
            H.coeffRef(base + 1, base    ) += off;
        }
    }
}


void Simulation::shearFlowProfile(F y,
    F& v_fx,      //  sin(π y/L)
    F& dvf_dy,    //  (π/L) cos(π y/L)
    F& d2vf_dy2)  // −(π/L)^2 sin(π y/L)
{
    const F k = M_PI / L;            // π / L
    v_fx     = V0 * std::sin(k * y);
    F c      = std::cos(k * y);
    F s      = std::sin(k * y);
    dvf_dy   = V0 * k * c;
    d2vf_dy2 = -V0 * k * k * s;
}

F Simulation::wrapX(F x) const {
    const F W = maxX - minX;
    // std::cout << W << std::endl;
    // If the domain has not been initialised yet (W == 0)
    // simply return the original coordinate to avoid Inf / NaN.
    if (W == 0 || !periodicX)
        return x;

    return x - std::floor((x - minX) / W) * W;
}

F Simulation::periodicDx(F x1,int ix1, F x2,int ix2) const
{
    // const F W = maxX - minX;
    // F dx = x1 - x2;
    // dx = dx - std::round(dx / W) * W;   // W = maxX-minX
    // return dx;

    if (!periodicX) return x1 - x2;

    const F W = maxX - minX;
    if (W <= std::numeric_limits<F>::epsilon())
        return x1 - x2;

    /* true signed separation in the covering space */
    F dx = (x1 - x2) + (F)(ix1 - ix2) * W;

    /* wrap onto (‑½W , ½W]  –  round() gives the nearest integer */
    dx -= std::round(dx / W) * W;
    return dx;
}

void Simulation::renormalise()
{
    const F W = maxX - minX;
    const I S = stride();
    for (int i = 0; i < (int)particles2D.size(); ++i) {
        auto &p = particles2D[i];

        /* old wrap counter */
        int old_ix = p.ix;

        /* keep pos.x in [minX,maxX) and update p.ix */
        while (p.pos(0) <  minX) { p.pos(0) += W; --p.ix; }
        while (p.pos(0) >= maxX) { p.pos(0) -= W; ++p.ix; }

        /* Δix since last frame */
        int dix = p.ix - old_ix;
        if (dix != 0) {
            F shift = (F)dix * W;
            globalState_1(S*i    ) += shift;
            globalState_2(S*i    ) += shift;
        }
    }
}

void Simulation::updateEffectiveNeighborCounts() {
    // Loop over each particle.
    for (size_t i = 0; i < particles2D.size(); i++) {
        F rawEffectiveCount = 0;
        const int S = stride();
        // Fetch mod‑position + wrap counter for convenience
        const F xi   = globalPositions(S*i    );  // already ∈ [minX,maxX)
        const int ixi = particles2D[i].ix;

        // Compute the raw effective count using a soft kernel w(d) = exp(−d²/σ²)
        for (int j : particles2D[i].neighborIndices) {
            // x‑difference with correct periodic image
            const F xj   = globalPositions(S*j    );
            const int ixj = particles2D[j].ix;
            F dx = periodicDx(xi, ixi, xj, ixj);

            // y is non‑periodic
            F dy = globalPositions(S*i + 1)
                 - globalPositions(S*j + 1);

            F d2 = dx*dx + dy*dy;  // ignore θ
            rawEffectiveCount += std::exp(-d2 / (kernelSigma * kernelSigma));
        }

        // Store the raw count for this frame
        particles2D[i].effectiveCountCurrent = rawEffectiveCount;
        
        // Temporal filtering
        F filteredCount = alpha * rawEffectiveCount
                        + (1 - alpha) * particles2D[i].prevEffectiveCount;
        // Clamp to avoid runaway stiffness
        filteredCount = std::min(filteredCount, F(1.0));

        // **Remember to write it back** so next frame can filter against it
        particles2D[i].prevEffectiveCount = filteredCount;
    }
}

void Simulation::updateEffectiveNeighborCountsFinal() {
    // Loop over each particle.
    for (size_t i = 0; i < particles2D.size(); i++) {
        F rawEffectiveCount = 0;

        const int S = stride();
        // Cached mod‐position + wrap counter
        const F  xi   = globalPositions(S*i    );  // ∈ [minX,maxX)
        const int ixi = particles2D[i].ix;

        // Soft‐kernel sum over neighbours
        for (int j : particles2D[i].neighborIndices) {
            const F  xj   = globalPositions(S*j    );
            const int ixj = particles2D[j].ix;

            F dx = periodicDx(xi, ixi, xj, ixj);
            F dy = globalPositions(S*i + 1)
                 - globalPositions(S*j + 1);

            F d2 = dx*dx + dy*dy;
            rawEffectiveCount += std::exp(-d2 / (kernelSigma * kernelSigma));
        }

        // Store the raw count
        particles2D[i].effectiveCountCurrent = rawEffectiveCount;
        
        // Temporal filtering
        F filteredCount = alpha * rawEffectiveCount
                        + (1 - alpha) * particles2D[i].prevEffectiveCount;
        filteredCount = std::min(filteredCount, F(1.0));

        // Write it back for the next frame
        particles2D[i].prevEffectiveCount = filteredCount;
    }
}

std::vector<Particle2D> Simulation::createRandomParticles2D()
{
    /* ------------------------------------------------------------ 0. Parameters */
    constexpr F SHRINK        = 0.60f;      // initial radius scale (0–1)
    constexpr int MAX_RELAX   = 180;        // number of growth sweeps
    constexpr F STEP_FRAC     = 0.52f;      // ≥0.5 → resolve pair overlap in one push
    constexpr F TOL           = 1e-10f;     // stop when max overlap < TOL·R_mean
    constexpr F JITTER        = 0.05f;      // random ±5 % R_mean offset
    constexpr F MAX_PACK_FRAC = 0.88f;      // safety: abort if requested area exceeds 88 % of domain

    /* ------------------------------------------------------------ 1. Bounds */
    F xMin = -0.8f, xMax = 0.8f;
    F yMin = -0.5f, yMax = 0.5f;

    if (!scenarioObjects.empty()) {
        if (auto* tun = dynamic_cast<Tunnel2D*>(scenarioObjects[0].get())) {
            const BoundingBox& BB = tun->getBoundingBox();
            xMin = BB.min_x; 
            xMax = BB.max_x;
            yMin = BB.min_y; 
            yMax = BB.max_y;
        }
    }

    /* domain area (needed for safety check) */
    const F domainArea = (xMax - xMin) * (yMax - yMin);

    /* ------------------------------------------------------------ 2. RNG */
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<F>  rDist(radiusMean, radiusStd);
    std::uniform_real_distribution<F> uni(-1.0f, 1.0f);

    /* quick over‑crowding check for RSA branch (approximate) */
    {
        const F probeArea = static_cast<F>(M_PI) * radiusMean * radiusMean;
        if (!densify && static_cast<F>(numParticles) * probeArea > MAX_PACK_FRAC * domainArea) {
            std::cerr << "[RSA] Too many particles for the available space (area ratio > "
                      << MAX_PACK_FRAC << "). Aborting.\n";
            return {};
        }
    }

    /* plain RSA branch (unchanged) */
    if (!densify) {
        std::vector<Particle2D> particles;
        particles.reserve(numParticles);
        std::uniform_real_distribution<F> distX(xMin, xMax);

        for (int i = 0; i < numParticles; ++i) {
            bool ok = false;
            int attempts = 0;
            Particle2D cand(0.0f, *this);

            while (!ok && attempts < maxAttemptsPerParticle) {
                ++attempts;
                F r = std::max<F>(F(0.01f), rDist(gen));
                cand = Particle2D(r, *this);
                cand.ix = 0;

                std::uniform_real_distribution<F> distY(yMin + r, yMax - r);
                cand.pos << distX(gen), distY(gen);
                cand.pos(0) = std::clamp(cand.pos(0), xMin, std::nextafter(xMax, xMin));

                ok = true;
                for (auto& ex : particles) {
                    F dx = periodicDx(cand.pos(0), cand.ix, ex.pos(0), ex.ix);
                    F dy = cand.pos(1) - ex.pos(1);
                    if (std::sqrt(dx*dx + dy*dy) < cand.radius + ex.radius) {
                        ok = false;
                        break;
                    }
                }
            }
            if (ok) {
                cand.X0 = cand.pos;
                particles.push_back(std::move(cand));
            } else {
                std::cerr << "[RSA] could not place #" << i + 1 << "\n";
            }
        }
        return particles;
    }
    // std::cout << "Dense" << std::endl;
    /* dense mode (densify == true) */
    auto clampX = [&](Particle2D& p) {
        if (periodicX) {
            p.pos(0) = wrapX(p.pos(0));
        } else {
            F r = p.radius;
            p.pos(0) = std::clamp(p.pos(0), xMin + r, xMax - r);
        }
    };

    /* sample target radii */
    std::vector<F> R_target(numParticles);
    for (F& r : R_target) {
        r = std::max<F>(F(0.01f), rDist(gen));
    }

    /* safety: total area check */
    F totalArea = 0.0f;
    const F PI = std::acos(F(-1));
    for (const F r : R_target) {
        totalArea += PI * r * r;
    }
    if (totalArea > MAX_PACK_FRAC * domainArea) {
        std::cerr << "[densify] Requested " << numParticles
                  << " particles cannot fit (area ratio = " << (totalArea / domainArea)
                  << " > " << MAX_PACK_FRAC << "). Aborting.\n";
        return {};
    }

    /* hex seed */
    const F a  = 2.0f * radiusMean * SHRINK;
    const F ay = a * std::sqrt(3.0f) / 2.0f;
    std::vector<Particle2D> discs;
    discs.reserve(numParticles);

    int row = 0;
    for (F y = yMin + ay; y < yMax - ay && discs.size() < static_cast<size_t>(numParticles);
         y += ay, ++row)
    {
        bool odd = row & 1;
        for (F x = xMin + (odd ? a * 0.5f : a);
             x < xMax - a && discs.size() < static_cast<size_t>(numParticles);
             x += a)
        {
            Particle2D p(R_target[discs.size()] * SHRINK, *this);
            p.ix = 0;
            p.pos << x, y;
            p.X0 = p.pos;
            clampX(p);
            discs.push_back(std::move(p));
        }
    }

    /* top‑up if lattice too small */
    std::uniform_real_distribution<F> distX(xMin, xMax);
    int topupAttempts = 0;
    const int MAX_TOPUP = 50 * numParticles;
    while (discs.size() < static_cast<size_t>(numParticles) &&
           topupAttempts < MAX_TOPUP)
    {
        ++topupAttempts;
        F r = R_target[discs.size()] * SHRINK;
        Particle2D p(r, *this);
        std::uniform_real_distribution<F> distY(yMin + r, yMax - r);
        p.ix = 0;
        p.pos << distX(gen), distY(gen);
        p.X0 = p.pos;
        clampX(p);
        discs.push_back(std::move(p));
        // std::cout << p.pos(0) << std::endl;
    }
    if (discs.size() < static_cast<size_t>(numParticles)) {
        std::cerr << "[densify] Could seed only " << discs.size()
                  << " / " << numParticles << " discs. Aborting.\n";
        return {};
    }
    /* jitter */
    for (auto& p : discs) {
        p.pos(0) += uni(gen) * JITTER * radiusMean;
        p.pos(1) += uni(gen) * JITTER * radiusMean;
        p.pos(1) = std::clamp(p.pos(1), yMin + p.radius, yMax - p.radius);
        clampX(p);
        // std::cout << p.pos(0) << std::endl;
    }

    /* wall relax helper */
    auto wallRelax = [&](std::vector<Particle2D>& ps) {
        for (auto& p : ps) {
            F pen = (yMin + p.radius) - p.pos(1);
            if (pen > 0) p.pos(1) += pen + TOL * radiusMean;
            pen = p.pos(1) - (yMax - p.radius);
            if (pen > 0) p.pos(1) -= pen + TOL * radiusMean;
            clampX(p);
        }
    };

    /* growth + relax */
    for (int sweep = 0; sweep < MAX_RELAX; ++sweep) {
        F g = std::pow((F)(sweep + 1) / (F)MAX_RELAX, 3.0f);
        for (size_t i = 0; i < discs.size(); ++i) {
            discs[i].radius = SHRINK * (1.0f - g) * R_target[i] + g * R_target[i];
        }
        wallRelax(discs);

        F maxOv;
        int innerIter = 0;
        const int INNER_LIMIT = 10000;
        do {
            if (++innerIter > INNER_LIMIT) {
                std::cerr << "[densify] Relaxation stalled (" << innerIter << "). Aborting.\n";
                return {};
            }
            maxOv = 0.0f;
            std::shuffle(discs.begin(), discs.end(), gen);
            for (size_t i = 0; i < discs.size(); ++i)
            for (size_t j = i + 1; j < discs.size(); ++j) {
                auto& a = discs[i];
                auto& b = discs[j];
                F dx = periodicDx(a.pos(0), a.ix, b.pos(0), b.ix);
                F dy = a.pos(1) - b.pos(1);
                F d2 = dx*dx + dy*dy;
                F rSum = a.radius + b.radius;
                if (d2 >= rSum*rSum || d2 < 1e-12f) continue;
                F d = std::sqrt(d2);
                F ov = rSum - d;
                maxOv = std::max(maxOv, ov);
                F push = STEP_FRAC * ov / d;
                a.pos(0) += dx * push;  a.pos(1) += dy * push;
                b.pos(0) -= dx * push;  b.pos(1) -= dy * push;
                clampX(a); clampX(b);
                // std::cout << a.pos(0) << std::endl;
            }
            wallRelax(discs);
        } while (maxOv > TOL * radiusMean);
    }

    for (auto& p : discs) p.X0 = p.pos;
    return discs;
}



void Simulation::colorParticleRed(int particleID) {
    if (particleID < 0 || particleID >= static_cast<int>(particles2D.size())) {
        std::cerr << "Error: Invalid particle ID: " << particleID << std::endl;
        return;
    }
    // Set the particle's color to red (RGB: 1, 0, 0)
    particles2D[particleID].color = Color(1.0f, 0.0f, 0.0f);
}


void Simulation::detectBoundaryCollision2D(Particle2D &p) const {
    p.BoundaryCollision = 0;
    // Loop over each scenario object and use the virtual function.
    for (const auto &so : scenarioObjects) {
        int collision = so->detectCollision(p);
        if (collision != 0) {
            p.BoundaryCollision = collision;
            // Optionally break out if you only need to detect the first collision.
            break;
        }
    }
}

void Simulation::startScenarioAnimation() {
    if (!scenarioObjects.empty()) {
        animationTimer = 0;
        // Store the current position of the first scenario object as the starting point.
        initialScenarioPosition = scenarioObjects[0]->position;
        animateScenario = true;
    }
}

void Simulation::updateScenarioAnimation(F dt) {
    if (!animateScenario || scenarioObjects.empty())
        return;

    // Increment the timer.
    animationTimer += dt;
    if (animationTimer > animationDuration) {
        animationTimer = animationDuration;
        animateScenario = false; // Stop the animation once complete.
    }

    // Compute progress as a fraction between 0 and 1.
    F t = animationTimer / animationDuration;

    // Start from the initial position.
    Vector3F newPos = initialScenarioPosition;
    
    // Apply linear motion along the x-axis.
    // newPos(0) += t * animationDistance;
    
    // Superimpose a shaking motion along the y-axis.
    // Here, the y-offset is given by 0.5*sin(4PI*t) (t goes from 0 to 1).
    newPos(0) += animationDistance * sin(8.0 * M_PI * t);

    // Update each scenario object (regenerate its vertices, and update auxiliary structures).
    for (auto &so : scenarioObjects) {
         so->position = newPos;
         so->generateVertices();
         buildGridDataStructure();
         updateAuxiliaryStructures();
    }
}

int Simulation::stride() const {          // (= DOF per particle)
    return boolSoftDEM ? DOF + 1   // x , y , epsV
                       : DOF;      // x , y
}

void Simulation::computeCFL() 
{
    F advectiveDisplacement = V0 * timeStep;

std::cout 
    << "Advective CFL check: V0*Δt = " << advectiveDisplacement
    << "  ;  minimum allowed = " << minParticleDiam
    << "  →  " 
    << (advectiveDisplacement < minParticleDiam ? "OK\n" : "TOO LARGE!\n");
}

/* ================================================================== */
/*  BUILD MASS MATRIX (diagonal)                                      */
/* ================================================================== */
void Simulation::buildMassMatrix(SparseMatrixF& M) const
{
    const int S = stride();
    const int n = int(particles2D.size());
    std::vector<Trip> T; T.reserve(S*n);

    for (int i = 0; i < n; ++i) {
        int base = S * i;
        F m = particles2D[i].mass;
        T.emplace_back(base    , base    , m);   // x
        T.emplace_back(base + 1, base + 1, m);   // y
        if (boolSoftDEM)
            T.emplace_back(base + 2, base + 2,
                           0.25 * m * particles2D[i].radius * particles2D[i].radius); // ε
    }
    
    M.resize(S*n,S*n);
    M.setFromTriplets(T.begin(),T.end());
}


Square::Square(F halfLength, F halfWidth, const Vector3F& pos)
    : halfLength(halfLength), halfWidth(halfWidth)
{
    position = pos;
    generateVertices();
}

void Square::generateVertices() {
    vertices.clear();
    // Compute bounds.
    min_x = position(0) - halfLength;
    max_x = position(0) + halfLength;
    min_y = position(1) - halfWidth;
    max_y = position(1) + halfWidth;
    
    // Update the cached bounding box.
    BB.min_x = min_x;
    BB.max_x = max_x;
    BB.min_y = min_y;
    BB.max_y = max_y;
    
    // Compute the vertices.
    vertices.push_back(Vector3F(min_x, min_y, position(2)));
    vertices.push_back(Vector3F(min_x, max_y, position(2)));
    vertices.push_back(Vector3F(max_x, max_y, position(2)));
    vertices.push_back(Vector3F(max_x, min_y, position(2)));
}

int Square::detectCollision(const Particle2D &p) const {
    // Use the cached bounding box (BB) for collision detection.
    if ((p.pos(0) - p.effectiveRadius()) < BB.min_x || (p.pos(0) + p.effectiveRadius()) > BB.max_x ||
        (p.pos(1) - p.effectiveRadius()) < BB.min_y || (p.pos(1) + p.effectiveRadius()) > BB.max_y)
    {
        return 2; // Collision with square boundary.
    }
    return 0;
}

// Implementation for Circle.

Circle::Circle(F radius, I numSegments, const Vector3F& pos)
    : radius(radius), numSegments(numSegments)
{
    position = pos;
    generateVertices();
}

void Circle::generateVertices() {
    vertices.clear();
    // Generate vertices approximating the circle.
    for (I i = 0; i < numSegments; i++) {
        F theta = 2.0f * F(M_PI) * F(i) / F(numSegments);
        F x = radius * std::cos(theta);
        F y = radius * std::sin(theta);
        vertices.push_back(Vector3F(position(0) + x, position(1) + y, position(2)));
    }
    
    // Update the bounding box for the circle.
    BB.min_x = position(0) - radius;
    BB.max_x = position(0) + radius;
    BB.min_y = position(1) - radius;
    BB.max_y = position(1) + radius;
}

int Circle::detectCollision(const Particle2D &p) const {
    // Use the circle's center and radius to detect collision.
    F dx = p.pos(0) - position(0);
    F dy = p.pos(1) - position(1);
    F distSq = dx * dx + dy * dy;
    F limit = radius - p.effectiveRadius();
    if (distSq > (limit * limit))
        return 1; // Collision with circle boundary.
    return 0;
}

Tunnel2D::Tunnel2D(F halfLen, F halfWid, const Vector3F& pos)
    : halfWidth(halfWid), halfLength(halfLen)
{
    position = pos;
    generateVertices();
}

void Tunnel2D::generateVertices()
{
    vertices.clear();

    /* true physical Y–walls */
    BB.min_y = position(1) - halfWidth;
    BB.max_y = position(1) + halfWidth;

    /* we expose an *artificially long* X span so the viewer can draw it.
       It *doesn’t* constrain particles – X is handled by periodic wrapping
       inside Simulation. */
    BB.min_x = position(0) - halfLength;
    BB.max_x = position(0) + halfLength;

    /* rectangle mesh (purely visual) */
    vertices.emplace_back(BB.min_x, BB.min_y, position(2));
    vertices.emplace_back(BB.min_x, BB.max_y, position(2));
    vertices.emplace_back(BB.max_x, BB.max_y, position(2));
    vertices.emplace_back(BB.max_x, BB.min_y, position(2));
}

int Tunnel2D::detectCollision(const Particle2D& p) const
{
    /* only top / bottom walls act as barriers */
    if (p.pos(1) - p.effectiveRadius() < BB.min_y || p.pos(1) + p.effectiveRadius() > BB.max_y)
        return 3;                 // same collision code used by Square
    return 0;
}

Particle2D::Particle2D(F radius, const Simulation& sim)
    : pos(Vector2F::Zero()), vel(Vector2F::Zero()), acc(Vector2F::Zero()), radius(radius)
{
    // For a 2D disc, mass = area * density.
    mass = M_PI * radius * radius * sim.density;
    // Moment of inertia for a uniform disc about its center: I = 1/2 * m * r^2.
    inertia = 0.5 * mass * radius * radius;

    X0      = pos;                     // take spawn position as “rest”
    F_prev.setIdentity();
    D.setZero();
}


Particle3D::Particle3D(F radius, const Simulation& sim)
    : pos(Vector3F::Zero()), vel(Vector3F::Zero()), acc(Vector3F::Zero()), radius(radius)
{
    // For a sphere, mass = volume * density.
    mass = (4.0 / 3.0) * M_PI * std::pow(radius, 3) * sim.density;
    // Moment of inertia for a solid sphere: I = 2/5 * m * r^2.
    inertia = (2.0 / 5.0) * mass * radius * radius;
}

F Particle2D::effectiveRadius() const
{
    return radius * (1.0f + epsV);
}