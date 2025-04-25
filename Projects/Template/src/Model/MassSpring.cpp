#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

#include "Projects/Template/include/Model/MassSpring.h"
#include "CRLHelper/MapleHelper.h"

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <unordered_set>

/* ------------------------------------------------------------------ */
/*  Local helpers / constants                                          */
/* ------------------------------------------------------------------ */
namespace {

using Vec2  = Eigen::Matrix<F,2,1>;
using Vec3  = Eigen::Matrix<F,3,1>;
using Trip  = Eigen::Triplet<F>;

inline F cube(F x)                   { return x * x * x; }
constexpr F oneThird   = static_cast<F>(1.0 / 3.0);
constexpr F half       = static_cast<F>(0.5);
constexpr F two        = static_cast<F>(2.0);
constexpr F tiny       = static_cast<F>(1e-12);

/* uniform wrapper for ImGui scalar inputs -------------------------------- */
void scalarInput(const char* label, F& var, F step, const char* fmt = "%.3f")
{ ImGui::InputDouble(label, &var, step, step*5, fmt); }

} // anonymous namespace

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
        F x0 = p.pos(0) - p.radius;
        F x1 = p.pos(0) + p.radius;
        F y0 = p.pos(1) - p.radius;
        F y1 = p.pos(1) + p.radius;

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
        int cxMinRaw = (int)std::floor((cxWrapped - p.radius - minX) / cellSize);
        int cxMaxRaw = (int)std::floor((cxWrapped + p.radius - minX) / cellSize);
        int cyMin    = std::max(0,  (int)std::floor((p.pos(1)-p.radius - minY)/cellSize));
        int cyMax    = std::min(numCellsY-1,
                                (int)std::floor((p.pos(1)+p.radius - minY)/cellSize));

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
            if (std::sqrt(dx*dx + dy*dy) < p.radius + q.radius)
                p.neighborIndices.push_back(j);
        }
    }
}


void Simulation::updateGlobalPositions() {
    int n = 0;
    if (use3D) {
        n = static_cast<int>(particles3D.size());
        globalPositions.resize(DOF_FULL * n);
        for (int i = 0; i < n; i++) {
            globalPositions.segment<DOF_FULL>(DOF_FULL*i) = particles3D[i].pos;
        }
    }else{
        n = static_cast<int>(particles2D.size());
        globalPositions.resize(DOF * n);
        for (int i = 0; i < n; i++) {
            globalPositions.segment<DOF>(DOF*i) = particles2D[i].pos;
        }
    }
}

// Similarly, when updating particles from the global state:
void Simulation::applyGlobalPositions(const VectorXF &positions) {
    if (use3D) {
        int n = static_cast<int>(particles3D.size());
        for (int i = 0; i < n; i++) {
            particles3D[i].pos = positions.segment<DOF_FULL>(DOF_FULL*i);
        }
    } else {
        int n = (int)particles2D.size();
        for (int i = 0; i < n; ++i) {
            particles2D[i].pos = positions.segment<DOF>(DOF*i);

            if (periodicX)   // keep it canonical
                particles2D[i].pos(0) = wrapX(particles2D[i].pos(0));
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
    // renormalise(); 
}

void Simulation::updateAuxiliaryStructures() {
    // Build the grid and update neighbor lists based on the current particle positions.
    insertParticlesIntoGrid();
    updateNeighborLists();
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
        if (experiment == Experiment::Default)
        {
            const F y0 = scenarioObjects.empty()
                        ? F(0)                       /* fallback datum       */
                        : scenarioObjects[0]->getBoundingBox().min_y;

            const F mgy = (dynamic ? p.mass : F(1)) * gravity(1)
                        * (p.pos(1) - y0);

            value += mgy;
        }

        /* -------------------- boundary objects --------------------------- */
        if (p.BoundaryCollision == 1)                       /* CIRCLE wall  */
        {
            const F px = p.pos(0),  py = p.pos(1),  r = p.radius;

            for (const auto &so : scenarioObjects)
                if (auto *c = dynamic_cast<Circle*>(so.get()))
                {
                    const F dx   = px - c->position(0);
                    const F dy   = py - c->position(1);
                    const F dist = std::sqrt(dx*dx + dy*dy);

                    F depth = dist + r - c->radius;          // δ
                    if (depth <= F(0)) continue;

                    value += oneThird * overlapParam * cube(depth);
                }
        }
        else if (p.BoundaryCollision == 2)                   /* SQUARE wall */
        {
            const F px = p.pos(0),  py = p.pos(1),  r = p.radius;

            for (const auto &so : scenarioObjects)
                if (auto *sq = dynamic_cast<Square*>(so.get()))
                {
                    /* face–wise overlaps (δ ≥ 0) -------------------------- */
                    const F δxL = sq->min_x - (px - r);
                    const F δxR = (px + r) - sq->max_x;
                    const F δyB = sq->min_y - (py - r);
                    const F δyT = (py + r) - sq->max_y;

                    const F δx = (δxR > 0 ? δxR : (δxL > 0 ? δxL : F(0)));
                    const F δy = (δyT > 0 ? δyT : (δyB > 0 ? δyB : F(0)));

                    if (δx==0 && δy==0) continue;

                    /* single face  →  δ = n    (because n=δ)               */
                    /* corner       →  n = √(δx²+δy²)                       */
                    const F n = std::sqrt(δx*δx + δy*δy);

                    value += oneThird * overlapParam * cube(n);
                }
        }
        else if (p.BoundaryCollision == 3)                   /* TUNNEL top/bot */
        {
            const F py = p.pos(1),  r = p.radius;

            for (const auto &so : scenarioObjects)
                if (auto *tu = dynamic_cast<Tunnel2D*>(so.get()))
                {
                    F δ = F(0);
                    if (py - r < tu->BB.min_y)       δ = tu->BB.min_y - (py - r);
                    else if (py + r > tu->BB.max_y)  δ = (py + r) - tu->BB.max_y;

                    if (δ > 0)
                        value += oneThird * overlapParam * cube(δ);
                }
        }
    }

    /* ------------------- inter–particle contacts ------------------------- */
    for (std::size_t i = 0; i < particles2D.size(); ++i)
    {
        const Particle2D &p = particles2D[i];

        for (int j : p.neighborIndices)           /* j > i by construction */
        {
            const Particle2D &q = particles2D[j];

            /* periodic dx so cell lists work                                    */
            const F dx = periodicDx(p.pos(0), p.ix, q.pos(0), q.ix);
            const F dy = p.pos(1) - q.pos(1);
            const F d  = std::sqrt(dx*dx + dy*dy);

            const F δ = (p.radius + q.radius) - d;         /* overlap */
            if (δ <= F(0)) continue;

            value += oneThird * interactionParam * cube(δ);
        }
    }

    /* ------------- pin spring (prevent rigid motion) -------------------- */
    if (periodicX && !particles2D.empty())
    {
        const F dx0 = particles2D[0].pos(0) - pinXref;
        value += F(0.5) * pinK * dx0 * dx0;
    }
}

/* -------------------------------------------------------------------------- */
/*  GRADIENT                                                                  */
/* -------------------------------------------------------------------------- */
void Simulation::compute_gradient(VectorXF &g) const
{
    g.resize(DOF * particles2D.size());
    g.setZero(); 

    for (std::size_t ii = 0; ii < particles2D.size(); ++ii)
    {
        Particle2D p = particles2D[ii];
        detectBoundaryCollision2D(p);

        const int xIdx = DOF*ii;
        const int yIdx = xIdx + 1;

        /* gravity ---------------------------------------------------------- */
        if (experiment == Experiment::Default)
            g[yIdx] += (dynamic ? p.mass : F(1)) * gravity(1);

        /* ----------------------------- CIRCLE wall ----------------------- */
        if (p.BoundaryCollision == 1)
        {
            const F px = p.pos(0), py = p.pos(1), r = p.radius;

            for (const auto &so : scenarioObjects)
                if (auto *c = dynamic_cast<Circle*>(so.get()))
                {
                    const F dx = px - c->position(0);
                    const F dy = py - c->position(1);
                    const F dist2 = dx*dx + dy*dy;
                    if (dist2 == F(0)) continue;

                    const F dist  = std::sqrt(dist2);
                    const F depth = dist + r - c->radius;          /* δ */
                    if (depth <= F(0)) continue;

                    const F kΔ2_over_d = overlapParam * depth * depth / dist;

                    g[xIdx] += kΔ2_over_d * dx;
                    g[yIdx] += kΔ2_over_d * dy;
                }
        }
        /* ----------------------------- SQUARE wall ----------------------- */
        else if (p.BoundaryCollision == 2)
        {
            const F px = p.pos(0), py = p.pos(1), r = p.radius;

            for (const auto &so : scenarioObjects)
                if (auto *sq = dynamic_cast<Square*>(so.get()))
                {
                    const F δxL = sq->min_x - (px - r);
                    const F δxR = (px + r) - sq->max_x;
                    const F δyB = sq->min_y - (py - r);
                    const F δyT = (py + r) - sq->max_y;

                    F δx = F(0), δy = F(0);
                    int sx = 0, sy = 0;

                    if (δxR > 0)      { δx = δxR;  sx = +1; }
                    else if (δxL > 0) { δx = δxL;  sx = -1; }

                    if (δyT > 0)      { δy = δyT;  sy = +1; }
                    else if (δyB > 0) { δy = δyB;  sy = -1; }

                    if (δx==0 && δy==0) continue;

                    if (δx>0 && δy>0)                        /* CORNER */
                    {
                        const F n = std::sqrt(δx*δx + δy*δy);
                        const F coef = overlapParam * n;

                        g[xIdx] += coef * sx * δx;
                        g[yIdx] += coef * sy * δy;
                    }
                    else                                     /* SINGLE */
                    {
                        if (δx > 0) g[xIdx] += sx * overlapParam * δx * δx;
                        if (δy > 0) g[yIdx] += sy * overlapParam * δy * δy;
                    }
                }
        }
        /* ----------------------------- TUNNEL ---------------------------- */
        else if (p.BoundaryCollision == 3)
        {
            const F py = p.pos(1), r = p.radius;

            for (const auto &so : scenarioObjects)
                if (auto *tu = dynamic_cast<Tunnel2D*>(so.get()))
                {
                    F δ = F(0);  int sy = 0;
                    if (py - r < tu->BB.min_y)      { δ = tu->BB.min_y - (py - r); sy = -1; }
                    else if (py + r > tu->BB.max_y) { δ = (py + r) - tu->BB.max_y; sy = +1; }

                    if (δ > 0)
                        g[yIdx] += overlapParam * sy * δ * δ;
                }
        }
    }

    /* --------------------- inter-particle overlaps ----------------------- */
    for (std::size_t i = 0; i < particles2D.size(); ++i)
    {
        const Particle2D &p = particles2D[i];
        const int ix = DOF*i, iy = ix+1;

        for (int j : p.neighborIndices)
        {
            if (j <= static_cast<int>(i)) continue;
            const Particle2D &q = particles2D[j];

            const F dx = periodicDx(p.pos(0), p.ix, q.pos(0), q.ix);
            const F dy = p.pos(1) - q.pos(1);
            const F dist2 = dx*dx + dy*dy;
            if (dist2 == F(0)) continue;

            const F dist = std::sqrt(dist2);
            const F δ    = (p.radius + q.radius) - dist;
            if (δ <= F(0)) continue;

            const F kΔ2_over_d = interactionParam * δ * δ / dist;

            const F gx = kΔ2_over_d * dx;
            const F gy = kΔ2_over_d * dy;

            g[ix] -= gx;  g[iy] -= gy;          /* particle i  */
            g[DOF*j] += gx; g[DOF*j+1] += gy;   /* particle j  */
        }
    }

    /* spring pin in x ------------------------------------------------------ */
    if (periodicX && !particles2D.empty())
        g[0] += pinK * (particles2D[0].pos(0) - pinXref);
}

/* -------------------------------------------------------------------------- */
/*  HESSIAN                                                                   */
/*  (sparse COO additions – same pattern everywhere)                          */
/* -------------------------------------------------------------------------- */
void Simulation::compute_hessian(SparseMatrixF &H) const
{
    const int N = DOF * particles2D.size();
    H.resize(N, N);  H.setZero();

    /* lambda to accumulate entries ---------------------------------------- */
    auto add = [&](int r,int c,F v){ H.coeffRef(r,c) += v; };

    /* --------------------------- boundary walls -------------------------- */
    for (std::size_t ii = 0; ii < particles2D.size(); ++ii)
    {
        Particle2D p = particles2D[ii];
        detectBoundaryCollision2D(p);

        const int xIdx = DOF*ii;
        const int yIdx = xIdx+1;

        /* --------------------------- CIRCLE ----------------------------- */
        if (p.BoundaryCollision == 1)
        {
            const F px = p.pos(0), py = p.pos(1), r = p.radius;
        
            for (const auto &so : scenarioObjects)
                if (auto *c = dynamic_cast<Circle *>(so.get()))
                {
                    const F dx   = px - c->position(0);
                    const F dy   = py - c->position(1);
                    const F d2   = dx*dx + dy*dy;
                    if (d2 <= tiny) continue;                 // centres coincide
        
                    const F  d   = std::sqrt(d2);
                    const F  δ   = d + r - c->radius;         // penetration
                    if (δ <= F(0)) continue;
        
                    const F kΔ   = overlapParam * δ;          // k Δ
                    const F kΔ2  = kΔ * δ;                    // k Δ²
        
                    const F invD2 = F(1) / d2;
                    const F invD3 = invD2 / d;
                    const F kΔ2_over_d = kΔ2 / d;             // common c₁ = kΔ²/d
        
                    /* pre-squared helpers */
                    const F dx2 = dx*dx, dy2 = dy*dy, dxy = dx*dy;
        
                    /* Hessian entries (see Maple) */
                    const F Hxx = 2.0 * dx2 * invD2 * kΔ    //  + 2 kΔ (dx/d)²
                                -       dx2 * invD3 * kΔ2   //  − kΔ² dx² / d³
                                + kΔ2_over_d;               //  + kΔ² / d
        
                    const F Hyy = 2.0 * dy2 * invD2 * kΔ
                                -       dy2 * invD3 * kΔ2
                                + kΔ2_over_d;
        
                    const F Hxy = 2.0 * dxy * invD2 * kΔ
                                -       dxy * invD3 * kΔ2;
        
                    add(xIdx, xIdx, Hxx);  add(yIdx, yIdx, Hyy);
                    add(xIdx, yIdx, Hxy);  add(yIdx, xIdx, Hxy);   // symmetry
                }
        }
        /* --------------------------- SQUARE ----------------------------- */
        else if (p.BoundaryCollision == 2)
        {
            const F px = p.pos(0), py = p.pos(1), r = p.radius;

            for (const auto &so : scenarioObjects)
                if (auto *sq = dynamic_cast<Square*>(so.get()))
                {
                    const F δxL = sq->min_x - (px - r);
                    const F δxR = (px + r) - sq->max_x;
                    const F δyB = sq->min_y - (py - r);
                    const F δyT = (py + r) - sq->max_y;

                    F δx = F(0), δy = F(0);
                    if (δxR > 0)      δx = δxR;
                    else if (δxL > 0) δx = δxL;

                    if (δyT > 0)      δy = δyT;
                    else if (δyB > 0) δy = δyB;

                    if (δx==0 && δy==0) continue;

                    const F k = overlapParam;

                    if (δx>0 && δy>0)                            /* CORNER */
                    {
                        const F n     = std::sqrt(δx*δx + δy*δy);
                        const F inv_n = F(1) / n;
                        const F Hxx = k * ( n + δx*δx*inv_n );
                        const F Hyy = k * ( n + δy*δy*inv_n );
                        const F Hxy = k * ( δx*δy*inv_n );

                        add(xIdx,xIdx,Hxx); add(yIdx,yIdx,Hyy);
                        add(xIdx,yIdx,Hxy); add(yIdx,xIdx,Hxy);
                    }
                    else                                          /* FACE  */
                    {
                        if (δx>0) add(xIdx,xIdx, 2.0 * k * δx);
                        if (δy>0) add(yIdx,yIdx, 2.0 * k * δy);
                    }
                }
        }
        /* --------------------------- TUNNEL ----------------------------- */
        else if (p.BoundaryCollision == 3)
        {
            const F py = p.pos(1), r = p.radius;

            for (const auto &so : scenarioObjects)
                if (auto *tu = dynamic_cast<Tunnel2D*>(so.get()))
                {
                    F δ = F(0);
                    if (py - r < tu->BB.min_y)      δ = tu->BB.min_y - (py - r);
                    else if (py + r > tu->BB.max_y) δ = (py + r) - tu->BB.max_y;

                    if (δ > 0) add(yIdx,yIdx, 2.0 * overlapParam * δ);
                }
        }
    }

    /* ----------------------- inter-particle blocks ----------------------- */
    for (std::size_t i = 0; i < particles2D.size(); ++i)
    {
        const Particle2D &p = particles2D[i];
        const int ix = DOF*i, iy = ix+1;

        for (int j : p.neighborIndices)
        {
            if (j <= static_cast<int>(i)) continue;
            const Particle2D &q = particles2D[j];

            const F dx = periodicDx(p.pos(0), p.ix, q.pos(0), q.ix);
            const F dy = p.pos(1) - q.pos(1);
            const F d2 = dx*dx + dy*dy;
            if (d2 == F(0)) continue;

            const F d  = std::sqrt(d2);
            const F δ  = (p.radius + q.radius) - d;
            if (δ <= F(0)) continue;

            const F kΔ  = interactionParam * δ;
            const F kΔ2 = kΔ * δ;

            const F invD2 = F(1)/d2;
            const F invD3 = invD2 / d;

            const F fac1 = 2.0 * kΔ  * invD2;
            const F fac2 =       kΔ2 * invD3;
            const F T27  =       kΔ2 /  d;

            const F dx2=dx*dx, dy2=dy*dy, dxy=dx*dy;

            const F Hii_xx = fac1*dx2 + fac2*dx2 - T27;
            const F Hii_yy = fac1*dy2 + fac2*dy2 - T27;
            const F Hii_xy = fac1*dxy + fac2*dxy;

            const F Hij_xx = -Hii_xx;
            const F Hij_yy = -Hii_yy;
            const F Hij_xy = -Hii_xy;

            const int jx = DOF*j, jy = jx+1;

            /* i-block */
            add(ix,ix,Hii_xx); add(iy,iy,Hii_yy);
            add(ix,iy,Hii_xy); add(iy,ix,Hii_xy);

            /* j-block */
            add(jx,jx,Hii_xx); add(jy,jy,Hii_yy);
            add(jx,jy,Hii_xy); add(jy,jx,Hii_xy);

            /* i–j off-diagonal */
            add(ix,jx,Hij_xx); add(jx,ix,Hij_xx);
            add(ix,jy,Hij_xy); add(jy,ix,Hij_xy);
            add(iy,jx,Hij_xy); add(jx,iy,Hij_xy);
            add(iy,jy,Hij_yy); add(jy,iy,Hij_yy);
        }
    }

    /* spring pin in x ------------------------------------------------------ */
    if (periodicX && !particles2D.empty())
        add(0,0,pinK);
}

void Simulation::compute_energy_dyn(F &value) {

    compute_energy(value);
    /* make sure previous states exist -------------------------------------- */
    if (globalState_1.size() != globalPositions.size()) {
        globalState_1 = globalPositions;
        globalState_2 = globalPositions;
    }

    /* 3) Δ²-like difference  a = xⁿ⁺¹ − 2xⁿ + xⁿ⁻¹ ----------------------- */
    const VectorXF a = globalPositions
                     - 2 * globalState_1
                     + globalState_2;

    const F inv_h2 = 1.0 / (timeStep * timeStep);

    /* ½/h² · aᵀ M a ------------------------------------------------------- */
    value += 0.5 * inv_h2 * a.dot(M * a);

    if (viscosity) {
        for (size_t i = 0; i < particles2D.size(); i++) {
            // Use the filtered effective neighbor count stored in the particle.
            F effectiveCountFiltered = particles2D[i].prevEffectiveCount;
            VectorXF diff = globalPositions.segment(DOF * i, DOF) - 
                            globalState_1.segment(DOF * i, DOF);
            value += 0.5 * viscosityCoeff * effectiveCountFiltered / (timeStep * timeStep) * diff.squaredNorm();
        }
    }
    if (experiment == Experiment::ShearFlow)
    {
        for (int i = 0; i < (int)particles2D.size(); ++i)
        {
            const int base = DOF * i;

            /* current & previous positions */
            F x1 = globalPositions(base    );
            F y1 = globalPositions(base + 1);
            F x0 = globalState_1  (base    );
            F y0 = globalState_1  (base + 1);

            /* fluid velocity */
            F vfx , dvf_dy , dummy;
            shearFlowProfile(y1, vfx, dvf_dy, dummy);

            /* particle velocity components */
            F dvx = (x1 - x0)/timeStep - vfx;
            F dvy = (y1 - y0)/timeStep;

            value += 0.5 * fluidViscosity *
                     (dvx*dvx + dvy*dvy);
        }
    }
    
}

void Simulation::compute_gradient_dyn(VectorXF &gradient) {
    
    compute_gradient(gradient);
    if (globalState_1.size() != globalPositions.size()) {
        globalState_1 = globalPositions;
        globalState_2 = globalPositions;
    }

    const VectorXF a      = globalPositions
                          - 2 * globalState_1
                          + globalState_2;
    const F        inv_h2 = 1.0 / (timeStep * timeStep);

    /* ∇E_dyn = 1/h² · M a -------------------------------------------------- */
    gradient += inv_h2 * (M * a);
    
    if (viscosity) {
        for (size_t i = 0; i < particles2D.size(); i++) {
            F effectiveCountFiltered = particles2D[i].prevEffectiveCount;
            gradient.segment(DOF * i, DOF) += viscosityCoeff * effectiveCountFiltered / (timeStep * timeStep) *
                                          (globalPositions.segment(DOF * i, DOF) - 
                                           globalState_1.segment(DOF * i, DOF));
        }
    }
    if (experiment == Experiment::ShearFlow)
    {
        const F nu      = fluidViscosity;
        const F inv_dt  = 1.0 / timeStep;
        const F nu_over_dt  = nu * inv_dt;
        const F nu_over_dt2 = nu * inv_dt * inv_dt;

        for (int i = 0; i < (int)particles2D.size(); ++i)
        {
            const int base = DOF * i;

            /* positions */
            F x1 = globalPositions(base    );
            F y1 = globalPositions(base + 1);
            F x0 = globalState_1  (base    );
            F y0 = globalState_1  (base + 1);

            /* fluid profile & derivative */
            F vfx , dvf_dy , dummy;
            shearFlowProfile(y1, vfx, dvf_dy, dummy);

            /* velocity differences */
            F dvx = (x1 - x0)*inv_dt - vfx;   // (v_p − v_f)_x
            F dvy = (y1 - y0)*inv_dt;         // v_p,y

            // -------- grad_x --------
            gradient(base) += nu_over_dt * dvx;

            // -------- grad_y --------
            gradient(base + 1) += nu_over_dt * dvy  // from v_p,y
                                - nu * dvx * dvf_dy; // from v_f(y)
        }
    }
    
    
}

void Simulation::compute_hessian_dyn(SparseMatrixF &hessian) {

    compute_hessian(hessian);
    const F inv_h2 = 1.0 / (timeStep * timeStep);

    /* H += 1/h² · M  (diagonal, so this is cheap) ------------------------- */
    hessian += inv_h2 * M; 
    
    if (viscosity) {
        for (size_t i = 0; i < particles2D.size(); i++) {
            F effectiveCountFiltered = particles2D[i].prevEffectiveCount;
            for (int j = 0; j < DOF; j++) {
                int index = static_cast<int>(DOF * i + j);
                hessian.coeffRef(index, index) += viscosityCoeff * effectiveCountFiltered / (timeStep * timeStep);
            }
        }
    }
    if (experiment == Experiment::ShearFlow)
    {
        const F nu       = fluidViscosity;
        const F inv_dt   = 1.0 / timeStep;
        const F nu_dt2   = nu * inv_dt * inv_dt;
        const F k        = M_PI / L;           // π/L

        for (int i = 0; i < (int)particles2D.size(); ++i)
        {
            const int base = DOF * i;

            /* positions & velocities */
            F x1 = globalPositions(base    );
            F y1 = globalPositions(base + 1);
            F x0 = globalState_1  (base    );

            F sin_k_y , cos_k_y;
            {
                F vfx , dvf_dy , d2vf_dy2;
                shearFlowProfile(y1, vfx, dvf_dy, d2vf_dy2);
                sin_k_y = std::sin(k*y1);
                cos_k_y = std::cos(k*y1);
            }

            /* pre‑compute   dvx   */
            F dvx = (x1 - x0)*inv_dt - V0*sin_k_y;

            /* ===== diagonal blocks ===== */
            hessian.coeffRef(base    , base    ) += nu_dt2;   // d²E/dx²

            /* d²E/dy² :
               nu_dt2  from v_p,y term
             + nu * (V0 k)^2 * cos²(·)   from ∂dvx/∂y
             + nu * dvx * V0 * k^2 * sin(·)   from ∂²v_f/∂y²
            */
            F term1 = nu_dt2;
            F term2 = nu * (V0*V0) * k*k * cos_k_y*cos_k_y;
            F term3 = nu * dvx * V0 * k*k * sin_k_y;
            hessian.coeffRef(base + 1, base + 1) += term1 + term2 + term3;

            /* ===== off‑diagonal  d²E/dxdy  (symmetric) =====
               ∂grad_x/∂y = nu * (-V0 k cos) / dt
            */
            F off = -nu * V0 * k * cos_k_y * inv_dt;
            hessian.coeffRef(base    , base + 1) += off;
            hessian.coeffRef(base + 1, base    ) += off;
        }
    }
    
    
}

inline void Simulation::shearFlowProfile(F y,
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

    // If the domain has not been initialised yet (W == 0)
    // simply return the original coordinate to avoid Inf / NaN.
    if (W == 0 || !periodicX)
        return x;

    return x - std::floor((x - minX) / W) * W;
}

inline F Simulation::periodicDx(F x1,int ix1, F x2,int ix2) const
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
            globalState_1(DOF*i    ) += shift;
            globalState_2(DOF*i    ) += shift;
        }
    }
}

void Simulation::updateEffectiveNeighborCounts() {
    // Loop over each particle.
    for (size_t i = 0; i < particles2D.size(); i++) {
        F rawEffectiveCount = 0;

        // Fetch mod‑position + wrap counter for convenience
        const F xi   = globalPositions(DOF*i    );  // already ∈ [minX,maxX)
        const int ixi = particles2D[i].ix;

        // Compute the raw effective count using a soft kernel w(d) = exp(−d²/σ²)
        for (int j : particles2D[i].neighborIndices) {
            // x‑difference with correct periodic image
            const F xj   = globalPositions(DOF*j    );
            const int ixj = particles2D[j].ix;
            F dx = periodicDx(xi, ixi, xj, ixj);

            // y is non‑periodic
            F dy = globalPositions(DOF*i + 1)
                 - globalPositions(DOF*j + 1);

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

        // Cached mod‐position + wrap counter
        const F  xi   = globalPositions(DOF*i    );  // ∈ [minX,maxX)
        const int ixi = particles2D[i].ix;

        // Soft‐kernel sum over neighbours
        for (int j : particles2D[i].neighborIndices) {
            const F  xj   = globalPositions(DOF*j    );
            const int ixj = particles2D[j].ix;

            F dx = periodicDx(xi, ixi, xj, ixj);
            F dy = globalPositions(DOF*i + 1)
                 - globalPositions(DOF*j + 1);

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


Particle2D::Particle2D(F radius, const Simulation& simParams)
    : pos(Vector2F::Zero()), vel(Vector2F::Zero()), acc(Vector2F::Zero()), radius(radius)
{
    // For a 2D disc, mass = area * density.
    mass = M_PI * radius * radius * simParams.density;
    // Moment of inertia for a uniform disc about its center: I = 1/2 * m * r^2.
    inertia = 0.5 * mass * radius * radius;
}

Particle3D::Particle3D(F radius, const Simulation& sim)
    : pos(Vector3F::Zero()), vel(Vector3F::Zero()), acc(Vector3F::Zero()), radius(radius)
{
    // For a sphere, mass = volume * density.
    mass = (4.0 / 3.0) * M_PI * std::pow(radius, 3) * sim.density;
    // Moment of inertia for a solid sphere: I = 2/5 * m * r^2.
    inertia = (2.0 / 5.0) * mass * radius * radius;
}

std::vector<Particle2D> Simulation::createRandomParticles2D() {
    std::vector<Particle2D> particles;
    particles.reserve(numParticles);

    // ------------------------------------------------------------
    // ❶  Find the area where we may drop particles
    // ------------------------------------------------------------
    F xMin = -0.8 , xMax = 0.8;          // ← fall‑back values
    F yMin = -0.5 , yMax = 0.5;

    if (!scenarioObjects.empty()) {
        if (auto* tun = dynamic_cast<Tunnel2D*>(scenarioObjects[0].get())) {
            const BoundingBox& BB = tun->getBoundingBox();
            xMin = BB.min_x;
            xMax = BB.max_x;
            yMin = BB.min_y;
            yMax = BB.max_y;
        }
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<F> radiusDist(radiusMean, radiusStd);   
    std::uniform_real_distribution<F> distX(xMin, xMax);

    // ------------------------------------------------------------
    // ❷  Rejection sampling with exact periodic distance
    // ------------------------------------------------------------
    for (int i = 0; i < numParticles; ++i) {
        bool ok = false;
        int attempts = 0;
        Particle2D cand(0.0, *this);
        cand.ix = 0;  // start in the base cell

        while (!ok && attempts < maxAttemptsPerParticle) {
            ++attempts;

            // sample radius
            F r = std::max<F>(F(0.01), radiusDist(gen));
            cand = Particle2D(r, *this);
            cand.ix = 0;

            // sample position in the fundamental domain
            std::uniform_real_distribution<F> distY(yMin + r, yMax - r);
            cand.pos(0) = distX(gen);
            cand.pos(1) = distY(gen);

            // enforce pos.x ∈ [xMin, xMax)
            // (distX already does this, but just to be safe)
            cand.pos(0) = std::clamp(cand.pos(0), xMin, std::nextafter(xMax, xMin));

            // overlap test using the new 4‑arg periodicDx
            ok = true;
            for (auto& ex : particles) {
                F dx = periodicDx(
                    cand.pos(0), cand.ix,
                    ex.pos(0),   ex.ix
                );
                F dy = cand.pos(1) - ex.pos(1);
                F dist = std::sqrt(dx*dx + dy*dy);
                if (dist < cand.radius + ex.radius) {
                    ok = false;
                    break;
                }
            }
        }

        if (ok) {
            particles.push_back(std::move(cand));
        } else {
            std::cerr << "createRandomParticles2D: could not place particle "
                      << (i+1) << " after " << maxAttemptsPerParticle << " tries.\n";
        }
    }

    // All new particles have ix == 0 and pos.x in [xMin, xMax).
    return particles;
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
    if ((p.pos(0) - p.radius) < BB.min_x || (p.pos(0) + p.radius) > BB.max_x ||
        (p.pos(1) - p.radius) < BB.min_y || (p.pos(1) + p.radius) > BB.max_y)
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
    F limit = radius - p.radius;
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
    if (p.pos(1) - p.radius < BB.min_y || p.pos(1) + p.radius > BB.max_y)
        return 3;                 // same collision code used by Square
    return 0;
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
    const int n = int(particles2D.size());
    std::vector<Trip> T; T.reserve(2*n);

    for (int i=0;i<n;++i)
    {
        const F m = particles2D[i].mass;
        T.emplace_back(DOF*i    ,DOF*i    ,m);
        T.emplace_back(DOF*i + 1,DOF*i + 1,m);
    }
    M.resize(DOF*n,DOF*n);
    M.setFromTriplets(T.begin(),T.end());
}