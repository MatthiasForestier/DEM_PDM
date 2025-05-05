#include "Projects/Template/include/Model/EnergyFunctions.h"
#include "CRLHelper/MapleHelper.h"
#include <cmath>

inline F cube(F x)                   { return x * x * x; }
inline F sqr(F x)                   { return x * x; }
constexpr F oneThird   = static_cast<F>(1.0 / 3.0);
inline F segmentLength(F Rbar, F delta)                // eq. (55) Mollon
{   return 2.0f * std::sqrt(Rbar * delta); }

namespace EnergyFunctions {

    F gravity(const Particle2D& p, const Simulation& sim, F y0) {
        F m = sim.dynamic ? p.mass : F(1);
        return m * sim.gravity(1) * (p.pos(1) - y0);
    }

    F volumetricStrain(const Particle2D& p, const Simulation& sim) {
        if (!sim.boolSoftDEM) return F(0);
        F r = p.radius;
        return F(0.5) * sim.K * M_PI * sqr(r) * sqr(p.epsV);
    }

    F boundaryEnergy(const Particle2D& p, const Simulation& sim) {
        F E = F(0);
        if (p.BoundaryCollision==1) {
        // circle wall
        for (auto &so : sim.scenarioObjects)
            if (auto *c = dynamic_cast<Circle*>(so.get())) {
                F dx = p.pos(0)-c->position(0);
                F dy = p.pos(1)-c->position(1);
                F dist = std::sqrt(dx*dx+dy*dy);
                F δ = dist + p.effectiveRadius() - c->radius;
                if (δ>0) E += oneThird * sim.overlapParam * cube(δ);
            }
        } else if (p.BoundaryCollision==2) {
        // square wall (same as original)
        for (auto &so : sim.scenarioObjects)
            if (auto *sq = dynamic_cast<Square*>(so.get())) {
                F px = p.pos(0), py = p.pos(1), r = p.effectiveRadius();
                F δxL = sq->min_x - (px - r);
                F δxR = (px + r) - sq->max_x;
                F δyB = sq->min_y - (py - r);
                F δyT = (py + r) - sq->max_y;
                F δx = δxR>0?δxR:(δxL>0?δxL:0);
                F δy = δyT>0?δyT:(δyB>0?δyB:0);
                if (δx||δy) {
                    F n = std::sqrt(δx*δx + δy*δy);
                    E += oneThird * sim.overlapParam * cube(n);
                }
            }
        } else if (p.BoundaryCollision==3) {
        // tunnel top/bot
        for (auto &so : sim.scenarioObjects)
            if (auto *tu = dynamic_cast<Tunnel2D*>(so.get())) {
                F py=p.pos(1), r=p.effectiveRadius(), δ=0;
                if (py-r<tu->BB.min_y) δ=tu->BB.min_y-(py-p.effectiveRadius());
                else if (py+r>tu->BB.max_y) δ=(py+p.effectiveRadius())-tu->BB.max_y;
                if (δ>0) E += oneThird * sim.overlapParam * cube(δ);
            }
        }
        return E;
    }

    F interParticleOld(const Particle2D& p, const Particle2D& q, const Simulation& sim) {
        F dx = sim.periodicDx(p.pos(0),p.ix,q.pos(0),q.ix);
        F dy = p.pos(1)-q.pos(1);
        F d = std::sqrt(dx*dx+dy*dy);
        F δ = (p.effectiveRadius()+q.effectiveRadius()) - d;
        return δ>0 ? oneThird * sim.interactionParam * cube(δ) : F(0);
    }

    F interParticleSoft(const Particle2D& p,
    const Particle2D& q,
    const Simulation& sim)
    {
        // 1) compute the two relative coords
        F dx = sim.periodicDx(p.pos(0),p.ix, q.pos(0),q.ix);
        F dy = p.pos(1) - q.pos(1);

        // 2) pack inputs exactly as Maple expects:
        //    x1, y1, eps1, x2, y2, eps2
        //    we choose (dx,dy) vs (0,0) so that x1-x2 = dx, y1-y2 = dy
        F inputs[6] = {
        dx,
        dy,
        p.epsV,
        F(0),
        F(0),
        q.epsV
        };

        // 3) invoke the Maple‐dumped code for U(δ,ε1,ε2)
        // clang-format off
        F x1    = inputs[0];
        F y1    = inputs[1];
        F epsV1 = inputs[2];
        F x2    = inputs[3];
        F y2    = inputs[4];
        F epsV2 = inputs[5];

        F t1  = x1 * x1;
        F t4  = x2 * x2;
        F t6  = pow(y1 - y2, 0.2e1);
        F t8  = sqrt(-0.2e1 * x1 * x2 + t1 + t4 + t6);
        F t10 = p.radius * (0.1e1 + epsV1);
        F t12 = q.radius * (0.1e1 + epsV2);
        F t13 = -t8 + t10 + t12;
        if (t13 <= 1e-6)                      // full overlap or ε < –1
            return 0;
        F t14 = t13 * t13;
        F t17 = sqrt(0.2e1);
        F t21 = sqrt(t13 * (t10 + t12));

        F unknown = 0.1e1 / (t21 * t17 * sim.a + sim.b * t13) * sim.Young * t14 * t13 / 0.3e1;

        // clang-format on
        return unknown;
    }

    F pinSpring(const Simulation& sim) {
        if (!sim.periodicX || sim.particles2D.empty()) return F(0);
        F dx0 = sim.particles2D[0].pos(0) - sim.pinXref;
        return F(0.5) * sim.pinK * dx0 * dx0;
    }

}

namespace GradientFunctions {

    F gravityGradient(const Particle2D& p, const Simulation& sim) {
    // ∂U/∂y = (dynamic?m:1) * g_y
        F m = sim.dynamic ? p.mass : F(1);
        return m * sim.gravity(1);
    }

    F volumetricStrainGradient(const Particle2D& p, const Simulation& sim) {
        if (!sim.boolSoftDEM) return F(0);
        F eps = p.epsV;
        return sim.K * M_PI * sqr(p.radius) * eps;
    }

    void boundaryGradient(const Particle2D& p,
    const Simulation& sim,
    F& fx, F& fy, F& feps)
    {
        fx = fy = feps = F(0);

        // 1) Circle wall
        if (p.BoundaryCollision == 1) {
            for (auto& so : sim.scenarioObjects) {
                if (auto* c = dynamic_cast<Circle*>(so.get())) {
                    F dx = p.pos(0) - c->position(0);
                    F dy = p.pos(1) - c->position(1);
                    F d2 = dx*dx + dy*dy;
                    if (d2 == F(0)) continue;
                    F d = std::sqrt(d2);
                    F δ = d + p.effectiveRadius() - c->radius;
                    if (δ <= F(0)) continue;

                    F k = sim.overlapParam;
                    // ∂E/∂x = k δ² (dx/d)
                    F coeff = k * δ*δ / d;
                    fx += coeff * dx;
                    fy += coeff * dy;

                    // soft‐DEM coupling to ε
                    if (sim.boolSoftDEM)
                        feps += k * δ*δ * p.radius;
                }
            }
        }
        // 2) Square wall
        else if (p.BoundaryCollision == 2) {
            const F px = p.pos(0), py = p.pos(1), r = p.effectiveRadius();
            for (auto& so : sim.scenarioObjects) {
                if (auto* sq = dynamic_cast<Square*>(so.get())) {
                    F δxL = sq->min_x - (px - r);
                    F δxR = (px + r) - sq->max_x;
                    F δyB = sq->min_y - (py - r);
                    F δyT = (py + r) - sq->max_y;

                    F δx = F(0), δy = F(0);
                    int sx = 0, sy = 0;

                    if (δxR > 0)      { δx = (px + p.effectiveRadius()) - sq->max_x; sx = +1; }
                    else if (δxL > 0) { δx = sq->min_x - (px - p.effectiveRadius()); sx = -1; }

                    if (δyT > 0)      { δy = (py + p.effectiveRadius()) - sq->max_y; sy = +1; }
                    else if (δyB > 0) { δy = sq->min_y - (py - p.effectiveRadius()); sy = -1; }

                    if (δx==0 && δy==0) continue;

                    F k = sim.overlapParam;
                    // corner
                    if (δx > 0 && δy > 0) {
                        F n = std::sqrt(δx*δx + δy*δy);
                        F coef = k * n;
                        fx += coef * sx * δx;
                        fy += coef * sy * δy;
                        if (sim.boolSoftDEM)
                            feps += k * n * (δx * r + δy * r);  
                    }
                    // single face
                    else {
                        if (δx > 0){
                            fx += sx * k * δx * δx;
                            if(sim.boolSoftDEM)
                                feps += k * δx * δx * r;
                        }
                        if (δy > 0){
                            fy += sy * k * δy * δy;
                            if(sim.boolSoftDEM)
                                feps += k * δy * δy * r;
                        }
                    }
                }
            }
        }
        // 3) Tunnel top/bottom
        else if (p.BoundaryCollision == 3) {
            const F py = p.pos(1), r = p.radius;
            for (auto& so : sim.scenarioObjects) {
                if (auto* tu = dynamic_cast<Tunnel2D*>(so.get())) {
                    F δ = F(0); int sy = 0;
                    if (py - r < tu->BB.min_y)      { δ = tu->BB.min_y - (py - p.effectiveRadius()); sy = -1; }
                    else if (py + r > tu->BB.max_y) { δ = (py + p.effectiveRadius()) - tu->BB.max_y; sy = +1; }

                    if (δ > 0) {
                        F k = sim.overlapParam;
                        fy += k * sy * δ * δ;
                        if (sim.boolSoftDEM)
                            feps += k * δ * δ * r;
                    }
                }
            }
        }
    }

    std::pair<F,F> interParticleOldGradient(const Particle2D& p,
                                        const Particle2D& q,
                                        const Simulation& sim) {
        F dx = sim.periodicDx(p.pos(0),p.ix,q.pos(0),q.ix);
        F dy = p.pos(1)-q.pos(1);
        F d = std::sqrt(dx*dx+dy*dy);
        F δ = (p.effectiveRadius()+q.effectiveRadius()) - d;
        if (δ<=0) return {F(0),F(0)};
        F Fn = sim.interactionParam * δ*δ;
        return { -Fn*dx/d, -Fn*dy/d };
        }

    std::tuple<F,F,F,F> interParticleSoftGradient(
    const Particle2D& p,
    const Particle2D& q,
    const Simulation& sim)
        {
        // 1) relative coords
        F dx = sim.periodicDx(p.pos(0),p.ix, q.pos(0),q.ix);
        F dy = p.pos(1) - q.pos(1);

        // 2) pack exactly for Maple
        F inputs[6] = {
        dx,
        dy,
        p.epsV,
        F(0),
        F(0),
        q.epsV
        };

        // 3) call Maple‐dumped gradient (6×1)
        F gradU[6];
        {
        // clang-format off
        F x1    = inputs[0];
        F y1    = inputs[1];
        F epsV1 = inputs[2];
        F x2    = inputs[3];
        F y2    = inputs[4];
        F epsV2 = inputs[5];

        F t1 = x1 * x1;
        F t4 = x2 * x2;
        F t5 = y1 - y2;
        F t6 = t5 * t5;
        F t8 = sqrt(-0.2e1 * x1 * x2 + t1 + t4 + t6);
        if (t8 < 1e-14)                     // same particle / identical pos
            return {0,0,0,0};
        F t10 = p.radius * (0.1e1 + epsV1);
        F t12 = q.radius * (0.1e1 + epsV2);
        F t13 = -t8 + t10 + t12;
        if (t13 <= 0)                     // already separated – no Soft‑DEM
            return {0,0,0,0};
        F t14 = t13 * t13;
        if (t14 <= 1e-12)
            return {0,0,0,0};
        F t15 = sim.Young * t14;
        F t16 = sqrt(0.2e1);
        F t17 = t16 * sim.a;
        F t18 = t10 + t12;
        F t20 = sqrt(t13 * t18);
        F t23 = sim.b * t13 + t20 * t17;
        if (std::abs(t23) < 1e-10)        // ill‑conditioned
            return {0,0,0,0};
        F t24 = 0.1e1 / t23 / 0.3e1;
        F t25 = 0.1e1 / t8;
        F t26 = t25 * t24;
        F t27 = x1 - x2;
        F t32 = sim.Young * t14 * t13;
        F t34 = pow(t23, -0.2e1) / 0.9e1;
        F t35 = 0.1e1 / t20;
        F t36 = t35 * t17;
        F t37 = t25 * t18;

        F unknown[6];

        unknown[0] = -0.3e1 * t27 * t26 * t15
                    - (
                        -0.3e1 / 0.2e1 * t27 * t37 * t36
                        - 0.3e1 * sim.b * t27 * t25
                    ) * t34 * t32;
        unknown[1] = -0.3e1 * t5  * t26 * t15
                    - (
                        -0.3e1 / 0.2e1 * t5  * t37 * t36
                        - 0.3e1 * sim.b * t5  * t25
                    ) * t34 * t32;
        unknown[2] =  0.3e1 * p.radius * t24 * t15
                    - (
                        0.3e1 / 0.2e1 * (t13 * p.radius + p.radius * t18) * t35 * t17
                        + 0.3e1 * p.radius * sim.b
                    ) * t34 * t32;
        unknown[3] =  0.3e1 * t27 * t26 * t15
                    - (
                        0.3e1 / 0.2e1 * t27 * t37 * t36
                        + 0.3e1 * sim.b * t27 * t25
                    ) * t34 * t32;
        unknown[4] =  0.3e1 * t5  * t26 * t15
                    - (
                        0.3e1 / 0.2e1 * t5  * t37 * t36
                        + 0.3e1 * sim.b * t5  * t25
                    ) * t34 * t32;
        unknown[5] =  0.3e1 * q.radius * t24 * t15
                    - (
                        0.3e1 / 0.2e1 * (t13 * q.radius + q.radius * t18) * t35 * t17
                        + 0.3e1 * q.radius * sim.b
                    ) * t34 * t32;


        VectorXF gradU(6);
        processMapleOutput(reinterpret_cast<F*>(unknown),
                        /*→ vector overload needs MatrixXF →*/ gradU,
                        /*rows*/6, /*cols*/1);

        // now extract
        F dUx = gradU(0);
        F dUy = gradU(1);
        F dUde1 = gradU(2);
        F dUde2 = gradU(5);

        // **do not** negate these here
        return { dUx, dUy, dUde1, dUde2 };
    }

    // 4) gradU now holds [∂U/∂x1, ∂U/∂y1, ∂U/∂ε1, ∂U/∂x2, ∂U/∂y2, ∂U/∂ε2]
    //    We want forces = -dU/dx, and volumetric couplings = +dU/dε
    F fx    = -gradU[0];
    F fy    = -gradU[1];
    F deps1 =  gradU[2];
    F deps2 =  gradU[5];

    return {fx, fy, deps1, deps2};
    }

    F pinSpringGradient(const Simulation& sim) {
        if (!sim.periodicX || sim.particles2D.empty()) return F(0);
            return sim.pinK * (sim.particles2D[0].pos(0) - sim.pinXref);
    }

}

namespace HessianFunctions {

    void boundaryHessian(
        const Particle2D& p,
        const Simulation& sim,
        int               i,
        AddFunc           add)
    {
        const int S    = sim.stride();
        const int xIdx = S*i+0;
        const int yIdx = S*i+1;
        const int eIdx = S*i+2;

        // 0) volumetric‐strain diag (Soft-DEM only)
        if (sim.boolSoftDEM) {
            F eps  = p.epsV;
            F diag = sim.K * M_PI * sqr(p.radius);
            add(eIdx,eIdx, diag);
        }

        // 1) Circle wall
        if (p.BoundaryCollision == 1) {
            for (auto& so : sim.scenarioObjects) if (auto* c = dynamic_cast<Circle*>(so.get())) {
                F dx = p.pos(0)-c->position(0);
                F dy = p.pos(1)-c->position(1);
                F d2 = dx*dx+dy*dy; if (d2==F(0)) continue;
                F d    = std::sqrt(d2);
                F δ    = d + p.effectiveRadius() - c->radius;
                if (δ <= F(0)) continue;

                F k    = sim.overlapParam;
                F kΔ  = k * δ;
                F kΔ2 = kΔ * δ;
                F invD2 = F(1)/d2;
                F invD3 = invD2/d;

                // translational block
                F dx2 = dx*dx, dy2 = dy*dy, dxy=dx*dy;
                F Hxx = 2*dx2*invD2*kΔ - dx2*invD3*kΔ2 + kΔ2/d;
                F Hyy = 2*dy2*invD2*kΔ - dy2*invD3*kΔ2 + kΔ2/d;
                F Hxy = 2*dxy*invD2*kΔ - dxy*invD3*kΔ2;
                add(xIdx,xIdx,Hxx);
                add(yIdx,yIdx,Hyy);
                add(xIdx,yIdx,Hxy);
                add(yIdx,xIdx,Hxy);

                // ε couplings (Soft-DEM)
                if (sim.boolSoftDEM) {
                    F Hxe = 2 * k * p.radius * dx / d * δ;
                    F Hye = 2 * k * p.radius * dy / d * δ;
                    F Hee = 2 * k * p.radius*p.radius * δ;
                    add(xIdx,eIdx,Hxe); add(eIdx,xIdx,Hxe);
                    add(yIdx,eIdx,Hye); add(eIdx,yIdx,Hye);
                    add(eIdx,eIdx,Hee);
                }
            }
        }
        // 2) Square wall
        else if (p.BoundaryCollision == 2) {
            F px = p.pos(0), py = p.pos(1), r = p.radius;
            for (auto& so : sim.scenarioObjects) if (auto* sq = dynamic_cast<Square*>(so.get())) {
                F δxL = sq->min_x - (px - r);
                F δxR = (px + r) - sq->max_x;
                F δyB = sq->min_y - (py - r);
                F δyT = (py + r) - sq->max_y;

                F δx = F(0), δy = F(0);
                int sx = 0, sy = 0;

                if (δxR > 0)      { δx = (px + p.effectiveRadius()) - sq->max_x; sx = +1; }
                else if (δxL > 0) { δx = sq->min_x - (px - p.effectiveRadius()); sx = -1; }

                if (δyT > 0)      { δy = (py + p.effectiveRadius()) - sq->max_y; sy = +1; }
                else if (δyB > 0) { δy = sq->min_y - (py - p.effectiveRadius()); sy = -1; }
                if (δx==0 && δy==0) continue;

                F k = sim.overlapParam;
                // corner
                if (!sim.boolSoftDEM){
                    if (δx>0 && δy>0) {
                        F n    = std::sqrt(δx*δx + δy*δy);
                        F invn = F(1)/n;
                        F Hxx  = k*(n + δx*δx*invn);
                        F Hyy  = k*(n + δy*δy*invn);
                        F Hxy  = k*(δx*δy*invn*sx*sy);
                        add(xIdx,xIdx,Hxx); add(yIdx,yIdx,Hyy);
                        add(xIdx,yIdx,Hxy); add(yIdx,xIdx,Hxy);
                    }
                    // face
                    else {
                        if (δx>0) add(xIdx,xIdx,2*k*δx);
                        if (δy>0) add(yIdx,yIdx,2*k*δy);
                    }
                }else{
                    if (δx>0 && δy>0) {
                        F n    = std::sqrt(δx*δx + δy*δy);
                        F invn = F(1)/n;
                        F Hxx  = k*(n + δx*δx*invn);
                        F Hyy  = k*(n + δy*δy*invn);
                        F Hxy  = k*(δx*δy*invn*sx*sy);
                        F Hxe = k * (r * (n*sx + (δx + δy)*δx*invn));
                        F Hye = k * (r * (n*sy + (δx + δy)*δy*invn));
                        F Hee = k * sqr(r) * (2 * n + sqr(δx + δy)*invn);
                        add(xIdx,eIdx,Hxe); add(eIdx,xIdx,Hxe);
                        add(yIdx,eIdx,Hye); add(eIdx,yIdx,Hye);
                        add(eIdx,eIdx,Hee);
                        add(xIdx,xIdx,Hxx); add(yIdx,yIdx,Hyy);
                        add(xIdx,yIdx,Hxy); add(yIdx,xIdx,Hxy);
                    }
                    // face
                    else {
                        if (δx>0){
                            add(xIdx,xIdx,2*k*δx);
                            add(xIdx,eIdx,2*k*δx*sx*r);
                            add(eIdx,xIdx,2*k*δx*sx*r);
                            add(eIdx,eIdx,2*k*δx*sqr(r));
                        }
                        if (δy>0){ 
                            add(yIdx,yIdx,2*k*δy);
                            add(yIdx,eIdx,2*k*δy*sy*r);
                            add(eIdx,yIdx,2*k*δy*sy*r);
                            add(eIdx,eIdx,2*k*δy*sqr(r));
                        }
                    }
                }
            }
        }
        // 3) Tunnel wall
        else if (p.BoundaryCollision == 3) {
            F py = p.pos(1), r = p.radius;
            for (auto& so : sim.scenarioObjects) if (auto* tu = dynamic_cast<Tunnel2D*>(so.get())) {
            F δ=0; int sy=0;
            if      (py - r < tu->BB.min_y) { δ = tu->BB.min_y - (py - p.effectiveRadius()); sy = -1; }
            else if (py + r > tu->BB.max_y) { δ = (py + p.effectiveRadius()) - tu->BB.max_y; sy = +1; }
            if (δ>0) add(yIdx,yIdx, 2*sim.overlapParam*δ);
            if (sim.boolSoftDEM) {
                F Hye = 2 * sy * sim.overlapParam * r * δ;
                F Hee = 2 * sim.overlapParam * r * r * δ;
                add(yIdx,eIdx,Hye); add(eIdx,yIdx,Hye);
                add(eIdx,eIdx,Hee);
                }
            }
        }
    }

    void interParticleOldHessian(
        const Particle2D& p,
        const Particle2D& q,
        const Simulation& sim,
        int               i,
        int               j,
        AddFunc           add)
    {
        const int S    = sim.stride();
        const int ix   = S*i+0, iy = S*i+1;
        const int jx   = S*j+0, jy = S*j+1;

        F dx = sim.periodicDx(p.pos(0),p.ix,q.pos(0),q.ix);
        F dy = p.pos(1)-q.pos(1);
        F d2 = dx*dx + dy*dy;
        if (d2 == F(0)) return;
        F d  = std::sqrt(d2);
        F δ  = (p.effectiveRadius()+q.effectiveRadius()) - d;
        if (δ <= 0) return;

        F k    = sim.interactionParam;
        F kΔ   = k * δ;
        F kΔ2  = kΔ * δ;
        F invD2 = F(1)/d2;
        F invD3 = invD2/d;
        F fac1 = 2*kΔ*invD2;
        F fac2 = kΔ2*invD3;
        F T27  = kΔ2 / d;

        F dx2=dx*dx, dy2=dy*dy, dxy=dx*dy;

        // diagonal blocks
        F Hxx = fac1*dx2 + fac2*dx2 - T27;
        F Hyy = fac1*dy2 + fac2*dy2 - T27;
        F Hxy = fac1*dxy + fac2*dxy;
        // off‐diagonal blocks
        F Jxx = -Hxx, Jyy = -Hyy, Jxy = -Hxy;

        // i‐block
        add(ix, ix, Hxx);
        add(iy, iy, Hyy);
        add(ix, iy, Hxy);
        add(iy, ix, Hxy);
        // j‐block
        add(jx, jx, Hxx);
        add(jy, jy, Hyy);
        add(jx, jy, Hxy);
        add(jy, jx, Hxy);
        // cross‐couplings
        add(ix, jx, Jxx); add(jx, ix, Jxx);
        add(ix, jy, Jxy); add(jy, ix, Jxy);
        add(iy, jx, Jxy); add(jx, iy, Jxy);
        add(iy, jy, Jyy); add(jy, iy, Jyy);
    }

    void interParticleSoftHessian(
        const Particle2D& p,
        const Particle2D& q,
        const Simulation& sim,
        int               i,
        int               j,
        AddFunc           add)
    {
        const int S    = sim.stride();
        const int ix   = S*i+0, iy = S*i+1, ie = S*i+2;
        const int jx   = S*j+0, jy = S*j+1, je = S*j+2;

        F dx = sim.periodicDx(p.pos(0), p.ix, q.pos(0), q.ix);
        F dy = p.pos(1) - q.pos(1);
        F d2 = dx*dx+dy*dy; if (d2==F(0)) return;
        F d  = std::sqrt(d2);
        F δ  = (p.effectiveRadius()+q.effectiveRadius()) - d;
        if (δ <= F(0)) return;


        F inputs[6] = {
        // pick x₁ - x₂ = dx  and  y₁ - y₂ = dy
        // we can choose e.g. x₁ = dx/2, x₂ = -dx/2, likewise for y
        dx * F(0.5),
        dy * F(0.5),
        p.epsV,
        -dx * F(0.5),
        -dy * F(0.5),
        q.epsV
        };

        MatrixXF Hloc;
        Hloc.resize(6,6);

        // --- 2) local 6×6 Soft-DEM Hessian block ---

        // clang-format off
        F x1 = inputs[0];
        F y1 = inputs[1];
        F epsV1 = inputs[2];
        F x2 = inputs[3];
        F y2 = inputs[4];
        F epsV2 = inputs[5];

        F t1 = x1 * x1;
        F t4 = x2 * x2;
        F t5 = y1 - y2;
        F t6 = t5 * t5;
        F t7 = -0.2e1 * x1 * x2 + t1 + t4 + t6;
        F t8 = sqrt(t7);
        if (t8 < 1e-14)                     // same particle / identical pos
            return;
        F t10 = p.radius * (0.1e1 + epsV1);
        F t12 = q.radius * (0.1e1 + epsV2);
        F t13 = -t8 + t10 + t12;
        if (t13 <= 0)                     // already separated – no Soft‑DEM
            return;
        F t14 = sim.Young * t13;
        F t15 = sqrt(0.2e1);
        F t16 = t15 * sim.a;
        F t17 = t10 + t12;
        F t18 = t13 * t17;
        F t19 = sqrt(t18);
        F t22 = sim.b * t13 + t19 * t16;
        if (std::abs(t22) < 1e-10)        // ill‑conditioned
            return;
        F t23 = 0.1e1 / t22 / 0.3e1;
        F t24 = 0.1e1 / t7;
        F t25 = t24 * t23;
        F t26 = x1 - x2;
        F t27 = 0.4e1 * t26 * t26;
        F t30 = 0.3e1 / 0.2e1 * t27 * t25 * t14;
        F t31 = t13 * t13;
        F t32 = sim.Young * t31;
        F t33 = 0.9e1 * t22 * t22;
        F t34 = 0.1e1 / t33;
        F t35 = t34 * t32;
        F t36 = 0.1e1 / t8;
        F t37 = 0.2e1 * t26 * t36;
        F t38 = 0.1e1 / t19;
        F t39 = t38 * t16;
        F t40 = t36 * t17;
        F t46 = -0.3e1 / 0.2e1 * t26 * t40 * t39 - 0.3e1 / 0.2e1 * sim.b * t37;
        F t51 = 0.1e1 / t8 / t7;
        F t52 = t51 * t23;
        F t55 = 0.3e1 / 0.4e1 * t27 * t52 * t32;
        F t58 = 0.3e1 * t36 * t23 * t32;
        F t60 = sim.Young * t31 * t13;
        F t63 = 0.1e1 / t33 / t22 / 0.3e1;
        F t64 = t46 * t46;
        F t69 = 0.1e1 / t19 / t18;
        F t70 = t69 * t16;
        F t71 = t17 * t17;
        F t72 = t24 * t71;
        F t76 = t51 * t17;
        F t83 = 0.3e1 / 0.2e1 * t36 * t17 * t38 * t16;
        F t88 = 0.3e1 * sim.b * t36;
        F t91 = (-0.3e1 / 0.16e2 * t27 * t72 * t70 + 0.3e1 / 0.8e1 * t27 * t76 * t39 - t83 + 0.3e1 / 0.4e1 * sim.b * t27 * t51 - t88) * t34 * t60;
        F t93 = t23 * t14;
        F t94 = 0.2e1 * t26 * t24;
        F t101 = 0.2e1 * t5 * t36;
        F t104 = -0.3e1 / 0.2e1 * t5 * t40 * t39 - 0.3e1 / 0.2e1 * sim.b * t101;
        F t108 = t23 * t32;
        F t109 = 0.2e1 * t26 * t51;
        F t113 = t36 * t46;
        F t117 = t46 * t63;
        F t121 = 0.4e1 * t5 * t26;
        F t124 = 0.3e1 / 0.16e2 * t121 * t72 * t70;
        F t127 = 0.3e1 / 0.8e1 * t121 * t76 * t39;
        F t134 = 0.3e1 * t5 * t94 * t93 + 0.3e1 / 0.2e1 * t104 * t37 * t35 + 0.3e1 / 0.2e1 * t5 * t109 * t108 + 0.3e1 * t5 * t113 * t35 + 0.2e1 * t104 * t117 * t60 - (-t124 + t127 + 0.3e1 / 0.2e1 * t5 * sim.b * t109) * t34 * t60;
        F t135 = p.radius * t37;
        F t140 = t13 * p.radius + p.radius * t17;
        F t146 = 0.3e1 / 0.2e1 * t140 * t38 * t16 + 0.3e1 * p.radius * sim.b;
        F t150 = t46 * t34;
        F t166 = -0.3e1 * t135 * t93 + 0.3e1 / 0.2e1 * t146 * t37 * t35 - 0.3e1 * p.radius * t150 * t32 + 0.2e1 * t146 * t117 * t60 - (0.3e1 / 0.4e1 * t140 * t26 * t40 * t70 - 0.3e1 / 0.4e1 * t135 * t39) * t34 * t60;
        F t173 = -0.2e1 * t26 * t36;
        F t176 = 0.3e1 / 0.2e1 * t26 * t40 * t39 - 0.3e1 / 0.2e1 * sim.b * t173;
        F t189 = -0.4e1 * t26 * t26;
        F t196 = -0.2e1 * t26 * sim.b;
        F t202 = -0.3e1 * t26 * t94 * t93 + 0.3e1 / 0.2e1 * t176 * t37 * t35 - 0.3e1 / 0.2e1 * t26 * t109 * t108 + t58 - 0.3e1 * t26 * t113 * t35 + 0.2e1 * t176 * t117 * t60 - (-0.3e1 / 0.16e2 * t189 * t72 * t70 + 0.3e1 / 0.8e1 * t189 * t76 * t39 + t83 + 0.3e1 / 0.4e1 * t196 * t109 + t88) * t34 * t60;
        F t209 = -0.2e1 * t5 * t36;
        F t212 = 0.3e1 / 0.2e1 * t5 * t40 * t39 - 0.3e1 / 0.2e1 * sim.b * t209;
        F t225 = -0.4e1 * t5 * t26;
        F t228 = 0.3e1 / 0.16e2 * t225 * t72 * t70;
        F t231 = 0.3e1 / 0.8e1 * t225 * t76 * t39;
        F t232 = -0.2e1 * t5 * sim.b;
        F t238 = -0.3e1 * t5 * t94 * t93 + 0.3e1 / 0.2e1 * t212 * t37 * t35 - 0.3e1 / 0.2e1 * t5 * t109 * t108 - 0.3e1 * t5 * t113 * t35 + 0.2e1 * t212 * t117 * t60 - (-t228 + t231 + 0.3e1 / 0.4e1 * t232 * t109) * t34 * t60;
        F t239 = q.radius * t37;
        F t244 = t13 * q.radius + q.radius * t17;
        F t250 = 0.3e1 / 0.2e1 * t244 * t38 * t16 + 0.3e1 * q.radius * sim.b;
        F t269 = -0.3e1 * t239 * t93 + 0.3e1 / 0.2e1 * t250 * t37 * t35 - 0.3e1 * q.radius * t150 * t32 + 0.2e1 * t250 * t117 * t60 - (0.3e1 / 0.4e1 * t244 * t26 * t40 * t70 - 0.3e1 / 0.4e1 * t239 * t39) * t34 * t60;
        F t270 = 0.4e1 * t5 * t5;
        F t273 = 0.3e1 / 0.2e1 * t270 * t25 * t14;
        F t279 = 0.3e1 / 0.4e1 * t270 * t52 * t32;
        F t280 = t104 * t104;
        F t295 = (-0.3e1 / 0.16e2 * t270 * t72 * t70 + 0.3e1 / 0.8e1 * t270 * t76 * t39 - t83 + 0.3e1 / 0.4e1 * sim.b * t270 * t51 - t88) * t34 * t60;
        F t297 = p.radius * t101;
        F t303 = t104 * t34;
        F t307 = t104 * t63;
        F t320 = -0.3e1 * t297 * t93 + 0.3e1 / 0.2e1 * t146 * t101 * t35 - 0.3e1 * p.radius * t303 * t32 + 0.2e1 * t146 * t307 * t60 - (0.3e1 / 0.4e1 * t140 * t5 * t40 * t70 - 0.3e1 / 0.4e1 * t297 * t39) * t34 * t60;
        F t321 = 0.2e1 * t5 * t24;
        F t328 = 0.2e1 * t5 * t51;
        F t332 = t36 * t104;
        F t344 = -0.3e1 * t26 * t321 * t93 + 0.3e1 / 0.2e1 * t176 * t101 * t35 - 0.3e1 / 0.2e1 * t26 * t328 * t108 - 0.3e1 * t26 * t332 * t35 + 0.2e1 * t176 * t307 * t60 - (-t228 + t231 + 0.3e1 / 0.4e1 * t196 * t328) * t34 * t60;
        F t360 = -0.4e1 * t5 * t5;
        F t372 = -0.3e1 * t5 * t321 * t93 + 0.3e1 / 0.2e1 * t212 * t101 * t35 - 0.3e1 / 0.2e1 * t5 * t328 * t108 + t58 - 0.3e1 * t5 * t332 * t35 + 0.2e1 * t212 * t307 * t60 - (-0.3e1 / 0.16e2 * t360 * t72 * t70 + 0.3e1 / 0.8e1 * t360 * t76 * t39 + t83 + 0.3e1 / 0.4e1 * t232 * t328 + t88) * t34 * t60;
        F t373 = q.radius * t101;
        F t394 = -0.3e1 * t373 * t93 + 0.3e1 / 0.2e1 * t250 * t101 * t35 - 0.3e1 * q.radius * t303 * t32 + 0.2e1 * t250 * t307 * t60 - (0.3e1 / 0.4e1 * t244 * t5 * t40 * t70 - 0.3e1 / 0.4e1 * t373 * t39) * t34 * t60;
        F t395 = p.radius * p.radius;
        F t399 = p.radius * t34;
        F t403 = t146 * t146;
        F t407 = t140 * t140;
        F t418 = t36 * p.radius;
        F t419 = -0.2e1 * t26 * t418;
        F t425 = t36 * t146;
        F t429 = t146 * t63;
        F t433 = t17 * t140;
        F t442 = -0.3e1 * t419 * t93 - 0.3e1 * t176 * t399 * t32 - 0.3e1 * t26 * t425 * t35 + 0.2e1 * t176 * t429 * t60 - (0.3e1 / 0.8e1 * t173 * t433 * t70 - 0.3e1 / 0.4e1 * t419 * t39) * t34 * t60;
        F t443 = -0.2e1 * t5 * t418;
        F t463 = -0.3e1 * t443 * t93 - 0.3e1 * t212 * t399 * t32 - 0.3e1 * t5 * t425 * t35 + 0.2e1 * t212 * t429 * t60 - (0.3e1 / 0.8e1 * t209 * t433 * t70 - 0.3e1 / 0.4e1 * t443 * t39) * t34 * t60;
        F t489 = 0.6e1 * q.radius * p.radius * t23 * t14 - 0.3e1 * t250 * t399 * t32 - 0.3e1 * q.radius * t146 * t34 * t32 + 0.2e1 * t250 * t429 * t60 - (-0.3e1 / 0.4e1 * t244 * t140 * t69 * t16 + 0.3e1 * q.radius * p.radius * t38 * t16) * t34 * t60;
        F t493 = t176 * t176;
        F t505 = -0.2e1 * t26 * t51;
        F t513 = t176 * t63;
        F t522 = 0.6e1 * t5 * t26 * t24 * t93 + 0.3e1 / 0.2e1 * t212 * t173 * t35 - 0.3e1 / 0.2e1 * t5 * t505 * t108 - 0.3e1 * t5 * t36 * t176 * t35 + 0.2e1 * t212 * t513 * t60 - (-t124 + t127 + 0.3e1 / 0.4e1 * t232 * t505) * t34 * t60;
        F t523 = q.radius * t173;
        F t545 = -0.3e1 * t523 * t93 + 0.3e1 / 0.2e1 * t250 * t173 * t35 - 0.3e1 * q.radius * t176 * t34 * t32 + 0.2e1 * t250 * t513 * t60 - (-0.3e1 / 0.4e1 * t244 * t26 * t40 * t70 - 0.3e1 / 0.4e1 * t523 * t39) * t34 * t60;
        F t549 = t212 * t212;
        F t554 = q.radius * t209;
        F t577 = -0.3e1 * t554 * t93 + 0.3e1 / 0.2e1 * t250 * t209 * t35 - 0.3e1 * q.radius * t212 * t34 * t32 + 0.2e1 * t250 * t212 * t63 * t60 - (-0.3e1 / 0.4e1 * t244 * t5 * t40 * t70 - 0.3e1 / 0.4e1 * t554 * t39) * t34 * t60;
        F t578 = q.radius * q.radius;
        F t586 = t250 * t250;
        F t590 = t244 * t244;
    
        F unknown[6][6];
    
        unknown[0][0] = 0.3e1 * t46 * t37 * t35 + 0.2e1 * t64 * t63 * t60 + t30 + t55 - t58 - t91;
        unknown[0][1] = t134;
        unknown[0][2] = t166;
        unknown[0][3] = t202;
        unknown[0][4] = t238;
        unknown[0][5] = t269;
        unknown[1][0] = t134;
        unknown[1][1] = 0.3e1 * t104 * t101 * t35 + 0.2e1 * t280 * t63 * t60 + t273 + t279 - t295 - t58;
        unknown[1][2] = t320;
        unknown[1][3] = t344;
        unknown[1][4] = t372;
        unknown[1][5] = t394;
        unknown[2][0] = t166;
        unknown[2][1] = t320;
        unknown[2][2] = 0.6e1 * t395 * t23 * t14 - 0.6e1 * t146 * t399 * t32 + 0.2e1 * t403 * t63 * t60 - (-0.3e1 / 0.4e1 * t407 * t69 * t16 + 0.3e1 * t395 * t38 * t16) * t34 * t60;
        unknown[2][3] = t442;
        unknown[2][4] = t463;
        unknown[2][5] = t489;
        unknown[3][0] = t202;
        unknown[3][1] = t344;
        unknown[3][2] = t442;
        unknown[3][3] = 0.3e1 * t176 * t173 * t35 + 0.2e1 * t493 * t63 * t60 + t30 + t55 - t58 - t91;
        unknown[3][4] = t522;
        unknown[3][5] = t545;
        unknown[4][0] = t238;
        unknown[4][1] = t372;
        unknown[4][2] = t463;
        unknown[4][3] = t522;
        unknown[4][4] = 0.3e1 * t212 * t209 * t35 + 0.2e1 * t549 * t63 * t60 + t273 + t279 - t295 - t58;
        unknown[4][5] = t577;
        unknown[5][0] = t269;
        unknown[5][1] = t394;
        unknown[5][2] = t489;
        unknown[5][3] = t545;
        unknown[5][4] = t577;
        unknown[5][5] = 0.6e1 * t578 * t23 * t14 - 0.6e1 * t250 * q.radius * t34 * t32 + 0.2e1 * t586 * t63 * t60 - (-0.3e1 / 0.4e1 * t590 * t69 * t16 + 0.3e1 * t578 * t38 * t16) * t34 * t60;
    


        processMapleOutput(reinterpret_cast<F*>(unknown), Hloc, 6, 6);
        // clang-format on

        // --- 3) scatter the 6×6 block into the global Hessian ---
        const int idx[6] = { ix, iy, ie, jx, jy, je };
        for(int a=0; a<6; ++a)
            for(int b=0; b<6; ++b)
                add(idx[a], idx[b], Hloc(a,b));

     // clang-format off
    }


    void pinSpringHessian(const Simulation& sim, AddFunc add) {
        if (sim.periodicX && !sim.particles2D.empty()) {
        add(0,0, sim.pinK);
        }
    }

}

namespace EnergyDynFunctions {

    F deformationViscous(const Particle2D& p, const Simulation& sim) {
        
        // F x1 = inputs[0];
        // F y1 = inputs[1];
        // F epsV1 = inputs[2];
        // F x2 = inputs[3];
        // F y2 = inputs[4];

        // clang-format off

        // F t2 = pow(0.1e1 + epsV1, 0.2e1);
        // F t4 = r0 * r0;
        // F t6 = Fprev11 * Fprev11;
        // F t8 = Fprev12 * Fprev12;
        // F t9 = Fprev12 * Fprev21;
        // F t11 = Fprev21 * Fprev21;
        // F t12 = Fprev22 * Fprev22;
        // F t14 = 0.2e1 * t6 + t8 + 0.2e1 * t9 + t11 + 0.2e1 * t12;
        // F t15 = X1 * X1;
        // F t16 = t15 * t15;
        // F t20 = y2 - y1;
        // F t21 = Fprev12 * t20;
        // F t22 = Fprev21 * t20;
        // F t23 = -x2 + x1;
        // F t25 = 0.2e1 * t23 * Fprev11;
        // F t30 = X2 * X2;
        // F t35 = -Fprev12 * t20 - Fprev21 * t20 + t25;
        // F t38 = Y1 * Y1;
        // F t40 = 0.2e1 * t38 * t14;
        // F t41 = -Y2 * t14;
        // F t43 = -Fprev12 * t23;
        // F t44 = -Fprev21 * t23;
        // F t46 = -0.2e1 * t20 * Fprev22;
        // F t49 = 0.2e1 * Y1 * (0.2e1 * t41 + t43 + t44 - t46);
        // F t50 = Y2 * Y2;
        // F t52 = 0.2e1 * t50 * t14;
        // F t53 = Fprev12 * t23;
        // F t54 = Fprev21 * t23;
        // F t56 = Y2 * (t53 + t54 + t46);
        // F t57 = 0.2e1 * t56;
        // F t58 = x1 * x1;
        // F t59 = 0.2e1 * t58;
        // F t60 = x1 * x2;
        // F t61 = 0.4e1 * t60;
        // F t62 = x2 * x2;
        // F t63 = 0.2e1 * t62;
        // F t64 = t20 * t20;
        // F t67 = t30 * X2;
        // F t70 = t21 + t22 - t25;
        // F t87 = Y1 - Y2;
        // F t92 = (Y1 * t35 + Y2 * t70 + t20 * t23) * t87;
        // F t96 = t30 * t30;
        // F t106 = t6 + t8 / 0.2e1 + t9 + t11 / 0.2e1 + t12;
        // F t114 = t87 * t87;
        // F t122 = pow(-0.2e1 * X2 * X1 + t114 + t15 + t30, 0.2e1);
        // F t124 = dt * dt;

        // F unknown[1];

        // unknown[0] = 0.1e1 / t124 / t122 * 0.3141592654e1 * (t16 * t14 + 0.2e1 * t15 * X1 * (-0.2e1 * X2 * t14 + t21 + t22 - t25) + t15 * (0.6e1 * X2 * t35 + 0.6e1 * t14 * t30 + t40 + t49 + t52 + t57 + t59 - t61 + t63 + t64) + 0.2e1 * X1 * (-0.2e1 * t67 * t14 + 0.3e1 * t30 * t70 + X2 * (-0.2e1 * t38 * t14 + 0.2e1 * Y1 * (0.2e1 * Y2 * t14 + t46 + t53 + t54) - 0.2e1 * t50 * t14 + 0.2e1 * Y2 * (t43 + t44 - t46) - t59 + t61 - t63 - t64) - t92) + t96 * t14 + 0.2e1 * t67 * t35 + t30 * (t40 + t49 + t52 + t57 + t59 - t61 + t63 + t64) + 0.2e1 * X2 * t92 + 0.2e1 * t114 * (t38 * t106 + Y1 * (t41 + t43 + t44 - t46) + t50 * t106 + t56 + t58 / 0.2e1 - t60 + t62 / 0.2e1 + t64)) * t4 * eta * t2 / 0.4e1;

        
        return 0; // unknown[0];
    }
}

namespace GradientDynFunctions {

    void deformationViscous(const Particle2D& p,
        const Simulation& sim,
        F &fx, F &fy, F &fe) 
    {
        // F x1 = inputs[0];
        // F y1 = inputs[1];
        // F epsV1 = inputs[2];
        // F x2 = inputs[3];
        // F y2 = inputs[4];

        // F t1 = 0.1e1 + epsV1;
        // F t2 = t1 * t1;
        // F t4 = r0 * r0;
        // F t5 = t4 * eta * t2;
        // F t6 = X1 * X1;
        // F t7 = t6 * X1;
        // F t12 = -Fprev12 - Fprev21;
        // F t13 = Y1 * t12;
        // F t14 = 0.2e1 * t13;
        // F t15 = -Y2 * t12;
        // F t16 = 0.2e1 * t15;
        // F t17 = 0.4e1 * x1;
        // F t18 = 0.4e1 * x2;
        // F t21 = X2 * X2;
        // F t24 = -Y1 * t12;
        // F t26 = Y2 * t12;
        // F t30 = Y1 - Y2;
        // F t36 = (0.2e1 * Fprev11 * Y1 - 0.2e1 * Fprev11 * Y2 - y1 + y2) * t30;
        // F t40 = t21 * X2;
        // F t48 = t30 * t30;
        // F t56 = pow(-0.2e1 * X2 * X1 + t21 + t48 + t6, 0.2e1);
        // F t58 = dt * dt;
        // F t60 = 0.1e1 / t58 / t56;
        // F t68 = Fprev22 * Y1;
        // F t69 = 0.4e1 * t68;
        // F t70 = Fprev22 * Y2;
        // F t71 = 0.4e1 * t70;
        // F t72 = 0.2e1 * y2;
        // F t73 = 0.2e1 * y1;
        // F t78 = t69 - t71 - t73 + t72;
        // F t81 = (t24 + t26 + x2 - x1) * t30;
        // F t100 = Fprev11 * Fprev11;
        // F t102 = Fprev12 * Fprev12;
        // F t103 = Fprev12 * Fprev21;
        // F t105 = Fprev21 * Fprev21;
        // F t106 = Fprev22 * Fprev22;
        // F t108 = 0.2e1 * t100 + t102 + 0.2e1 * t103 + t105 + 0.2e1 * t106;
        // F t109 = t6 * t6;
        // F t113 = y2 - y1;
        // F t114 = Fprev12 * t113;
        // F t115 = Fprev21 * t113;
        // F t116 = -x2 + x1;
        // F t118 = 0.2e1 * t116 * Fprev11;
        // F t126 = -Fprev12 * t113 - Fprev21 * t113 + t118;
        // F t129 = Y1 * Y1;
        // F t131 = 0.2e1 * t129 * t108;
        // F t132 = -Y2 * t108;
        // F t134 = -Fprev12 * t116;
        // F t135 = -Fprev21 * t116;
        // F t137 = -0.2e1 * t113 * Fprev22;
        // F t140 = 0.2e1 * Y1 * (0.2e1 * t132 + t134 + t135 - t137);
        // F t141 = Y2 * Y2;
        // F t143 = 0.2e1 * t141 * t108;
        // F t144 = Fprev12 * t116;
        // F t145 = Fprev21 * t116;
        // F t147 = Y2 * (t144 + t145 + t137);
        // F t148 = 0.2e1 * t147;
        // F t149 = x1 * x1;
        // F t150 = 0.2e1 * t149;
        // F t151 = x1 * x2;
        // F t152 = 0.4e1 * t151;
        // F t153 = x2 * x2;
        // F t154 = 0.2e1 * t153;
        // F t155 = t113 * t113;
        // F t160 = t114 + t115 - t118;
        // F t181 = (Y1 * t126 + Y2 * t160 + t113 * t116) * t30;
        // F t185 = t21 * t21;
        // F t195 = t100 + t102 / 0.2e1 + t103 + t105 / 0.2e1 + t106;

        // F unknown[3];

        // unknown[0] = t60 * 0.3141592654e1 * (-0.4e1 * t7 * Fprev11 + t6 * (0.12e2 * Fprev11 * X2 + t14 + t16 + t17 - t18) + 0.2e1 * X1 * (-0.6e1 * t21 * Fprev11 + X2 * (0.2e1 * t24 + 0.2e1 * t26 - t17 + t18) - t36) + 0.4e1 * t40 * Fprev11 + t21 * (t14 + t16 + t17 - t18) + 0.2e1 * X2 * t36 + 0.2e1 * t48 * (t13 + t15 + x1 - x2)) * t5 / 0.4e1;
        // unknown[1] = t60 * 0.3141592654e1 * (0.2e1 * t7 * t12 + t6 * (-0.6e1 * X2 * t12 - t69 + t71 - t72 + t73) + 0.2e1 * X1 * (X2 * t78 + 0.3e1 * t21 * t12 - t81) - 0.2e1 * t40 * t12 - t21 * t78 + 0.2e1 * X2 * t81 + 0.4e1 * t48 * (-t68 + t70 + y1 - y2)) * t5 / 0.4e1;
        // unknown[2] = t60 * 0.3141592654e1 * (t109 * t108 + 0.2e1 * t7 * (-0.2e1 * X2 * t108 + t114 + t115 - t118) + t6 * (0.6e1 * X2 * t126 + 0.6e1 * t21 * t108 + t131 + t140 + t143 + t148 + t150 - t152 + t154 + t155) + 0.2e1 * X1 * (-0.2e1 * t40 * t108 + 0.3e1 * t21 * t160 + X2 * (-0.2e1 * t129 * t108 + 0.2e1 * Y1 * (0.2e1 * Y2 * t108 + t137 + t144 + t145) - 0.2e1 * t141 * t108 + 0.2e1 * Y2 * (t134 + t135 - t137) - t150 + t152 - t154 - t155) - t181) + t185 * t108 + 0.2e1 * t40 * t126 + t21 * (t131 + t140 + t143 + t148 + t150 - t152 + t154 + t155) + 0.2e1 * X2 * t181 + 0.2e1 * t48 * (t129 * t195 + Y1 * (t132 + t134 + t135 - t137) + t141 * t195 + t147 + t149 / 0.2e1 - t151 + t153 / 0.2e1 + t155)) * t4 * eta * t1 / 0.2e1;

        // fx =  unknown[0];
        // fy =  unknown[1];
        // fe =  unknown[2];
    }

}

namespace HessianDynFunctions {

    void deformationViscous(const Particle2D& p,
        const Simulation& sim,
        int i,
        AddFunc add)
    {
        // const int S   = sim.stride();
        // const int xIx = S*i+0, yIx = S*i+1, eIx = S*i+2;

        // F x1 = inputs[0];
        // F y1 = inputs[1];
        // F epsV1 = inputs[2];
        // F x2 = inputs[3];
        // F y2 = inputs[4];
    
        // F t1 = 0.1e1 + epsV1;
        // F t2 = t1 * t1;
        // F t4 = r0 * r0;
        // F t5 = t4 * eta * t2;
        // F t6 = X1 * X1;
        // F t8 = X2 * X1;
        // F t10 = X2 * X2;
        // F t12 = Y1 - Y2;
        // F t13 = t12 * t12;
        // F t19 = pow(t6 - 0.2e1 * t8 + t10 + t13, 0.2e1);
        // F t20 = 0.1e1 / t19;
        // F t21 = dt * dt;
        // F t22 = 0.1e1 / t21;
        // F t23 = t22 * t20;
        // F t33 = t23 * 0.3141592654e1 * (X1 * t12 - X2 * t12) * t5 / 0.2e1;
        // F t35 = t4 * eta * t1;
        // F t36 = t6 * X1;
        // F t41 = -Fprev12 - Fprev21;
        // F t42 = Y1 * t41;
        // F t43 = 0.2e1 * t42;
        // F t44 = -Y2 * t41;
        // F t45 = 0.2e1 * t44;
        // F t46 = 0.4e1 * x1;
        // F t47 = 0.4e1 * x2;
        // F t52 = -Y1 * t41;
        // F t54 = Y2 * t41;
        // F t63 = (0.2e1 * Fprev11 * Y1 - 0.2e1 * Fprev11 * Y2 - y1 + y2) * t12;
        // F t67 = t10 * X2;
        // F t81 = t23 * 0.3141592654e1 * (-0.4e1 * t36 * Fprev11 + t6 * (0.12e2 * Fprev11 * X2 + t43 + t45 + t46 - t47) + 0.2e1 * X1 * (-0.6e1 * t10 * Fprev11 + X2 * (0.2e1 * t52 + 0.2e1 * t54 - t46 + t47) - t63) + 0.4e1 * t67 * Fprev11 + t10 * (t43 + t45 + t46 - t47) + 0.2e1 * X2 * t63 + 0.2e1 * t13 * (t42 + t44 + x1 - x2)) * t35 / 0.2e1;
        // F t95 = Fprev22 * Y1;
        // F t96 = 0.4e1 * t95;
        // F t97 = Fprev22 * Y2;
        // F t98 = 0.4e1 * t97;
        // F t99 = 0.2e1 * y2;
        // F t100 = 0.2e1 * y1;
        // F t105 = t96 - t98 - t100 + t99;
        // F t108 = (t52 + t54 + x2 - x1) * t12;
        // F t124 = t23 * 0.3141592654e1 * (0.2e1 * t36 * t41 + t6 * (-0.6e1 * X2 * t41 + t100 - t96 + t98 - t99) + 0.2e1 * X1 * (X2 * t105 + 0.3e1 * t10 * t41 - t108) - 0.2e1 * t67 * t41 - t10 * t105 + 0.2e1 * X2 * t108 + 0.4e1 * t13 * (-t95 + t97 + y1 - y2)) * t35 / 0.2e1;
        // F t126 = Fprev11 * Fprev11;
        // F t128 = Fprev12 * Fprev12;
        // F t129 = Fprev12 * Fprev21;
        // F t131 = Fprev21 * Fprev21;
        // F t132 = Fprev22 * Fprev22;
        // F t134 = 0.2e1 * t126 + t128 + 0.2e1 * t129 + t131 + 0.2e1 * t132;
        // F t135 = t6 * t6;
        // F t139 = y2 - y1;
        // F t140 = Fprev12 * t139;
        // F t141 = Fprev21 * t139;
        // F t142 = -x2 + x1;
        // F t144 = 0.2e1 * t142 * Fprev11;
        // F t152 = -Fprev12 * t139 - Fprev21 * t139 + t144;
        // F t155 = Y1 * Y1;
        // F t157 = 0.2e1 * t155 * t134;
        // F t158 = -Y2 * t134;
        // F t160 = -Fprev12 * t142;
        // F t161 = -Fprev21 * t142;
        // F t163 = -0.2e1 * t139 * Fprev22;
        // F t166 = 0.2e1 * Y1 * (0.2e1 * t158 + t160 + t161 - t163);
        // F t167 = Y2 * Y2;
        // F t169 = 0.2e1 * t167 * t134;
        // F t170 = Fprev12 * t142;
        // F t171 = Fprev21 * t142;
        // F t173 = Y2 * (t170 + t171 + t163);
        // F t174 = 0.2e1 * t173;
        // F t175 = x1 * x1;
        // F t176 = 0.2e1 * t175;
        // F t177 = x1 * x2;
        // F t178 = 0.4e1 * t177;
        // F t179 = x2 * x2;
        // F t180 = 0.2e1 * t179;
        // F t181 = t139 * t139;
        // F t186 = t140 + t141 - t144;
        // F t207 = (Y1 * t152 + Y2 * t186 + t139 * t142) * t12;
        // F t211 = t10 * t10;
        // F t221 = t126 + t128 / 0.2e1 + t129 + t131 / 0.2e1 + t132;
    
        // F unknown[3][3];
    
        // unknown[0][0] = t23 * 0.3141592654e1 * (0.4e1 * t6 - 0.8e1 * t8 + 0.4e1 * t10 + 0.2e1 * t13) * t5 / 0.4e1;
        // unknown[0][1] = t33;
        // unknown[0][2] = t81;
        // unknown[1][0] = t33;
        // unknown[1][1] = t23 * 0.3141592654e1 * (0.2e1 * t6 - 0.4e1 * t8 + 0.2e1 * t10 + 0.4e1 * t13) * t5 / 0.4e1;
        // unknown[1][2] = t124;
        // unknown[2][0] = t81;
        // unknown[2][1] = t124;
        // unknown[2][2] = t22 * t20 * 0.3141592654e1 * (t135 * t134 + 0.2e1 * t36 * (-0.2e1 * X2 * t134 + t140 + t141 - t144) + t6 * (0.6e1 * X2 * t152 + 0.6e1 * t10 * t134 + t157 + t166 + t169 + t174 + t176 - t178 + t180 + t181) + 0.2e1 * X1 * (-0.2e1 * t67 * t134 + 0.3e1 * t10 * t186 + X2 * (-0.2e1 * t155 * t134 + 0.2e1 * Y1 * (0.2e1 * Y2 * t134 + t163 + t170 + t171) - 0.2e1 * t167 * t134 + 0.2e1 * Y2 * (t160 + t161 - t163) - t176 + t178 - t180 - t181) - t207) + t211 * t134 + 0.2e1 * t67 * t152 + t10 * (t157 + t166 + t169 + t174 + t176 - t178 + t180 + t181) + 0.2e1 * X2 * t207 + 0.2e1 * t13 * (t155 * t221 + Y1 * (t158 + t160 + t161 - t163) + t167 * t221 + t173 + t175 / 0.2e1 - t177 + t179 / 0.2e1 + t181)) * t4 * eta / 0.2e1;
    
        // processMapleOutput(reinterpret_cast<F *>(unknown), value.hessian, 3, 3);
    
    

        // // now scatter into the big H:
        // for(int a=0; a<3; ++a)
        //     for(int b=0; b<3; ++b)
        //         add(  /* idx[a] = { xIx,yIx,eIx } */[a],
        //         /* idx[b] */      [b],
        //         Hloc[a][b]);
    }
}