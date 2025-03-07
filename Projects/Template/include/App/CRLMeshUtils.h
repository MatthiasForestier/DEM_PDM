#pragma once
#include "CRLHelper/VecMatDef.h"
#include "Projects/Template/include/Model/MassSpring.h"
#include "polyscope/surface_mesh.h"
#include <vector>
#include <string>

namespace CRLMeshUtils {
    void registerDiskMeshForParticle(const Particle2D& particle, const std::string& name);
}
