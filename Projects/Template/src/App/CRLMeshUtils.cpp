#include "Projects/Template/include/App/CRLMeshUtils.h"
#include "polyscope/surface_mesh.h"
#include <cmath>

namespace CRLMeshUtils {

void registerDiskMeshForParticle(const Particle2D& particle, const std::string& name) {
    std::vector<Vector3F> diskVertices;
    std::vector<Vector3I> diskFaces;
    int numSegments = 32;
    F radius = particle.radius;
    Vector3F center = particle.pos;

    // Add the center vertex.
    diskVertices.push_back(center);

    // Generate perimeter vertices using polar coordinates.
    for (int i = 0; i < numSegments; ++i) {
        F theta = 2.0 * M_PI * i / numSegments;
        Vector3F vertex = center + Vector3F(radius * cos(theta), radius * sin(theta), 0.0);
        diskVertices.push_back(vertex);
    }

    // Create faces (triangles) using a triangle fan from the center.
    for (int i = 1; i < numSegments; ++i) {
        diskFaces.push_back(Vector3I(0, i, i + 1));
    }
    // Close the circle.
    diskFaces.push_back(Vector3I(0, numSegments, 1));

    // Register the mesh with Polyscope.
    polyscope::registerSurfaceMesh(name, diskVertices, diskFaces);
}

} // namespace CRLMeshUtils
