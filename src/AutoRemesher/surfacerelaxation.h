#ifndef AUTO_REMESHER_SURFACE_RELAXATION_H
#define AUTO_REMESHER_SURFACE_RELAXATION_H
#include <AutoRemesher/Vector2>
#include <AutoRemesher/Vector3>
#include <array>
#include <cstddef>
#include <map>
#include <vector>

namespace AutoRemesher {
struct SurfaceRelaxationFace {
    Vector3 direction;
    double spacingU = 1, spacingV = 1, anisotropy = 0;
    bool strictCurvature = false;
};
struct SurfaceRelaxationStencil {
    bool valid = false;
    double blend = 0;
    std::array<size_t, 4> neighbors {};
    std::array<Vector3, 4> outward;
    std::array<double, 4> minimum {}, maximum {};
    Vector3 target(const std::vector<Vector3>& points, size_t vertex) const;
};
// Samples the initial output once using nearest prepared faces and chart axes.
// Invalid charts and neighborhoods other than a closed valence-four fan retain
// the caller's existing target. Returned stencils own all their data.
std::vector<SurfaceRelaxationStencil> buildSurfaceRelaxationStencils(
    const std::vector<Vector3>& preparedVertices,
    const std::vector<std::vector<size_t>>& preparedTriangles,
    const std::vector<std::vector<Vector2>>& finalTriangleUvs,
    const std::vector<SurfaceRelaxationFace>& faces,
    const std::vector<Vector3>& outputVertices,
    const std::vector<std::vector<size_t>>& outputPolygons,
    const std::map<std::pair<size_t, size_t>, size_t>* edgeUseCount = nullptr);
}
#endif
