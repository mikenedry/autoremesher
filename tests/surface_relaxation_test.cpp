#include <AutoRemesher/SurfaceRelaxation>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace AutoRemesher;
using Polygons = std::vector<std::vector<size_t>>;
static void require(bool ok, const char* message)
{
    if (!ok)
        throw std::runtime_error(message);
}
struct Ridge {
    std::vector<Vector3> prepared, output;
    Polygons triangles, quads;
    std::vector<std::vector<Vector2>> uv;
    std::vector<SurfaceRelaxationFace> metric;
    Ridge()
    {
        // Cylinder cross-section with unequal angular and longitudinal samples.
        for (double y : { -1., 0., .6 })
            for (double angle : { -.4, 0., .7 })
                output.emplace_back(std::sin(angle), y, std::cos(angle));
        for (size_t y = 0; y < 2; ++y)
            for (size_t x = 0; x < 2; ++x)
                quads.push_back({ y * 3 + x, y * 3 + x + 1, (y + 1) * 3 + x + 1, (y + 1) * 3 + x });
        for (size_t y = 0; y < 3; ++y)
            for (size_t x = 0; x <= 220; ++x) {
                const double angle = -.4 + .005 * x;
                prepared.emplace_back(std::sin(angle), y == 0 ? -1 : y == 1 ? 0
                                                                            : .6,
                    std::cos(angle));
            }
        for (size_t y = 0; y < 2; ++y)
            for (size_t x = 0; x < 220; ++x) {
                size_t a = y * 221 + x, b = a + 1, c = a + 222, d = a + 221;
                triangles.push_back({ a, b, c });
                triangles.push_back({ a, c, d });
            }
        for (size_t y = 0; y < 3; ++y)
            for (size_t x = 0; x < 3; ++x)
                output[y * 3 + x] = prepared[y * 221 + (x == 0 ? 0 : x == 1 ? 80
                                                                            : 220)];
        for (const auto& triangle : triangles) {
            std::vector<Vector2> chart;
            for (size_t v : triangle)
                chart.emplace_back(prepared[v].y(), std::atan2(prepared[v].x(), prepared[v].z()));
            uv.push_back(chart);
            SurfaceRelaxationFace face;
            face.direction = Vector3(0, 1, 0);
            face.spacingU = 1;
            face.spacingV = .7;
            face.anisotropy = 1;
            face.strictCurvature = true;
            metric.push_back(face);
        }
    }
    std::vector<SurfaceRelaxationStencil> build(const std::map<std::pair<size_t, size_t>, size_t>* uses = nullptr) const
    {
        return buildSurfaceRelaxationStencils(prepared, triangles, uv, metric, output, quads, uses);
    }
};
static void sharedEdgeOwnership()
{
    // The two sides of a continuous chart intentionally carry unlike metrics.
    const Polygons triangles { { 0, 1, 2 }, { 0, 2, 3 } };
    const Polygons quads { { 0, 1, 4, 3 }, { 1, 2, 5, 4 }, { 3, 4, 7, 6 }, { 4, 5, 8, 7 } };
    const std::vector<std::vector<Vector2>> uv { { { -1, -1 }, { 1, -1 }, { 1, 1 } }, { { -1, -1 }, { 1, 1 }, { -1, 1 } } };
    for (double scale : { .25, 1., 7. })
        for (bool translated : { false, true }) {
            const auto transform = [&](const Vector3& p) {
                return scale * (Vector3(p.z(), p.x(), p.y()) + (translated ? Vector3(13, -9, 4) : Vector3()));
            };
            std::vector<Vector3> source { { -1, -1, 0 }, { 1, -1, 0 }, { 1, 1, 0 }, { -1, 1, 0 } };
            for (auto& p : source)
                p = transform(p);
            std::vector<SurfaceRelaxationFace> metric(2);
            for (auto& face : metric) {
                face.direction = Vector3(0, 1, 0);
                face.spacingU = face.spacingV = .4 * scale;
            }
            metric[0].strictCurvature = true;
            metric[0].anisotropy = 1;
            for (bool nearby : { true, false }) {
                Vector3 targets[2];
                double blends[2];
                for (size_t side = 0; side < 2; ++side) {
                    std::vector<Vector3> points { { -.6, -.6, 0 }, { 0, -.6, 0 }, { .6, -.6, 0 },
                        { -.2, 0, 0 }, { 0, 0, 0 }, { .7, 0, 0 }, { -.6, .6, 0 }, { 0, .3, 0 }, { .6, .6, 0 } };
                    const double step = nearby ? 4 * std::numeric_limits<double>::epsilon() * (translated ? 13 : 1) : .001;
                    const double offset = side ? step : -step;
                    points[4] = Vector3(offset, -offset, 0);
                    for (auto& p : points)
                        p = transform(p);
                    const auto stencil = buildSurfaceRelaxationStencils(source, triangles, uv, metric, points, quads);
                    require(stencil[4].valid, "shared-edge fan has no stencil");
                    targets[side] = stencil[4].target(points, 4);
                    blends[side] = stencil[4].blend;
                }
                if (nearby)
                    require((targets[0] - targets[1]).length() < 1e-10 * scale, "shared-edge roundoff changed metric target");
                else
                    require(blends[0] == 0 && blends[1] == 1, "distinct face interiors lost their metrics");
            }
        }
}

int main()
{
    try {
        sharedEdgeOwnership();
        Ridge ridge;
        const auto original = ridge.build();
        require(original.size() == ridge.output.size() && original[4].valid, "closed ridge fan has no stencil");
        require(std::count_if(original.begin(), original.end(), [](const auto& s) { return s.valid; }) == 1, "boundary vertex acquired a metric stencil");
        const Vector3 target = original[4].target(ridge.output, 4);
        const Vector3 uniform = (ridge.output[1] + ridge.output[3] + ridge.output[5] + ridge.output[7]) * .25;
        require(std::fabs(target.z() - 1) < .005, "metric target loses curved ridge height");
        require(target.z() > uniform.z() + .03, "curved ridge fixture does not discriminate uniform smoothing");
        std::cout << "ridge_target_z=" << target.z() << " uniform_z=" << uniform.z() << '\n';
        std::map<std::pair<size_t, size_t>, size_t> uses;
        for (const auto& f : ridge.quads)
            for (size_t k = 0; k < f.size(); ++k)
                ++uses[std::minmax(f[k], f[(k + 1) % f.size()])];
        require((ridge.build(&uses)[4].target(ridge.output, 4) - target).length() < 1e-12, "reusing output edge counts changes target");
        // A height-preserving no-op must fail: rebalance a skewed planar fan.
        std::vector<Vector3> plane;
        for (double y : { -1., 0., 1. })
            for (double x : { -1., 0., 1. })
                plane.emplace_back(x, y, 0);
        Polygons planeTriangles;
        std::vector<std::vector<Vector2>> planeUvs;
        std::vector<SurfaceRelaxationFace> planeMetric;
        for (const auto& q : ridge.quads) {
            planeTriangles.push_back({ q[0], q[1], q[2] });
            planeTriangles.push_back({ q[0], q[2], q[3] });
        }
        for (const auto& t : planeTriangles) {
            std::vector<Vector2> uv;
            for (size_t v : t)
                uv.emplace_back(plane[v].x(), plane[v].y());
            planeUvs.push_back(uv);
            SurfaceRelaxationFace face;
            face.direction = Vector3(1, 0, 0);
            face.anisotropy = 1;
            face.strictCurvature = true;
            planeMetric.push_back(face);
        }
        auto skewed = plane;
        skewed[4] = Vector3(0, .3, 0);
        const auto redistribute = buildSurfaceRelaxationStencils(plane, planeTriangles, planeUvs, planeMetric, skewed, ridge.quads);
        require(redistribute[4].valid, "skewed planar fan has no stencil");
        const Vector3 balanced = redistribute[4].target(skewed, 4);
        const double beforeImbalance = std::fabs((skewed[7] - skewed[4]).length() - (skewed[1] - skewed[4]).length());
        const double afterImbalance = std::fabs((skewed[7] - balanced).length() - (skewed[1] - balanced).length());
        require(balanced.y() < .2 && afterImbalance < .75 * beforeImbalance, "metric target did not redistribute skewed spacing");
        require(std::fabs(balanced.x()) < 1e-12 && std::fabs(balanced.z()) < 1e-12, "planar redistribution left its axis or plane");
        std::cout << "skewed_target_y=" << balanced.y() << " before_imbalance=" << beforeImbalance << " after_imbalance=" << afterImbalance << '\n';
        Ridge chartTurn = ridge;
        for (auto& face : chartTurn.uv)
            for (auto& uv : face)
                uv = Vector2(-uv.y(), uv.x());
        require((chartTurn.build()[4].target(ridge.output, 4) - target).length() < 1e-12, "chart quarter-turn changes target");
        Ridge fieldTurn = ridge;
        for (size_t f = 0; f < fieldTurn.metric.size(); ++f) {
            const auto& t = ridge.triangles[f];
            const Vector3 normal = Vector3::crossProduct(ridge.prepared[t[1]] - ridge.prepared[t[0]], ridge.prepared[t[2]] - ridge.prepared[t[0]]).normalized();
            fieldTurn.metric[f].direction = Vector3::crossProduct(normal, ridge.metric[f].direction);
            std::swap(fieldTurn.metric[f].spacingU, fieldTurn.metric[f].spacingV);
        }
        require((fieldTurn.build()[4].target(ridge.output, 4) - target).length() < 1e-12, "field quarter-turn with spacing swap changes target");
        for (double scale : { .001, 1000. }) {
            Ridge scaled = ridge;
            for (auto& p : scaled.prepared)
                p *= scale;
            for (auto& p : scaled.output)
                p *= scale;
            for (auto& f : scaled.metric) {
                f.spacingU *= scale;
                f.spacingV *= scale;
            }
            const auto stencil = scaled.build();
            require(stencil[4].valid && (stencil[4].target(scaled.output, 4) - target * scale).length() < 1e-11 * scale, "physical unit rescaling changes target");
        }
        Ridge reversed = ridge;
        std::reverse(reversed.quads.begin(), reversed.quads.end());
        for (auto& f : reversed.quads)
            std::reverse(f.begin(), f.end());
        require((reversed.build()[4].target(ridge.output, 4) - target).length() < 1e-12, "fan ordering changes target");
        Ridge boundary = ridge;
        boundary.quads.pop_back();
        require(!boundary.build()[4].valid, "four-valence open fan accepted");
        Ridge triangleStar = ridge;
        triangleStar.quads = { { 4, 5, 7 }, { 4, 7, 3 }, { 4, 3, 1 }, { 4, 1, 5 } };
        require(triangleStar.build()[4].valid, "closed four-triangle star rejected");
        require((triangleStar.build()[4].target(ridge.output, 4) - target).length() < 1e-12, "closed four-triangle star changes target");
        Ridge disconnected = ridge;
        disconnected.quads = { { 4, 1, 3 }, { 4, 3, 1 }, { 4, 5, 7 }, { 4, 7, 5 } };
        require(!disconnected.build()[4].valid, "disconnected four-valence fan accepted");
        Ridge badUv = ridge;
        for (auto& f : badUv.uv)
            for (auto& p : f)
                p = Vector2();
        const auto invalid = badUv.build();
        require(!invalid[4].valid, "singular UV should use old generic fallback");
        require((invalid[4].target(ridge.output, 4) - ridge.output[4]).length() == 0, "invalid stencil does not return current point");
        Ridge badMetric = ridge;
        for (auto& f : badMetric.metric)
            f.spacingU = std::numeric_limits<double>::quiet_NaN();
        require(!badMetric.build()[4].valid, "nonfinite spacing accepted");
        Ridge badDomain = ridge;
        badDomain.uv.pop_back();
        require(!badDomain.build()[4].valid, "mismatched face records accepted");
        std::cout << "curved ridge, chart/field gauge, physical unit rescaling, fan topology and invalid-input controls passed\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
