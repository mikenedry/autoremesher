// Geometry regressions and an optional headless OBJ driver.
#include <AutoRemesher/AutoRemesher>
#include <AutoRemesher/FrameField>
#include <AutoRemesher/IsotropicRemesher>
#include <AutoRemesher/MeshSeparator>
#include <AutoRemesher/MixedIntegerLeastSquares>
#include <AutoRemesher/QuadExtractor>
#include <AutoRemesher/Parameterizer>
#include <AutoRemesher/QuadParameterizer>
#include <AutoRemesher/SurfaceAnalysis>
#include <isotropichalfedgemesh.h>
#include <isotropicremesher.h>
#include <map>
#include <set>
#define TINYOBJLOADER_IMPLEMENTATION
#include "../src/tiny_obj_loader.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <type_traits>

using Remesher = AutoRemesher::AutoRemesher;
using Vec3 = AutoRemesher::Vector3;
using Faces = std::vector<std::vector<size_t>>;
using Feature = Remesher::ReferenceSurface::EdgeFeature;
using UV = AutoRemesher::Vector2;

static void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

namespace AutoRemesher {
// Test-only access to cycle extraction, before smoothing or hole repair.
struct QuadExtractorTestAccess {
    static const SurfaceAnalysis& workflow(bool open)
    {
        const std::vector<Vec3> p { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
        static const SurfaceAnalysis sheet(SurfaceMesh({ p[0], p[1], p[2] }, { { 0, 1, 2 } }),
            .1, 90, 0, 0, false, false, false, nullptr, true);
        static const SurfaceAnalysis solid(SurfaceMesh(p,
            { { 0, 2, 1 }, { 0, 1, 3 }, { 1, 2, 3 }, { 2, 0, 3 } }), .1, 90, 0, 0, false, false);
        return open ? sheet : solid;
    }
    static std::vector<Vec3> connections(const std::vector<Vec3>& points, const Faces& triangles,
        const std::vector<std::vector<UV>>& uv, bool cloth = true)
    {
        QuadExtractor extractor(&points, &triangles, &uv);
        extractor.setSurfaceAnalysis(&workflow(cloth));
        std::vector<Vec3> crossings;
        std::vector<size_t> owners;
        std::set<std::pair<size_t, size_t>> edges;
        extractor.extractConnections(&crossings, &owners, &edges);
        return crossings;
    }
    static Faces finish(const std::vector<Vec3>& points, const Faces& faces, bool cloth = true)
    {
        QuadExtractor extractor(&points, &faces, nullptr);
        extractor.setSurfaceAnalysis(&workflow(cloth));
        extractor.m_remeshedVertices = points;
        extractor.m_remeshedPolygons = faces;
        extractor.restoreBoundary();
        return extractor.m_remeshedPolygons;
    }
    static std::pair<std::vector<Vec3>, Faces> mergePentagons(const std::vector<Vec3>& source, const Faces& triangles,
        const SurfaceAnalysis& analysis, const std::vector<Vec3>& points, const Faces& faces)
    {
        QuadExtractor extractor(&source, &triangles, nullptr);
        extractor.setSurfaceAnalysis(&analysis);
        extractor.m_remeshedVertices = points;
        extractor.m_remeshedPolygons = faces;
        extractor.mergeSharedFiveEdgeFaces();
        return { extractor.m_remeshedVertices, extractor.m_remeshedPolygons };
    }

    static Faces splitSix(const std::vector<Vec3>& points, const Faces& faces, bool open = true)
    {
        std::vector<Vec3> source { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
        const Faces surface = open ? Faces { { 0, 1, 2 } }
                                   : Faces { { 0, 2, 1 }, { 0, 1, 3 }, { 1, 2, 3 }, { 2, 0, 3 } };
        if (open)
            source.resize(3);
        SurfaceAnalysis analysis(SurfaceMesh(source, surface), .1, 90, 0, 0, false, false, false, nullptr, true);
        require(analysis.openFoldWorkflow() == open, "hex fixture has the wrong source workflow");
        QuadExtractor extractor(&points, &faces, nullptr);
        extractor.setSurfaceAnalysis(&analysis);
        extractor.m_remeshedVertices = points;
        extractor.m_remeshedPolygons = faces;
        extractor.splitSixEdgeFaces();
        require(extractor.m_remeshedVertices.size() == points.size(), "hex split changed vertex count");
        for (size_t i = 0; i < points.size(); ++i)
            require((extractor.m_remeshedVertices[i] - points[i]).lengthSquared() == 0, "hex split moved a vertex");
        return extractor.m_remeshedPolygons;
    }
    using Graph = std::unordered_map<size_t, std::unordered_set<size_t>>;
    static Faces extract(std::vector<Vec3> points, const std::vector<size_t>& owners, Graph graph, bool open = true)
    {
        const std::vector<Vec3> source { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
        const Faces triangles { { 0, 2, 3 }, { 0, 3, 1 }, { 0, 1, 2 }, { 0, 2, 1 } };
        const std::vector<std::vector<UV>> uv;
        const Faces surface = open ? Faces { { 0, 1, 2 } }
                                   : Faces { { 0, 2, 1 }, { 0, 1, 3 }, { 1, 2, 3 }, { 2, 0, 3 } };
        auto support = source;
        if (open)
            support.resize(3);
        SurfaceAnalysis analysis(SurfaceMesh(support, surface), .1, 90, 0, 0, false, false, false, nullptr, true);
        QuadExtractor extractor(&source, &triangles, &uv);
        extractor.setSurfaceAnalysis(&analysis);
        Faces result;
        extractor.extractMesh(points, owners, graph, &result);
        return result;
    }
};
}

static void constantCoordinateIntersections()
{
    const std::vector<Vec3> points { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 } };
    const Faces triangles { { 0, 1, 2 } };
    for (bool swap : { false, true })
        for (double sign : { -1., 1. }) {
            std::vector<std::vector<UV>> uv { { { sign * 25, 32.769546213747894 },
                { sign * 25, 29.02313352630401 }, { sign * 25, 36.5159034436706 } } };
            if (swap)
                for (auto& p : uv[0])
                    p = UV(p.y(), p.x());
            const auto crossings = AutoRemesher::QuadExtractorTestAccess::connections(points, triangles, uv);
            // Seven integer levels cross the perimeter twice, plus three corners.
            require(crossings.size() == 17, "constant UV line gained a spurious crossing");
            for (const auto& p : crossings)
                require(p.x() < 1e-12 || p.y() < 1e-12 || p.x() + p.y() > 1 - 1e-12,
                    "constant UV line created an interior intersection");
        }
    const auto ordinary = AutoRemesher::QuadExtractorTestAccess::connections(points, triangles,
        { { { 0, 0 }, { 3, 0 }, { 0, 3 } } });
    require(std::any_of(ordinary.begin(), ordinary.end(), [](const Vec3& p) {
        return p.x() > .3 && p.y() > .3 && p.x() + p.y() < .7;
    }), "ordinary UV triangle lost its interior grid crossing");
}

static void curveConstrainedPentagonMerge()
{
    // A hard fold supports one endpoint of the pentagons' shared edge.
    // The ordinary midpoint leaves that fold, even though the merge stays manifold.
    const std::vector<Vec3> closedPoints { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
    const Faces closedTriangles { { 0, 2, 1 }, { 0, 1, 3 }, { 1, 2, 3 }, { 2, 0, 3 } };
    const AutoRemesher::SurfaceAnalysis closed(AutoRemesher::SurfaceMesh(closedPoints, closedTriangles), 1, 60, 0, 0, false, false);
    for (double scale : { .25, 1., 7. })
        for (bool rotate : { false, true })
            for (bool originalOpen : { false, true }) {
                const double angle = 75 * M_PI / 180;
                const auto point = [&](double x, double y) {
                    const Vec3 p = x <= 0 ? Vec3(x, y, 0) : Vec3(x * std::cos(angle), y, x * std::sin(angle));
                    return (rotate ? Vec3(p.z(), p.x(), p.y()) : p) * scale;
                };
                std::vector<Vec3> source;
                Faces triangles;
                for (int y = -1; y <= 1; ++y)
                    for (int x = -1; x <= 1; ++x)
                        source.push_back(point(2 * x, 2 * y));
                for (size_t y = 0; y < 2; ++y)
                    for (size_t x = 0; x < 2; ++x) {
                        const size_t a = 3 * y + x;
                        triangles.push_back({ a, a + 1, a + 4 });
                        triangles.push_back({ a, a + 4, a + 3 });
                    }
                const AutoRemesher::SurfaceAnalysis analysis(AutoRemesher::SurfaceMesh(source, triangles), scale, 60, 0, 0, true, false, false, originalOpen ? nullptr : &closed, true);
                const std::vector<Vec3> points { point(0, 0), point(.2, 0), point(.8, .6), point(.8, 1.2),
                    point(-.8, 1.2), point(-.8, -1.2), point(.8, -1.2), point(.8, -.6) };
                require(!analysis.supportsRimConstraints(), "crease merge fixture acquired a protected rim");
                require(analysis.bindCurve(points[0], .02 * scale).chain != AutoRemesher::SurfaceMesh::npos
                        && analysis.bindCurve(points[1], .02 * scale).chain == AutoRemesher::SurfaceMesh::npos,
                    "crease merge fixture lacks one curve-bound endpoint");
                const auto result = AutoRemesher::QuadExtractorTestAccess::mergePentagons(source, triangles, analysis,
                    points, { { 0, 1, 2, 3, 4 }, { 1, 0, 5, 6, 7 } });
                require(result.first.size() == 7 && result.second.size() == 2
                        && result.second[0].size() == 4 && result.second[1].size() == 4,
                    "valid crease-adjacent pentagons were not merged");
                double nearest = 1e100;
                for (const auto& p : result.first)
                    nearest = std::min(nearest, (p - point(0, 0)).length());
                require(originalOpen ? nearest < 1e-10 * scale : nearest > .05 * scale,
                    "pentagon merge changed curve support or legacy midpoint behavior");
            }
}

static void safeHexagonSplits()
{
    const auto checkSplit = [](const std::vector<Vec3>& points, const std::vector<size_t>& hex) {
        const auto children = AutoRemesher::QuadExtractorTestAccess::splitSix(points, { hex });
        require(children.size() == 2, "splittable hexagon lost its valid alternative");
        Vec3 normal;
        for (size_t i = 1; i + 1 < hex.size(); ++i)
            normal += Vec3::crossProduct(points[hex[i]] - points[hex[0]], points[hex[i + 1]] - points[hex[0]]);
        const auto turn = [&](size_t a, size_t b, size_t c) {
            return Vec3::dotProduct(Vec3::crossProduct(points[b] - points[a], points[c] - points[a]), normal);
        };
        for (const auto& quad : children) {
            require(quad.size() == 4, "hex split did not produce two quads");
            for (size_t i = 0; i < 2; ++i) {
                const size_t a = quad[i], b = quad[i + 1], c = quad[i + 2], d = quad[(i + 3) % 4];
                require(!(turn(a, b, c) * turn(a, b, d) < 0 && turn(c, d, a) * turn(c, d, b) < 0),
                    "hex split introduced a crossing quad perimeter");
            }
        }
    };
    // The best corner score chooses diagonal 2--5, which crosses this simple hexagon.
    // Either of the other diagonals can produce two valid quads without moving a vertex.
    const std::vector<Vec3> dent { { 8, 2, 0 }, { 5, 1.90, 0 }, { 0, 2, 0 },
        { 0, 0, 0 }, { 8.2, 0, 0 }, { 8.2, 1.9, 0 } };
    require(AutoRemesher::QuadExtractorTestAccess::splitSix(dent, { { 0, 1, 2, 3, 4, 5 } }, false)
            == Faces({ { 2, 3, 4, 5 }, { 5, 0, 1, 2 } }),
        "closed source changed the legacy hexagon score choice");
    for (double scale : { .25, 1.0, 7.0 })
        for (bool tilted : { false, true })
            for (bool reversed : { false, true })
                for (size_t start = 0; start < 6; ++start) {
                    auto points = dent;
                    for (auto& p : points)
                        p = (tilted ? Vec3(.8 * p.x(), p.y(), .6 * p.x()) : p) * scale;
                    std::vector<size_t> hex { 0, 1, 2, 3, 4, 5 };
                    std::rotate(hex.begin(), hex.begin() + start, hex.end());
                    if (reversed)
                        std::reverse(hex.begin(), hex.end());
                    checkSplit(points, hex);
                }
    std::vector<Vec3> convex;
    for (size_t i = 0; i < 6; ++i)
        convex.push_back({ std::cos(i * M_PI / 3), std::sin(i * M_PI / 3), 0 });
    checkSplit(convex, { 0, 1, 2, 3, 4, 5 });
    const std::vector<Vec3> line { { 0, 0, 0 }, { 1, 0, 0 }, { 2, 0, 0 }, { 3, 0, 0 }, { 4, 0, 0 }, { 5, 0, 0 } };
    const Faces unsplittable { { 0, 1, 2, 3, 4, 5 } };
    require(AutoRemesher::QuadExtractorTestAccess::splitSix(line, unsplittable) == unsplittable,
        "hexagon without a valid split was changed");
}

static void occupiedFaceInteriors()
{
    using Access = AutoRemesher::QuadExtractorTestAccess;
    const auto canonical = [](Faces faces) {
        for (auto& face : faces)
            std::rotate(face.begin(), std::min_element(face.begin(), face.end()), face.end());
        std::sort(faces.begin(), faces.end());
        return faces;
    };
    const auto extract = [&](const std::vector<Vec3>& points, const Faces& loops) {
        std::set<std::pair<size_t, size_t>> edges;
        for (const auto& face : loops)
            for (size_t i = 0; i < face.size(); ++i)
                edges.insert(std::minmax(face[i], face[(i + 1) % face.size()]));
        Access::Graph graph;
        for (const auto& edge : edges) {
            graph[edge.first].insert(edge.second);
            graph[edge.second].insert(edge.first);
        }
        return canonical(Access::extract(points, std::vector<size_t>(points.size(), 2), std::move(graph)));
    };
    // The nonplanar larger cycle cuts through the planar accepted quad
    // along edge3--4, creating an overlapping extra fan.
    const std::vector<Vec3> points { { -1, 1, 0 }, { -1, 0, 0 }, { -1, -.02, -.1 },
        { 0, -.1, 0 }, { 0, .1, 0 }, { 1, 0, 0 }, { 1, 1, 0 } };
    require(extract(points, { { 3, 5, 4, 1 }, { 0, 1, 2, 3, 4, 5, 6 } })
            == canonical({ { 3, 5, 4, 1 } }),
        "larger cycle overlapped an already extracted cell");
    // A triangle and a concave quad legitimately share a two-edge path.
    require(extract({ { 0, 0, 0 }, { .2, .5, 0 }, { 0, 1, 0 }, { 1, .5, 0 } },
                { { 0, 1, 2, 3 }, { 0, 1, 2 } })
            == canonical({ { 3, 2, 1, 0 }, { 0, 1, 2 } }),
        "connected shared boundary path lost a cell");
    // An exterior diagonal of the concave cell must remain harmless.
    require(extract({ { 0, 0, 0 }, { 2, 0, 0 }, { .9, .9, 0 }, { 0, 2, 0 }, { 2, 2, 0 } },
                { { 0, 1, 2, 3 }, { 1, 4, 3 } })
            == canonical({ { 4, 3, 2, 1 }, { 2, 3, 0, 1 } }),
        "exterior graph chord removed a valid concave cell");
    // A cell can surround an earlier quad without crossing its interior.
    require(extract({ { -2, 2, 0 }, { 2, 2, 0 }, { 2, -2, 0 }, { 1, -1, 0 },
                { .5, 0, 0 }, { -2, -2, 0 }, { 1, 1, 0 } },
                { { 6, 1, 2, 3 }, { 1, 6, 4, 3, 2, 5, 0 } })
            == canonical({ { 2, 1, 6, 3 }, { 4, 6, 1, 0, 5, 2, 3 }, { 4, 3, 6 } }),
        "surrounding cell with no interior chord was removed");
}

static void subdividedExtractionCorners()
{
    using Access = AutoRemesher::QuadExtractorTestAccess;
    const auto canonical = [](Faces faces) {
        for (auto& face : faces)
            std::rotate(face.begin(), std::min_element(face.begin(), face.end()), face.end());
        std::sort(faces.begin(), faces.end());
        return faces; // Retain duplicates and winding so either causes comparison failure.
    };
    const auto graphFor = [](const Faces& faces) {
        std::set<std::pair<size_t, size_t>> edges;
        for (const auto& face : faces)
            for (size_t i = 0; i < face.size(); ++i)
                edges.insert(std::minmax(face[i], face[(i + 1) % face.size()]));
        Access::Graph graph;
        for (const auto& edge : edges) {
            graph[edge.first].insert(edge.second);
            graph[edge.second].insert(edge.first);
        }
        return graph;
    };

    // Five cells meet at 0. Nodes 1 and 2 only subdivide two shared edges.
    // Reserving either undirected corner wrongly omits one of its incident cells.
    std::vector<Vec3> points { { 0, 0, 0 }, { -.25, 0, 0 }, { .75, 0, 0 },
        { 0, 0, -1 }, { 0, 1, 0 }, { 1, 0, 0 }, { 0, -1, 0 }, { -1, 0, 0 },
        { 0, 1, -1 }, { -1, -1, 0 }, { -1, 1, 0 }, { 1, -1, 0 }, { 1, 0, -1 } };
    // Source-triangle owners provide +X, +Y and +Z normals at the sharp corner.
    std::vector<size_t> owners { 2, 2, 2, 1, 2, 1, 2, 2, 0, 2, 2, 2, 1 };
    const Faces expected { { 4, 0, 3, 8 }, { 11, 5, 2, 0, 6 }, { 10, 7, 1, 0, 4 },
        { 0, 1, 7, 9, 6 }, { 0, 2, 5, 12, 3 } };
    auto graph = graphFor(expected);
    require(canonical(Access::extract(points, owners, graph)) == canonical(expected),
        "subdivided extraction edge lost a cell or emitted a duplicate/reversed cell");
    require(Access::extract(points, owners, graph, false).size() == 3,
        "open-fold recovery changed closed-source corner ownership");

    // A true junction still reserves its undirected corner.
    for (size_t node : { size_t(1), size_t(2) }) {
        graph[node].insert(points.size());
        graph[points.size()].insert(node);
        points.push_back(points[node] + Vec3(0, 0, 1));
        owners.push_back(2);
    }
    require(Access::extract(points, owners, graph).size() == 3,
        "higher-degree extraction corner reservation was bypassed");

    // An all-degree-two loop still emits exactly one source-facing cell.
    const std::vector<Vec3> square { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 } };
    for (size_t normal : { size_t(2), size_t(3) }) {
        Faces cell { { 0, 1, 2, 3 } };
        if (normal == 3)
            std::reverse(cell[0].begin(), cell[0].end());
        require(canonical(Access::extract(square, std::vector<size_t>(4, normal), graphFor(cell))) == canonical(cell),
            "degree-two loop duplicated a face or crossed the source sheet");
    }
}

static void collapsedOutputFaces()
{
    for (double scale : { 1., 7. }) {
        std::vector<Vec3> points { { 0, 0, 0 }, { 1, 0, 0 }, { 2, 0, 0 },
            { 3, 0, 0 }, { 4, 0, 0 }, { 1, 1e-180, 0 } };
        for (auto& p : points)
            p = scale * Vec3(p.z(), p.x(), p.y());
        const Faces valid { { 0, 1, 5 } };
        require(AutoRemesher::QuadExtractorTestAccess::finish(points,
                    { { 0, 1, 2, 3, 4 }, valid[0] }) == valid,
            "collapsed output face survived or a thin valid face was removed");
        const Faces input { { 0, 1, 2, 3, 4 }, valid[0] };
        require(AutoRemesher::QuadExtractorTestAccess::finish(points, input, false) == input,
            "cloth cleanup changed legacy closed-face retention");
    }
}

static void tube(size_t sides, size_t rings, std::vector<Vec3>& p, Faces& t,
    std::vector<std::vector<UV>>& uv)
{
    for (size_t j = 0; j < rings; ++j)
        for (size_t i = 0; i < sides; ++i)
            p.emplace_back(std::cos(2 * M_PI * i / sides), std::sin(2 * M_PI * i / sides), double(j));
    for (size_t j = 0; j + 1 < rings; ++j)
        for (size_t i = 0; i < sides; ++i) {
            const size_t a = j * sides + i, b = j * sides + (i + 1) % sides, c = b + sides, d = a + sides;
            t.push_back({ a, b, c });
            uv.push_back({ { double(i), double(j) }, { double(i + 1), double(j) }, { double(i + 1), double(j + 1) } });
            t.push_back({ a, c, d });
            uv.push_back({ { double(i), double(j) }, { double(i + 1), double(j + 1) }, { double(i), double(j + 1) } });
        }
}

static std::map<std::pair<size_t, size_t>, size_t> edgeUses(const Faces& faces)
{
    std::map<std::pair<size_t, size_t>, size_t> uses;
    for (const auto& f : faces)
        for (size_t k = 0; k < f.size(); ++k)
            ++uses[std::minmax(f[k], f[(k + 1) % f.size()])];
    return uses;
}

static void protectedPreparationRims()
{
    // A subdivided 90-degree rim must survive classification, flips and relaxation.
    for (double scale : { 1., .01 }) {
        const size_t sides = 32;
        std::vector<::Vector3> points;
        Faces triangles;
        for (double y : { -1., 1. })
            for (size_t i = 0; i < sides; ++i)
                points.emplace_back(scale * std::cos(2 * M_PI * i / sides), scale * y,
                    scale * std::sin(2 * M_PI * i / sides));
        for (size_t i = 0; i < sides; ++i) {
            const size_t j = (i + 1) % sides;
            triangles.push_back({ i, i + sides, j + sides });
            triangles.push_back({ i, j + sides, j });
        }
        for (size_t i = 1; i + 1 < sides; ++i) {
            triangles.push_back({ 0, i, i + 1 });
            triangles.push_back({ sides, sides + i + 1, sides + i });
        }
        ::IsotropicRemesher remesher(&points, &triangles);
        remesher.setTargetEdgeLength(.1 * scale);
        remesher.setSharpEdgeIncludedAngle(90);
        remesher.remesh(3);
        auto* mesh = remesher.remeshedHalfedgeMesh();
        for (auto* v = mesh->moveToNextVertex(nullptr); v; v = mesh->moveToNextVertex(v))
            require(!v->position.containsNan() && !v->position.containsInf(), "nonfinite prepared rim");
        // Every original rim segment must still be covered by working edges.
        for (size_t i = 0; i < points.size(); ++i) {
            const auto a = points[i], edge = points[i / sides * sides + (i + 1) % sides] - a;
            const auto onSegment = [&](const ::Vector3& p) {
                const double t = ::Vector3::dotProduct(p - a, edge) / edge.lengthSquared();
                return t >= -1e-9 && t <= 1 + 1e-9 && (p - a - edge * t).length() < 1e-9 * scale;
            };
            double covered = 0;
            for (auto* f = mesh->moveToNextFace(nullptr); f; f = mesh->moveToNextFace(f)) {
                auto* h = f->halfedge;
                for (size_t k = 0; k < 3; ++k, h = h->nextHalfedge)
                    if (onSegment(h->startVertex->position) && onSegment(h->nextHalfedge->startVertex->position))
                        covered += .5 * (h->startVertex->position - h->nextHalfedge->startVertex->position).length();
            }
            require(std::fabs(covered - edge.length()) < 1e-7 * scale, "preparation erased part of a sharp rim");
        }
    }
}

static void balancedRefinement()
{
    // A folded strip with alternating diagonals has avoidable 4/8-valence
    // fans. Refinement may balance its planar patches, preserving the crease.
    std::vector<Vec3> vertices;
    Faces triangles;
    for (size_t y = 0; y < 5; ++y)
        for (size_t x = 0; x < 7; ++x)
            vertices.emplace_back(double(x), double(std::min(y, size_t(2))), double(y > 2 ? y - 2 : 0));
    for (size_t y = 0; y < 4; ++y)
        for (size_t x = 0; x < 6; ++x) {
            const size_t a = y * 7 + x, b = a + 1, c = a + 7, d = c + 1;
            if ((x + y) % 2) {
                triangles.push_back({ a, b, c });
                triangles.push_back({ b, d, c });
            } else {
                triangles.push_back({ a, b, d });
                triangles.push_back({ a, d, c });
            }
        }
    const auto defect = [&](const Faces& faces) {
        const auto uses = edgeUses(faces);
        std::vector<size_t> valence(vertices.size());
        std::set<size_t> boundary;
        for (const auto& e : uses) {
            ++valence[e.first.first];
            ++valence[e.first.second];
            if (e.second == 1) {
                boundary.insert(e.first.first);
                boundary.insert(e.first.second);
            }
        }
        double result = 0;
        for (size_t v = 0; v < valence.size(); ++v)
            result += std::pow(double(valence[v]) - (boundary.count(v) ? 4 : 6), 2);
        return result;
    };
    AutoRemesher::IsotropicRemesher remesher(vertices, triangles);
    remesher.setRefineOnly(true);
    remesher.setBalanceDiagonals(true);
    remesher.setTargetEdgeLength(10);
    remesher.setSharpEdgeDegrees(90);
    require(remesher.remesh(), "strip refinement failed");
    const auto& out = remesher.remeshedVertices();
    const auto& faces = remesher.remeshedTriangles();
    require(out.size() == vertices.size() && faces.size() == triangles.size(), "diagonal balancing changed counts");
    for (size_t v = 0; v < out.size(); ++v)
        require((out[v] - vertices[v]).lengthSquared() == 0, "diagonal balancing moved a point");
    require(defect(faces) < defect(triangles), "refinement left avoidable diagonal fans");
    const auto before = edgeUses(triangles), after = edgeUses(faces);
    for (const auto& e : before)
        if (e.second == 1 || (e.first.first / 7 == 2 && e.first.second / 7 == 2)) {
            const auto found = after.find(e.first);
            require(found != after.end() && found->second == e.second, "balancing changed boundary or crease");
        }
    double area = 0;
    for (const auto& f : faces) {
        const Vec3 n = Vec3::crossProduct(out[f[1]] - out[f[0]], out[f[2]] - out[f[0]]);
        require(n.lengthSquared() > 0, "balancing created degenerate triangle");
        require((n.x() == 0 && n.y() == 0 && n.z() > 0) || (n.x() == 0 && n.y() < 0 && n.z() == 0), "balancing crossed folded sheets");
        area += .5 * n.length();
    }
    require(std::fabs(area - 24) < 1e-12, "balancing changed source area");
}

static void checkReferences(Remesher& remesher, const std::vector<Vec3>& vertices,
    const Faces& faces)
{
    std::vector<bool> found(faces.size(), false);
    size_t vertexOffset = 0, triangleOffset = 0;
    for (const auto& prepared : remesher.preparedIslands()) {
        require(bool(prepared.reference), "missing reference");
        const auto& ref = *prepared.reference;
        require(bool(prepared.analysis), "missing persistent source analysis");
        require(prepared.analysis->faces().size() == ref.triangles.size(), "source metric domain");
        for (const auto& chain : prepared.analysis->chains()) {
            require(chain.first < ref.vertices.size() && chain.last < ref.vertices.size(), "chain endpoint provenance");
            for (size_t corner : chain.corners)
                require(corner < 3 * ref.triangles.size(), "chain corner provenance");
        }
        static_assert(std::is_const<typename std::remove_reference<decltype(ref)>::type>::value,
            "reference surface must be immutable");
        require(ref.vertices.size() == ref.sourceVertexIds.size(), "vertex provenance size");
        require(ref.triangles.size() == ref.sourceTriangleIds.size(), "triangle provenance size");
        require(ref.edgeFeatures.size() == 3 * ref.triangles.size(), "feature provenance size");
        for (size_t i = 0; i < ref.vertices.size(); ++i) {
            require(ref.sourceVertexIds[i] < vertices.size(), "source vertex domain");
            require((ref.vertices[i] - vertices[ref.sourceVertexIds[i]]).lengthSquared() == 0, "reference position changed");
        }
        for (size_t i = 0; i < ref.triangles.size(); ++i) {
            const size_t source = ref.sourceTriangleIds[i];
            require(source < faces.size() && !found[source], "lost or duplicated source triangle");
            found[source] = true;
            for (size_t corner = 0; corner < 3; ++corner) {
                const size_t local = ref.triangles[i][corner];
                require(local < ref.vertices.size(), "reference vertex domain");
                require(ref.sourceVertexIds[local] == faces[source][corner], "source corner order changed");
            }
        }
        require(prepared.vertexOffset == vertexOffset && prepared.triangleOffset == triangleOffset,
            "prepared island offsets");
        vertexOffset += prepared.vertexCount;
        triangleOffset += prepared.triangleCount;
        require(triangleOffset <= remesher.isotropicTriangles().size(), "prepared triangle range");
        for (size_t i = prepared.triangleOffset; i < triangleOffset; ++i)
            for (size_t v : remesher.isotropicTriangles()[i])
                require(v >= prepared.vertexOffset && v < vertexOffset, "working triangle crosses reference islands");
    }
    for (bool seen : found)
        require(seen, "source triangle missing after preparation");
    require(vertexOffset == remesher.isotropicVertices().size(), "unmapped working vertices");
    require(triangleOffset == remesher.isotropicTriangles().size(), "unmapped working triangles");
    require(remesher.isotropicTriangleUvs().size() == triangleOffset, "selected candidate UV ranges do not match preparation");
}

static void fixtures()
{
    // Two coincident, disconnected folded strips with different input indices.
    // Both contain a 90-degree crease and four open boundary edges.
    std::vector<Vec3> vertices { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 },
        { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
    const auto original = vertices;
    const Faces faces { { 0, 1, 2 }, { 1, 0, 3 }, { 4, 5, 6 }, { 5, 4, 7 } };
    Remesher remesher(vertices, faces);
    vertices[0] = { 99, 99, 99 }; // The caller's array is not the reference owner.
    remesher.setTargetTriangleCount(100);
    remesher.setSharpEdgeDegrees(30);
    remesher.remesh();
    checkReferences(remesher, original, faces);
    require(remesher.preparedIslands().size() == 2, "coincident islands aliased");
    for (const auto& prepared : remesher.preparedIslands()) {
        const auto& ref = *prepared.reference;
        size_t sharp = 0, boundary = 0;
        for (Feature edge : ref.edgeFeatures) {
            sharp += edge == Feature::Sharp;
            boundary += edge == Feature::Boundary;
        }
        require(sharp == 2 && boundary == 4, "original crease or boundary classification lost");
        require(ref.sharpEdgeDegrees == 30, "feature threshold provenance");
        require(prepared.vertexCount != ref.vertices.size(), "fixture did not resample");
    }
    const auto retained = remesher.preparedIslands().front().reference;
    remesher.setSharpEdgeDegrees(120);
    remesher.remesh();
    checkReferences(remesher, original, faces);
    require(retained->sharpEdgeDegrees == 30, "previous snapshot mutated on rerun");
    for (Feature edge : remesher.preparedIslands().front().reference->edgeFeatures)
        require(edge != Feature::Sharp, "new threshold not captured");
    remesher.setTargetTriangleCount(0);
    require(!remesher.remesh() && remesher.preparedIslands().empty(), "failed run exposes stale preparation");
    require(remesher.remeshedVertices().empty() && remesher.remeshedQuads().empty() && remesher.isotropicVertices().empty() && remesher.isotropicTriangles().empty() && remesher.isotropicTriangleUvs().empty(), "failed run exposes stale delivery");
    require(retained->vertices[0][0] == original[0][0], "retained reference expired");

    // Equal face values must retain distinct input row IDs.
    const Faces duplicates { { 0, 1, 2 }, { 0, 1, 2 } };
    std::vector<Faces> islands;
    std::vector<std::vector<size_t>> ids;
    AutoRemesher::MeshSeparator::splitToIslands(duplicates, islands, &ids);
    require(ids.size() == 2 && ids[0][0] == 0 && ids[1][0] == 1, "duplicate source face identity");
}

static void tinyIsland()
{
    const std::vector<Vec3> cube { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 },
        { 0, 0, 1 }, { 1, 0, 1 }, { 1, 1, 1 }, { 0, 1, 1 } };
    const Faces quads { { 0, 3, 2, 1 }, { 4, 5, 6, 7 }, { 0, 1, 5, 4 }, { 1, 2, 6, 5 }, { 2, 3, 7, 6 }, { 3, 0, 4, 7 } };
    std::vector<Vec3> vertices = cube;
    Faces triangles;
    for (const auto& p : cube)
        vertices.push_back(p * .01 + Vec3(3, 0, 0));
    for (size_t offset : { size_t(0), size_t(8) })
        for (const auto& q : quads) {
            triangles.push_back({ q[0] + offset, q[1] + offset, q[2] + offset });
            triangles.push_back({ q[0] + offset, q[2] + offset, q[3] + offset });
        }
    Remesher remesher(vertices, triangles);
    remesher.setTargetTriangleCount(2000);
    remesher.setSharpEdgeDegrees(80);
    require(remesher.remesh(), "tiny island regression failed to remesh");
    bool large = false, small = false;
    for (const auto& q : remesher.remeshedQuads()) {
        const auto& p = remesher.remeshedVertices()[q[0]];
        large = large || p.x() < 2;
        small = small || p.x() > 2;
    }
    require(large && small, "global sizing removed the tiny disconnected part");
}

static void thinTubeCleanup()
{
    // Six-sided capped tube: valence-only cleanup used to unzip it from an end.
    std::vector<Vec3> p;
    Faces t;
    std::vector<std::vector<UV>> uv;
    const size_t rings = 25;
    for (size_t j = 0; j < rings; ++j)
        for (size_t i = 0; i < 6; ++i)
            p.emplace_back(double(j), .1 * std::cos(i * M_PI / 3), .1 * std::sin(i * M_PI / 3));
    const auto quad = [&](const std::vector<size_t>& q, const std::vector<UV>& u) {
        for (size_t k = 1; k < 3; ++k) {
            t.push_back({ q[0], q[k], q[k + 1] });
            uv.push_back({ u[0], u[k], u[k + 1] });
        }
    };
    for (size_t j = 0; j + 1 < rings; ++j)
        for (size_t i = 0; i < 6; ++i)
            quad({ 6 * j + i, 6 * j + (i + 1) % 6, 6 * (j + 1) + (i + 1) % 6, 6 * (j + 1) + i },
                { { double(i), double(j) }, { double(i + 1), double(j) }, { double(i + 1), double(j + 1) }, { double(i), double(j + 1) } });
    quad({ 5, 4, 1, 0 }, { { 0, 0 }, { 1, 0 }, { 1, 1 }, { 0, 1 } });
    quad({ 4, 3, 2, 1 }, { { 1, 0 }, { 2, 0 }, { 2, 1 }, { 1, 1 } });
    size_t o = 6 * (rings - 1);
    quad({ o, o + 1, o + 4, o + 5 }, { { 0, 0 }, { 1, 0 }, { 1, 1 }, { 0, 1 } });
    quad({ o + 1, o + 2, o + 3, o + 4 }, { { 1, 0 }, { 2, 0 }, { 2, 1 }, { 1, 1 } });
    AutoRemesher::QuadExtractor e(&p, &t, &uv);
    require(e.extract(), "thin tube extraction failed");
    double area = 0;
    const auto& v = e.remeshedVertices();
    for (const auto& f : e.remeshedQuads())
        for (size_t k = 1; k + 1 < f.size(); ++k)
            area += .5 * Vec3::crossProduct(v[f[k]] - v[f[0]], v[f[k + 1]] - v[f[0]]).length();
    require(area > .8 * .6 * (rings - 1), "valence cleanup flattened the thin tube");
    // Multiple extracted patches with source ownership must both survive cleanup.
    const size_t vertices = p.size(), triangles = t.size();
    for (size_t i = 0; i < vertices; ++i)
        p.push_back(p[i] + Vec3(0, 1, 0));
    for (size_t i = 0; i < triangles; ++i) {
        auto face = t[i];
        for (auto& v : face)
            v += vertices;
        t.push_back(face);
        uv.push_back(uv[i]);
    }
    AutoRemesher::SurfaceAnalysis analysis(AutoRemesher::SurfaceMesh(p, t), .1, 90, 1, 1);
    AutoRemesher::QuadExtractor both(&p, &t, &uv);
    both.setSurfaceAnalysis(&analysis);
    require(both.extract(), "fragment regression failed to extract");
    bool lower = false, upper = false;
    for (const auto& f : both.remeshedQuads()) {
        const double y = both.remeshedVertices()[f[0]].y();
        lower = lower || y < .3;
        upper = upper || y > .7;
    }
    require(lower && upper, "cleanup discarded an extracted source patch");
}

static void telescopingTube()
{
    // A four-node circumference is a graph cycle, but is not a surface cap.
    // The former extractor capped two shoulder rings and split this single tube.
    std::vector<Vec3> p;
    Faces f;
    std::vector<std::vector<UV>> uv;
    const double x[] = { 0, 1, 1.02, 2, 2.02, 3 }, r[] = { .4, .4, .3, .3, .2, .2 };
    for (size_t j = 0; j < 6; ++j)
        for (size_t i = 0; i < 4; ++i)
            p.push_back({ x[j], r[j] * std::cos(i * M_PI / 2), r[j] * std::sin(i * M_PI / 2) });
    const auto quad = [&](const std::vector<size_t>& q, const std::vector<UV>& u) {
        for (size_t k = 1; k < 3; ++k) {
            f.push_back({ q[0], q[k], q[k + 1] });
            uv.push_back({ u[0], u[k], u[k + 1] });
        }
    };
    // Shoulders occur first so physical shoulder normals seed shared source corners.
    for (size_t j : { 1, 3, 0, 2, 4 })
        for (size_t i = 0; i < 4; ++i)
            quad({ 4 * j + i, 4 * j + (i + 1) % 4, 4 * (j + 1) + (i + 1) % 4, 4 * (j + 1) + i },
                { { double(i), double(j) }, { double(i + 1), double(j) }, { double(i + 1), double(j + 1) }, { double(i), double(j + 1) } });
    quad({ 3, 2, 1, 0 }, { { 0, 0 }, { 1, 0 }, { 1, 1 }, { 0, 1 } });
    quad({ 20, 21, 22, 23 }, { { 0, 0 }, { 1, 0 }, { 1, 1 }, { 0, 1 } });
    AutoRemesher::SurfaceAnalysis a(AutoRemesher::SurfaceMesh(p, f), .5, 90, 0, 0, true, false);
    AutoRemesher::QuadExtractor q(&p, &f, &uv);
    q.setSurfaceAnalysis(&a);
    require(q.extract(true), "telescoping tube extraction failed");
    std::vector<Faces> components;
    AutoRemesher::MeshSeparator::splitToIslands(q.remeshedQuads(), components);
    require(components.size() == 1, "source-unsupported shoulder caps severed a telescoping tube");
    for (const auto& face : q.remeshedQuads())
        for (size_t k = 0; k < face.size(); ++k)
            require((q.remeshedVertices()[face[k]] - q.remeshedVertices()[face[(k + 1) % face.size()]]).lengthSquared() > 0,
                "telescoping tube contains a collapsed edge");
}

static void disconnectedFans()
{
    // Two affine 3x3 sheets touch only at the origin. Graph simplification
    // may trim their outer corners, but cleanup must retain the shared-origin
    // neighborhood on both sheets and give each a separate boundary fan.
    std::vector<Vec3> p { { 0, 0, 0 } };
    Faces t;
    std::vector<std::vector<UV>> uv;
    for (int sign : { -1, 1 }) {
        size_t grid[4][4];
        for (size_t j = 0; j < 4; ++j)
            for (size_t i = 0; i < 4; ++i) {
                grid[j][i] = i || j ? p.size() : 0;
                if (i || j)
                    p.push_back({ sign * double(i), sign * double(j), 0 });
            }
        for (size_t j = 0; j < 3; ++j)
            for (size_t i = 0; i < 3; ++i) {
                const size_t a = grid[j][i], b = grid[j][i + 1], c = grid[j + 1][i + 1], d = grid[j + 1][i];
                t.push_back({ a, b, c });
                uv.push_back({ { double(i), double(j) }, { double(i + 1), double(j) }, { double(i + 1), double(j + 1) } });
                t.push_back({ a, c, d });
                uv.push_back({ { double(i), double(j) }, { double(i + 1), double(j + 1) }, { double(i), double(j + 1) } });
            }
    }
    AutoRemesher::SurfaceAnalysis analysis(AutoRemesher::SurfaceMesh(p, t), 1, 90, 0, 0);
    AutoRemesher::QuadExtractor extractor(&p, &t, &uv);
    extractor.setSurfaceAnalysis(&analysis);
    require(extractor.extract(), "disconnected fans failed to extract");
    const auto& points = extractor.remeshedVertices();
    const auto& faces = extractor.remeshedQuads();
    std::map<size_t, unsigned> originSides;
    std::map<std::pair<size_t, size_t>, size_t> edges;
    for (const auto& face : faces) {
        require(face.size() == 4, "fan repair did not produce quads");
        Vec3 center;
        for (size_t v : face)
            center += points[v];
        require(center.x() * center.y() > 0, "fan repair connected different source sheets");
        for (size_t k = 0; k < face.size(); ++k) {
            const size_t v = face[k];
            if (points[v].lengthSquared() < 1e-18)
                originSides[v] |= center.x() > 0 ? 1 : 2;
            ++edges[std::minmax(v, face[(k + 1) % face.size()])];
        }
    }
    require(originSides.size() == 2, "fan repair discarded the shared-origin neighborhoods");
    unsigned sides = 0;
    for (const auto& v : originSides) {
        require(v.second == 1 || v.second == 2, "fan repair left the source sheets joined");
        sides |= v.second;
    }
    require(sides == 3, "fan repair discarded one source sheet");
    std::map<size_t, size_t> boundaryDegree;
    for (const auto& edge : edges) {
        require(edge.second <= 2, "fan repair created a nonmanifold edge");
        if (edge.second == 1) {
            ++boundaryDegree[edge.first.first];
            ++boundaryDegree[edge.first.second];
        }
    }
    for (const auto& v : boundaryDegree)
        require(v.second == 2, "fan repair left touching boundary fans");
    for (const auto& v : originSides)
        require(boundaryDegree[v.first] == 2, "fan repair capped a source boundary");
}

static void rotatedBoundaryCover()
{
    const std::vector<Vec3> points = { { 0, 0, 0 }, { 4, 0, 0 }, { 4, 2, 0 }, { 0, 2, 0 }, { .07, .73, 0 } };
    const Faces triangles = { { 0, 1, 4 }, { 1, 2, 4 }, { 2, 3, 4 }, { 3, 0, 4 } };
    const AutoRemesher::SurfaceMesh mesh(points, triangles);
    require(!AutoRemesher::SurfaceAnalysis(mesh, .5, 90, 0, 0).supportsRimConstraints(),
        "a shallow patch acquired hard rim constraints");
    for (double scale : { 1., 1e-9 }) {
        auto shell = points;
        shell.back()[2] = -1.;
        for (auto& p : shell)
            p *= scale;
        const AutoRemesher::SurfaceMesh deep(shell, triangles);
        require(AutoRemesher::SurfaceAnalysis(deep, .5 * scale, 90, 0, 0).supportsRimConstraints(),
            "a resolved planar rim lost its constraints");
        require(!AutoRemesher::SurfaceAnalysis(deep, 4. * scale, 90, 0, 0).supportsRimConstraints(),
            "a subcell opening acquired hard rim constraints");
        shell[2][2] = scale;
        require(!AutoRemesher::SurfaceAnalysis(AutoRemesher::SurfaceMesh(shell, triangles), .5 * scale, 90, 0, 0).supportsRimConstraints(),
            "a spatial rim acquired planar constraints");
    }
    // Authored boundaries remain integer isolines even when the field misses their tangent.
    for (double angle : { M_PI / 9, M_PI / 4 }) {
        const std::vector<Vec3> field(triangles.size(), Vec3(std::cos(angle), std::sin(angle), 0));
        for (bool featureLayout : { false, true })
            for (bool preserveBoundary : { false, true }) {
                AutoRemesher::QuadParameterizer::Result result;
                require(AutoRemesher::QuadParameterizer::parameterize(points, triangles, &field,
                            .5 / mesh.averageEdgeLength(), 90, &result, nullptr, nullptr, nullptr, nullptr, nullptr, featureLayout, nullptr, preserveBoundary),
                    "rotated boundary cover failed");
                double area = 0;
                for (const auto& t : result.triangleUvs)
                    area += std::fabs((t[1].x() - t[0].x()) * (t[2].y() - t[0].y()) - (t[1].y() - t[0].y()) * (t[2].x() - t[0].x()));
                require(area > 1e-6, "boundary constraints collapsed the cover");
                if (!preserveBoundary)
                    require(std::fabs(area - 64.) < 1e-5, "unconstrained cover changed its affine area");
                for (size_t c = 0; c < mesh.cornerCount(); ++c)
                    if (preserveBoundary && mesh.isBoundaryCorner(c)) {
                        const auto& a = result.triangleUvs[c / 3][c % 3];
                        const auto& b = result.triangleUvs[c / 3][(c + 1) % 3];
                        bool isoline = false;
                        for (size_t axis = 0; axis < 2; ++axis)
                            isoline |= std::fabs(a[axis] - b[axis]) < 1e-6 && std::fabs(a[axis] - std::round(a[axis])) < 1e-6
                                && std::fabs(a[1 - axis] - b[1 - axis]) > 1e-6;
                        require(isoline, "authored boundary missed a noncollapsed integer isoline");
                    }
            }
    }
}

static void directionalQuadCover()
{
    // A very uneven triangulation must reproduce an affine rectangular grid,
    // including anisotropic spacing and a change of world units.
    const Faces triangles = { { 0, 1, 4 }, { 1, 2, 4 }, { 2, 3, 4 }, { 3, 0, 4 } };
    for (double worldScale : { 1., 7. })
        for (double uSpacing : { 1., 2. })
            for (bool adaptiveSizing : { false, true }) {
                std::vector<Vec3> positions = { { 0, 0, 0 }, { 4, 0, 0 }, { 4, 2, 0 }, { 0, 2, 0 }, { .07, .73, 0 } };
                for (auto& p : positions)
                    p = p * worldScale;
                const AutoRemesher::SurfaceMesh mesh(positions, triangles);
                const std::vector<Vec3> field(4, Vec3(1, 0, 0));
                const std::vector<double> u(4, uSpacing), v(4, 1);
                // Local detail may resolve below the component's nominal spacing floor.
                // Constant physical targets must survive the final transported sizing
                // stage, including unequal triangle areas and a change of world units.
                AutoRemesher::SurfaceGuidance sizing;
                sizing.faces.resize(4);
                sizing.spacingUpper = 10;
                sizing.spacingLower = 2.5;
                sizing.spacingAspect = 2.3;
                for (auto& face : sizing.faces)
                    face.adaptiveWeight = true;
                const std::vector<double> scalar(4, 1);
                AutoRemesher::QuadParameterizer::Result result;
                if (worldScale == 1 && uSpacing == 1 && !adaptiveSizing) {
                    const std::vector<char> shortCorners;
                    const std::vector<double> shortScale;
                    AutoRemesher::SurfaceGuidance shortGuidance;
                    std::vector<Vec3> outputField;
                    require(!AutoRemesher::QuadParameterizer::parameterize(positions, triangles, &field, 1, 90, &result,
                                nullptr, nullptr, nullptr, nullptr, &shortCorners, true),
                        "short feature array accepted");
                    require(!AutoRemesher::QuadParameterizer::parameterize(positions, triangles, &field, 1, 90, &result,
                                &shortScale, nullptr, nullptr, nullptr, nullptr, true),
                        "short sizing array accepted");
                    require(!AutoRemesher::FrameField::create(mesh, 90, &outputField, &shortGuidance), "short field guidance accepted");
                }
                require(AutoRemesher::QuadParameterizer::parameterize(positions, triangles, &field,
                            worldScale / mesh.averageEdgeLength(), 90, &result, adaptiveSizing ? &scalar : nullptr, &u, &v,
                            nullptr, nullptr, adaptiveSizing, adaptiveSizing ? &sizing : nullptr),
                    "directional quad cover failed");
                for (size_t f = 0; f < triangles.size(); ++f)
                    for (size_t k = 0; k < 3; ++k) {
                        const auto& p = positions[triangles[f][k]];
                        const auto& uv = result.triangleUvs[f][k];
                        require(std::fabs(uv.x() - p.x() / (worldScale * uSpacing)) < 1e-5 && std::fabs(uv.y() - p.y() / worldScale) < 1e-5, adaptiveSizing ? "final spacing changed a resolved affine grid" : "quad cover distorted an affine grid");
                    }
                AutoRemesher::SurfaceAnalysis analysis(mesh, worldScale, 90, 0, 0);
                AutoRemesher::QuadExtractor extractor(&positions, &triangles, &result.triangleUvs);
                extractor.setSurfaceAnalysis(&analysis);
                require(extractor.extract(), "affine rectangle extraction failed");
                const auto& points = extractor.remeshedVertices();
                const auto& faces = extractor.remeshedQuads();
                require(faces.size() == size_t(8 / uSpacing), "affine rectangle cell count changed");
                std::set<size_t> used;
                double area = 0;
                for (const auto& face : faces) {
                    require(face.size() == 4, "affine grid lost a quad");
                    used.insert(face.begin(), face.end());
                    for (size_t k = 1; k + 1 < face.size(); ++k)
                        area += .5 * Vec3::crossProduct(points[face[k]] - points[face[0]], points[face[k + 1]] - points[face[0]]).length();
                }
                for (size_t c = 0; c < 4; ++c) {
                    bool found = false;
                    for (size_t v : used)
                        found |= (points[v] - positions[c]).length() < 1e-9 * worldScale;
                    require(found, "graph simplification erased an authored boundary corner");
                }
                require(std::fabs(area - 8 * worldScale * worldScale) < 1e-8 * worldScale * worldScale, "affine rectangle lost area");
            }
}

static void staggeredFeatureBand()
{
    // Two staggered, nearly parallel feature segments on an anisotropic patch.
    // Their closest tips have the opposite transverse order from the continuous
    // cover. Opening the band must preserve that cover order and one cell of width.
    for (double scale : { .25, 1., 7. })
        for (bool rotate : { false, true }) {
            std::vector<Vec3> p;
            Faces t;
            for (size_t y = 0; y < 4; ++y)
                for (size_t x = 0; x < 3; ++x) {
                    double px = double(x) - 1, py = y < 2 ? double(y) - 1 : double(y) - 2 + .1;
                    if (x == 1)
                        px = y < 2 ? .1 * (1 - double(y)) : .002 - .1 * (double(y) - 2);
                    p.push_back((rotate ? Vec3(0, px, py) : Vec3(px, py, 0)) * scale);
                }
            for (size_t y = 0; y < 3; ++y)
                for (size_t x = 0; x < 2; ++x) {
                    const size_t a = 3 * y + x;
                    t.push_back({ a, a + 1, a + 4 });
                    t.push_back({ a, a + 4, a + 3 });
                }
            const AutoRemesher::SurfaceMesh mesh(p, t);
            std::vector<char> features(mesh.cornerCount(), 0);
            for (size_t c = 0; c < mesh.cornerCount(); ++c) {
                const auto a = mesh.cornerVertex(c), b = mesh.cornerVertex(mesh.nextCorner(c));
                features[c] = (std::min(a, b) == 1 && std::max(a, b) == 4)
                    || (std::min(a, b) == 7 && std::max(a, b) == 10);
            }
            const std::vector<Vec3> field(t.size(), rotate ? Vec3(0, 0, 1) : Vec3(0, 1, 0));
            const std::vector<double> u(t.size(), 1), v(t.size(), .05);
            AutoRemesher::QuadParameterizer::Result result;
            require(AutoRemesher::QuadParameterizer::parameterize(p, t, &field, .4 * scale / mesh.averageEdgeLength(),
                        180, &result, nullptr, &u, &v, nullptr, &features, true, nullptr, true, true),
                "staggered feature cover failed");
            require(std::fabs(result.triangleUvs[8][1].y() - result.triangleUvs[0][2].y() - 1) < 1e-7,
                "feature separation reversed or collapsed the staggered band");
            AutoRemesher::QuadParameterizer::Result legacy;
            require(AutoRemesher::QuadParameterizer::parameterize(p, t, &field, .4 * scale / mesh.averageEdgeLength(),
                        180, &legacy, nullptr, &u, &v, nullptr, &features, true, nullptr, true),
                "legacy staggered feature cover failed");
            require(std::fabs(legacy.triangleUvs[8][1].y() - legacy.triangleUvs[0][2].y() + 1) < 1e-7,
                "default feature separation changed legacy ordering");
            for (const auto& uv : result.triangleUvs)
                require((uv[1].x() - uv[0].x()) * (uv[2].y() - uv[0].y())
                        - (uv[1].y() - uv[0].y()) * (uv[2].x() - uv[0].x()) > 0,
                    "staggered feature cover folded a triangle");
            for (size_t c = 0; c < features.size(); ++c)
                if (features[c]) {
                    const auto& a = result.triangleUvs[c / 3][c % 3];
                    const auto& b = result.triangleUvs[c / 3][(c + 1) % 3];
                    require(std::fabs(a.y() - b.y()) < 1e-7 && std::fabs(a.y() - std::round(a.y())) < 1e-7
                            && std::fabs(a.x() - b.x()) > 1,
                        "staggered feature lost its noncollapsed integer isoline");
                }
        }
}

static void oddTubePeriod()
{
    // Five cells around a tube require an odd chart translation. An even-only
    // lattice changes the circumference and can collapse narrow curved tubes.
    std::vector<Vec3> p;
    Faces t;
    std::vector<std::vector<UV>> uv;
    const size_t sides = 16;
    tube(sides, 7, p, t, uv);
    AutoRemesher::SurfaceMesh mesh(p, t);
    std::vector<Vec3> field(t.size(), Vec3(0, 0, 1));
    std::vector<double> u(t.size(), 1), v(t.size(), 2 * sides * std::sin(M_PI / sides) / 5);
    AutoRemesher::QuadParameterizer::Result result;
    require(AutoRemesher::QuadParameterizer::parameterize(p, t, &field, 1 / mesh.averageEdgeLength(),
                90, &result, nullptr, &u, &v, nullptr, nullptr, true),
        "odd tube cover failed");
    double circumference = 0;
    for (size_t i = 0; i < sides; ++i) {
        const auto d = result.triangleUvs[2 * i][1] - result.triangleUvs[2 * i][0];
        circumference += std::hypot(d.x(), d.y());
    }
    require(std::fabs(circumference - 5) < 1e-5, "tube circumference was forced onto an even lattice");
}

static void featureIntegerLattice()
{
    using System = AutoRemesher::MixedIntegerLeastSquares;
    // Independent near-half-integer targets admit only one variable per pass.
    // The cover's iteration allowance must scale beyond its old 100-pass cap.
    System slow(201);
    for (size_t i = 0; i < 201; ++i) {
        slow.setVariablePeriod(i, 1);
        slow.addEnergy({ { i, 1 } }, .49);
    }
    slow.finalizeConstraints();
    const size_t limit = slow.integerKernelVariableCount() + 1;
    for (size_t i = 0; i < limit && !slow.converged(); ++i)
        require(slow.solveIteration(true, true), "bounded integer rounding failed");
    require(slow.converged(), "integer rounding stopped before admitting every variable");
    for (bool star : { false, true }) {
        System system(6);
        for (size_t i = 0; i < 3; ++i) {
            system.setVariablePeriod(i, 1);
            system.setVariablePeriod(i + 3, 2);
            system.addConstraint(i, 1, i + 3, -1);
            system.addEnergy({ { i, 1 } }, .1 * i);
        }
        system.finalizeConstraints();
        require(system.solveIteration(false), "continuous lattice solve");
        require(system.separateIntegerCoordinates(0, 1, 1), "first feature separation rejected");
        require(system.separateIntegerCoordinates(star ? 0 : 1, 2, 1), "adjacent feature separation rejected after elimination");
        for (size_t i = 0; i < 100 && !system.converged(); ++i)
            require(system.solveIteration(true, true), "budgeted lattice solve failed");
        require(system.converged(), "budgeted rounding failed to converge");
        require(std::fabs(system.value(1) - system.value(0) - 2) < 1e-8 && std::fabs(system.value(2) - system.value(star ? 0 : 1) - 2) < 1e-8, "positive feature lattice spacing lost");
        for (size_t i = 0; i < 6; ++i)
            require(std::fabs(system.value(i) / 2 - std::round(system.value(i) / 2)) < 1e-8, "period-two cover lost integrality");
        require(!system.separateIntegerCoordinates(99, 0, 1), "out of range feature accepted");
    }
    // A sparse solve can land on either side of a half-integer by roundoff.
    // Original stored-float decisions must use the same tie for scoring and fixing.
    for (double target : { 1. - 1e-12, 1. + 1e-12 }) {
        System tie(1);
        tie.setVariablePeriod(0, 2);
        tie.addEnergy({ { 0, 1 } }, target);
        tie.finalizeConstraints();
        require(tie.solveIteration(false) && tie.solveIteration(true, true) && tie.converged(), "half-integer solve failed");
        require(std::fabs(tie.value(0) - 2) < 1e-8, "half-integer decision changed with solver roundoff");
    }
}

static void tubeOpenings()
{
    for (size_t sides : { 4, 12 }) {
        std::vector<Vec3> p;
        Faces t;
        std::vector<std::vector<UV>> uv;
        tube(sides, sides == 4 ? 3 : 5, p, t, uv);
        AutoRemesher::SurfaceAnalysis analysis(AutoRemesher::SurfaceMesh(p, t), .5, 90, 0, 0);
        AutoRemesher::QuadExtractor extractor(&p, &t, &uv);
        if (sides == 12)
            extractor.setSurfaceAnalysis(&analysis);
        require(extractor.extract(), "tube opening extraction failed");
        const auto edges = edgeUses(extractor.remeshedQuads());
        require(!edges.empty(), "hole closure discarded the tube");
        if (sides == 4) {
            // Without source boundaries, fill each small hole exactly once.
            for (const auto& edge : edges)
                require(edge.second == 2, "four-edge hole was filled more than once");
        } else {
            // Authored source rims must remain open.
            require(analysis.onSourceBoundary(p[0]) && !analysis.onSourceBoundary(p[2 * sides]), "source rim confused with interior");
            size_t lower = 0, upper = 0;
            for (const auto& edge : edges)
                if (edge.second == 1) {
                    const auto& a = extractor.remeshedVertices()[edge.first.first];
                    const auto& b = extractor.remeshedVertices()[edge.first.second];
                    lower += std::fabs(a.z()) < 1e-8 && std::fabs(b.z()) < 1e-8;
                    upper += std::fabs(a.z() - 4) < 1e-8 && std::fabs(b.z() - 4) < 1e-8;
                }
            require(lower == sides && upper == sides, "authored source opening was capped");
        }
    }
}

static void radialCover()
{
    std::vector<Vec3> p { { 0, 0, 0 } };
    Faces t;
    const size_t sides = 32;
    for (double r : { .35, .4, .45, 1. })
        for (size_t i = 0; i < sides; ++i)
            p.emplace_back(r * std::cos(2 * M_PI * i / sides), r * std::sin(2 * M_PI * i / sides), 0);
    for (size_t i = 0; i < sides; ++i)
        t.push_back({ 0, 1 + i, 1 + (i + 1) % sides });
    for (size_t j = 0; j < 3; ++j)
        for (size_t i = 0; i < sides; ++i) {
            const size_t a = 1 + j * sides + i, b = 1 + j * sides + (i + 1) % sides, c = b + sides, d = a + sides;
            t.push_back({ a, d, c });
            t.push_back({ a, c, b });
        }
    std::vector<Vec3> field;
    std::vector<char> features;
    for (const auto& f : t) {
        field.push_back((p[f[0]] + p[f[1]] + p[f[2]]).normalized());
        for (size_t k = 0; k < 3; ++k)
            features.push_back(std::fabs(p[f[k]].length() - p[f[(k + 1) % 3]].length()) < 1e-8);
    }
    AutoRemesher::SurfaceMesh mesh(p, t);
    AutoRemesher::QuadParameterizer::Result result;
    require(AutoRemesher::QuadParameterizer::parameterize(p, t, &field, .25 / mesh.averageEdgeLength(), 90, &result,
                nullptr, nullptr, nullptr, nullptr, &features, true),
        "radial cover failed");
    require(result.fullTurnVertices == std::vector<size_t> { 0 }, "radial full-turn center was lost modulo four");
    double area = 0;
    for (const auto& q : result.triangleUvs) {
        const double a = (q[1].x() - q[0].x()) * (q[2].y() - q[0].y()) - (q[1].y() - q[0].y()) * (q[2].x() - q[0].x());
        require(a >= -1e-6, "radial cover folded");
        area += a;
    }
    require(area > 10, "radial grid collapsed to a line");

    // Join two radial disks into a thin closed shell. Pole ownership must follow
    // source connectivity even when both poles satisfy a loop's plane tolerance.
    const size_t offset = p.size(), diskFaces = t.size();
    for (size_t v = 0; v < offset; ++v)
        p.push_back(p[v] - Vec3(0, 0, .002));
    for (size_t f = 0; f < diskFaces; ++f) {
        t.push_back({ t[f][2] + offset, t[f][1] + offset, t[f][0] + offset });
        field.push_back(field[f]);
        for (size_t k : { size_t(1), size_t(0), size_t(2) })
            features.push_back(features[3 * f + k]);
    }
    for (size_t i = 0; i < sides; ++i) {
        const size_t a = 1 + 3 * sides + i, b = 1 + 3 * sides + (i + 1) % sides;
        t.push_back({ a, a + offset, b + offset });
        t.push_back({ a, b + offset, b });
        field.insert(field.end(), 2, Vec3(0, 0, 1));
        features.insert(features.end(), 6, 0);
    }
    AutoRemesher::SurfaceMesh shell(p, t);
    require(AutoRemesher::QuadParameterizer::parameterize(p, t, &field, .25 / shell.averageEdgeLength(), 90, &result,
                nullptr, nullptr, nullptr, nullptr, &features, true),
        "closed radial cover failed");
    require(result.fullTurnVertices == std::vector<size_t>({ 0, offset }), "closed shell lost a radial pole");
    AutoRemesher::SurfaceAnalysis analysis(shell, .25, 90, 0, 0, true, false);
    Faces firstFaces;
    std::vector<Vec3> firstPoints;
    for (size_t order = 0; order < 2; ++order) {
        AutoRemesher::QuadExtractor q(&p, &t, &result.triangleUvs);
        q.setSurfaceAnalysis(&analysis);
        q.setFullTurnVertices(&result.fullTurnVertices);
        require(q.extract(), "closed radial extraction failed");
        size_t fans[2] = { 0, 0 };
        for (const auto& face : q.remeshedQuads())
            if (face.size() == 3)
                for (size_t v : face) {
                    const auto& center = q.remeshedVertices()[v];
                    if (center.x() != 0 || center.y() != 0)
                        continue;
                    ++fans[center.z() < -.001];
                    for (size_t corner : face)
                        require(std::fabs(q.remeshedVertices()[corner].z() - center.z()) < 1e-10,
                            "pole fan crossed to another source sheet");
                }
        require(fans[0] > 0 && fans[1] > 0, "pole ownership discarded a radial fan");
        if (order) {
            require(q.remeshedQuads() == firstFaces && q.remeshedVertices().size() == firstPoints.size(), "pole order changed topology");
            for (size_t v = 0; v < firstPoints.size(); ++v)
                require((q.remeshedVertices()[v] - firstPoints[v]).lengthSquared() == 0, "pole order changed geometry");
        } else {
            firstFaces = q.remeshedQuads();
            firstPoints = q.remeshedVertices();
        }
        std::reverse(result.fullTurnVertices.begin(), result.fullTurnVertices.end());
    }
}

static void boundedRefinement()
{
    std::vector<::Vector3> points { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 } };
    Faces triangles { { 0, 1, 2 }, { 0, 2, 3 } };
    for (size_t limit : { size_t(3), size_t(80) }) {
        ::IsotropicRemesher remesher(&points, &triangles);
        remesher.setTargetEdgeLength(.001);
        remesher.setRefinementVertexLimit(limit);
        remesher.remesh(0);
        auto* mesh = remesher.remeshedHalfedgeMesh();
        size_t count = 0;
        double area = 0;
        for (auto* v = mesh->moveToNextVertex(nullptr); v; v = mesh->moveToNextVertex(v))
            ++count;
        for (auto* f = mesh->moveToNextFace(nullptr); f; f = mesh->moveToNextFace(f)) {
            auto* e = f->halfedge;
            area += ::Vector3::area(e->startVertex->position, e->nextHalfedge->startVertex->position, e->previousHalfedge->startVertex->position);
        }
        require(count == std::max(points.size(), limit), "refinement exceeded its vertex budget or deleted original geometry");
        require(std::fabs(area - 1) < 1e-12, "budget stop lost source coverage");
    }
}

static void originalOpenWorkflow()
{
    // A weakly creased cube is closed; removing its warped top is preparation only.
    // Reversing that relationship checks that capping does not erase provenance.
    for (double scale : { 1., 7. }) {
        std::vector<Vec3> points = { { -1, -1, -1 }, { 1, -1, -1 }, { 1, 1, -1 }, { -1, 1, -1 },
            { -1, -1, 1.4 }, { 1, -1, 1 }, { 1, 1, 1 }, { -1, 1, 1 } };
        for (auto& p : points)
            p = scale * Vec3(p.z(), p.x(), p.y());
        const Faces quads = { { 0, 3, 2, 1 }, { 4, 5, 6, 7 }, { 0, 1, 5, 4 },
            { 1, 2, 6, 5 }, { 2, 3, 7, 6 }, { 3, 0, 4, 7 } };
        Faces closed, open;
        for (size_t f = 0; f < quads.size(); ++f)
            for (size_t k = 1; k < 3; ++k) {
                closed.push_back({ quads[f][0], quads[f][k], quads[f][k + 1] });
                if (f != 1)
                    open.push_back(closed.back());
            }
        const auto hasCurvature = [](const AutoRemesher::SurfaceAnalysis& analysis) {
            return std::any_of(analysis.faces().begin(), analysis.faces().end(),
                [](const AutoRemesher::SurfaceGuidance::Face& face) { return face.major > 1e-6; });
        };
        Faces inconsistent = closed;
        std::reverse(inconsistent[0].begin(), inconsistent[0].end());
        AutoRemesher::SurfaceAnalysis windingGap(AutoRemesher::SurfaceMesh(points, inconsistent), .2 * scale, 100, 1, 1, true, true, false, nullptr, true);
        require(!hasCurvature(windingGap), "two-sided incidence was mistaken for an authored opening");
        for (bool originallyOpen : { false, true }) {
            const Faces& source = originallyOpen ? open : closed;
            const Faces& prepared = originallyOpen ? closed : open;
            AutoRemesher::SurfaceAnalysis original(AutoRemesher::SurfaceMesh(points, source), .2 * scale, 100, 1, 1, true, true, false, nullptr, true);
            AutoRemesher::SurfaceAnalysis measured(AutoRemesher::SurfaceMesh(points, prepared), .2 * scale, 100, 1, 1,
                true, true, false, &original);
            require(hasCurvature(original) == originallyOpen, "closed weak creases changed their legacy curvature exclusion");
            require(hasCurvature(measured) == originallyOpen, "prepared hole or cap changed the original curvature policy");
            AutoRemesher::Parameterizer parameterizer(&points, &prepared, nullptr);
            parameterizer.setSurfaceAnalysis(&original);
            parameterizer.setSharpEdgeDegrees(100);
            parameterizer.setGradientAdaptivity(1);
            parameterizer.setAnisotropy(1);
            require(parameterizer.parameterize(true), "source-provenance cover failed");
            require(!parameterizer.relaxationMetric().empty() == originallyOpen,
                "prepared hole or cap changed the original relaxation metric policy");
        }
    }
}

static void curvatureFieldSupport()
{
    // A marginal curvature estimate should yield to the adjacent authored axis.
    // Well-supported curvature, explicit features and legacy solves stay fixed.
    for (double scale : { 1., 7. })
        for (bool taper : { true, false })
            for (bool feature : { false, true })
                for (double support : { .28, .8 }) {
                    const std::vector<Vec3> points { { 0, 0, 0 }, { scale, 0, 0 },
                        { scale, scale, 0 }, { 0, scale, 0 } };
                    const Faces faces { { 0, 1, 2 }, { 0, 2, 3 } };
                    AutoRemesher::SurfaceMesh mesh(points, faces);
                    AutoRemesher::SurfaceGuidance guidance;
                    guidance.faces.resize(2);
                    guidance.featureCorners.assign(6, 0);
                    guidance.featureCorners[0] = 1;
                    guidance.featureCorners[4] = feature;
                    const double angle = M_PI / 8;
                    guidance.faces[0].direction = Vec3(1, 0, 0);
                    guidance.faces[1].direction = Vec3(std::cos(angle), std::sin(angle), 0);
                    for (auto& face : guidance.faces) {
                        face.confidence = 1;
                        face.major = support / scale;
                        face.radius = scale;
                    }
                    std::vector<Vec3> field;
                    require(AutoRemesher::FrameField::create(mesh, 90, &field, &guidance, true, taper),
                        "curvature support field solve failed");
                    const auto difference = [](const Vec3& a, const Vec3& b) {
                        const double cosine = std::min(1., std::fabs(Vec3::dotProduct(a, b)));
                        const double radians = std::acos(cosine);
                        return std::min(radians, M_PI / 2 - radians);
                    };
                    require(difference(field[0], guidance.faces[0].direction) < 1e-7,
                        "curvature support changed an authored field axis");
                    if (taper && !feature && support < .3)
                        require(difference(field[1], field[0]) < .01,
                            "marginal curvature overrode the adjacent authored axis");
                    else
                        require(difference(field[1], guidance.faces[1].direction) < 1e-7,
                            "curvature support weakened a reliable, explicit or legacy constraint");
                }
}

int main(int argc, char** argv)
{
    try {
        if (argc == 1) {
            collapsedOutputFaces();
            safeHexagonSplits();
            curveConstrainedPentagonMerge();
            curvatureFieldSupport();
            originalOpenWorkflow();
            occupiedFaceInteriors();
            constantCoordinateIntersections();
            subdividedExtractionCorners();
            protectedPreparationRims();
            boundedRefinement();
            balancedRefinement();
            fixtures();
            tinyIsland();
            thinTubeCleanup();
            telescopingTube();
            disconnectedFans();
            directionalQuadCover();
            rotatedBoundaryCover();
            staggeredFeatureBand();
            oddTubePeriod();
            featureIntegerLattice();
            tubeOpenings();
            radialCover();
        } else {
            require(argc == 4 || (argc == 5 && std::string(argv[4]) == "--cloth"),
                "usage: reference_surface_test [input.obj target_quads output.obj [--cloth]]");
            tinyobj::attrib_t attributes;
            std::vector<tinyobj::shape_t> shapes;
            std::vector<tinyobj::material_t> materials;
            std::string warning, error;
            require(tinyobj::LoadObj(&attributes, &shapes, &materials, &warning, &error, argv[1]), "OBJ load failed");
            std::vector<Vec3> vertices;
            for (size_t i = 0; i < attributes.vertices.size(); i += 3)
                vertices.push_back({ attributes.vertices[i], attributes.vertices[i + 1], attributes.vertices[i + 2] });
            Faces triangles;
            for (const auto& shape : shapes)
                for (size_t i = 0; i < shape.mesh.indices.size(); i += 3)
                    triangles.push_back({ size_t(shape.mesh.indices[i].vertex_index),
                        size_t(shape.mesh.indices[i + 1].vertex_index), size_t(shape.mesh.indices[i + 2].vertex_index) });
            Remesher remesher(vertices, triangles);
            remesher.setClothFoldGuidance(argc == 5);
            remesher.setTargetTriangleCount(std::stoull(argv[2]) * 2);
            remesher.setScaling(1.0);
            remesher.setSharpEdgeDegrees(90.0);
            require(remesher.remesh(), "remesh failed");
            checkReferences(remesher, vertices, triangles);
            require(!remesher.remeshedQuads().empty(), "empty result");
            for (const auto& p : remesher.remeshedVertices())
                for (size_t axis = 0; axis < 3; ++axis)
                    require(std::isfinite(p[axis]), "non-finite output vertex");
            for (const auto& face : remesher.remeshedQuads()) {
                require(face.size() >= 3, "short output face");
                for (size_t v : face)
                    require(v < remesher.remeshedVertices().size(), "output face index out of range");
            }
            std::ofstream out(argv[3]);
            out << std::setprecision(17);
            for (const auto& p : remesher.remeshedVertices())
                out << "v " << p.x() << ' ' << p.y() << ' ' << p.z() << '\n';
            for (const auto& f : remesher.remeshedQuads()) {
                out << 'f';
                for (size_t v : f)
                    out << ' ' << v + 1;
                out << '\n';
            }
            require(out.good(), "OBJ write failed");
            size_t chains = 0, selected = 0;
            for (const auto& island : remesher.preparedIslands())
                for (const auto& chain : island.analysis->chains()) {
                    ++chains;
                    selected += chain.strength > 0;
                }
            std::cout << "Feature chains: " << selected << " selected of " << chains << '\n';
            std::cout << "Reference islands: " << remesher.preparedIslands().size()
                      << ", decimation used: " << remesher.decimated() << '\n';
        }
        std::cout << "Reference checks passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
