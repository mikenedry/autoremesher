#include <AutoRemesher/SurfaceAnalysis>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <type_traits>
static_assert(!std::is_move_constructible<AutoRemesher::SurfaceAnalysis>::value, "analysis tree storage must stay put");
using namespace AutoRemesher;
using V = AutoRemesher::Vector3;
static void require(bool ok, const char* message)
{
    if (!ok)
        throw std::runtime_error(message);
}

static void checkMetric(const SurfaceAnalysis& analysis)
{
    for (const auto& metric : analysis.faces()) {
        require(std::isfinite(metric.scale) && metric.scale > 0 && std::isfinite(metric.ratio), "nonfinite sizing");
        require(std::isfinite(metric.confidence) && metric.confidence >= 0 && metric.confidence <= 1 + 1e-8, "invalid curvature confidence");
        require(metric.ratio >= 1 / 2.3 - 1e-8 && metric.ratio <= 2.3 + 1e-8, "aspect ratio exceeded limit");
        require(metric.scale * std::sqrt(metric.ratio) <= 1 + 1e-12 && metric.scale / std::sqrt(metric.ratio) <= 1 + 1e-12, "curvature enlarged the base spacing");
    }
}

static void cylinder(bool featureLayout = false)
{
    std::vector<V> p;
    std::vector<std::vector<size_t>> t;
    const size_t sides = 48, rings = 13;
    for (size_t z = 0; z < rings; ++z)
        for (size_t x = 0; x < sides; ++x) {
            const double a = 2 * M_PI * x / sides;
            p.emplace_back(std::cos(a), std::sin(a), double(z) / 6);
        }
    for (size_t z = 0; z + 1 < rings; ++z)
        for (size_t x = 0; x < sides; ++x) {
            size_t a = z * sides + x, b = z * sides + (x + 1) % sides, c = a + sides, d = b + sides;
            t.push_back({ a, b, d });
            t.push_back({ a, d, c });
        }
    SurfaceMesh mesh(p, t);
    SurfaceAnalysis analysis(mesh, .25, 90, 1, 1, featureLayout);
    checkMetric(analysis);
    size_t reliable = 0;
    for (size_t i = 0; i < t.size(); ++i) {
        const auto& f = analysis.faces()[i];
        const auto& v = p[t[i][0]];
        const V tangent = V(-v.y(), v.x(), 0).normalized();
        if (f.major > .5 && f.major < 2 && f.minor < .1 && std::fabs(V::dotProduct(f.direction, tangent)) > .95)
            ++reliable;
    }
    require(reliable > t.size() / 2, "cylinder curvature/directional sizing is incorrect");
    SurfaceAnalysis coarse(mesh, 1, 90, 1, 1, featureLayout);
    size_t directional = 0;
    for (const auto& f : coarse.faces())
        directional += f.ratio < .9 && f.scale < .8;
    require(directional > t.size() / 2, "physical curvature did not contract the circumference direction");
    // A uniformly thin island must gain resolution even at a fixed global budget.
    auto thin = p;
    for (auto& v : thin)
        v = v * .01;
    SurfaceAnalysis wire(SurfaceMesh(thin, t), .25, 90, 1, 1, featureLayout);
    size_t resolved = 0;
    for (const auto& f : wire.faces())
        resolved += f.scale < .15 && .25 * f.scale * std::sqrt(f.ratio) < .025;
    require(wire.scalarSize(thin[0]) >= .25 * .25, "triangle sampling exceeded its density budget");
    require(resolved > t.size() / 2, "thin island refinement was normalized away or stopped at the ordinary floor");
    // Physical scale changes curvature inversely, and leave relative sizes alone.
    for (auto& v : p)
        v = v * 7;
    SurfaceAnalysis scaled(SurfaceMesh(p, t), 1.75, 90, 1, 1, featureLayout);
    for (size_t f = 0; f < t.size(); ++f) {
        require(std::fabs(analysis.faces()[f].major - 7 * scaled.faces()[f].major) < 1e-7, "curvature is not scale covariant");
        require(std::fabs(analysis.faces()[f].scale - scaled.faces()[f].scale) < 1e-7, "sizes changed under uniform scaling");
    }
    // A coarse rim chord belongs to its source boundary even when its midpoint
    // is inside the circle; endpoints on distinct rims must not be conflated.
    require(scaled.onSourceBoundary(p[0], p[2]), "coarse rim chord lost boundary ownership");
    require(!scaled.onSourceBoundary(p[0], p.back()), "separate source rims were conflated");
    // Every round boundary is one closed chain, not many edge fragments.
    require(analysis.chains().size() == 2 && analysis.chains()[0].closed && analysis.chains()[1].closed,
        "closed boundary chains were not recognized");
}

static void creaseAndTransfer()
{
    // 45-degree crease: automatic at a 90-degree hard threshold.
    std::vector<V> p;
    std::vector<std::vector<size_t>> t;
    for (size_t i = 0; i < 9; ++i) {
        const double x = double(i) / 4;
        p.emplace_back(x, -1, 0);
        p.emplace_back(x, 0, 0);
        p.emplace_back(x, 1, 1);
    }
    for (size_t i = 0; i < 8; ++i)
        for (size_t j = 0; j < 2; ++j) {
            const size_t a = 3 * i + j, b = a + 3;
            t.push_back({ a, b, b + 1 });
            t.push_back({ a, b + 1, a + 1 });
        }
    SurfaceMesh mesh(p, t);
    SurfaceAnalysis analysis(mesh, .2, 90, 1, 1);
    for (const auto& face : analysis.faces())
        require(face.major < 1e-12 && face.minor < 1e-12 && face.scale == 1, "a crease injected curvature into its planar patches");
    // A weak alignment guide is not a sharp-crease exclusion for vertex tensors.
    SurfaceAnalysis aligned(mesh, .2, 90, 1, 1, true, true, false, nullptr, true);
    require(std::count_if(aligned.faces().begin(), aligned.faces().end(),
                [](const SurfaceGuidance::Face& face) { return face.major > .1; }) > t.size() / 2,
        "weak alignment guide erased the fold's curvature samples");
    const auto alignedGuidance = aligned.transfer(mesh, true);
    size_t alignedCorners = 0;
    for (size_t c = 0; c < mesh.cornerCount(); ++c)
        alignedCorners += !mesh.isBoundaryCorner(c) && alignedGuidance.featureCorners[c];
    require(alignedCorners == 16,
        "curvature sampling removed the alignment constraints");
    auto narrow = p;
    for (auto& point : narrow)
        point = V(point.x(), .01 * point.y(), .01 * point.z());
    SurfaceAnalysis band(SurfaceMesh(narrow, t), .2, 90, 1, 1);
    require(std::any_of(band.faces().begin(), band.faces().end(), [](const SurfaceGuidance::Face& face) { return face.major > 1 && face.scale < 1; }), "narrow planar bands lost their curvature density safeguard");
    require(analysis.onSourceBoundary(V(.005, 0, 0)), "interior crease hid a nearby authored boundary");
    require(!analysis.onSourceBoundary(V(1, 0, 0)), "interior crease was mistaken for a boundary");
    const auto exact = analysis.transfer(mesh, true);
    std::vector<char> expected(mesh.cornerCount(), 0);
    for (const auto& chain : analysis.chains())
        if (chain.strength > 0)
            for (size_t c : chain.corners) {
                expected[c] = 1;
                if (mesh.oppositeCorner(c) != SurfaceMesh::npos)
                    expected[mesh.oppositeCorner(c)] = 1;
            }
    require(exact.featureCorners == expected, "identity transfer changed the feature corner IDs");
    checkMetric(analysis);
    bool supported = false;
    for (const auto& c : analysis.chains())
        if (c.strength > .5 && c.strength < 1)
            supported = true;
    require(supported, "orthogonal endpoint context did not strengthen the crease");
    const auto binding = analysis.bindCurve(V(1, .01, 0), .05);
    require(binding.chain != SurfaceMesh::npos, "source crease binding failed");
    const auto slide = analysis.projectCurve(binding, V(1, -.99, 0));
    require(std::fabs(slide.y()) < 1e-12 && std::fabs(slide.z()) < 1e-12, "bound vertex hopped to another curve");
    auto working = p;
    for (size_t i = 0; i < 9; ++i)
        working[3 * i + 1] = working[3 * i + 1] + V(0, .01, 0);
    const auto before = working;
    require(analysis.finishCurves(working, t) > 0, "final curve constraints were not exercised");
    for (size_t i = 0; i < 9; ++i)
        require(std::fabs(working[3 * i + 1].y()) < 1e-10 && std::fabs(working[3 * i + 1].z()) < 1e-10, "finishing drifted away from the original crease");
    require((working[1] - p[1]).length() < 1e-10 && (working[25] - p[25]).length() < 1e-10, "junction endpoints moved");
    for (const auto& f : t) {
        const V a = V::crossProduct(before[f[1]] - before[f[0]], before[f[2]] - before[f[0]]);
        const V b = V::crossProduct(working[f[1]] - working[f[0]], working[f[2]] - working[f[0]]);
        require(V::dotProduct(a, b) > 0, "curve finishing flipped an incident triangle");
    }

    SurfaceAnalysis supportedLink(mesh, 1.2, 90, 1, 1);
    bool rescued = false;
    for (const auto& c : supportedLink.chains())
        rescued = rescued || (c.strength > 0 && c.strength < 1);
    require(rescued, "short link supported at both endpoints was not rescued");
    SurfaceAnalysis shortChain(mesh, 3, 90, 1, 1);
    for (const auto& c : shortChain.chains())
        require(c.strength == 0 || c.strength == 1, "short automatic fragment was admitted");
    const auto guidance = analysis.transfer(mesh);
    size_t selected = 0;
    for (size_t c = 0; c < mesh.cornerCount(); ++c)
        if (!mesh.isBoundaryCorner(c) && guidance.featureCorners[c])
            ++selected;
    require(selected == 16, "crease provenance did not transfer to working corners");
    // With both sizing controls off, the shared analysis gives a uniform metric.
    SurfaceAnalysis uniform(mesh, .2, 90, 0, 0);
    for (const auto& f : uniform.faces())
        require(std::fabs(f.scale - 1) < 1e-12 && f.ratio == 1, "disabled sizing still adapted");
    SurfaceAnalysis empty(SurfaceMesh({}, {}), .2, 90, 1, 1);
    require(empty.scalarSize(V()) == .2, "empty reference fallback failed");
}

static void connectedRelaxation()
{
    // Nearby disconnected sheets must not attract a vertex off its original sheet.
    std::vector<V> p;
    std::vector<std::vector<size_t>> t;
    for (double z : { 0., .02 }) {
        const size_t offset = p.size();
        for (size_t x = 0; x < 5; ++x) {
            p.emplace_back(x, -1, z);
            p.emplace_back(x, 1, z);
        }
        for (size_t x = 0; x < 4; ++x) {
            size_t a = offset + 2 * x;
            t.push_back({ a, a + 2, a + 3 });
            t.push_back({ a, a + 3, a + 1 });
        }
    }
    SurfaceAnalysis analysis(SurfaceMesh(p, t), .2, 90, 1, 1);
    std::vector<V> output = { V(.25, 0, 0), V(3.5, 0, .08) };
    std::vector<std::unordered_set<size_t>> neighbors = { { 1 }, { 0 } };
    analysis.relaxSurface(output, neighbors, { false, true }, {}, 8);
    require(std::fabs(output[0].z()) < 1e-12, "relaxation jumped to a nearby sheet");
    require(output[0].x() > 3, "source traversal stalled at a triangle edge");
    require((output[1] - V(3.5, 0, .08)).length() < 1e-12, "relaxation moved a locked vertex");
    // A captured boundary cannot slide because of an off-curve neighbor.
    output = { V(.5, -.99, 0), V(2, 0, 0) };
    analysis.relaxSurface(output, neighbors, { false, true }, {}, 8);
    require(std::fabs(output[0].y() + 1) < 1e-12 && std::fabs(output[0].z()) < 1e-12,
        "relaxation lost its source curve binding");
    require(std::fabs(output[0].x() - .5) < 1e-12, "off-curve neighbor moved a bound vertex");
    output = { V(.5, -.99, 0), V(1, -1, 0), V(3, -1, 0) };
    analysis.relaxSurface(output, { { 1, 2 }, { 0 }, { 0 } }, { false, true, true }, {}, 8);
    require(output[0].x() > 1.9 && std::fabs(output[0].y() + 1) < 1e-12, "same-curve neighbors could not slide a bound vertex");
    // Snapping the triangle apex across its opposite edge would flip this face.
    output = { V(.5, -.995, 0), V(.7, -.999, 0), V(.3, -.999, 0) };
    const auto apex = output[0];
    neighbors = { { 1, 2 }, { 0, 2 }, { 0, 1 } };
    analysis.relaxSurface(output, neighbors, { false, true, true }, { { 0, 1, 2 } }, 8);
    require((output[0] - apex).length() < 1e-12, "curve capture folded an incident face");
    // Nearby interior rows and cross-curve neighbors must not pull a feature row.
    SurfaceAnalysis stencil(SurfaceMesh(p, t), .5, 90, 0, 0, true);
    output = { { .5, -.8, 0 }, { 0, -.8, 0 }, { 1, -.8, 0 } };
    neighbors = { { 1, 2 }, { 0 }, { 0 } };
    stencil.relaxSurface(output, neighbors, { false, true, true }, {}, 2);
    require((output[0] - V(.5, -.8, 0)).length() < 1e-10, "interior vertex captured a nearby curve");
    output = { { .5, -1, 0 }, { .25, -1, 0 }, { .75, -1, 0 }, { 2, 0, 0 } };
    neighbors = { { 1, 2, 3 }, { 0 }, { 0 }, { 0 } };
    stencil.relaxSurface(output, neighbors, { false, true, true, true }, {}, 2);
    require((output[0] - V(.5, -1, 0)).length() < 1e-10, "cross-curve neighbor disturbed a regular feature row");
    // Crossing a planar fan requires visiting faces tied at its center.
    p = { { 0, 0, 0 } };
    t.clear();
    for (size_t i = 0; i < 8; ++i)
        p.emplace_back(std::cos(i * M_PI / 4), std::sin(i * M_PI / 4), 0);
    for (size_t i = 0; i < 8; ++i)
        t.push_back({ 0, 1 + i, 1 + (i + 1) % 8 });
    SurfaceAnalysis fan(SurfaceMesh(p, t), 1, 90, 0, 0, false, false);
    output = { { .2, .05, 0 }, { -.7, -.2, 0 }, { -.7, -.1, 0 } };
    fan.relaxSurface(output, { { 1, 2 }, { 0, 2 }, { 0, 1 } }, { false, true, true }, { { 0, 1, 2 } }, 1);
    require((output[0] - V(-.25, -.05, 0)).length() < 1e-12, "source traversal pinned a vertex at a planar fan center");
}

static void curvedProtectedTube()
{
    std::vector<V> p;
    std::vector<std::vector<size_t>> t;
    for (size_t j = 0; j <= 12; ++j)
        for (size_t i = 0; i < 10; ++i) {
            const double a = j * M_PI / 24, b = i * M_PI / 5, r = 2 + .1 * std::cos(b);
            p.emplace_back(r * std::cos(a), r * std::sin(a), .1 * std::sin(b));
        }
    for (size_t j = 0; j < 12; ++j)
        for (size_t i = 0; i < 10; ++i) {
            const size_t a = 10 * j + i, b = 10 * j + (i + 1) % 10;
            t.push_back({ a, b, b + 10 });
            t.push_back({ a, b + 10, a + 10 });
        }
    SurfaceMesh mesh(p, t);
    SurfaceAnalysis analysis(mesh, .1, 90, 1, 1, true);
    const auto guidance = analysis.transfer(mesh, true);
    size_t protectedCurves = 0;
    for (const auto& chain : analysis.chains())
        if (chain.strength > 0 && !chain.directional) {
            ++protectedCurves;
            const size_t c = chain.corners[chain.corners.size() / 2];
            require(!guidance.featureCorners[c], "a shallow curved tube crease forced a grid axis");
            const auto binding = analysis.bindCurve(mesh.position(mesh.cornerVertex(c)), .001);
            require(binding.chain != SurfaceMesh::npos, "protected tube curve lost its finishing provenance");
        }
    require(protectedCurves > 0, "curved tube retained faceting as grid constraints");
}

static void protectedJunctions()
{
    // Five-way icosahedron junctions cannot all be quad-field directions.
    std::vector<V> p;
    std::vector<std::vector<size_t>> triangles;
    const double phi = (1 + std::sqrt(5.)) / 2;
    for (double a : { -1., 1. })
        for (double b : { -1., 1. }) {
            p.emplace_back(0, a, b * phi);
            p.emplace_back(a, b * phi, 0);
            p.emplace_back(b * phi, 0, a);
        }
    for (size_t a = 0; a < p.size(); ++a)
        for (size_t b = a + 1; b < p.size(); ++b)
            for (size_t c = b + 1; c < p.size(); ++c)
                if (std::fabs((p[a] - p[b]).lengthSquared() - 4) < 1e-10 && std::fabs((p[b] - p[c]).lengthSquared() - 4) < 1e-10 && std::fabs((p[c] - p[a]).lengthSquared() - 4) < 1e-10) {
                    if (V::dotProduct(V::crossProduct(p[b] - p[a], p[c] - p[a]), p[a]) > 0)
                        triangles.push_back({ a, b, c });
                    else
                        triangles.push_back({ a, c, b });
                }
    SurfaceMesh mesh(p, triangles);
    SurfaceAnalysis reference(mesh, .4, 90, 1, 1, true, false);
    require(triangles.size() == 20 && reference.chains().size() == 30, "invalid junction fixture");
    const auto guidance = reference.transfer(mesh, true);
    for (char corner : guidance.featureCorners)
        require(!corner, "five-way protected junction forced a quad axis");
    const V midpoint = (p[triangles[0][0]] + p[triangles[0][1]]) * .5;
    const auto binding = reference.bindCurve(midpoint, .05);
    require(binding.chain != SurfaceMesh::npos, "protected geometry lost its original curve");
    require((reference.projectCurve(binding, midpoint) - midpoint).length() < 1e-12, "protected curve finishing moved the source edge");
}

static void missingSurface()
{
    // Every point of a surviving half lies on the source, despite losing half
    // the sheet. Reverse distance must distinguish it from a complete quad.
    const std::vector<V> p = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 } };
    SurfaceAnalysis reference(SurfaceMesh(p, { { 0, 1, 2 }, { 0, 2, 3 } }), 1, 90, 0, 0);
    require(reference.surfaceDistanceSquared(V(.75, .25, 0)) < 1e-20, "coplanar fragment left its source");
    require(reference.missingSurfaceError(p, { { 0, 1, 2 } }) > .01, "missing half was scored as a complete sheet");
    require(reference.missingSurfaceError(p, { { 0, 1, 2, 3 } }) < 1e-20, "complete quad has reverse fitting error");
    require(std::isinf(reference.missingSurfaceError({}, {})), "empty output has finite fitting error");
}

static void rimWorkflow()
{
    // Duplicate adjacency is not an authored opening, but must retain the
    // established rim constraint independently of fold-workflow classification.
    const SurfaceMesh duplicatedMesh(
        { { 0, 0, 0 }, { 4, 0, 0 }, { 0, 4, 0 }, { 0, 0, 4 } },
        { { 0, 2, 1 }, { 0, 1, 3 }, { 0, 3, 2 }, { 1, 2, 3 }, { 0, 2, 1 } });
    const SurfaceAnalysis duplicated(duplicatedMesh, 1, 90, 0, 0, false, false);
    require(!duplicated.openFoldWorkflow() && duplicated.supportsRimConstraints(),
        "authored-opening classification changed the existing rim workflow");
    const SurfaceAnalysis clothDuplicated(duplicatedMesh, 1, 90, 0, 0, false, false, false, nullptr, true);
    require(!clothDuplicated.openFoldWorkflow() && clothDuplicated.supportsRimConstraints(),
        "cloth mode treated duplicate adjacency as an authored opening");
    // A cut solid has one planar opening; bending that opening makes it a sheet.
    const std::vector<std::vector<size_t>> triangles = {
        { 0, 1, 4 }, { 1, 2, 4 }, { 2, 3, 4 }, { 3, 0, 4 }
    };
    for (double scale : { 1., 7. }) {
        auto pose = [&](const V& p) { return scale * V(p.z(), p.x(), p.y()); };
        std::vector<V> points = { V(-1, -1, 0), V(1, -1, 0), V(1, 1, 0), V(-1, 1, 0), V(0, 0, 1) };
        for (auto& p : points)
            p = pose(p);
        SurfaceMesh cut(points, triangles);
        SurfaceAnalysis solid(cut, .1 * scale, 90, 0, 0, true, false, false, nullptr, true);
        require(solid.supportsRimConstraints() && !solid.openFoldWorkflow(), "planar cut lost its established rim workflow");
        points[0] += pose(V(0, 0, 1));
        const SurfaceMesh sheetMesh(points, triangles);
        SurfaceAnalysis legacySheet(sheetMesh, .1 * scale, 90, 0, 0, true, false);
        require(!legacySheet.supportsRimConstraints() && !legacySheet.openFoldWorkflow(),
            "default analysis enabled cloth guidance for a spatial opening");
        SurfaceAnalysis sheet(sheetMesh, .1 * scale, 90, 0, 0, true, false, false, nullptr, true);
        require(!sheet.supportsRimConstraints() && sheet.openFoldWorkflow(), "spatial opening lost fold guidance");
        // Preparation may flatten a rim; its local shape must not change the source policy.
        SurfaceAnalysis prepared(cut, .1 * scale, 90, 0, 0, true, false, false, &sheet);
        require(prepared.supportsRimConstraints() && prepared.openFoldWorkflow(), "prepared rim replaced the source workflow");
        SurfaceAnalysis preparedLegacy(sheetMesh, .1 * scale, 90, 0, 0, true, false, false, &legacySheet, true);
        require(!preparedLegacy.openFoldWorkflow(), "preparation enabled cloth guidance for a legacy source");
    }
}

static void roundRimRecovery()
{
    const size_t sides = 64, rings = 8;
    std::vector<V> source { V(0, 0, 0) };
    std::vector<std::vector<size_t>> triangles, quads;
    for (size_t r = 1; r <= rings; ++r)
        for (size_t i = 0; i < sides; ++i) {
            const double angle = 2 * M_PI * i / sides;
            source.emplace_back(double(r) / rings * std::cos(angle), double(r) / rings * std::sin(angle), 0);
        }
    for (size_t i = 0; i < sides; ++i) {
        quads.push_back({ 0, 1 + i, 1 + (i + 1) % sides });
        triangles.push_back(quads.back());
    }
    for (size_t r = 1; r < rings; ++r)
        for (size_t i = 0; i < sides; ++i) {
            const size_t a = 1 + (r - 1) * sides + i, b = 1 + (r - 1) * sides + (i + 1) % sides;
            quads.push_back({ a, a + sides, b + sides, b });
            triangles.push_back({ a, a + sides, b + sides });
            triangles.push_back({ a, b + sides, b });
        }
    auto output = source;
    for (size_t i = 0; i < sides; ++i)
        output[1 + (rings - 1) * sides + i] *= i % 2 ? .84 : .98;
    const auto dented = output;
    SurfaceAnalysis round(SurfaceMesh(source, triangles), .08, 90, 0, 0, false, false);
    require(round.restoreRoundBoundary(output, quads) == sides, "round rim recovery declined a dented disk");
    for (size_t i = 0; i < sides; ++i)
        require(std::fabs(output[1 + (rings - 1) * sides + i].length() - 1) < 1e-8, "round rim retained a dent");
    require(output[0].lengthSquared() == 0, "rim recovery moved the fixed interior");
    auto tilted = source;
    output = dented;
    const auto pose = [](V v) { return V(3 + 7 * v.z(), 2 + 7 * v.x(), 1 + 7 * v.y()); };
    for (auto& v : tilted) v = pose(v);
    for (auto& v : output) v = pose(v);
    SurfaceAnalysis rotated(SurfaceMesh(tilted, triangles), .56, 90, 0, 0, false, false);
    require(rotated.restoreRoundBoundary(output, quads) == sides, "rim recovery depends on axis or scale");
    for (size_t i = 0; i < sides; ++i)
        require(std::fabs((output[1 + (rings - 1) * sides + i] - pose(V())).length() - 7) < 1e-7, "posed rim retained a dent");
    for (auto& v : source) v = V(v.x(), .7 * v.y(), v.z());
    SurfaceAnalysis ellipse(SurfaceMesh(source, triangles), .08, 90, 0, 0, false, false);
    output = dented;
    require(ellipse.restoreRoundBoundary(output, quads) == 0, "elliptical opening was forced round");
    for (size_t i = 0; i < output.size(); ++i)
        require((output[i] - dented[i]).lengthSquared() == 0, "rejected rim recovery was not atomic");
}

static void residualHoles()
{
    // A concave prism: removing the top requires more than a convex fan.
    std::vector<V> points;
    for (double z : { 0., 1. })
        for (auto xy : { V(0, 0, 0), V(2, 0, 0), V(2, 1, 0), V(1, 1, 0), V(1, 2, 0), V(0, 2, 0) })
            points.emplace_back(xy.x(), xy.y(), z);
    std::vector<std::vector<size_t>> source;
    for (size_t i = 0; i < 6; ++i) {
        const size_t j = (i + 1) % 6;
        source.push_back({ i, j, j + 6 });
        source.push_back({ i, j + 6, i + 6 });
    }
    const std::vector<std::vector<size_t>> cap { { 0, 1, 3 }, { 1, 2, 3 }, { 0, 3, 5 }, { 3, 4, 5 } };
    for (auto face : cap) {
        std::reverse(face.begin(), face.end());
        source.push_back(face);
    }
    const auto open = source;
    for (auto face : cap) {
        for (auto& v : face) v += 6;
        source.push_back(face);
    }
    SurfaceAnalysis closed(SurfaceMesh(points, source), .2, 90, 0, 0, false, false);
    const auto watertight = [](const std::vector<std::vector<size_t>>& faces) {
        std::map<std::pair<size_t, size_t>, size_t> edges;
        for (const auto& f : faces)
            for (size_t i = 0; i < f.size(); ++i)
                ++edges[{ f[i], f[(i + 1) % f.size()] }];
        for (const auto& e : edges)
            require(e.second == 1 && edges[{ e.first.second, e.first.first }] == 1, "repair left an opening or duplicated a directed edge");
    };
    auto output = open;
    require(closed.closeBoundaryHoles(points, output) == 1, "concave hole was not filled");
    require(std::equal(open.begin(), open.end(), output.begin()), "hole repair changed existing faces");
    watertight(output);
    require(closed.closeBoundaryHoles(points, output) == 0, "closed mesh was modified again");
    SurfaceAnalysis intentional(SurfaceMesh(points, open), .2, 90, 0, 0, false, false);
    output = open;
    require(intentional.closeBoundaryHoles(points, output) == 0 && output == open, "source opening was capped");
    // A source orientation seam is not an open border.
    auto inconsistent = source;
    std::reverse(inconsistent[0].begin(), inconsistent[0].end());
    SurfaceAnalysis seam(SurfaceMesh(points, inconsistent), .2, 90, 0, 0, false, false);
    output = open;
    require(seam.closeBoundaryHoles(points, output) == 1, "source orientation seam hid an accidental hole");
    watertight(output);
    // Split one corner of an existing triangle without changing its position.
    points.push_back(points[0]);
    output = source;
    output[0][0] = points.size() - 1;
    require(closed.closeBoundaryHoles(points, output) == 1, "coincident seam was not sewn");
    require(output.size() == source.size(), "seam repair introduced a zero-area cap");
    watertight(output);
}

static void competingFold()
{
    // A two-edge curved crease meets a short, nearly parallel 78-degree
    // branch. The branch lacks geometric strength but still conflicts with
    // fixing the crease as a grid axis. Coordinates form a manifold patch.
    const std::vector<V> p = {
        { -1.1739, 1.5743, -1.2032 }, { -0.9719, 1.5606, -1.3530 },
        { -0.2848, 0.1854, 0.1866 }, { -0.7125, 0.8861, -0.5067 },
        { -0.0240, 0.0766, 0.2073 }, { -0.5075, 0.8640, -0.5964 },
        { 0.0128, -1.9678, 1.1396 }, { -1.1317, 1.5229, -1.5903 },
        { 0.5004, -2.3070, 1.4073 }, { 0.4967, -0.9420, 0.7634 },
        { 0.8051, -1.4823, 1.0759 }, { 0.8672, -1.3979, 0.9903 },
        { 0.5962, -2.5440, 1.2930 }, { 0.5539, -1.2169, 0.7005 },
        { 0.1299, -0.6350, 0.2416 }, { 0.5210, -0.6479, 0.5285 },
        { 0.0000, 0.0000, 0.0000 }, { -0.6024, 0.7218, -0.7605 },
        { -0.6768, 0.6547, -0.8524 },
    };
    const std::vector<std::vector<size_t>> t = {
        { 3, 2, 4 }, { 4, 5, 3 }, { 0, 3, 5 }, { 5, 1, 0 },
        { 6, 4, 2 }, { 9, 4, 6 }, { 6, 8, 9 }, { 10, 9, 8 },
        { 14, 13, 12 }, { 16, 4, 9 }, { 9, 15, 16 }, { 15, 9, 10 },
        { 10, 11, 15 }, { 13, 15, 11 }, { 14, 16, 15 }, { 15, 13, 14 },
        { 7, 17, 18 }, { 17, 7, 1 }, { 1, 5, 17 }, { 16, 17, 5 },
        { 5, 4, 16 }, { 14, 18, 17 }, { 17, 16, 14 },
    };
    for (double scale : { 1., 7. }) {
        std::vector<V> q;
        for (const auto& v : p)
            q.push_back(scale * (scale == 1 ? v : V(v.z(), v.x(), v.y())));
        SurfaceMesh mesh(q, t);
        SurfaceAnalysis analysis(mesh, 1.4 * scale, 90, 1, 1, true, true, false, nullptr, true);
        const auto guidance = analysis.transfer(mesh, true);
        const auto chainOnEdge = [&](size_t a, size_t b) -> const SurfaceAnalysis::Chain* {
            for (const auto& chain : analysis.chains())
                for (size_t c : chain.corners) {
                    const size_t x = mesh.cornerVertex(c), y = mesh.cornerVertex(mesh.nextCorner(c));
                    if ((x == a && y == b) || (x == b && y == a))
                        return &chain;
                }
            return nullptr;
        };
        const auto* fold = chainOnEdge(15, 16);
        const auto* branch = chainOnEdge(9, 16);
        require(fold && fold->strength == 1 && branch && branch->strength == 0,
            "invalid strong-fold/weak-competitor surface");
        for (size_t c : fold->corners)
            require(!guidance.featureCorners[c], "zero-strength competitor was ignored when fixing a fold axis");
        for (size_t c = 0; c < mesh.cornerCount(); ++c)
            if (mesh.isBoundaryCorner(c))
                require(guidance.featureCorners[c], "directional conflict removed an authored border axis");
        const V midpoint = .5 * (q[15] + q[16]);
        const V displaced = midpoint + .02 * scale * mesh.faceNormal(fold->corners.front() / 3);
        const auto binding = analysis.bindCurve(displaced, .05 * scale);
        require(binding.chain != SurfaceMesh::npos && analysis.chains()[binding.chain].strength == 1,
            "directional conflict removed geometric crease protection");
        require((analysis.projectCurve(binding, displaced) - midpoint).length() < 1e-10 * scale,
            "directional conflict changed crease projection");
    }
}

static void foldedStrip(bool curved, bool closed = false, bool varying = false)
{
    // The constant control is a uniform 120-degree crease. The varying fold
    // changes its cross-section smoothly from 50 to 130 degrees along the arc.
    // Only the middle row lies off the authored open border.
    std::vector<V> p;
    std::vector<std::vector<size_t>> t;
    for (size_t i = 0; i <= 12; ++i)
        for (int j = -1; j <= 1; ++j) {
            const double a = i * M_PI / 18, r = 2 + .15 * j;
            const double phi = (varying ? 50 + 80 * double(i) / 12 : 120) * M_PI / 180;
            const double height = .15 * std::tan(phi / 2) * std::abs(j);
            p.push_back(curved ? V(r * std::cos(a), r * std::sin(a), height)
                               : V(2 * a, .15 * j, height));
        }
    for (size_t i = 0; i < 12; ++i) {
        for (size_t j = 0; j < 2; ++j) {
            const size_t a = 3 * i + j, b = a + 3;
            t.push_back({ a, b, b + 1 });
            t.push_back({ a, b + 1, a + 1 });
        }
        if (closed) {
            const size_t a = 3 * i;
            t.push_back({ a + 2, a + 5, a + 3 });
            t.push_back({ a + 2, a + 3, a });
        }
    }
    if (closed) {
        t.push_back({ 0, 1, 2 });
        t.push_back({ 38, 37, 36 });
    }
    SurfaceMesh mesh(p, t);
    SurfaceAnalysis analysis(mesh, .2, 90, 1, 1, true, true, false, nullptr, true);
    const auto guidance = analysis.transfer(mesh, true);
    size_t foldEdges = 0, foldAxes = 0, borderEdges = 0, physical = 0;
    for (size_t c = 0; c < mesh.cornerCount(); ++c) {
        if (mesh.isBoundaryCorner(c)) {
            ++borderEdges;
            require(guidance.featureCorners[c], "curved authored border lost its grid constraint");
        }
        if (mesh.cornerVertex(c) % 3 == 1 && mesh.cornerVertex(mesh.nextCorner(c)) % 3 == 1) {
            ++foldEdges;
            foldAxes += !!guidance.featureCorners[c];
        }
    }
    for (const auto& f : analysis.faces()) physical += f.major > .1;
    const V midpoint = (p[16] + p[19]) * .5, displaced = midpoint + V(0, 0, .02);
    const auto binding = analysis.bindCurve(displaced, .05);
    require(binding.chain != SurfaceMesh::npos, "fold lost its protected geometry");
    require((analysis.projectCurve(binding, displaced) - midpoint).length() < 1e-12,
        "protected fold no longer restores displaced vertices");
    require(foldEdges == 24 && borderEdges == (closed ? 0 : 28), "invalid folded-strip topology");
    const bool curvatureDriven = curved && !closed && varying;
    require(curvatureDriven ? physical > t.size() / 2 : physical == 0,
        "curved cloth fold lost physical curvature samples, or protected control gained curvature");
    require(foldAxes == (curvatureDriven ? 0 : foldEdges),
        "cloth fold forced a grid axis, or protected control lost its axis");
    if (curvatureDriven) {
        // A nearby interior sample is not curve-bound, but its source walk
        // must still respect a protected fold that does not force a grid axis.
        const V inner = (p[15] + p[18]) * .5, outer = (p[17] + p[20]) * .5;
        // Barycentric in triangle (15, 19, 16), even when the band twists.
        std::vector<V> output = { .8 * midpoint + .2 * p[15], p[17], p[20] };
        const auto before = output;
        require(analysis.bindCurve(output[0], .02).chain == SurfaceMesh::npos,
            "fold relaxation sample was unexpectedly curve-bound");
        analysis.relaxSurface(output, { { 1, 2 }, { 0, 2 }, { 0, 1 } },
            { false, true, true }, { { 0, 1, 2 } }, 5);
        require(V::dotProduct(output[0] - midpoint, outer - inner) <= 1e-12,
            "relaxation crossed a protected nondirectional fold");
        require(analysis.surfaceDistanceSquared(output[0]) < 1e-20,
            "protected-fold relaxation left the source surface");
        require(V::dotProduct(V::crossProduct(before[1] - before[0], before[2] - before[0]),
                    V::crossProduct(output[1] - output[0], output[2] - output[0])) > 0,
            "protected-fold relaxation reversed an incident face");
        require((output[1] - before[1]).length() < 1e-12 && (output[2] - before[2]).length() < 1e-12,
            "protected-fold relaxation moved locked neighbors");
    }
}

static void supportedBend()
{
    // A gently curved strip approaches a smooth 103-degree bend. The bend
    // first enters the widest support at row 48: finer radii see only the
    // gentle background. That supported sample should affect spacing while
    // remaining a soft guide, since adjacent scales do not confirm it yet.
    const size_t columns = 101, rows = 21, shoulder = 48 * 40 + 20;
    for (double scale : { 1., 7. }) {
        std::vector<std::vector<size_t>> triangles;
        for (size_t i = 0; i + 1 < columns; ++i)
            for (size_t j = 0; j + 1 < rows; ++j) {
                const size_t a = i * rows + j, b = a + rows;
                triangles.push_back({ a, b, b + 1 });
                triangles.push_back({ a, b + 1, a + 1 });
            }
        const auto strip = [&](bool bend) {
            std::vector<V> points;
            double x = 0, z = 0, previous = 0;
            for (size_t i = 0; i < columns; ++i) {
                const double s = -.2 + .004 * i;
                const double angle = 2 * s + (bend ? .9 * (1 + std::tanh((s - .06) / .01)) : 0);
                if (i) {
                    x += .004 * std::cos((angle + previous) * .5);
                    z += .004 * std::sin((angle + previous) * .5);
                }
                previous = angle;
                for (size_t j = 0; j < rows; ++j) {
                    const double y = .1 * (double(j) - 10);
                    points.push_back(scale * (scale == 1 ? V(x, y, z) : V(x, z, -y)));
                }
            }
            return points;
        };
        SurfaceMesh curvedMesh(strip(true), triangles), gentleMesh(strip(false), triangles);
        SurfaceAnalysis curved(curvedMesh, .05 * scale, 90, 1, 1, true, true, true, nullptr, true);
        SurfaceAnalysis gentle(gentleMesh, .05 * scale, 90, 1, 1, true, true, true, nullptr, true);
        const auto& guide = curved.faces()[shoulder];
        const auto& background = gentle.faces()[shoulder];
        require(background.scale == 1 && background.major * scale < 3,
            "gentle control acquired bend density");
        require(guide.major > 3 * background.major && guide.scale < .98 && guide.ratio < .98,
            "supported bend was discarded for under-resolved background curvature");
        require(guide.confidence > 0 && guide.confidence < 1,
            "one supported scale became a fixed curvature root");
        require(curved.faces()[65 * 40 + 20].confidence == 1,
            "bend interior lost its confirmed curvature root");
    }
}

static void sourceFlow()
{
    const std::vector<std::vector<size_t>> triangles { { 0, 1, 2 }, { 0, 2, 3 } };
    for (double scale : { 1., 7. }) {
        const auto transform = [&](V v) { return scale * (scale == 1 ? v : V(v.z(), v.x(), v.y())); };
        std::vector<V> source;
        for (V v : { V(-1, -1, 0), V(1, -1, 0), V(1, 1, 0), V(-1, 1, 0) })
            source.push_back(transform(v));
        SurfaceAnalysis reference(SurfaceMesh(source, triangles), .1 * scale, 90, 0, 0, true);
        std::vector<SurfaceGuidance::Face> flow(2);
        for (auto& face : flow) {
            face.direction = transform(V(1, 0, 0));
            face.confidence = 1;
        }
        struct Grid { std::vector<V> points; std::vector<std::vector<size_t>> faces; };
        const auto grid = [&](size_t nx, size_t ny, double angle, double right = 2.) {
            Grid result;
            for (size_t j = 0; j <= ny; ++j)
                for (size_t i = 0; i <= nx; ++i) {
                    const double x = -2 + (right + 2) * i / nx, y = -2 + 4. * j / ny;
                    result.points.push_back(transform(V(std::cos(angle) * x - std::sin(angle) * y,
                        std::sin(angle) * x + std::cos(angle) * y, 0)));
                }
            for (size_t j = 0; j < ny; ++j)
                for (size_t i = 0; i < nx; ++i) {
                    const size_t a = j * (nx + 1) + i;
                    result.faces.push_back({ a, a + 1, a + nx + 2, a + nx + 1 });
                }
            return result;
        };
        // Analytic directions isolate layout scoring from curvature estimation.
        // Rotating/subdividing the same regular grid must not buy better flow.
        for (double angle : { 0., M_PI / 4, M_PI / 15 })
            for (size_t nx : { size_t(1), size_t(8), size_t(3) }) {
                const auto output = grid(nx, nx == 3 ? 12 : nx, angle);
                const auto fit = reference.measureLayoutFit(output.points, output.faces, &flow);
                require(std::fabs(fit.flowError - std::pow(std::sin(2 * angle), 2)) < 1e-12,
                    "flow depends on subdivision, aspect, world rotation or scale");
                require(std::fabs(fit.flowWeight - 4 * scale * scale) < 1e-10,
                    "candidate density changed source flow weight");
            }
        auto output = grid(8, 8, M_PI / 15);
        for (auto& face : output.faces) std::reverse(face.begin(), face.end());
        for (V axis : { V(-1, 0, 0), V(0, 1, 0) }) {
            for (auto& face : flow) face.direction = transform(axis);
            require(std::fabs(reference.measureLayoutFit(output.points, output.faces, &flow).flowError - std::pow(std::sin(2 * M_PI / 15), 2)) < 1e-12,
                "axis sign, quarter turn or polygon winding changed cross flow");
        }
        const auto missing = grid(4, 4, 0, 0);
        const auto fit = reference.measureLayoutFit(missing.points, missing.faces, &flow);
        require(std::fabs(fit.flowWeight - 4 * scale * scale) < 1e-10 && fit.flowError < 1e-12,
            "missing aligned strip changed source weighting");
        require(std::fabs(fit.missingError - scale * scale / 18) < 1e-10,
            "flow scoring hid a missing strip from the existing fitting measure");
        require(fit.missingError == reference.missingSurfaceError(missing.points, missing.faces),
            "optional flow changed geometric fitting");
        const auto flat = reference.measureLayoutFit(output.points, output.faces, &reference.faces());
        require(flat.flowWeight == 0 && flat.flowError == 0, "flat source invented directional evidence");
        for (auto& face : flow) face.direction = V();
        const auto zero = reference.measureLayoutFit(output.points, output.faces, &flow);
        require(zero.flowWeight == 0 && zero.flowError == 0, "zero direction did not fall back to geometry");
    }
}

static void invalidFlowFit()
{
    const std::vector<V> source { V(-1, -1, 0), V(1, -1, 0), V(1, 1, 0), V(-1, 1, 0) };
    SurfaceAnalysis reference(SurfaceMesh(source, { { 0, 1, 2 }, { 0, 2, 3 } }), .1, 90, 0, 0, false, false);
    const std::vector<std::vector<size_t>> quad { { 0, 1, 2, 3 } };
    std::vector<SurfaceGuidance::Face> flow(2);
    for (auto& face : flow) { face.direction = V(1, 0, 0); face.confidence = 1; }
    auto far = source;
    for (auto& v : far) v += V(1e200, 1e200, 1e200);
    require(std::isinf(reference.measureLayoutFit(far, quad, &flow).flowError),
        "failed nearest-face lookup did not fail closed");
    auto nonfinite = source;
    nonfinite[0] = V(std::numeric_limits<double>::quiet_NaN(), 0, 0);
    require(std::isinf(reference.measureLayoutFit(nonfinite, quad, &flow).flowError),
        "nonfinite output geometry has a finite flow score");
    require(std::isinf(reference.measureLayoutFit(source, { { 0, 1, 2, 4 } }, &flow).flowError),
        "invalid output indices reached the candidate BVH");
}

static void coherentQuadRelaxation()
{
    // A thin warped quad with a closed one-ring around its only movable vertex.
    for (double scale : { .25, 1., 7. })
        for (bool rotate : { false, true })
            for (double fraction : { .1, .3, 1. }) {
                std::vector<V> points = {
                    { .08610773016252324, -.43752701009304035, .24165775902321926 },
                    { -.11330633235041171, .39788264620953756, -.23873349134840416 },
                    { -.12217656198252361, .4308693265025071, -.25167236535143317 },
                    { .14937516417041208, -.39122496261897316, .24874809767661804 }
                };
                const V goal(-.12647046828443562, .4284101238856443, -.2570446132709799);
                const V normal(.7874016444877138, -.029775616421176657, -.615720767007889);
                const V tangent = (goal - points[2]).normalized();
                const V across = V::crossProduct(normal, tangent).normalized();
                const V wanted = points[2] + fraction * (goal - points[2]);
                // Supporting neighbors prescribe the undamped mean; projection stays on one sheet.
                V center = ((2 * wanted - points[2]) * 4 - points[1] - points[3]) * .5;
                center = center - normal * V::dotProduct(center - points[2], normal);
                points.push_back(center - .5 * across);
                points.push_back(center + .5 * across);
                std::vector<V> source = { points[2] - 10 * tangent - 10 * across,
                    points[2] + 20 * tangent - 10 * across, points[2] - 10 * tangent + 20 * across };
                auto transform = [&](const V& p) { return scale * (rotate ? V(p.z(), p.x(), p.y()) : p); };
                for (auto& p : points)
                    p = transform(p);
                for (auto& p : source)
                    p = transform(p);
                const std::vector<std::vector<size_t>> polygons = { { 0, 1, 2, 3 }, { 2, 1, 4 }, { 2, 4, 5 }, { 2, 5, 3 } };
                std::vector<std::unordered_set<size_t>> neighbors(points.size());
                for (const auto& polygon : polygons)
                    for (size_t i = 0; i < polygon.size(); ++i) {
                        size_t a = polygon[i], b = polygon[(i + 1) % polygon.size()];
                        neighbors[a].insert(b);
                        neighbors[b].insert(a);
                    }
                const auto before = points;
                SurfaceAnalysis analysis(SurfaceMesh(source, { { 0, 1, 2 } }), scale, 90, 0, 0, false, false, false, nullptr, true);
                for (size_t start = 0; start < 4; ++start) {
                    auto cyclic = polygons;
                    std::rotate(cyclic[0].begin(), cyclic[0].begin() + start, cyclic[0].end());
                    auto output = before;
                    analysis.relaxSurface(output, neighbors, { true, true, false, true, true, true }, cyclic, 1);
                    const V expected = fraction < 1 ? transform(wanted) : before[2];
                    require((output[2] - expected).length() < 1e-10 * scale,
                        "quad start changed acceptance of a valid concave or invalid move");
                    for (size_t i = 0; i < output.size(); ++i)
                        if (i != 2)
                            require((output[i] - before[i]).length() == 0, "surface relaxation moved a locked support");
                    if (fraction < 1)
                        require((output[2] - before[2]).length() > 1e-5 * scale, "quad safeguard disabled valid smoothing");
                }
            }
}



static void semanticCurvePaths()
{
    // Geometric corner splits must not hide a weak open path's short bridge.
    // A loop assembled from open pieces still has no semantic endpoints.
    for (bool closed : { false, true })
        for (double degrees : { 38., 100. })
            for (double scale : { 1., 7. }) {
                const double width = closed ? .2 : .01;
                const double height = width * std::tan(degrees * M_PI / 360);
                const std::vector<V> centers = closed
                    ? std::vector<V>{ V(-1, -1, 0), V(1, -1, 0), V(1, 1, 0), V(-1, 1, 0) }
                    : std::vector<V>{ V(0, 0, 0), V(1, 0, 0), V(1, .1, 0), V(0, .1, 0) };
                const std::vector<V> offsets = closed ? centers
                    : std::vector<V>{ V(0, 1, 0), V(-1, 1, 0), V(-1, -1, 0), V(0, -1, 0) };
                auto transform = [&](const V& v) { return scale * (scale == 1 ? v : V(v.z(), v.x(), v.y())); };
                std::vector<V> points;
                std::vector<std::vector<size_t>> triangles;
                for (size_t i = 0; i < centers.size(); ++i)
                    for (int j = -1; j <= 1; ++j)
                        points.push_back(transform(centers[i] + width * j * offsets[i] + V(0, 0, height * std::abs(j))));
                for (size_t i = 0; i < (closed ? 4 : 3); ++i)
                    for (size_t j = 0; j < 2; ++j) {
                        const size_t a = 3 * i + j, b = 3 * ((i + 1) % 4) + j;
                        triangles.push_back({ a, b, b + 1 });
                        triangles.push_back({ a, b + 1, a + 1 });
                    }
                SurfaceMesh mesh(points, triangles);
                SurfaceAnalysis analysis(mesh, .2 * scale, 90, 1, 1, true, true, false, nullptr, true);
                const auto guidance = analysis.transfer(mesh, true);
                size_t weak = 0;
                for (const auto& chain : analysis.chains()) {
                    if (chain.strength > 0 && chain.strength < 1) {
                        ++weak;
                        require(chain.directional == closed, "semantic path direction ignored its open/closed topology");
                    }
                    if (chain.strength >= 1)
                        require(chain.directional, "semantic grouping removed a hard crease or border direction");
                    for (size_t c : chain.corners) {
                        const bool selected = chain.strength > 0 && chain.directional;
                        require(bool(guidance.featureCorners[c]) == selected, "semantic direction did not reach corner guidance");
                        if (mesh.oppositeCorner(c) != SurfaceMesh::npos)
                            require(guidance.featureCorners[c] == guidance.featureCorners[mesh.oppositeCorner(c)], "semantic path made asymmetric edge constraints");
                    }
                }
                require(weak == (degrees > 90 ? 0 : closed ? 4 : 2), "semantic grouping changed geometric feature retention");
                const V midpoint = closed ? V(0, -1, 0) : V(.5, 0, 0);
                const auto binding = analysis.bindCurve(transform(midpoint), .001 * scale);
                require(binding.chain != SurfaceMesh::npos, "semantic grouping removed geometric curve binding");
                const V along = midpoint + V(.1, 0, 0);
                require((analysis.projectCurve(binding, transform(along + V(0, .001, .001))) - transform(along)).length() < 1e-10 * scale,
                    "semantic grouping changed curve projection ownership");
                if (!closed && degrees < 90)
                    require(analysis.bindCurve(transform(V(1, .05, 0)), .001 * scale).chain == SurfaceMesh::npos,
                        "semantic grouping reactivated the unsupported geometric bridge");
            }
}


static void semanticPathGuards()
{
    // Two long automatic arms meet across two unsupported short side branches.
    // Varying the arm dihedrals distinguishes a protected bend from a strong crease.
    for (double degrees : { 75., 100. }) {
        const double h = .01 * std::tan(degrees * M_PI / 360);
        const V a(-1, 0, 0), b(std::cos(50 * M_PI / 180), std::sin(50 * M_PI / 180), 0);
        const V c(0, .01, h), d(0, -.01, h);
        const std::vector<V> points = { V(), c, d, a, b, a + c, a + d, b + c, b + d };
        const std::vector<std::vector<size_t>> triangles = {
            { 3, 0, 1 }, { 3, 1, 5 }, { 3, 2, 0 }, { 3, 6, 2 },
            { 0, 4, 1 }, { 4, 7, 1 }, { 0, 2, 4 }, { 4, 2, 8 }
        };
        SurfaceMesh mesh(points, triangles);
        SurfaceAnalysis analysis(mesh, .2, 150, 1, 1, true, true, false, nullptr, true);
        size_t arms = 0, branches = 0;
        for (const auto& chain : analysis.chains())
            if (!mesh.isBoundaryCorner(chain.corners.front())) {
                if (chain.strength > 0 && chain.strength < 1) {
                    ++arms;
                    require(chain.directional == (degrees > 90), "semantic grouping confused a strong crease with a protected bend");
                } else if (chain.strength == 0) {
                    ++branches;
                }
            }
        require(arms == 2 && branches == 2, "semantic grouping changed arm or unsupported branch retention");
    }
}

int main()
{
    try {
        rimWorkflow();
        semanticPathGuards();
        semanticCurvePaths();
        coherentQuadRelaxation();
        sourceFlow();
        invalidFlowFit();
        competingFold();
        supportedBend();
        foldedStrip(false);
        foldedStrip(true, true);
        foldedStrip(true);
        foldedStrip(true, false, true);
        residualHoles();
        roundRimRecovery();
        missingSurface();
        curvedProtectedTube();
        protectedJunctions();
        cylinder();
        cylinder(true);
        creaseAndTransfer();
        connectedRelaxation();
        std::cout << "Surface analysis checks passed\n";
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
