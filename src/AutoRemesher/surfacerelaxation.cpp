#include "surfacerelaxation.h"
#include <algorithm>
#include <axisalignedboundingboxtree.h>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <numeric>

namespace AutoRemesher {
namespace {
    double dot(const Vector3& a, const Vector3& b) { return Vector3::dotProduct(a, b); }
    bool finite(const Vector3& p) { return std::isfinite(p.x()) && std::isfinite(p.y()) && std::isfinite(p.z()); }
    double triangleDistance2(const Vector3& p, const Vector3& a, const Vector3& b, const Vector3& c)
    {
        if ((p - a).lengthSquared() == 0 || (p - b).lengthSquared() == 0 || (p - c).lengthSquared() == 0)
            return 0;
        const Vector3 n = Vector3::crossProduct(b - a, c - a);
        if (n.lengthSquared() > 0) {
            const Vector3 q = p - n * (dot(p - a, n) / n.lengthSquared());
            if (dot(Vector3::crossProduct(b - a, q - a), n) >= 0 && dot(Vector3::crossProduct(c - b, q - b), n) >= 0 && dot(Vector3::crossProduct(a - c, q - c), n) >= 0)
                return (p - q).lengthSquared();
        }
        double best = std::numeric_limits<double>::infinity();
        for (const auto& edge : { std::make_pair(a, b), std::make_pair(b, c), std::make_pair(c, a) }) {
            const Vector3 d = edge.second - edge.first;
            const double t = d.lengthSquared() > 0 ? std::max(0., std::min(1., dot(p - edge.first, d) / d.lengthSquared())) : 0;
            best = std::min(best, (p - edge.first - d * t).lengthSquared());
        }
        return best;
    }
    struct Frame {
        std::array<Vector3, 2> axis, gradient;
        std::array<double, 2> size {};
        bool valid = false;
    };
    struct Edge {
        size_t uses = 0;
        Vector3 direction;
        double minimum = 0, maximum = 0;
        bool valid = false;
    };
    struct Endpoint {
        Vector3 direction;
        double minimum = 0, maximum = 0;
        bool valid = false;
    };
}

Vector3 SurfaceRelaxationStencil::target(const std::vector<Vector3>& points, size_t vertex) const
{
    if (vertex >= points.size())
        return Vector3();
    const Vector3 center = points[vertex];
    if (!valid)
        return center;
    Vector3 directional, metric;
    for (size_t axis = 0; axis < 2; ++axis) {
        const size_t a = axis, b = axis + 2;
        if (neighbors[a] >= points.size() || neighbors[b] >= points.size())
            return center;
        const Vector3 first = points[neighbors[a]], second = points[neighbors[b]];
        const double l0 = std::max(minimum[a], std::min(maximum[a], (first - center).length()));
        const double l1 = std::max(minimum[b], std::min(maximum[b], (second - center).length()));
        if (!(l0 + l1 > 0) || !std::isfinite(l0 + l1))
            return center;
        directional += (first - outward[a] * l0 + second - outward[b] * l1) * .5;
        metric += first * (l1 / (l0 + l1)) + second * (l0 / (l0 + l1));
    }
    const Vector3 result = (directional * blend + metric * (1 - blend)) * .5;
    return finite(result) ? result : center;
}

std::vector<SurfaceRelaxationStencil> buildSurfaceRelaxationStencils(
    const std::vector<Vector3>& preparedVertices,
    const std::vector<std::vector<size_t>>& preparedTriangles,
    const std::vector<std::vector<Vector2>>& finalTriangleUvs,
    const std::vector<SurfaceRelaxationFace>& faces,
    const std::vector<Vector3>& outputVertices,
    const std::vector<std::vector<size_t>>& outputPolygons,
    const std::map<std::pair<size_t, size_t>, size_t>* edgeUseCount)
{
    std::vector<SurfaceRelaxationStencil> result(outputVertices.size());
    const size_t count = preparedTriangles.size(), none = std::numeric_limits<size_t>::max();
    if (!count || count != faces.size() || count != finalTriangleUvs.size())
        return result;
    for (const auto& p : preparedVertices)
        if (!finite(p))
            return result;
    for (const auto& p : outputVertices)
        if (!finite(p))
            return result;
    std::vector<Frame> frames(count);
    std::vector<AxisAlignedBoudingBox> boxes(count);
    AxisAlignedBoudingBox bounds;
    for (size_t f = 0; f < count; ++f) {
        const auto& t = preparedTriangles[f];
        if (t.size() != 3 || *std::max_element(t.begin(), t.end()) >= preparedVertices.size())
            return result;
        for (size_t v : t) {
            const auto& p = preparedVertices[v];
            const ::Vector3 q(p.x(), p.y(), p.z());
            boxes[f].update(q);
            bounds.update(q);
        }
        boxes[f].updateCenter();
        const auto& uv = finalTriangleUvs[f];
        const auto& data = faces[f];
        auto& frame = frames[f];
        if (uv.size() != 3 || !finite(data.direction) || !(data.spacingU > 0) || !(data.spacingV > 0) || !std::isfinite(data.spacingU + data.spacingV + data.anisotropy))
            continue;
        const Vector3 a = preparedVertices[t[1]] - preparedVertices[t[0]], b = preparedVertices[t[2]] - preparedVertices[t[0]];
        const Vector3 n = Vector3::crossProduct(a, b);
        const double area2 = n.lengthSquared();
        const Vector2 u = uv[1] - uv[0], v = uv[2] - uv[0];
        const double determinant = u.x() * v.y() - u.y() * v.x();
        const double uvScale = std::hypot(u.x(), u.y()) * std::hypot(v.x(), v.y());
        if (!(area2 > 0) || !std::isfinite(determinant) || std::fabs(determinant) <= 1e-12 * uvScale)
            continue;
        const Vector3 normal = n.normalized();
        Vector3 first = data.direction - normal * dot(data.direction, normal);
        if (!(first.length() > 1e-12))
            continue;
        first = first.normalized();
        const Vector3 second = Vector3::crossProduct(normal, first);
        const Vector3 inverseU = (a * v.y() - b * u.y()) / determinant, inverseV = (b * u.x() - a * v.x()) / determinant;
        const bool swap = std::fabs(dot(inverseU, first)) < std::fabs(dot(inverseU, second));
        frame.axis = { swap ? second : first, swap ? first : second };
        frame.size = { swap ? data.spacingV : data.spacingU, swap ? data.spacingU : data.spacingV };
        if (dot(frame.axis[0], inverseU) < 0)
            frame.axis[0] = -frame.axis[0];
        if (dot(frame.axis[1], inverseV) < 0)
            frame.axis[1] = -frame.axis[1];
        const Vector3 ga = Vector3::crossProduct(b, n) / area2, gb = Vector3::crossProduct(n, a) / area2;
        frame.gradient = { ga * u.x() + gb * v.x(), ga * u.y() + gb * v.y() };
        frame.valid = finite(frame.gradient[0]) && finite(frame.gradient[1]);
    }
    // Geometry-only index: no feature classification or curvature is repeated.
    bounds.updateCenter();
    std::vector<size_t> indices(count);
    std::iota(indices.begin(), indices.end(), 0);
    AxisAlignedBoudingBoxTree tree(&boxes, indices, bounds);
    std::vector<size_t> mapped(outputVertices.size(), none);
    double coordinateScale = 0;
    for (size_t k = 0; k < 3; ++k)
        coordinateScale = std::max({ coordinateScale, std::fabs(bounds.lowerBound()[k]), std::fabs(bounds.upperBound()[k]) });
    std::vector<std::pair<double, size_t>> candidates;
    for (size_t i = 0; i < outputVertices.size(); ++i) {
        candidates.clear();
        const Vector3 p = outputVertices[i];
        double best = std::numeric_limits<double>::infinity();
        // Distance arithmetic at a shared edge differs by coordinate roundoff.
        // Retain all numerical ties, then choose an owner against the true minimum.
        constexpr double roundoff = 32 * std::numeric_limits<double>::epsilon();
        const double coordinateError = roundoff * std::max({ coordinateScale, std::fabs(p.x()), std::fabs(p.y()), std::fabs(p.z()) });
        const auto tolerance = [&](double squaredDistance) {
            // Convert the coordinate allowance to squared-distance units, including
            // the first-order term for points away from the source surface.
            return std::isfinite(squaredDistance)
                ? coordinateError * (2 * std::sqrt(squaredDistance) + coordinateError) + roundoff * squaredDistance
                : std::numeric_limits<double>::infinity();
        };
        const auto distance = [&](const AxisAlignedBoudingBox& box) {
            double sum = 0;
            for (size_t k = 0; k < 3; ++k) {
                const double d = std::max({ box.lowerBound()[k] - p[k], 0., p[k] - box.upperBound()[k] });
                sum += d * d;
            }
            return sum;
        };
        std::function<void(const AxisAlignedBoudingBoxTree::Node*)> visit = [&](const AxisAlignedBoudingBoxTree::Node* node) {
            if (!node || distance(node->boundingBox) > best + tolerance(best))
                return;
            if (node->isLeaf())
                for (size_t f : node->boxIndices) {
                    const auto& t = preparedTriangles[f];
                    const double d = triangleDistance2(p, preparedVertices[t[0]], preparedVertices[t[1]], preparedVertices[t[2]]);
                    best = std::min(best, d);
                    if (d <= best + tolerance(best))
                        candidates.emplace_back(d, f);
                }
            else {
                const bool left = distance(node->left->boundingBox) < distance(node->right->boundingBox);
                visit(left ? node->left : node->right);
                visit(left ? node->right : node->left);
            }
        };
        visit(tree.root());
        for (const auto& candidate : candidates)
            if (candidate.first <= best + tolerance(best))
                mapped[i] = std::min(mapped[i], candidate.second);
    }
    using Key = std::pair<size_t, size_t>;
    std::map<Key, Edge> edges;
    if (edgeUseCount)
        for (const auto& edge : *edgeUseCount) {
            if (edge.first.first >= outputVertices.size() || edge.first.second >= outputVertices.size())
                return result;
            edges[edge.first].uses = edge.second;
        }
    std::vector<std::vector<Key>> links(outputVertices.size());
    for (const auto& polygon : outputPolygons) {
        if (polygon.size() < 3 || *std::max_element(polygon.begin(), polygon.end()) >= outputVertices.size())
            return std::vector<SurfaceRelaxationStencil>(outputVertices.size());
        for (size_t k = 0; k < polygon.size(); ++k) {
            const size_t a = polygon[k], b = polygon[(k + 1) % polygon.size()], previous = polygon[(k + polygon.size() - 1) % polygon.size()];
            auto& edge = edges[std::minmax(a, b)];
            if (!edgeUseCount)
                ++edge.uses;
            links[a].push_back({ previous, b });
        }
    }
    const auto endpoint = [&](size_t a, size_t b) {
        Endpoint value;
        const size_t f = mapped[a];
        if (f == none || !frames[f].valid || !(outputVertices[b] - outputVertices[a]).lengthSquared())
            return value;
        const auto& frame = frames[f];
        const Vector3 outgoing = outputVertices[b] - outputVertices[a];
        const double u = dot(frame.gradient[0], outgoing), v = dot(frame.gradient[1], outgoing);
        const size_t axis = std::fabs(v) > std::fabs(u);
        value.direction = frame.axis[axis];
        if (dot(value.direction, outgoing) < 0)
            value.direction = -value.direction;
        value.minimum = .5 * frame.size[axis];
        value.maximum = (faces[f].strictCurvature ? 1.2 : 2.5) * frame.size[axis];
        value.valid = std::isfinite(value.maximum);
        return value;
    };
    for (auto& item : edges) {
        auto& edge = item.second;
        const auto a = endpoint(item.first.first, item.first.second), b = endpoint(item.first.second, item.first.first);
        const Vector3 direction = a.direction - b.direction;
        if (!a.valid || !b.valid || !(direction.length() > 1e-12))
            continue;
        edge.direction = direction.normalized();
        edge.minimum = .5 * (a.minimum + b.minimum);
        edge.maximum = .5 * (a.maximum + b.maximum);
        edge.valid = std::isfinite(edge.minimum + edge.maximum);
    }
    for (size_t i = 0; i < result.size(); ++i) {
        if (links[i].size() != 4 || mapped[i] == none)
            continue;
        std::map<size_t, std::vector<size_t>> ring;
        for (const auto& link : links[i]) {
            ring[link.first].push_back(link.second);
            ring[link.second].push_back(link.first);
        }
        bool valid = ring.size() == 4;
        for (const auto& link : ring) {
            const auto& edge = edges.at(std::minmax(i, link.first));
            valid &= link.second.size() == 2 && link.second[0] != link.second[1] && edge.uses == 2 && edge.valid;
        }
        if (!valid)
            continue;
        auto& stencil = result[i];
        stencil.neighbors[0] = ring.begin()->first;
        for (size_t k = 1; k < 4; ++k) {
            const auto& next = ring.at(stencil.neighbors[k - 1]);
            stencil.neighbors[k] = k == 1 ? std::min(next[0], next[1]) : next[0] == stencil.neighbors[k - 2] ? next[1]
                                                                                                             : next[0];
            valid &= std::find(stencil.neighbors.begin(), stencil.neighbors.begin() + k, stencil.neighbors[k]) == stencil.neighbors.begin() + k;
        }
        const auto& last = ring.at(stencil.neighbors[3]);
        valid &= last[0] == stencil.neighbors[0] || last[1] == stencil.neighbors[0];
        if (!valid)
            continue;
        for (size_t k = 0; k < 4; ++k) {
            const auto& edge = edges.at(std::minmax(i, stencil.neighbors[k]));
            stencil.outward[k] = i < stencil.neighbors[k] ? edge.direction : -edge.direction;
            stencil.minimum[k] = edge.minimum;
            stencil.maximum[k] = edge.maximum;
        }
        const auto& face = faces[mapped[i]];
        const double anisotropy = face.strictCurvature ? std::max(0., std::min(1., face.anisotropy)) : 0;
        stencil.blend = anisotropy < .2 ? 0 : anisotropy > .85 ? 1
                                                               : anisotropy;
        stencil.valid = true;
    }
    return result;
}
}
