/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Wednesday, September 20, 2017 - 12:04:15
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <labstreamlines/integrator.h>

namespace inviwo {

void Integrator::drawPoint(const dvec2& p, const vec4& color, IndexBufferRAM* indexBuffer,
                           std::vector<BasicMesh::Vertex>& vertices) {
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(p[0], p[1], 0), vec3(0, 0, 1), vec3(p[0], p[1], 0), color});
}

// Alias for draw point
void Integrator::drawNextPointInPolyline(const dvec2& p, const vec4& color,
                                         IndexBufferRAM* indexBuffer,
                                         std::vector<BasicMesh::Vertex>& vertices) {
    Integrator::drawPoint(p, color, indexBuffer, vertices);
}

void Integrator::drawLineSegment(const dvec2& v1, const dvec2& v2, const vec4& color,
                                 IndexBufferRAM* indexBuffer,
                                 std::vector<BasicMesh::Vertex>& vertices) {
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(v1[0], v1[1], 0), vec3(0, 0, 1), vec3(v1[0], v1[1], 0), color});
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(v2[0], v2[1], 0), vec3(0, 0, 1), vec3(v2[0], v2[1], 0), color});
}

// TASK 4.1 IMPLEMENTATION BEGIN

dvec2 Integrator::RK4(const VectorField2& vectorField, const dvec2& position, const float& stepSize) {
    // Obtain the four intermediate points
	dvec2 v1 = vectorField.interpolate(position);
    dvec2 v2 = vectorField.interpolate(position + Multiply(v1, stepSize / 2));
    dvec2 v3 = vectorField.interpolate(position + Multiply(v2, stepSize / 2));
    dvec2 v4 = vectorField.interpolate(position + Multiply(v3, stepSize));
    
	// Divide the intermediate points according to the RK4 formula
    v1 = Divide(v1, 6);
    v2 = Divide(v2, 3);
    v3 = Divide(v3, 3);
    v4 = Divide(v4, 6);

	return position + Multiply((v1 + v2 + v3 + v4), stepSize);
}

dvec2 Integrator::Euler(const VectorField2& vectorField, const dvec2& position, const float& stepSize) {
    dvec2 value = vectorField.interpolate(position);
    return position + Multiply(value, stepSize);
}

dvec2 Integrator::Multiply(const dvec2& vector, const float& factor) { 
	dvec2 v = vector;
    v.x *= factor;
    v.y *= factor;
    return v;
}

dvec2 Integrator::Divide(const dvec2& vector, const float& divisor) {
    dvec2 v = vector;
    v.x /= divisor;
    v.y /= divisor;
    return v;
}

// TASK 4.1 IMPLEMENTATION END

} // namespace inviwo
