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

dvec2 Integrator::Euler(const VectorField2& vectorField, const dvec2& position,
                        const float& stepSize) {
    dvec2 value = vectorField.interpolate(position);
    return position + Multiply(value, stepSize);
}

dvec2 Integrator::RK4(const VectorField2& vectorField, const dvec2& position,
                      const float& stepSize) {
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

std::string Integrator::EulerLine(const VectorField2& vectorField, const dvec2& start, dvec2& end,
                                  double& arcLength, double stepSize, float minVelocity,
                                  float maxArchLength, bool normalize, bool inverted) {

    // Interpolate vector field from our current point
    dvec2 vecValue = vectorField.interpolate(start);
    dvec2 bbmin = vectorField.getBBoxMin();
    dvec2 bbmax = vectorField.getBBoxMax();

    // Compute the euclidean norm and stop integration on very low / zero velocity
    double vecNorm = Magnitude(vecValue);
    if (vecNorm == 0) {
        return "Stopping integration due to zeros in the vector field.";
    } else if (vecNorm < minVelocity) {
        return "Stopping integration due to slow velocity.";
    }
    // Invert the direction of the vector
    if (inverted) vecValue *= -1;

    // obtain the normalized vector field
    if (normalize) {
        vecValue /= vecNorm;
        vecNorm = 1.0;
    }
    // Keep track of the accumulating arc length
    arcLength += vecNorm;
    end = start + ((double)stepSize * vecValue);

    if (end.x < bbmin.x || end.x > bbmax.x || end.y < bbmin.y || end.y > bbmax.y) {
        return "Stopping integration at the domain's boundary.";
    } else if (arcLength > maxArchLength) {
        return "Stopping integration due to exceeded arc length.";
    }

    return "";
}

std::string Integrator::RK4line(const VectorField2& vectorField, const dvec2& start, dvec2& end,
                                double& arcLength, double stepSize, float minVelocity,
                                float maxArchLength, bool normalize, bool inverted) {
    dvec2 bbmin = vectorField.getBBoxMin();
    dvec2 bbmax = vectorField.getBBoxMax();

    // Obtain the four intermediate points
    dvec2 v1 = vectorField.interpolate(start);
    if (Magnitude(v1) == 0) return "Stopping integration due to zeros in the vector field.";
    if (inverted) v1 *= -1;
    if (normalize) v1 /= Magnitude(v1);
    dvec2 v2 = vectorField.interpolate(start + (v1 * (stepSize / 2)));
    if (Magnitude(v2) == 0) return "Stopping integration due to zeros in the vector field.";
    if (inverted) v2 *= -1;
    if (normalize) v2 /= Magnitude(v2);
    dvec2 v3 = vectorField.interpolate(start + (v2 * (stepSize / 2)));
    if (Magnitude(v3) == 0) return "Stopping integration due to zeros in the vector field.";
    if (inverted) v3 *= -1;
    if (normalize) v3 /= Magnitude(v3);
    dvec2 v4 = vectorField.interpolate(start + (v3 * stepSize));
    if (Magnitude(v4) == 0) return "Stopping integration due to zeros in the vector field.";
    if (inverted) v4 *= -1;
    if (normalize) v4 /= Magnitude(v4);

    // Divide the intermediate points according to the RK4 formula
    v1 /= 6;
    v2 /= 3;
    v3 /= 3;
    v4 /= 6;

    // Obtain final vector
    dvec2 vecValue = v1 + v2 + v3 + v4;
    double vecNorm = Magnitude(vecValue);
    if (vecNorm < minVelocity) {
        return "Stopping integration due to slow velocity.";
    }

    // Keep track of the accumulating arc length
    arcLength += vecNorm;
    end = start + (stepSize * vecValue);

    if (end.x < bbmin.x || end.x > bbmax.x || end.y < bbmin.y || end.y > bbmax.y) {
        return "Stopping integration at the domain's boundary.";
    } else if (arcLength > maxArchLength) {
        return "Stopping integration due to exceeded arc length.";
    }

    return "";
}

std::string Integrator::EulerLoop(const VectorField2& vectorField, const dvec2& start,
                                  std::shared_ptr<inviwo::BasicMesh>& mesh,
                                  std::vector<BasicMesh::Vertex>& vertices, int& stepsTaken,
                                  float stepSize, float minVelocity, float maxArchLength,
                                  bool normalize, const vec4& color, int steps, bool showSteps,
                                  bool inverted) {
    std::string msg = "";
    // Bypass entire function if the Euler color alpha is zero
    if (color.a == 0) return "";

    auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);
    auto indexBufferLine = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
    // Draw start point
    Integrator::drawPoint(start, vec4(0, 0, 0, 1), indexBufferPoints.get(), vertices);

    // Create one stream line from the given start point
    dvec2 current = start;
    double arcLength = 0;
    for (stepsTaken = 0; stepsTaken < steps; stepsTaken++) {

        dvec2 next;
        msg = EulerLine(vectorField, current, next, arcLength, stepSize, minVelocity, maxArchLength,
                        normalize, inverted);
        if (msg != "") break;

        drawLineSegment(current, next, color, indexBufferLine.get(), vertices);
        if (showSteps) drawPoint(next, color, indexBufferPoints.get(), vertices);
        current = next;
    }
    return msg;
}

std::string Integrator::RK4Loop(const VectorField2& vectorField, const dvec2& start,
                                std::shared_ptr<inviwo::BasicMesh>& mesh,
                                std::vector<BasicMesh::Vertex>& vertices, int& stepsTaken,
                                float stepSize, float minVelocity, float maxArchLength,
                                bool normalize, const vec4& color, int steps, bool showSteps,
                                bool inverted) {
    std::string msg = "";
    // Bypass entire function if the RK4 color alpha is zero
    if (color.a == 0) return "";

    auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);
    auto indexBufferLine = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
    // Draw start point
    Integrator::drawPoint(start, vec4(0, 0, 0, 1), indexBufferPoints.get(), vertices);

    // Create one stream line from the given start point
    dvec2 current = start;
    double arcLength = 0;
    for (stepsTaken = 0; stepsTaken < steps; stepsTaken++) {

        dvec2 next;
        msg = RK4line(vectorField, current, next, arcLength, stepSize, minVelocity, maxArchLength,
                      normalize, inverted);
        if (msg != "") break;

        Integrator::drawLineSegment(current, next, color, indexBufferLine.get(), vertices);
        if (showSteps) Integrator::drawPoint(next, color, indexBufferPoints.get(), vertices);
        current = next;
    }
    return msg;
}

bool inviwo::Integrator::RK4Lite(const VectorField2& vectorField, const dvec2& start,
                                        dvec2& end, double stepSize, bool normalize, bool inverted) {
    ++RK4CallCounter;
	
	end = start;
	
	dvec2 bbmin = vectorField.getBBoxMin();
    dvec2 bbmax = vectorField.getBBoxMax();

    // Obtain the four intermediate points
    dvec2 v1 = vectorField.interpolate(start);
    if (Magnitude(v1) < 0.001) return false;
    if (inverted) v1 *= -1;
    if (normalize) v1 /= Magnitude(v1);
    dvec2 v2 = vectorField.interpolate(start + (v1 * (stepSize / 2)));
    if (Magnitude(v2) < 0.001) return false;
    if (inverted) v2 *= -1;
    if (normalize) v2 /= Magnitude(v2);
    dvec2 v3 = vectorField.interpolate(start + (v2 * (stepSize / 2)));
    if (Magnitude(v3) < 0.001) return false;
    if (inverted) v3 *= -1;
    if (normalize) v3 /= Magnitude(v3);
    dvec2 v4 = vectorField.interpolate(start + (v3 * stepSize));
    if (Magnitude(v4) < 0.001) return false;
    if (inverted) v4 *= -1;
    if (normalize) v4 /= Magnitude(v4);

    // Divide the intermediate points according to the RK4 formula
    v1 /= 6;
    v2 /= 3;
    v3 /= 3;
    v4 /= 6;

    // Obtain final vector
    dvec2 vecValue = v1 + v2 + v3 + v4;
    double vecNorm = Magnitude(vecValue);

    end = start + (stepSize * vecValue);

	if (end.x < bbmin.x || end.x > bbmax.x || end.y < bbmin.y || end.y > bbmax.y) return false;
	return true;
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

double Integrator::Magnitude(const dvec2& vec) { return sqrt(pow(vec.x, 2) + pow(vec.y, 2)); }

// TASK 4.1 IMPLEMENTATION END

}  // namespace inviwo
