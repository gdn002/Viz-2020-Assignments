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

std::string Integrator::EulerLine(const VectorField2& vectorField, const dvec2& start, dvec2& end,
                          double& arcLength, float stepSize, float minVelocity, float maxArchLength,
                          bool normalize, bool inverted) {

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

    return 0;
}

std::string Integrator::RK4line(const VectorField2& vectorField, const dvec2& start,
                                           dvec2& end, double& arcLength, bool inverted) {
    double stepSize = propStepSize.get();

    // Obtain the four intermediate points
    dvec2 v1 = vectorField.interpolate(start);
    dvec2 v2 = vectorField.interpolate(start + (v1 * (stepSize / 2)));
    dvec2 v3 = vectorField.interpolate(start + (v2 * (stepSize / 2)));
    dvec2 v4 = vectorField.interpolate(start + (v3 * stepSize));

    if (Magnitude(v1) == 0) {
        LogProcessorInfo("Stopping integration due to zeros in the vector field.");
        return false;
    }

    // Invert the direction of the vectors
    if (inverted) {
        v1 *= -1;
        v2 *= -1;
        v3 *= -1;
        v4 *= -1;
    }

    // Obtain normalized vectors
    if (propNormalizedField.get()) {
        v1 /= Magnitude(v1);
        v2 /= Magnitude(v2);
        v3 /= Magnitude(v3);
        v4 /= Magnitude(v4);
    }

    // Divide the intermediate points according to the RK4 formula
    v1 /= 6;
    v2 /= 3;
    v3 /= 3;
    v4 /= 6;

    // Obtain final vector
    dvec2 vecValue = v1 + v2 + v3 + v4;
    double vecNorm = Magnitude(vecValue);
    if (vecNorm < propMinimumVelocity.get()) {
        LogProcessorInfo("Stopping integration due to slow velocity.");
        return false;
    }

    // Keep track of the accumulating arc length
    arcLength += vecNorm;
    end = start + (stepSize * vecValue);

    if (end.x < BBoxMin_.x || end.x > BBoxMax_.x || end.y < BBoxMin_.y || end.y > BBoxMax_.y) {
        LogProcessorInfo("Stopping integration at the domain's boundary.");
        return false;
    } else if (arcLength > propMaximumArcLength.get()) {
        LogProcessorInfo("Stopping integration due to exceeded arc length.");
        return false;
    }

    return true;
}

std::string Integrator::EulerLoop(const VectorField2& vectorField, const dvec2& start,
                          std::shared_ptr<inviwo::BasicMesh>& mesh,
                          std::vector<BasicMesh::Vertex>& vertices, int& stepsTaken, float stepSize,
                          float minVelocity, float maxArchLength, bool normalize,
                          const vec4& color, int steps, bool showSteps, bool inverted) {
    // Bypass entire function if the Euler color alpha is zero
    if (color.a == 0) return;

    auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);
    auto indexBufferLine = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
    // Draw start point
    Integrator::drawPoint(start, vec4(0, 0, 0, 1), indexBufferPoints.get(), vertices);

    // Create one stream line from the given start point
    dvec2 current = start;
    double arcLength = 0;
    std::string msg;
    for (stepsTaken = 0; stepsTaken < steps; stepsTaken++) {

        dvec2 next;
        msg = EulerLine(vectorField, current, next, arcLength, stepSize, minVelocity,
                             maxArchLength, normalize, inverted);
        if (msg != 0) break;

        drawLineSegment(current, next, color, indexBufferLine.get(),
                                    vertices);
        if (showSteps)
            drawPoint(next, color, indexBufferPoints.get(), vertices);
        current = next;
    }
    return msg;
}

std::string Integrator::RK4Loop(const VectorField2& vectorField, const dvec2& start,
                                           std::shared_ptr<inviwo::BasicMesh>& mesh,
                                           std::vector<BasicMesh::Vertex>& vertices,
                                           bool inverted) {
    // Bypass entire function if the RK4 color alpha is zero
    if (propRK4Color.get().a == 0) return;

    auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);
    auto indexBufferLine = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
    // Draw start point
    Integrator::drawPoint(start, vec4(0, 0, 0, 1), indexBufferPoints.get(), vertices);

    // Create one stream line from the given start point
    dvec2 current = start;
    double arcLength = 0;
    int num_steps;
    for (num_steps = 0; num_steps < propNumberSteps.get(); num_steps++) {

        dvec2 next;
        if (!RK4(vectorField, current, next, arcLength, inverted)) break;

        Integrator::drawLineSegment(current, next, propRK4Color.get(), indexBufferLine.get(),
                                    vertices);
        if (propShowPoints.get())
            Integrator::drawPoint(next, propRK4Color.get(), indexBufferPoints.get(), vertices);
        current = next;
    }

    // Use the propNumStepsTaken property to show how many steps have actually been
    // integrated This could be different from the desired number of steps due to stopping
    // conditions (too slow, boundary, ...)
    return num_steps;
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

double Integrator::Magnitude(const dvec2& vec) {
    return sqrt(pow(vec.x, 2) + pow(vec.y, 2));
}

// TASK 4.1 IMPLEMENTATION END

} // namespace inviwo
