/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Tuesday, September 19, 2017 - 15:08:33
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <inviwo/core/interaction/events/mouseevent.h>
#include <inviwo/core/util/utilities.h>
#include <labstreamlines/integrator.h>
#include <labstreamlines/streamlineintegrator.h>
#include <labutils/scalarvectorfield.h>
#include <stdlib.h>
#include <time.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming
// scheme
const ProcessorInfo StreamlineIntegrator::processorInfo_{
    "org.inviwo.StreamlineIntegrator",  // Class identifier
    "Streamline Integrator",            // Display name
    "KTH Lab",                          // Category
    CodeState::Experimental,            // Code state
    Tags::None,                         // Tags
};

const ProcessorInfo StreamlineIntegrator::getProcessorInfo() const { return processorInfo_; }

StreamlineIntegrator::StreamlineIntegrator()
    : Processor()
    , inData("volIn")
    , meshOut("meshOut")
    , meshBBoxOut("meshBBoxOut")
    , propStartPoint("startPoint", "Start Point", vec2(0.5f, 0.5f), vec2(-1.f), vec2(1.f),
                     vec2(0.1))
    , propSeedMode("seedMode", "Seeds")
    , propNumStepsTaken("numstepstaken", "Number of actual steps", 0, 0, 100000)
    , propEulerColor("eulerColor", "Euler Stream Line Color", vec4(1.0f, 0.0f, 0.0f, 1.0f),
                     vec4(0.0f), vec4(1.0f), vec4(0.1f), InvalidationLevel::InvalidOutput,
                     PropertySemantics::Color)
    , propRK4Color("rk4Color", "RK4 Stream Line Color", vec4(0.0f, 0.0f, 1.0f, 1.0f), vec4(0.0f),
                   vec4(1.0f), vec4(0.1f), InvalidationLevel::InvalidOutput,
                   PropertySemantics::Color)
    , propShowPoints("showPoints", "Show Points")
    , propNumberSteps("numberSteps", "Number of Steps", 5, 1, 50, 1)
    , propStepSize("stepSize", "Step Size", 1.0f, 0, 5.0f, 0.1f)
    , propForwardDirection("forwardDirection", "Forward Direction", true)
    , propBackwardDirection("backwardDirection", "Backward Direction", false)
    , propNormalizedField("normalizedField", "Normalized Vector Field")
    , propMinimumVelocity("minimumVelocity", "Minimum Velocity", 0.1f, 0, 5, 0.1f)
    , propMaximumArcLength("maximumArcLength", "Maximum Arc Length", 5, 0, 20, 0.5f)
    , mouseMoveStart(
          "mouseMoveStart", "Move Start", [this](Event* e) { eventMoveStart(e); },
          MouseButton::Left, MouseState::Press | MouseState::Move)
    // TODO: Initialize additional properties
    // propertyName("propertyIdentifier", "Display Name of the Propery",
    // default value (optional), minimum value (optional), maximum value (optional),
    // increment (optional)); propertyIdentifier cannot have spaces
    , propNumStreamLines("numStreamLines", "Stream Lines", 1, 0, 10, 1) {
    // Register Ports
    addPort(inData);
    addPort(meshOut);
    addPort(meshBBoxOut);

    // Register Properties
    propSeedMode.addOption("one", "Single Start Point", 0);
    propSeedMode.addOption("multiple", "Multiple Seeds", 1);
    addProperty(propSeedMode);
    addProperty(propStartPoint);
    addProperty(propNumStepsTaken);
    propNumStepsTaken.setReadOnly(true);
    propNumStepsTaken.setSemantics(PropertySemantics::Text);
    addProperty(mouseMoveStart);
    addProperty(propNumStreamLines);

    // Register additional properties
    addProperty(propEulerColor);
    addProperty(propRK4Color);
    addProperty(propShowPoints);
    addProperty(propNumberSteps);
    addProperty(propStepSize);
    addProperty(propForwardDirection);
    addProperty(propBackwardDirection);
    addProperty(propNormalizedField);
    addProperty(propMinimumVelocity);
    addProperty(propMaximumArcLength);

    // Show properties for a single seed and hide properties for multiple seeds
    // (TODO)
    propSeedMode.onChange([this]() {
        if (propSeedMode.get() == 0) {
            util::show(propStartPoint, mouseMoveStart, propNumStepsTaken);
            util::hide(propNumStreamLines);
        } else {
            util::hide(propStartPoint, mouseMoveStart, propNumStepsTaken);
            util::show(propNumStreamLines);
        }
    });
}

// TASK 4.x FUNCTIONS IMPLEMENTATION START
vec2 StreamlineIntegrator::getRandomPoint(vec2 min, vec2 max) {
    float r;
    vec2 point;
    for (int i = 0; i < 2; i++) {
        r = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);  // random between 0.0
                                                                        // and 1.0
        r = r * (max[i] - min[i]) + min[i];
        point[i] = r;
    }
    return point;
}

/*int drawStreamLine(dvec2 startPoint, double stepSize, int nSteps, float minVel, bool forward, bool
normalize, const VectorField2& vectorField, IndexBufferRAM* indexBufferLine,
std::vector<BasicMesh::Vertex>& vertices){
                                   
    dvec2 BBoxMin_ = vectorField.getBBoxMin();
    dvec2 BBoxMax_ = vectorField.getBBoxMax();
    dvec2 current = startPoint;
    double arcLength = 0;
    int num_steps;
    for (num_steps = 0; num_steps < nSteps; num_steps++) {
        // Interpolate vector field from our current point
        dvec2 vecValue = vectorField.interpolate(current);

        // Compute the euclidean norm and stop integration on very low / zero velocity
        double vecNorm = sqrt(pow(vecValue.x, 2) + pow(vecValue.y, 2));
        if (vecNorm == 0) {
            //LogProcessorInfo("Stopping integration due to zeros in the vector field.");
            break;
        } else if (vecNorm < minVel) {
            //LogProcessorInfo("Stopping integration due to slow velocity.");
            break;
        }
        // Invert the direction of the vector
        if (!forward) vecValue *= -1;

        // obtain the normalized vector field
        if (normalize) {
            vecValue /= vecNorm;
            vecNorm = 1.0;
        }
        // Keep track of the accumulating arc length
        arcLength += vecNorm;
        dvec2 next = current + stepSize * vecValue;

        if (next.x < BBoxMin_.x || next.x > BBoxMax_.x ||
            next.y < BBoxMin_.y || next.y > BBoxMax_.y) {
            //LogProcessorInfo("Stopping integration at the domain's boundary.");
            break;
        } else if (arcLength > propMaximumArcLength.get()) {
            //LogProcessorInfo("Stopping integration due to exceeded arc length.");
            break;
        }

        Integrator::drawLineSegment(current, next, propLineColor.get(), indexBufferLine.get(),
vertices); if (propShowPoints.get()) Integrator::drawPoint(next, propLineColor.get(),
indexBufferPoints.get(), vertices); current = next;
    }
    return num_steps;
}*/
// TASK 4.x FUNCTIONS IMPLEMENTATION END

void StreamlineIntegrator::eventMoveStart(Event* event) {
    if (!inData.hasData()) return;
    auto mouseEvent = static_cast<MouseEvent*>(event);
    vec2 mousePos = mouseEvent->posNormalized();

    // Map to bounding box range
    mousePos[0] *= static_cast<float>(BBoxMax_[0] - BBoxMin_[0]);
    mousePos[1] *= static_cast<float>(BBoxMax_[1] - BBoxMin_[1]);
    mousePos += static_cast<vec2>(BBoxMin_);

    // Update starting point
    propStartPoint.set(mousePos);
    event->markAsUsed();
}

void StreamlineIntegrator::process() {
    // Get input
    if (!inData.hasData()) {
        return;
    }
    auto vol = inData.getData();

    srand(time(NULL));

    // Retreive data in a form that we can access it
    auto vectorField = VectorField2::createFieldFromVolume(vol);
    BBoxMin_ = vectorField.getBBoxMin();
    BBoxMax_ = vectorField.getBBoxMax();

    auto bboxMesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> bboxVertices;

    // Make bounding box without vertex duplication, instead of line segments which duplicate
    // vertices, create line segments between each added points with connectivity type of the index
    // buffer
    auto indexBufferBBox = bboxMesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
    // Bounding Box vertex 0
    vec4 black = vec4(0, 0, 0, 1);
    Integrator::drawNextPointInPolyline(BBoxMin_, black, indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(vec2(BBoxMin_[0], BBoxMax_[1]), black,
                                        indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(BBoxMax_, black, indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(vec2(BBoxMax_[0], BBoxMin_[1]), black,
                                        indexBufferBBox.get(), bboxVertices);
    // Connect back to the first point, to make a full rectangle
    indexBufferBBox->add(static_cast<std::uint32_t>(0));
    bboxMesh->addVertices(bboxVertices);
    meshBBoxOut.setData(bboxMesh);

    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;

    double stepSize = propStepSize.get();
    int steps = propNumberSteps.get();
    bool inverted = propBackwardDirection.get();
    if (propSeedMode.get() == 0) {

        dvec2 startingPoint = propStartPoint.get();

        if (propForwardDirection.get()) {
            EulerLoop(vectorField, startingPoint, mesh, vertices);
            RK4Loop(vectorField, startingPoint, mesh, vertices);
		}
        if (propBackwardDirection.get()) {
            EulerLoop(vectorField, startingPoint, mesh, vertices, true);
            RK4Loop(vectorField, startingPoint, mesh, vertices, true);
        }

    } else {
        // TODO: Seed multiple stream lines either randomly or using a uniform grid
        // (TODO: Bonus, sample randomly according to magnitude of the vector field)
        // TASK 4.3 IMPLEMENTATION BEGIN

        // Draw start point

        auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);
        vec2 point;
        for (int i = 0; i < propNumStreamLines.get(); i++) {
            point = getRandomPoint(BBoxMin_, BBoxMax_);
            EulerLoop(vectorField, point, mesh, vertices, inverted);
        }
        // TASK 4.3 IMPLEMENTATION END
    }

    mesh->addVertices(vertices);
    meshOut.setData(mesh);
}  // process

void inviwo::StreamlineIntegrator::EulerLoop(const VectorField2& vectorField, const dvec2& start,
                                             std::shared_ptr<inviwo::BasicMesh>& mesh,
                                             std::vector<BasicMesh::Vertex>& vertices,
                                             bool inverted) {
    // Bypass entire function if the Euler color alpha is zero
    if (propEulerColor.get().a == 0) return;

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
        if (!Euler(vectorField, current, next, arcLength, inverted)) break;

        Integrator::drawLineSegment(current, next, propEulerColor.get(), indexBufferLine.get(),
                                    vertices);
        if (propShowPoints.get())
            Integrator::drawPoint(next, propEulerColor.get(), indexBufferPoints.get(), vertices);
        current = next;
    }

    // Use the propNumStepsTaken property to show how many steps have actually been
    // integrated This could be different from the desired number of steps due to stopping
    // conditions (too slow, boundary, ...)
    propNumStepsTaken.set(num_steps);
}

bool inviwo::StreamlineIntegrator::Euler(const VectorField2& vectorField, const dvec2& start,
                                         dvec2& end, double& arcLength, bool inverted) {

    // Interpolate vector field from our current point
    dvec2 vecValue = vectorField.interpolate(start);

    // Compute the euclidean norm and stop integration on very low / zero velocity
    double vecNorm = Magnitude(vecValue);
    if (vecNorm == 0) {
        LogProcessorInfo("Stopping integration due to zeros in the vector field.");
        return false;
    } else if (vecNorm < propMinimumVelocity.get()) {
        LogProcessorInfo("Stopping integration due to slow velocity.");
        return false;
    }
    // Invert the direction of the vector
    if (inverted) vecValue *= -1;

    // obtain the normalized vector field
    if (propNormalizedField.get()) {
        vecValue /= vecNorm;
        vecNorm = 1.0;
    }
    // Keep track of the accumulating arc length
    arcLength += vecNorm;
    end = start + ((double)propStepSize.get() * vecValue);

    if (end.x < BBoxMin_.x || end.x > BBoxMax_.x || end.y < BBoxMin_.y || end.y > BBoxMax_.y) {
        LogProcessorInfo("Stopping integration at the domain's boundary.");
        return false;
    } else if (arcLength > propMaximumArcLength.get()) {
        LogProcessorInfo("Stopping integration due to exceeded arc length.");
        return false;
    }

    return true;
}

void inviwo::StreamlineIntegrator::RK4Loop(const VectorField2& vectorField, const dvec2& start,
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
    propNumStepsTaken.set(num_steps);
}

bool inviwo::StreamlineIntegrator::RK4(const VectorField2& vectorField, const dvec2& start,
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

double inviwo::StreamlineIntegrator::Magnitude(const dvec2& vec) {
    return sqrt(pow(vec.x, 2) + pow(vec.y, 2));
}

}  // namespace inviwo
