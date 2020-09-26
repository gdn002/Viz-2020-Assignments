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
    , propRandom("random", "Random Placement", true)
    , propNumStreamLines("numStreamLines", "Stream Lines", 1, 0, 100, 1)
    , propGridPointX("gridPointX", "X", 1, 0, 100, 1)
    , propGridPointY("gridPointY", "Y", 1, 0, 100, 1) {
    
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
    addProperty(propRandom);
    addProperty(propNumStreamLines);
    addProperty(propGridPointX);
    addProperty(propGridPointY);

    // Show properties for a single seed and hide properties for multiple seeds
    propSeedMode.onChange([this]() {
        if (propSeedMode.get() == 0) {
            util::show(propStartPoint, mouseMoveStart, propNumStepsTaken);
            util::hide(propRandom, propNumStreamLines, propGridPointX, propGridPointY);
        } else {
            util::hide(propStartPoint, mouseMoveStart, propNumStepsTaken);
            util::show(propRandom, propNumStreamLines, propGridPointX, propGridPointY);
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
    int steps = propNumberSteps.get(), num_steps, msg;
    dvec2 startingPoint = propStartPoint.get();
    vec2 point;
    if (propSeedMode.get() == 0) {
        if (propForwardDirection.get()) {
            LogProcessorInfo(Integrator::EulerLoop(
                vectorField, point, mesh, vertices, num_steps, propStepSize.get(),
                propMinimumVelocity.get(), propMaximumArcLength.get(), propNormalizedField.get(),
                propEulerColor.get(), propNumberSteps.get(), propShowPoints.get()));
            LogProcessorInfo(Integrator::RK4Loop(
                vectorField, point, mesh, vertices, num_steps, propStepSize.get(),
                propMinimumVelocity.get(), propMaximumArcLength.get(), propNormalizedField.get(),
                propRK4Color.get(), propNumberSteps.get(), propShowPoints.get()));
        }
        if (propBackwardDirection.get()) {
            LogProcessorInfo(Integrator::EulerLoop(
                vectorField, point, mesh, vertices, num_steps, propStepSize.get(),
                propMinimumVelocity.get(), propMaximumArcLength.get(), propNormalizedField.get(),
                propEulerColor.get(), propNumberSteps.get(), propShowPoints.get(), true));
            LogProcessorInfo(Integrator::RK4Loop(
                vectorField, point, mesh, vertices, num_steps, propStepSize.get(),
                propMinimumVelocity.get(), propMaximumArcLength.get(), propNormalizedField.get(),
                propRK4Color.get(), propNumberSteps.get(), propShowPoints.get(), true));
        }
    } // TASK 4.3 IMPLEMENTATION BEGIN
    else{
      if(propRandom.get()){
        for (int i = 0; i < propNumStreamLines.get(); i++) {
            point = getRandomPoint(BBoxMin_, BBoxMax_);
            if (propForwardDirection.get()) {
                num_steps = Integrator::EulerLoop(vectorField, point, mesh, vertices);
                num_steps = Integrator::RK4Loop(vectorField, point, mesh, vertices);
            }
            if (propBackwardDirection.get()) {
                num_steps = Integrator::EulerLoop(vectorField, point, mesh, vertices, true);
                num_steps = Integrator::RK4Loop(vectorField, point, mesh, vertices, true);
            }
        }
      }
      else{
        int nX = propGridPointX.get();
        int nY = propGridPointY.get();
        float unitX = (BBoxMax_.x - BBoxMin_.x)/nX;
        float unitY = (BBoxMax_.y - BBoxMin_.y)/nY;
        for (float j = BBoxMin_.y; j<BBoxMax_.y; j+= unitY){
          for(float i = BBoxMin_.x; i<BBoxMax_.x; i+= unitX){
            point = {i,j};
            if (propForwardDirection.get()) {
                num_steps = Integrator::EulerLoop(vectorField, point, mesh, vertices);
                num_steps = Integrator::RK4Loop(vectorField, point, mesh, vertices);
            }
            if (propBackwardDirection.get()) {
                num_steps = Integrator::EulerLoop(vectorField, point, mesh, vertices, true);
                num_steps = Integrator::RK4Loop(vectorField, point, mesh, vertices, true);
            }
          }
        }
      }
    } // TASK 4.3 IMPLEMENTATION END

    // Use the propNumStepsTaken property to show how many steps have actually been
    // integrated This could be different from the desired number of steps due to stopping
    // conditions (too slow, boundary, ...)
    propNumStepsTaken.set(num_steps);

    mesh->addVertices(vertices);
    meshOut.setData(mesh);
}  // process









}  // namespace inviwo
