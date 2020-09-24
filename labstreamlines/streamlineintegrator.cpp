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
    , propStartPoint("startPoint", "Start Point", vec2(0.5f, 0.5f), vec2(-1.f), vec2(1.f), vec2(0.1))
    , propSeedMode("seedMode", "Seeds")
    , propNumStepsTaken("numstepstaken", "Number of actual steps", 0, 0, 100000)
    , propLineColor("lineColor", "Stream Line Color", vec4(1.0f, 0.0f, 0.0f, 1.0f), vec4(0.0f), vec4(1.0f),
					vec4(0.1f), InvalidationLevel::InvalidOutput, PropertySemantics::Color)
    , propShowPoints("showPoints", "Show Points")
    , propNumberSteps("numberSteps", "Number of Steps", 5, 1, 50, 1)
    , propStepSize("stepSize", "Step Size", 1.0f, 0, 5.0f, 0.1f)
    , propForwardDirection("forwardDirection", "Forward Direction", true)
    , propNormalizedField("normalizedField", "Normalized Vector Field")
    , propMinimumVelocity("minimumVelocity", "Minimum Velocity", 0.1f, 0, 5, 0.1f)
    , propMaximumArcLength("maximumArcLength", "Maximum Arc Length", 5, 0, 10, 0.5f)
    , mouseMoveStart("mouseMoveStart", "Move Start", [this](Event* e) { eventMoveStart(e); },
                     MouseButton::Left, MouseState::Press | MouseState::Move)
// TODO: Initialize additional properties
// propertyName("propertyIdentifier", "Display Name of the Propery",
// default value (optional), minimum value (optional), maximum value (optional),
// increment (optional)); propertyIdentifier cannot have spaces
    , propNumStreamLines("numStreamLines", "Stream Lines", 1, 0, 10, 1)
{
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
    addProperty(propLineColor);
    addProperty(propShowPoints);
    addProperty(propNumberSteps);
    addProperty(propStepSize);
    addProperty(propForwardDirection);
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

//TASK 4.x FUNCTIONS IMPLEMENTATION START
vec2 StreamlineIntegrator::getRandomPoint(){
  float r;
  vec2 point;
  for(int i=0; i<2; i++){
    r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);  //random between 0.0 and 1.0
    r = r*2-1; //remap to interval [-1;1] 
    point[i] = r;
  }
  return point;
}

int drawStreamLine(dvec2 startPoint, double stepSize, int nSteps, float minVel, bool forward, bool normalize, const VectorField2& vectorField, IndexBufferRAM* indexBufferLine, std::vector<BasicMesh::Vertex>& vertices){
                                   
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

        Integrator::drawLineSegment(current, next, propLineColor.get(), indexBufferLine.get(), vertices);
        if (propShowPoints.get()) Integrator::drawPoint(next, propLineColor.get(), indexBufferPoints.get(), vertices);
        current = next;
    }
    return num_steps;
}
//TASK 4.x FUNCTIONS IMPLEMENTATION END

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
    if (propSeedMode.get() == 0) {
        auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);
        auto indexBufferLine = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
        vec2 startPoint = propStartPoint.get();
        // Draw start point
        Integrator::drawPoint(startPoint, vec4(0, 0, 0, 1), indexBufferPoints.get(), vertices);

        // Create one stream line from the given start point
        int num_steps = drawStreamLine(startPoint, stepSize, steps, vectorField, indexBufferLine, vertices);

        // Use the propNumStepsTaken property to show how many steps have actually been
        // integrated This could be different from the desired number of steps due to stopping
        // conditions (too slow, boundary, ...)
        propNumStepsTaken.set(num_steps);

    } else {
        // TODO: Seed multiple stream lines either randomly or using a uniform grid
        // (TODO: Bonus, sample randomly according to magnitude of the vector field)
        // TASK 4.3 IMPLEMENTATION BEGIN

        // Draw start point
     
        auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);
        vec2 point;
        for(int i=0; i<propNumStreamLines.get(); i++){
          point = getRandomPoint();
          drawStreamLine(startPoint, stepSize, steps, vectorField, indexBufferLine, vertices);
        }
        // TASK 4.3 IMPLEMENTATION END
    }

    mesh->addVertices(vertices);
    meshOut.setData(mesh);
}  // process

}  // namespace inviwo
