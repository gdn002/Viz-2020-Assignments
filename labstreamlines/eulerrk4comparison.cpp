/*********************************************************************
 *  Author  : Himangshu Saikia, Wiebke Koepp
 *  Init    : Tuesday, September 19, 2017 - 15:08:24
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <inviwo/core/datastructures/geometry/basicmesh.h>
#include <inviwo/core/interaction/events/mouseevent.h>
#include <labstreamlines/eulerrk4comparison.h>
#include <labstreamlines/integrator.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo EulerRK4Comparison::processorInfo_{
    "org.inviwo.EulerRK4Comparison",  // Class identifier
    "Euler RK4 Comparison",           // Display name
    "KTH Lab",                        // Category
    CodeState::Experimental,          // Code state
    Tags::None,                       // Tags
};

const ProcessorInfo EulerRK4Comparison::getProcessorInfo() const { return processorInfo_; }

EulerRK4Comparison::EulerRK4Comparison()
    : Processor()
    , inData("inData")
    , meshOut("meshOut")
    , meshBBoxOut("meshBBoxOut")
    , propStartPoint("startPoint", "Start Point", vec2(0.5, 0.5), vec2(0), vec2(1024), vec2(0.5))
    , propEulerColor("eulerColor", "Euler Color", vec4(1.0f, 0.0f, 0.0f, 1.0f), vec4(0.0f), vec4(1.0f), 
					vec4(0.1f), InvalidationLevel::InvalidOutput, PropertySemantics::Color)
    , propEulerShowPoints("eulerShowPoints", "Show Euler Points")
    , propEulerNumberSteps("eulerNumberSteps", "Euler Number of Steps", 10, 1, 500, 1)
    , propEulerStepSize("eulerStepSize", "Euler Step Size", 1.0f, 0, 5.0f, 0.1f)
    , propRK4Color("rk4Color", "RK4 Color", vec4(0.0f, 0.0f, 1.0f, 1.0f), vec4(0.0f), vec4(1.0f),
					vec4(0.1f), InvalidationLevel::InvalidOutput, PropertySemantics::Color)
    , propRK4ShowPoints("rk4ShowPoints", "Show RK4 Points")
    , propRK4NumberSteps("rk4NumberSteps", "RK4 Number of Steps", 10, 1, 500, 1)
    , propRK4StepSize("rk4StepSize", "RK4 Step Size", 1.0f, 0, 5.0f, 0.1f)
    , mouseMoveStart(
          "mouseMoveStart", "Move Start", [this](Event* e) { eventMoveStart(e); },
          MouseButton::Left, MouseState::Press | MouseState::Move)
// TODO: Initialize additional properties
// propertyName("propertyIdentifier", "Display Name of the Propery",
// default value (optional), minimum value (optional), maximum value (optional), increment
// (optional)); propertyIdentifier cannot have spaces
{
    // Register Ports
    addPort(meshOut);
    addPort(meshBBoxOut);
    addPort(inData);

    // Register Properties
    addProperty(propStartPoint);

    addProperty(propEulerColor);
    addProperty(propEulerShowPoints);
    addProperty(propEulerNumberSteps);
    addProperty(propEulerStepSize);

    addProperty(propRK4Color);
    addProperty(propRK4ShowPoints);
    addProperty(propRK4NumberSteps);
    addProperty(propRK4StepSize);

    addProperty(mouseMoveStart);

    // TODO: Register additional properties
    // addProperty(propertyName);
}

void EulerRK4Comparison::eventMoveStart(Event* event) {
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

void EulerRK4Comparison::process() {
    // Get input
    if (!inData.hasData()) {
        return;
    }
    auto vol = inData.getData();

    // Retreive data in a form that we can access it
    const VectorField2 vectorField = VectorField2::createFieldFromVolume(vol);
    BBoxMin_ = vectorField.getBBoxMin();
    BBoxMax_ = vectorField.getBBoxMax();

    // The start point should be inside the volume (set maximum to the upper right corner)
    propStartPoint.setMinValue(BBoxMin_ - dvec2(1, 1));
    propStartPoint.setMaxValue(BBoxMax_ + dvec2(1, 1));

    // Initialize mesh, vertices and index buffers for the two streamlines and the points
    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;

    auto indexBufferEuler = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
    auto indexBufferRK = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
    auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);
    auto indexBufferEulerPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);
    auto indexBufferRKPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);

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

    // Draw start point
    dvec2 startPoint = propStartPoint.get();
    Integrator::drawPoint(startPoint, black, indexBufferPoints.get(), vertices);

	// TASK 4.1 IMPLEMENTATION BEGIN

	// Run Euler integration
    dvec2 current = startPoint;
	for (int i = 0; i < propEulerNumberSteps.get(); i++) {
        dvec2 next = Integrator::Euler(vectorField, current, propEulerStepSize.get());
        Integrator::drawLineSegment(current, next, propEulerColor.get(), indexBufferEuler.get(), vertices);
        if (propEulerShowPoints) Integrator::drawPoint(next, propEulerColor.get(), indexBufferEulerPoints.get(), vertices);
        current = next;
    }

	// Run RK4 integration
	current = startPoint;
	for (int i = 0; i < propRK4NumberSteps.get(); i++) {
        dvec2 next = Integrator::RK4(vectorField, current, propRK4StepSize.get());
        Integrator::drawLineSegment(current, next, propRK4Color.get(), indexBufferRK.get(), vertices);
        if (propRK4ShowPoints) Integrator::drawPoint(next, propRK4Color.get(), indexBufferRKPoints.get(), vertices);
        current = next;
    }

	// TASK 4.1 IMPLEMENTATION END

    mesh->addVertices(vertices);
    meshOut.setData(mesh);
}

}  // namespace inviwo
