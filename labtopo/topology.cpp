/*********************************************************************
 *  Author  : Anke Friederici
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 **********************************************************************/

#include <inviwo/core/datastructures/geometry/basicmesh.h>
#include <inviwo/core/datastructures/volume/volumeram.h>
#include <labstreamlines/integrator.h>
#include <labutils/scalarvectorfield.h>
#include <labtopo/topology.h>
#include <labtopo/utils/gradients.h>

namespace inviwo {

const vec4 Topology::ColorsCP[6] = {
    vec4(1, 1, 0, 1),    // Saddle - Yellow
    vec4(1, 0, 0, 1),    // AttractingNode - Red
    vec4(0, 0, 1, 1),    // RepellingNode - Blue
    vec4(0.5, 0, 1, 1),  // AttractingFocus - Purple
    vec4(1, 0.5, 0, 1),  // RepellingFocus - Orange
    vec4(0, 1, 0, 1)     // Center - Green
};

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo Topology::processorInfo_{
    "org.inviwo.Topology",    // Class identifier
    "Vector Field Topology",  // Display name
    "KTH Lab",                // Category
    CodeState::Experimental,  // Code state
    Tags::None,               // Tags
};

const ProcessorInfo Topology::getProcessorInfo() const { return processorInfo_; }

Topology::Topology()
    : Processor()
    , inData("inData")
    , outMesh("meshOut")
    , meshBBoxOut("meshBBoxOut")
// TODO: Initialize additional properties
// propertyName("propertyIdentifier", "Display Name of the Propery",
// default value (optional), minimum value (optional), maximum value (optional), increment
// (optional)); propertyIdentifier cannot have spaces
    , propMinLength("minimumLength", "Minimum Cell Length", 1.0f, 0.000001f, 5, 0.01f)
	, propDrawSeparatrices("drawSeparatrices", "Draw Separatrices", false)
	, propMaxSeparatrices("maxSeparatrices", "Maximum Separatrice Cycles", 100, 10, 1000)
	, propSeparatriceOffset("separatriceOffset", "Separatrice Offset", 0.001, 0.001, 1)
{
    // Register Ports
    addPort(outMesh);
    addPort(inData);
    addPort(meshBBoxOut);

    // TODO: Register additional properties
    // addProperty(propertyName);
    addProperty(propMinLength);
    addProperty(propDrawSeparatrices);
    addProperty(propMaxSeparatrices);
    addProperty(propSeparatriceOffset);
}

void Topology::process() {
    // Get input
    if (!inData.hasData()) {
        return;
    }
    auto vol = inData.getData();

    // Retreive data in a form that we can access it
    const VectorField2 vectorField = VectorField2::createFieldFromVolume(vol);

    // Add a bounding box to the mesh
    const dvec2& BBoxMin = vectorField.getBBoxMin();
    const dvec2& BBoxMax = vectorField.getBBoxMax();
    auto bboxMesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> bboxVertices;
    auto indexBufferBBox = bboxMesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
    // Bounding Box vertex 0
    vec4 black = vec4(0, 0, 0, 1);
    Integrator::drawNextPointInPolyline(BBoxMin, black, indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(vec2(BBoxMin[0], BBoxMax[1]), black, indexBufferBBox.get(),
                                        bboxVertices);
    Integrator::drawNextPointInPolyline(BBoxMax, black, indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(vec2(BBoxMax[0], BBoxMin[1]), black, indexBufferBBox.get(),
                                        bboxVertices);
    // Connect back to the first point, to make a full rectangle
    indexBufferBBox->add(static_cast<std::uint32_t>(0));
    bboxMesh->addVertices(bboxVertices);
    meshBBoxOut.setData(bboxMesh);

    // Initialize mesh, vertices and index buffers for seperatrices
    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;
    // Either add all line segments to this index buffer (one large buffer, two consecutive points
    // make up one line), or use several index buffers with connectivity type strip.
    auto indexBufferSeparatrices = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);
    // auto indexBufferSeparatrices = mesh->addIndexBuffer(DrawType::Lines,
    // ConnectivityType::Strip);

    auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);

    // TODO: Compute the topological skeleton of the input vector field.
    // Find the critical points and color them according to their type.
    // Integrate all separatrices.

    size2_t dims = vectorField.getNumVerticesPerDim();

    // Looping through all cells in the vector field.
    dvec2 p_00, p_01, p_10, p_11;
    for (size_t j = 0; j < dims[1]-1; ++j) {
        for (size_t i = 0; i < dims[0]-1; ++i) {
            p_00 = vectorField.getPositionAtVertex(size2_t(i, j));
            p_10 = vectorField.getPositionAtVertex(size2_t(i+1, j));
            p_01 = vectorField.getPositionAtVertex(size2_t(i, j+1));
            p_11 = vectorField.getPositionAtVertex(size2_t(i+1, j+1));

            checkChangeOfSign(vectorField, indexBufferPoints.get(),
                              indexBufferSeparatrices.get(), vertices, p_00, p_10, p_01, p_11,
                              p_10.x - p_00.x, p_01.y - p_00.y);
        }
    }

    // Loop through horizontal (X-axis) boundary lines
    dvec2 p0, p1;
    for (size_t i = 0; i < dims[0]-1; ++i) {
        // Find boundary switch points on the lower line
        p0 = vectorField.getPositionAtVertex(size2_t(i, 0));
        p1 = vectorField.getPositionAtVertex(size2_t(i+1, 0));
        checkChangeOfSignX(vectorField, indexBufferPoints.get(),
                           indexBufferSeparatrices.get(), vertices, p0, p1, p1.x - p0.x);

       // Find boundary switch points on the upper line
       p0 = vectorField.getPositionAtVertex(size2_t(i, dims[1]-1));
       p1 = vectorField.getPositionAtVertex(size2_t(i+1, dims[1]-1));
       checkChangeOfSignX(vectorField, indexBufferPoints.get(),
                          indexBufferSeparatrices.get(), vertices, p0, p1, p1.x - p0.x);

    }

    // Loop through vertical (Y-axis) boundary lines
    for (size_t j = 0; j < dims[1]-1; ++j) {
        // Find boundary switch points on the left line
        p0 = vectorField.getPositionAtVertex(size2_t(0, j));
        p1 = vectorField.getPositionAtVertex(size2_t(0, j+1));
        checkChangeOfSignY(vectorField, indexBufferPoints.get(),
                           indexBufferSeparatrices.get(), vertices, p0, p1, p1.y - p0.y);

       // Find boundary switch points on the right line
       p0 = vectorField.getPositionAtVertex(size2_t(dims[0]-1, j));
       p1 = vectorField.getPositionAtVertex(size2_t(dims[0]-1, j+1));
       checkChangeOfSignY(vectorField, indexBufferPoints.get(),
                          indexBufferSeparatrices.get(), vertices, p0, p1, p1.y - p0.y);

    }

    // Other helpful functions
    // dvec2 pos = vectorField.getPositionAtVertex(size2_t(i, j));
    // Computing the jacobian at a position
    // dmat2 jacobian = vectorField.derive(pos);
    // Doing the eigen analysis
    // auto eigenResult = util::eigenAnalysis(jacobian);
    // The result of the eigen analysis has attributed eigenvaluesRe eigenvaluesIm and
    // eigenvectors

    // Accessing the colors
    vec4 colorCenter = ColorsCP[static_cast<int>(TypeCP::Center)];

    mesh->addVertices(vertices);
    outMesh.setData(mesh);
}

void Topology::checkChangeOfSign(const VectorField2& vectorField,
    IndexBufferRAM* indexBufferPoints, IndexBufferRAM* indexBufferLines, std::vector<BasicMesh::Vertex>& vertices, dvec2
    pos00,dvec2 pos10,dvec2 pos01,dvec2 pos11, float lengthX, float lengthY){
   dvec2 v_00 = vectorField.interpolate(pos00);
   dvec2 v_10 = vectorField.interpolate(pos10);
   dvec2 v_01 = vectorField.interpolate(pos01);
   dvec2 v_11 = vectorField.interpolate(pos11);

    //if(v_00.x == 0 && v_00.y == 0){
    //    // this is a critical point
    //    Integrator::drawPoint(pos00, {0,0,0,1}, indexBuffer, vertices);
    //}
    //if(v_01.x == 0 && v_01.y == 0){
    //    // this is a critical point
    //    Integrator::drawPoint(pos00, {0,0,0,1}, indexBuffer, vertices);
    //}
    //if(v_10.x == 0 && v_10.y == 0){
    //    // this is a critical point
    //    Integrator::drawPoint(pos00, {0,0,0,1}, indexBuffer, vertices);
    //}
    //if(v_11.x == 0 && v_11.y == 0){
    //    // this is a critical point
    //    Integrator::drawPoint(pos00, {0,0,0,1}, indexBuffer, vertices);
    //}
    //else
      if ((v_00.x>=0 && v_01.x>=0 && v_10.x>=0 && v_11.x>=0) || //if not all x have the same sign
              (v_00.x<0 && v_01.x<0 && v_10.x<0 && v_11.x<0) ||
          (v_00.y>=0 && v_01.y>=0 && v_10.y>=0 && v_11.y>=0) || //and not all y have the same sign
              (v_00.y<0 && v_01.y<0 && v_10.y<0 && v_11.y<0)){}
      else{
        // this may be a critical point
        //stop condition
        if(lengthX < propMinLength.get() || lengthY < propMinLength.get()){
            auto t = analyzeCriticalPoint(vectorField, pos00);

			if (propDrawSeparatrices.get() && t == TypeCP::Saddle) {
                drawSeparatrices(vectorField, indexBufferLines, vertices, pos00);

            }
            Integrator::drawPoint(pos00, ColorsCP[(int)t], indexBufferPoints, vertices);
            return;
        }

        lengthX /= 2;
        lengthY /= 2;

       dvec2 midD = {pos00.x+lengthX, pos00.y};
       dvec2 midL = {pos00.x, pos00.y+lengthY};
       dvec2 midC = {pos00.x+lengthX, pos00.y+lengthY};
       dvec2 midR = {pos10.x, pos10.y+lengthY};
       dvec2 midU = {pos01.x+lengthX, pos01.y};
        checkChangeOfSign(vectorField,indexBufferPoints,indexBufferLines,vertices,pos00,midD,midL,midC,lengthX,lengthY);
        checkChangeOfSign(vectorField,indexBufferPoints,indexBufferLines,vertices,midD,pos10,midC,midR,lengthX,lengthY);
        checkChangeOfSign(vectorField,indexBufferPoints,indexBufferLines,vertices,midL,midC,pos01,midU,lengthX,lengthY);
        checkChangeOfSign(vectorField,indexBufferPoints,indexBufferLines,vertices,midC,midR,midU,pos11,lengthX,lengthY);

    }

}

void Topology::checkChangeOfSignX(const VectorField2& vectorField,
    IndexBufferRAM* indexBufferPoints, IndexBufferRAM* indexBufferLines,
    std::vector<BasicMesh::Vertex>& vertices, dvec2 p0, dvec2 p1, float lengthX) {
    dvec2 v0 = vectorField.interpolate(p0);
    dvec2 v1 = vectorField.interpolate(p1);

    if ((v0.y>=0 && v1.y>=0) || (v0.y<0 && v1.y<0)){}
    else {
        // This may be a boundary switch point
        if(lengthX < propMinLength.get()) {
            // Integrate and draw the separatrices (can use the switch point
            // as seed since this boundary point is not a critical point)
            if (propDrawSeparatrices.get()) {
                drawSingleSeparatrice(vectorField, indexBufferLines, vertices, p0, true);
                drawSingleSeparatrice(vectorField, indexBufferLines, vertices, p0, false);
            }

            vec4 color = {1.0, 1.0, 1.0, 1.0};
            Integrator::drawPoint(p0, color, indexBufferPoints, vertices);
            return;
        }

        lengthX /= 2;
        dvec2 midPoint = {p0.x+lengthX, p0.y};

        checkChangeOfSignX(vectorField, indexBufferPoints, indexBufferLines,
                           vertices, p0, midPoint, lengthX);
        checkChangeOfSignX(vectorField, indexBufferPoints, indexBufferLines,
                           vertices, midPoint, p1, lengthX);
    }
}

void Topology::checkChangeOfSignY(const VectorField2& vectorField,
    IndexBufferRAM* indexBufferPoints, IndexBufferRAM* indexBufferLines,
    std::vector<BasicMesh::Vertex>& vertices, dvec2 p0, dvec2 p1, float lengthY) {
    dvec2 v0 = vectorField.interpolate(p0);
    dvec2 v1 = vectorField.interpolate(p1);

    if ((v0.x>=0 && v1.x>=0) || (v0.x<0 && v1.x<0)){}
    else {
        // This may be a boundary switch point
        if(lengthY < propMinLength.get()) {
            // Integrate and draw the separatrices (can use the switch point
            // as seed since this boundary point is not a critical point)
            if (propDrawSeparatrices.get()) {
                drawSingleSeparatrice(vectorField, indexBufferLines, vertices, p0, true);
                drawSingleSeparatrice(vectorField, indexBufferLines, vertices, p0, false);
            }

            vec4 color = {1.0, 1.0, 1.0, 1.0};
            Integrator::drawPoint(p0, color, indexBufferPoints, vertices);
            return;
        }

        lengthY /= 2;
        dvec2 midPoint = {p0.x, p0.y+lengthY};

        checkChangeOfSignY(vectorField, indexBufferPoints, indexBufferLines,
                           vertices, p0, midPoint, lengthY);
        checkChangeOfSignY(vectorField, indexBufferPoints, indexBufferLines,
                           vertices, midPoint, p1, lengthY);
    }
}

void Topology::drawLineSegment(const dvec2& v1, const dvec2& v2, const vec4& color,
                               IndexBufferRAM* indexBuffer,
                               std::vector<BasicMesh::Vertex>& vertices) {
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(v1[0], v1[1], 0), vec3(0, 0, 1), vec3(v1[0], v1[1], 0), color});
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(v2[0], v2[1], 0), vec3(0, 0, 1), vec3(v2[0], v2[1], 0), color});
}

Topology::TypeCP inviwo::Topology::analyzeCriticalPoint(const VectorField2& vectorField,
                                                        const dvec2& criticalPoint) {

    // The critical point MUST be an isolated critical point
    auto eigenValues = util::eigenAnalysis(vectorField.derive(criticalPoint));

    // Imaginary values are zero (Saddle, Attracting Node, Repelling Node)
    if (eigenValues.eigenvaluesIm[0] == 0 && eigenValues.eigenvaluesIm[1] == 0) {

        // Both real values are greater than zero (Repelling Node)
        if (eigenValues.eigenvaluesRe[0] > 0 && eigenValues.eigenvaluesRe[1] > 0)
            return TypeCP::RepellingNode;

        // Both real values are lesser than zero (Attracting Node)
        if (eigenValues.eigenvaluesRe[0] < 0 && eigenValues.eigenvaluesRe[1] < 0)
            return TypeCP::AttractingNode;

        // One real value is greater than zero and the other is lesser than zero (Saddle)
        return TypeCP::Saddle;
    }
    // Imaginary values are nonzero (Center, Attracting Focus, Repelling Focus)
    else {

        // Both real values are greater than zero (Repelling Focus)
        if (eigenValues.eigenvaluesRe[0] > 0 && eigenValues.eigenvaluesRe[1] > 0)
            return TypeCP::RepellingFocus;

        // Both real values are lesser than zero (Attracting Focus)
        if (eigenValues.eigenvaluesRe[0] < 0 && eigenValues.eigenvaluesRe[1] < 0)
            return TypeCP::AttractingFocus;

        // Both real values are zero (Center)
        return TypeCP::Center;
    }
}

void inviwo::Topology::getSeedingPointsFromSaddle(const VectorField2& vectorField,
                                                  const dvec2& saddlePoint, dvec2 (&out)[2],
                                                  dvec2 (&in)[2]) {

    auto eigenValues = util::eigenAnalysis(vectorField.derive(saddlePoint));

    // The initial step from the saddle point at which the line is seeded
    double firstStep = propSeparatriceOffset.get();


    // Check for outbound seeding points (positive EigenValue)
    // These must be integrated FORWARD
    for (size_t i = 0; i < 2; i++) {
        if (eigenValues.eigenvaluesRe[i] > 0) {
            // Two points, one in each direction
            out[0] = Integrator::Add(saddlePoint, (eigenValues.eigenvectors[i] * firstStep));
            out[1] = Integrator::Add(saddlePoint, (eigenValues.eigenvectors[i] * -firstStep));
            break;
        }
    }

    // Check for inbound seeding points (negative EigenValue)
    // These must be integrated BACKWARD
    for (size_t i = 0; i < 2; i++) {
        if (eigenValues.eigenvaluesRe[i] < 0) {
            // Two points, one in each direction
            in[0] = Integrator::Add(saddlePoint, (eigenValues.eigenvectors[i] * firstStep));
            in[1] = Integrator::Add(saddlePoint, (eigenValues.eigenvectors[i] * -firstStep));
            break;
        }
    }
}

void inviwo::Topology::drawSeparatrices(const VectorField2& vectorField,
                                        IndexBufferRAM* indexBuffer,
                                        std::vector<BasicMesh::Vertex>& vertices,
                                        const dvec2& saddlePoint) {
    dvec2 inbound[2];
    dvec2 outbound[2];

	getSeedingPointsFromSaddle(vectorField, saddlePoint, outbound, inbound);

	for (size_t i = 0; i < 2; i++) {
		drawSingleSeparatrice(vectorField, indexBuffer, vertices, inbound[i], true);
		drawSingleSeparatrice(vectorField, indexBuffer, vertices, outbound[i], false);
    }
}

void inviwo::Topology::drawSingleSeparatrice(const VectorField2& vectorField,
                                             IndexBufferRAM* indexBuffer,
                                             std::vector<BasicMesh::Vertex>& vertices,
                                             const dvec2& seedingPoint, bool inverted) {
    dvec2 current = seedingPoint;
    dvec2 next;

	// TODO: Make this a property
	double stepSize = 0.1;

	// Integrate line until it hits a sink or a boundary
    for (size_t i = 0; i < propMaxSeparatrices.get(); i++) {
		if (!Integrator::RK4Lite(vectorField, current, next, stepSize, false, inverted)) return;
        Integrator::drawLineSegment(current, next, {1, 1, 1, 1}, indexBuffer, vertices);
        current = next;
    }
}

}  // namespace inviwo
