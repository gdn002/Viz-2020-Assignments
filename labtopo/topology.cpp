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
    , propMinLength("minimumLength", "Minimum Cell Length", 1.0f, 0, 5, 0.01f)
{
    // Register Ports
    addPort(outMesh);
    addPort(inData);
    addPort(meshBBoxOut);

    // TODO: Register additional properties
    // addProperty(propertyName);
    addProperty(propMinLength);
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
            p_01 = vectorField.getPositionAtVertex(size2_t(i, j+1));
            p_10 = vectorField.getPositionAtVertex(size2_t(i+1, j));
            p_11 = vectorField.getPositionAtVertex(size2_t(i+1, j+1));
            
            checkChangeOfSign(vectorField, indexBufferPoints.get(), vertices, p_00, p_01, p_10, p_11, p_10.x-p_00.x, p_01.y-p_00.y, propMinLength.get());
        }
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
    IndexBufferRAM* indexBuffer, std::vector<BasicMesh::Vertex>& vertices, dvec2
    pos00,dvec2 pos01,dvec2 pos10,dvec2 pos11, float lengthX, float lengthY,
    float minLength){
   dvec2 v_00 = vectorField.interpolate(pos00);
   dvec2 v_01 = vectorField.interpolate(pos01);
   dvec2 v_10 = vectorField.interpolate(pos10);
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
      if ((!(v_00.x>0 && v_01.x>0 && v_10.x>0 && v_11.x>0) && //if not all x have the same sign
              !(v_00.x<0 && v_01.x<0 && v_10.x<0 && v_11.x<0)) &&
            (!(v_00.y>0 && v_01.y>0 && v_10.y>0 && v_11.y>0) && //and not all y have the same sign
             !(v_00.y<0 && v_01.y<0 && v_10.y<0 && v_11.y<0))){
        // this may be a critical point
        //stop condition
        if(lengthX < minLength || lengthY < minLength){
            Integrator::drawPoint(pos00, {0,0,0,1}, indexBuffer, vertices);
            return;
        }

        lengthX /= 2;
        lengthY /= 2;
    
       dvec2 midD = {pos00.x+lengthX, pos00.y};
       dvec2 midL = {pos00.x, pos00.y+lengthY};
       dvec2 midC = {pos00.x+lengthX, pos00.y+lengthY};
       dvec2 midR = {pos10.x, pos10.y+lengthY};
       dvec2 midU = {pos01.x+lengthX, pos01.y};
        checkChangeOfSign(vectorField,indexBuffer,vertices,pos00,midD,midL,midC,lengthX,lengthY,minLength);
        checkChangeOfSign(vectorField,indexBuffer,vertices,midD,pos01,midC,midR,lengthX,lengthY,minLength);
        checkChangeOfSign(vectorField,indexBuffer,vertices,midL,midC,pos10,midU,lengthX,lengthY,minLength);
        checkChangeOfSign(vectorField,indexBuffer,vertices,midC,midR,midU,pos11,lengthX,lengthY,minLength);
        
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

}  // namespace inviwo
