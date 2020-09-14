/*********************************************************************
 *  Author  : Himangshu Saikia, Wiebke Koepp, Anke Friederici
 *  Init    : Monday, September 11, 2017 - 12:58:42
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <labmarchingsquares/marchingsquares.h>
#include <inviwo/core/util/utilities.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo MarchingSquares::processorInfo_{
    "org.inviwo.MarchingSquares",  // Class identifier
    "Marching Squares",            // Display name
    "KTH Lab",                     // Category
    CodeState::Experimental,       // Code state
    Tags::None,                    // Tags
};

const ProcessorInfo MarchingSquares::getProcessorInfo() const { return processorInfo_; }

MarchingSquares::MarchingSquares()
    : Processor()
    , inData("volumeIn")
    , meshIsoOut("meshIsoOut")
    , meshGridOut("meshGridOut")
    , propShowGrid("showGrid", "Show Grid")
    , propGridColor("gridColor", "Grid Lines Color", vec4(0.0f, 0.0f, 0.0f, 1.0f), vec4(0.0f),
                    vec4(1.0f), vec4(0.1f), InvalidationLevel::InvalidOutput,
                    PropertySemantics::Color)
    , propDeciderType("deciderType", "Decider Type")
    , propRandomSeed("seed", "Random Seed", 0, 0, std::mt19937::max())
    , propMultiple("multiple", "Iso Levels")
    , propIsoValue("isovalue", "Iso Value")
    , propIsoColor("isoColor", "Color", vec4(0.0f, 0.0f, 1.0f, 1.0f), vec4(0.0f), vec4(1.0f),
                   vec4(0.1f), InvalidationLevel::InvalidOutput, PropertySemantics::Color)
    , propNumContours("numContours", "Number of Contours", 1, 1, 50, 1)
    , propIsoTransferFunc("isoTransferFunc", "Colors", &inData) {
    // Register ports
    addPort(inData);
    addPort(meshIsoOut);
    addPort(meshGridOut);

    // Register properties
    addProperty(propShowGrid);
    addProperty(propGridColor);

    addProperty(propDeciderType);
    propDeciderType.addOption("asymptotic", "Asymptotic", 0);
    propDeciderType.addOption("random", "Random", 1);

    addProperty(propRandomSeed);
    propRandomSeed.setSemantics(PropertySemantics::Text);

    addProperty(propMultiple);

    propMultiple.addOption("single", "Single", 0);
    addProperty(propIsoValue);
    addProperty(propIsoColor);

    propMultiple.addOption("multiple", "Multiple", 1);
    addProperty(propNumContours);
    addProperty(propIsoTransferFunc);

    // The default transfer function has just two blue points
    propIsoTransferFunc.get().clear();
    propIsoTransferFunc.get().add(0.0f, vec4(0.0f, 0.0f, 1.0f, 1.0f));
    propIsoTransferFunc.get().add(1.0f, vec4(0.0f, 0.0f, 1.0f, 1.0f));
    propIsoTransferFunc.setCurrentStateAsDefault();

    util::hide(propGridColor, propRandomSeed, propNumContours, propIsoTransferFunc);

    propDeciderType.onChange([this]() {
        if (propDeciderType.get() == 1) {
            util::show(propRandomSeed);
        } else {
            util::hide(propRandomSeed);
        }
    });

    // Show the grid color property only if grid is actually displayed
    propShowGrid.onChange([this]() {
        if (propShowGrid.get()) {
            util::show(propGridColor);
        } else {
            util::hide(propGridColor);
        }
    });

    // Show options based on display of one or multiple iso contours
    propMultiple.onChange([this]() {
        if (propMultiple.get() == 0) {
            util::show(propIsoValue, propIsoColor);
            util::hide(propNumContours, propIsoTransferFunc);
        } else {
            util::hide(propIsoValue);
            util::show(propIsoColor, propNumContours);

            // TODO (Bonus): Comment out above if you are using the transfer function
            // and comment in below instead
            // util::hide(propIsoValue, propIsoColor);
            // util::show(propNumContours, propIsoTransferFunc);
        }
    });
}

void MarchingSquares::process() {
    if (!inData.hasData()) {
        return;
    }

    // Create a structured grid from the input volume
    auto vol = inData.getData();
    auto grid = ScalarField2::createFieldFromVolume(vol);

    // Extract the minimum and maximum value from the input data
    const double minValue = grid.getMinValue();
    const double maxValue = grid.getMaxValue();

    // Set the range for the isovalue to that minimum and maximum
    propIsoValue.setMinValue(minValue);
    propIsoValue.setMaxValue(maxValue);

    // You can print to the Inviwo console with Log-commands:
    LogProcessorInfo("This scalar field contains values between " << minValue << " and " << maxValue
                                                                  << ".");
    // You can also inform about errors and warnings:
    // LogProcessorWarn("I am warning about something"); // Will print warning message in yellow
    // LogProcessorError("I am letting you know about an error"); // Will print error message in red
    // (There is also LogNetwork...() and just Log...(), these display a different source,
    // LogProcessor...() for example displays the name of the processor in the workspace while
    // Log...() displays the identifier of the processor (thus with multiple processors of the
    // same kind you would not know which one the information is coming from

    // Get the definition of our structured grid with
    // - number of vertices in each dimension {nx, ny}
    const ivec2 nVertPerDim = grid.getNumVerticesPerDim();
    // - bounding box {xmin, ymin} - {xmax, ymax}
    const dvec2 bBoxMin = grid.getBBoxMin();
    const dvec2 bBoxMax = grid.getBBoxMax();
    // - cell size {dx, dy}
    const dvec2 cellSize = grid.getCellSize();

    // Values at the vertex positions can be accessed by the indices of the vertex
    // with index i ranging between [0, nx-1] and j in [0, ny-1]
    ivec2 ij = {0, 0};
    double valueAt00 = grid.getValueAtVertex(ij);
    LogProcessorInfo("The value at (0,0) is: " << valueAt00 << ".");

    // Initialize the output: mesh and vertices for the grid and bounding box
    auto gridmesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> gridvertices;

    auto indexBufferBBox = gridmesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);
    // bottomLeft to topLeft
    drawLineSegment(bBoxMin, vec2(bBoxMin[0], bBoxMax[1]), propGridColor.get(),
                    indexBufferBBox.get(), gridvertices);
    // topLeft to topRight
    drawLineSegment(vec2(bBoxMin[0], bBoxMax[1]), bBoxMax, propGridColor.get(),
                    indexBufferBBox.get(), gridvertices);
    // topRight to bottomRight
    drawLineSegment(bBoxMax, vec2(bBoxMax[0], bBoxMin[1]), propGridColor.get(),
                    indexBufferBBox.get(), gridvertices);
    // bottomRight to bottomLeft
    drawLineSegment(vec2(bBoxMax[0], bBoxMin[1]), bBoxMin, propGridColor.get(),
                    indexBufferBBox.get(), gridvertices);

    // Set the random seed to the one selected in the interface
    randGenerator.seed(static_cast<std::mt19937::result_type>(propRandomSeed.get()));
    // You can create a random sample between min and max with
    float minRand = 0.0;
    float maxRand = 1.0;
    float rand = randomValue(minRand, maxRand);
    LogProcessorInfo("The first random sample for seed " << propRandomSeed.get() << " between "
                                                         << minRand << " and " << maxRand << " is "
                                                         << rand << ".");

    // Properties are accessed with propertyName.get()
    if (propShowGrid.get()) {
        // The function drawLineSegments creates two vertices at the specified positions,
        // that are placed into the Vertex vector defining our mesh.
        // An index buffer specifies which of those vertices should be grouped into to make up
        // lines/trianges/quads. Here two vertices make up a line segment.
        auto indexBufferGrid = gridmesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);

        // ** TASK 3.1.A BEGIN **

        // Draw vertical grid lines
        float step = bBoxMin.x;
        for (size_t i = 0; i < grid.getNumVerticesPerDim().x; i++) {
            vec2 v1 = vec2(step, bBoxMin.y);
            vec2 v2 = vec2(step, bBoxMax.y);
            drawLineSegment(v1, v2, propGridColor.get(), indexBufferGrid.get(), gridvertices);

            step += grid.getCellSize().x;
        }

        // Draw horizontal grid lines
        step = bBoxMin.y;
        for (size_t i = 0; i < grid.getNumVerticesPerDim().y; i++) {
            // Draw a line segment from v1 to v2 with a the given color for the grid
            vec2 v1 = vec2(bBoxMin.x, step);
            vec2 v2 = vec2(bBoxMax.x, step);
            drawLineSegment(v1, v2, propGridColor.get(), indexBufferGrid.get(), gridvertices);

            step += grid.getCellSize().y;
        }

        // ** TASK 3.1.A END **
    }

    // Set the created grid mesh as output
    // (this was pushed further down in the function)
    // gridmesh->addVertices(gridvertices);
    // meshGridOut.setData(gridmesh);

    // TODO (Bonus) Gaussian filter
    // Our input is const (i.e. cannot be altered), but you need to compute smoothed data and write
    // it somewhere
    // Create an editable structured grid with ScalarField2 smoothedField =
    // ScalarField2(nVertPerDim, bBoxMin, bBoxMax - bBoxMin); Values can be set with
    // smoothedField.setValueAtVertex({0, 0}, 4.2);
    // and read again in the same way as before
    // smoothedField.getValueAtVertex(ij);
    // Initialize the output: mesh and vertices
    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;

    if (propMultiple.get() == 0) {

        // ** TASK 3.1.B BEGIN **

        auto indexBufferIso = gridmesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);

        int cellCountX = grid.getNumVerticesPerDim().x - 1;
        int cellCountY = grid.getNumVerticesPerDim().y - 1;

        // March through all grid cells
        for (size_t i = 0; i < cellCountX; i++) {
            for (size_t j = 0; j < cellCountY; j++) {
                // Get values from grid points surrounding this cell
                auto f00 = grid.getValueAtVertex({i, j});
                auto f10 = grid.getValueAtVertex({i + 1, j});
                auto f01 = grid.getValueAtVertex({i, j + 1});
                auto f11 = grid.getValueAtVertex({i + 1, j + 1});

                // Mark vertices with true for above the isovalue and false for below the isovalue
                bool b00 = f00 >= propIsoValue;
                bool b10 = f10 >= propIsoValue;
                bool b01 = f01 >= propIsoValue;
                bool b11 = f11 >= propIsoValue;

                // Disregard this cell if all values have the same signal
                if (b00 == b10 && b01 == b11 && b00 == b11) continue;

                // Otherwise, proceed with finding the intersection points
                int intersectionCount = 0;
                vec2 intersections[4];

                // We will need coordinates for each grid point
                vec2 v00 = grid.getPositionAtVertex({i, j});
                vec2 v10 = grid.getPositionAtVertex({i + 1, j});
                vec2 v01 = grid.getPositionAtVertex({i, j + 1});
                vec2 v11 = grid.getPositionAtVertex({i + 1, j + 1});

                // We verify each line and add an intersection point whenever one is found
                if (b00 != b10) {
                    // We calculate the relative position along the line by dividing the slope at
                    // the isovalue by the total slope of the line This value 't' will be in the
                    // range [0-1]
                    float t = (propIsoValue - f00) / (f10 - f00);

                    // Then we interpolate between the coordinates of the two points using 't'
                    // (Note: interpolating in the Y axis here is redundant since both points are in
                    // the same Y coordinate)
                    intersections[intersectionCount].x = ((1 - t) * v00.x) + (t * v10.x);
                    intersections[intersectionCount].y = v00.y;

                    // Increment the intersection count
                    intersectionCount++;
                }
                // The same snippet is repeated for the second line...
                if (b01 != b11) {
                    float t = (propIsoValue - f01) / (f11 - f01);
                    intersections[intersectionCount].x = ((1 - t) * v01.x) + (t * v11.x);
                    intersections[intersectionCount].y = v01.y;
                    intersectionCount++;
                }
                // ...and the third line...
                if (b00 != b01) {
                    float t = (propIsoValue - f00) / (f01 - f00);
                    intersections[intersectionCount].x = v00.x;
                    intersections[intersectionCount].y = ((1 - t) * v00.y) + (t * v01.y);
                    intersectionCount++;
                }
                // ...and the last line
                if (b10 != b11) {
                    float t = (propIsoValue - f10) / (f11 - f10);
                    intersections[intersectionCount].x = v10.x;
                    intersections[intersectionCount].y = ((1 - t) * v10.y) + (t * v11.y);
                    intersectionCount++;
                }

                // Now we should have either 2 or 4 intersection points, if we have two then we can
                // simply connect them
                if (intersectionCount == 2) {
                    drawLineSegment(intersections[0], intersections[1], propIsoColor.get(),
                                    indexBufferIso.get(), gridvertices);
                    continue;
                }

                // If we have 4 intersections, we need to run a decider
                if (intersectionCount == 4) {
                    switch (propDeciderType) {
                        case 0: {
                            // Asymptotic Decider
                            // We will be using the "sort by axis" technique
                            // Bubble sort on X axis coordinate
                            bool isSorted = false;
                            while (!isSorted) {
                                isSorted = true;
                                for (size_t k = 0; k < 3; k++) {
                                    if (intersections[k].x > intersections[k + 1].x) {
                                        isSorted = false;
                                        vec2 v1 = intersections[k];
                                        vec2 v2 = intersections[k + 1];
                                        intersections[k] = v2;
                                        intersections[k + 1] = v1;
                                    }
                                }
                            }

                            // Connect 0-1 and 2-3
                            drawLineSegment(intersections[0], intersections[1], propIsoColor.get(),
                                            indexBufferIso.get(), gridvertices);
                            drawLineSegment(intersections[2], intersections[3], propIsoColor.get(),
                                            indexBufferIso.get(), gridvertices);
						} 
						break;
                            
						case 1: {
							// ** TASK 3.1.C BEGIN **

							// Random decider 
							// NOTE: Because of the order the lines are evaluated, we avoid the 0-1, 2-3 combination to prevent the lines
							// from crossing each other
							// Valid combinations are: 0-2, 1-3 and 0-3, 1-2
							// We give it a 50/50 chance of either
							
							vec2 combo1 = {0, 2};     
							vec2 combo2 = {1, 3};  

							if (randomValue(0, 1) > 0.5f) {
								combo1 = {0, 3};
								combo2 = {1, 2};                    
							}

							// Connect lines
                            drawLineSegment(intersections[(int)combo1.x], intersections[(int)combo1.y], propIsoColor.get(),
                                            indexBufferIso.get(), gridvertices);
                            drawLineSegment(intersections[(int)combo2.x], intersections[(int)combo2.y], propIsoColor.get(),
                                            indexBufferIso.get(), gridvertices);

							// ** TASK 3.1.C END **
						}
                    }

                    
                    continue;
                }
            }
        }
        // ** TASK 3.1.B END **
    }

    else {
        // TODO: Draw the given number (propNumContours) of isolines between
        // the minimum and maximum value

        // TODO (Bonus): Use the transfer function property to assign a color
        // The transfer function normalizes the input data and sampling colors
        // from the transfer function assumes normalized input, that means
        // vec4 color = propIsoTransferFunc.get().sample(0.0f);
        // is the color for the minimum value in the data
        // vec4 color = propIsoTransferFunc.get().sample(1.0f);
        // is the color for the maximum value in the data
    }

    gridmesh->addVertices(gridvertices);
    meshGridOut.setData(gridmesh);

    // Note: It is possible to add multiple index buffers to the same mesh,
    // thus you could for example add one for the grid lines and one for
    // each isoline
    // Also, consider to write helper functions to avoid code duplication
    // e.g. for the computation of a single iso contour

    mesh->addVertices(vertices);
    meshIsoOut.setData(mesh);
}

float MarchingSquares::randomValue(const float min, const float max) const {
    return min + uniformReal(randGenerator) * (max - min);
}

void MarchingSquares::drawLineSegment(const vec2& v1, const vec2& v2, const vec4& color,
                                      IndexBufferRAM* indexBuffer,
                                      std::vector<BasicMesh::Vertex>& vertices) {
    // Add first vertex
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    // A vertex has a position, a normal, a texture coordinate and a color
    // we do not use normal or texture coordinate, but still have to specify them
    vertices.push_back({vec3(v1[0], v1[1], 0), vec3(0, 0, 1), vec3(v1[0], v1[1], 0), color});
    // Add second vertex
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(v2[0], v2[1], 0), vec3(0, 0, 1), vec3(v2[0], v2[1], 0), color});
}

}  // namespace inviwo