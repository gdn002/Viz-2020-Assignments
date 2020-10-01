/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Wednesday, September 20, 2017 - 12:04:15
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#pragma once

#include <inviwo/core/common/inviwo.h>
#include <inviwo/core/datastructures/geometry/basicmesh.h>
#include <inviwo/core/datastructures/volume/volumeram.h>
#include <labstreamlines/labstreamlinesmoduledefine.h>
#include <labutils/scalarvectorfield.h>
#include <string>

namespace inviwo {

class IVW_MODULE_LABSTREAMLINES_API Integrator {

public:
    // Construction / Deconstruction
public:
    Integrator() {}
    virtual ~Integrator() = default;

	inline static std::atomic_ulong RK4CallCounter;

    // Methods
public:
    // Add a point to a mesh
    static void drawPoint(const dvec2& p, const vec4& color, IndexBufferRAM* indexBuffer,
                          std::vector<BasicMesh::Vertex>& vertices);
    // Add a point to a polyline, assumes that the indexBuffer uses Strip Connectivity
    static void drawNextPointInPolyline(const dvec2& p, const vec4& color,
                                        IndexBufferRAM* indexBuffer,
                                        std::vector<BasicMesh::Vertex>& vertices);
    // Add a line segment to a mesh, assuming no connectivity in the indexBuffer
    static void drawLineSegment(const dvec2& v1, const dvec2& v2, const vec4& color,
                                IndexBufferRAM* indexBuffer,
                                std::vector<BasicMesh::Vertex>& vertices);
    // TODO: Implement the methods below (one integration step with either Euler or
    // Runge-Kutte of 4th order integration method)
    // Pass any other properties that influence the integration process
    // Examples would be the stepsize, inegreation direction, ...
    static dvec2 RK4(const VectorField2& vectorField, const dvec2& position, const float& stepSize);
    static dvec2 Euler(const VectorField2& vectorField, const dvec2& position,
                       const float& stepSize);

    static std::string EulerLine(const VectorField2& vectorField, const dvec2& start, dvec2& end,
								double& arcLength, double stepSize, float minVelocity,
								float maxArchLength, bool normalize, bool inverted = false);
    static std::string RK4line(const VectorField2& vectorField, const dvec2& start, dvec2& end,
                               double& arcLength, double stepSize, float minVelocity,
                               float maxArchLength, bool normalize, bool inverted = false);

    static std::string EulerLoop(const VectorField2& vectorField, const dvec2& start,
								std::shared_ptr<inviwo::BasicMesh>& mesh,
								std::vector<BasicMesh::Vertex>& vertices, int& stepsTaken, float stepSize,
								float minVelocity, float maxArchLength, bool normalize,
								const vec4& color = {0, 0, 0, 255}, int steps = 1, bool showSteps = false,
								bool inverted = false);
    static std::string RK4Loop(const VectorField2& vectorField, const dvec2& start,
                               std::shared_ptr<inviwo::BasicMesh>& mesh,
                               std::vector<BasicMesh::Vertex>& vertices, int& stepsTaken,
                               float stepSize, float minVelocity, float maxArchLength,
                               bool normalize, const vec4& color = {0, 0, 0, 255}, int steps = 1,
                               bool showSteps = false, bool inverted = false);
                              
	static bool RK4Lite(const VectorField2& vectorField, const dvec2& start, dvec2& end,
                               double stepSize, bool normalize, bool inverted);

    static dvec2 Multiply(const dvec2& vector, const float& factor);
    static dvec2 Divide(const dvec2& vector, const float& divisor);

    static double Magnitude(const dvec2& vec);
};

}  // namespace inviwo
