/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Tuesday, September 19, 2017 - 15:08:33
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
#include <inviwo/core/ports/meshport.h>
#include <inviwo/core/ports/volumeport.h>
#include <inviwo/core/processors/processor.h>
#include <inviwo/core/properties/boolproperty.h>
#include <inviwo/core/properties/compositeproperty.h>
#include <inviwo/core/properties/eventproperty.h>
#include <inviwo/core/properties/optionproperty.h>
#include <inviwo/core/properties/ordinalproperty.h>
#include <inviwo/core/properties/cameraproperty.h>
#include <labstreamlines/labstreamlinesmoduledefine.h>
#include <labutils/scalarvectorfield.h>

namespace inviwo {

/** \docpage{org.inviwo.StreamlineIntegrator, Streamline Integrator}
    ![](org.inviwo.StreamlineIntegrator.png?classIdentifier=org.inviwo.StreamlineIntegrator)

    Processor to integrate streamlines.

    ### Inports
    * __data__ The input here is 2-dimensional vector field (with vectors of
    two components thus two values within each voxel) but it is represented
    by a 3-dimensional volume.
    This processor deals with 2-dimensional data only, therefore it is assumed
    the z-dimension will have size 1 otherwise the 0th slice of the volume
    will be processed.

    ### Outports
    * __meshout__ The output mesh contains linesegments making up either a single or
    multiple stream lines
    * __meshBBoxOut__ Mesh with boundling box

    ### Properties
    * __propSeedMode__ Mode for the number of seeds, either a single start point
   or multiple
    * __propStartPoint__ Location of the start point
    * __mouseMoveStart__ Move the start point when a selected mouse button is
    * __numStepsTaken__ Number of steps actually taken for a single streamline
   pressed (default left)
*/

class IVW_MODULE_LABSTREAMLINES_API StreamlineIntegrator : public Processor {
    // Construction / Deconstruction
public:
    StreamlineIntegrator();
    virtual ~StreamlineIntegrator() = default;

    // Methods
public:
    virtual const ProcessorInfo getProcessorInfo() const override;
    static const ProcessorInfo processorInfo_;

protected:
    /// Our main computation function
    virtual void process() override;

    /// Function to handle mouse interaction for a single streamline
    void eventMoveStart(Event* event);

    // (TODO: You could define some helper functions here,
    // e.g. a function creating a single streamline from one seed point)

    //TASK 4.3 (a) returns a random point from the domain
    vec2 getRandomPoint(vec2 min, vec2 max);
    void EulerLoop(const VectorField2& vectorField, const dvec2& start,
                   std::shared_ptr<inviwo::BasicMesh>& mesh,
                   std::vector<BasicMesh::Vertex>& vertices, bool inverted = false);
    bool Euler(const VectorField2& vectorField, const dvec2& start, dvec2& end, double& arcLength,
               bool inverted = false);

	void RK4Loop(const VectorField2& vectorField, const dvec2& start,
                   std::shared_ptr<inviwo::BasicMesh>& mesh,
                   std::vector<BasicMesh::Vertex>& vertices, bool inverted = false);
    bool RK4(const VectorField2& vectorField, const dvec2& start, dvec2& end, double& arcLength, 
			    bool inverted = false);

	double Magnitude(const dvec2& vec);

    // Ports
public:
    // Input Vector Field
    VolumeInport inData;

    // Output mesh
    MeshOutport meshOut;

    // Output mesh for bounding box and gridlines
    MeshOutport meshBBoxOut;

    // Properties
public:

    FloatVec2Property propStartPoint;
    TemplateOptionProperty<int> propSeedMode;

    IntProperty propNumStepsTaken;
    EventProperty mouseMoveStart;

    // Declare additional properties
    FloatVec4Property propEulerColor;
    FloatVec4Property propRK4Color;
    BoolProperty propShowPoints;
    BoolProperty propForwardDirection;
    BoolProperty propBackwardDirection;
    BoolProperty propNormalizedField;
    IntProperty propNumberSteps;
    FloatProperty propStepSize;
    FloatProperty propMinimumVelocity;
    FloatProperty propMaximumArcLength;

    IntProperty propNumStreamLines;

    // Attributes
private:
    dvec2 BBoxMin_{0, 0};
    dvec2 BBoxMax_{0, 0};
};

}  // namespace inviwo
