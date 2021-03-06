/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Monday, October 02, 2017 - 13:31:17
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#pragma once

#include <inviwo/core/common/inviwo.h>
#include <inviwo/core/datastructures/image/imageram.h>
#include <inviwo/core/ports/imageport.h>
#include <inviwo/core/ports/volumeport.h>
#include <inviwo/core/processors/processor.h>
#include <inviwo/core/properties/boolproperty.h>
#include <inviwo/core/properties/ordinalproperty.h>
#include <inviwo/core/util/colorconversion.h>
#include <lablic/lablicmoduledefine.h>
#include <labutils/scalarvectorfield.h>
#include <labutils/rgbaimage.h>

namespace inviwo {

/** \docpage{org.inviwo.LICProcessor, LICProcessor}
    ![](org.inviwo.LICProcessor.png?classIdentifier=org.inviwo.LICProcessor)

    Line Integral Convolution with a box kernel.

    ### Inports
      * __vectorField__ 2-dimensional vector field (with vectors of
      two components thus two values within each voxel)
      This processor deals with 2-dimensional data only, therefore it is assumed
      the z-dimension will have size 1 otherwise the 0th slice of the volume
      will be processed.
      * __texture__ Texture to be convolved along the streamlines.

    ### Outports
      * __image__ The image resulting from smearing the given texture
      the streamlines of the given vector field.
*/
class IVW_MODULE_LABLIC_API LICProcessor : public Processor {
    // Friends
    // Types
public:
    // Construction / Deconstruction
public:
    LICProcessor();
    virtual ~LICProcessor() = default;

    // Methods
public:
    virtual const ProcessorInfo getProcessorInfo() const override;
    static const ProcessorInfo processorInfo_;

protected:
    /// Our main computation function
    virtual void process() override;

    // Helper functions to be defined here and then implemented in the .cpp
    std::string standardLIC(const VectorField2 &vectorField, const RGBAImage &inTex,
                                   RGBAImage &outImg);
    std::string parallelLIC(const VectorField2 &vectorField, const RGBAImage &inTex,
                            RGBAImage &outImg);
    std::string fastLIC(const VectorField2 &vectorField, const RGBAImage &inTex,
                                   RGBAImage &outImg);

    void enhanceContrast(RGBAImage& img, const float& desiredMu,
                         const float& desiredSigma);
    void applyColor(const VectorField2 &vectorField, RGBAImage &img);

    dvec2 PixelToGrid(const VectorField2 &vectorField, const size2_t &pixel);
    size2_t GridToPixel(const VectorField2 &vectorField, const dvec2 &grid);

    // Ports
public:
    // Input vector field
    VolumeInport volumeIn_;

    // Input texture
    ImageInport noiseTexIn_;

    // Output image
    ImageOutport licOut_;

    // Properties
public:
    // Declare properties
    BoolProperty propColoredTexture;
    IntProperty propKernelRadius;
    DoubleProperty propStepSize;
    BoolProperty propEnhancedContrast;
    FloatProperty propDesiredMu;
    FloatProperty propDesiredSigma;
    TemplateOptionProperty<int> propLICType;
    BoolProperty propInvalidate;

    // Attributes
private:
    size3_t vectorFieldDims_;
    size2_t texDims_;
};

}  // namespace inviwo
