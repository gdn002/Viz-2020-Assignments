/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Monday, October 02, 2017 - 13:31:17
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <inviwo/core/datastructures/volume/volumeram.h>
#include <lablic/licprocessor.h>
#include <labstreamlines/integrator.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo LICProcessor::processorInfo_{
    "org.inviwo.LICProcessor",  // Class identifier
    "LICProcessor",             // Display name
    "KTH Labs",                 // Category
    CodeState::Experimental,    // Code state
    Tags::None,                 // Tags
};

const ProcessorInfo LICProcessor::getProcessorInfo() const { return processorInfo_; }

LICProcessor::LICProcessor()
    : Processor()
    , volumeIn_("volIn")
    , noiseTexIn_("noiseTexIn")
    , licOut_("licOut")
// TODO: Register additional properties
{
    // Register ports
    addPort(volumeIn_);
    addPort(noiseTexIn_);
    addPort(licOut_);

    // Register properties
    // TODO: Register additional properties
}

void LICProcessor::process() {
    // Get input
    if (!volumeIn_.hasData()) {
        return;
    }

    if (!noiseTexIn_.hasData()) {
        return;
    }

    auto vol = volumeIn_.getData();
    const VectorField2 vectorField = VectorField2::createFieldFromVolume(vol);
    vectorFieldDims_ = vol->getDimensions();

    auto tex = noiseTexIn_.getData();
    const RGBAImage texture = RGBAImage::createFromImage(tex);
    texDims_ = tex->getDimensions();

    double value = texture.readPixelGrayScale(size2_t(0, 0));

    LogProcessorInfo(value);

    // Prepare the output, it has the same dimensions as the texture and rgba values in [0,255]
    auto outImage = std::make_shared<Image>(texDims_, DataVec4UInt8::get());
    RGBAImage licImage(outImage);

    std::vector<std::vector<double>> licTexture(texDims_.x, std::vector<double>(texDims_.y, 100.0));

    // Hint: Output an image showing which pixels you have visited for debugging
    std::vector<std::vector<int>> visited(texDims_.x, std::vector<int>(texDims_.y, 0));
    // TODO: Implement LIC and FastLIC
    // This code instead sets all pixels to the same gray value
    for (size_t j = 0; j < texDims_.y; j++) {
        for (size_t i = 0; i < texDims_.x; i++) {
            int val = int(licTexture[i][j]);
            licImage.setPixel(size2_t(i, j), dvec4(val, val, val, 255));
            // or
            licImage.setPixelGrayScale(size2_t(i, j), val);
        }
    }

    standardLIC(vectorField, texture, licImage);

    licOut_.setData(outImage);
}

std::string LICProcessor::standardLIC(const VectorField2& vectorField, const RGBAImage& inTex,
                                      RGBAImage& outImg) {
    std::string msg;
    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;
    std::vector<dvec2> points;
    int num_steps, L = 10;
    double sum;
    for (size_t j = 0; j < texDims_.y; j++) {
      for (size_t i = 0; i < texDims_.x; i++) {
        msg = Integrator::RK4LoopV2(vectorField, {i, j}, mesh, vertices, points, num_steps, 0.1f, 0.01f,
                                      2.0f, true, {255, 255, 255, 255}, L);
        sum = 0.0f;
        for (int k=0; k<points.size(); k++){
          sum += inTex.readPixelGrayScale(points[k]);
          sum /= points.size();
        }
        outImg.setPixelGrayScale(size2_t(i,j), sum);
      }
    }
    return msg;
}


}  // namespace inviwo
