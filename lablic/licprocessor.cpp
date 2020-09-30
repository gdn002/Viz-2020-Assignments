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
#include <limits.h>
#include <algorithm>

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

    std::vector<std::vector<double>> licTexture(texDims_.x, std::vector<double>(texDims_.y, -1.0));

    // Hint: Output an image showing which pixels you have visited for debugging
    std::vector<std::vector<int>> visited(texDims_.x, std::vector<int>(texDims_.y, 0));

	// LIC runs here
    //LogProcessorInfo(standardLIC(vectorField, texture, licImage));
    LogProcessorInfo(fastLIC(vectorField, texture, licImage, licTexture, 50));

    licOut_.setData(outImage);
}

std::string LICProcessor::fastLIC(const VectorField2& vectorField, const RGBAImage& inTex,
    RGBAImage& outImg, std::vector<std::vector<double>> &visited, int kernel_size) {

    int steps = INT_MAX;
    double stepSize = 0.001;

    int counter;
    double sum = 0.0;
    
    time_t startTime;
    time(&startTime);

    for (size_t j = 0; j < texDims_.y; j++) {
        for (size_t i = 0; i < texDims_.x; i++) {

          if(visited[i][j] == -1){
            dvec2 gridCenter = PixelToGrid(vectorField, {i, j});

            // Forward integration
            std::vector<size2_t> forwardSteps;
            dvec2 current = gridCenter;
            forwardSteps.push_back(GridToPixel(vectorField, current)); //the center is at the beginning of the forward integration
            for (int k = 0; k < steps; k++) {
                dvec2 next;
                if (!Integrator::RK4Lite(vectorField, current, next, stepSize, true, false)) break;
                if (k%5==0)
                  forwardSteps.push_back(GridToPixel(vectorField, next));
                current = next;
            }

            // Backward integration
            std::vector<size2_t> backwardSteps;
            current = gridCenter;
            for (int k = 0; k < steps; k++) {
                dvec2 next;
                if (!Integrator::RK4Lite(vectorField, current, next, stepSize, true, true)) break;
                if (k%5==0)
                  backwardSteps.push_back(GridToPixel(vectorField, next));
                current = next;
            }

            //resort and merge the 2 vectors
            std::reverse(backwardSteps.begin(),backwardSteps.end());
            backwardSteps.reserve(backwardSteps.size()+forwardSteps.size());
            backwardSteps.insert(backwardSteps.end(), forwardSteps.begin(), forwardSteps.end());

            //iterators based on the kernel
            std::vector<size2_t>::iterator k_start;
            std::vector<size2_t>::iterator k_center = k_start+kernel_size/2;
            std::vector<size2_t>::iterator k_end = k_start+kernel_size;

            // calculate the first kernel
            for (k_start = backwardSteps.begin(); k_start < k_end; k_start++) {
              sum += inTex.readPixelGrayScale(*k_start);
            }
            visited[k_center->x][k_center->y] = sum/kernel_size;
            outImg.setPixelGrayScale(size2_t(k_center->x, k_center->y), sum/kernel_size);
            k_center++;
            
            // Iterate along the streamline
            for (k_start = backwardSteps.begin(); k_end < backwardSteps.end()-1; k_start++, k_center++, k_end++) {
              sum -= inTex.readPixelGrayScale(*k_start);
              sum += inTex.readPixelGrayScale(*k_end);
              visited[k_center->x][k_center->y] = sum/kernel_size;
              outImg.setPixelGrayScale(size2_t(k_center->x, k_center->y), sum/kernel_size);
            }
          }
        }

    }

	time_t endTime;
    time(&endTime);

    return "Standard LIC completed in " + toString(difftime(endTime, startTime)) + " seconds.";
}
std::string LICProcessor::standardLIC(const VectorField2& vectorField, const RGBAImage& inTex,
                                      RGBAImage& outImg) {
    int steps = 50;
    double stepSize = 0.001;

    int counter;
    double sum;
    
    time_t startTime;
    time(&startTime);

    #pragma omp parallel
    #pragma omp for
    for (size_t j = 0; j < texDims_.y; j++) {
        for (size_t i = 0; i < texDims_.x; i++) {

            dvec2 gridCenter = PixelToGrid(vectorField, {i, j});

            // Forward integration
            std::vector<size2_t> forwardSteps;
            dvec2 current = gridCenter;
            for (int k = 0; k < steps; k++) {
                dvec2 next;
                if (!Integrator::RK4Lite(vectorField, current, next, stepSize, true, false)) break;
                if (k%5==0)
                  forwardSteps.push_back(GridToPixel(vectorField, next));
                current = next;
            }

            // Backward integration
            std::vector<size2_t> backwardSteps;
            current = gridCenter;
            for (int k = 0; k < steps; k++) {
                dvec2 next;
                if (!Integrator::RK4Lite(vectorField, current, next, stepSize, true, true)) break;
                if (k%5==0)
                  backwardSteps.push_back(GridToPixel(vectorField, next));
                current = next;
            }

            // Sum all values from kernel

            // Central pixel
            int counter = 1;
            double sum = inTex.readPixelGrayScale(size2_t(i, j));

            // Forward pixels
            for (int k = 0; k < forwardSteps.size(); k++) {
                sum += inTex.readPixelGrayScale(forwardSteps[k]);
                counter++;
            }

            // Backward pixels
            for (int k = 0; k < backwardSteps.size(); k++) {
                sum += inTex.readPixelGrayScale(backwardSteps[k]);
                counter++;
            }

            if (counter > 1) {
                sum /= counter; // Average grayscale values along the sampled pixels
            } else {
                sum = 0; // If only the center pixel was sampled, it means it is over a sink
            }
            outImg.setPixelGrayScale(size2_t(i, j), sum);
        }
    }

    time_t endTime;
    time(&endTime);

    return "Standard LIC completed in " + toString(difftime(endTime, startTime)) + " seconds.";
}

dvec2 inviwo::LICProcessor::PixelToGrid(const VectorField2& vectorField, const size2_t& pixel) {
    double relationX = (double)pixel.x / texDims_.x;
    double relationY = (double)pixel.y / texDims_.y;

    auto bBoxMin = vectorField.getBBoxMin();
    auto bBoxMax = vectorField.getBBoxMax();

    // For this I will assume pixels are counted Left-Right, Top-Bottom
    double x = (bBoxMin.x * (1 - relationX)) + (bBoxMax.x * relationX);
    double y = (bBoxMax.y * (1 - relationY)) + (bBoxMin.y * relationY);

    return {x, y};
}

size2_t inviwo::LICProcessor::GridToPixel(const VectorField2& vectorField, const dvec2& grid) {
    auto bBoxMin = vectorField.getBBoxMin();
    auto bBoxMax = vectorField.getBBoxMax();

    double relationX = (grid.x - bBoxMin.x) / (bBoxMax.x - bBoxMin.x);
    double relationY = (grid.y - bBoxMin.y) / (bBoxMax.y - bBoxMin.y);

    // For this I will assume pixels are counted Left-Right, Top-Bottom
    double x = texDims_.x * relationX;
    double y = texDims_.y * (1 - relationY);

    // Clamp
    if (x < 0) x = 0;
    if (x >= texDims_.x) x = texDims_.x - 1;
    if (y < 0) y = 0;
    if (y >= texDims_.y) y = texDims_.y - 1;

    return {x, y};
}

}  // namespace inviwo
