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

    // LIC runs here
    LogProcessorInfo(standardLIC(vectorField, texture, licImage));
    //LogProcessorInfo(parallelLIC(vectorField, texture, licImage));
    //LogProcessorInfo(fastLIC(vectorField, texture, licImage));

    licOut_.setData(outImage);
}

std::string LICProcessor::fastLIC(const VectorField2& vectorField, const RGBAImage& inTex,
                                  RGBAImage& outImg) {
    Integrator::RK4CallCounter = 0;

    int steps = 60;
    double stepSize = 0.003;

    int counter;
    double sum;

	unsigned int skipCount = 0;

    time_t startTime;
    time(&startTime);

	dvec2 gridCenter;
    dvec2 current;

    std::list<size2_t> forwardSteps;
    std::list<size2_t> backwardSteps;

    std::vector<std::vector<unsigned char>> buffer(texDims_.x, std::vector<unsigned char>(texDims_.y, 0));

    for (size_t j = 0; j < texDims_.y; j++) {
        for (size_t i = 0; i < texDims_.x; i++) {

            size2_t center = {i, j};
            if (buffer[center.x][center.y] == 0) {

                gridCenter = PixelToGrid(vectorField, center);

                // Backward integration
                backwardSteps.clear();
                current = gridCenter;
                for (int k = 0; k < steps; k++) {
                    dvec2 next;
                    if (!Integrator::RK4Lite(vectorField, current, next, stepSize, true, true))
                        break;
                    backwardSteps.push_back(GridToPixel(vectorField, next));
                    current = next;
                }

				// Forward integration (done in this order because we need to keep the last "current" value)
                forwardSteps.clear();
                current = gridCenter;
                for (int k = 0; k < steps; k++) {
                    dvec2 next;
                    if (!Integrator::RK4Lite(vectorField, current, next, stepSize, true, false))
                        break;
                    forwardSteps.push_back(GridToPixel(vectorField, next));
                    current = next;
                }

                // Sum all values from kernel

                // Central pixel
                counter = 1;
                sum = inTex.readPixelGrayScale(center);

                // Forward pixels
                for (auto it = forwardSteps.begin(); it != forwardSteps.end(); it++) {
                    sum += inTex.readPixelGrayScale(*it);
                    counter++;
                }

                // Backward pixels
                for (auto it = backwardSteps.begin(); it != backwardSteps.end(); it++) {
                    sum += inTex.readPixelGrayScale(*it);
                    counter++;
                }

                if (counter == 1) {
                    sum = 0;  // If only the center pixel was sampled, it means it is over a sink
                }
                
				// Average grayscale values along the sampled pixels
				outImg.setPixelGrayScale(center, sum / counter);

				// Mark as visited
                buffer[center.x][center.y] = 1;

				bool reachedEnd = false;

                // Begin iterating forward
                while (!forwardSteps.empty()) {
                    // Pop last pixel in backward vector (if the backward vector is at maximum
                    // length)
                    if (backwardSteps.size() == steps) {
                        sum -= inTex.readPixelGrayScale(backwardSteps.back());
                        backwardSteps.pop_back();
                        counter--;
                    }

                    // Move center pixel to the front of backward vector
                    backwardSteps.push_front(center);

                    // Move first pixel from forward vector to center
                    center = forwardSteps.front();
                    forwardSteps.pop_front();

					// Check for loops
                    if (buffer[center.x][center.y] != 0) {
                        skipCount++;
                        break;
                    }

					// Integrate next step and put at end of forward vector (if there is a next step)
                    if (!reachedEnd) {
                        dvec2 next;
						if (Integrator::RK4Lite(vectorField, current, next, stepSize, true, false)) {
							forwardSteps.push_back(GridToPixel(vectorField, next));
							sum += inTex.readPixelGrayScale(forwardSteps.back());
							current = next;
							counter++;
						}
						else {
							reachedEnd = true; // No more forward steps possible, let forward list run empty
						}
                    }

					// Recalculate the average for the new center pixel based on the new kernel values and mark it as visited
                    outImg.setPixelGrayScale(center, sum / counter);
                    buffer[center.x][center.y] = 1;
                }
            } else {
                skipCount++;
            }
        }
    }

    time_t endTime;
    time(&endTime);

    return "Fast LIC completed in " + toString(difftime(endTime, startTime)) + " seconds, with " 
		+ toString(Integrator::RK4CallCounter) + " RK4 calls and " + toString(skipCount) + " skips.";
}

std::string LICProcessor::standardLIC(const VectorField2& vectorField, const RGBAImage& inTex,
                                      RGBAImage& outImg) {

    Integrator::RK4CallCounter = 0;

    int steps = 60;
    double stepSize = 0.003;

    int counter;
    double sum;

	dvec2 current;

	std::vector<size2_t> forwardSteps;
    std::vector<size2_t> backwardSteps;

	forwardSteps.reserve(steps);
    backwardSteps.reserve(steps);

    time_t startTime;
    time(&startTime);

    for (size_t j = 0; j < texDims_.y; j++) {
        for (size_t i = 0; i < texDims_.x; i++) {

            dvec2 gridCenter = PixelToGrid(vectorField, {i, j});

            // Forward integration
            forwardSteps.clear();
            current = gridCenter;
            for (int k = 0; k < steps; k++) {
                dvec2 next;
                if (!Integrator::RK4Lite(vectorField, current, next, stepSize, true, false)) break;
                forwardSteps.push_back(GridToPixel(vectorField, next));
                current = next;
            }

            // Backward integration
            backwardSteps.clear();
            current = gridCenter;
            for (int k = 0; k < steps; k++) {
                dvec2 next;
                if (!Integrator::RK4Lite(vectorField, current, next, stepSize, true, true)) break;
                backwardSteps.push_back(GridToPixel(vectorField, next));
                current = next;
            }

            // Sum all values from kernel

            // Central pixel
            counter = 1;
            sum = inTex.readPixelGrayScale(size2_t(i, j));

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
                sum /= counter;  // Average grayscale values along the sampled pixels
            } else {
                sum = 0;  // If only the center pixel was sampled, it means it is over a sink
            }
            outImg.setPixelGrayScale(size2_t(i, j), sum);
        }
    }

    time_t endTime;
    time(&endTime);

    return "Standard LIC completed in " + toString(difftime(endTime, startTime)) + " seconds, with " +
           toString(Integrator::RK4CallCounter) + " RK4 calls.";
}

std::string LICProcessor::parallelLIC(const VectorField2& vectorField, const RGBAImage& inTex,
                                      RGBAImage& outImg) {
    int steps = 60;
    double stepSize = 0.008;

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
                forwardSteps.push_back(GridToPixel(vectorField, next));
                current = next;
            }

            // Backward integration
            std::vector<size2_t> backwardSteps;
            current = gridCenter;
            for (int k = 0; k < steps; k++) {
                dvec2 next;
                if (!Integrator::RK4Lite(vectorField, current, next, stepSize, true, true)) break;
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
                sum /= counter;  // Average grayscale values along the sampled pixels
            } else {
                sum = 0;  // If only the center pixel was sampled, it means it is over a sink
            }
            outImg.setPixelGrayScale(size2_t(i, j), sum);
        }
    }

    time_t endTime;
    time(&endTime);

	// RK4 counter is not thread-safe, so this information is ommited here
    return "Parallel LIC completed in " + toString(difftime(endTime, startTime)) + ".";
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
