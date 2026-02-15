#include "tft/StftMoment.h"
#include <cstring>
#include <iostream>

namespace TFT {

StftMoment::StftMoment(unsigned int momentInTime,
                       unsigned int windowCenterOffset,
                       double overlap,
                       const std::shared_ptr<FFTParam> h, std::shared_ptr<std::vector<float> > _window) :
    window(_window),
    windowEnergy(0),
    windowLength((_window.get())->size()),
    duration((_window.get())->size() * (1.0 - overlap / 100)),
    bandwidth(1.0/(_window.get())->size()),
    moment(momentInTime),
    windowOffset(windowCenterOffset),
    resultZ((_window.get())->size()),
    handle(h),
    momentSignalBuffer(0)
{
    assert(windowLength && (windowLength & (windowLength - 1)) == 0);        ///< Must be power-of-two
    assert(overlap < 100);
    assert(h->Points == windowLength / 2);
    assert(windowOffset < windowLength);

    // Calculate window energy
    for (int i = 0; i < windowLength; i++) {
        TF_DATA_TYPE val = (*window.get())[i];
        windowEnergy += val * val;
    }
}

StftMoment::~StftMoment() {
}

void StftMoment::transform(const TF_DATA_TYPE * pSamples) {
    // Construct a buffer for our windowed signal
    // Note that window is assumed to provice proper scaling
    std::vector<TF_DATA_TYPE> scratch(windowLength);
    pSamples += getConstructedMomentOffset();
    for (int i = 0; i < windowLength; i++) {
        scratch[i] = (*window.get())[i] * *pSamples++;
    }

    RealFFTf(&scratch[0], handle.get());

    // Do the bit-reverse into result. Package DC and N/2 into first complex value
    resultZ[0] = scratch[0];
    resultZ[1] = scratch[1];
    for(int i = 1; i < handle->Points; i++) {
        const int index = handle->BitReversed[i];
        resultZ[2 * i    ] = scratch[index];
        resultZ[2 * i + 1] = scratch[index + 1];
    }
}

std::vector<TF_DATA_TYPE> StftMoment::constructMomentSignal(const std::unique_ptr<IRegion> & region) const {
    // If someone has already constructed a result, we will be using that
    if (!momentSignalBuffer.empty()) {
        auto rval = momentSignalBuffer;
        momentSignalBuffer.clear();
        return rval;
    }

    std::vector<TF_DATA_TYPE> rval(2 * handle->Points);
    std::vector<TF_DATA_TYPE> scratch(windowLength);
    std::memcpy(&scratch[0], &resultZ[0], sizeof(TF_DATA_TYPE) * windowLength);
    // Now is the time for inverse FFT and proper normalisation
    // Firstly, we need to map out all those values that are outside our region
    for (int i = 0; i <= handle->Points; i++) {
        float frequency = (0.5 * i) / handle->Points;
        if (!region->isWithin(moment, frequency)) {
            // Zero this
            if (i == 0) {
                scratch[0] = 0;
            } else if (i == handle->Points) {
                scratch[1] = 0;
            } else {
                scratch[2 * i] = scratch[2 * i + 1] = 0;
            }
        }
    }

    // Next, we transform back
    InverseRealFFTf(&scratch[0], handle.get());

    double normFactor = 2 * bandwidth * duration / windowEnergy;
    // And do the bitreverse thing while normalizing. Not sure about factor "N". Need some testing here...
    for(int i = 0; i < handle->Points; i++) {
        const int index = handle->BitReversed[i];
        rval[2 * i    ] = scratch[index    ] * normFactor;
        rval[2 * i + 1] = scratch[index + 1] * normFactor;
    }

    return rval;
}

void StftMoment::constructMomentSignalBuffer(const std::unique_ptr<IRegion> & region) const {
    momentSignalBuffer = constructMomentSignal(region);
}

} // namespace TFT

