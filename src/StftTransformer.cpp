#include <tft/StftTransformer.h>
#include <cstring>

namespace TFT {

/**
     * @brief StftTransformer Construct a STFT transformer
     * @param nOctaves Number of octaves to cover
     * @param fmax maximum frequency to cover
     * @param Q Q factor of implied Gaussian wavelet
     * @param overlapPercentage Overlap of transform. Applies to both time and frequency domains. More overlap means more work, but improves quality. 75% is a good default value
     */
StftTransformer::StftTransformer(unsigned int windowLength,
                                 float overlapPercentage) : ConfinedGaussian(sigma),
    overlap(overlapPercentage), nSamples(0), windowLength(windowLength), momentStep((unsigned int)((1.0 - overlapPercentage / 100.0) * windowLength + 0.5)),
    pre(0), post(0)
{
    if (momentStep <= 0) {
        momentStep = 1;
    }
    assert(overlap < 100);
    assert(windowLength && (windowLength & (windowLength - 1)) == 0);        ///< Must be power-of-two

    // Adjust the overlap
    overlap = 100.0 * (windowLength - momentStep) / (float) windowLength;

    // Create the FFT handle we are going to use along with the particular window function
    pWindow = std::make_shared<std::vector<TF_DATA_TYPE>>(windowLength, 0);   // Simple square window. Note that average must be two due to one-sided spectrum
    windowCenterOffset = windowLength / 2 - 1;

    // Now, construct a confined Gaussian window
    setLength(windowLength - 1);
    double sum = 0;
    for (int i = 0; i < windowLength - 1; i++) {
        (*pWindow)[i] = approximateConfinedGaussian(i);
        sum += (*pWindow)[i];
    }
    TF_DATA_TYPE factor = 1.0 / sum;
    for (int i = 0; i < windowLength - 1; i++) {
        (*pWindow)[i] *= factor;
    }

    hFFT = GetFFT(windowLength);
}

StftTransformer::~StftTransformer() {

}

/**
     * @brief prepare See parent interface interface ITimeFrequencyTransformer
     * @param nSamples
     * @param nPre
     * @param nPost
     */
void  StftTransformer::prepare(unsigned int _nSamples, unsigned int & nPre, unsigned int & nPost) {
    assert(_nSamples > 0);
    if (nSamples != _nSamples) {
        // re-allocate the input
        input.clear();
        stftMoments.clear();

        nSamples = _nSamples;
        pre = nPre = windowCenterOffset;

        int nMoments = 1 + (nSamples - 1) / momentStep;
        int left = nSamples - (1 + (nMoments - 1) * momentStep);
        assert(left >= 0);
        if (2 * left >= momentStep) {
            nMoments++;
        }
        post = nPost = 1 + (nMoments - 1) * momentStep - nSamples + windowCenterOffset;
        // allocate the input and moments
        input.resize(pre + nSamples + post);
        for (int i = 0; i < nMoments; i++) {
            stftMoments.push_back(StftMoment(i * momentStep, windowCenterOffset, overlap, hFFT, pWindow ));
        }
    } else {
        assert(input.size());
        assert(stftMoments.size());
        nPre = pre;
        nPost = post;
    }
}

/**
         * @brief prepareParallelSequences
         * @param pSamples
         * @param nSamples
         * @param nValidSamplesBefore
         * @param nValidSamplesAfter
         * @return Number of sequences that can be invoked in parallel
         */
int StftTransformer::prepareParallelForwardSequences(const TF_DATA_TYPE* pSamples, unsigned int nSamples, unsigned int nValidSamplesBefore, unsigned int nValidSamplesAfter) {
    prepare(nSamples, pre, post);

    // Take in samples
    int missingLeft = pre - nValidSamplesBefore;
    int missingRight = post - nValidSamplesAfter;
    if (missingLeft < 0)
        missingLeft = 0;
    if (missingRight < 0)
        missingRight = 0;

    // Pad left
    memset(&input[0], 0, sizeof(TF_DATA_TYPE) * missingLeft);
    // Pad right
    memset(&input[pre + nSamples + post - missingRight], 0, sizeof(TF_DATA_TYPE) * missingRight);
    // Fill in samples
    memcpy(&input[missingLeft], pSamples - pre + missingLeft, (pre + nSamples + post - missingLeft - missingRight) * sizeof(TF_DATA_TYPE));

    return stftMoments.size();
}

/**
        Execute a given sequence as prepared above. Sequences can be executed in any order, even parallel
        @param iSequence ranges from 0 to number of sequences minus 1
        */
void StftTransformer::executeForwardSequence(int iSequence) {
    stftMoments[iSequence].transform(&input[pre]);
}

/**
     * @brief forwardTransform See parent interface ITimeFrequencyTransformer
     * @param pSamples
     * @param nSamples
     * @param nValidSamplesBefore
     * @param nValidSamplesAfter
     * @return
     */
unsigned int StftTransformer::forwardTransform(const TF_DATA_TYPE* pSamples, unsigned int nSamples, unsigned int nValidSamplesBefore, unsigned int nValidSamplesAfter) {
    prepareParallelForwardSequences(pSamples, nSamples, nValidSamplesBefore, nValidSamplesAfter);

    for (int i = 0; i < stftMoments.size(); i++) {
        executeForwardSequence(i);
    }

    return stftMoments.size() * hFFT->Points;
}

/**
     * @brief backwardTransform See parent interface interface ITimeFrequencyTransformer
     * @param signal
     * @return
     */
unsigned int StftTransformer::backwardTransform(std::vector<TF_DATA_TYPE> & signal) const {
    signal = std::vector<TF_DATA_TYPE>(nSamples, 0);
    for (auto iter = stftMoments.cbegin(); iter != stftMoments.cend(); iter ++) {
        // Accumulate into signal
        auto partialSignal = iter->constructMomentSignal(this->cloneRegion());
        int partial_begins_at = iter->getConstructedMomentOffset();
        int partial_length = partialSignal.size();
        int partial_inx = 0;
        if (partial_begins_at < 0) {
            partial_inx += -partial_begins_at;
            partial_length -= -partial_begins_at;
            partial_begins_at = 0;
        }
        if (partial_begins_at + partial_length  > nSamples) {
            partial_length -= (partial_begins_at + partial_length  - nSamples);
        }

        TF_DATA_TYPE * dst = &signal        [partial_begins_at];
        TF_DATA_TYPE * src = &partialSignal [partial_inx];
        for (int i = partial_length; i--; ) {
            *dst++ += *src++;
        }
    }
    return signal.size();
}

/**
         * @brief prepareParallelSequences
         * @return Number of sequences that can be invoked in parallel
         */
int StftTransformer::prepareParallelBackwardSequences() const {
    return stftMoments.size();
}

/**
        Execute a given sequence as prepared above. Sequences can be executed in any order, even parallel
        @param iSequence ranges from 0 to number of sequences minus 1
        */
void StftTransformer::executeBackwardSequence(int iSequence) const {
    stftMoments[iSequence].constructMomentSignalBuffer(this->cloneRegion());
}

/**
     * @brief setPolygonRegion Set a region in time/frequency plane. Only points inside the region will be used for backwardTransform
     * @param region : Sequence of inter-connected polygon points. Last point is connected to first point
     * @return void
     */
void StftTransformer::setPolygonRegion(const std::vector<std::pair<float, float>> & region) {
    setCorners(region);
}





}
