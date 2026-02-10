#ifndef STFTTRANSFORMER_H
#define STFTTRANSFORMER_H

#include "PolygonRegionVertical.h"
#include "tft.h"
#include "tft/ConfinedGaussian.h"
#include "tft/StftMoment.h"

namespace TFT {

class StftTransformer : public ITimeFrequencyTransformer, protected PolygonRegionVertical, protected ConfinedGaussian
{
public:
    /**
     * @brief StftTransformer Construct a STFT transformer
     * @param nOctaves Number of octaves to cover
     * @param fmax maximum frequency to cover
     * @param Q Q factor of implied Gaussian wavelet
     * @param overlapPercentage Overlap of transform. Applies to both time and frequency domains. More overlap means more work, but improves quality. 75% is a good default value
     */
    StftTransformer(unsigned int windowLength,
                       float overlapPercentage);
    ~StftTransformer();

    /**
     * @brief prepare See parent interface interface ITimeFrequencyTransformer
     * @param nSamples
     * @param nPre
     * @param nPost
     */
    virtual void  prepare(unsigned int nSamples, unsigned int & nPre, unsigned int & nPost) override;

    /**
         * @brief prepareParallelSequences
         * @param pSamples
         * @param nSamples
         * @param nValidSamplesBefore
         * @param nValidSamplesAfter
         * @return Number of sequences that can be invoked in parallel
         */
    virtual int prepareParallelForwardSequences(const TF_DATA_TYPE* pSamples, unsigned int nSamples, unsigned int nValidSamplesBefore, unsigned int nValidSamplesAfter) override;

    /**
        Execute a given sequence as prepared above. Sequences can be executed in any order, even parallel
        @param iSequence ranges from 0 to number of sequences minus 1
        */
    virtual void executeForwardSequence(int iSequence) override;

    /**
     * @brief forwardTransform See parent interface ITimeFrequencyTransformer
     * @param pSamples
     * @param nSamples
     * @param nValidSamplesBefore
     * @param nValidSamplesAfter
     * @return
     */
    virtual unsigned int forwardTransform(const TF_DATA_TYPE* pSamples, unsigned int nSamples, unsigned int nValidSamplesBefore, unsigned int nValidSamplesAfter) override;

    /**
     * @brief backwardTransform See parent interface interface ITimeFrequencyTransformer
     * @param signal
     * @param fromSample
     * @param toSample
     * @return
     */
    virtual unsigned int backwardTransform(std::vector<TF_DATA_TYPE> & signal) const override;

    /**
         * @brief prepareParallelSequences
         * @return Number of sequences that can be invoked in parallel
         */
    virtual int prepareParallelBackwardSequences() const override;

    /**
        Execute a given sequence as prepared above. Sequences can be executed in any order, even parallel
        @param iSequence ranges from 0 to number of sequences minus 1
        */
    virtual void executeBackwardSequence(int iSequence) const override;

    /**
     * @brief setPolygonRegion Set a region in time/frequency plane. Only points inside the region will be used for backwardTransform
     * @param region : Sequence of inter-connected polygon points. Last point is connected to first point
     * @return void
     */
    virtual void setPolygonRegion(const std::vector<std::pair<float, float>> & region) override;

protected:
    std::vector<StftMoment> stftMoments;
    float overlap;
    unsigned int nSamples;
    unsigned int windowLength;
    unsigned int windowCenterOffset;
    std::shared_ptr<std::vector<TF_DATA_TYPE>> pWindow;
    std::shared_ptr<FFTParam> hFFT;
    unsigned int momentStep;
    std::vector<TF_DATA_TYPE> input;
    unsigned int pre;
    unsigned int post;

private:
    static constexpr float sigma = 0.1;
};

}
#endif // STFTTRANSFORMER_H
