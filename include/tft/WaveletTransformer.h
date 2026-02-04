/**
Time Frequency Calculator Library
Copyright (C) 2025  Klaus Gram-Hansen

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#pragma once
#include <stddef.h>
#include "WaveletContainer.h"
#include "TimeFrequencyTransformer.h"
#include "PolygonRegion.h"

namespace TFT {


/**
     */
class WaveletTransformer : public WaveletContainer, public ITimeFrequencyTransformer, protected PolygonRegion
{
public:
    /**
     * @brief WaveletTransformer Construct a wavelet transformer
     * @param nOctaves Number of octaves to cover
     * @param fmax maximum frequency to cover
     * @param Q Q factor of implied Gaussian wavelet
     * @param overlapPercentage Overlap of transform. Applies to both time and frequency domains. More overlap means more work, but improves quality. 75% is a good default value
     */
    WaveletTransformer( unsigned int nOctaves,
                      double fmax,
                      float Q,
                      float overlapPercentage);
    ~WaveletTransformer();

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
     * @return true if region could be interpreted
     */
    virtual void setPolygonRegion(const std::vector<std::pair<float, float>> & region) override;
};
}
