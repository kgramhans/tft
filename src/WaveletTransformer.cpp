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
#include <assert.h>
#include <cstring>
#include "tft/WaveletTransformer.h"

TFT::WaveletTransformer::WaveletTransformer(unsigned int nOctaves,
                                     double fmax,
                                     float Q,
                                     float overlapPercentage) : WaveletContainer(nOctaves, fmax, Q, overlapPercentage)
{
}

TFT::WaveletTransformer::~WaveletTransformer()
{
}


 /**
  Reset any internal state to that of a newly created object. This may be a time-saver as opposed to free/allocate af new object
 */
void  TFT::WaveletTransformer::prepare(unsigned int n_samples, unsigned int & nPre, unsigned int & nPost)
{
    nSamples = n_samples;
    for (auto iter = waveletVoices.begin(); iter != waveletVoices.end(); iter++)
    {
       (*iter)->allocateResult(nSamples, 0, true);
    }

    // Then do allocation of dyadic filter since we now do know requirements
    getRequiredPaddingSamples(nPre, nPost);
    dyadicFilter.doAllocation(nSamples, nPre, nPost);
}

unsigned int TFT::WaveletTransformer::forwardTransform(const TF_DATA_TYPE* pSamples, unsigned int nSamples, unsigned int nValidSamplesBefore, unsigned int nValidSamplesAfter) {
    unsigned int pre, post;
    prepare(nSamples, pre, post);
    // Feed data into our dyadic filter
    dyadicFilter.filterSamples(pSamples, nSamples, std::min(nValidSamplesBefore, pre), std::min(nValidSamplesAfter, post));

    // And then ask all voices to co-operate
    unsigned int rval = 0;
    for (auto iter = waveletVoices.begin(); iter != waveletVoices.end(); iter++)
    {
        rval += (*iter)->transform();
    }
    return rval;
}

int TFT::WaveletTransformer::prepareParallelBackwardSequences() const
{
    return waveletVoices.size();
}

void TFT::WaveletTransformer::executeBackwardSequence(int iSequence) const
{
    assert(iSequence >= 0);
    assert(iSequence < waveletVoices.size());
    waveletVoices[iSequence]->constructVoiceSignalBuffer(*this);
}

unsigned int TFT::WaveletTransformer::backwardTransform(std::vector<TF_DATA_TYPE> & signal) const {
    signal = std::vector<TF_DATA_TYPE>(nSamples, 0);

    // Iterate all voices and ask for contribution
    for (auto iter = waveletVoices.begin(); iter != waveletVoices.end(); iter++)
    {
        auto v = (*iter)->constructVoiceSignal(*this);
        for (int inx = 0; inx < nSamples; inx++) {
            signal[inx] += v[inx];
        }
    }
    return signal.size();
}

/**
 * @brief prepareParallelSequences
 * @param pSamples
 * @param nSamples
 * @param nValidSamplesBefore
 * @param nValidSamplesAfter
 * @return Number of sequences that can be invoked in parallel
 */
int TFT::WaveletTransformer::prepareParallelForwardSequences(const TF_DATA_TYPE* pSamples, unsigned int nSamples, unsigned int nValidSamplesBefore, unsigned int nValidSamplesAfter) {
    // Process samples into dyadic filter - otherwise parallel work is not possible
    unsigned int pre, post;
    prepare(nSamples, pre, post);
    // Feed data into our dyadic filter
    dyadicFilter.filterSamples(pSamples, nSamples, nValidSamplesBefore, nValidSamplesAfter);

    return (int) waveletVoices.size();  // All must do some work
}

/**
    Execute a given sequence as prepared above. Sequences can be executed in any order, even parallel
    @param iSequence ranges from 0 to number of sequences minus 1
    */
void TFT::WaveletTransformer::executeForwardSequence(int iSequence) {
    assert(iSequence >= 0 && iSequence < waveletVoices.size());
    waveletVoices[iSequence]->transform();
}

/**
 * @brief setRegion
 * @param region
 */
void TFT::WaveletTransformer::setPolygonRegion(const std::vector<std::pair<float, float>> & region) {
    setCorners(region);
}





