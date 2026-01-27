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
#include <vector>
#include <utility>
#include "version.h"

namespace TFT {


    /**
     * An generic interface for performing transformation, modification and inverse transformation on 1D signals
     * As opposed to TimeFrequencyCalculator, this class buffers transformed data internally
     *
     * Known subclasses: WaveletTransformer
     */
    class ITimeFrequencyTransformer
    {
    public:
        /**
         * @brief prepare for transformation: Tell how much data is needed before and after a given signal
         * @param nSamples length of signal to-be-transformed
         * @param nPre number of samples needed prior to signal
         * @param nPost number of samples needed after signal
         */
        virtual void  prepare(unsigned int nSamples, unsigned int & nPre, unsigned int & nPost) = 0;

        /**
         * @brief forwardTransform
         * @param pSamples Points to samples to-be-transformed
         * @param nSamples Number of samples to be transformed
         * @param nValidSamplesBefore Available samples before signal start. Zero-padding will occur if lower than nPre
         * @param nValidSamplesAfter Available samples after signal end. Zero-padding will occur if lower than nPost
         * @return
         */
        virtual unsigned int forwardTransform(const TF_DATA_TYPE* pSamples, unsigned int nSamples, unsigned int nValidSamplesBefore, unsigned int nValidSamplesAfter) = 0;

        /**
         * @brief backwardTransform. Bring back the signal in the time-domain
         * @param signal This is whereto the signal will be stored
         * @param fromSample starting sample of transform (< nSamples from forward transform)
         * @param toSample ending sample of transform (<= nSamples from forward transform)
         * @return length of signal upon return
         */
        virtual unsigned int backwardTransform(std::vector<TF_DATA_TYPE> & signal, unsigned int fromSample, unsigned int toSample) const = 0;

        /**
         * @brief backwardTransform short form
         * @param signal allocated to hold generated signal
         * @return  length of signal upon return
         */
        virtual unsigned int backwardTransform(std::vector<TF_DATA_TYPE> & signal) const { return backwardTransform(signal, 0, signal.size());}

        /**
         * @brief setPolygonRegion Set a region in time/frequency plane. Only points inside the region will be used for backwardTransform
         * @param region : Sequence of inter-connected polygon points. Last point is connected to first point
         * @return true if region could be interpreted
         */
        virtual bool setPolygonRegion(const std::vector<std::pair<float, float>> & region) = 0;

        /**
        Virtual destructor does nothing - but must be defined
        */
        virtual     ~ITimeFrequencyTransformer() {}
    };
}
