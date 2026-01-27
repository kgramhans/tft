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

//
//  DyadicFilter.hpp
//  tfq
//
//  Created by Klaus Gram-Hansen on 12/11/2025.
//

#ifndef DyadicFilter_h
#define DyadicFilter_h

#include <utility>
#include <vector>
#include <limits>
#include "version.h"

namespace TFT {
    class DyadicFilter
    {
    public:
        /**
        * @brief DyadicFilter
        * @param nOctaves is the number of octaves covered by the filter
        */
        DyadicFilter(const unsigned int nOctaves);
        ~DyadicFilter();

        /**
         * @brief getExtraSamples needed for context on either side of a signal to be filtered. This static member is for information purposes only
         * @param nOctaves
         * @return number of extra samples needed before and after a signal to be filtered
         */
        static unsigned int  getExtraSamples(unsigned int nOctaves);

       /**
        Find lowest octave that still contains frequencies up to fmax
        @param fmax is a positive number
        @return index of octave containing frequencies. If fmax > 0.5, return is still 0 for robustness even though frequencies above 0.5 make no sense
        */
       unsigned int findOctave(const float fmax) const;

        /**
        * @brief filterSamples Perform filtering operation on a signal of samples
        * @param pSamples Points to start of signal buffer
        * @param n_samples Number of samples in scope
        * @param n_samples_before Number of samples available before pSamples
        * @param n_samples_after Number of samples available after pSamples + n_samples - 1
        */
        void filterSamples(const TF_DATA_TYPE * pSamples, const unsigned int n_samples, unsigned int n_samples_before, unsigned int n_samples_after);

        /**
         * @brief getSamples
         * @param octave
         * @param fromSample Integer specifying offset from pSamples passed previously to filterSamples
         * @return a pair of pointer to data and length of array
         */
        const std::pair<TF_DATA_TYPE *, unsigned int> getSamples(const unsigned int octave, const int fromSample) const;

        /**
         * @brief doAllocation allocates internal buffers BEFORE calling filterSamples
         * @param nSamples
         * @param nSamplesBefore
         * @param nSamplesAfter
         */
        void doAllocation(unsigned int nSamples, unsigned int nSamplesBefore, unsigned int nSamplesAfter);

        /**
         * @brief upsample (interpolate) a signal by a factor 2
         * @param pSamples points to start of data
         * @param n_samples length of signal
         * @param n_samples_before valid samples before pSamples. Must equal or exceed getUpsamplingPaddingSize()
         * @param n_samples_after valid samples after pSamples + n_samples - 1. Must equal or exceed getUpsamplingPaddingSize()
         * @return the upsampled signal padded with getUpsamplingPaddingSize() samples in either end
         */
        static std::vector<TF_DATA_TYPE> upsample(const TF_DATA_TYPE * pSamples, unsigned int n_samples, unsigned int n_samples_before, unsigned int n_samples_after);
        /**
         * @brief upsample
         * @param vSamples vector of samples to be padded with getUpsamplingPaddingSize() samples in either end
         * @return the upsampled signal padded with getUpsamplingPaddingSize() samples in either end
         */
        static std::vector<TF_DATA_TYPE> upsample(const std::vector<TF_DATA_TYPE> vSamples);

        /**
        * @brief getUpsamplingPaddingSize
        * @return number of samples needed to be padded in either end of a signal when upsampling
        */
        static int getUpsamplingPaddingSize() { return cstFilterTaps >> 1;}
    private:
        std::vector<TF_DATA_TYPE *> vBufferBegin;      // Point to start of allocated buffer
        std::vector<TF_DATA_TYPE *> vBufferTimeZero;   // Point to location of time zero in allocated buffer
        std::vector<unsigned int> vBufferLengths;            // Point to total length of allocated buffer (starting from BufferBegin)
        void verify();
        bool isNaN(TF_DATA_TYPE val) {
            return val != val;
        }

       /*

       FIR filter designed with
       http://t-filter.engineerjs.com

       sampling frequency: 1000 Hz

       * 0 Hz - 200 Hz
         gain = 1
         desired ripple = 0.1 dB
         actual ripple = 0.07053769787516298 dB

       * 300 Hz - 500 Hz
         gain = 0
         desired attenuation = -90 dB
         actual attenuation = -90.38473284370106 dB

       */

        static constexpr unsigned int cstFilterTaps = 39;
        static constexpr float cstFreqLimit = 0.2;
        static constexpr TF_DATA_TYPE cstMagic = std::numeric_limits<TF_DATA_TYPE>::signaling_NaN();
        static const TF_DATA_TYPE filter_taps[cstFilterTaps];
    };
}
#endif /* DyadicFilter_h */
