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
#include "DyadicFilter.h"
#include "WaveletVoice.h"

namespace TFT {


/**
     */
class WaveletContainer
{
public:
protected:
    WaveletContainer( unsigned int nOctaves,
                     double fmax,
                     float Q,
                     float overlapPercentage);
    ~WaveletContainer();
    DyadicFilter dyadicFilter;
    std::vector<WaveletVoice *> waveletVoices;
    unsigned int nSamples;

    /**
        Tell caller about the amount of padding required (left,right) when making calls to doTransform.
        This is essential information to the caller when seeking to avoid artificial transients at the edges of signal
        being investigated
        */
    virtual void getRequiredPaddingSamples(unsigned int &nPre, unsigned int & nPost) const;

};
}
