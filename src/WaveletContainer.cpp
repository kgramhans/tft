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
#include <cmath>
#include <cstring>
#include "tft/WaveletContainer.h"
#include "tft/ConfinedGaussianWaveletVoice.h"

TFT::WaveletContainer::WaveletContainer(unsigned int nOctaves,
                                     double fmax,
                                     float Q,
                                     float overlapPercentage) : dyadicFilter(nOctaves), nSamples(0)
{
   assert(nOctaves > 0);
   assert(fmax > 0 && fmax <= 0.5);
   assert(Q > 0.5);
   assert(overlapPercentage < 100); // negative is ok, though unusual
   
   // Calculate increment
   float increment = (Q - 0.5) / (Q + 0.5 - overlapPercentage/100) ;
   assert(increment < 1); //Avoid recursion
   
   double fCenter = fmax;  // Place first filter as high as possible
   double flo = fCenter * std::sqrt(increment);
   double fhi = flo / increment;
   int nVoices = (int)(nOctaves * std::log(2) /( -std::log(increment)));
   for (int i = nVoices; i--; fCenter *= increment, fhi = flo, flo *= increment)
   {
      assert(fhi > fCenter);
      assert(fCenter > flo);
      waveletVoices.push_back(new ConfinedGaussianWaveletVoice(fCenter, flo, fhi, Q, overlapPercentage, &dyadicFilter));
   }
}

TFT::WaveletContainer::~WaveletContainer()
{
   while (!waveletVoices.empty())
   {
      delete waveletVoices.back();
      waveletVoices.pop_back();
   }
}

/**
 Tell caller about the amount of padding required (left,right) when making calls to doTransform.
 This is essential information to the caller when seeking to avoid artificial transients at the edges of signal
 being investigated
 */
void TFT::WaveletContainer::getRequiredPaddingSamples(unsigned int &nPre, unsigned int & nPost) const
{
    unsigned int pre(0);
    unsigned int post(0);
    nPre = 0;
    nPost = 0;
    for (auto iter = waveletVoices.begin(); iter != waveletVoices.end(); iter++)
    {
        (*iter)->getRequiredPaddingSamples(pre, post);
        nPre = std::max(pre, nPre);
        nPost = std::max(post, nPost);
    }
}




