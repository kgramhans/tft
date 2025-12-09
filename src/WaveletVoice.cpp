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
//  WaveletVoice.cpp
//  test
//
//  Created by Klaus Gram-Hansen on 17/11/2025.
//

#include "WaveletVoice.h"
#include <iostream>
#include <assert.h>

WaveletVoice::WaveletVoice(const float overlapPercentage,
                           const DyadicFilter * dFilter,
                           const double fCenter) :
   dyadicFilter(dFilter),
   octave(0),
   resultLen(0),
   overlap(overlapPercentage),
   waveletHalfLen(0),
   duration(0),
   transformLength(0),
   resultStep(0),
   frequency(fCenter)
{}

void WaveletVoice::dump()
{
   for (int i = 0; i <= waveletHalfLen; i++)
   {
//      cout << i << " : " << waveletRe[i] << ", " << waveletIm[i] << endl;
   }
}

WaveletVoice::~WaveletVoice()
{
   if(waveletHalfLen)
   {
      delete[]waveletRe;
      delete[]waveletIm;
      waveletHalfLen=0;
   }
   if(resultLen)
   {
      delete[]resultRe;
      delete[]resultIm;
      resultLen=0;
   }
}
   
void WaveletVoice::getRequiredPaddingSamples(unsigned int & pre, unsigned int & post) const
{
   // We need WaveletHalfLen samples from the LP filter before and WaveletHalfLen+1 after
   // These numbers must be multiplied by 1<<octave. Additionally, we must add requirement from LP filter
   assert(waveletHalfLen); // Must be initialised by derived constructor
   assert(transformLength);
   assert(resultLen);
   assert(resultStep > 0);
   pre = (waveletHalfLen << octave) + dyadicFilter->getExtraSamples(octave);
   post = pre + resultLen * resultStep - transformLength;
}

pair<unsigned int, unsigned int> WaveletVoice::calculateResultLenAndStep(unsigned int _resolution) const
{
   assert(transformLength);
   unsigned int resultStep = std::max(transformStep, (_resolution / 2) >> octave ) << octave; // Be more conservative internally with resolution
   unsigned int rval = 1 + (transformLength - 1) / resultStep;
   int remaining = (transformLength - 1) - (rval -1) * resultStep;
   while (2 * remaining >= (long) resultStep)
   {
      remaining -= resultStep;
      assert(2 * remaining < (long) resultStep);
      rval++;
   }
   return pair<unsigned int, unsigned int>(rval, resultStep);
}

void WaveletVoice::allocateResult(unsigned int nSamples, unsigned int _resolution)
{
   if (resultLen)
   {
      // Allocation already done. Get rid of that if we have a change in length
      if (nSamples == transformLength)
      {
         // Check if allocation will result in same step
         if (make_pair(resultLen, resultStep) == calculateResultLenAndStep(_resolution))
         {
            return; // We are good
         }
      }
      delete[]resultRe;
      delete[]resultIm;
      resultLen=0;
      transformLength = 0;
      resultStep = 0;
   }

   assert(resultLen == 0);
   transformLength = nSamples;
   tie(resultLen, resultStep) = calculateResultLenAndStep(_resolution);
   resultRe = new TF_DATA_TYPE[resultLen];
   resultIm = new TF_DATA_TYPE[resultLen];
}
   
int WaveletVoice::transform()
{
   assert(resultLen);
   assert(resultStep);
   
   // Simply step through data until result is full!
   pair<TF_DATA_TYPE *, unsigned int> rval = dyadicFilter->getSamples(octave, -(waveletHalfLen << octave));
   if (rval.second == 0)
   {
      for (int i = resultLen; i--;)
      {
         resultRe[i] = resultIm[i] = 0.0;
      }
      return;
   }

   TF_DATA_TYPE * ptRe = resultRe;
   TF_DATA_TYPE * ptIm = resultIm;
   TF_DATA_TYPE * ptS = rval.first + waveletHalfLen;
   size_t Sstep = resultStep >> octave;
   for (int inx = 0; inx < resultLen; inx++)
   {
      TF_DATA_TYPE * ptR = ptS + 1;
      TF_DATA_TYPE * ptL = ptS - 1;
      TF_DATA_TYPE * ptWRe = waveletRe + 1;
      TF_DATA_TYPE * ptWIm = waveletIm + 1;
      TF_DATA_TYPE sumRe = *waveletRe * *ptS;
      TF_DATA_TYPE sumIm = 0;
      for (int j = waveletHalfLen; j--;)
      {
         sumRe += (*ptL   + *ptR  ) * *ptWRe++;
         sumIm += (*ptL-- - *ptR++) * *ptWIm++;
      }
      *ptRe++ = sumRe;
      *ptIm++ = sumIm;
      ptS +=  Sstep;
      
      // Verify consistency in allocations
      if (inx == 0)
      {
         assert(ptL - rval.first >= -1);
      }
      if (inx == resultLen -1)
      {
         assert(ptR - rval.first <= rval.second);
      }
   }
   return resultLen;
}
   
TF_DATA_TYPE WaveletVoice::get(double timestamp) const
{
   assert(timestamp >= 0 && timestamp < transformLength);
   assert(resultLen);
   unsigned int inx = (unsigned int) (timestamp / resultStep + 0.5);
   assert(inx < resultLen);
      
   // No interpolation, no rounding - just plain get it as easy as possible
   return resultRe[inx] * resultRe[inx] + resultIm[inx] * resultIm[inx];
}
   
