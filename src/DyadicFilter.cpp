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
//  DyadicFilter.cpp
//  tfq
//
//  Created by Klaus Gram-Hansen on 12/11/2025.
//

#include "tft/DyadicFilter.h"
#include <assert.h>
#include <cstring>

const TF_DATA_TYPE TFT::DyadicFilter::filter_taps[DyadicFilter::cstFilterTaps] = {
    0.0003447458782673763,
    0.0006103374853720806,
    -0.0008549219193086368,
    -0.0031822698064080607,
    -0.0012819812679645117,
    0.00459280853085552,
    0.003372387915406647,
    -0.00750276202736566,
    -0.008046301348208626,
    0.010593379830783777,
    0.015720981218986844,
    -0.013907244312892805,
    -0.028387944283777083,
    0.017004235983861287,
    0.05037323690390934,
    -0.01955196418023554,
    -0.09766172462122302,
    0.021225674955957806,
    0.3154139711915178,
    0.4781902537699047,
    0.3154139711915178,
    0.021225674955957806,
    -0.09766172462122302,
    -0.01955196418023554,
    0.05037323690390934,
    0.017004235983861287,
    -0.028387944283777083,
    -0.013907244312892805,
    0.015720981218986844,
    0.010593379830783777,
    -0.008046301348208626,
    -0.00750276202736566,
    0.003372387915406647,
    0.00459280853085552,
    -0.0012819812679645117,
    -0.0031822698064080607,
    -0.0008549219193086368,
    0.0006103374853720806,
    0.0003447458782673763
};

TFT::DyadicFilter::DyadicFilter(const unsigned int nOctaves) : vBufferBegin(nOctaves, NULL), vBufferTimeZero(nOctaves, NULL), vBufferLengths(nOctaves, 0)
{
   assert(nOctaves);
   static_assert(cstFilterTaps & 1, "Filter must have odd length");
}

void TFT::DyadicFilter::doAllocation(unsigned int nSamples, unsigned int nSamplesBefore, unsigned int nSamplesAfter)
{
   // Allocation Already done?
   if (vBufferBegin[0])
   {
      // Same allocation happening?
      unsigned int total_samples_before = vBufferLengths[0];
      long pre_samples_before = vBufferTimeZero[0] - vBufferBegin[0];
      
      if (total_samples_before != nSamplesBefore + nSamples + nSamplesAfter)
      {
         // We need to re-allocate
         // Discard the current buffers
         for (std::vector<TF_DATA_TYPE *>::iterator iter = vBufferBegin.begin(); iter != vBufferBegin.end(); iter++)
         {
            if (*iter)
            {
               delete [] (*iter - 1);
               *iter = NULL;
            }
         }
      }
      else
      {
         // We here assume that nb of pre samples has not changed. Enforce that assumption
         assert(pre_samples_before == nSamplesBefore);
         return;
      }
   }

   // Verify enough samples
   unsigned int padding = getExtraSamples((unsigned int)vBufferBegin.size() - 1);
   assert(nSamplesBefore >= padding);
   assert(nSamplesAfter >= padding);
   // Define some iterators
   std::vector<unsigned int>::iterator iterBufferLength;
   std::vector<TF_DATA_TYPE *>::iterator iterBufferBegin;
   std::vector<TF_DATA_TYPE *>::iterator iterBufferTimeZero;
   for (iterBufferBegin = vBufferBegin.begin(),
        iterBufferLength = vBufferLengths.begin(),
        iterBufferTimeZero = vBufferTimeZero.begin();
        iterBufferBegin != vBufferBegin.end();
        iterBufferBegin++,
        iterBufferLength++,
        iterBufferTimeZero++)
   {
      assert((int)nSamplesBefore >= 0);
      assert((int)nSamplesAfter >= 0);

      int allocationSize = nSamplesBefore + nSamples + nSamplesAfter;
      if (nSamples == 0 || allocationSize <= 0) break;
      *iterBufferBegin = (new TF_DATA_TYPE[allocationSize + 2]) + 1;
      *iterBufferLength = allocationSize;
      *iterBufferTimeZero = *iterBufferBegin + nSamplesBefore;
      
      // Prepare for next octave
      nSamplesBefore = (nSamplesBefore - (cstFilterTaps>> 1)) >> 1;
      nSamplesAfter = (nSamples + nSamplesAfter - (cstFilterTaps>> 1)) >> 1;
      nSamples >>= 1;
      nSamplesAfter -= nSamples;
      
      (*iterBufferBegin)[-1] = cstMagic;
      (*iterBufferBegin)[allocationSize] = cstMagic;
      memset(*iterBufferBegin, 0, sizeof(TF_DATA_TYPE) * allocationSize);
      verify();
   }
}

void TFT::DyadicFilter::verify()
{
   // Define some iterators
   std::vector<unsigned int>::iterator iterBufferLength;
   std::vector<TF_DATA_TYPE *>::iterator iterBufferBegin;
   for (iterBufferBegin = vBufferBegin.begin(),
        iterBufferLength = vBufferLengths.begin();
        iterBufferBegin != vBufferBegin.end();
        iterBufferBegin++,
        iterBufferLength++)
   {
      if (*iterBufferBegin)
      {
         assert(isNaN((*iterBufferBegin)[-1]));
         assert(isNaN((*iterBufferBegin)[*iterBufferLength]));
      }
   }
}

/**
 Find lowest octave that still contains frequencies up to fmax
 */
unsigned int TFT::DyadicFilter::findOctave(const float fmax) const
{
   assert(fmax >= 0);
   float belowF = cstFreqLimit;
   unsigned int octave = 0;
   while (fmax <= belowF && octave < vBufferBegin.size() - 1 )
   {
      octave++;   
      belowF /= 2.0;
   }
   return octave;
}

TFT::DyadicFilter::~DyadicFilter()
{
   for (std::vector<TF_DATA_TYPE *>::iterator iter = vBufferBegin.begin(); iter != vBufferBegin.end(); iter++)
   {
      if (*iter)
      {
         delete [] (*iter - 1);
         *iter = NULL;
      }
   }
}

unsigned int  TFT::DyadicFilter::getExtraSamples(unsigned int octaveIndex)
{

   unsigned int xtra = 0;

   // Iterate octaves
   while (octaveIndex--)
   {
      xtra = 2 * xtra + (cstFilterTaps >> 1);
   }
   return xtra;
}

void TFT::DyadicFilter::filterSamples(const TF_DATA_TYPE * pSamples, const unsigned int n_samples, unsigned int n_samples_before, unsigned int n_samples_after)
{
   // We must have room
   assert(vBufferBegin.size() && vBufferBegin[0]);
   assert(n_samples + n_samples_before + n_samples_after <= vBufferLengths[0]);
   assert(n_samples_before  <= vBufferTimeZero[0] - vBufferBegin[0]);
   assert(n_samples + n_samples_after <= vBufferLengths[0] - (vBufferTimeZero[0] - vBufferBegin[0]));
   
   // Take it all in
   // First octave will be initially zeroed in order to provide zeropadding
   memset(vBufferBegin[0],0, sizeof(TF_DATA_TYPE) * vBufferLengths[0]);
   std::vector<TF_DATA_TYPE *>::iterator dst = vBufferBegin.begin();
   std::vector<TF_DATA_TYPE *>::iterator dst0 = vBufferTimeZero.begin();
   std::vector<unsigned int>::iterator countIter = vBufferLengths.begin();
   memcpy(*dst0 - n_samples_before, pSamples - n_samples_before, sizeof(TF_DATA_TYPE) * (n_samples_before + n_samples + n_samples_after));

   TF_DATA_TYPE * src = *dst;
   // Adjust for even number of pre-samples
   static_assert((cstFilterTaps >> 1) & 1, "Odd number of half-filter taps expected"); // If this was not the case, we would have to check for odd number of pre-samples (also below: two places)
   if (((*dst0 - *dst) & 1) == 0)
   {
      src++;
   }
   while (++dst != vBufferBegin.end() && *dst)
   {
      dst0++;
      countIter++;
      TF_DATA_TYPE * ptrIn = src;
      TF_DATA_TYPE * ptrOut = *dst;
      if (*countIter == 0) break;
      for (unsigned int i = *countIter; i--;)
      {
         TF_DATA_TYPE * ptrSignal = ptrIn;
         const TF_DATA_TYPE * ptrFilter = filter_taps;
         double sum = 0.0;
         for (int j = cstFilterTaps; j--;)
         {
            sum += *ptrSignal++ * *ptrFilter++;
         }
         *ptrOut++ = sum;
         ptrIn+=2;
      }
      src = *dst;
      // Adjust for even number of pre-samples
      if (((*dst0 - *dst) & 1) == 0)
      {
         src++;
      }
      verify();
   }
}

/**
 @param nSamples will be filled with number of samples in range from last "insertSamples" operation for the requested octave
 Note that it is ok to request from a negative sample or even a sample beyound number of samples fed in
 We will always get the sample which is clostst to the requested one
 */
const std::pair<TF_DATA_TYPE *, unsigned int> TFT::DyadicFilter::getSamples(const unsigned int octave, const int fromSample) const
{
   assert(octave < vBufferBegin.size());
   if (vBufferBegin[octave] == NULL)
   {
      return std::pair<TF_DATA_TYPE *, unsigned int>(NULL, 0);
   }
   int offsetFromZero = (octave == 0) ? fromSample : (((fromSample >> (octave - 1)) + 1) >> 1); // This is rounding!
   TF_DATA_TYPE * from = vBufferTimeZero[octave] + offsetFromZero;
   long offsetFromStart = from - vBufferBegin[octave];
   unsigned long nb = vBufferLengths[octave] - offsetFromStart;
   assert(from >= vBufferBegin[octave]);
   assert(from < vBufferBegin[octave] + vBufferLengths[octave]);
   assert(nb);
   return std::pair<TF_DATA_TYPE *, unsigned int>(from, nb);
}

