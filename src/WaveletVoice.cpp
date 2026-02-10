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

#include "tft/WaveletVoice.h"
#include <assert.h>
#include <tuple>
#include <cmath>

TFT::WaveletVoice::WaveletVoice(const DyadicFilter * dFilter,
                           double fCenter,
                           double flow,
                           double fhigh) :
    dyadicFilter(dFilter),
    octave(0),
    waveletHalfLength(0),
    resultLen(0),
    duration(0),
    resultStep(0),
    transformLength(0),
    frequency(fCenter),
    fLow(flow),
    fHigh(fhigh),
    resultZ(0),
    voiceSignalBuffer(0)
{
}

void TFT::WaveletVoice::dump() const
{
}

TFT::WaveletVoice::~WaveletVoice()
{
   if (waveletHalfLength)
   {
      delete[]waveletRe;
      delete[]waveletIm;
      delete[]waveletEnvelope;
      waveletHalfLength=0;
   }
   if (hasResultBuffer())
   {
       resultZ.clear();
       resultLen=0;
   }

}
   
void TFT::WaveletVoice::getRequiredPaddingSamples(unsigned int & pre, unsigned int & post) const
{
   // We need WaveletHalfLen samples from the LP filter before and WaveletHalfLen+1 after
   // These numbers must be multiplied by 1<<octave. Additionally, we must add requirement from LP filter
   assert(waveletHalfLength); // Must be initialised by derived constructor
   assert(transformLength);
   assert(resultLen);
   assert(resultStep > 0);
   pre = (waveletHalfLength << octave) + dyadicFilter->getExtraSamples(octave);
   post = pre + resultLen * resultStep - transformLength;
}

std::pair<unsigned int, unsigned int> TFT::WaveletVoice::calculateResultLenAndStep(unsigned int _resolution) const
{
   assert(transformLength);
   unsigned int resultStep = std::max(transformStep, (_resolution / 2) >> octave ) << octave; // Be more conservative internally with resolution
   unsigned int rval = 1 + (transformLength - 1) / resultStep;
   int remaining = (transformLength - 1) - (rval -1) * resultStep;
   while (remaining > 0)
   {
      remaining -= resultStep;
      assert(remaining <= 0);
      rval++;
   }
   return std::pair<unsigned int, unsigned int>(rval, resultStep);
}

void TFT::WaveletVoice::allocateResult(unsigned int nSamples, unsigned int _resolution, bool withBuffer)
{
    valueCache.invalidate();
    // Currently does not support change of buffering scheme for existing object. Not that it could not be done, but this does not match usage patterns
    if (resultLen) {
        // 2nd allocation must not change behavior
        assert(withBuffer == hasResultBuffer());

        // Allocation already done. Get rid of that if we have a change in length
        if (nSamples == transformLength)
        {
            // Check if allocation will result in same step
            if (std::make_pair(resultLen, resultStep) == calculateResultLenAndStep(_resolution))
            {
                return; // We are good
            }
        }
    }

    if (hasResultBuffer())
    {
        resultZ.clear();
    }
    transformLength = nSamples;
    std::tie(resultLen, resultStep) = calculateResultLenAndStep(_resolution);
    if (withBuffer) {
        resultZ.resize(2 * resultLen);
    }
}
   
int TFT::WaveletVoice::transform()
{
    if (hasResultBuffer())
    {
        assert(resultLen);
        assert(resultStep);

        // Simply step through data until result is full!
        std::pair<TF_DATA_TYPE *, unsigned int> rval = dyadicFilter->getSamples(octave, -(waveletHalfLength << octave));
        if (rval.second == 0)
        {
            for (int i = resultLen; i--;)
            {
                resultZ[2 * i] = 0.0;
                resultZ[2 * i + 1] = 0.0;
            }
            return resultLen;
        }

        TF_DATA_TYPE * ptZ = &resultZ[0];
        TF_DATA_TYPE * ptS = rval.first + waveletHalfLength;
        size_t Sstep = resultStep >> octave;
        for (int inx = 0; inx < resultLen; inx++)
        {
            TF_DATA_TYPE * ptR = ptS + 1;
            TF_DATA_TYPE * ptL = ptS - 1;
            TF_DATA_TYPE * ptWRe = waveletRe + 1;
            TF_DATA_TYPE * ptWIm = waveletIm + 1;
            TF_DATA_TYPE sumRe = *waveletRe * *ptS;
            TF_DATA_TYPE sumIm = 0;
            for (int j = waveletHalfLength; j--;)
            {
                sumRe += (*ptL   + *ptR  ) * *ptWRe++;
                sumIm += (*ptL-- - *ptR++) * *ptWIm++;
            }
            *ptZ++ = sumRe;
            *ptZ++ = sumIm;
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
    } else {
       return 0; // We do not do anything when unbuffered
    }
}

TF_DATA_TYPE TFT::WaveletVoice::getBuffered(double timestamp) const
{
    assert(timestamp >= 0 && timestamp < transformLength);
    assert(resultLen);
    assert(hasResultBuffer());
    unsigned int inx = (unsigned int) (timestamp / resultStep + 0.5);
    assert(inx < resultLen);

    // No interpolation, no rounding - just plain get it as easy as possible
    TF_DATA_TYPE val;
    const TF_DATA_TYPE * pZ = &resultZ[2 * inx];
    if (valueCache.lookup(pZ, val))
    {
        return val;
    }
    val = *pZ * *pZ + *(pZ + 1) * *(pZ + 1);
    return valueCache.set(pZ, val);
}

TF_DATA_TYPE TFT::WaveletVoice::getUnbuffered(double timestamp) const
{
   // During getting, we do the calculations
   assert(timestamp >= 0 && timestamp < transformLength);
   assert(resultLen);
   assert(resultStep);
   
   // Simply step through data until result is full!
   std::pair<TF_DATA_TYPE *, unsigned int> rval = dyadicFilter->getSamples(octave, (unsigned int)(timestamp + 0.5) -(waveletHalfLength << octave));
   if (rval.second == 0)
   {
      return 0;
   }
   
   TF_DATA_TYPE val;
   if (valueCache.lookup(rval.first, val))
   {
      return val;
   }

   TF_DATA_TYPE * ptS = rval.first + waveletHalfLength;
   TF_DATA_TYPE * ptR = ptS + 1;
   TF_DATA_TYPE * ptL = ptS - 1;
   TF_DATA_TYPE * ptWRe = waveletRe + 1;
   TF_DATA_TYPE * ptWIm = waveletIm + 1;
   TF_DATA_TYPE sumRe = *waveletRe * *ptS;
   TF_DATA_TYPE sumIm = 0;
   for (int j = waveletHalfLength; j--;)
   {
      sumRe += (*ptL   + *ptR  ) * *ptWRe++;
      sumIm += (*ptL-- - *ptR++) * *ptWIm++;
   }
   assert(ptL - rval.first >= -1);
   assert(ptR - rval.first <= rval.second);
   val = sumRe * sumRe + sumIm * sumIm;
   return valueCache.set(rval.first, val);
}


void TFT::WaveletVoice::executeSequence(int freqStride, int timeStride, TF_DATA_TYPE * out, std::vector<double>::const_iterator timeIterBegin, std::vector<double>::const_iterator timeIterEnd, bool transpose)
{
   transform();
   valueCache.invalidate();
   // We now iterate all timestamps in the sequence and place values in the result buffer. A single timestamp may fill one or more frequencies
   // If not transposed, "out" will contain values for first frequency, then values for second frequency ... finally values for "freqStride'th frequency
   // If transposed, "out" will contain values for first timestamp (length freqStride), then values for second timestamp (length freqStrid) etc
   // Finally, iterate all timestamps
   
   int firstFreqInx = (int)std::ceil(2 * fLow * freqStride);
   int lastFreqInx = (int)std::floor(2 * fHigh * freqStride);
   if (lastFreqInx >= freqStride)
      lastFreqInx = freqStride - 1;
   if (lastFreqInx < firstFreqInx) return;
   if (transpose)
   {
      TF_DATA_TYPE * pOut = out;
      for (auto iterTime = timeIterBegin; iterTime != timeIterEnd; iterTime++)
      {
         assert(*iterTime <= transformLength - 1);
         TF_DATA_TYPE val = getUnbuffered(*iterTime);
         for (auto inx = firstFreqInx; inx <= lastFreqInx; inx++ )
         {
            pOut[inx] = val;
         }
         pOut += freqStride;
      }
   }
   else
   {
      TF_DATA_TYPE * pOut = out;
      for (auto iterTime = timeIterBegin; iterTime != timeIterEnd; iterTime++)
      {
         assert(*iterTime <= transformLength - 1);
         TF_DATA_TYPE val = getUnbuffered(*iterTime);
         for (auto inx = firstFreqInx; inx <= lastFreqInx; inx++ )
         {
            pOut[inx * timeStride] = val;
         }
         pOut++;
      }
   }
}

std::vector<TF_DATA_TYPE> TFT::WaveletVoice::constructVoiceSignal(const std::unique_ptr<IRegion> & region) const {
    // If someone has already constructed a result, we will be using that
    if (!voiceSignalBuffer.empty()) {
        auto rval = voiceSignalBuffer;
        voiceSignalBuffer.clear();
        return rval;
    }

    assert(hasResultBuffer());
    assert(resultStep >> octave == transformStep); // Otherwise we do not fit nicely to dyadic grid
    std::vector<TF_DATA_TYPE> vRe(waveletHalfLength + (resultLen - 1) * transformStep + 1 + waveletHalfLength, 0);
    std::vector<TF_DATA_TYPE> vIm(waveletHalfLength + (resultLen - 1) * transformStep + 1 + waveletHalfLength, 0);
    std::vector<TF_DATA_TYPE> vRePart(waveletHalfLength + 1 + waveletHalfLength);
    std::vector<TF_DATA_TYPE> vImPart(waveletHalfLength + 1 + waveletHalfLength);

    TF_DATA_TYPE omega = 2 * pi * getDecimatedFrequency();

    for (int inx = 0; inx < resultLen; inx++) {
        // Translate coefficients down to DC in order to simplify interpolation
        int iTime = inx * transformStep;

        // Check if we are within region or not
        if (!region->isWithin(iTime << octave, getUndecimatedFrequency())) {
            continue;   // Exclude from summation
        }

        // Phase adjust to current origin of time
        TF_DATA_TYPE delayRe = cos(-omega * iTime);
        TF_DATA_TYPE delayIm = sin(-omega * iTime);

        // Current transform value
        TF_DATA_TYPE valRe = resultZ[2 * inx];
        TF_DATA_TYPE valIm = resultZ[2 * inx + 1];

        // Current factor (complex number)
        TF_DATA_TYPE facRe = valRe * delayRe - delayIm * valIm;
        TF_DATA_TYPE facIm = valRe * delayIm + delayRe * valIm;

        // Fill current signal contribution
        int wInx = waveletHalfLength;
        vRePart[wInx] = waveletEnvelope[0] * facRe;
        vImPart[wInx] = waveletEnvelope[0] * facIm;

        for ( ; wInx; wInx--) {
            vRePart[waveletHalfLength - wInx] = vRePart[waveletHalfLength + wInx] = waveletEnvelope[wInx] * facRe;
            vImPart[waveletHalfLength - wInx] = vImPart[waveletHalfLength + wInx] = waveletEnvelope[wInx] * facIm;
        }

        // And add it to the accumulator
        for (int i = 2 * waveletHalfLength + 1; i--; ) {
            vRe[iTime + i] += vRePart[i];
            vIm[iTime + i] += vImPart[i];
        }
    }

    // extract slice with feed-in
    int too_much = waveletHalfLength - DyadicFilter::getUpsamplingPaddingSize();
    assert(too_much >= 0);
    std::vector<TF_DATA_TYPE> signalRe(vRe.cbegin() + too_much, vRe.cend() - too_much);
    std::vector<TF_DATA_TYPE> signalIm(vIm.cbegin() + too_much, vIm.cend() - too_much);

    // Now we upsample to current sampling rate
    for (int ioctave = octave; ioctave--; ) {
        signalRe = DyadicFilter::upsample(signalRe);
        signalIm = DyadicFilter::upsample(signalIm);
        omega *= 0.5;
    }

    // And finally translate to relevant frequency and form real part
    assert(transformLength <= signalRe.size() - 2 * DyadicFilter::getUpsamplingPaddingSize());
    std::vector<TF_DATA_TYPE> signal(transformLength);
    // The area that we occupy in t/f plane is bw * step
    double normFactor = 2 * (fHigh - fLow) * resultStep / getWaveletEnergy(); // Factor 2 seems reasonable compared to full bandwidth
    for (int i = 0; i < signal.size(); i++) {
        signal[i] = normFactor * (signalRe[DyadicFilter::getUpsamplingPaddingSize() + i] * cos(omega * i) - signalIm[DyadicFilter::getUpsamplingPaddingSize() + i] * sin(omega * i));
    }

    return signal; // Voila!
}

void TFT::WaveletVoice::constructVoiceSignalBuffer(const std::unique_ptr<IRegion> & region) const
{
    voiceSignalBuffer = constructVoiceSignal(region);
}

std::vector<TF_DATA_TYPE> TFT::WaveletVoice::getWavelet() const {
    std::vector<TF_DATA_TYPE> w(DyadicFilter::getUpsamplingPaddingSize() * 2 + waveletHalfLength * 2 + 1, 0);
    for (int i = 0; i <= waveletHalfLength; i++) {
        w[DyadicFilter::getUpsamplingPaddingSize() + waveletHalfLength + i] = w[waveletHalfLength + DyadicFilter::getUpsamplingPaddingSize() - i] = waveletRe[i];
    }

    // Now, interpolate as required
    for (int i = octave; i--; ) {
        w = DyadicFilter::upsample(w);
    }

    return w;
}
