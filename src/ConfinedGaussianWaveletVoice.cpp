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
//  WaveletVoice.h
//  test
//
//  Created by Klaus Gram-Hansen on 17/11/2025.
//

#include <assert.h>
#include <cmath>
#include <iostream>
#include "tft/ConfinedGaussianWaveletVoice.h"
using namespace std;

TFT::ConfinedGaussianWaveletVoice::ConfinedGaussianWaveletVoice(double fCenter,
                double fDelay,
                double flow,
                double fhigh,
                double Q,
                const float overlapPercentage,
                const DyadicFilter * dFilter) : WaveletVoice(dFilter, fCenter, flow, fhigh), q(Q), ConfinedGaussian(sigma)
{
   // When constructing, firstly let us know about the upper frequency limit by the wavelet to be defined
   assert(fCenter > 0 && fCenter <= 0.5);
   assert(dFilter);
   assert(overlapPercentage < 100.0);
   assert(Q > 1); // Otherwise wavelet simply gets too degenerated
   float deltaF = fCenter / Q / 2.0;
   float fmax = bw70Ratio * deltaF + fCenter;
   
   // From fmax, locate the octave we will be using
   octave = dFilter->findOctave(fmax);
   
   // recalculate fCenter relative to actual octave
   fCenter *= (1 << octave);
   fDelay /= (1 << octave);
   
   // Calculate Len of wavelet, halfLen = (1 / (2 * pi * sigma * fCenter / Q)) / 2
   waveletHalfLength = (unsigned int) ((Q / (2.0 * pi * sigma * fCenter)) / 2 );
   unsigned int waveletLength = 2 * waveletHalfLength + 1;
   setLength(waveletLength);
   waveletRe = new TF_DATA_TYPE[waveletHalfLength + 1];
   waveletIm = new TF_DATA_TYPE[waveletHalfLength + 1];
   waveletEnvelope = new TF_DATA_TYPE[waveletHalfLength + 1];

   TF_DATA_TYPE sum = waveletEnvelope[0] = approximateConfinedGaussian(waveletHalfLength); // Assure we always have the same SUM (fDelay parameter is only in use for normalisation)
   waveletEnergy = waveletEnvelope[0] * waveletEnvelope[0];
   TF_DATA_TYPE argument = 2 * pi * fCenter * fDelay;
   waveletRe[0] = waveletEnvelope[0] * cos(argument);
   waveletIm[0] = waveletEnvelope[0] * sin(argument);
   
   for (int i = 1; i <= waveletHalfLength; i++)
   {
      //     wavelet = np.append(wavelet, env * math.cos(2 * math.pi * fc * xs[inx]))
      argument = 2 * pi * fCenter * (i + fDelay);
      TF_DATA_TYPE envelope = approximateConfinedGaussian(waveletHalfLength + i);
      sum += 2 * envelope;
      waveletEnergy += 2 * envelope * envelope;
      waveletRe[i] = envelope * cos(argument);
      waveletIm[i] = envelope * sin(argument);
      waveletEnvelope[i] = envelope;
   }
   TF_DATA_TYPE factor = 2 / sum; // Factor "2" bcs of one-sided spectrum
   waveletEnergy *= factor * factor;
   
   for (int i = 0; i <= waveletHalfLength; i++)
   {
      waveletRe[i] *= factor;
      waveletIm[i] *= factor;
      waveletEnvelope[i] *= factor;
   }
   
   // Calculate duration (measured in samples at current rate
   // Approximate temporal resolution is dT =  (2 * halfLen + 1) * sigma
   duration = 2.0 * waveletLength *sigma;
   float step = duration * (100-overlapPercentage) / 100.0;

   // Select power-of-two just larger. Not chosing any integer is due to inverse transform
   int exponent = (int)floor(log2(step));
   if (exponent < 0)
   {
       exponent = 0;
   }
   transformStep = 1 << exponent;
}


void TFT::ConfinedGaussianWaveletVoice::dump() const
{
   cout << "Approximate Confined Gaussian Wavelet in octave " << octave << ", duration " << duration << ", step " << transformStep << "[ " << fLow << "; " << getUndecimatedFrequency() << "; " << fHigh << " ]" << endl;
   WaveletVoice::dump();
}

TFT::ConfinedGaussianWaveletVoice::~ConfinedGaussianWaveletVoice()
{
   
}
