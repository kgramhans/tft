# _tft - Time/Frequency Transform library_
tft is a library for performing transformation of signals (1D) into the time/frequency domain (2D). Also, backward transform from time/frequency domain to pure time domain is offered - allowing for masking out parts of the signal in the time/frequency domain.
An infinite number of possible time/frequency transformations exists due to the inherent time/frequency duality. It is the aim of the current library to provide a uniform interface for implementing any such possible transform.
Focus is on implementing time+FREQUENCY representations of signals as applicable to acoustic signals, as opposed to time+SCALE or time+AR (AR: Auto regressive) representations. This is so because frequency represents an easily interpretable quantity in acoustics.
Examples of Time-Frequency representations are:
- Short-Time Fourier Transform (STFT). Also known as a Spectrogram
- Continuous Wavelet transform
- Wigner-Ville distribution

## Features
- Generic interfaces
- 2^n transform length STFT with Gaussian Window
- Continuous Wavelet Transform with selectable Q factor
- Designed for parallel processing
- Control of transform overlap
- Specific interface optimised for coarse representations (aka Time/Frequency Calculator)
- Specific interface optimised for time/frequency modifications and inverse transform (aka Time/Frequency Transformer)

## Generic Interfaces
### TFT::ITimeFrequencyCalculator
This interface is meant to allow for efficiency in computation (only calculate what is necessary) and memory usage (calculations happen when extracting values)
| Member | Explanation |
|--------|-------|
|prepare()| Interrogate transform in order to obtain information about required context|
|doTransform()| Feed samples into the calculator|
|extractFrequencySlices()|Extract values from the transform|
|prepareSequences()|Set up Calculator for parallel sequences of doTransform()+extractFrequencySlices() operations|
|executeSequence()|Execute such sequences set up above (parallel execution)

### TFT::ITimeFrequencyTransformer
This interface is meant to actually transform a time signal into a time/frequency representation and back. When transforming back, filtering in the time/frequency plan can take place. This feature allows for extraction of signatures or cancelling of unwanted signatures. Compared to the TFT::ITimeFrequencyCalculator interface more memory is needed since a proper representation of the signal in the time/frequency plane is required.
| Member | Explanation |
|--------|-------|
|prepare()| Interrogate transform in order to obtain information about required context|
|forwardTransform()| Generate a time/frequency representation of a signal |
|setPolygonRegion()|Define a polygon to be used by backwardTransform()|
|backwardTransform()|Transform the time/frequency representation back to a signal. The filter region will be applied in the process|

As for the TFT::ITimeFrequencyCalculator interface, prepare/execute members also do exist that allow for parallelisation of both the forward and backward transforms:
| Member | Explanation |
|--------|-------|
|prepareParallelForwardSequences()| Set up for parallel execution of forward Transformation|
|executeForwardSequence()| Execute such sequences set up above |
|prepareParallelBackwardSequences()|Set up for parallel execution of backward Transformation|
|executeBackwardSequence()|Execute such sequences set up above|
## Specific implementations of Interfaces
### WaveletCalculator
Implements TFT::ITimeFrequencyCalculator Interface, allowing for a configurable overlap, a number of frequency octaves and a Q factor of the analysis. The analysing kernel for the transform is a so-called Confined Gaussian Wavelet based on a windowing function described [here](https://en.wikipedia.org/wiki/Window_function#Confined_Gaussian_window). The Confined Gaussian is the most favourable function when it comes to simultaneous time and frequency resolution. The WaveletCalculator will produce a time/frequency representation with constant Q: Short duration/broad bandwidth at high frequencies combined with Long duration/narrow bandwidth at low frequencies
### StftCalculator
Implements TFT::ITimeFrequencyCalculator Interface, allowing for a configurable overlap and a fixed window length (2^n). The analysing kernel for the transform defined by an external window function passed to the implementation. The StftCalculator will produce a time/frequency representation with constant bandwidth and constant time resolution. A multitude of variants are made possible via the external windowing function
### WaveletTransformer
Implements TFT::ITimeFrequencyTransformer Interface with support for region filtering in the time/frequency plane. Like the WaveletCalculator, the implementation is based on the Confined Gaussian Wavelet with constant Q
### StftTransformer
Implements TFT::ITimeFrequencyTransformer Interface with support for region filtering in the time/frequency plane. The implementation is based on the Confined Gaussian Window with constant duration (no external window)
