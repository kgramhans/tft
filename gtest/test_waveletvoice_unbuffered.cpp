#include <gtest/gtest.h>
#include <tft/tft.h>
#include <tft/RealFFTf.h>
#include <tft/ConfinedGaussianWaveletVoice.h>
#include <random>
#include <numeric>

using namespace TFT;
using namespace std;


// Test WaveletVoice
TEST(WaveletVoiceUnbuffered, ConfinedGaussianConstruction) {
    double Q = 10.0;
    double overlap = 50.0;
    int n_octaves = 10;
    DyadicFilter filter(n_octaves);

    double frequency;

    frequency = -1.0;
    EXPECT_DEBUG_DEATH(ConfinedGaussianWaveletVoice w(frequency, frequency - 0.1, frequency + 0.2, Q, overlap, &filter), "Assertion") << "Freq " << frequency;

    frequency = 0.0;
    EXPECT_DEBUG_DEATH(ConfinedGaussianWaveletVoice w(frequency, frequency - 0.1, frequency + 0.2, Q, overlap, &filter), "Assertion") << "Freq " << frequency;

    frequency = 0.6;
    EXPECT_DEBUG_DEATH(ConfinedGaussianWaveletVoice w(frequency, frequency - 0.1, frequency + 0.2, Q, overlap, &filter), "Assertion") << "Freq " << frequency;

    frequency = 0.25;
    EXPECT_DEBUG_DEATH(ConfinedGaussianWaveletVoice w(frequency, frequency - 0.1, frequency + 0.2, Q, overlap, NULL), "Assertion") << "Filter is NULL";

    overlap = 110;
    EXPECT_DEBUG_DEATH(ConfinedGaussianWaveletVoice w(frequency, frequency - 0.1, frequency + 0.2, Q, overlap, &filter), "Assertion") << "Overlap is " << overlap;

    overlap = 50;
    Q = 0;
    EXPECT_DEBUG_DEATH(ConfinedGaussianWaveletVoice w(frequency, frequency - 0.1, frequency + 0.2, Q, overlap, &filter), "Assertion") << "Q is zero";

    // Verify that we can create across a broad range of frequencies (just hammer it with randomness)
    // Only criterion here is that we do not die
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist2500(1,2500);
    std::uniform_int_distribution<std::mt19937::result_type> dist100(1,100);

    int runs = 10000;
    while (runs--) {
        frequency = dist2500(rng) / 10000.0;
        Q = dist100(rng);
        if (Q > 1) {
            EXPECT_NO_THROW(
                ConfinedGaussianWaveletVoice w(frequency, frequency - 0.1, frequency + 0.2, Q, overlap, &filter);
                ) << "Q: " << Q << ", f; " << frequency;
        } else {
            EXPECT_DEBUG_DEATH(
                ConfinedGaussianWaveletVoice w(frequency, frequency - 0.1, frequency + 0.2, Q, overlap, &filter), "Assertion"
                ) << "Q: " << Q << ", f; " << frequency;
        }
    }

}

TEST(WaveletVoiceUnbuffered, GeneralConfiguration) {
    double Q = 10.0;
    double overlap = 99.999; // For this get as close as possible to continuous transform
    int n_octaves = 10;
    DyadicFilter filter(n_octaves);
    unsigned int pre, post;
    double frequency = 0.25;
    unsigned int n_samples = 1 << (n_octaves + 4);

    // Verify that we cannot proceed before "allocation" was done
    EXPECT_DEBUG_DEATH(
    {
        frequency = 0.2;
        ConfinedGaussianWaveletVoice w(frequency, frequency - 0.1, frequency + 0.2, Q, overlap, &filter);
        w.getRequiredPaddingSamples(pre, post);
    }
        , "Assertion")<< frequency << ", " << pre << ", " << post << endl;

    // Verify that padding samples grows approximately linearly with factor between 2 and 3 with lower frequency
    {
        frequency = 0.2;
        ConfinedGaussianWaveletVoice w(frequency, frequency - 0.1, frequency + 0.2, Q, overlap, &filter);
        w.allocateResult(n_samples, 0, false);
        w.getRequiredPaddingSamples(pre, post);
    }
    for (int i = n_octaves - 1; i--;)
    {
        unsigned int _pre, _post;
        frequency /= 2;
        ConfinedGaussianWaveletVoice w(frequency, frequency - 0.1, frequency + 0.2, Q, overlap, &filter);
        w.allocateResult(n_samples, 0, false);
        w.getRequiredPaddingSamples(_pre, _post);
        EXPECT_GT(_pre, 2 * pre) << "Frequency : " << frequency;
        EXPECT_LT(_pre, 3 * pre) << "Frequency : " << frequency;

        EXPECT_GT(_post, 2 * post) << "Frequency : " << frequency;
        EXPECT_LT(_post, 3 * post) << "Frequency : " << frequency;
        pre = _pre;
        post = _post;
    }

    // Verify that padding samples grows linearly with Q for a wavelet in first octave
    {
        frequency = 0.2;
        Q = 2;
        ConfinedGaussianWaveletVoice w(frequency, frequency - 0.1, frequency + 0.2, Q, overlap, &filter);
        w.allocateResult(n_samples, 0, false);
        w.getRequiredPaddingSamples(pre, post);
    }
    for (int i = 10; i--;)
    {
        unsigned int _pre, _post;
        Q *= 2;
        ConfinedGaussianWaveletVoice w(frequency, frequency - 0.1, frequency + 0.2, Q, overlap, &filter);
        w.allocateResult(n_samples, 0, false);
        w.getRequiredPaddingSamples(_pre, _post);
        EXPECT_GE(_pre, 2 * pre - 1) << "Q : " << Q;  // Take into account rounding effects
        EXPECT_LE(_pre, 2 * pre + 1) << "Q : " << Q;

        EXPECT_GE(_post, 2 * post - 1) << "Q : " << Q; // Take into account rounding effects
        EXPECT_LE(_post, 2 * post + 1) << "Q : " << Q;
        pre = _pre;
        post = _post;
    }

    // Test that frequency contained is remembered
    frequency = 0.2;
    Q = 5;
    {
        double delta = 0.001;
        ConfinedGaussianWaveletVoice w(frequency, frequency - 0.1, frequency + 0.2, Q, overlap, &filter);
        EXPECT_FALSE(w.containsFrequency(frequency - 0.100 - delta));
        EXPECT_TRUE(w.containsFrequency(frequency - 0.100)); // Lower bound included
        EXPECT_TRUE(w.containsFrequency(frequency - 0.100 + delta));
        EXPECT_TRUE(w.containsFrequency(frequency));
        EXPECT_TRUE(w.containsFrequency(frequency + 0.200 - delta));
        EXPECT_FALSE(w.containsFrequency(frequency + 0.200));
        EXPECT_FALSE(w.containsFrequency(frequency + 0.200 + delta));
    }

    // Test allocation
    // Nothing much we can do here except test that illegal value is rejected (just a single illegal value)
    EXPECT_DEBUG_DEATH(
    {
        ConfinedGaussianWaveletVoice w(frequency, frequency - 0.1, frequency + 0.2, Q, overlap, &filter);
        w.allocateResult(0, 1, false);
    }
        , "Assertion");
}

TEST(WaveletVoiceUnbuffered, Transformation) {
    double Q = 10.0;
    double overlap = 99.999; // For this get as close as possible to continuous transform
    int n_octaves = 1;       // For simplicity have just a single octave
    DyadicFilter filter(n_octaves);
    double frequency = 0.2;

    ConfinedGaussianWaveletVoice w(frequency, 0, 0.5, Q, overlap, &filter);

    // Depending on Baseclass, we test differently:
    // * WaveletVoiceUnbuffered only supports sequenced operation: Transform does nothing
    // * WaveletVoice does support buffered transform
    // Currently only support for the former has been implemented test-wise

    // We proceed with test case by supressing a delta function
    // This should return the impulse response of the wavelet that we can then use for testing properties of the wavelet
    int n_samples = 1024; // Assumed sufficiently large to contain the wavelet duration
    unsigned int pre, post;

    // Mimic allocation scheme already implemented in WaveletCalculator
    w.allocateResult(n_samples, 1, false);
    w.getRequiredPaddingSamples(pre, post);
    filter.doAllocation(n_samples, pre, post);

    // Allocate our signal array and fill in delta
    vector<TF_DATA_TYPE> signal(pre + n_samples + post, 0);
    signal[pre + n_samples / 2] = 1;

    // Allocate output
    vector<TF_DATA_TYPE> out(n_samples);

    // Dyadic filtration makes us ready
    filter.filterSamples(&signal[pre], n_samples, pre, post);

    // After all this synopsis we ar now ready for the actual wavelet filtering operation that should reveal its filter response
    // We re-use the time signal array herefor
    vector<double> vTime(n_samples);
    iota(vTime.begin(), vTime.end(), 0);
    w.executeSequence(1, 1, &out[0], vTime.cbegin(), vTime.cend(), false);

    // Having gotten access to wavelet impulse response we can now
    // Verify in time domain
    double duration_sum = 0;
    double weight_sum = 0;
    for (int i = 0; i < n_samples; i++) {
        duration_sum += i * sqrt(out[i]);
        weight_sum += sqrt(out[i]);
    }
    double center_time = duration_sum / weight_sum;
    EXPECT_EQ((int)(center_time + 0.5), n_samples / 2);

    double duration_sqr_sum = 0;
    weight_sum = 0.0;
    for (int i = 0; i < n_samples; i++) {
        duration_sqr_sum += (i - center_time) *  (i - center_time) * out[i];
        weight_sum += out[i];
    }
    double duration = sqrt(duration_sqr_sum / weight_sum);
    double theoretical_duration = Q / 2 / 4 / atan(1) / frequency;    // Theoretical RMS duration is: 1/Q = 2*dF/relF = 1/2/pi/dT/relF ==> dT = Q/(2*pi*relF) = 10 / (2*pi*0.2) = 7.96
    double sigma = 0.1;
    double precision = 1/(theoretical_duration / sigma);
    EXPECT_LE(abs(duration - theoretical_duration) / theoretical_duration, precision) << "RMS Duration " << duration << " (" << theoretical_duration << ")";


    // Verify in frequency domain
    for (int i = n_samples ; i --; ) {
        out[i] = sqrt(out[i]);  // Prepare for FFT of absolute values
    }

    // Do the DFT thing
    auto handle = TFT::GetFFT(n_samples);
    TFT::RealFFTf(&out[0], handle.get());

    // Convert to magnitude
    vector<TF_DATA_TYPE> power(handle->Points);
    power[0] = out[0] * out[0];

    for (int i = 1; i < handle->Points; i++) {
        const int index = handle->BitReversed[i];
        const float re = out[index], im = out[index + 1];
        power[i] = re * re + im * im;
    }

    double bw_sqr_sum = 0;
    weight_sum = 0;
    for (int i = 0; i < handle->Points; i++) {
        bw_sqr_sum += i * i * power[i];
        weight_sum += power[i];
    }
    double bw = sqrt(bw_sqr_sum / weight_sum) / n_samples;
    double theoretical_bw = frequency / Q / 2;
    EXPECT_LE(abs(bw - theoretical_bw) / theoretical_bw, precision) << "bw " << bw << ", theoretical bw " << theoretical_bw;

    EXPECT_LE(bw * duration * 4 * 4 * atan(1), 1.01);

    // Verify also, that level is below 70dB of the level at 0 when we are more than 6.5 times away
    auto from_i = (int)ceil(6.5 * theoretical_bw * n_samples);
    for (auto i = from_i; i < n_samples / 2; i++) {
        EXPECT_LE(10 * log10(power[i]), 10 * log10(power[0]) - 70) << i << from_i << " " << 6.5 * theoretical_bw * n_samples ;
    }


}
