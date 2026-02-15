#include <gtest/gtest.h>
#include <tft/tft.h>
#include <tft/RealFFTf.h>
#include <tft/ConfinedGaussianWaveletVoice.h>
#include <random>
#include <numeric>

using namespace TFT;
using namespace std;


// Test WaveletVoice
TEST(WaveletVoiceBuffered, GeneralConfiguration) {
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
        w.allocateResult(n_samples, 0, true);
        w.getRequiredPaddingSamples(pre, post);
    }
    for (int i = n_octaves - 1; i--;)
    {
        unsigned int _pre, _post;
        frequency /= 2;
        ConfinedGaussianWaveletVoice w(frequency, frequency - 0.1, frequency + 0.2, Q, overlap, &filter);
        w.allocateResult(n_samples, 0, true);
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
        w.allocateResult(n_samples, 0, true);
        w.getRequiredPaddingSamples(pre, post);
    }
    for (int i = 10; i--;)
    {
        unsigned int _pre, _post;
        Q *= 2;
        ConfinedGaussianWaveletVoice w(frequency, frequency - 0.1, frequency + 0.2, Q, overlap, &filter);
        w.allocateResult(n_samples, 0, true);
        w.getRequiredPaddingSamples(_pre, _post);
        EXPECT_GE(_pre, 2 * pre) << "Q : " << Q;
        EXPECT_LE(_pre, 2 * pre + 1) << "Q : " << Q;

        EXPECT_GE(_post, 2 * post) << "Q : " << Q;
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
        w.allocateResult(0, 1, true);
    }
        , "Assertion");
}

TEST(WaveletVoiceBuffered, Calculation) {
    double Q = 10.0;
    double overlap = 99.9; // For this get as close as possible to continuous transform
    int n_octaves = 1;       // For simplicity have just a single octave
    DyadicFilter filter(n_octaves);
    double frequency = 0.2;

    ConfinedGaussianWaveletVoice w(frequency, 0, 0.5, Q, overlap, &filter);

    // We proceed with test case by supressing a delta function
    // This should return the impulse response of the wavelet that we can then use for testing properties of the wavelet
    int n_samples = 1024; // Assumed sufficiently large to contain the wavelet duration
    unsigned int pre, post;

    // Mimic allocation scheme already implemented in WaveletCalculator
    w.allocateResult(n_samples, 1, true);
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
    w.transform();
    vector<double> vTime(n_samples);
    iota(vTime.begin(), vTime.end(), 0);
    int inx = 0;
    for (auto iter = vTime.begin(); iter != vTime.end(); iter++) {
        out[inx++] = w.get(*iter);
    }

    // Having gotten access to wavelet impulse response we can now
    // Verify in time domain
    double duration_sum = 0;
    double weight_sum = 0;
    for (int i = 0; i < n_samples; i++) {
        duration_sum += i * sqrt(out[i]);
        weight_sum += sqrt(out[i]);
    }
    double center_time = duration_sum / weight_sum;
    EXPECT_EQ((int)(center_time + 0.5), n_samples / 2) << center_time;

    double duration_sqr_sum = 0;
    weight_sum = 0.0;
    for (int i = 0; i < n_samples; i++) {
        duration_sqr_sum += (i - center_time) *  (i - center_time) * out[i];
        weight_sum += out[i];
    }
    double duration = sqrt(duration_sqr_sum / weight_sum);
    double theoretical_duration = Q / 2 / 4 / atan(1) / frequency;    // Theoretical RMS duration is: 1/Q = 2*dF/relF = 1/2/pi/dT/relF ==> dT = Q/(2*pi*relF) = 10 / (2*pi*0.2) = 7.96

    EXPECT_LE(abs(duration - theoretical_duration) / theoretical_duration, 0.01) << "RMS Duration " << duration << " (" << theoretical_duration << ")";


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
    EXPECT_LE(abs(bw - theoretical_bw) / theoretical_bw, 0.015) << "bw " << bw << ", theoretical bw " << theoretical_bw;

    EXPECT_LE(bw * duration * 4 * 4 * atan(1), 1.01);

    // Verify also, that level is below 70dB of the level at 0 when we are more than 6.5 times away
    auto from_i = (int)ceil(6.5 * theoretical_bw * n_samples);
    for (auto i = from_i; i < n_samples / 2; i++) {
        EXPECT_LE(10 * log10(power[i]), 10 * log10(power[0]) - 70) << i << from_i << " " << 6.5 * theoretical_bw * n_samples ;
    }
}

TEST(WaveletVoiceBuffered, Transformation) {
    double Q = 20.0; // Relatively narrow band
    double overlap = 75; // For this get as close as possible to continuous transform
    int n_octaves = 10;
    DyadicFilter filter(n_octaves);
    int n_samples = 1024; // Assumed sufficiently large to contain the wavelet duration
    double frequency = 100.0 / n_samples;

    ConfinedGaussianWaveletVoice w(frequency, 0, 0.5, Q, overlap, &filter);

    // We proceed with test case by supressing a delta function
    // This should return the prolonged impulse response of the wavelet when restoring
    unsigned int pre = 10000, post = 10000;

    // Mimic allocation scheme already implemented in WaveletCalculator
    w.allocateResult(n_samples, 1, true);
    filter.doAllocation(n_samples, pre, post);

    // Allocate our signal array and fill in delta
    vector<TF_DATA_TYPE> signal(pre + n_samples + post, 0);
    signal[pre + n_samples / 2] = 1e6;

    // Dyadic filtration makes us ready
    filter.filterSamples(&signal[pre], n_samples, pre, post);

    // After all this synopsis we ar now ready for the actual wavelet filtering operation that should reveal its filter response
    // We re-use the time signal array herefor
    w.transform();

    auto v = w.constructVoiceSignal(std::make_unique<PolygonRegionHorizontal>());

    // We expect a signal centered around sample 512 and frequency 0.1 and which has BWrms =  0.1 / 20 / 2 = 0.0025 and rms duration exceeding 1 / (4 * pi * BWrms) = 31.8 (by how much?)
    // First time domain
    double bwRms = frequency / Q / 2;
    double durRms = 1.0 / 4.0 / (4.0*atan(1)) / bwRms;
    double duration_sum = 0;
    double weight_sum = 0;
    for (int i = 0; i < n_samples; i++) {
        duration_sum += i * abs(v[i]);
        weight_sum += abs(v[i]);
    }
    double center_time = duration_sum / weight_sum;
    EXPECT_EQ((int)(center_time + 0.5), n_samples / 2) << center_time;

    double duration_sqr_sum = 0;
    weight_sum = 0.0;
    for (int i = 0; i < n_samples; i++) {
        duration_sqr_sum += (i - center_time) *  (i - center_time) * v[i] * v[i];
        weight_sum += v[i] * v[i];
    }
    double duration = sqrt(duration_sqr_sum / weight_sum);

    EXPECT_GT(duration, durRms);
    EXPECT_LT(duration, 2 * durRms);

    // Verify in frequency domain
    // Do the DFT thing
    auto handle = TFT::GetFFT(n_samples);
    TFT::RealFFTf(&v[0], handle.get());

    // Convert to magnitude
    vector<TF_DATA_TYPE> mag(handle->Points);
    mag[0] = sqrt(v[0] * v[0]);

    for (int i = 1; i < handle->Points; i++) {
        const int index = handle->BitReversed[i];
        const float re = v[index], im = v[index + 1];
        mag[i] = sqrt(re * re + im * im);
    }

    double bw_sum = 0;
    weight_sum = 0;
    for (int i = 0; i < handle->Points; i++) {
        bw_sum += i * mag[i];
        weight_sum += mag[i];
    }
    double center_freq = bw_sum / weight_sum;
    EXPECT_EQ((int)(center_freq + 0.5), (int)(n_samples * frequency + 0.5)) << center_freq;

    double bw_sqr_sum = 0;
    weight_sum = 0;
    for (int i = 0; i < handle->Points; i++) {
        bw_sqr_sum += (i - center_freq) * (i - center_freq) * mag[i] * mag[i];
        weight_sum += mag[i] * mag[i];
    }
    double bw = sqrt(bw_sqr_sum / weight_sum) / n_samples;
    double theoretical_bw = frequency / Q / 2;
    EXPECT_LT(bw, theoretical_bw); // Sharper bcs it is longer

    double heisenbergProduct = theoretical_bw * duration * 4 * 4 * atan(1);

    // This is reasonably located in t/f if equal to sqrt(2)
    EXPECT_LT(heisenbergProduct, sqrt(2) * 1.0001);
}
