#include <gtest/gtest.h>
#include <tft/tft.h>
#include <iomanip>

using namespace TFT;
using namespace std;

// Test Dyadic filter
TEST(DyadicFilter, BasicConfiguration) {
    EXPECT_DEBUG_DEATH(DyadicFilter filter(0), "Assertion");
    constexpr int n_octaves = 10;
    DyadicFilter filter(n_octaves);

    EXPECT_DEBUG_DEATH(filter.findOctave(-1.0), "Assertion");
    EXPECT_EQ(0, filter.findOctave(0.5));
    EXPECT_EQ(0, filter.findOctave(1.0));
    EXPECT_EQ(0, filter.findOctave(0.21));
    EXPECT_EQ(1, filter.findOctave(0.19));
    EXPECT_EQ(1, filter.findOctave(0.21/2));
    EXPECT_EQ(2, filter.findOctave(0.19/2));
    EXPECT_EQ(n_octaves - 1, filter.findOctave(0.00001));
    EXPECT_EQ(n_octaves - 1, filter.findOctave(0.0));

    EXPECT_EQ(0, DyadicFilter::getExtraSamples(0));
    EXPECT_LT(0, DyadicFilter::getExtraSamples(1));
    EXPECT_EQ(DyadicFilter::getExtraSamples(2), 3*DyadicFilter::getExtraSamples(1));
    EXPECT_EQ(DyadicFilter::getExtraSamples(3), 7*DyadicFilter::getExtraSamples(1));
    EXPECT_EQ(DyadicFilter::getExtraSamples(4), 15*DyadicFilter::getExtraSamples(1));
    EXPECT_EQ(DyadicFilter::getExtraSamples(5), 31*DyadicFilter::getExtraSamples(1));
    EXPECT_EQ(DyadicFilter::getExtraSamples(20), ((1<<20) - 1)*DyadicFilter::getExtraSamples(1));

    int xtra = DyadicFilter::getExtraSamples(n_octaves - 1);
    int n_samples = 1000 + 2 * xtra;
    vector<TF_DATA_TYPE> vInput(n_samples, 0);
    // Cannot filter without allocation
    EXPECT_DEBUG_DEATH(filter.filterSamples(&vInput[xtra], n_samples, xtra, xtra), "Assertion");

    // But can filter once allocatd
    EXPECT_DEBUG_DEATH(filter.doAllocation(n_samples, 0, 0), "Assertion");       // This would fail as long as we do not provide enough samples
    EXPECT_DEBUG_DEATH(filter.doAllocation(n_samples, xtra - 1, xtra - 1), "Assertion");       // This would fail as long as we do not provide enough samples

    // Test that we can allocat
    filter.doAllocation(n_samples, xtra, xtra);
    filter.doAllocation(n_samples, xtra, xtra);
    filter.doAllocation(n_samples, 2*xtra, 2*xtra);

    // Test that filtering succeeds. Do not care about result
    filter.filterSamples(&vInput[xtra], n_samples, xtra, xtra);
}

// Test shape of impulse response qualitatively
TEST(DyadicFilter, ImpulseResponse) {
    constexpr int n_octaves = 10;
    DyadicFilter filter(n_octaves);
    int xtra = DyadicFilter::getExtraSamples(n_octaves - 1);
    int n_samples = DyadicFilter::getExtraSamples(1)* (1 << (n_octaves - 1));
    vector<TF_DATA_TYPE> vInput(n_samples + 2 * xtra, 0);
    vInput[xtra] = 1;  // Impulse response
    filter.doAllocation(n_samples, xtra, xtra);
    filter.filterSamples(&vInput[xtra], n_samples, xtra, xtra);

    // In verifying the result of filtering a delta function we check
    // 1st sample will drop by a factor 0.45 and 0.5 for each octave
    // Every octave output will have 4 fast low-level (<10%) diminusing oscillations

    TF_DATA_TYPE level = 1;
    for (int i = 0; i < n_octaves; i++)
    {
        TF_DATA_TYPE * ptr;
        int nb;
        tie(ptr, nb) = filter.getSamples(i, 0);

        EXPECT_GE(nb << i, n_samples); // Must give us enough as configured


        EXPECT_GT(*ptr, 0.90 * level);
        EXPECT_LT(*ptr, 1.02 * level);

        TF_DATA_TYPE power = *ptr * *ptr;
        double sum_duration = 0;
        for (int i = 1; i < nb; i++)
        {
            power += 2 * ptr[i] * ptr[i];
            sum_duration += 2 * i * i * ptr[i] * ptr[i];
        }
        double duration = sqrt(sum_duration / power);
        EXPECT_LE(duration, 0.5); // verify response width in order to assure it does not spread out due to dispersion

        EXPECT_GT(sqrt(power) * (1<<i), 0.93); // Verify that signal level approximately still corresponds to filtered bandwidth
        EXPECT_LE(sqrt(power) * (1<<i), 1); // Verify that signal level approximately still corresponds to filtered bandwidth

        TF_DATA_TYPE side_lobe_level = 0.1 * level;
        for (int side_lobe_inx = 1; side_lobe_inx <= 4; side_lobe_inx++)
        {
            EXPECT_GE(ptr[2 * side_lobe_inx - 1], 0);
            EXPECT_LE(ptr[2 * side_lobe_inx    ], 0);

            EXPECT_LT(ptr[2 * side_lobe_inx - 1], side_lobe_level);
            EXPECT_GT(ptr[2 * side_lobe_inx    ], -side_lobe_level);

        }

        level /= 2;
    }
}

// Test that a sinewave gets filtered out
TEST(DyadicFilter, SineWave) {
    constexpr int n_octaves = 10;
    double PI = 4.0 * atan(1);
    constexpr double test_frequency = 1.0 / 20.0; // 0.05. OK in octave 0 (< 0.4), 1 (< 0.2), 2 (< 0.1), 3 (> 0.05). OK here means within 0.1 dB . Gone in Octave 4 and below. Gone here means level attenuated by 90dB
    DyadicFilter filter(n_octaves);
    EXPECT_EQ(filter.findOctave(test_frequency), 3);
    int xtra = DyadicFilter::getExtraSamples(n_octaves - 1);
    int n_samples = 1 / test_frequency * 100; // Integer number of periods to make easy calculations
    vector<TF_DATA_TYPE> vInput(n_samples + 2 * xtra, 0);

    // Fill vector with sine wave
    for (int i = 0; i < vInput.size(); i++) {
        vInput[i] = sin(2 * PI * test_frequency * i);
    }
    filter.doAllocation(n_samples, xtra, xtra);
    filter.filterSamples(&vInput[xtra], n_samples, xtra, xtra);

    int nb = n_samples;
    for (int i = 0; i < n_octaves; i++)
    {
        TF_DATA_TYPE * ptr;
        int _nb;
        tie(ptr, _nb) = filter.getSamples(i, 0);
        ASSERT_GE(_nb, nb); // Logically wrong test if not fulfilled

        TF_DATA_TYPE power = *ptr * *ptr;
        for (int i = 1; i < nb; i++)
        {
            power += ptr[i] * ptr[i];
        }

        double db = 10 * log(power / nb * 2);
        if (i <= filter.findOctave(test_frequency)) {
            EXPECT_LE(abs(db), 0.1);
        } else {
            EXPECT_LT(db, -90);
        }
        nb /= 2;
    }
}
