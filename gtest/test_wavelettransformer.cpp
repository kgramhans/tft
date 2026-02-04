#include "tft/ConfinedGaussianWaveletVoice.h"
#include <gtest/gtest.h>
#include <iomanip>
#include <tft/tft.h>
#include <random>
#include <thread>

using namespace TFT;
using namespace std;

// Test WaveletVoice with a simple forward/backward transform on a time/frequency localized signal
TEST(WaveletTransformer, ForwardBackward) {
    int octaves = 10;
    double fmax = 0.4;
    float Q = 10.0;
    float overlap = 75.0;
    float nSamples = 1024;
    vector<TF_DATA_TYPE> vKernel(nSamples, 0);

    // Use a confined Gaussian as our test signal (limited in time and frequency)
    DyadicFilter dummyFilter(10);
    ConfinedGaussianWaveletVoice w(0.1, 0, 0.5, Q, 0, &dummyFilter);
    auto vw = w.getWavelet();
    int halfSize = vw.size() >> 1;
    for (int i = 0; i < vw.size(); i++) {
        vKernel[nSamples / 2 - halfSize + i] = vw[i];
    }

    ITimeFrequencyTransformer * tft = new WaveletTransformer(octaves, fmax, Q, overlap);
    int nb = tft->forwardTransform(&vKernel[0], vKernel.size(), 0, 0);

    EXPECT_GT(nb, 0);

    vector<TF_DATA_TYPE> vSignal(vKernel.size());
    nb = tft->backwardTransform(vSignal);
    EXPECT_EQ(nb, vKernel.size());

    // Let's have a look at the diff
    // Calculate diff
    double sumDiffSqr = 0;
    double sumOriginal = 0;

    for (int i = 0; i < vKernel.size(); i++) {
        sumOriginal += vKernel[i] * vKernel[i];
        sumDiffSqr += (vKernel[i] - vSignal[i]) * (vKernel[i] - vSignal[i]);
    }
    double SNR = 10 * log(sumDiffSqr/sumOriginal);
    EXPECT_LT(SNR, -60);

    delete tft;
}

// Test that we get equal results for sequenced and non-sequenced transforms
TEST(WaveletTransformer, Sequenced) {
    int octaves = 10;
    double fmax = 0.4;
    float Q = 10.0;
    float overlap = 75.0;
    float nSamples = 1024;
    vector<TF_DATA_TYPE> vKernel(nSamples, 0);
    queue<thread> threadPool;

    // Use a confined Gaussian as our test signal (limited in time and frequency)
    DyadicFilter dummyFilter(10);
    ConfinedGaussianWaveletVoice w(0.1, 0, 0.5, Q, 0, &dummyFilter);
    auto vw = w.getWavelet();
    int halfSize = vw.size() >> 1;
    for (int i = 0; i < vw.size(); i++) {
        vKernel[nSamples / 2 - halfSize + i] = vw[i];
    }

    ITimeFrequencyTransformer * tft1 = new WaveletTransformer(octaves, fmax, Q, overlap);
    int nb = tft1->forwardTransform(&vKernel[0], vKernel.size(), 0, 0);
    EXPECT_GT(nb, 0);

    WaveletTransformer * tft2 = new WaveletTransformer(octaves, fmax, Q, overlap);
    int sequences = tft2->prepareParallelForwardSequences(&vKernel[0], vKernel.size(), 0, 0);
    for (int i = sequences; i--;) {
        if (threadPool.size() == thread::hardware_concurrency()) {
            threadPool.front().join();
            threadPool.pop();
        }
        threadPool.push(thread(&WaveletTransformer::executeForwardSequence, tft2, i)); //tft2->executeForwardSequence(i);
    }

    // Join all remaining threads
    while (!threadPool.empty()) {
        threadPool.front().join();
        threadPool.pop();
    }

    vector<TF_DATA_TYPE> vSignal1(vKernel.size());
    vector<TF_DATA_TYPE> vSignal2(vKernel.size());
    vector<TF_DATA_TYPE> vSignal3(vKernel.size());
    int nb1 = tft1->backwardTransform(vSignal1);
    int nb2 = tft2->backwardTransform(vSignal2);
    sequences = tft2->prepareParallelBackwardSequences();

    // In below loop we do not fully exploit concurrency since it is not necessarily the front thread which finishes first
    // But at least we get some degree of parallell threads
    // Doing more would need some degree of sync in order to monitor threads
    for (int i = sequences; i--;) {
        if (threadPool.size() == thread::hardware_concurrency()) {
            threadPool.front().join();
            threadPool.pop();
        }
        threadPool.push(thread(&WaveletTransformer::executeBackwardSequence, tft2, i));
    }

    // Join all remaining threads
    while (!threadPool.empty()) {
        threadPool.front().join();
        threadPool.pop();
    }

    int nb3 = tft2->backwardTransform(vSignal3);
    EXPECT_EQ(nb1, nb2);
    EXPECT_EQ(nb1, nb3);

    for (int i = 0; i < nb1; i++) {
        EXPECT_EQ(vSignal1[i], vSignal2[i]);
        EXPECT_EQ(vSignal1[i], vSignal3[i]);
    }

    delete tft1;
    delete tft2;
}


TEST(WaveletTransformer, TimingSingleThreaded) {
    int octaves = 10;
    double fmax = 0.4;
    float Q = 10.0;
    float overlap = 75.0;
    float nSamples = 44100;
    vector<TF_DATA_TYPE> vKernel(nSamples, 0);

    // Use a delta function
    vKernel[nSamples/2] = 1;

    ITimeFrequencyTransformer * tft = new WaveletTransformer(octaves, fmax, Q, overlap);
    int nb = tft->forwardTransform(&vKernel[0], vKernel.size(), 0, 0);
    EXPECT_GT(nb, 0);

    vector<TF_DATA_TYPE> vSignal(vKernel.size());
    nb = tft->backwardTransform(vSignal);
    EXPECT_EQ(nb, vKernel.size());

    delete tft;
}

TEST(WaveletTransformer, TimingMultiThreaded) {
    int octaves = 10;
    double fmax = 0.4;
    float Q = 10.0;
    float overlap = 75.0;
    float nSamples = 44100;
    vector<TF_DATA_TYPE> vKernel(nSamples, 0);
    queue<thread> threadPool;

    // Use a delta function
    vKernel[nSamples/2] = 1;

    WaveletTransformer * tft = new WaveletTransformer(octaves, fmax, Q, overlap);
    int sequences = tft->prepareParallelForwardSequences(&vKernel[0], vKernel.size(), 0, 0);
    for (int i = sequences; i--;) {
        if (threadPool.size() == thread::hardware_concurrency()) {
            threadPool.front().join();
            threadPool.pop();
        }
        threadPool.push(thread(&WaveletTransformer::executeForwardSequence, tft, i)); //tft2->executeForwardSequence(i);
    }
    // Join all remaining threads
    while (!threadPool.empty()) {
        threadPool.front().join();
        threadPool.pop();
    }


    vector<TF_DATA_TYPE> vSignal(vKernel.size());
    sequences = tft->prepareParallelBackwardSequences();
    for (int i = sequences; i--;) {
        if (threadPool.size() == thread::hardware_concurrency()) {
            threadPool.front().join();
            threadPool.pop();
        }
        threadPool.push(thread(&WaveletTransformer::executeBackwardSequence, tft, i));
    }

    // Join all remaining threads
    while (!threadPool.empty()) {
        threadPool.front().join();
        threadPool.pop();
    }

    int nb = tft->backwardTransform(vSignal);
    EXPECT_EQ(nb, vKernel.size());

    delete tft;
}

// Test WaveletVoice with a simple forward/backward transform on a time/frequency localized signal
TEST(WaveletTransformer, Region) {
    int octaves = 10;
    double fmax = 0.41;
    float Q = 10.0;
    float overlap = 75.0;
    float nSamples = 1024;
    float kernelFreq = 0.1;
    vector<TF_DATA_TYPE> vKernel(nSamples, 0);

    // Use a confined Gaussian as our test signal (limited in time and frequency)
    DyadicFilter dummyFilter(10);
    ConfinedGaussianWaveletVoice w(kernelFreq, 0, 0.5, Q, 0, &dummyFilter);
    auto vw = w.getWavelet();
    int halfSize = vw.size() >> 1;
    for (int i = 0; i < vw.size(); i++) {
        vKernel[nSamples / 2 - halfSize + i] = vw[i];
    }

    ITimeFrequencyTransformer * tft = new WaveletTransformer(octaves, fmax, Q, overlap);
    int nb = tft->forwardTransform(&vKernel[0], vKernel.size(), 0, 0);

    EXPECT_GT(nb, 0);

    vector<TF_DATA_TYPE> vSignalRef(vKernel.size());
    vector<TF_DATA_TYPE> vSignal(vKernel.size());
    nb = tft->backwardTransform(vSignalRef);
    EXPECT_EQ(nb, vKernel.size());

    vector<pair<float, float>> region;

    // Apply an invalid region and let's see that it generates an empty output
    region.clear();
    region.push_back(make_pair(nSamples / 2.0, 0.0));
    region.push_back(make_pair(nSamples / 2.0, 0.5));
    tft->setPolygonRegion(region);
    nb = tft->backwardTransform(vSignal);
    EXPECT_EQ(nb, vKernel.size());
    for (int i = 0; i < vKernel.size(); i++) {
        EXPECT_EQ(vSignal[i], 0);
    }

    // Apply an empty region. Then we should get a zero signal
    region.clear();
    region.push_back(make_pair(nSamples + 1 , 0.0));
    region.push_back(make_pair(nSamples + 1, 0.25));
    region.push_back(make_pair(nSamples + 1, 0.5));
    tft->setPolygonRegion(region);
    nb = tft->backwardTransform(vSignal);
    EXPECT_EQ(nb, vKernel.size());
    for (int i = 0; i < vKernel.size(); i++) {
        EXPECT_EQ(vSignal[i], 0);
    }
    region.clear();
    region.push_back(make_pair(nSamples + 1 , 0.0));
    region.push_back(make_pair(nSamples + 2, 0.25));
    region.push_back(make_pair(nSamples + 3, 0.5));
    tft->setPolygonRegion(region);
    nb = tft->backwardTransform(vSignal);
    EXPECT_EQ(nb, vKernel.size());
    for (int i = 0; i < vKernel.size(); i++) {
        EXPECT_EQ(vSignal[i], 0);
    }

    // Apply a region which includes everything. Then we should get equal signals
    region.clear();
    region.push_back(make_pair(-1, 0.0));
    region.push_back(make_pair(-1, 0.5));
    region.push_back(make_pair(nSamples, 0.5));
    region.push_back(make_pair(nSamples, 0.0));
    tft->setPolygonRegion(region);
    nb = tft->backwardTransform(vSignal);
    EXPECT_EQ(nb, vKernel.size());
    for (int i = 0; i < vKernel.size(); i++) {
        EXPECT_EQ(vSignal[i], vSignalRef[i]);
    }

    // Apply half the region. Verify plausible energy
    region.clear();
    region.push_back(make_pair(-1, 0.0));
    region.push_back(make_pair(-1, 0.5));
    region.push_back(make_pair(nSamples / 2, 0.5));
    region.push_back(make_pair(nSamples / 2, 0.0));
    tft->setPolygonRegion(region);
    nb = tft->backwardTransform(vSignal);
    EXPECT_EQ(nb, vKernel.size());
    double powerRef = 0;
    double powerLeft = 0;
    for (int i = 0; i < vKernel.size(); i++) {
        powerLeft    += vSignal   [i] * vSignal   [i];
        powerRef += vSignalRef[i] * vSignalRef[i];
    }
    EXPECT_LT(powerLeft, 0.5 * powerRef);

    region.clear();
    region.push_back(make_pair(2 * nSamples, 0.0));
    region.push_back(make_pair(2 * nSamples, 0.5));
    region.push_back(make_pair(nSamples / 2, 0.5));
    region.push_back(make_pair(nSamples / 2, 0.0));
    tft->setPolygonRegion(region);
    nb = tft->backwardTransform(vSignal);
    EXPECT_EQ(nb, vKernel.size());
    double powerRight = 0;
    for (int i = 0; i < vKernel.size(); i++) {
        powerRight    += vSignal   [i] * vSignal   [i];
    }
    EXPECT_LT(powerRight, 0.5 * powerRef);
    EXPECT_LT(powerRight + powerLeft, powerRef); //(a+b)(a+b) > a*a + b*b
    EXPECT_GT(powerRight + powerLeft, 0.5 * powerRef); //(a+b)(a+b) > a*a + b*b

    region.clear();
    region.push_back(make_pair(-1, 0.0));
    region.push_back(make_pair(-1, kernelFreq));
    region.push_back(make_pair(2 * nSamples , kernelFreq));
    region.push_back(make_pair(2 * nSamples, 0.0));
    tft->setPolygonRegion(region);
    nb = tft->backwardTransform(vSignal);
    EXPECT_EQ(nb, vKernel.size());
    double powerLF = 0;
    for (int i = 0; i < vKernel.size(); i++) {
        powerLF    += vSignal   [i] * vSignal   [i];
    }
    EXPECT_LT(powerLF, 0.5 * powerRef);

    region.clear();
    region.push_back(make_pair(-1, 1.0));
    region.push_back(make_pair(-1, kernelFreq));
    region.push_back(make_pair(2 * nSamples , kernelFreq));
    region.push_back(make_pair(2 * nSamples, 1.0));
    tft->setPolygonRegion(region);
    nb = tft->backwardTransform(vSignal);
    EXPECT_EQ(nb, vKernel.size());
    double powerHF = 0;
    for (int i = 0; i < vKernel.size(); i++) {
        powerHF    += vSignal   [i] * vSignal   [i];
    }
    EXPECT_LT(powerHF, 0.5 * powerRef);
    EXPECT_LT(powerLF + powerHF, powerRef); //(a+b)(a+b) > a*a + b*b
    EXPECT_GT(powerLF + powerHF, 0.5 * powerRef); //(a+b)(a+b) > a*a + b*b


    // Apply an invalid region again and let's see that we again have zero output
    region.clear();
    region.push_back(make_pair(nSamples / 2.0, 0.0));
    region.push_back(make_pair(nSamples / 2.0, 0.5));
    tft->setPolygonRegion(region);
    nb = tft->backwardTransform(vSignal);
    EXPECT_EQ(nb, vKernel.size());
    for (int i = 0; i < vKernel.size(); i++) {
        EXPECT_EQ(vSignal[i], 0);
    }

    delete tft;
}

// Test normalisation of wavelet transform under verying analysis methods and stimulation
TEST(WaveletTransformer, Normalisation) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist100(0,100);

    int octaves = 10;
    double fmax = 0.4;
    float nSamples = 1024;
    vector<TF_DATA_TYPE> vKernel(nSamples, 0);

    int repeat = 10;
    while (repeat--) {
        float frequency = 0.125 + 0.00025 * dist100(rng); // Vary in range 0.1 .. 0.15
        float Q = 10.0 + 0.2 * dist100(rng); // Vary in range 10..30
        float overlap = 90.0 - 0.4 * dist100(rng); // Vary in range 50..90
        float nSamples = 1024;
        float fDelay = -1 + 0.02 * dist100(rng); // Vary in range -1 .. 1

        // Use a confined Gaussian as our test signal (limited in time and frequency)
        DyadicFilter dummyFilter(10);
        ConfinedGaussianWaveletVoice w(frequency, fDelay, 0, 0.5, 10, 0, &dummyFilter);
        auto vw = w.getWavelet();
        int halfSize = vw.size() >> 1;
        for (int i = 0; i < vw.size(); i++) {
            vKernel[nSamples / 2 - halfSize + i] = vw[i];
        }

        ITimeFrequencyTransformer * tft = new WaveletTransformer(octaves, fmax, Q, overlap);
        int nb = tft->forwardTransform(&vKernel[0], vKernel.size(), 0, 0);
        vector<TF_DATA_TYPE> vSignal(vKernel.size());
        tft->backwardTransform(vSignal);

        // Generate the two sums and diff power
        double sumOriginal = 0;
        double sumRestored = 0;
        double sumDiffSqr = 0;

        for (int i = 0; i < vKernel.size(); i++) {
            sumOriginal += vKernel[i] * vKernel[i];
            sumRestored += vSignal[i] * vSignal[i];
        }

        double ratio = sqrt(sumOriginal / sumRestored);

        // Calculate diff
        for (int i = 0; i < vKernel.size(); i++) {
            sumDiffSqr += (vKernel[i] - vSignal[i]) * (vKernel[i] - vSignal[i]);
        }
        double SNR = 10 * log(sumDiffSqr/sumOriginal);

        EXPECT_LT(SNR, -32)  << std::fixed << std::setw(10) << nb << "\t" << frequency << "\t" << Q << "\t" << overlap << "\t" << fDelay << "\t" << ratio << "\t" << (int)SNR << "\t" << sqrt(sumOriginal) << "\t" << sqrt(sumRestored)<< "\t";
        delete tft;
    }
}

// Test restore with a random noise signal
// test signal is upsampled in order to avoid aliasing effects close to Fs/2
TEST(WaveletTransformer, Noise) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist100(0,100);
    std::uniform_real_distribution<> destUniform(-1,1);

    int octaves = 10;
    double fmax = 0.4;
    float nSamples = 1024;
    vector<TF_DATA_TYPE> vKernel(nSamples + DyadicFilter::getUpsamplingPaddingSize() * 2, 0);

    int repeat = 10;
    while (repeat--) {
        float Q = 10.0 + 0.2 * dist100(rng); // Vary in range 10..30
        float overlap = 87.5;
        float nSamples = 1000;


        // Use a 10% cosine tapered random signal
        for (int i = 0; i < nSamples; i++) {
            TF_DATA_TYPE factor = 1.0;
            if (i < nSamples / 10) {
                factor = 0.5 * (1 - cos(4 * atan(1) * i * 10.0 / nSamples));
            }
            if ((nSamples - i - 1) < nSamples / 10) {
                factor = 0.5 * (1 - cos(4 * atan(1) * (nSamples - i - 1) * 10.0 / nSamples));
            }
            vKernel[i + DyadicFilter::getUpsamplingPaddingSize()] = factor * destUniform(rng);
        }
        vKernel = DyadicFilter::upsample(vKernel);

        ITimeFrequencyTransformer * tft = new WaveletTransformer(octaves, fmax, Q, overlap);
        int nb = tft->forwardTransform(&vKernel[DyadicFilter::getUpsamplingPaddingSize()], 2 * nSamples, DyadicFilter::getUpsamplingPaddingSize(), DyadicFilter::getUpsamplingPaddingSize());
        vector<TF_DATA_TYPE> vSignal(nSamples * 2);
        tft->backwardTransform(vSignal);

        // Generate the two sums and diff power
        double sumOriginal = 0;
        double sumDiffSqr = 0;

        for (int i = 0; i < nSamples * 2; i++) {
            sumOriginal += vKernel[i + DyadicFilter::getUpsamplingPaddingSize()] * vKernel[i + DyadicFilter::getUpsamplingPaddingSize()];
        }

        // Calculate diff
        for (int i = 0; i < nSamples * 2; i++) {
            sumDiffSqr += (vKernel[i + DyadicFilter::getUpsamplingPaddingSize()] - vSignal[i]) * (vKernel[i + DyadicFilter::getUpsamplingPaddingSize()] - vSignal[i]);
        }
        double SNR = 10 * log(sumDiffSqr/sumOriginal);

        EXPECT_LT(SNR, -40);
        delete tft;
    }
}

