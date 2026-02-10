#include "tft/PolygonRegionVertical.h"
#include <gtest/gtest.h>
#include <tft/tft.h>
#include <tft/StftMoment.h>
#include <random>

using namespace TFT;
using namespace std;


// Test WaveletVoice
TEST(StftMoment, GeneralConfiguration) {
    // Verify general construction safeguards against illegal parameters


    // Illegal overlap
    EXPECT_DEBUG_DEATH(
    {
        double when = 20;
        unsigned int windowLen = 1024;
        auto pWindow = std::make_shared<std::vector<TF_DATA_TYPE>>(windowLen, 1);   // Simple square window
        unsigned int windowCenter = windowLen / 2;
        float overlap = 100; // Illegal
        auto h = GetFFT(windowLen);

        StftMoment moment(when, windowCenter, overlap, h, pWindow);
    }, "Assertion");

    // Illegal windowLen
    EXPECT_DEBUG_DEATH(
        {
            double when = 20;
            unsigned int windowLen = 1023;
            auto pWindow = std::make_shared<std::vector<TF_DATA_TYPE>>(windowLen, 1);   // Simple square window
            unsigned int windowCenter = windowLen / 2;
            float overlap = 90;
            auto h = GetFFT(windowLen);

            StftMoment moment(when, windowCenter, overlap, h, pWindow);
        }
        , "Assertion");

    // Illegal FFT handle
    EXPECT_DEBUG_DEATH(
        {
            double when = 20;
            unsigned int windowLen = 1024;
            auto pWindow = std::make_shared<std::vector<TF_DATA_TYPE>>(windowLen, 1);   // Simple square window
            unsigned int windowCenter = windowLen / 2;
            float overlap = 90;
            auto h = GetFFT(windowLen * 2);

            StftMoment moment(when, windowCenter, overlap, h, pWindow);
        }
        , "Assertion");

    // Illegal window center
    EXPECT_DEBUG_DEATH(
        {
            double when = 20;
            unsigned int windowLen = 1024;
            auto pWindow = std::make_shared<std::vector<TF_DATA_TYPE>>(windowLen, 1);   // Simple square window
            unsigned int windowCenter = windowLen * 2;
            float overlap = 90;
            auto h = GetFFT(windowLen);

            StftMoment moment(when, windowCenter, overlap, h, pWindow);
        }
        , "Assertion");

    // Everything OK
        {
            double when = 20;
            unsigned int windowLen = 1024;
            auto pWindow = std::make_shared<std::vector<TF_DATA_TYPE>>(windowLen, 1);   // Simple square window
            unsigned int windowCenter = windowLen / 2;
            float overlap = 90;
            auto h = GetFFT(windowLen);

            StftMoment moment(when, windowCenter, overlap, h, pWindow);
        }

}

TEST(StftMoment, TransformNoOverlap) {
    double when = 1024;
    unsigned int windowLen = 1024;
    auto pWindow = std::make_shared<std::vector<TF_DATA_TYPE>>(windowLen, 2.0 / windowLen);   // Simple square window. Note that average must be two due to one-sided spectrum
    unsigned int windowCenter = windowLen / 2;
    float overlap = 0;
    auto h = GetFFT(windowLen);

    StftMoment moment(when, windowCenter, overlap, h, pWindow);

    // Do a transform of a simple signal
    std::vector<TF_DATA_TYPE> vSignal(2 * windowLen, 0);
    vSignal[windowLen] = 1;

    moment.transform(&vSignal[0]);

    // Transform back and verify result
    PolygonRegionVertical region;
    auto signal = moment.constructMomentSignal(region.cloneRegion());

    // Verify equality
    EXPECT_EQ(signal.size(), windowLen);
    for (int i = 0; i < windowLen; i++) {
        EXPECT_EQ(signal[i], vSignal[i + moment.getConstructedMomentOffset()]) << i;
    }
}

TEST(StftMoment, TransformOverlap90) {
    unsigned int windowLen = 1024;
    auto pWindow = std::make_shared<std::vector<TF_DATA_TYPE>>(windowLen, 2.0 / (windowLen - 1));   // Simple square window. Note that average must be two due to one-sided spectrum
    pWindow->back() = 0;
    unsigned int windowCenter = windowLen / 2 - 1;
    float overlap = 90;
    auto h = GetFFT(windowLen);
    unsigned int step = (unsigned int)((100-overlap)/100 * 1024.0 + 0.5);
    unsigned int transforms = windowLen / step;
    float actual_overlap = (1.0 - step/1024.0) * 100.0;

    // Do a transform of a simple signal
    unsigned int deltaPosition = 1024;
    std::vector<TF_DATA_TYPE> vSignal(2 * windowLen, 0);
    vSignal[deltaPosition] = 1;

    // Place first transform such that it will just include the delta
    std::vector<double> vWhen(transforms, 0);
    for (int i = 0; i < transforms; i++) {
        vWhen[i] = deltaPosition - windowCenter + i * step;
    }

    // Create an array of moments
    std::vector<StftMoment> vMoment;
    for (int i = 0; i < transforms; i++) {
        vMoment.push_back(StftMoment(vWhen[i], windowCenter, actual_overlap, h, pWindow));
    }

    EXPECT_EQ(pWindow.use_count(), transforms + 1);

    // transform moments
    for (int i = 0; i < transforms; i++) {
        vMoment[i].transform(&vSignal[0]);
    }

    // Transform back and verify result
    PolygonRegionVertical region;
    std::vector<TF_DATA_TYPE> constructedSignal(2 * windowLen, 0);
    for (int i = 0; i < transforms; i++) {
        auto signal = vMoment[i].constructMomentSignal(region.cloneRegion());
        for (int inx = 0; inx < signal.size(); inx ++) {
            int ii = inx + vMoment[i].getConstructedMomentOffset();
            EXPECT_GE(ii, 0); // Sanity check
            EXPECT_LT(ii, 2 * windowLen); // Sanity check
            constructedSignal[ii] += signal[inx];
        }
    }

    // Verify almost-equality: Extract the two and verify the noise
    double err_sqr_sum = 0;
    for (int i = 0; i < 2 * windowLen; i++) {
        float delta = constructedSignal[i] - vSignal[i];
        err_sqr_sum += delta * delta;
    }
    float rms = sqrt(err_sqr_sum / (2 * windowLen));
    float SNR = - 20 * log(rms);
    EXPECT_GT(SNR, 100);

    // Test some cleaning up. Do not see why it should not work? But anyhow test it
    vMoment.clear();
    EXPECT_EQ(pWindow.use_count(), 1);
}
