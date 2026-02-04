#include <gtest/gtest.h>
#include <tft/PolygonRegion.h>

using namespace TFT;
using namespace std;

// Test WaveletVoice with a simple forward/backward transform on a time/frequency localized signal
TEST(Region, Polygon) {
    PolygonRegion region;

    // Test no region is always within
    EXPECT_TRUE(region.isWithin(0,0));
    EXPECT_TRUE(region.isWithin(1,1));
    EXPECT_TRUE(region.isWithin(0.5,0.5));
    EXPECT_TRUE(region.isWithin(0.4,0.5));
    EXPECT_TRUE(region.isWithin(0.6,0.5));
    EXPECT_TRUE(region.isWithin(0.5,0.4));
    EXPECT_TRUE(region.isWithin(0.5,0.6));
    EXPECT_TRUE(region.isWithin(0,1));
    EXPECT_TRUE(region.isWithin(1,0));
    EXPECT_TRUE(region.isWithin(0.5,0.5));

    // Test invalid region is never within
    region.setCorners ({make_pair(0,0), make_pair(1,1)});
    EXPECT_FALSE(region.isWithin(0,0));
    EXPECT_FALSE(region.isWithin(1,1));
    EXPECT_FALSE(region.isWithin(0.5,0.5));
    EXPECT_FALSE(region.isWithin(0.4,0.5));
    EXPECT_FALSE(region.isWithin(0.6,0.5));
    EXPECT_FALSE(region.isWithin(0.5,0.4));
    EXPECT_FALSE(region.isWithin(0.5,0.6));
    EXPECT_FALSE(region.isWithin(0,1));
    EXPECT_FALSE(region.isWithin(1,0));

    // Test for a square
    region.setCorners ({make_pair(0.1,0.1), make_pair(0.1,0.9), make_pair(0.9,0.9), make_pair(0.9,0.1)});
    EXPECT_FALSE(region.isWithin(0,0));
    EXPECT_FALSE(region.isWithin(1,1));
    EXPECT_TRUE(region.isWithin(0.5,0.5));
    EXPECT_TRUE(region.isWithin(0.4,0.5));
    EXPECT_TRUE(region.isWithin(0.6,0.5));
    EXPECT_TRUE(region.isWithin(0.5,0.4));
    EXPECT_TRUE(region.isWithin(0.5,0.6));
    EXPECT_FALSE(region.isWithin(0,1));
    EXPECT_FALSE(region.isWithin(1,0));

    // Test for a more detailed square
    region.setCorners ({make_pair(0.1,0.1), make_pair(0.1,0.5), make_pair(0.1,0.9), make_pair(0.5,0.9), make_pair(0.9,0.9), make_pair(0.9,0.5), make_pair(0.9,0.1), make_pair(0.5,0.1)});
    EXPECT_FALSE(region.isWithin(0,0));
    EXPECT_FALSE(region.isWithin(1,1));
    EXPECT_TRUE(region.isWithin(0.5,0.5));
    EXPECT_TRUE(region.isWithin(0.4,0.5));
    EXPECT_TRUE(region.isWithin(0.6,0.5));
    EXPECT_TRUE(region.isWithin(0.5,0.4));
    EXPECT_TRUE(region.isWithin(0.5,0.6));
    EXPECT_FALSE(region.isWithin(0,1));
    EXPECT_FALSE(region.isWithin(1,0));

    // Test for a weird square
    region.setCorners ({make_pair(0.1,0.1), make_pair(0.1,0.9), make_pair(0.9,0.9), make_pair(0.9,0.9), make_pair(0.9,0.9), make_pair(0.9,0.1), make_pair(0.1,0.1)});
    EXPECT_FALSE(region.isWithin(0,0));
    EXPECT_FALSE(region.isWithin(1,1));
    EXPECT_TRUE(region.isWithin(0.5,0.5));
    EXPECT_TRUE(region.isWithin(0.4,0.5));
    EXPECT_TRUE(region.isWithin(0.6,0.5));
    EXPECT_TRUE(region.isWithin(0.5,0.4));
    EXPECT_TRUE(region.isWithin(0.5,0.6));
    EXPECT_FALSE(region.isWithin(0,1));
    EXPECT_FALSE(region.isWithin(1,0));

    // Test for a butterfly
    region.setCorners ({make_pair(0.1,0.1), make_pair(0.9,0.9), make_pair(0.1,0.9), make_pair(0.9,0.1)});
    EXPECT_FALSE(region.isWithin(0,0));
    EXPECT_FALSE(region.isWithin(1,1));
    EXPECT_FALSE(region.isWithin(0.4,0.5));
    EXPECT_FALSE(region.isWithin(0.6,0.5));
    EXPECT_TRUE(region.isWithin(0.5,0.4));
    EXPECT_TRUE(region.isWithin(0.5,0.6));
    EXPECT_FALSE(region.isWithin(0,1));
    EXPECT_FALSE(region.isWithin(1,0));

    // Test for a butterfly the other way round
    region.setCorners ({make_pair(0.1,0.1), make_pair(0.9,0.9), make_pair(0.9,0.1), make_pair(0.1,0.9)});
    EXPECT_FALSE(region.isWithin(0,0));
    EXPECT_FALSE(region.isWithin(1,1));
    EXPECT_TRUE(region.isWithin(0.4,0.5));
    EXPECT_TRUE(region.isWithin(0.6,0.5));
    EXPECT_FALSE(region.isWithin(0.5,0.4));
    EXPECT_FALSE(region.isWithin(0.5,0.6));
    EXPECT_FALSE(region.isWithin(0,1));
    EXPECT_FALSE(region.isWithin(1,0));

    // Test for a double butterfly the other way round. It does make sense that we will never be within
    region.setCorners ({make_pair(0.1,0.1), make_pair(0.9,0.9), make_pair(0.9,0.1), make_pair(0.1,0.9), make_pair(0.1,0.1), make_pair(0.9,0.9), make_pair(0.9,0.1), make_pair(0.1,0.9)});
    EXPECT_FALSE(region.isWithin(0,0));
    EXPECT_FALSE(region.isWithin(1,1));
    EXPECT_FALSE(region.isWithin(0.4,0.5));
    EXPECT_FALSE(region.isWithin(0.6,0.5));
    EXPECT_FALSE(region.isWithin(0.5,0.4));
    EXPECT_FALSE(region.isWithin(0.5,0.6));
    EXPECT_FALSE(region.isWithin(0,1));
    EXPECT_FALSE(region.isWithin(1,0));

    // Test for another square
    region.setCorners ({make_pair(0.0,0.5), make_pair(0.5,1.0), make_pair(1.0,0.5), make_pair(0.5,0.0)});
    EXPECT_FALSE(region.isWithin(0,0));
    EXPECT_FALSE(region.isWithin(1,1));
    EXPECT_TRUE(region.isWithin(0.5,0.5));
    EXPECT_TRUE(region.isWithin(0.4,0.5));
    EXPECT_TRUE(region.isWithin(0.6,0.5));
    EXPECT_TRUE(region.isWithin(0.5,0.4));
    EXPECT_TRUE(region.isWithin(0.5,0.6));
    EXPECT_FALSE(region.isWithin(0,1));
    EXPECT_FALSE(region.isWithin(1,0));

}
