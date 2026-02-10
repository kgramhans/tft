#include <tft/PolygonRegionVertical.h>

namespace TFT {
PolygonRegionVertical::PolygonRegionVertical() {}

PolygonRegionVertical::PolygonRegionVertical(const std::vector<std::pair<float, float>> & region) : PolygonRegionHorizontal(region) {
}

/**
         * @brief isWithin
         * @param abscissa
         * @param ordinate
         * @return
         */
bool PolygonRegionVertical::isWithin(float abscissa, float ordinate) const {
    // An empty region is the same as everything within
    if (corners.empty()) {
        return true;
    }
    std::vector<float> crossings = ordinatesOf(abscissa);

    auto lower =
        std::lower_bound(crossings.cbegin(), crossings.cend(), ordinate);
    return (lower == crossings.cend())
               ? false
               : (bool)(0x1 &
                         std::distance(
                             crossings.cbegin(),
                             lower)); // Uneven number of crossings means within
}

/**
         * @brief ordinatesOf
         * @param abscissa
         * @return An ordered list of crossings of "abscissa". Crossing also if region just touches. If region segment is parallel to ordinate axis, no crossing exists
         */
std::vector<float> PolygonRegionVertical::ordinatesOf(float abscissa) const {
    std::vector<float> crossings;
    if (crossingsCache.lookup(abscissa, crossings)) {
        return  crossings;
    }

    // Iterate all line segments and add crossings with abscissa
    // In this process, p2 is trailing p1
    // We end iteration when p1 reaches the end
    for (auto p1 = corners.cbegin(), p2 = corners.cbegin(); ++p1 != corners.cend(); p2++) {
        bool hasCross(false);
        hasCross = p1->first != p2->first && (abscissa >= p1->first && abscissa < p2->first || abscissa > p2->first && abscissa <= p1->first);
        if (hasCross) {
            float a = (p2->second - p1->second) / (p2->first - p1->first);
            float b = p1->second - a * p1->first;
            crossings.push_back(a * abscissa + b);
        }
    }

    std::sort(crossings.begin(), crossings.end());
    return crossingsCache.set(abscissa, crossings);
}

}
