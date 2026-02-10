#include <tft/PolygonRegion.h>
#include <iterator>
#include <algorithm>

TFT::PolygonRegionHorizontal::PolygonRegionHorizontal() {}
TFT::PolygonRegionHorizontal::PolygonRegionHorizontal(const std::vector<std::pair<float, float>> & region) {
    setCorners(region);
}

/**
 * @brief setCorners
 * @param region
 * @return Validity of region
 */
void TFT::PolygonRegionHorizontal::setCorners(const std::vector<std::pair<float, float>> & region) {
    corners = region;

    // add first point in the end in order to assure we close the polygon
    if (!corners.empty()) {
        corners.push_back(corners[0]);
    }
    crossingsCache.invalidate();

    // Remove  double points since this is confusing
    auto iter = corners.begin();
    std::pair<float, float> p;
    while(iter != corners.end()) {
        p = *iter;
        iter++;
        while (iter != corners.end() && p == *iter) {
            // remove duplicate
            iter = corners.erase(iter);
        }
    }
}

/**
 * @brief abscissassOf
 * @param ordinate
 * @return An ordered list of crossings of "ordinate". Crossing also if region just touches. If region segment is parallel to abscissa axis, no crossing exists
 */
std::vector<float> TFT::PolygonRegionHorizontal::abscissasOf(float ordinate) const {
    std::vector<float> crossings;
    if (crossingsCache.lookup(ordinate, crossings)) {
        return  crossings;
    }

    // Iterate all line segments and add crossings with ordinate
    // In this process, p2 is trailing p1
    // We end iteration when p1 reaches the end
    for (auto p1 = corners.cbegin(), p2 = corners.cbegin(); ++p1 != corners.cend(); p2++) {
        bool hasCross(false);
        hasCross = p1->second != p2->second && (ordinate >= p1->second && ordinate < p2->second || ordinate > p2->second && ordinate <= p1->second);
        if (hasCross) {
            if (p1->first == p2->first) {
                crossings.push_back(p1->first);
            } else {
                float a = (p2->second - p1->second) / (p2->first - p1->first);
                float b = p1->second - a * p1->first;
                crossings.push_back((ordinate - b) / a);
            }
        }
    }

    std::sort(crossings.begin(), crossings.end());
    return crossingsCache.set(ordinate, crossings);
}

bool TFT::PolygonRegionHorizontal::isWithin(float abscissa, float ordinate) const {
    // An empty region is the same as everything within
    if (corners.empty()) {
        return true;
    }
    std::vector<float> crossings = abscissasOf(ordinate);

    auto lower =
        std::lower_bound(crossings.cbegin(), crossings.cend(), abscissa);
    return (lower == crossings.cend())
               ? false
               : (bool)(0x1 &
                         std::distance(
                             crossings.cbegin(),
                             lower)); // Uneven number of crossings means within
}

