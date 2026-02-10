#ifndef POLYGONREGIONVERTICAL_H
#define POLYGONREGIONVERTICAL_H

#include "PolygonRegion.h"

namespace TFT {
class PolygonRegionVertical : public PolygonRegionHorizontal
{
public:
    PolygonRegionVertical();
    PolygonRegionVertical(const std::vector<std::pair<float, float>> & region);

    /**
         * @brief isWithin
         * @param abscissa
         * @param ordinate
         * @return
         */
    bool isWithin(float abscissa, float ordinate) const override;

    /**
         * @brief clone
         * @return pointer to yet another instance
         */
    std::unique_ptr<IRegion> cloneRegion() const override {return std::make_unique<PolygonRegionVertical>(*this);}

    ~PolygonRegionVertical() {}

private:
    /**
         * @brief ordinatesOf
         * @param abscissa
         * @return An ordered list of crossings of "abscissa". Crossing also if region just touches. If region segment is parallel to ordinate axis, no crossing exists
         */
    std::vector<float> ordinatesOf(float abscissa) const;
};
}
#endif // POLYGONREGIONVERTICAL_H
