#ifndef POLYGONREGION_H
#define POLYGONREGION_H

#include <vector>
#include <utility>
#include <limits>
#include "Region.h"


namespace TFT {
    class PolygonRegion : public IRegion {
    public:
        PolygonRegion();
        PolygonRegion(const std::vector<std::pair<float, float>> & region);

        /**
         * @brief setCorners
         * @param region
         * @return Validity of region
         */
        void setCorners(const std::vector<std::pair<float, float>> & region) override;

        /**
         * @brief isWithin
         * @param abscissa
         * @param ordinate
         * @return
         */
        bool isWithin(float abscissa, float ordinate) const override;

    private:
        /**
         * @brief abscissassOf
         * @param ordinate
         * @return An ordered list of crossings of "ordinate". Crossing also if region just touches. If region segment is parallel to abscissa axis, no crossing exists
         */
        std::vector<float> abscissasOf(float ordinate) const;

        std::vector<std::pair<float, float>> corners;

        mutable struct {
            float key;
            std::vector<float> values;
            void invalidate()
            {
                key = std::numeric_limits<float>::signaling_NaN();
                values.clear();
            }
            std::vector<float> set(float k, const std::vector<float> vals)
            {
                key = k;
                values = vals;
                return values;
            }
            bool lookup(float k, std::vector<float> & rval) const
            {
                if (k != key) return false;
                rval = values;
                return true;
            }
        } crossingsCache;
    };
}

#endif // POLYGONREGION_H
