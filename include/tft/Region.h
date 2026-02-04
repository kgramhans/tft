#ifndef REGION_H
#define REGION_H

#include <vector>
#include <utility>


namespace TFT {
class IRegion {
public:
    /**
         * @brief setCorners
         * @param region
         * @return Validity of region
         */
    virtual void setCorners(const std::vector<std::pair<float, float>> & region) = 0;

    /**
         * @brief isWithin
         * @param abscissa
         * @param ordinate
         * @return
         */
    virtual bool isWithin(float abscissa, float ordinate) const = 0;

};
}

#endif // REGION_H
