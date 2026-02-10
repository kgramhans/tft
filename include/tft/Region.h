#ifndef REGION_H
#define REGION_H

#include <vector>
#include <utility>
#include <memory>


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

    /**
     * @brief clone
     * @return pointer to yet another instance
     */
    virtual std::unique_ptr<IRegion> cloneRegion() const = 0;

    virtual ~IRegion() {}
};
}

#endif // REGION_H
