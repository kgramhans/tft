#ifndef STFTMOMENT_H
#define STFTMOMENT_H

#include "tft/RealFFTf.h"
#include "tft/Region.h"
#include "tft/version.h"
#include <cassert>
#include <vector>
#include <memory>

namespace TFT {

class StftMoment
{
public:
    StftMoment(unsigned int momentInTime,
               unsigned int windowCenterOffset,
               double overlap,
               const std::shared_ptr<FFTParam> h,
               std::shared_ptr<std::vector<TF_DATA_TYPE>> _window);
    virtual ~StftMoment();
    virtual void transform(const TF_DATA_TYPE * pSamples);
    std::vector<TF_DATA_TYPE> constructMomentSignal(const std::unique_ptr<IRegion> & region) const;
    void constructMomentSignalBuffer(const std::unique_ptr<IRegion> &region) const;
    int getConstructedMomentOffset() const {return moment - windowOffset;}

protected:
    double getWindowEnergy() const {assert(windowEnergy);return windowEnergy;}

    std::shared_ptr<std::vector<TF_DATA_TYPE>> window;
    double windowEnergy;
    unsigned int windowLength;
    float duration;
    float bandwidth;
    unsigned int moment;
    unsigned int windowOffset;
    mutable struct
    {
        const void * key;
        TF_DATA_TYPE value;
        void invalidate()
        {
            key = NULL;
        }
        TF_DATA_TYPE set(const void * k, TF_DATA_TYPE val)
        {
            key = k;
            value = val;
            return value;
        }
        bool lookup(const void * k, TF_DATA_TYPE & rval) const
        {
            if (k != key) return false;
            rval = value;
            return true;
        }
    } valueCache;
    std::vector<TF_DATA_TYPE> resultZ; ///< This member points to interleaved (Re,Im) data
    mutable std::vector<TF_DATA_TYPE> momentSignalBuffer;    ///< Used for storage of interim results during parallel execution
    const std::shared_ptr<FFTParam> handle;
    static constexpr float pi = 3.14159265359;
};

} // namespace TFT

#endif // STFTMOMENT_H
