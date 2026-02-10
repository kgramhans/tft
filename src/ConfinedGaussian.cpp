#include <cassert>
#include <tft/ConfinedGaussian.h>
#include <cmath>

namespace TFT {

ConfinedGaussian::ConfinedGaussian(double _sigma, unsigned int _length)
    : sigma(_sigma) {
    setLength(_length);
}

ConfinedGaussian::ConfinedGaussian(double _sigma)
    : sigma(_sigma), Length(0), Center(0) {
}

void ConfinedGaussian::setLength(unsigned int _length) {
    assert(_length & 1);
    Length = _length;
    Center = (Length - 1) / 2;
}

float TFT::ConfinedGaussian::gaussian(float x)
{
    float argument = (x - Center) / 2 / Length / sigma;
    return exp(- argument * argument);
}

float TFT::ConfinedGaussian::approximateConfinedGaussian(float x)
{
    assert(Length > 0);
    if (x < 0.5 || x > Length - 0.5) {
        return 0.0;
    }
    return gaussian(x) - gaussian(-0.5) * (gaussian(x + Length) + gaussian(x - Length)) / (gaussian(-0.5 + Length) + gaussian(-0.5 - Length));
}

std::vector<TF_DATA_TYPE> TFT::ConfinedGaussian::getWindow(unsigned int N, unsigned int & center) {
    std::vector<TF_DATA_TYPE> rval(N, 0);
    for (int i = 0; i < N; i++) {
        rval[i] = approximateConfinedGaussian(i);
    }
    center = Center;
    return rval;
}

} // namespace TFT
