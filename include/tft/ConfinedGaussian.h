#ifndef CONFINEDGAUSSIAN_H
#define CONFINEDGAUSSIAN_H

#include <cassert>
#include <vector>
#include <tft/version.h>

namespace TFT {

/***
 * DEFINITIONS:
 *
 * Compute a Confined Gaussian Window function with continuous length L:
 * The continuous Confined Guassian ( cg(x) ) is defined by:
 * cg(x) = 0 for x < -0.5
 * cg(x) = 0 for x > L - 0.5
 * cg(x) = max for x = (L - 1)/2
 *
 * Thus, If we want max to occur at an integer number, it follows that L must be an odd integer and that
 *
 * cg(n) != 0 for n=0,1,2, ... ,L-1, totalling to L non-zero values
 *
 * Note that it is thus not possible to center the window on a sample without making it odd
 * Specifically, this must be taken into account when applying it to STFT which uses integer (power-of-two) window lengths
 * *****/

class ConfinedGaussian
{
protected:
    ConfinedGaussian(double _sigma);
    ConfinedGaussian(double _sigma, unsigned int _length);
    void setLength(unsigned int _length);
    float approximateConfinedGaussian(float x);  // x ranging from 0 to Length
    std::vector<TF_DATA_TYPE> getWindow(unsigned int N, unsigned int & center);
private:
    float gaussian(float x);
    double sigma;
    int Center;     ///< TOffset of center
    int Length;      ///< Note that continuous Confined Gaussian window extends half a sample on either side of halfLen: Length = 2 * (HalfLen + 0.5), which also equals the number of non-zero samples - which is odd
};

} // namespace TFT

#endif // CONFINEDGAUSSIAN_H
