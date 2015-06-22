#ifndef NORMALRAND_H
#define NORMALRAND_H

#include <RandomVariable.h>
#include "UniformRand.h"

class RANDLIBSHARED_EXPORT NormalRand : public ContinuousRand
{
    double mu, sigma;
    double sigmaSqrt2Inv; /// 1 / (sigma * sqrt(2))

    //TODO: find a way to initialize them only once!!!
    /// Tables for ziggurat
    unsigned long kn[128];
    double wn[128], fn[128];
    UniformRand U;
    /**
     * @brief ziggurat
     * @return standard normal variable: X ~ N(0,1)
     */
    double ziggurat();

public:
    NormalRand(double mean, double sigma);

    void setMean(double mean);
    void setSigma(double rootVar);
    double getSigma() { return sigma; }

    virtual double pdf(double x);
    virtual double cdf(double x);
    virtual double value();

    inline double M() { return mu; }
    inline double Var() { return sigma * sigma; }

    inline double Median() { return mu; }
    inline double Mode() { return mu; }
    inline double Skewness() { return 0; }
    inline double ExcessKurtosis() { return 0; }
};


/*

double erfinv(double p)
{
    if (p < 0 || p > 1)
        return NAN;
    double t = M_SQRT2;
    t *= ((p < 0.5) ? sqrt(-log(p)) : sqrt(-log(1 - p)));
    double c[] = {2.515517, 0.802853, 0.010328};
    double d[] = {1.432788, 0.189269, 0.001308};
    double res = t - ((c[2] * t + c[1]) * t + c[0]) /
               (((d[2] * t + d[1]) * t + d[0]) * t + 1.0);
    return ((p < 0.5) ? -res : res);
}

double erfcinv(double p)
{
    if (p < 0 || p > 1)
        return NAN;
    return erfinv(1 - p);
}

 * */

#endif // NORMALRAND_H
