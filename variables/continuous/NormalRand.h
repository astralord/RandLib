#ifndef NORMALRAND_H
#define NORMALRAND_H

#include "ContinuousRand.h"
#include "UniformRand.h"
#include "../BasicRandGenerator.h"

/**
 * @brief The NormalRand class
 */
class RANDLIBSHARED_EXPORT NormalRand : public ContinuousRand
{
    double mu, sigma;
    double sigmaSqrt2Inv; /// 1 / (sigma * sqrt(2))

    //TODO: find a way to initialize them without dummy
    /// Tables for ziggurat
    static unsigned long kn[128];
    static double wn[128], fn[128];
    static const bool dummy;
    /// Set up ziggurat tables
    static bool setupTables();

    UniformRand U;
    BasicRandGenerator B;

    /**
     * @brief ziggurat
     * @return standard normal variable: X ~ N(0,1)
     */
    double ziggurat();

public:
    NormalRand(double mean = 0, double var = 1);

    void setMean(double mean);
    void setSigma(double rootVar);
    void setVar(double var) { setSigma(std::sqrt(std::max(var, MIN_POSITIVE))); }
    double getSigma() const { return sigma; }

    virtual double f(double x) const;
    virtual double F(double x) const;
    virtual double value();

    double E() const { return mu; }
    double Var() const { return sigma * sigma; }

    inline double Median() const { return mu; }
    inline double Mode() const { return mu; }
    static constexpr double Skewness() { return 0; }
    static constexpr double ExcessKurtosis() { return 0; }

    bool fitToData(const QVector<double> &sample);
};

#endif // NORMALRAND_H
