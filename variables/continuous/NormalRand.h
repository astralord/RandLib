#ifndef NORMALRAND_H
#define NORMALRAND_H

#include "ContinuousRand.h"
#include "UniformRand.h"
#include "ExponentialRand.h"
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
    static double stairWidth[257], stairHeight[256];
    static constexpr double x1 = 3.6541528853610088;
    static const bool dummy;
    /// Set up ziggurat tables
    static bool setupTables();

public:
    NormalRand(double mean = 0, double var = 1);
    virtual std::string name() override;

    void setMean(double mean);
    void setSigma(double rootVar);
    void setVar(double var) { setSigma(std::sqrt(std::max(var, MIN_POSITIVE))); }
    double getMean() const { return mu; }
    double getSigma() const { return sigma; }
    double getVar() const { return sigma * sigma; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    static double variate(double mean, double rootVar);
    static double standardVariate();

    double E() const { return mu; }
    double Var() const { return sigma * sigma; }

    double quantile(double p);

    inline double Median() const { return mu; }
    inline double Mode() const { return mu; }
    static constexpr double Skewness() { return 0; }
    static constexpr double ExcessKurtosis() { return 0; }

    bool fitToData(const QVector<double> &sample);
};

#endif // NORMALRAND_H
