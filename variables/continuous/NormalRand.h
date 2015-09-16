#ifndef NORMALRAND_H
#define NORMALRAND_H

#include "ContinuousRand.h"

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
    std::string name() override;

    void setMean(double mean);
    void setSigma(double rootVar);
    void setVar(double var);
    inline double getMean() const { return mu; }
    inline double getSigma() const { return sigma; }
    inline double getVar() const { return sigma * sigma; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    static double variate(double mean, double rootVar);
    static double standardVariate();

    double Mean() const { return mu; }
    double Variance() const { return sigma * sigma; }

    std::complex<double> CF(double t) const override;
    double Quantile(double p) const override;

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    bool fitToData(const QVector<double> &sample);
};

#endif // NORMALRAND_H
