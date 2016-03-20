#ifndef GEOMETRICRAND_H
#define GEOMETRICRAND_H

#include "NegativeBinomialRand.h"
#include "../continuous/ExponentialRand.h"
#include "../continuous/BetaRand.h"

/**
 * @brief The GeometricRand class
 * Geometric distribution
 * X ~ Geometric(p)
 *
 * P(X = k) = p (1 - p)^k
 *
 * X ~ NB(1, p)
 */
class RANDLIBSHARED_EXPORT GeometricRand : public PascalRand
{
public:
    explicit GeometricRand(double probability = 0.5);
    std::string name() override;

protected:
    using NegativeBinomialRand::setParameters;

public:
    void setProbability(double probability);

    double P(int k) const override;
    double F(double x) const override;
    double variate() const override;
    static double variate(double probability);

    void sample(QVector<double> &outputData) const override;

    double Median() const override;

    double Entropy() const;
    
    bool checkValidity(const QVector<double> &sample);

    bool fit_MLE(const QVector<double> &sample);
    bool fit_MM(const QVector<double> &sample);
    bool fitProbability_Bayes(const QVector<double> &sample, BetaRand &priorDistribution);
};

#endif // GEOMETRICRAND_H
