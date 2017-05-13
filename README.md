# RandLib

[![Build Status](https://travis-ci.org/StochasticEngineer/RandLib.svg?branch=master)](https://travis-ci.org/StochasticEngineer/RandLib)
<a href="https://scan.coverity.com/projects/randlib">
  <img alt="Coverity Scan Build Status"
       src="https://scan.coverity.com/projects/12703/badge.svg"/>
</a>

With RandLib one can easily work with probability distributions.
* Fast sampling. For instance, generate million variates from standard normal distribution:
```c++
NormalRand distribution(0, 1);
std::vector<double> data(1e6);
distribution.Sample(data);
```
![alt tag](https://github.com/StochasticEngineer/RandLib/blob/master/images/standardNormal.png)

* Calculate moments and other properties:
```c++
LogNormalRand distribution(1, 1);
std::cout << " Mean = " << distribution.Mean()
          << " and Variance = " << distribution.Variance()
          << "\n Median = " << distribution.Median()
          << " and Mode = " << distribution.Mode()
          << "\n Skewness = " << distribution.Skewness()
          << " and Excess kurtosis = " << distribution.ExcessKurtosis();
```
![alt tag](https://github.com/StochasticEngineer/RandLib/blob/master/images/lognormal11.png)
```
Mean = 4.48169 and Variance = 34.5126
Median = 2.71828 and Mode = 1
Skewness = 6.18488 and Excess Kurtosis = 110.936
```
* Fitting parameters using different estimators:
```c++
using std::cout;

NormalRand X(0, 1);
std::vector<double> data(10);
X.Sample(data);
cout << "True distribution: " << X.Name() << "\n";
cout << "Sample: ";
for (double var : data)
    cout << var << "  ";
cout << "\n";

/// Bayesian estimation
NormalInverseGammaRand prior(0, 1, 1, 1);
NormalInverseGammaRand posterior = X.FitMeanAndVarianceBayes(data, prior);
cout << "Bayesian estimator: " << X.Name() << "\n";
cout << "(Posterior distribution: " << posterior.Name() << ")\n";

/// Uniformly minimum variance unbiased estimator
X.FitMeanAndVarianceUMVU(data);
cout << "UMVU estimator: " << X.Name() << "\n";

/// Maximum-likelihood estimator
X.FitMeanAndVarianceMLE(data);
cout << "Maximum-likelihood estimator: " << X.Name() << "\n";
```
![alt tag](https://github.com/StochasticEngineer/RandLib/blob/master/images/normalFit.png)
```
True distribution: Normal(0, 1)
Sample: -0.328154  0.709122  -0.607214  1.11472  -1.23726  -0.123584  0.59374  -1.20573  -0.397376  -1.63173
Bayesian estimator: Normal(-0.283042, 0.951348)
(Posterior distribution: Normal-Inverse-Gamma(-0.283042, 11, 6, 4.75674))
UMVU estimator: Normal(-0.311347, 0.82504)
Maximum-likelihood estimator: Normal(-0.311347, 0.742536)
```

Implemented distributions (under titles special cases are listed):
* Continuous
  * Beta distribution
    * Arcsine distribution
    * Balding-Nichols distribution
    * Uniform distribution
  * Beta-prime distribution (Beta of the second kind)
  * Degenerate distribution
  * Exponentially-modified Gaussian distribution (EMG)
  * F-distribution (Fisher-Snedecor)
  * Gamma distribution
    * Chi-Squared distribution
    * Erlang distribution
    * Exponential distribution
  * Generalised extreme value distribution
    * Gumbel distribution
    * Frechet distribution
    * Weibull distribution
  * Geometric Stable distribution
    * Laplace distribution
  * Inverse-Gamma distribution
  * Inverse-Gaussian distribution
  * Irwin-Hall distribution
  * Kolmogorov-Smirnov distribution
  * Kumaraswamy distribution
  * Logistic distribution
  * Log-Normal distribution
  * Nakagami distribution
    * Chi distribution
    * Maxwell-Boltzmann distribution
    * Rayleigh distribution
  * Noncentral Chi-Squared distribution
  * Pareto distribution
  * Planck distribution
  * Raised-cosine distribution
    * Raab-Green distribution
  * Sech distribution
  * Stable distribution
    * Cauchy distribution
    * Holtsmark distribution
    * Landau distribution
    * Levy distribution
    * Normal distribution
  * t-distribution
  * Triangular distribution
  * Von-Mises distribution
  * Wigner-Semicircle distribution
* Disrete:
  * Beta-binomial distribution
  * Binomial distribution
    * Bernoulli distribution
  * Categorical distribution
  * Hypergeometric distribution
  * Logarithmic distribution
  * Negative-binomial (Polya) distribution
    * Pascal distribution
    * Geometric distribution
  * Negative hypergeometric distribution
  * Poisson distribution
  * Rademacher distribution
  * Skellam distribution
  * Uniform discrete distribution
  * Yule distribution
  * Zeta distribution
  * Zipf distribution
* Singular:
  * Cantor distribution
* Bivariate:
  * Bivariate Normal distribution
  * Normal-inverse-Gamma distribution
  
  
