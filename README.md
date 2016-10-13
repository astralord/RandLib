# RandLib
Stochastic calculus

With RandLib one can work with probability distributions.
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
* Fitting parameters:
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

/// Bayesian
NormalInverseGammaRand prior(0, 1, 1, 1);
X.FitBayes(data, prior);
cout << "Bayesian estimator: " << X.Name() << "\n";
cout << "(Posterior distribution: " << prior.Name() << ")\n";

/// UMVU
X.FitUMVU(data);
cout << "UMVU estimator: " << X.Name() << "\n";

/// Maximum-likelihood
X.FitMLE(data);
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

List of available distributions:

:white_check_mark: - ok

:warning: - numerically unstable and/or can take too much time in some cases

:x: - not yet implemented



Continuous distributions:

|    Title     |     F(x)     |     f(x)     |   variate    |   CF(t)    |
| ------------ | ------------ | ------------ | ------------ | ------------ |
|    Arcsine   | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark: |
|     Beta     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark: |
|     Beta Prime     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:warning:|
|     Cauchy     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Chi    | :white_check_mark: | :white_check_mark: | :white_check_mark: |:warning:|
|     Chi-squared     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Erlang     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Exponential     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Exponential-Normal (EMG)     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     F    | :white_check_mark: | :white_check_mark: | :white_check_mark: |:warning:|
|     Frechet    | :white_check_mark: | :white_check_mark: | :white_check_mark: |:warning:|
|     Gamma     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Geometric Stable     | :warning: | :warning: | :white_check_mark: |:white_check_mark:|
|     Gumbel     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:warning:|
|     Irwin-Hall     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Laplace     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Levy     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Log-Cauchy    | :white_check_mark: | :white_check_mark: | :white_check_mark: |:warning:|
|     Log-normal     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:warning:|
|     Logistic     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Maxwell-Boltzmann     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:warning:|
|     Nakagami     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:warning:|
|     Non-central chi-squared | :warning: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Normal     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Pareto     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:warning:|
|     Planck     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:warning:|
|     Raised cosine     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Rayleigh     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:warning:|
|     Sech (Hyperbolic secant)    | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Stable     | :warning: | :warning: | :white_check_mark: |:white_check_mark:|
|     Student's t     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:warning:|
|     Triangular     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Uniform     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     von Mises     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:warning:|
|     Wald (Inverse Gaussian)     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Weibull     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:warning:|
|     Wigner Semicircle     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|

Discrete distributions:

|    Title     |     F(x)     |     P(X = x)     |   variate    |   CF(t)    |
| ------------ | ------------ | ------------ | ------------ | ------------ |
|     Bernoulli     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Beta-Binomial     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Binomial     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Geometric    | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Hypergeometric     | :white_check_mark: | :white_check_mark: | :warning: |:white_check_mark:|
|     Logarithmic     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Negative binomial     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Poisson     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Rademacher     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Skellam    | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Uniform     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Yule     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:warning:|
|     Zeta     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:warning:|
|     Zipf     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|

Singular distributions:

|    Title     |     F(x)     |  variate    |   CF(t)    |
| ------------ | ------------ | ------------ | ------------ |
|     Cantor     | :white_check_mark: | :white_check_mark: | :white_check_mark: |
