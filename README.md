# RandLib
Stochastic calculus

With RandLib one can work with univariate distributions.
* Fast sampling. For instance, generate ten thousand variates from standard normal distribution:
```c++
NormalRand randomVariable(0, 1);
std::vector<double> data(10000);
randomVariable.sample(data);
```
* Calculate moments and other properties:
```c++
LaplaceRand randomVariable(2, 3);
std::cout << " Mean = " << randomVariable.Mean()
          << " Variance = " << randomVariable.Variance()
          << " Median = " << randomVariable.Median()
          << " Mode = " << randomVariable.Mode()
          << " Skewness = " << randomVariable.Skewness()
          << " Ex. kurtosis = " << randomVariable.ExcessKurtosis();
```
* Calculate probabilities for discrete distributions and probability density functions for continuous:
```c++
GeometricRand geometricRandomVariable(4);
std::cout << "Probability to get 5 for Geometric(4) is " << geometricRandomVariable.P(5);
ExponentialRand expRandomVariable(4);
std::cout << "Density function at point 5 for Exponential(4) is " << expRandomVariable.f(5);
```
* Get cumulative density function for random variables with sophisticated distribution:
```c++
BetaRand randomVariable(6, 7);
int size = 100;
std::vector<double> x(size), y(size);
double sizem1Inv = 1.0 / (size - 1);
for (int i = 0; i != size; ++i)
    x[i] = i * sizem1Inv;
randomVariable.cdf(x, y);
for (int i = 0; i != size; ++i)
    std::cout << "P(X < " << x[i] << ") = " << y[i];
```

List of implemented distributions:

Continuous:

|    Title     |     F(x)     |     f(x)     |   variate    |   CF(t)    |
| ------------ | ------------ | ------------ | ------------ | ------------ |
|    Arcsine   | :white_check_mark: | :white_check_mark: | :white_check_mark: |:x:|
|     Beta     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:x:|
|     Beta Prime     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:x:|
|     Cauchy     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Chi-squared     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Erlang     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Exponential     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     F    | :white_check_mark: | :white_check_mark: | :white_check_mark: |:x:|
|     Gamma     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Geometric Stable     | :x: | :x: | :white_check_mark: |:white_check_mark:|
|     Gumbel     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:x:|
|     Laplace     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Levy     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Log-normal     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:x:|
|     Logistic     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:x:|
|     Maxwell-Boltzmann     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:x:|
|     Nakagami     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:x:|
|     Normal     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Pareto     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:x:|
|     Raised cosine     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Rayleigh     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:x:|
|     Sech (Hyperbolic secant)    | :white_check_mark: | :white_check_mark: | :white_check_mark: |:x:|
|     Stable     | :heavy_exclamation_mark: | :heavy_exclamation_mark: | :white_check_mark: |:white_check_mark:|
|     Student's t     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:x:|
|     Triangular     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Uniform     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     von Mises     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:x:|
|     Wald (Inverse Gaussian)     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Weibull     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:x:|
|     von Wigner Semicircle     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:x:|

Discrete:

|    Title     |     F(x)     |     P(X = x)     |   variate    |   CF(t)    |
| ------------ | ------------ | ------------ | ------------ | ------------ |
|     Bernoulli     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Binomial     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Geometric    | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Hypergeometric     | :white_check_mark: | :white_check_mark: | :heavy_exclamation_mark: |:x:|
|     Logarithmic     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Negative binomial     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Poisson     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Rademacher     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Uniform     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:white_check_mark:|
|     Yule     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:x:|
|     Zeta     | :white_check_mark: | :white_check_mark: | :white_check_mark: |:x:|
|     Zipf     | :white_check_mark: | :white_check_mark: | :x: |:x:|

Singular:

|    Title     |     F(x)     |  variate    |   CF(t)    |
| ------------ | ------------ | ------------ | ------------ |
|     Cantor     | :white_check_mark: | :white_check_mark: | :white_check_mark: |
