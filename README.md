# RandLib
Stochastic calculus

With RandLib one can work with univariate distributions.
* Fast sampling. For instance, generate ten thousand variates from standard normal distribution:
```c++
NormalRand randomVariable(0, 1);
std::vector<double> data(10000);
randomVariable.sample(data);
```
* Calculate moments:
```c++
LaplaceRand randomVariable(2, 3);
std::cout << " Mean = " << randomVariable.E()
          << " Variance = " << randomVariable.Var()
          << " Skewness = " << randomVariable.Skewness()
          << " Ex. kurtosis = " << randomVariable.ExcessKurtosis();
```
* Calculate probabilities for discrete distributions and density functions for continuous:
```c++
GeometricRand geometricRandomVariable(4);
std::cout << "Probability to get 5 for Geometric(4) is " << geometricRandomVariable.P(5);
ExponentialRand expRandomVariable(4);
std::cout << "Probability density function at point 5 for Exponential(4) is " << expRandomVariable.f(5);
```
* Get cumulative density function for random variables with sophisticated distribution:
```c++
BetaRand randomVariable(6, 7);
std::vector<double> data(100);
randomVariable.cdf(data);
for (const auto& i : data)
  std::cout << i << " ";
```
