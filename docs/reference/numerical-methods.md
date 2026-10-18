# Shared numerical methods {#numerical_methods_reference}

Distribution entries describe their specialized algorithms. This page records the common kernels
used across the library, so an entry can state which kernel it uses without repeating its full
derivation.

## Stable evaluation of probabilities

RandLib evaluates logarithms before exponentiating whenever a density, mass or likelihood contains
large products or gamma-function ratios. Expressions such as \f$1-e^x\f$ and \f$\log(1-e^x)\f$ use
`expm1`, `log1p` and a dedicated `log1mexp` split to avoid cancellation. Regularized incomplete beta
and gamma functions return lower and upper probabilities together, including their logarithms, so a
small upper tail is not obtained by subtracting a rounded lower tail from one.

The same rule applies at the distribution layer: `S(x)` is a first-class probability function.
Quantile searches choose `F` or `S` according to the requested probability and use logarithmic tail
comparisons where an ordinary probability would underflow.

## Finite and infinite quadrature

The general integrator is a global adaptive 21-point Gauss--Kronrod rule with its embedded 10-point
Gauss estimate. On a segment \f$[a,b]\f$, the Kronrod approximation \f$K_{21}\f$ supplies the value
and

\f[
 \varepsilon_{a,b}=\lvert K_{21}-G_{10}\rvert
\f]

supplies the local error indicator. A heap repeatedly bisects the segment with the largest indicator.
Iteration stops when the accumulated estimate satisfies the requested relative tolerance or reaches
a roundoff floor proportional to the estimated integral of \f$|f|\f$.

Semi-infinite intervals use the QUADPACK-style rational map

\f[
 x=a+\frac{1-t}{t},\qquad dx=-\frac{dt}{t^2},\qquad 0<t\le 1,
\f]

with the analogous reflection for a finite upper limit. A doubly infinite integral evaluates the
positive and negative branches together after the same mapping. The rule never samples an infinite
endpoint.

## Endpoint singularities

Finite integrals with algebraic or logarithmic endpoint singularities use double-exponential
tanh--sinh quadrature. With midpoint \f$c=(a+b)/2\f$ and half-width \f$h=(b-a)/2\f$, the transform is

\f[
 x(u)=c+h\tanh\!\left(\frac{\pi}{2}\sinh u\right).
\f]

The transformed weights decay double exponentially. Successive step halvings are compared against
the integral of the absolute transformed integrand. Specialized callables may also receive the
accurately computed distance to the nearer endpoint, preventing cancellation when `x` itself rounds
to that endpoint. Stable and log-stable probability integrals use this facility.

## Oscillatory Fourier integrals

Characteristic functions without a practical closed form use a QAWF-style semi-infinite Fourier
integrator. For positive frequency \f$\omega\f$, the domain is split into half-periods of length
\f$\pi/\omega\f$. Each cosine and sine contribution is evaluated by adaptive Gauss--Kronrod, and
the partial sums

\f[
 s_n=\sum_{k=0}^{n}\int_{a+k\pi/\omega}^{a+(k+1)\pi/\omega}
       g(x)e^{i\omega x}\,dx
\f]

are accelerated with Wynn's epsilon algorithm. Three negligible consecutive cycles terminate a
directly convergent tail; otherwise the best stable even-order epsilon estimate is returned. This is
the production path for several heavy-tailed characteristic functions.

## Quantiles and scalar roots

The generic continuous quantile implementation first brackets the requested probability by expanding
away from the median or another distribution-specific seed. It then applies Alefeld--Potra--Shi ACM
TOMS Algorithm 748. The method maintains a sign-changing bracket and combines inverse-cubic,
Newton-quadratic and secant interpolation. Bisection is forced whenever interpolation does not shrink
the bracket sufficiently, so the search retains bracketing robustness while converging faster on
smooth probability functions.

Specialized parameter estimators also use safeguarded first- and second-order Newton iterations. They
accept an update only while it remains finite and inside the parameter domain; otherwise they
backtrack or report non-convergence rather than evaluating a likelihood at an invalid parameter.

## Optimization and fitting

One-dimensional smooth objectives use Brent minimization: parabolic interpolation is accepted when
it lies safely inside the current bracket, otherwise a golden-section step is taken.

The generic constrained MLE maps bounded parameters to unconstrained coordinates, then minimizes the
negative log-likelihood. Its derivative-free fallback is the Gao--Han dimension-adaptive
Nelder--Mead simplex. BFGS is available for smooth objectives and uses a central-difference gradient,
a backtracking Armijo line search and the standard inverse-Hessian rank-two update. Numerical
gradients use steps proportional to \f$\epsilon^{1/3}\f$; observed-information Hessians use symmetric
central differences with steps proportional to \f$\epsilon^{1/4}\f$.

For fitted parameters \f$\hat\theta\f$, the covariance estimate is the inverse observed information

\f[
 \widehat{\operatorname{Cov}}[\hat\theta]
   =\left[-\nabla^2\ell(\hat\theta)\right]^{-1},
\f]

when that matrix is finite and invertible. Closed-form and distribution-specific solvers override the
generic path where they are more accurate or substantially faster.

## Standard normal and exponential variates

The shared standard-normal and unit-exponential samplers use 256-layer Marsaglia--Tsang ziggurats.
The normal sampler consumes an explicit sign bit and uses Marsaglia's exponential-tail rejection
method below the base strip. The exponential base strip recurses beyond its table boundary. Rejection
loops have a defensive iteration bound that signals an implementation failure; they do not truncate
the mathematical distribution.

All random bits come from the caller's standard uniform random bit generator. Distribution objects
therefore contain no hidden mutable engine or sampling cache, and reproducibility is controlled by
the caller-owned engine state.

## Online sample moments

The common sample-statistics helper uses the one-pass Pébay/Welford recurrences. For the \f$n\f$th
observation \f$x\f$, define \f$\delta=x-M_1\f$, \f$\delta_n=\delta/n\f$ and
\f$t=\delta\delta_n(n-1)\f$. The unnormalized central sums update as

\f[
\begin{aligned}
 M_1' &= M_1+\delta_n,\\
 M_4' &= M_4+t\delta_n^2(n^2-3n+3)+6\delta_n^2M_2-4\delta_nM_3,\\
 M_3' &= M_3+t\delta_n(n-2)-3\delta_nM_2,\\
 M_2' &= M_2+t.
\end{aligned}
\f]

The implementation periodically merges blocks with the corresponding parallel central-moment
formulas, limiting loss of precision for long samples. Reported population moments are
\f$M_r/n\f$; standardized skewness and excess kurtosis are derived from those central sums.
