# Continuous distributions {#continuous_distribution_reference}

This page documents every concrete continuous distribution in RandLib. Related special cases share
one section so that their parameter transformations and common numerical implementation remain
visible. The detailed entries are being migrated from the archived notes and checked against the v2
headers; no legacy accuracy statement is retained without current code or test evidence.

## Coverage

The continuous reference covers the following public classes:

`ArcsineRand`, `AsymmetricLaplaceRand`, `BaldingNicholsRand`, `BetaPrimeRand`, `BetaRand`,
`CauchyRand`, `CenteredAsymmetricLaplaceRand`, `ChiRand`, `ChiSquaredRand`, `DegenerateRand`,
`ErlangRand`, `ExponentialRand`, `ExponentiallyModifiedGaussianRand`, `FisherFRand`, `FrechetRand`,
`GammaRand`, `GeneralizedExtremeValueRand`, `GeneralizedParetoRand`, `GeometricStableRand`,
`GumbelRand`, `HoltsmarkRand`, `InverseGammaRand`, `InverseGaussianRand`, `IrwinHallRand`,
`KolmogorovSmirnovRand`, `LandauRand`, `LaplaceRand`, `LevyRand`, `LogNormalRand`, `LogStableRand`,
`LogisticRand`, `MarchenkoPasturRand`, `MaxwellBoltzmannRand`, `NakagamiRand`,
`NoncentralChiSquaredRand`, `NoncentralTRand`, `NormalRand`, `ParetoRand`, `PlanckRand`,
`RaabGreenRand`, `RaisedCosineRand`, `RayleighRand`, `RiceRand`, `SechRand`, `StableRand`,
`StudentTRand`, `TriangularRand`, `UniformRand`, `WeibullRand` and `WignerSemicircleRand`.

## Beta family

### Beta distribution (`BetaRand`)

`BetaRand(alpha, beta, a, b)` uses shapes \f$\alpha,\beta>0\f$ and finite support
\f$[a,b]\f$, where \f$a<b\f$. Write

\f[
 d=b-a,\qquad y=\frac{x-a}{d}.
\f]

On the support its density and distribution functions are

\f[
 f(x)=\frac{y^{\alpha-1}(1-y)^{\beta-1}}{dB(\alpha,\beta)},\qquad
 F(x)=I_y(\alpha,\beta),\qquad
 S(x)=I_{1-y}(\beta,\alpha),
\f]

where \f$B\f$ is the beta function and \f$I\f$ is the regularized incomplete beta function. The
direct expression for \f$S\f$ is also the implementation path; the upper tail is not formed as
`1 - F(x)`. The characteristic function has the closed form

\f[
 \phi(t)=e^{ita}\,{}_1F_1\!\left(\alpha;\alpha+\beta;itd\right).
\f]

RandLib evaluates it by quadrature instead of evaluating the confluent hypergeometric function. For
\f$\alpha,\beta\geq1\f$ it integrates \f$e^{itx}f(x)\f$ directly. When either endpoint is singular,
the implementation subtracts its limiting trigonometric value before integrating. For example, with
\f$z=td\f$ the standardized real part is evaluated from

\f[
 \operatorname{Re}\phi_0(z)
 =1+\frac{1}{B(\alpha,\beta)}
 \left[
 \int_0^1 (1-x)^{\beta-1}
 \left(x^{\alpha-1}(\cos zx-1)-(\cos z-1)\right)dx
 +\frac{\cos z-1}{\beta}
 \right].
\f]

This algebraically removes both endpoint singularities before the shared adaptive
Gauss--Kronrod integrator is called. The imaginary part uses the analogous subtraction with
\f$\sin(zx)\f$ and \f$\sin z\f$, after which multiplication by \f$e^{ita}\f$ restores the location.

The principal moments are

\f[
 \operatorname{E}X=a+d\frac{\alpha}{\alpha+\beta},\qquad
 \operatorname{Var}X=d^2\frac{\alpha\beta}
 {\left(\alpha+\beta\right)^2\left(\alpha+\beta+1\right)},
\f]

\f[
 \operatorname{Skew}X=
 \frac{2(\beta-\alpha)\sqrt{\alpha+\beta+1}}
 {(\alpha+\beta+2)\sqrt{\alpha\beta}},
\f]

and

\f[
 \operatorname{ExKurt}X=
 \frac{6\left[(\alpha-\beta)^2(\alpha+\beta+1)
 -\alpha\beta(\alpha+\beta+2)\right]}
 {\alpha\beta(\alpha+\beta+2)(\alpha+\beta+3)}.
\f]

The mean absolute deviation about the mean is available as

\f[
 \operatorname{E}|X-\operatorname{E}X|
 =\frac{2d\,\alpha^\alpha\beta^\beta}
 {B(\alpha,\beta)(\alpha+\beta)^{\alpha+\beta+1}}.
\f]

`GeometricMean()` and `GeometricVariance()` retain historical API names but return the log-moments
of the standardized variable \f$Y=(X-a)/d\f$:

\f[
 \operatorname{E}\log Y=\psi(\alpha)-\psi(\alpha+\beta),\qquad
 \operatorname{Var}(\log Y)=\psi_1(\alpha)-\psi_1(\alpha+\beta).
\f]

For \f$\alpha,\beta>1\f$ the unique mode is
\f$a+d(\alpha-1)/(\alpha+\beta-2)\f$. If exactly one shape is at most one, the corresponding
endpoint is a mode. If both are below one, both endpoints are locally modal; `Mode()` returns one
representative endpoint (the upper endpoint when \f$\alpha>\beta\f$, otherwise the lower endpoint).
For \f$\alpha=\beta=1\f$ every point is a mode and the API likewise returns the lower endpoint.

The median is analytic for the symmetric case, for \f$\alpha=1\f$, and for \f$\beta=1\f$:

\f[
 m=\frac{a+b}{2}\quad(\alpha=\beta),\qquad
 m=a+d\left(1-2^{-1/\beta}\right)\quad(\alpha=1),\qquad
 m=a+d\,2^{-1/\alpha}\quad(\beta=1).
\f]

Other medians and quantiles use the generic bracketed TOMS 748 inversion. When both shapes are at
least one, the median search starts at the approximation

\f[
 a+d\frac{\alpha-1/3}{\alpha+\beta-2/3}.
\f]

Quantiles themselves have direct formulas for the symmetric arcsine and uniform cases and whenever
one shape is one.

Sampling is selected by the shape region. RandLib uses direct uniform and arcsine transforms for
\f$(1,1)\f$ and \f$(1/2,1/2)\f$; specialized uniform rejection for symmetric shapes between one and
two; normal rejection for larger symmetric shapes; the Cheng, Johnk, and Atkinson--Whittaker
algorithms in their applicable regions; and otherwise the ratio
\f$G_\alpha/(G_\alpha+G_\beta)\f$ of independent gamma variates. The result is then shifted and
scaled to \f$[a,b]\f$.

With support fixed at \f$[0,1]\f$, `fit<BetaRand>` performs maximum likelihood for the two shapes.
For a standardized sample \f$y_i\f$, the score equations are

\f[
 \overline{\log y}+\psi(\alpha+\beta)-\psi(\alpha)=0,
 \qquad
 \overline{\log(1-y)}+\psi(\alpha+\beta)-\psi(\beta)=0.
\f]

Their two-dimensional Newton solve is initialized by the method of moments, falling back to
\f$(0.001,0.001)\f$ when the sample moments do not yield positive shapes. The class also exposes the
two-parameter exponential-family representation with sufficient statistic
\f$(\log Y,\log(1-Y))\f$ and natural parameter \f$(\alpha-1,\beta-1)\f$.

### Arcsine subfamily (`ArcsineRand`)

`ArcsineRand(shape, a, b)` accepts \f$0<q<1\f$ and is the Beta subfamily

\f[
 X\sim\operatorname{Beta}(1-q,q,a,b).
\f]

Thus, for \f$0<y<1\f$,

\f[
 f(x)=\frac{\sin(\pi q)}{\pi d}\,y^{-q}(1-y)^{q-1},\qquad
 F(x)=I_y(1-q,q),
\f]

\f[
 \phi(t)=e^{ita}\,{}_1F_1(1-q;1;itd),\qquad
 \operatorname{E}X=a+d(1-q),\qquad
 \operatorname{Var}X=\frac{d^2q(1-q)}{2}.
\f]

Only the default \f$q=1/2\f$ is the classical symmetric arcsine distribution. In that case,

\f[
 F(x)=\frac{2}{\pi}\arcsin\sqrt y,\qquad
 S(x)=\frac{2}{\pi}\arccos\sqrt y,\qquad
 Q(p)=a+d\sin^2\!\left(\frac{\pi p}{2}\right),
\f]

and RandLib uses these formulas directly, including the inverse-CDF sampler. Other shapes inherit
the Beta probability, quadrature, quantile, and sampling paths.

`fit<ArcsineRand>` estimates \f$q\f$ by maximum likelihood on the standard support. If
\f$L_0=\overline{\log Y}\f$ and \f$L_1=\overline{\log(1-Y)}\f$, the score equation is

\f[
 L_1-L_0+\pi\cot(\pi q)=0,
\f]

which the implementation solves in closed form with an `atan` expression and the appropriate branch
in \f$(0,1)\f$.

### Balding--Nichols distribution (`BaldingNicholsRand`)

`BaldingNicholsRand(fixationIndex, frequency)` accepts \f$F,p\in(0,1)\f$ and is the standard-support
Beta law

\f[
 X\sim\operatorname{Beta}\!\left(p\frac{1-F}{F},
 (1-p)\frac{1-F}{F}\right).
\f]

Consequently \f$\operatorname{E}X=p\f$ and \f$\operatorname{Var}X=Fp(1-p)\f$; its density, CDF,
characteristic function, higher moments, quantiles, and sampler are the corresponding Beta formulas
and algorithms above. The constructor order is \f$(F,p)\f$, despite the conventional notation often
listing frequency first. There is currently no `fit<BaldingNicholsRand>` specialization.

### Uniform distribution (`UniformRand`)

`UniformRand(a, b)` is `BetaRand(1, 1, a, b)`. On \f$[a,b]\f$,

\f[
 f(x)=\frac1d,\qquad F(x)=\frac{x-a}{d},\qquad S(x)=\frac{b-x}{d},
\f]

with the usual constant continuation outside the support, and

\f[
 \phi(t)=
 \begin{cases}
  \dfrac{e^{itb}-e^{ita}}{itd},&t\ne0,\\
  1,&t=0.
 \end{cases}
\f]

Its mean, variance, median, skewness, excess kurtosis, and differential entropy are

\f[
 \operatorname{E}X=\operatorname{Median}X=\frac{a+b}{2},\qquad
 \operatorname{Var}X=\frac{d^2}{12},\qquad
 \operatorname{Skew}X=0,
\f]

\f[
 \operatorname{ExKurt}X=-\frac65,\qquad h(X)=\log d.
\f]

Every support point is a mode; `Mode()` returns the midpoint as a representative value. Sampling is
the affine transform \f$a+dU\f$ of the caller-engine uniform generator. `StandardVariate`,
`StandardVariateClosed`, and `StandardVariateHalfClosed` expose respectively the open, closed, and
left-closed/right-open unit-interval conventions.

`fit<UniformRand>` supports maximum likelihood and UMVU estimation. With order statistics
\f$x_{(1)}\f$ and \f$x_{(n)}\f$, the MLE is

\f[
 \hat a=x_{(1)},\qquad \hat b=x_{(n)},
\f]

while the unbiased estimators are

\f[
 \widetilde a=\frac{nx_{(1)}-x_{(n)}}{n-1},\qquad
 \widetilde b=\frac{nx_{(n)}-x_{(1)}}{n-1}.
\f]

When \f$a\f$ is fixed and the width \f$\theta=b-a\f$ has prior
\f$\operatorname{Pareto}(\alpha,\sigma)\f$, `fit_bayes` returns the conjugate posterior

\f[
 \theta\mid x\sim\operatorname{Pareto}\!\left(
 \alpha+n,\max\{\sigma,x_{(n)}-a\}\right).
\f]

The posterior mean or mode supplies the Bayesian or MAP point estimate, respectively.

## Gamma family

### Gamma distribution (`GammaRand`)

`GammaRand(shape, rate)` uses the rate parameterization
\f$\alpha>0,\beta>0\f$. The named factories `FromShapeRate` and `FromShapeScale` make the
rate-versus-scale convention explicit; the scale is \f$\theta=1/\beta\f$. For \f$x\geq0\f$,

\f[
 f(x)=\frac{\beta^\alpha}{\Gamma(\alpha)}x^{\alpha-1}e^{-\beta x},\qquad
 F(x)=P(\alpha,\beta x),\qquad
 S(x)=Q(\alpha,\beta x),
\f]

where \f$P\f$ and \f$Q\f$ are the regularized lower and upper incomplete gamma functions. RandLib
has direct implementations of all four of \f$F,S,\log F,\log S\f$; small tails therefore do not
depend on subtracting from one. The characteristic function is

\f[
 \phi(t)=\left(1-\frac{it}{\beta}\right)^{-\alpha}.
\f]

The moments and principal summaries are

\f[
 \operatorname{E}X=\frac{\alpha}{\beta},\qquad
 \operatorname{Var}X=\frac{\alpha}{\beta^2},\qquad
 \operatorname{Mode}X=\frac{\max\{\alpha-1,0\}}{\beta},
\f]

\f[
 \operatorname{Skew}X=\frac{2}{\sqrt\alpha},\qquad
 \operatorname{ExKurt}X=\frac{6}{\alpha}.
\f]

The historical methods `GeometricMean()` and `GeometricVariance()` return log-moments:

\f[
 \operatorname{E}\log X=\psi(\alpha)-\log\beta,
 \qquad \operatorname{Var}(\log X)=\psi_1(\alpha).
\f]

The differential entropy exposed through the exponential-family API is

\f[
 h(X)=\alpha-\log\beta+\log\Gamma(\alpha)+(1-\alpha)\psi(\alpha).
\f]

Except at \f$\alpha=1\f$, where the exponential formula applies, medians and quantiles are found by
inverting the appropriate incomplete-gamma tail. Initial approximations follow Gil, Segura, and
Temme: a power series seeds small lower-tail probabilities, a Lambert-\f$W_{-1}\f$ asymptotic seeds
large values, and an inverse-complementary-error-function plus Lambert-\f$W\f$ approximation handles
large shapes. The result is refined by the shared bracketed TOMS 748 solver. Upper quantiles use
\f$Q\f$ and \f$\log Q\f$ directly.

The sampler dispatches by shape:

- integer shapes one through three are sums of unit exponential variates;
- \f$\alpha=3/2\f$ is an exponential plus half the square of a standard normal;
- very small shapes use Best's RGS algorithm in a configured object and Ahrens--Dieter GS in the
  static `StandardVariate` path;
- \f$1<\alpha<1.2\f$ uses Fishman's rejection algorithm;
- the remaining region uses Marsaglia--Tsang rejection.

The accepted unit-rate variate is divided by \f$\beta\f$. All rejection algorithms consume the
caller-owned engine and enforce the library's defensive iteration limit.

`fit<GammaRand>` performs maximum likelihood for shape and rate on positive samples. If
\f$\bar x\f$ and \f$\overline{\log x}\f$ are the sample means, the profiled rate is
\f$\hat\beta=\hat\alpha/\bar x\f$ and the remaining equation is

\f[
 \log\hat\alpha-\psi(\hat\alpha)
 =\log\bar x-\overline{\log x}.
\f]

RandLib initializes \f$\hat\alpha\f$ with a closed approximation and solves this equation by
safeguarded Newton iteration using \f$1/\alpha-\psi_1(\alpha)\f$. With a known shape, the rate MLE
is \f$\alpha/\bar x\f$ and its UMVU alternative is
\f$(\alpha-1/n)/\bar x\f$.

For fixed shape \f$\alpha\f$, a Gamma prior
\f$\beta\sim\operatorname{Gamma}(\kappa,\gamma)\f$ is conjugate to the rate:

\f[
 \beta\mid x\sim\operatorname{Gamma}\!\left(
 \alpha n+\kappa,\gamma+\sum_i x_i\right).
\f]

`fit_bayes<ExponentialRand>` exposes the \f$\alpha=1\f$ case. The class also exposes the
two-parameter exponential-family representation with sufficient statistic
\f$(\log X,X)\f$ and natural parameter \f$(\alpha-1,-\beta)\f$.

### Chi-squared distribution (`ChiSquaredRand`)

`ChiSquaredRand(degree)` takes a positive integer \f$k\f$ and is

\f[
 X\sim\chi_k^2=\operatorname{Gamma}\!\left(\frac{k}{2},\frac12\right).
\f]

Its density, CDF, survival function, characteristic function, quantiles, and sampler are inherited
from the Gamma implementation. In particular,

\f[
 \phi(t)=(1-2it)^{-k/2},\qquad
 \operatorname{E}X=k,\qquad
 \operatorname{Var}X=2k,
\f]

\f[
 \operatorname{Mode}X=\max\{k-2,0\},\qquad
 \operatorname{Skew}X=\sqrt{\frac8k},\qquad
 \operatorname{ExKurt}X=\frac{12}{k}.
\f]

It is also the distribution of \f$\sum_{j=1}^{k}Z_j^2\f$ for independent standard normals.
`fit<ChiSquaredRand>` treats the degree as discrete: it evaluates the log likelihood over positive
integer degrees up to an adaptive, bounded search limit and returns the maximizer. A fixed degree
bypasses the search.

### Erlang distribution (`ErlangRand`)

`ErlangRand(shape, rate)` restricts the Gamma shape to a positive integer \f$k\f$:

\f[
 X\sim\operatorname{Erlang}(k,\beta)
 =\operatorname{Gamma}(k,\beta)
 \stackrel{d}{=}\sum_{j=1}^{k}E_j,
 \qquad E_j\sim\operatorname{Exponential}(\beta).
\f]

It otherwise uses all Gamma formulas and algorithms. For integer shape the sampler explicitly sums
\f$k\f$ exponential variates when \f$k\leq3\f$ and uses the general Gamma dispatch for larger
values. `fit<ErlangRand>` profiles the rate as \f$k/\bar x\f$ for every candidate integer shape and
selects the largest likelihood. The candidate ceiling is derived from the moment estimate
\f$\bar x^2/s^2\f$ and bounded by the implementation; fixed shape and rate constraints are honored.

### Exponential distribution (`ExponentialRand`)

`ExponentialRand(rate)` is the \f$\alpha=1\f$ Gamma law. The factory `FromScale` accepts its mean
\f$\theta=1/\beta\f$. For \f$x\geq0\f$,

\f[
 f(x)=\beta e^{-\beta x},\qquad
 F(x)=1-e^{-\beta x},\qquad
 S(x)=e^{-\beta x},
\f]

\f[
 Q(p)=-\frac{\log(1-p)}{\beta},\qquad
 \phi(t)=\frac{\beta}{\beta-it}.
\f]

RandLib evaluates \f$F\f$ with `expm1`, evaluates \f$S\f$ directly, and uses `log1p` for the lower
quantile so both tails retain precision. Sampling scales the shared 256-layer unit-exponential
ziggurat by \f$1/\beta\f$.

The mean, variance, median, mode, skewness, excess kurtosis, entropy, and raw integer moments are

\f[
 \operatorname{E}X=\frac1\beta,\quad
 \operatorname{Var}X=\frac1{\beta^2},\quad
 \operatorname{Median}X=\frac{\log2}{\beta},\quad
 \operatorname{Mode}X=0,
\f]

\f[
 \operatorname{Skew}X=2,\qquad
 \operatorname{ExKurt}X=6,\qquad
 h(X)=1-\log\beta,\qquad
 \operatorname{E}X^r=\frac{r!}{\beta^r}\quad(r=0,1,\ldots).
\f]

`fit<ExponentialRand>` supports the rate MLE \f$1/\bar x\f$ and its UMVU version
\f$(1-1/n)/\bar x\f$. A Gamma prior is conjugate, as described for the Gamma family above.

## Transformed beta and gamma laws

### Beta-prime distribution (`BetaPrimeRand`)

`BetaPrimeRand(alpha, beta)` accepts \f$\alpha,\beta>0\f$ and has support \f$[0,\infty)\f$. Its
density and tail probabilities are

\f[
 f(x)=\frac{x^{\alpha-1}(1+x)^{-\alpha-\beta}}{B(\alpha,\beta)},
 \qquad
 F(x)=I_{x/(1+x)}(\alpha,\beta),
\f]

\f[
 S(x)=I_{1/(1+x)}(\beta,\alpha).
\f]

The direct survival expression is important numerically: for large \f$x\f$, forming
\f$x/(1+x)\f$ can round to one. RandLib instead passes \f$1/(1+x)\f$ and separately computed
logarithms to the incomplete-beta implementation.

The transformation

\f[
 Y=\frac{X}{1+X}\sim\operatorname{Beta}(\alpha,\beta),
 \qquad
 X=\frac{Y}{1-Y}
\f]

drives probability evaluation and sampling. A Beta quantile transformed by \f$y/(1-y)\f$ supplies
the initial guess for a final TOMS 748 inversion; that refinement avoids losing upper-tail accuracy
in the transformation itself. Sampling applies the same transform to the selected Beta sampler.

The mathematical characteristic function can be written with Tricomi's confluent hypergeometric
function,

\f[
 \phi(t)=\frac{\Gamma(\alpha+\beta)}{\Gamma(\beta)}
 U(\alpha,1-\beta,-it).
\f]

RandLib does not evaluate \f$U\f$; it uses the shared half-period Fourier integration and Wynn
epsilon acceleration on the density over \f$[0,\infty)\f$.

The principal moments, with their existence conditions, are

\f[
 \operatorname{E}X=\frac{\alpha}{\beta-1}\quad(\beta>1),
\qquad
 \operatorname{Var}X=\frac{\alpha(\alpha+\beta-1)}
 { (\beta-2)(\beta-1)^2}\quad(\beta>2),
\f]

\f[
 \operatorname{Mode}X=\max\!\left\{\frac{\alpha-1}{\beta+1},0\right\},
\f]

\f[
 \operatorname{Skew}X=
 \frac{2(2\alpha+\beta-1)}{\beta-3}
 \sqrt{\frac{\beta-2}{\alpha(\alpha+\beta-1)}}
 \quad(\beta>3),
\f]

\f[
 \operatorname{ExKurt}X=
 \frac{6}{(\beta-3)(\beta-4)}
 \left[
 \frac{(\beta-1)^2(\beta-2)}{\alpha(\alpha+\beta-1)}+5\beta-11
 \right]
 \quad(\beta>4).
\f]

RandLib returns infinity when a listed moment does not exist. The median is one when
\f$\alpha=\beta\f$ and is otherwise found through the quantile path.

`fit<BetaPrimeRand>` performs MLE by transforming every observation to
\f$Y_i=X_i/(1+X_i)\f$ and applying the Beta score-equation solver to
\f$\overline{\log Y}\f$ and \f$\overline{\log(1-Y)}\f$.

### Fisher F distribution (`FisherFRand`)

`FisherFRand(degree1, degree2)` requires positive integer degrees \f$d_1,d_2\f$. It is a scaled
beta-prime law,

\f[
 \frac{d_1}{d_2}X\sim
 \operatorname{BetaPrime}\!\left(\frac{d_1}{2},\frac{d_2}{2}\right),
\f]

and therefore

\f[
 f(x)=\frac{(d_1/d_2)^{d_1/2}x^{d_1/2-1}}
 {B(d_1/2,d_2/2)\left(1+d_1x/d_2\right)^{(d_1+d_2)/2}},
\f]

\f[
 F(x)=I_{d_1x/(d_1x+d_2)}\!\left(\frac{d_1}{2},\frac{d_2}{2}\right),
 \qquad
 S(x)=I_{d_2/(d_1x+d_2)}\!\left(\frac{d_2}{2},\frac{d_1}{2}\right).
\f]

Its characteristic function is the beta-prime characteristic function evaluated at
\f$td_2/d_1\f$ and is computed by the same oscillatory integrator. Quantiles and samples are likewise
scaled beta-prime results.

The main summaries are

\f[
 \operatorname{E}X=\frac{d_2}{d_2-2}\quad(d_2>2),
\f]

\f[
 \operatorname{Var}X=
 \frac{2d_2^2(d_1+d_2-2)}{d_1(d_2-2)^2(d_2-4)}
 \quad(d_2>4),
\f]

\f[
 \operatorname{Mode}X=
 \begin{cases}
 \dfrac{d_2(d_1-2)}{d_1(d_2+2)},&d_1>2,\\
 0,&d_1\leq2,
 \end{cases}
\f]

\f[
 \operatorname{Skew}X=
 \frac{2d_1+d_2-2}{d_2-6}
 \sqrt{\frac{8(d_2-4)}{d_1(d_1+d_2-2)}}
 \quad(d_2>6),
\f]

\f[
 \operatorname{ExKurt}X=
 \frac{12}{(d_2-6)(d_2-8)}
 \left[
 \frac{(d_2-2)^2(d_2-4)}{d_1(d_1+d_2-2)}+5d_2-22
 \right]
 \quad(d_2>8).
\f]

Undefined higher moments are reported as infinity. `fit<FisherFRand>` treats both degrees as
discrete and searches the likelihood grid from one through 80 for each unfixed degree. This finite
implementation search bound is part of the fitting API behavior, not a property of the law.

### Inverse-gamma distribution (`InverseGammaRand`)

`InverseGammaRand(shape, rate)` uses \f$\alpha,\beta>0\f$ and the reciprocal relationship

\f[
 X=\frac1Y,\qquad Y\sim\operatorname{Gamma}(\alpha,\beta).
\f]

For \f$x>0\f$,

\f[
 f(x)=\frac{\beta^\alpha}{\Gamma(\alpha)}x^{-\alpha-1}e^{-\beta/x},
 \qquad
 F(x)=Q\!\left(\alpha,\frac\beta x\right),
 \qquad
 S(x)=P\!\left(\alpha,\frac\beta x\right).
\f]

RandLib evaluates those swapped Gamma tails directly, uses reciprocal Gamma quantiles, and samples
by taking the reciprocal of a Gamma variate. The mathematical characteristic function is

\f[
 \phi(t)=\frac{2}{\Gamma(\alpha)}
 \left(\sqrt{-i\beta t}\right)^\alpha
 K_\alpha\!\left(2\sqrt{-i\beta t}\right),
\f]

with the continuous value one at \f$t=0\f$. The implementation instead uses the shared oscillatory
Fourier integrator.

The moments and shape summaries are

\f[
 \operatorname{E}X=\frac\beta{\alpha-1}\quad(\alpha>1),
 \qquad
 \operatorname{Var}X=\frac{\beta^2}{(\alpha-1)^2(\alpha-2)}
 \quad(\alpha>2),
\f]

\f[
 \operatorname{Mode}X=\frac\beta{\alpha+1},qquad
 \operatorname{Median}X=\frac{1}{\operatorname{Median}(Y)},
\f]

\f[
 \operatorname{Skew}X=\frac{4\sqrt{\alpha-2}}{\alpha-3}
 \quad(\alpha>3),
 \qquad
 \operatorname{ExKurt}X=\frac{30\alpha-66}{(\alpha-3)(\alpha-4)}
 \quad(\alpha>4).
\f]

`fit<InverseGammaRand>` initializes \f$\alpha=2+\bar x^2/s^2\f$ and
\f$\beta=\bar x(\alpha-1)\f$, then runs the common constrained maximum-likelihood optimizer. The law
also serves as the conditional conjugate prior for a normal variance when the normal location is
fixed; that model is described with the normal family.

## Normal family

### Normal distribution (`NormalRand`)

`NormalRand(location, variance)` uses mean \f$\mu\in\mathbb R\f$ and variance
\f$\sigma^2>0\f$; `FromMeanStd` is the explicit standard-deviation factory. Its density, tails, and
characteristic function are

\f[
 f(x)=\frac{1}{\sigma\sqrt{2\pi}}
 \exp\!\left[-\frac{(x-\mu)^2}{2\sigma^2}\right],
\f]

\f[
 F(x)=\frac12\operatorname{erfc}\!\left(
 \frac{\mu-x}{\sigma\sqrt2}\right),
 \qquad
 S(x)=\frac12\operatorname{erfc}\!\left(
 \frac{x-\mu}{\sigma\sqrt2}\right),
\f]

\f[
 \phi(t)=\exp\!\left(i\mu t-\frac12\sigma^2t^2\right).
\f]

Both tails are evaluated directly with `erfc`. Quantiles use the corresponding direct inverse-tail
forms

\f[
 Q(p)=\mu-\sigma\sqrt2\,\operatorname{erfc}^{-1}(2p),
 \qquad
 Q(1-p)=\mu+\sigma\sqrt2\,\operatorname{erfc}^{-1}(2p).
\f]

RandLib's inverse error function uses piecewise Chebyshev approximations, with separate central,
intermediate-tail, and extreme-tail coefficient sets. Sampling uses the shared 256-layer normal
ziggurat and then applies \f$\mu+\sigma Z\f$.

The normal law has

\f[
 \operatorname{E}X=\operatorname{Median}X=\operatorname{Mode}X=\mu,
 \qquad
 \operatorname{Var}X=\sigma^2,
\f]

\f[
 \operatorname{Skew}X=0,qquad
 \operatorname{ExKurt}X=0,qquad
 h(X)=\frac12\log(2\pi e\sigma^2).
\f]

Raw integer moments are evaluated from

\f[
 \operatorname{E}X^n=
 \sum_{\substack{0\leq k\leq n\\k\ \mathrm{even}}}
 {n\choose k}\mu^{n-k}\sigma^k(k-1)!!,
\f]

where the \f$k=0\f$ term is one. In RandLib's stable-law convention the same law is
\f$S_1(2,0,\sigma/\sqrt2,\mu)\f$.

`fit<NormalRand>` supports MLE, method of moments, and UMVU estimation. With both parameters free,

\f[
 \hat\mu=\bar x,qquad
 \hat\sigma^2_{\mathrm{MLE}}=\frac1n\sum_i(x_i-\bar x)^2,qquad
 \widetilde{\sigma^2}=\frac1{n-1}\sum_i(x_i-\bar x)^2.
\f]

The legacy member API additionally exposes normal/Student-\f$t\f$ confidence intervals for the mean,
chi-squared confidence intervals for the variance, and the gamma-ratio correction that makes the
sample standard deviation unbiased when the mean is known.

The conjugate models exposed by `fit_bayes` are:

- for known \f$\sigma^2\f$ and prior
  \f$\mu\sim N(\mu_0,1/\tau_0)\f$,
  \f$\tau_n=\tau_0+n/\sigma^2\f$ and
  \f$\mu_n=(\tau_0\mu_0+\sum_i x_i/\sigma^2)/\tau_n\f$;
- for known \f$\mu\f$ and prior
  \f$\sigma^2\sim\operatorname{InvGamma}(\alpha_0,\beta_0)\f$,
  \f$\alpha_n=\alpha_0+n/2\f$ and
  \f$\beta_n=\beta_0+\frac12\sum_i(x_i-\mu)^2\f$;
- with both parameters unknown, a normal--inverse-gamma prior updates to

\f[
 \lambda_n=\lambda_0+n,qquad
 \mu_n=\frac{\lambda_0\mu_0+n\bar x}{\lambda_n},qquad
 \alpha_n=\alpha_0+\frac n2,
\f]

\f[
 \beta_n=\beta_0+\frac n2\left[
 s_n^2+\frac{\lambda_0}{\lambda_n}(\bar x-\mu_0)^2
 \right],qquad
 s_n^2=\frac1n\sum_i(x_i-\bar x)^2.
\f]

The joint prior itself is documented under `NormalInverseGammaRand` in the bivariate reference.

### Log-normal distribution (`LogNormalRand`)

`LogNormalRand(location, squaredScale)` parameterizes the logarithm:

\f[
 Y=\log X\sim N(\mu,\sigma^2).
\f]

The second constructor argument is the log-variance \f$\sigma^2\f$, while `GetScale()` returns the
log-standard deviation \f$\sigma\f$. The explicit factories `FromLogMeanLogStd` and
`FromLogMeanLogVariance` avoid that distinction at call sites.

For \f$x>0\f$,

\f[
 f(x)=\frac{1}{x\sigma\sqrt{2\pi}}
 \exp\!\left[-\frac{(\log x-\mu)^2}{2\sigma^2}\right],
 \qquad
 F(x)=\Phi\!\left(\frac{\log x-\mu}{\sigma}\right),
\f]

\f[
 S(x)=\Phi\!\left(\frac{\mu-\log x}{\sigma}\right),qquad
 Q(p)=\exp\!\left[\mu+\sigma\Phi^{-1}(p)\right].
\f]

Probability functions and both quantile tails delegate to the direct normal kernels after taking a
logarithm; sampling exponentiates a normal variate. There is no elementary characteristic function:

\f[
 \phi(t)=\frac{1}{\sigma\sqrt{2\pi}}
 \int_{-\infty}^{\infty}
 \exp\!\left(ite^y-\frac{(y-\mu)^2}{2\sigma^2}\right)dy.
\f]

For \f$\sigma^2<4\f$, RandLib evaluates the real and imaginary expectations with the ordinary
adaptive infinite-interval quadrature. For \f$\sigma^2\geq4\f$, it switches to the half-period
oscillatory integrator with Wynn epsilon acceleration.

Writing \f$v=e^{\sigma^2}\f$, the summaries are

\f[
 \operatorname{E}X=e^{\mu+\sigma^2/2},qquad
 \operatorname{Var}X=(v-1)e^{2\mu+\sigma^2},
\f]

\f[
 \operatorname{Median}X=e^\mu,qquad
 \operatorname{Mode}X=e^{\mu-\sigma^2},qquad
 \operatorname{Skew}X=(v+2)\sqrt{v-1},
\f]

\f[
 \operatorname{ExKurt}X=v^4+2v^3+3v^2-6.
\f]

`fit<LogNormalRand>` is the closed-form normal MLE applied to \f$\log x_i\f$: the fitted location is
the mean log and the fitted squared scale is the population variance of the logs. The same normal,
inverse-gamma, and normal--inverse-gamma conjugate models are available through `fit_bayes`, again
applied in log space.

### Exponentially modified Gaussian distribution (`ExponentiallyModifiedGaussianRand`)

`ExponentiallyModifiedGaussianRand(location, variance, rate)` represents the independent sum

\f[
 Z=X+Y,qquad X\sim N(\mu,\sigma^2),qquad
 Y\sim\operatorname{Exponential}(\lambda).
\f]

For all real \f$z\f$,

\f[
 f(z)=\frac\lambda2
 \exp\!\left[\frac\lambda2(2\mu+\lambda\sigma^2-2z)\right]
 \operatorname{erfc}\!\left(
 \frac{\mu+\lambda\sigma^2-z}{\sigma\sqrt2}\right).
\f]

With \f$w=(z-\mu)/\sigma\f$ and \f$v=\lambda\sigma\f$,

\f[
 F(z)=\Phi(w)-e^{-\lambda(z-\mu)+v^2/2}\Phi(w-v),
\f]

\f[
 S(z)=\Phi(-w)+e^{-\lambda(z-\mu)+v^2/2}\Phi(w-v).
\f]

The separately implemented survival expression avoids subtracting the complete CDF from one. The
characteristic function follows directly from independence:

\f[
 \phi(t)=\frac{\exp(i\mu t-\sigma^2t^2/2)}{1-it/\lambda}.
\f]

Sampling adds the library's normal and exponential variates. The mean, variance, skewness, and
excess kurtosis are

\f[
 \operatorname{E}Z=\mu+\frac1\lambda,qquad
 \operatorname{Var}Z=\sigma^2+\frac1{\lambda^2},
\f]

\f[
 \operatorname{Skew}Z=
 \frac{2\lambda^{-3}}{(\sigma^2+\lambda^{-2})^{3/2}},
 \qquad
 \operatorname{ExKurt}Z=
 \frac{6\lambda^{-4}}{(\sigma^2+\lambda^{-2})^2}.
\f]

No elementary median or mode is known. `Median()` uses the generic tail-aware TOMS 748 quantile
solver, and `Mode()` minimizes the negative log density with the shared Brent search, starting from
the mean. There is currently no `fit<ExponentiallyModifiedGaussianRand>` specialization.

## Nakagami family

### Nakagami distribution (`NakagamiRand`)

`NakagamiRand(shape, spread)` uses \f$\mu>0\f$ and \f$\omega>0\f$. Although the conventional
wireless-fading parameter range is often restricted to \f$\mu\geq1/2\f$, RandLib implements the
normalized law for every positive shape. Its defining transformation is

\f[
 X^2\sim\operatorname{Gamma}\!\left(\mu,\frac{\mu}{\omega}\right),
\f]

where the second Gamma parameter is a rate. Therefore, for \f$x\geq0\f$,

\f[
 f(x)=\frac{2\mu^\mu}{\Gamma(\mu)\omega^\mu}
 x^{2\mu-1}e^{-\mu x^2/\omega},
\f]

\f[
 F(x)=P\!\left(\mu,\frac{\mu x^2}{\omega}\right),
 \qquad
 S(x)=Q\!\left(\mu,\frac{\mu x^2}{\omega}\right).
\f]

The direct Gamma upper tail is used for \f$S\f$. Quantiles and samples are respectively the square
roots of Gamma quantiles and Gamma variates, so they inherit the Gamma family's tail inversion and
shape-dependent sampling dispatch.

All raw moments are summarized by

\f[
 \operatorname{E}X^r=
 \frac{\Gamma(\mu+r/2)}{\Gamma(\mu)}
 \left(\frac{\omega}{\mu}\right)^{r/2},
 \qquad r>-2\mu.
\f]

In particular,

\f[
 m=\operatorname{E}X=
 \frac{\Gamma(\mu+1/2)}{\Gamma(\mu)}
 \sqrt{\frac\omega\mu},qquad
 \operatorname{Var}X=\omega-m^2,
\f]

\f[
 \operatorname{Mode}X=
 \begin{cases}
 \sqrt{\omega(1-1/(2\mu))},&\mu>1/2,\\
 0,&0<\mu\leq1/2.
 \end{cases}
\f]

`Skewness()` and `ExcessKurtosis()` combine the same closed raw moments through order four. The
median is the square root of the corresponding Gamma median.

The characteristic function has the hypergeometric representation

\f[
 \begin{aligned}
 \phi(t)={}&{}_1F_1\!\left(\mu;\frac12;-\frac{\omega t^2}{4\mu}\right)\\
 &+it\sqrt{\frac\omega\mu}
 \frac{\Gamma(\mu+1/2)}{\Gamma(\mu)}
 {}_1F_1\!\left(\mu+\frac12;\frac32;-\frac{\omega t^2}{4\mu}\right).
 \end{aligned}
\f]

RandLib evaluates the defining real and imaginary expectations by adaptive infinite-interval
quadrature. When \f$\mu<1/2\f$, the density is singular at zero; the real part is then regularized as
\f$1+\int_0^\infty(\cos(tx)-1)f(x)\,dx\f$. This corrects the archived document's broader
\f$\mu<1\f$ threshold. There is currently no `fit<NakagamiRand>` specialization.

### Chi distribution (`ChiRand`)

`ChiRand(degree)` takes a positive integer \f$k\f$ and is

\f[
 X\sim\chi_k=\operatorname{Nakagami}\!\left(\frac{k}{2},k\right),
 \qquad X^2\sim\chi_k^2.
\f]

Thus

\f[
 f(x)=\frac{2^{1-k/2}}{\Gamma(k/2)}x^{k-1}e^{-x^2/2},
 \qquad
 F(x)=P\!\left(\frac{k}{2},\frac{x^2}{2}\right),
\f]

\f[
 \operatorname{E}X=\sqrt2\frac{\Gamma((k+1)/2)}{\Gamma(k/2)},qquad
 \operatorname{Var}X=k-(\operatorname{E}X)^2,qquad
 \operatorname{Mode}X=\sqrt{k-1}.
\f]

Its higher moments, characteristic function, direct upper tail, quantiles, and sampler are the
corresponding Nakagami operations above. There is currently no `fit<ChiRand>` specialization.

### Maxwell--Boltzmann distribution (`MaxwellBoltzmannRand`)

`MaxwellBoltzmannRand(scale)` takes \f$\sigma>0\f$ and represents

\f[
 \frac{X}{\sigma}\sim\chi_3,
 \qquad
 X\sim\operatorname{Nakagami}\!\left(\frac32,3\sigma^2\right).
\f]

The factor three in the Nakagami spread is required by the current parameterization and corrects
the missing factor in the archived document. RandLib specializes the density and tails as

\f[
 f(x)=\sqrt{\frac2\pi}\frac{x^2}{\sigma^3}e^{-x^2/(2\sigma^2)},
\f]

\f[
 F(x)=\operatorname{erf}\!\left(\frac{x}{\sigma\sqrt2}\right)
 -\sqrt{\frac2\pi}\frac{x}{\sigma}e^{-x^2/(2\sigma^2)},
\f]

\f[
 S(x)=\operatorname{erfc}\!\left(\frac{x}{\sigma\sqrt2}\right)
 +\sqrt{\frac2\pi}\frac{x}{\sigma}e^{-x^2/(2\sigma^2)}.
\f]

The summaries are

\f[
 \operatorname{E}X=2\sigma\sqrt{\frac2\pi},qquad
 \operatorname{Var}X=\sigma^2\left(3-\frac8\pi\right),qquad
 \operatorname{Mode}X=\sqrt2\,\sigma,
\f]

\f[
 \operatorname{Skew}X=
 \frac{2\sqrt2(16-5\pi)}{(3\pi-8)^{3/2}},
 \qquad
 \operatorname{ExKurt}X=
 \frac{4\left[\pi(40-3\pi)-96\right]}{(3\pi-8)^2}.
\f]

Sampling uses \f$\sigma\sqrt{2E+Z^2}\f$ with an independent unit exponential \f$E\f$ and standard
normal \f$Z\f$. Quantiles and the characteristic function use the general Nakagami paths. There is
currently no `fit<MaxwellBoltzmannRand>` specialization.

### Rayleigh distribution (`RayleighRand`)

`RayleighRand(scale)` takes \f$\sigma>0\f$ and is

\f[
 \frac{X}{\sigma}\sim\chi_2,
 \qquad
 X\sim\operatorname{Nakagami}(1,2\sigma^2).
\f]

Its direct probability and quantile formulas are

\f[
 f(x)=\frac{x}{\sigma^2}e^{-x^2/(2\sigma^2)},qquad
 F(x)=1-e^{-x^2/(2\sigma^2)},qquad
 S(x)=e^{-x^2/(2\sigma^2)},
\f]

\f[
 Q(p)=\sigma\sqrt{-2\log(1-p)},qquad
 Q(1-p)=\sigma\sqrt{-2\log p}.
\f]

The lower probability uses `expm1` and the lower quantile uses `log1p`. Sampling uses
\f$\sigma\sqrt{2E}\f$ with the shared unit-exponential ziggurat. The characteristic function is the
\f$\mu=1,\omega=2\sigma^2\f$ case of the Nakagami formula and is evaluated by its quadrature path.

The summaries are

\f[
 \operatorname{E}X=\sigma\sqrt{\frac\pi2},qquad
 \operatorname{Var}X=\frac{4-\pi}{2}\sigma^2,qquad
 \operatorname{Median}X=\sigma\sqrt{2\log2},qquad
 \operatorname{Mode}X=\sigma,
\f]

\f[
 \operatorname{Skew}X=
 \frac{2\sqrt\pi(\pi-3)}{(4-\pi)^{3/2}},
 \qquad
 \operatorname{ExKurt}X=
 -\frac{6\pi^2-24\pi+16}{(4-\pi)^2}.
\f]

`fit<RayleighRand>` supports MLE and the unbiased scale estimator. If
\f$R^2=(2n)^{-1}\sum_i x_i^2\f$, then

\f[
 \hat\sigma_{\mathrm{MLE}}=R,qquad
 \widetilde\sigma=
 \frac{\sqrt n\,\Gamma(n)}{\Gamma(n+1/2)}R.
\f]

The correction follows from
\f$\operatorname{E}R=\sigma\Gamma(n+1/2)/(\sqrt n\,\Gamma(n))\f$. RandLib evaluates the exact
factor in the log domain through \f$n=100\f$ and uses its log-gamma-ratio asymptotic series for
larger samples.

## Exponential-tail location families

### Asymmetric Laplace distributions (`AsymmetricLaplaceRand`, `CenteredAsymmetricLaplaceRand`)

`AsymmetricLaplaceRand(shift, scale, asymmetry)` uses
\f$m\in\mathbb R\f$, \f$\gamma>0\f$, and \f$\kappa>0\f$. With \f$z=x-m\f$,

\f[
 f(x)=\frac{\kappa}{\gamma(1+\kappa^2)}
 \begin{cases}
  e^{z/(\gamma\kappa)},&z<0,\\
  e^{-\kappa z/\gamma},&z\geq0,
 \end{cases}
\f]

\f[
 F(x)=
 \begin{cases}
  \dfrac{\kappa^2}{1+\kappa^2}e^{z/(\gamma\kappa)},&z<0,\\
  1-\dfrac{1}{1+\kappa^2}e^{-\kappa z/\gamma},&z\geq0,
 \end{cases}
\f]

and the direct survival implementation uses one minus the first branch below \f$m\f$ and the
exponential term in the second branch above it. The characteristic function is

\f[
 \phi(t)=\frac{e^{imt}}
 {(1-i\gamma t/\kappa)(1+i\gamma\kappa t)}.
\f]

This follows from the sampling representation

\f[
 X=m+\gamma\left(\frac{E_1}{\kappa}-\kappa E_2\right),
 \qquad E_1,E_2\stackrel{\mathrm{iid}}{\sim}\operatorname{Exponential}(1),
\f]

which is also RandLib's sampler. The summaries are

\f[
 \operatorname{E}X=m+\gamma(\kappa^{-1}-\kappa),
 \qquad
 \operatorname{Var}X=\gamma^2(\kappa^{-2}+\kappa^2),
 \qquad
 \operatorname{Mode}X=m,
\f]

\f[
 \operatorname{Skew}X=
 \frac{2(1-\kappa^6)}{(1+\kappa^4)^{3/2}},
 \qquad
 \operatorname{ExKurt}X=
 6-\frac{12}{(\kappa^2+\kappa^{-2})^2},
\f]

\f[
 h(X)=1+\log\!\left[\gamma(\kappa+\kappa^{-1})\right].
\f]

The last expression is the corrected B41 implementation; the inherited 1.x code had used an
unrelated `log1p` expression. The median is obtained without root finding:

\f[
 \operatorname{Median}X=
 \begin{cases}
  m+\gamma\kappa\log\!\left[(1+\kappa^{-2})/2\right],&\kappa>1,\\
  m-(\gamma/\kappa)\log\!\left[(1+\kappa^2)/2\right],&0<\kappa\leq1.
 \end{cases}
\f]

All other quantiles use the corresponding logarithm of the appropriate exponential branch.

`CenteredAsymmetricLaplaceRand(scale, asymmetry)` fixes \f$m=0\f$. “Centered” therefore refers to
the density's join and mode, not its mean unless \f$\kappa=1\f$. It exposes the exponential-family
sufficient statistic
\f$(\max(X,0),\max(-X,0))\f$ with natural parameters
\f$(-\kappa/\gamma,-1/(\gamma\kappa))\f$.

There is no free `fit` specialization for either asymmetric class. The deprecated member API can
fit the scale from the weighted absolute deviations, fit the shift as the sample quantile whose
level is \f$\kappa^2/(1+\kappa^2)\f$, and fit asymmetry from the resulting one-sided sums.

### Laplace distribution (`LaplaceRand`)

`LaplaceRand(shift, scale)` is the symmetric case \f$\kappa=1\f$:

\f[
 f(x)=\frac1{2\gamma}e^{-|x-m|/\gamma},
 \qquad
 \phi(t)=\frac{e^{imt}}{1+\gamma^2t^2}.
\f]

It uses the asymmetric family's direct piecewise CDF, survival function, logarithmic quantiles, and
difference-of-exponentials sampler; the single-variate fast path equivalently draws an exponential
magnitude and an independent sign bit. Its summaries reduce to

\f[
 \operatorname{E}X=\operatorname{Median}X=\operatorname{Mode}X=m,
 \qquad
 \operatorname{Var}X=2\gamma^2,
\f]

\f[
 \operatorname{Skew}X=0,qquad
 \operatorname{ExKurt}X=3,qquad
 h(X)=1+\log(2\gamma).
\f]

`fit<LaplaceRand>` performs the closed-form MLE: a sample median estimates \f$m\f$, and the mean
absolute deviation from that median estimates \f$\gamma\f$.

### Logistic distribution (`LogisticRand`)

`LogisticRand(location, scale)` uses \f$\mu\in\mathbb R\f$ and \f$s>0\f$. With
\f$z=(x-\mu)/s\f$,

\f[
 f(x)=\frac{e^{-z}}{s(1+e^{-z})^2}
 =\frac{e^{-|z|}}{s(1+e^{-|z|})^2},
\f]

\f[
 F(x)=\frac1{1+e^{-z}},qquad
 S(x)=\frac1{1+e^z},qquad
 Q(p)=\mu+s\log\frac{p}{1-p}.
\f]

The density uses the symmetric second form to avoid overflow, `logf` uses a softplus kernel, and
both quantile tails use `log1p`. Sampling applies the logit transform to an open-interval uniform
variate. The characteristic function and summaries are

\f[
 \phi(t)=e^{i\mu t}\frac{\pi st}{\sinh(\pi st)},
\f]

\f[
 \operatorname{E}X=\operatorname{Median}X=\operatorname{Mode}X=\mu,
 \qquad
 \operatorname{Var}X=\frac{\pi^2s^2}{3},
\f]

\f[
 \operatorname{Skew}X=0,qquad
 \operatorname{ExKurt}X=\frac65,qquad
 h(X)=2+\log s.
\f]

`fit<LogisticRand>` initializes location with the sample mean and scale with
\f$\sqrt{3s_n^2}/\pi\f$, then runs the common constrained MLE optimizer. With scale fixed, the
deprecated `FitLocation` member instead solves the location score by safeguarded Newton iteration.

### Hyperbolic secant distribution (`SechRand`)

`SechRand` is parameter-free and standardized to unit variance:

\f[
 f(x)=\frac1{2\cosh(\pi x/2)},qquad
 F(x)=\frac2\pi\arctan(e^{\pi x/2}),qquad
 S(x)=\frac2\pi\arctan(e^{-\pi x/2}),
\f]

\f[
 Q(p)=\frac2\pi\log\tan\frac{\pi p}{2},qquad
 \phi(t)=\frac1{\cosh t}.
\f]

The log density uses softplus rather than taking the logarithm of a potentially overflowing
hyperbolic cosine. Sampling uses
\f$(2/\pi)\log|C|\f$ for a standard Cauchy variate \f$C\f$. Finally,

\f[
 \operatorname{E}X=\operatorname{Median}X=\operatorname{Mode}X=0,qquad
 \operatorname{Var}X=1,qquad
 \operatorname{Skew}X=0,qquad
 \operatorname{ExKurt}X=2,
\f]

\f[
 h(X)=\frac{4G}{\pi},
\f]

where \f$G\f$ is Catalan's constant. There are no parameters to fit.

## Bounded algebraic distributions

### Triangular distribution (`TriangularRand`)

`TriangularRand(lower, mode, upper)` uses \f$a<b\f$ and \f$a\leq c\leq b\f$. Its density is

\f[
 f(x)=
 \begin{cases}
  \dfrac{2(x-a)}{(b-a)(c-a)},&a\leq x<c,\\[4pt]
  \dfrac{2}{b-a},&x=c,\\[4pt]
  \dfrac{2(b-x)}{(b-a)(b-c)},&c<x\leq b,\\[4pt]
  0,&\text{otherwise},
 \end{cases}
\f]

with the natural one-sided interpretation when \f$c=a\f$ or \f$c=b\f$. Its distribution and
survival functions are evaluated directly as

\f[
 F(x)=
 \begin{cases}
  0,&x<a,\\
  \dfrac{(x-a)^2}{(b-a)(c-a)},&a\leq x\leq c,\\[4pt]
  1-\dfrac{(b-x)^2}{(b-a)(b-c)},&c<x<b,\\[4pt]
  1,&x\geq b,
 \end{cases}
\f]

and \f$S\f$ uses the complementary square in each branch rather than subtracting a computed tail
from one. Inverting these branches gives both the quantile and sampling algorithm: for a uniform
\f$U\f$, use

\f[
 X=
 \begin{cases}
  a+\sqrt{U(b-a)(c-a)},&U<(c-a)/(b-a),\\
  b-\sqrt{(1-U)(b-a)(b-c)},&\text{otherwise}.
 \end{cases}
\f]

For an interior mode the characteristic function is

\f[
 \phi(t)=
 \frac{2\left[(b-a)e^{itc}-(b-c)e^{ita}-(c-a)e^{itb}\right]}
 {(b-a)(c-a)(b-c)t^2}.
\f]

RandLib uses separate limiting expressions for a mode at either endpoint and returns one directly
when cancellation makes the small-\f$t\f$ numerator indistinguishable from zero.

The principal summaries are

\f[
 \operatorname{E}X=\frac{a+b+c}{3},\qquad
 \operatorname{Var}X=\frac{a^2+b^2+c^2-ab-ac-bc}{18},
\f]

\f[
 \operatorname{Median}X=
 \begin{cases}
  a+\sqrt{(b-a)(c-a)/2},&2c>a+b,\\
  b-\sqrt{(b-a)(b-c)/2},&2c\leq a+b,
 \end{cases}
 \qquad \operatorname{Mode}X=c,
\f]

\f[
 \operatorname{Skew}X=
 \frac{\sqrt2(a+b-2c)(2a-b-c)(a-2b+c)}
 {5\left(a^2+b^2+c^2-ab-ac-bc\right)^{3/2}},
 \qquad \operatorname{ExKurt}X=-\frac35.
\f]

There is currently no `fit<TriangularRand>` specialization.

### Raised-cosine distributions (`RaisedCosineRand`, `RaabGreenRand`)

`RaisedCosineRand(location, scale)` takes \f$\mu\in\mathbb R\f$ and \f$s>0\f$, with support
\f$[\mu-s,\mu+s]\f$. Writing \f$y=(x-\mu)/s\f$,

\f[
 f(x)=\frac{1+\cos(\pi y)}{2s},\qquad
 F(x)=\frac12\left(1+y+\frac{\sin(\pi y)}{\pi}\right),
 \qquad |y|<1,
\f]

and the direct survival formula changes the two plus signs after one to minus signs. The
characteristic function is

\f[
 \phi(t)=e^{i\mu t}\frac{\pi^2\sin(st)}{st\,[\pi^2-(st)^2]}.
\f]

At \f$st=\pi\f$ the displayed expression is indeterminate. The implementation rewrites it with a
cardinal-sine factor, preserving the analytic limit instead of evaluating \f$0/0\f$ (B33).

The mean, median, and mode all equal \f$\mu\f$, and

\f[
 \operatorname{Var}X=s^2\left(\frac13-\frac{2}{\pi^2}\right),\qquad
 \operatorname{Skew}X=0,
\f]

\f[
 \operatorname{ExKurt}X=
 \frac{6(90-\pi^4)}{5(\pi^2-6)^2}.
\f]

The CDF has no elementary inverse, so quantiles use the shared bracketed TOMS 748 solver. Sampling
uses Devroye's alternating-series accept/reject method rather than numerically inverting the CDF.
`RaabGreenRand` is the fixed \f$\mu=0,s=\pi\f$ specialization. Neither class currently has a
free `fit` specialization.

### Wigner semicircle distribution (`WignerSemicircleRand`)

`WignerSemicircleRand(radius)` takes \f$R>0\f$ and has support \f$[-R,R]\f$:

\f[
 f(x)=\frac{2}{\pi R^2}\sqrt{R^2-x^2},
\f]

\f[
 F(x)=\frac12+\frac1\pi\left[
 \frac{x\sqrt{R^2-x^2}}{R^2}+\arcsin\frac{x}{R}
 \right],\qquad -R\leq x\leq R.
\f]

The implementation does not evaluate that lower tail by adding nearly cancelling terms. For
\f$x\leq0\f$ it sets

\f[
 \delta=1+x/R,\qquad \varphi=2\arcsin\sqrt{\delta/2},\qquad
 F(x)=\frac{2\varphi-\sin(2\varphi)}{2\pi},
\f]

using a local series for \f$z-\sin z\f$ at small \f$z\f$; positive arguments and survival
probabilities follow by symmetry. Quantiles use TOMS 748. Sampling exploits

\f[
 Y\sim\operatorname{Beta}(3/2,3/2),\qquad X=R(2Y-1).
\f]

The closed-form characteristic function and summaries are

\f[
 \phi(t)=\frac{2J_1(Rt)}{Rt},
\f]

\f[
 \operatorname{E}X=\operatorname{Median}X=\operatorname{Mode}X=0,
 \qquad \operatorname{Var}X=\frac{R^2}{4},
\f]

\f[
 \operatorname{Skew}X=0,\qquad
 \operatorname{ExKurt}X=-1,\qquad
 h(X)=\log(\pi R)-\frac12.
\f]

The class currently obtains \f$\phi\f$ through the common finite-support expectation quadrature
rather than a dedicated Bessel-function path. There is no `fit<WignerSemicircleRand>`
specialization.

### Irwin--Hall distribution (`IrwinHallRand`)

`IrwinHallRand(number)` is the sum of \f$n\geq1\f$ independent standard uniform variables. On
\f$0\leq x\leq n\f$,

\f[
 f(x)=\frac1{(n-1)!}\sum_{k=0}^{\lfloor x\rfloor}
 (-1)^k{n\choose k}(x-k)^{n-1},
\f]

\f[
 F(x)=\frac1{n!}\sum_{k=0}^{\lfloor x\rfloor}
 (-1)^k{n\choose k}(x-k)^n.
\f]

RandLib has direct branches for \f$n=1\f$ and \f$n=2\f$; the general finite sums use a cached
log-factorial and recursively updated binomial coefficients. The characteristic function is

\f[
 \phi(t)=\left(\frac{e^{it}-1}{it}\right)^n
 =e^{int/2}\left[\frac{\sin(t/2)}{t/2}\right]^n.
\f]

Sampling explicitly sums \f$n\f$ independent uniforms. Other summaries are

\f[
 \operatorname{E}X=\operatorname{Median}X=\frac n2,\qquad
 \operatorname{Var}X=\frac n{12},\qquad
 \operatorname{Skew}X=0,\qquad
 \operatorname{ExKurt}X=-\frac6{5n}.
\f]

The API returns \f$n/2\f$ as a representative mode; for \f$n=1\f$ the density is constant and the
mode is not unique. Quantiles use TOMS 748, and there is no `fit<IrwinHallRand>` specialization.

### Degenerate distribution (`DegenerateRand`)

`DegenerateRand(value)` represents the point mass \f$P(X=a)=1\f$. It is retained in the continuous
class hierarchy for API composition, but it is singular with respect to Lebesgue measure. Its
`f(a)` and `logf(a)` return positive infinity as a sentinel, while both return zero or negative
infinity, respectively, away from \f$a\f$. Its actual distribution functions are

\f[
 F(x)=\mathbf1_{[a,\infty)}(x),\qquad
 S(x)=\mathbf1_{(-\infty,a)}(x),\qquad
 Q(p)=a,
\f]

and

\f[
 \phi(t)=e^{iat},\qquad
 \operatorname{E}X=\operatorname{Median}X=\operatorname{Mode}X=a,
 \qquad \operatorname{Var}X=0,
\f]

\f[
 h(X)=0.
\f]

Skewness and excess kurtosis are undefined and return NaN. Sampling always returns \f$a\f$.
`fit<DegenerateRand>` supports MLE for a nonempty sample only when all observations are exactly
equal; the common observation becomes the fitted point mass.
