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
