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
