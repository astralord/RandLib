# Discrete distributions {#discrete_distribution_reference}

The discrete reference covers `BernoulliRand`, `BetaBinomialRand`, `BinomialRand`, `CategoricalRand`,
`GeometricRand`, `HyperGeometricRand`, `LogarithmicRand`, `NegativeBinomialRand`,
`NegativeHyperGeometricRand`, `PoissonRand`, `RademacherRand`, `SkellamRand`, `UniformDiscreteRand`,
`YuleRand`, `ZetaRand` and `ZipfRand`.

For every class, `F(k)` means \f$P(X\leq k)\f$ and `S(k)` means \f$P(X>k)\f$. Classes with a
dedicated `S` evaluate that upper tail directly; otherwise the shared implementation returns
\f$1-F(k)\f$. Discrete quantiles are the smallest integers whose CDF reaches the requested
probability. The generic characteristic-function path sums outward from the median and uses an
Euler-transformed remainder for unbounded slowly decaying tails.

## Binomial family

### Bernoulli distribution (`BernoulliRand`)

`BernoulliRand(probability)` accepts \f$0\leq p\leq1\f$ and has support \f$\{0,1\}\f$:

\f[
 P(X=k)=
 \begin{cases}
  1-p,&k=0,\\
  p,&k=1,\\
  0,&\text{otherwise}.
 \end{cases}
\f]

Consequently,

\f[
 F(k)=
 \begin{cases}
  0,&k<0,\\
  1-p,&0\leq k<1,\\
  1,&k\geq1,
 \end{cases}
 \qquad
 S(k)=
 \begin{cases}
  1,&k<0,\\
  p,&0\leq k<1,\\
  0,&k\geq1.
 \end{cases}
\f]

Its probability-generating and characteristic functions are

\f[
 G(z)=1-p+pz,\qquad
 \phi(t)=1-p+pe^{it}.
\f]

The summaries are

\f[
 \operatorname{E}X=p,\qquad
 \operatorname{Var}X=p(1-p),
\f]

\f[
 \operatorname{Skew}X=\frac{1-2p}{\sqrt{p(1-p)}},
 \qquad
 \operatorname{ExKurt}X=\frac{1}{p(1-p)}-6,
\f]

for non-degenerate \f$p\f$. The API median is zero when \f$p\leq1/2\f$ and one otherwise. Its
single reported mode is zero for \f$p<1/2\f$ and one for \f$p\geq1/2\f$; mathematically both
points are modes at equality. The entropy is

\f[
 H(X)=-p\log p-(1-p)\log(1-p),
\f]

with \f$0\log0=0\f$, so the deterministic endpoints return zero.

Sampling compares one uniform variate with \f$p\f$. Bulk sampling at \f$p=1/2\f$ instead consumes
one unbiased random bit per result. `fit<BernoulliRand>` supports closed-form MLE and
method-of-moments fitting, both with \f$\widehat p=\bar x\f$.

### Binomial distribution (`BinomialRand`)

`BinomialRand(number, probability)` uses integer \f$n>0\f$ and \f$0\leq p\leq1\f$. For
\f$k=0,\ldots,n\f$,

\f[
 P(X=k)={n\choose k}p^k(1-p)^{n-k}.
\f]

RandLib evaluates interior log probabilities with a saddle-point expression based on Stirling
errors and the deviance term `bd0`; this avoids cancellation in the direct logarithm of the
binomial coefficient for large \f$n\f$. The two tails are

\f[
 F(k)=I_{1-p}(n-k,k+1),
 \qquad
 S(k)=I_p(k+1,n-k),
\f]

for \f$0\leq k<n\f$, with exact support-end values outside that range. Both use direct regularized
incomplete-beta evaluations.

The generating functions are

\f[
 G(z)=(1-p+pz)^n,\qquad
 \phi(t)=(1-p+pe^{it})^n.
\f]

Writing \f$q=1-p\f$,

\f[
 \operatorname{E}X=np,\qquad
 \operatorname{Var}X=npq,
\f]

\f[
 \operatorname{Skew}X=\frac{q-p}{\sqrt{npq}},
 \qquad
 \operatorname{ExKurt}X=\frac{1/(pq)-6}{n}.
\f]

The standardized moments apply for \f$0<p<1\f$.
The mode returned by the API is
\f$\min\{n,\lfloor(n+1)p\rfloor\}\f$; when \f$(n+1)p\f$ is an integer in the interior,
the preceding integer is also a mathematical mode. The median is selected from
\f$\lfloor np\rfloor\f$ and \f$\lceil np\rceil\f$ by evaluating the direct CDF.

The cached sampler chooses among three algorithms:

- direct Bernoulli summation for tiny \f$n\f$ or moderate \f$n\f$ close to \f$p=1/2\f$;
- a geometric waiting-time construction when \f$n\min(p,q)\f$ is small;
- the Devroye--Naderlsamanl rejection method otherwise.

The rejection branch first samples at the nearby lattice probability
\f$p_0=\lfloor n\min(p,q)\rfloor/n\f$, adds an independent residual binomial by waiting times when
needed, and reflects the result for \f$p>1/2\f$. Exact \f$p=0\f$ and \f$p=1\f$ bulk samples are
filled without consuming the engine.

`fit<BinomialRand>` currently requires \f$n\f$ as a fixed parameter and returns the closed-form
MLE

\f[
 \widehat p=\frac{1}{mn}\sum_{j=1}^m x_j.
\f]

With a \f$\operatorname{Beta}(a,b)\f$ prior, the conditional conjugate helper returns

\f[
 p\mid x_1,\ldots,x_m\sim
 \operatorname{Beta}\!\left(a+\sum_jx_j,\ b+mn-\sum_jx_j\right).
\f]

### Beta-binomial distribution (`BetaBinomialRand`)

`BetaBinomialRand(number, shape1, shape2)` uses \f$n\in\mathbb N\f$ and
\f$\alpha,\beta>0\f$. It is the beta mixture

\f[
 P\sim\operatorname{Beta}(\alpha,\beta),
 \qquad
 X\mid P\sim\operatorname{Binomial}(n,P),
\f]

and therefore

\f[
 P(X=k)={n\choose k}
 \frac{B(k+\alpha,n-k+\beta)}{B(\alpha,\beta)},
 \qquad k=0,\ldots,n.
\f]

The log PMF is assembled from `lgamma` and log-factorials. The current CDF explicitly sums the
PMF from zero through \f$k\f$; survival and quantiles use the shared discrete paths. The
probability-generating function is

\f[
 G(z)={}_2F_1(-n,\alpha;\alpha+\beta;1-z),
 \qquad
 \phi(t)=G(e^{it}).
\f]

The first two moments are

\f[
 \operatorname{E}X=\frac{n\alpha}{\alpha+\beta},
\f]

\f[
 \operatorname{Var}X=
 \frac{n\alpha\beta(\alpha+\beta+n)}
 {(\alpha+\beta)^2(\alpha+\beta+1)}.
\f]

Its skewness is

\f[
 \operatorname{Skew}X=
 \frac{(\alpha+\beta+2n)(\beta-\alpha)}
 {\alpha+\beta+2}
 \sqrt{\frac{\alpha+\beta+1}
 {n\alpha\beta(\alpha+\beta+n)}}.
\f]

For completeness, with \f$s=\alpha+\beta\f$, the implemented excess kurtosis is

\f[
 \operatorname{ExKurt}X=
 \frac{s^2(s+1)}
 {n\alpha\beta(s+2)(s+3)(s+n)}
 \left\{
 s(s-1+6n)+3\left[
 \alpha\beta(n-2)+2n^2
 -\frac{\alpha\beta n(6-n)}s
 -\frac{6\alpha\beta n^2}{s^2}
 \right]\right\}-3.
\f]

The mode
starts from the beta mixing mode scaled by \f$n\f$ and compares the adjacent log probabilities;
the median and characteristic function use the generic finite-support algorithms.

Sampling draws \f$P\f$ from the beta law and then a conditional binomial variate. When
\f$\alpha=\beta=1\f$, the beta-binomial law is exactly uniform on
\f$\{0,\ldots,n\}\f$, so the implementation samples that integer directly. There is currently no
`fit<BetaBinomialRand>` specialization or dedicated entropy member.

## Finite categorical laws

### Categorical distribution (`CategoricalRand`)

`CategoricalRand(probabilities)` stores a nonempty vector
\f$(p_0,\ldots,p_{K-1})\f$ with nonnegative entries summing to one. Its support is
\f$\{0,\ldots,K-1\}\f$ and

\f[
 P(X=k)=p_k,\qquad
 F(k)=\sum_{j=0}^{k}p_j.
\f]

The implementation sums the shorter index side of the CDF: forward from zero when \f$k\f$ is in
the lower half, or backward as one minus the upper probabilities otherwise. Its generating
functions are finite polynomials,

\f[
 G(z)=\sum_{k=0}^{K-1}p_kz^k,
 \qquad
 \phi(t)=\sum_{k=0}^{K-1}p_ke^{itk}.
\f]

The first two moments are

\f[
 \operatorname{E}X=\sum_{k=0}^{K-1}kp_k,
 \qquad
 \operatorname{Var}X=
 \sum_{k=0}^{K-1}k^2p_k-\left(\sum_{k=0}^{K-1}kp_k\right)^2.
\f]

Higher standardized moments use the shared exact finite-support expectation sums. `Mode()`
returns the first index attaining \f$\max_kp_k\f$, and the median is the first index whose
cumulative probability reaches one half. The mathematical entropy is

\f[
 H(X)=-\sum_{k=0}^{K-1}p_k\log p_k,
 \qquad 0\log0=0,
\f]

although the class does not expose a dedicated entropy member.

Sampling is inverse transform: one uniform draw is accumulated through the stored probabilities
until the corresponding category is reached. This is linear in \f$K\f$ and does not construct an
alias table. There is currently no `fit<CategoricalRand>` specialization.

### Discrete-uniform distribution (`UniformDiscreteRand`)

`UniformDiscreteRand(minimum, maximum)` uses integer bounds \f$a\leq b\f$. Set
\f$n=b-a+1\f$. On \f$\{a,\ldots,b\}\f$,

\f[
 P(X=k)=\frac1n,
 \qquad
 F(k)=\frac{k-a+1}{n}\quad(a\leq k\leq b).
\f]

The probability-generating function and characteristic function are

\f[
 G(z)=\frac{z^a(1-z^n)}{n(1-z)},
\f]

\f[
 \phi(t)=
 e^{i(a+b)t/2}\frac{\sin(nt/2)}{n\sin(t/2)}.
\f]

At integer multiples of \f$2\pi\f$, the continuous limit is exactly one. RandLib first reduces the
angle modulo \f$2\pi\f$, returns that exact limit at a lattice period, and otherwise evaluates the
centered sinc form. This avoids both the former \f$0/0\f$ and precision loss from subtracting nearby
trigonometric values.

For \f$n>1\f$,

\f[
 \operatorname{E}X=\frac{a+b}{2},
 \qquad
 \operatorname{Var}X=\frac{n^2-1}{12},
\f]

\f[
 \operatorname{Skew}X=0,
 \qquad
 \operatorname{ExKurt}X=
 -\frac{6(n^2+1)}{5(n^2-1)}.
\f]

The entropy is \f$\log n\f$. Every support point is a mathematical mode; the API reports an
integer midpoint, and reports the lower integer midpoint as its median when two central points
exist. The \f$n=1\f$ law is degenerate, with zero variance and undefined standardized fourth
moment.

Sampling uses rejection-based `uniformBelow`, so reducing raw engine output modulo \f$n\f$ cannot
introduce bias. `fit<UniformDiscreteRand>` provides the closed-form MLE

\f[
 \widehat a=\min_jx_j,\qquad
 \widehat b=\max_jx_j.
\f]

### Rademacher distribution (`RademacherRand`)

The parameter-free Rademacher law has

\f[
 P(X=-1)=P(X=1)=\frac12,
\f]

\f[
 F(k)=
 \begin{cases}
  0,&k<-1,\\
  1/2,&-1\leq k<1,\\
  1,&k\geq1.
 \end{cases}
\f]

Its generating functions and summaries are

\f[
 G(z)=\frac{z+z^{-1}}2,
 \qquad
 \phi(t)=\cos t,
\f]

\f[
 \operatorname{E}X=0,\qquad
 \operatorname{Var}X=1,\qquad
 \operatorname{Skew}X=0,\qquad
 \operatorname{ExKurt}X=-2.
\f]

Both support points are modes and medians under common set-valued definitions; RandLib consistently
returns the lower one, \f$-1\f$. The entropy is \f$\log2\f$. Sampling maps one unbiased engine bit
to \f$-1\f$ or \f$1\f$. There are no parameters to fit.

## Poisson family

### Poisson distribution (`PoissonRand`)

`PoissonRand(rate)` uses \f$\lambda>0\f$ and support
\f$\{0,1,\ldots\}\f$:

\f[
 P(X=k)=e^{-\lambda}\frac{\lambda^k}{k!}.
\f]

The interior log PMF uses Loader's saddle-point form,

\f[
 \log P(X=k)=
 -\operatorname{stirlerr}(k)-\operatorname{bd0}(k,\lambda)
 -\frac12\log(2\pi k),
\f]

instead of subtracting three large terms in
\f$k\log\lambda-\lambda-\log(k!)\f$. The direct tails are regularized incomplete-gamma ratios:

\f[
 F(k)=Q(k+1,\lambda),
 \qquad
 S(k)=P(k+1,\lambda),
 \qquad k\geq0.
\f]

Here \f$P(a,x)=\gamma(a,x)/\Gamma(a)\f$ and
\f$Q(a,x)=\Gamma(a,x)/\Gamma(a)\f$.

The probability-generating and characteristic functions are

\f[
 G(z)=e^{\lambda(z-1)},
 \qquad
 \phi(t)=e^{\lambda(e^{it}-1)}.
\f]

All cumulants equal \f$\lambda\f$, giving

\f[
 \operatorname{E}X=\operatorname{Var}X=\lambda,
 \qquad
 \operatorname{Skew}X=\lambda^{-1/2},
 \qquad
 \operatorname{ExKurt}X=\lambda^{-1}.
\f]

The mode is \f$\lfloor\lambda\rfloor\f$; for positive integer \f$\lambda\f$,
\f$\lambda-1\f$ is also a mode and the API reports \f$\lambda\f$. `Median()` returns the
standard accurate approximation

\f[
 \max\!\left(0,\left\lfloor\lambda+\frac13-\frac{0.02}{\lambda}\right\rfloor\right)
\f]

rather than running a probability search.

Sampling uses inversion around \f$\lfloor\lambda\rfloor\f$ for \f$\lambda<10\f$. For larger rates
it uses a multi-region normal/exponential rejection envelope and accepts against the exact Poisson
log-PMF ratio. Constructed objects cache the envelope constants; the static convenience function
builds the same constants for its supplied rate and engine, so both entry points sample the exact
law.

`fit<PoissonRand>` uses the shared closed-form MLE and method-of-moments estimate

\f[
 \widehat\lambda=\bar x.
\f]

For a shape--rate prior \f$\lambda\sim\operatorname{Gamma}(a,b)\f$, `fit_bayes` returns

\f[
 \lambda\mid x_1,\ldots,x_m\sim
 \operatorname{Gamma}\!\left(a+\sum_jx_j,\ b+m\right).
\f]

The class exposes exponential-family adjusted cross-entropy terms for KL calculations, but no
dedicated Shannon-entropy member.

### Skellam distribution (`SkellamRand`)

If \f$Y\sim\operatorname{Poisson}(\mu_1)\f$ and
\f$Z\sim\operatorname{Poisson}(\mu_2)\f$ are independent, then
\f$X=Y-Z\f$ has the Skellam law. RandLib requires \f$\mu_1,\mu_2>0\f$ and supports every integer.
Its PMF is

\f[
 P(X=k)=
 e^{-(\mu_1+\mu_2)}
 \left(\frac{\mu_1}{\mu_2}\right)^{k/2}
 I_{|k|}\!\left(2\sqrt{\mu_1\mu_2}\right).
\f]

`logP` evaluates this with the scaled log-Bessel kernel. With
\f$Q_M(a,b)\f$ the generalized Marcum Q function and
\f$P_M(a,b)=1-Q_M(a,b)\f$, the two direct tails are

\f[
 F(k)=
 \begin{cases}
 P_{-k}(\sqrt{2\mu_1},\sqrt{2\mu_2}),&k<0,\\
 Q_{k+1}(\sqrt{2\mu_2},\sqrt{2\mu_1}),&k\geq0,
 \end{cases}
\f]

\f[
 S(k)=
 \begin{cases}
 Q_{-k}(\sqrt{2\mu_1},\sqrt{2\mu_2}),&k<0,\\
 P_{k+1}(\sqrt{2\mu_2},\sqrt{2\mu_1}),&k\geq0.
 \end{cases}
\f]

Thus neither tail is obtained by subtracting the other from one.

The generating functions are

\f[
 G(z)=
 \exp\!\left[\mu_1(z-1)+\mu_2(z^{-1}-1)\right],
\f]

\f[
 \phi(t)=
 \exp\!\left[
  (\mu_1+\mu_2)(\cos t-1)
  +i(\mu_1-\mu_2)\sin t
 \right].
\f]

The principal summaries are

\f[
 \operatorname{E}X=\mu_1-\mu_2,
 \qquad
 \operatorname{Var}X=\mu_1+\mu_2,
\f]

\f[
 \operatorname{Skew}X=
 \frac{\mu_1-\mu_2}{(\mu_1+\mu_2)^{3/2}},
 \qquad
 \operatorname{ExKurt}X=\frac1{\mu_1+\mu_2}.
\f]

The median is found by discrete probability search initialized at the mean. The mode search starts
at that median and walks through adjacent log probabilities until it reaches a maximum; merely
truncating \f$\mu_1-\mu_2\f$ is not generally correct.

Sampling subtracts independent draws from two cached `PoissonRand` objects. There is currently no
`fit<SkellamRand>` specialization or dedicated entropy member.
