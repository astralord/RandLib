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

## Geometric and logarithmic-series families

### Geometric distribution (`GeometricRand`)

`GeometricRand(p)` counts failures before the first success. With \f$0<p\leq1\f$, \f$q=1-p\f$,
and \f$k\in\{0,1,\ldots\}\f$,

\f[
 P(X=k)=p q^k,
 \qquad
 F(k)=1-q^{k+1},
 \qquad
 S(k)=q^{k+1}.
\f]

The implementation evaluates \f$F\f$ with `expm1` and evaluates \f$S\f$ directly, so neither short
tail is formed by subtracting a rounded probability from one. The probability-generating and
characteristic functions are

\f[
 G(z)=\frac{p}{1-qz},
 \qquad
 \phi(t)=\frac{p}{1-qe^{it}}.
\f]

Its principal summaries are

\f[
 \operatorname{E}X=\frac qp,
 \qquad
 \operatorname{Var}X=\frac q{p^2},
 \qquad
 \operatorname{Skew}X=\frac{1+q}{\sqrt q},
 \qquad
 \operatorname{ExKurt}X=6+\frac{p^2}{q},
\f]

with the standardized moments undefined at the deterministic endpoint \f$p=1\f$. The mode is zero
and the lower median is

\f[
 \left\lceil\frac{-\log2}{\log q}\right\rceil-1,
\f]

interpreted as zero when \f$p=1\f$. In natural units its Shannon entropy is

\f[
 H(X)=-\log p-\frac qp\log q
     =\frac{-p\log p-q\log q}{p},
\f]

with the usual \f$0\log0=0\f$ limit.

Sampling uses the identity

\f[
 X=\left\lfloor\frac{E}{-\log q}\right\rfloor,
 \qquad E\sim\operatorname{Exp}(1),
\f]

for small \f$p\f$. Otherwise a constructed distribution inverts a cached 16-entry CDF table,
recursively shifting by 16 when a uniform falls in the remaining tail. The static convenience
sampler uses the exponential transform below a slightly lower crossover and ordinary sequential
CDF inversion above it. All three paths sample the exact law.

Both `fit<GeometricRand>` and `fit_moments<GeometricRand>` use the closed-form estimate

\f[
 \widehat p=\frac1{1+\bar x}.
\f]

The inherited legacy beta-prior member gives
\f$\operatorname{Beta}(a+m,b+\sum_jx_j)\f$ for \f$m\f$ observations, but the unified
`fit_bayes<GeometricRand>` specialization is not currently provided.

### Negative-binomial distribution (`NegativeBinomialRand`, `PascalRand`, `PolyaRand`)

RandLib counts failures before the \f$r\f$th success. `PascalRand<Int>` is the integer-shape alias;
`PolyaRand<Int>` admits any real \f$r>0\f$. For \f$0<p\leq1\f$, \f$q=1-p\f$, and
\f$k\in\{0,1,\ldots\}\f$,

\f[
 P(X=k)=
 \frac{\Gamma(k+r)}{\Gamma(r)\,\Gamma(k+1)}p^r q^k.
\f]

The two tails are evaluated as regularized incomplete-beta functions,

\f[
 F(k)=I_p(r,k+1),
 \qquad
 S(k)=I_q(k+1,r),
\f]

and the log PMF uses `lgamma` and `lfact`. Its generating functions are

\f[
 G(z)=\left(\frac{p}{1-qz}\right)^r,
 \qquad
 \phi(t)=\left(\frac{p}{1-qe^{it}}\right)^r.
\f]

The principal summaries are

\f[
 \operatorname{E}X=\frac{rq}{p},
 \qquad
 \operatorname{Var}X=\frac{rq}{p^2},
\f]

\f[
 \operatorname{Skew}X=\frac{1+q}{\sqrt{rq}},
 \qquad
 \operatorname{ExKurt}X=\frac6r+\frac{p^2}{rq}.
\f]

The reported mode is

\f[
 \max\!\left(0,\left\lfloor\frac{(r-1)q}{p}\right\rfloor\right),
\f]

and the median is found by the generic discrete quantile search.

For integer \f$r\f$, sampling sums \f$r\f$ geometric variates, using the table or exponential
method described above. For noninteger \f$r\f$ and the remaining parameter region it uses the exact
gamma--Poisson mixture

\f[
 \Lambda\sim\operatorname{Gamma}\!\left(r,\frac p q\right),
 \qquad
 X\mid\Lambda\sim\operatorname{Poisson}(\Lambda),
\f]

where the gamma parameter shown is a rate. Overflowing geometric sums saturate at the integer
result type's maximum.

`fit<PolyaRand>` starts from the method-of-moments shape

\f[
 r_0=\frac{\bar x^2}{s^2-\bar x},
 \qquad
 p_0=\frac{r_0}{r_0+\bar x},
\f]

then solves the likelihood score for \f$\log r\f$ with guarded Newton steps and sets
\f$\widehat p=\widehat r/(\widehat r+\bar x)\f$. It requires \f$s^2>\bar x\f$.
With fixed \f$r\f$ and prior \f$p\sim\operatorname{Beta}(a,b)\f$, `fit_bayes` returns

\f[
 p\mid x_1,\ldots,x_m
 \sim\operatorname{Beta}\!\left(a+rm,\ b+\sum_jx_j\right).
\f]

There is no dedicated Shannon-entropy member.

### Logarithmic-series distribution (`LogarithmicRand`)

For \f$0<p<1\f$ and \f$k\in\{1,2,\ldots\}\f$,

\f[
 P(X=k)=-\frac{p^k}{k\log(1-p)}.
\f]

Writing

\f[
 B_p(a,0)=\int_0^p\frac{t^{a-1}}{1-t}\,dt
         =p^a\sum_{j=0}^{\infty}\frac{p^j}{a+j},
\f]

the implemented tails are

\f[
 F(k)=1+\frac{B_p(k+1,0)}{\log(1-p)},
 \qquad
 S(k)=-\frac{B_p(k+1,0)}{\log(1-p)}.
\f]

The series routine writes out its first three terms and continues until the next term is negligible;
the survival probability is evaluated directly from that positive tail.

The probability-generating and characteristic functions are

\f[
 G(z)=\frac{\log(1-pz)}{\log(1-p)},
 \qquad
 \phi(t)=\frac{\log(1-pe^{it})}{\log(1-p)}.
\f]

With \f$q=1-p\f$, the implemented moments are

\f[
 \operatorname{E}X=-\frac{p}{q\log q},
 \qquad
 \operatorname{Var}X=-\frac{p(p+\log q)}{q^2(\log q)^2}.
\f]

The mode is one and the median is obtained by the generic discrete quantile search. Sampling uses
Kemp's second accelerated generator: two uniforms select the atom at one, the atom at two, or an
inverse logarithmic tail in constant expected time. There is currently no fitting specialization or
dedicated entropy member.

## Sampling without replacement

### Hypergeometric distribution (`HyperGeometricRand`)

`HyperGeometricRand(N, n, K)` counts successes in \f$n\f$ draws without replacement from a
population of size \f$N\f$ containing \f$K\f$ successes. RandLib requires
\f$N>0\f$, \f$1\leq n\leq N\f$, and \f$1\leq K\leq N\f$. Its support is

\f[
 L=\max(0,n+K-N)\leq k\leq U=\min(n,K),
\f]

and its PMF is

\f[
 P(X=k)=\frac{\binom Kk\binom{N-K}{n-k}}{\binom Nn}.
\f]

The implementation evaluates the log PMF with log-factorials. For \f$F(k)\f$, it compares the
number of support points below and above \f$k\f$, then either sums forward from the true lower
support limit \f$L\f$ or subtracts the shorter upper sum from one. The generating functions are the
finite sums

\f[
 G(z)=\frac1{\binom Nn}
 \sum_{j=L}^{U}\binom Kj\binom{N-K}{n-j}z^j,
 \qquad
 \phi(t)=G(e^{it}).
\f]

The generic finite-support characteristic-function routine evaluates this sum numerically. The
principal moments are

\f[
 \operatorname{E}X=\frac{nK}{N},
 \qquad
 \operatorname{Var}X=\frac{nK(N-K)(N-n)}{N^2(N-1)},
\f]

\f[
 \operatorname{Skew}X=
 \frac{(N-2K)(N-2n)\sqrt{N-1}}
 {(N-2)\sqrt{nK(N-K)(N-n)}}.
\f]

Writing \f$A=nK(N-K)(N-n)\f$, the implemented excess kurtosis is

\f[
 \operatorname{ExKurt}X=
 \frac{
  N^2(N-1)\{N(N+1)-6[K(N-K)+n(N-n)]\}
  +6A(5N-6)}
 {A(N-2)(N-3)}.
\f]

The displayed skewness formula applies to nondegenerate populations with \f$N>2\f$, and the
displayed kurtosis formula applies when \f$N>3\f$. For the nondegenerate \f$N=2\f$ and \f$N=3\f$
cases, RandLib evaluates the relevant finite central moment directly; standardized moments are
undefined in deterministic cases. The reported mode is

\f[
 \left\lfloor\frac{(n+1)(K+1)}{N+2}\right\rfloor,
\f]

and the median uses the generic discrete quantile search.

Sampling performs the draws sequentially. After \f$i\f$ objects have been considered and \f$s\f$
successes observed, the next Bernoulli probability is \f$(K-s)/(N-i)\f$; the routine stops early
once every population success has already been selected. Thus the algorithm is exact and uses at
most \f$n\f$ engine draws. There is currently no fitting specialization or dedicated entropy
member.

### Negative-hypergeometric distribution (`NegativeHyperGeometricRand`)

`NegativeHyperGeometricRand(N, M, m)` counts failures observed before the \f$m\f$th success when
sampling without replacement from \f$N\f$ objects containing \f$M\f$ successes. The constraints
are \f$N>0\f$ and \f$1\leq m\leq M\leq N\f$; the support is
\f$k\in\{0,\ldots,N-M\}\f$. Its PMF is

\f[
 P(X=k)=
 \frac{\binom{k+m-1}{k}\binom{N-m-k}{M-m}}{\binom NM}.
\f]

As in the ordinary hypergeometric implementation, the log PMF uses log-factorials and the CDF
chooses between a forward sum and subtracting the shorter upper sum from one. Its generating
functions are the finite sums

\f[
 G(z)=\frac1{\binom NM}
 \sum_{j=0}^{N-M}
 \binom{j+m-1}{j}\binom{N-m-j}{M-m}z^j,
 \qquad
 \phi(t)=G(e^{it}),
\f]

which the generic finite-support routine evaluates numerically. The implemented moments are

\f[
 \operatorname{E}X=\frac{m(N-M)}{M+1},
\f]

\f[
 \operatorname{Var}X=
 \frac{m(N+1)(N-M)(M-m+1)}{(M+1)^2(M+2)}.
\f]

The median and mode use the generic discrete probability searches; skewness, kurtosis, entropy,
and fitting do not have dedicated implementations.

Sampling updates the remaining-success fraction after each draw and stops at the \f$m\f$th
success. If every failure must occur before that stopping success, it returns \f$N-M\f$ early.
The method is exact and requires at most \f$N\f$ sequential Bernoulli draws.

## Discrete power-law families

### Zeta distribution (`ZetaRand`)

For exponent \f$s>1\f$ and \f$k\in\{1,2,\ldots\}\f$,

\f[
 P(X=k)=\frac{k^{-s}}{\zeta(s)},
\f]

where \f$\zeta\f$ is the Riemann zeta function. With the generalized harmonic number
\f$H_{k,s}=\sum_{j=1}^k j^{-s}\f$,

\f[
 F(k)=\frac{H_{k,s}}{\zeta(s)},
 \qquad
 S(k)=\frac{\zeta(s,k+1)}{\zeta(s)},
\f]

where the second numerator is the Hurwitz zeta tail. Up to \f$k=10^5\f$, the implementation forms
\f$F\f$ from the harmonic sum and \f$S=1-F\f$. Beyond that threshold it evaluates \f$S\f$
directly with a ten-correction Euler--Maclaurin expansion and obtains \f$F=1-S\f$.

The generating functions are

\f[
 G(z)=\frac{\operatorname{Li}_s(z)}{\zeta(s)},
 \qquad
 \phi(t)=\frac{\operatorname{Li}_s(e^{it})}{\zeta(s)}.
\f]

RandLib's generic discrete CF evaluator sums outward from the median and applies a short Euler
transform to the oscillatory remainder; it does not call a complex polylogarithm. More generally,

\f[
 \operatorname{E}(X^m)=
 \begin{cases}
 \dfrac{\zeta(s-m)}{\zeta(s)},&s>m+1,\\
 +\infty,&s\leq m+1.
 \end{cases}
\f]

Thus

\f[
 \operatorname{E}X=\frac{\zeta(s-1)}{\zeta(s)}\quad(s>2),
\f]

\f[
 \operatorname{Var}X=
 \frac{\zeta(s-2)}{\zeta(s)}-
 \left(\frac{\zeta(s-1)}{\zeta(s)}\right)^2
 \quad(s>3),
\f]

and the absent moments are reported as infinity. Skewness and excess kurtosis are assembled from
the first four raw zeta-ratio moments and therefore require \f$s>4\f$ and \f$s>5\f$,
respectively. The mode is one; the median is found by the generic discrete quantile search.

Sampling uses Devroye's rejection method with a rounded Pareto proposal of shape \f$s-1\f$. For
moderate proposals the acceptance ratio is evaluated directly; above \f$10^4\f$ it switches to
`log1p` and `expm1` to avoid cancellation. Accepted values beyond the result type saturate at its
maximum. There is currently no fitting specialization or dedicated entropy member.

### Truncated Zipf distribution (`ZipfRand`)

`ZipfRand(s, n)` is the zeta power law truncated to \f$\{1,\ldots,n\}\f$. RandLib requires
\f$s>1\f$ and \f$n>0\f$. With \f$H_{n,s}=\sum_{j=1}^n j^{-s}\f$,

\f[
 P(X=k)=\frac{k^{-s}}{H_{n,s}},
 \qquad
 F(k)=\frac{H_{k,s}}{H_{n,s}}\quad(1\leq k<n).
\f]

Its generating functions are the finite sums

\f[
 G(z)=\frac1{H_{n,s}}\sum_{k=1}^n\frac{z^k}{k^s},
 \qquad
 \phi(t)=\frac1{H_{n,s}}\sum_{k=1}^n\frac{e^{ikt}}{k^s}.
\f]

The CF implementation evaluates the latter sum directly. All raw moments exist and satisfy

\f[
 \operatorname{E}(X^m)=\frac{H_{n,s-m}}{H_{n,s}}.
\f]

In particular,

\f[
 \operatorname{E}X=\frac{H_{n,s-1}}{H_{n,s}},
 \qquad
 \operatorname{Var}X=
 \frac{H_{n,s-2}}{H_{n,s}}-
 \left(\frac{H_{n,s-1}}{H_{n,s}}\right)^2.
\f]

Skewness and excess kurtosis are computed from these raw moments through order four. They are
undefined for the degenerate case \f$n=1\f$. The mode is one and the median uses the generic
discrete quantile search.

Sampling is exact CDF inversion. Construction caches the first \f$\min(n,16)\f$ cumulative masses;
a uniform in that prefix is resolved from the table, while a uniform in the remaining tail extends
the cumulative sum one integer at a time. There is currently no fitting specialization or dedicated
entropy member.

### Yule--Simon distribution (`YuleRand`)

For shape \f$\rho>0\f$ and \f$k\in\{1,2,\ldots\}\f$,

\f[
 P(X=k)=\rho B(k,\rho+1)
 =\rho\frac{\Gamma(\rho+1)\Gamma(k)}{\Gamma(k+\rho+1)}.
\f]

RandLib evaluates this expression in the log domain. Its direct tails are

\f[
 S(k)=kB(k,\rho+1),
 \qquad
 F(k)=1-kB(k,\rho+1),
\f]

with the CDF formed by `expm1` so a small lower tail is retained. The probability-generating and
characteristic functions are

\f[
 G(z)=\frac{\rho z}{\rho+1}
 {}_2F_1(1,1;\rho+2;z),
 \qquad
 \phi(t)=G(e^{it}).
\f]

The generic discrete CF evaluator uses probability summation with an accelerated oscillatory tail,
rather than evaluating the hypergeometric function. The principal summaries are

\f[
 \operatorname{E}X=
 \begin{cases}\dfrac{\rho}{\rho-1},&\rho>1,\\+\infty,&\rho\leq1,\end{cases}
\f]

\f[
 \operatorname{Var}X=
 \begin{cases}
 \dfrac{\rho^2}{(\rho-1)^2(\rho-2)},&\rho>2,\\
 +\infty,&\rho\leq2,
 \end{cases}
\f]

\f[
 \operatorname{Skew}X=
 \frac{(\rho+1)^2\sqrt{\rho-2}}{\rho(\rho-3)}
 \quad(\rho>3),
\f]

\f[
 \operatorname{ExKurt}X=
 \rho+3+
 \frac{\rho(11\rho^2-49)-22}{\rho(\rho-4)(\rho-3)}
 \quad(\rho>4).
\f]

Unavailable moments are reported as infinity. The mode is one, and the median uses the generic
discrete quantile search.

Sampling uses the exact mixture

\f[
 Y\sim\operatorname{Pareto}(\rho,1),
 \qquad
 X=1+\operatorname{Geometric}(1/Y).
\f]

Both constructed and static samplers implement that transformation. There is currently no fitting
specialization or dedicated entropy member.
