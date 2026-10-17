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
