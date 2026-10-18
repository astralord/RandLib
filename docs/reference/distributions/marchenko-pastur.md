# Marchenko--Pastur distribution {#marchenko_pastur_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Marchenko--Pastur distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Spectral and goodness-of-fit laws</span>
  <span role="listitem"><strong>C++ type</strong><code>MarchenkoPasturRand</code></span>
</div>

<div class="distribution-plots"><div class="distribution-plot"><img src="images/marchenko_pastur_pdf.png" loading="lazy" alt="Representative Marchenko--Pastur distribution density curves"><span class="distribution-plot__caption">Representative density curves.</span></div><div class="distribution-plot"><img src="images/marchenko_pastur_cdf.png" loading="lazy" alt="Representative Marchenko--Pastur distribution cumulative distribution curves"><span class="distribution-plot__caption">Corresponding cumulative distribution functions.</span></div></div>

<div class="plot-note"><span>The legends give the parameter sets used in the archived reference plots; formulas and current constructor conventions below are authoritative.</span></div>

## Definition and properties

`MarchenkoPasturRand(ratio, scale)` uses \f$\lambda,s>0\f$. Define

\f[
 a=(1-\sqrt\lambda)^2,\qquad b=(1+\sqrt\lambda)^2.
\f]

Its continuous density on \f$[sa,sb]\f$ is

\f[
 f(x)=\frac{\sqrt{(b-x/s)(x/s-a)}}{2\pi\lambda x}.
\f]

When \f$\lambda>1\f$ the distribution also has an atom

\f[
 P(X=0)=1-\frac1\lambda,
\f]

and the continuous component has total mass \f$1/\lambda\f$. `f(0)` returns infinity as an API
sentinel for that atom; it is not a Lebesgue density there. A compact definition of the CDF is

\f[
 F(x)=\left(1-\frac1\lambda\right)_+\mathbf1_{[0,\infty)}(x)
 +\int_{sa}^{\min(x,sb)} f(u)\,du.
\f]

## Numerical evaluation and quantiles

RandLib evaluates the integral's closed antiderivative with square roots and two branch-aware
`atan` terms. It has separate formulas for \f$\lambda\leq1\f$ and \f$\lambda>1\f$, computes
the latter's survival probability directly, and avoids cancellation at the lower edge.

The raw moments form a Narayana polynomial. For \f$n\geq1\f$,

\f[
 \operatorname{E}X^n=s^n
 \sum_{k=0}^{n-1}\frac1{k+1}{n\choose k}{n-1\choose k}\lambda^k.
\f]

In particular,

\f[
 \operatorname{E}X=s,\qquad
 \operatorname{Var}X=s^2\lambda,\qquad
 \operatorname{Skew}X=\sqrt\lambda,
 \qquad
 \operatorname{ExKurt}X=\lambda-1.
\f]

For \f$\lambda\leq1\f$ the continuous-density mode is

\f[
 \operatorname{Mode}X=s\frac{(\lambda-1)^2}{\lambda+1};
\f]

for \f$\lambda>1\f$ the atom makes zero the mode reported by the API.

## Sampling

Sampling uses beta rejection envelopes selected by ratio. For \f$\lambda<0.3\f$ or
\f$\lambda>3.3\f$, the proposal is a shifted/scaled Beta \f$(3/2,3/2)\f$ law on \f$[a,b]\f$;
intermediate ratios use Beta \f$(1/2,3/2)\f$ on \f$[0,b]\f$. When \f$\lambda>1\f$, an initial
uniform draw returns the zero atom with probability \f$1-1/\lambda\f$.

Quantiles account for the zero atom before applying TOMS 748 to the continuous component. The
characteristic function is evaluated from

\f[
 \phi(t)=\left(1-\frac1\lambda\right)_+
 +\int_{sa}^{sb}e^{itx}f(x)\,dx.
\f]

## Fitting

At \f$\lambda=1\f$ the real integrand is rewritten as
\f$1+(\cos(tx)-1)\f$ to remove the density singularity; for \f$\lambda>1\f$ the atom is added
explicitly. There is currently no `fit<MarchenkoPasturRand>` specialization or entropy member.
