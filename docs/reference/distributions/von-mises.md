# Von Mises distribution {#von_mises_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="other_distribution_reference.html">Circular and singular distributions</a>
  <span aria-hidden="true">/</span>
  <span>Von Mises distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Circular or singular</span>
  <span role="listitem"><strong>Family</strong>Circular distributions</span>
  <span role="listitem"><strong>C++ type</strong><code>VonMisesRand</code></span>
</div>

## Definition and properties

`VonMisesRand(mu, kappa)` requires \f$\kappa>0\f$ and represents angles on
\f$[\mu-\pi,\mu+\pi]\f$. Its density is

\f[
 f(x)=\frac{\exp\{\kappa\cos(x-\mu)\}}{2\pi I_0(\kappa)},
\f]

where \f$I_\nu\f$ is the modified Bessel function of the first kind. The CDF has no elementary
closed form. Below the implementation crossover \f$\kappa=10.5\f$, RandLib evaluates a backwards
trigonometric recurrence whose order is selected from \f$\kappa\f$. At larger concentration it
uses a normal/`erfc` approximation designed for the concentrated regime. The large-concentration
survival function has its own complementary-`erfc` branch; below the crossover it is \f$1-F\f$.

## Moments and shape

At integer harmonic \f$m\f$, the trigonometric moments are

\f[
 \mathbb{E}[e^{imX}]=
 e^{im\mu}\frac{I_m(\kappa)}{I_0(\kappa)}.
\f]

This Bessel ratio is not the ordinary characteristic function at noninteger arguments. Accordingly,
`CF(t)` uses it only when \f$t\f$ is an integer; otherwise the shared continuous-distribution path
numerically integrates \f$e^{itx}f(x)\f$ over the representative interval.

The circular summaries are

\f[
 \operatorname{CMean}[X]=\mu,
 \qquad
 \operatorname{CVar}[X]=1-\frac{I_1(\kappa)}{I_0(\kappa)}.
\f]

The linear representative's median and mode are both \f$\mu\f$.

## Sampling

Sampling uses the
ratio-of-uniforms rejection method with a concentration-dependent envelope scale and returns an
angle in \f$[\mu-\pi,\mu+\pi]\f$. There is currently no fitting specialization or dedicated
entropy member.
