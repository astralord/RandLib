# Triangular distribution {#triangular_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Triangular distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Bounded algebraic distributions</span>
  <span role="listitem"><strong>C++ type</strong><code>TriangularRand</code></span>
</div>

## Definition and properties

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

## Moments and shape

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

## Fitting

There is currently no `fit<TriangularRand>` specialization.
