"""High-precision reference implementations of Tier 1 distributions.

Every class follows RandLib's constructor parameterization, so that
``make("GammaRand", [0.3, 2])`` is the same law as ``GammaRand<>(0.3, 2)``.
Formulas come from the textbook definitions, not from RandLib's code.
"""

import mpmath as mp

INF = mp.inf


class Continuous:
    discrete = False
    lower = -INF
    upper = INF

    def logpdf(self, x):
        p = self.pdf(x)
        return mp.log(p) if p > 0 else -INF

    def sf(self, x):
        return 1 - self.cdf(x)

    def moments(self):
        return {}


class Discrete:
    discrete = True
    lower = 0
    upper = INF

    def logpmf(self, k):
        p = self.pmf(k)
        return mp.log(p) if p > 0 else -INF

    def moments(self):
        return {}


def standardized_moments(mean, var, third_central, fourth_central):
    return {
        "Mean": mean,
        "Variance": var,
        "Skewness": third_central / var ** mp.mpf(1.5),
        "ExcessKurtosis": fourth_central / var**2 - 3,
    }


class Normal(Continuous):
    def __init__(self, mu, var):
        self.mu, self.sigma = mp.mpf(mu), mp.sqrt(var)

    def pdf(self, x):
        return mp.npdf(x, self.mu, self.sigma)

    def logpdf(self, x):
        z = (x - self.mu) / self.sigma
        return -z * z / 2 - mp.log(self.sigma) - mp.log(2 * mp.pi) / 2

    def cdf(self, x):
        return mp.erfc((self.mu - x) / (self.sigma * mp.sqrt(2))) / 2

    def sf(self, x):
        return mp.erfc((x - self.mu) / (self.sigma * mp.sqrt(2))) / 2

    def cf(self, t):
        return mp.exp(1j * self.mu * t - (self.sigma * t) ** 2 / 2)

    def moments(self):
        return {"Mean": self.mu, "Variance": self.sigma**2, "Skewness": 0, "ExcessKurtosis": 0,
                "Median": self.mu, "Mode": self.mu}


class LogNormal(Continuous):
    lower = 0

    def __init__(self, mu, var):
        self.mu, self.sigma = mp.mpf(mu), mp.sqrt(var)
        self.normal = Normal(mu, var)

    def pdf(self, x):
        return self.normal.pdf(mp.log(x)) / x if x > 0 else mp.mpf(0)

    def logpdf(self, x):
        return self.normal.logpdf(mp.log(x)) - mp.log(x) if x > 0 else -INF

    def cdf(self, x):
        return self.normal.cdf(mp.log(x)) if x > 0 else mp.mpf(0)

    def sf(self, x):
        return self.normal.sf(mp.log(x)) if x > 0 else mp.mpf(1)

    def cf(self, t):
        """E[exp(itX)] after rotating the integration path onto the positive imaginary axis.

        For t > 0, int_0^inf exp(itx) f(x) dx = int_R exp(-t e^y) phi(y - mu + i pi / 2) dy,
        where phi is the N(0, sigma^2) density continued analytically. Unlike the integral
        along the real axis, this integrand does not oscillate rapidly.
        """
        t = mp.mpf(t)
        if t == 0:
            return mp.mpc(1)
        if t < 0:
            return mp.conj(self.cf(-t))
        mu, s = self.mu, self.sigma
        shift = 1j * mp.pi / 2

        def integrand(y):
            z = y - mu + shift
            return mp.exp(-t * mp.exp(y) - z * z / (2 * s * s)) / (s * mp.sqrt(2 * mp.pi))

        upper = mp.log((mp.mp.dps + 10) * mp.log(10) / t) + 2
        lower = mu - 40 * s
        points = [lower] + [mu + s * k for k in range(-12, 13) if lower < mu + s * k < upper] + [upper]
        return mp.quad(integrand, sorted(set(points)))

    def moments(self):
        s2 = self.sigma**2
        e = mp.exp(s2)
        return {"Mean": mp.exp(self.mu + s2 / 2), "Variance": (e - 1) * mp.exp(2 * self.mu + s2),
                "Skewness": (e + 2) * mp.sqrt(e - 1), "ExcessKurtosis": e**4 + 2 * e**3 + 3 * e**2 - 6,
                "Median": mp.exp(self.mu), "Mode": mp.exp(self.mu - s2)}


class Gamma(Continuous):
    lower = 0

    def __init__(self, shape, rate):
        self.k, self.rate = mp.mpf(shape), mp.mpf(rate)

    def logpdf(self, x):
        if x < 0:
            return -INF
        if x == 0:
            return mp.log(self.rate) if self.k == 1 else (INF if self.k < 1 else -INF)
        return self.k * mp.log(self.rate) + (self.k - 1) * mp.log(x) - self.rate * x - mp.loggamma(self.k)

    def pdf(self, x):
        lp = self.logpdf(x)
        return mp.exp(lp) if lp != -INF else mp.mpf(0)

    def cdf(self, x):
        return mp.gammainc(self.k, 0, self.rate * x, regularized=True) if x > 0 else mp.mpf(0)

    def sf(self, x):
        return mp.gammainc(self.k, self.rate * x, INF, regularized=True) if x > 0 else mp.mpf(1)

    def cf(self, t):
        return (1 - 1j * t / self.rate) ** (-self.k)

    def moments(self):
        k, r = self.k, self.rate
        m = {"Mean": k / r, "Variance": k / r**2, "Skewness": 2 / mp.sqrt(k), "ExcessKurtosis": 6 / k}
        if k >= 1:
            m["Mode"] = (k - 1) / r
        return m


class Beta(Continuous):
    def __init__(self, a, b, lo=0, hi=1):
        self.a, self.b = mp.mpf(a), mp.mpf(b)
        self.lower, self.upper = mp.mpf(lo), mp.mpf(hi)
        self.width = self.upper - self.lower

    def logpdf(self, x):
        if x < self.lower or x > self.upper:
            return -INF
        u = (x - self.lower) / self.width
        if u == 0 or u == 1:
            e = self.a if u == 0 else self.b
            if e == 1:
                return -mp.log(mp.beta(self.a, self.b) * self.width)
            return INF if e < 1 else -INF
        return (self.a - 1) * mp.log(u) + (self.b - 1) * mp.log1p(-u) - mp.log(mp.beta(self.a, self.b)) - mp.log(self.width)

    def pdf(self, x):
        lp = self.logpdf(x)
        return mp.exp(lp) if lp != -INF else mp.mpf(0)

    def cdf(self, x):
        if x <= self.lower:
            return mp.mpf(0)
        if x >= self.upper:
            return mp.mpf(1)
        return mp.betainc(self.a, self.b, 0, (x - self.lower) / self.width, regularized=True)

    def sf(self, x):
        if x <= self.lower:
            return mp.mpf(1)
        if x >= self.upper:
            return mp.mpf(0)
        return mp.betainc(self.b, self.a, 0, (self.upper - x) / self.width, regularized=True)

    def cf(self, t):
        t = mp.mpf(t)
        return mp.exp(1j * t * self.lower) * mp.hyp1f1(self.a, self.a + self.b, 1j * t * self.width)

    def moments(self):
        a, b, w = self.a, self.b, self.width
        mean = a / (a + b)
        var = a * b / ((a + b) ** 2 * (a + b + 1))
        skew = 2 * (b - a) * mp.sqrt(a + b + 1) / ((a + b + 2) * mp.sqrt(a * b))
        kurt = 6 * ((a - b) ** 2 * (a + b + 1) - a * b * (a + b + 2)) / (a * b * (a + b + 2) * (a + b + 3))
        m = {"Mean": self.lower + w * mean, "Variance": w * w * var, "Skewness": skew, "ExcessKurtosis": kurt}
        if a > 1 and b > 1:
            m["Mode"] = self.lower + w * (a - 1) / (a + b - 2)
        return m


class BetaPrime(Continuous):
    lower = 0

    def __init__(self, a, b):
        self.a, self.b = mp.mpf(a), mp.mpf(b)

    def logpdf(self, x):
        if x < 0:
            return -INF
        if x == 0:
            if self.a == 1:
                return -mp.log(mp.beta(self.a, self.b))
            return INF if self.a < 1 else -INF
        return (self.a - 1) * mp.log(x) - (self.a + self.b) * mp.log1p(x) - mp.log(mp.beta(self.a, self.b))

    def pdf(self, x):
        lp = self.logpdf(x)
        return mp.exp(lp) if lp != -INF else mp.mpf(0)

    def cdf(self, x):
        return mp.betainc(self.a, self.b, 0, x / (1 + x), regularized=True) if x > 0 else mp.mpf(0)

    def sf(self, x):
        return mp.betainc(self.b, self.a, 0, 1 / (1 + x), regularized=True) if x > 0 else mp.mpf(1)

    def cf(self, t):
        t = mp.mpf(t)
        if t == 0:
            return mp.mpc(1)
        if t < 0:
            return mp.conj(self.cf(-t))
        return mp.gamma(self.a + self.b) / mp.gamma(self.b) * mp.hyperu(self.a, 1 - self.b, -1j * t)

    def moments(self):
        a, b = self.a, self.b
        m = {}
        if b > 1:
            m["Mean"] = a / (b - 1)
        if b > 2:
            m["Variance"] = a * (a + b - 1) / ((b - 2) * (b - 1) ** 2)
        if b > 3:
            m["Skewness"] = 2 * (2 * a + b - 1) / (b - 3) * mp.sqrt((b - 2) / (a * (a + b - 1)))
        if a >= 1:
            m["Mode"] = (a - 1) / (b + 1)
        return m


class StudentT(Continuous):
    def __init__(self, nu, mu=0, sigma=1):
        self.nu, self.mu, self.sigma = mp.mpf(nu), mp.mpf(mu), mp.mpf(sigma)

    def logpdf(self, x):
        nu = self.nu
        z = (x - self.mu) / self.sigma
        return (mp.loggamma((nu + 1) / 2) - mp.loggamma(nu / 2) - mp.log(nu * mp.pi) / 2 - mp.log(self.sigma)
                - (nu + 1) / 2 * mp.log1p(z * z / nu))

    def pdf(self, x):
        return mp.exp(self.logpdf(x))

    def _tail(self, z):
        return mp.betainc(self.nu / 2, mp.mpf(1) / 2, 0, self.nu / (self.nu + z * z), regularized=True) / 2

    def cdf(self, x):
        z = (x - self.mu) / self.sigma
        return self._tail(z) if z < 0 else 1 - self._tail(z)

    def sf(self, x):
        z = (x - self.mu) / self.sigma
        return self._tail(z) if z > 0 else 1 - self._tail(z)

    def cf(self, t):
        z = mp.sqrt(self.nu) * self.sigma * abs(mp.mpf(t))
        if z == 0:
            return mp.mpc(1)
        half = self.nu / 2
        return mp.exp(1j * self.mu * t) * mp.besselk(half, z) * z**half / (mp.gamma(half) * 2 ** (half - 1))

    def moments(self):
        nu = self.nu
        m = {"Median": self.mu, "Mode": self.mu}
        if nu > 1:
            m["Mean"] = self.mu
        if nu > 2:
            m["Variance"] = self.sigma**2 * nu / (nu - 2)
        if nu > 3:
            m["Skewness"] = 0
        if nu > 4:
            m["ExcessKurtosis"] = 6 / (nu - 4)
        return m


class Cauchy(Continuous):
    def __init__(self, mu, gamma):
        self.mu, self.gamma = mp.mpf(mu), mp.mpf(gamma)

    def pdf(self, x):
        z = (x - self.mu) / self.gamma
        return 1 / (mp.pi * self.gamma * (1 + z * z))

    def cdf(self, x):
        return mp.mpf(1) / 2 + mp.atan((x - self.mu) / self.gamma) / mp.pi

    def sf(self, x):
        return mp.mpf(1) / 2 - mp.atan((x - self.mu) / self.gamma) / mp.pi

    def cf(self, t):
        return mp.exp(1j * self.mu * t - self.gamma * abs(t))

    def moments(self):
        return {"Median": self.mu, "Mode": self.mu}


class Levy(Continuous):
    def __init__(self, mu, gamma):
        self.mu, self.gamma = mp.mpf(mu), mp.mpf(gamma)
        self.lower = self.mu

    def logpdf(self, x):
        y = x - self.mu
        if y <= 0:
            return -INF
        return (mp.log(self.gamma / (2 * mp.pi)) / 2 - 3 * mp.log(y) / 2 - self.gamma / (2 * y))

    def pdf(self, x):
        lp = self.logpdf(x)
        return mp.exp(lp) if lp != -INF else mp.mpf(0)

    def cdf(self, x):
        y = x - self.mu
        return mp.erfc(mp.sqrt(self.gamma / (2 * y))) if y > 0 else mp.mpf(0)

    def sf(self, x):
        y = x - self.mu
        return mp.erf(mp.sqrt(self.gamma / (2 * y))) if y > 0 else mp.mpf(1)

    def cf(self, t):
        return mp.exp(1j * self.mu * t - mp.sqrt(-2j * self.gamma * t))

    def moments(self):
        return {"Mode": self.mu + self.gamma / 3, "Median": self.mu + self.gamma / (2 * mp.erfinv(mp.mpf(1) / 2) ** 2)}


class Stable(Continuous):
    """S1 parameterization (Samorodnitsky & Taqqu), defined through its characteristic function

    alpha != 1: phi(t) = exp(-|g t|^a (1 - i b tan(pi a / 2) sign t) + i mu t)
    alpha == 1: phi(t) = exp(-g |t| (1 + i b (2 / pi) sign t log|t|) + i mu t)

    Densities come from inverting phi numerically (alpha >= 1) or from the convergent
    series at infinity obtained by expanding phi term by term (alpha < 1). Distribution
    functions use the Gil-Pelaez formula (alpha >= 1) or the termwise integral of that
    series (alpha < 1).
    """

    def __init__(self, alpha, beta, gamma=1, mu=0):
        self.a, self.b = mp.mpf(alpha), mp.mpf(beta)
        self.g, self.mu = mp.mpf(gamma), mp.mpf(mu)
        if self.a < 1 and self.b == 1:
            self.lower = self.mu
        if self.a < 1 and self.b == -1:
            self.upper = self.mu

    def cf(self, t):
        t = mp.mpf(t)
        if t == 0:
            return mp.mpc(1)
        a, b, g = self.a, self.b, self.g
        s = mp.sign(t)
        if a == 1:
            psi = g * abs(t) * (1 + 1j * b * 2 / mp.pi * s * mp.log(abs(t)))
        else:
            psi = abs(g * t) ** a * (1 - 1j * b * mp.tan(mp.pi * a / 2) * s)
        return mp.exp(-psi + 1j * self.mu * t)

    def _standardize(self, x):
        """Z = (X - mu) / g, except that for alpha == 1 the S1 scale family also shifts:
        X = g Z + (2 / pi) b g log g + mu."""
        z = (mp.mpf(x) - self.mu) / self.g
        if self.a == 1:
            z -= 2 / mp.pi * self.b * mp.log(self.g)
        return z

    def _polar(self, beta):
        """1 - i beta tan(pi a / 2) = c exp(-i d)."""
        tau = beta * mp.tan(mp.pi * self.a / 2)
        return mp.sqrt(1 + tau * tau), mp.atan(tau)

    def _quad(self, integrand, z):
        cutoff = ((mp.mp.dps + 12) * mp.log(10)) ** (1 / self.a)
        period = 2 * mp.pi / max(abs(z), mp.mpf(1))
        pieces = int(min(cutoff / period, 20000)) + 1
        points = [mp.mpf(0)] + [cutoff * k / pieces for k in range(1, pieces + 1)]
        return mp.quad(integrand, points)

    def _pdf_cf(self, z):
        """f(z) = (1 / pi) int_0^inf Re[exp(-itz) phi(t)] dt."""
        if self.a == 1:
            b = self.b

            def integrand(t):
                return mp.exp(-t) * mp.cos(t * z + b * 2 / mp.pi * t * mp.log(t)) if t > 0 else mp.mpf(1)
        else:
            c, d = self._polar(self.b)

            def integrand(t):
                u = c * t ** self.a
                return mp.exp(-u * mp.cos(d)) * mp.cos(t * z - u * mp.sin(d))
        return self._quad(integrand, z) / mp.pi

    def _cdf_cf(self, z):
        """F(z) = 1/2 - (1 / pi) int_0^inf Im[exp(-itz) phi(t)] / t dt."""
        if self.a == 1:
            b = self.b

            def integrand(t):
                return -mp.exp(-t) * mp.sin(t * z + b * 2 / mp.pi * t * mp.log(t)) / t
        else:
            c, d = self._polar(self.b)

            def integrand(t):
                u = c * t ** self.a
                return mp.exp(-u * mp.cos(d)) * mp.sin(u * mp.sin(d) - t * z) / t
        return mp.mpf(1) / 2 - self._quad(integrand, z) / mp.pi

    def _series_at_infinity(self, z, beta, integrate):
        """Convergent series for alpha < 1 and z > 0 (density, or survival when integrate=True).

        Expanding phi(t) = sum_n (-c e^{-i d} t^a)^n / n! and integrating term by term against
        exp(-itz) gives f(z) = (1 / pi z) sum_n (-c)^n Gamma(n a + 1) / n! z^{-n a} cos(n d + pi (n a + 1) / 2).
        """
        a = self.a
        c, d = self._polar(beta)
        log_zc = mp.log(z) * a - mp.log(c)
        n_peak = self._series_peak(z, beta)
        total = mp.mpf(0)
        eps = mp.mpf(2) ** (-mp.mp.prec)
        n = 1
        while True:
            log_term = mp.loggamma(n * a + 1) - mp.loggamma(n + 1) - n * log_zc
            if integrate:
                log_term -= mp.log(n * a)
            phase = n * (mp.pi - d) - mp.pi * (n * a + 1) / 2
            total += mp.exp(log_term) * mp.cos(phase)
            if n > n_peak + 5 and mp.exp(log_term) < eps * abs(total):
                break
            n += 1
        return total / (mp.pi if integrate else mp.pi * z)

    def _series_peak(self, z, beta):
        c, _ = self._polar(beta)
        r = c * self.a**self.a * mp.mpf(z) ** (-self.a)
        return r ** (1 / (1 - self.a)) if r > 1 else mp.mpf(0)

    def _series(self, z, beta, integrate):
        """Series value checked at two working precisions; None if they disagree or are too costly."""
        extra = int(self._series_peak(z, beta) * (1 - self.a) / mp.log(10)) + 10
        if extra > 800:
            return None
        with mp.workdps(mp.mp.dps + extra):
            first = self._series_at_infinity(z, beta, integrate)
        with mp.workdps(mp.mp.dps + extra + 40):
            second = self._series_at_infinity(z, beta, integrate)
        if second == 0 or abs(first - second) > abs(second) * mp.mpf(10) ** (-mp.mp.dps + 5):
            return None
        return +second

    def _asymptotic(self, z, beta, integrate):
        """Series at infinity for alpha > 1, where it is asymptotic rather than convergent: sum
        while the terms decrease and accept only if the smallest term is negligible."""
        a = self.a
        c, d = self._polar(beta)
        log_zc = mp.log(z) * a - mp.log(c)
        total = mp.mpf(0)
        previous = mp.inf
        for n in range(1, 400):
            log_term = mp.loggamma(n * a + 1) - mp.loggamma(n + 1) - n * log_zc
            if integrate:
                log_term -= mp.log(n * a)
            size = mp.exp(log_term)
            if size >= previous:
                return None
            phase = n * (mp.pi - d) - mp.pi * (n * a + 1) / 2
            total += size * mp.cos(phase)
            if total != 0 and size < abs(total) * mp.mpf(10) ** (-mp.mp.dps - 2):
                return total / (mp.pi if integrate else mp.pi * z)
            previous = size
        return None

    def _rotated_unity(self, z, b, integrate):
        """alpha = 1, z > 0 and skewness b >= 0 towards z: rotating the inversion integral onto the
        imaginary axis gives f(z) = (1/pi) int_0^inf exp(-sz - (2b/pi) s log s) sin((1 + b) s) ds,
        and the survival function with an extra 1/s. For b = 0 this is the Cauchy density."""
        def integrand(s_):
            if s_ == 0:
                return (1 + b) if integrate else mp.mpf(0)
            v = mp.exp(-s_ * z - 2 * b / mp.pi * s_ * mp.log(s_)) * mp.sin((1 + b) * s_)
            return v / s_ if integrate else v
        scale = 1 / z
        return mp.quad(integrand, [0, scale / 10, scale, 10 * scale, 100 * scale, mp.inf]) / mp.pi

    def _tail(self, z, beta, integrate):
        """Density (or survival function when integrate) at z > 0 of the standardized law with
        skewness beta, by a method suited to the tails; None where none applies."""
        if self.a < 1:
            if beta == -1:
                return mp.mpf(0)
            return self._series(z, beta, integrate)
        if self.a > 1 and z > 5:
            return self._asymptotic(z, beta, integrate)
        if self.a == 1 and beta >= 0 and z > 5:
            return self._rotated_unity(z, beta, integrate)
        return None

    def _std_pdf(self, z):
        if self.a < 1 and abs(self.b) == 1 and z == 0:
            return mp.mpf(0)
        if self.a != 1 and z == 0:
            return self.pdf_at_zero()
        if z != 0:
            value = self._tail(abs(z), self.b if z > 0 else -self.b, integrate=False)
            if value is not None:
                return value
        return self._pdf_cf(z)

    def _std_tails(self, z):
        """(F(z), S(z)) for the standardized law, each computed directly where possible."""
        if z == 0 and self.a != 1:
            c, d = self._polar(self.b)
            if self.a < 1 and abs(self.b) == 1:
                return (mp.mpf(0), mp.mpf(1)) if self.b == 1 else (mp.mpf(1), mp.mpf(0))
            low = mp.mpf(1) / 2 - d / (mp.pi * self.a)
            return low, 1 - low
        if z > 0:
            upper = self._tail(z, self.b, integrate=True)
            if upper is None:
                low = self._cdf_cf(z)
                return low, 1 - low
            return 1 - upper, upper
        lower = self._tail(-z, -self.b, integrate=True)
        if lower is None:
            low = self._cdf_cf(z)
            return low, 1 - low
        return lower, 1 - lower

    def pdf(self, x):
        if x < self.lower or x > self.upper:
            return mp.mpf(0)
        return self._std_pdf(self._standardize(x)) / self.g

    def sf(self, x):
        if x <= self.lower:
            return mp.mpf(1)
        if x >= self.upper:
            return mp.mpf(0)
        return self._std_tails(self._standardize(x))[1]

    def cdf(self, x):
        if x <= self.lower:
            return mp.mpf(0)
        if x >= self.upper:
            return mp.mpf(1)
        return self._std_tails(self._standardize(x))[0]

    def moments(self):
        m = {}
        if self.a > 1:
            m["Mean"] = self.mu
        if self.a == 2:
            m["Variance"] = 2 * self.g**2
        if self.b == 0:
            m["Median"] = self.mu
            m["Mode"] = self.mu
        return m

    def pdf_at_zero(self):
        """Exact standardized density at 0 for alpha != 1: Gamma(1 + 1/a) c^(-1/a) cos(d / a) / pi."""
        c, d = self._polar(self.b)
        return mp.gamma(1 + 1 / self.a) * c ** (-1 / self.a) * mp.cos(d / self.a) / mp.pi


class Pareto(Continuous):
    def __init__(self, shape, scale):
        self.a, self.s = mp.mpf(shape), mp.mpf(scale)
        self.lower = self.s

    def logpdf(self, x):
        if x < self.s:
            return -INF
        return mp.log(self.a) + self.a * mp.log(self.s) - (self.a + 1) * mp.log(x)

    def pdf(self, x):
        return mp.exp(self.logpdf(x)) if x >= self.s else mp.mpf(0)

    def cdf(self, x):
        return 1 - (self.s / x) ** self.a if x > self.s else mp.mpf(0)

    def sf(self, x):
        return (self.s / x) ** self.a if x > self.s else mp.mpf(1)

    def cf(self, t):
        t = mp.mpf(t)
        if t == 0:
            return mp.mpc(1)
        z = -1j * self.s * t
        return self.a * z**self.a * mp.gammainc(-self.a, z)

    def moments(self):
        a, s = self.a, self.s
        m = {"Mode": s, "Median": s * mp.mpf(2) ** (1 / a)}
        if a > 1:
            m["Mean"] = a * s / (a - 1)
        if a > 2:
            m["Variance"] = s * s * a / ((a - 1) ** 2 * (a - 2))
        if a > 3:
            m["Skewness"] = 2 * (1 + a) / (a - 3) * mp.sqrt((a - 2) / a)
        if a > 4:
            m["ExcessKurtosis"] = 6 * (a**3 + a**2 - 6 * a - 2) / (a * (a - 3) * (a - 4))
        return m


class Weibull(Continuous):
    def cf(self, t):
        t = mp.mpf(t)
        if t == 0:
            return mp.mpc(1)
        if t < 0:
            return mp.conj(self.cf(-t))
        lam, k = self.lam, self.k
        if k <= 1:
            # rotate onto the positive imaginary axis: ∫ exp(itx) f(x) dx = i ∫ exp(-ts) f(is) ds,
            # valid since exp(-(x/λ)^k) decays in the whole first quadrant for k ≤ 1
            def rotated(s_):
                z = 1j * s_ / lam
                return 1j * mp.exp(-t * s_) * (k / lam) * z ** (k - 1) * mp.exp(-z**k)
            return mp.quad(rotated, [0, lam / 10, lam, 10 * lam, mp.inf])
        end = lam * (mp.log(10) * (mp.mp.dps + 20)) ** (1 / k)
        return head_cf(self.pdf, t, mp.mpf(0), end)

    lower = 0

    def __init__(self, scale, shape):
        self.lam, self.k = mp.mpf(scale), mp.mpf(shape)

    def logpdf(self, x):
        if x < 0:
            return -INF
        if x == 0:
            if self.k == 1:
                return -mp.log(self.lam)
            return INF if self.k < 1 else -INF
        u = x / self.lam
        return mp.log(self.k / self.lam) + (self.k - 1) * mp.log(u) - u**self.k

    def pdf(self, x):
        lp = self.logpdf(x)
        return mp.exp(lp) if lp != -INF else mp.mpf(0)

    def cdf(self, x):
        return -mp.expm1(-((x / self.lam) ** self.k)) if x > 0 else mp.mpf(0)

    def sf(self, x):
        return mp.exp(-((x / self.lam) ** self.k)) if x > 0 else mp.mpf(1)


    def moments(self):
        lam, k = self.lam, self.k
        g = [mp.gamma(1 + i / k) for i in range(5)]
        mean = lam * g[1]
        var = lam**2 * (g[2] - g[1] ** 2)
        third = lam**3 * (g[3] - 3 * g[1] * g[2] + 2 * g[1] ** 3)
        fourth = lam**4 * (g[4] - 4 * g[1] * g[3] + 6 * g[1] ** 2 * g[2] - 3 * g[1] ** 4)
        m = standardized_moments(mean, var, third, fourth)
        m["Median"] = lam * mp.log(2) ** (1 / k)
        m["Mode"] = lam * ((k - 1) / k) ** (1 / k) if k > 1 else mp.mpf(0)
        return m


class Frechet(Continuous):
    def cf(self, t):
        t = mp.mpf(t)
        if t == 0:
            return mp.mpc(1)
        if t < 0:
            return mp.conj(self.cf(-t))
        a, sc = self.a, self.s
        # head numerically; beyond Y expand exp(-(y/s)^-a) and integrate each power exactly
        end = max(sc * mp.mpf(10) ** (1 / a), 50 * 2 * mp.pi / t)
        f0 = lambda y: (a / sc) * (y / sc) ** (-1 - a) * mp.exp(-((y / sc) ** (-a))) if y > 0 else mp.mpf(0)
        value = head_cf(f0, t, mp.mpf(0), end)
        n = 0
        while True:
            term = (-1) ** n / mp.factorial(n) * a * sc ** (a * (n + 1)) * power_tail_cf(t, end, 1 + a * (n + 1))
            value += term
            if abs(term) < mp.mpf(10) ** (-mp.mp.dps - 5):
                break
            n += 1
        return mp.exp(1j * t * self.m) * value

    def __init__(self, shape, scale, loc):
        self.a, self.s, self.m = mp.mpf(shape), mp.mpf(scale), mp.mpf(loc)
        self.lower = self.m

    def logpdf(self, x):
        if x <= self.m:
            return -INF
        u = (x - self.m) / self.s
        return mp.log(self.a / self.s) - (1 + self.a) * mp.log(u) - u ** (-self.a)

    def pdf(self, x):
        lp = self.logpdf(x)
        return mp.exp(lp) if lp != -INF else mp.mpf(0)

    def cdf(self, x):
        return mp.exp(-(((x - self.m) / self.s) ** (-self.a))) if x > self.m else mp.mpf(0)

    def sf(self, x):
        return -mp.expm1(-(((x - self.m) / self.s) ** (-self.a))) if x > self.m else mp.mpf(1)


    def moments(self):
        a, s, loc = self.a, self.s, self.m
        m = {"Median": loc + s / mp.log(2) ** (1 / a), "Mode": loc + s * (a / (1 + a)) ** (1 / a)}
        g = lambda k: mp.gamma(1 - k / a)
        if a > 1:
            m["Mean"] = loc + s * g(1)
        if a > 2:
            m["Variance"] = s * s * (g(2) - g(1) ** 2)
        if a > 3:
            m["Skewness"] = (g(3) - 3 * g(2) * g(1) + 2 * g(1) ** 3) / (g(2) - g(1) ** 2) ** mp.mpf(1.5)
        if a > 4:
            m["ExcessKurtosis"] = (g(4) - 4 * g(3) * g(1) + 3 * g(2) ** 2) / (g(2) - g(1) ** 2) ** 2 - 6
        return m


class Gumbel(Continuous):
    def __init__(self, mu, beta):
        self.mu, self.beta = mp.mpf(mu), mp.mpf(beta)

    def logpdf(self, x):
        z = (x - self.mu) / self.beta
        return -mp.log(self.beta) - z - mp.exp(-z)

    def pdf(self, x):
        return mp.exp(self.logpdf(x))

    def cdf(self, x):
        return mp.exp(-mp.exp(-(x - self.mu) / self.beta))

    def sf(self, x):
        return -mp.expm1(-mp.exp(-(x - self.mu) / self.beta))

    def cf(self, t):
        return mp.gamma(1 - 1j * self.beta * t) * mp.exp(1j * self.mu * t)

    def moments(self):
        b = self.beta
        return {"Mean": self.mu + b * mp.euler, "Variance": mp.pi**2 * b * b / 6,
                "Skewness": 12 * mp.sqrt(6) * mp.zeta(3) / mp.pi**3, "ExcessKurtosis": mp.mpf(12) / 5,
                "Median": self.mu - b * mp.log(mp.log(2)), "Mode": self.mu}


class GeneralizedPareto(Continuous):
    def __init__(self, shape, scale, location):
        self.xi, self.scale, self.location = mp.mpf(shape), mp.mpf(scale), mp.mpf(location)
        self.lower = self.location
        if self.xi < 0:
            self.upper = self.location - self.scale / self.xi

    def logpdf(self, x):
        if x < self.lower or x > self.upper:
            return -INF
        z = (x - self.location) / self.scale
        if self.xi == 0:
            return -mp.log(self.scale) - z
        support = 1 + self.xi * z
        if support == 0:
            if self.xi < -1:
                return INF
            if self.xi == -1:
                return -mp.log(self.scale)
            return -INF
        return -mp.log(self.scale) - (1 / self.xi + 1) * mp.log(support)

    def pdf(self, x):
        return mp.exp(self.logpdf(x))

    def sf(self, x):
        if x <= self.lower:
            return mp.mpf(1)
        if x >= self.upper:
            return mp.mpf(0)
        z = (x - self.location) / self.scale
        return mp.exp(-z) if self.xi == 0 else (1 + self.xi * z) ** (-1 / self.xi)

    def cdf(self, x):
        return 1 - self.sf(x)

    def moments(self):
        xi, scale, location = self.xi, self.scale, self.location
        result = {"Median": self.quantile(mp.mpf(1) / 2),
                  "Mode": self.upper if xi < -1 else location}
        if xi < 1:
            result["Mean"] = location + scale / (1 - xi)
        if xi < mp.mpf(1) / 2:
            result["Variance"] = scale**2 / ((1 - xi) ** 2 * (1 - 2 * xi))
        if xi < mp.mpf(1) / 3:
            result["Skewness"] = 2 * (1 + xi) * mp.sqrt(1 - 2 * xi) / (1 - 3 * xi)
        if xi < mp.mpf(1) / 4:
            result["ExcessKurtosis"] = (3 * (1 - 2 * xi) * (2 * xi**2 + xi + 3)
                                          / ((1 - 3 * xi) * (1 - 4 * xi)) - 3)
        return result

    def quantile(self, p):
        exponential = -mp.log1p(-p)
        return (self.location + self.scale * exponential if self.xi == 0 else
                self.location + self.scale * mp.expm1(self.xi * exponential) / self.xi)


class GeneralizedExtremeValue(Continuous):
    def __init__(self, shape, scale, location):
        self.xi, self.scale, self.location = mp.mpf(shape), mp.mpf(scale), mp.mpf(location)
        if self.xi > 0:
            self.lower = self.location - self.scale / self.xi
        elif self.xi < 0:
            self.upper = self.location - self.scale / self.xi

    def _support(self, x):
        return 1 + self.xi * (x - self.location) / self.scale

    def logpdf(self, x):
        if x < self.lower or x > self.upper:
            return -INF
        z = (x - self.location) / self.scale
        if self.xi == 0:
            return -mp.log(self.scale) - z - mp.exp(-z)
        support = self._support(x)
        if support == 0:
            if self.xi < -1:
                return INF
            if self.xi == -1:
                return -mp.log(self.scale)
            return -INF
        return (-mp.log(self.scale) - (1 / self.xi + 1) * mp.log(support)
                - support ** (-1 / self.xi))

    def pdf(self, x):
        return mp.exp(self.logpdf(x))

    def cdf(self, x):
        if x <= self.lower:
            return mp.mpf(0)
        if x >= self.upper:
            return mp.mpf(1)
        exponent = mp.exp(-(x - self.location) / self.scale) if self.xi == 0 else self._support(x) ** (-1 / self.xi)
        return mp.exp(-exponent)

    def sf(self, x):
        if x <= self.lower:
            return mp.mpf(1)
        if x >= self.upper:
            return mp.mpf(0)
        exponent = mp.exp(-(x - self.location) / self.scale) if self.xi == 0 else self._support(x) ** (-1 / self.xi)
        return -mp.expm1(-exponent)

    def quantile(self, p):
        log_exponential = mp.log(-mp.log(p))
        return (self.location - self.scale * log_exponential if self.xi == 0 else
                self.location + self.scale * mp.expm1(-self.xi * log_exponential) / self.xi)

    def moments(self):
        xi, scale, location = self.xi, self.scale, self.location
        result = {"Median": self.quantile(mp.mpf(1) / 2)}
        if xi <= -1:
            result["Mode"] = self.upper
        elif xi == 0:
            result["Mode"] = location
        else:
            result["Mode"] = location + scale * mp.expm1(-xi * mp.log1p(xi)) / xi
        if xi == 0:
            result.update(Gumbel(location, scale).moments())
            return result
        gamma = lambda k: mp.gamma(1 - k * xi)
        if xi < 1:
            result["Mean"] = location + scale * (gamma(1) - 1) / xi
        if xi < mp.mpf(1) / 2:
            variance = gamma(2) - gamma(1) ** 2
            result["Variance"] = scale**2 * variance / xi**2
        if xi < mp.mpf(1) / 3:
            third = gamma(3) - 3 * gamma(1) * gamma(2) + 2 * gamma(1) ** 3
            result["Skewness"] = mp.sign(xi) * third / variance ** mp.mpf(1.5)
        if xi < mp.mpf(1) / 4:
            fourth = gamma(4) - 4 * gamma(1) * gamma(3) + 6 * gamma(1) ** 2 * gamma(2) - 3 * gamma(1) ** 4
            result["ExcessKurtosis"] = fourth / variance**2 - 3
        return result


class AsymmetricLaplace(Continuous):
    """f(x) = lam / (k + 1/k) * exp(-lam k (x - m)) for x >= m, exp(lam (x - m) / k) otherwise."""

    def __init__(self, m, lam, k=1):
        self.m, self.lam, self.k = mp.mpf(m), mp.mpf(lam), mp.mpf(k)

    def logpdf(self, x):
        y = x - self.m
        lp = mp.log(self.lam / (self.k + 1 / self.k))
        return lp - self.lam * self.k * y if y >= 0 else lp + self.lam * y / self.k

    def pdf(self, x):
        return mp.exp(self.logpdf(x))

    def cdf(self, x):
        y = x - self.m
        k2 = self.k**2
        if y <= 0:
            return k2 / (1 + k2) * mp.exp(self.lam * y / self.k)
        return 1 - mp.exp(-self.lam * self.k * y) / (1 + k2)

    def sf(self, x):
        y = x - self.m
        k2 = self.k**2
        if y <= 0:
            return 1 - k2 / (1 + k2) * mp.exp(self.lam * y / self.k)
        return mp.exp(-self.lam * self.k * y) / (1 + k2)

    def cf(self, t):
        return mp.exp(1j * self.m * t) / (1 + 1j * t / self.lam * (self.k - 1 / self.k) + t * t / self.lam**2)

    def moments(self):
        lam, k = self.lam, self.k
        mean = self.m + (1 - k * k) / (lam * k)
        var = (1 + k**4) / (lam**2 * k**2)
        skew = 2 * (1 - k**6) / (k**4 + 1) ** mp.mpf(1.5)
        kurt = 6 * (1 + k**8) / (1 + k**4) ** 2
        return {"Mean": mean, "Variance": var, "Skewness": skew, "ExcessKurtosis": kurt, "Mode": self.m}


class Logistic(Continuous):
    def __init__(self, mu, s):
        self.mu, self.s = mp.mpf(mu), mp.mpf(s)

    def logpdf(self, x):
        z = -abs(x - self.mu) / self.s
        return z - 2 * mp.log1p(mp.exp(z)) - mp.log(self.s)

    def pdf(self, x):
        return mp.exp(self.logpdf(x))

    def cdf(self, x):
        return 1 / (1 + mp.exp(-(x - self.mu) / self.s))

    def sf(self, x):
        return 1 / (1 + mp.exp((x - self.mu) / self.s))

    def cf(self, t):
        st = self.s * t
        return mp.exp(1j * self.mu * t) * (mp.pi * st / mp.sinh(mp.pi * st) if st != 0 else 1)

    def moments(self):
        return {"Mean": self.mu, "Variance": (self.s * mp.pi) ** 2 / 3, "Skewness": 0,
                "ExcessKurtosis": mp.mpf(6) / 5, "Median": self.mu, "Mode": self.mu}


class InverseGamma(Continuous):
    lower = 0

    def __init__(self, shape, rate):
        self.a, self.b = mp.mpf(shape), mp.mpf(rate)

    def logpdf(self, x):
        if x <= 0:
            return -INF
        return self.a * mp.log(self.b) - mp.loggamma(self.a) - (self.a + 1) * mp.log(x) - self.b / x

    def pdf(self, x):
        return mp.exp(self.logpdf(x)) if x > 0 else mp.mpf(0)

    def cdf(self, x):
        return mp.gammainc(self.a, self.b / x, INF, regularized=True) if x > 0 else mp.mpf(0)

    def sf(self, x):
        return mp.gammainc(self.a, 0, self.b / x, regularized=True) if x > 0 else mp.mpf(1)

    def cf(self, t):
        t = mp.mpf(t)
        if t == 0:
            return mp.mpc(1)
        z = -1j * self.b * t
        return 2 * z ** (self.a / 2) / mp.gamma(self.a) * mp.besselk(self.a, mp.sqrt(4 * z))

    def moments(self):
        a, b = self.a, self.b
        m = {"Mode": b / (a + 1)}
        if a > 1:
            m["Mean"] = b / (a - 1)
        if a > 2:
            m["Variance"] = b * b / ((a - 1) ** 2 * (a - 2))
        if a > 3:
            m["Skewness"] = 4 * mp.sqrt(a - 2) / (a - 3)
        if a > 4:
            m["ExcessKurtosis"] = 6 * (5 * a - 11) / ((a - 3) * (a - 4))
        return m


class InverseGaussian(Continuous):
    lower = 0

    def __init__(self, mean, shape):
        self.mu, self.lam = mp.mpf(mean), mp.mpf(shape)

    def logpdf(self, x):
        if x <= 0:
            return -INF
        return (mp.log(self.lam / (2 * mp.pi)) / 2 - 3 * mp.log(x) / 2
                - self.lam * (x - self.mu) ** 2 / (2 * self.mu**2 * x))

    def pdf(self, x):
        return mp.exp(self.logpdf(x)) if x > 0 else mp.mpf(0)

    def cdf(self, x):
        if x <= 0:
            return mp.mpf(0)
        r = mp.sqrt(self.lam / x)
        a = mp.erfc(-r * (x / self.mu - 1) / mp.sqrt(2)) / 2
        b = mp.exp(2 * self.lam / self.mu) * mp.erfc(r * (x / self.mu + 1) / mp.sqrt(2)) / 2
        return a + b

    def sf(self, x):
        if x <= 0:
            return mp.mpf(1)
        r = mp.sqrt(self.lam / x)
        a = mp.erfc(r * (x / self.mu - 1) / mp.sqrt(2)) / 2
        b = mp.exp(2 * self.lam / self.mu) * mp.erfc(r * (x / self.mu + 1) / mp.sqrt(2)) / 2
        return a - b

    def cf(self, t):
        return mp.exp(self.lam / self.mu * (1 - mp.sqrt(1 - 2j * self.mu**2 * t / self.lam)))

    def moments(self):
        mu, lam = self.mu, self.lam
        return {"Mean": mu, "Variance": mu**3 / lam, "Skewness": 3 * mp.sqrt(mu / lam),
                "ExcessKurtosis": 15 * mu / lam,
                "Mode": mu * (mp.sqrt(1 + 9 * mu * mu / (4 * lam * lam)) - 3 * mu / (2 * lam))}


class Uniform(Continuous):
    def __init__(self, a, b):
        self.lower, self.upper = mp.mpf(a), mp.mpf(b)

    def pdf(self, x):
        return 1 / (self.upper - self.lower) if self.lower <= x <= self.upper else mp.mpf(0)

    def cdf(self, x):
        return min(mp.mpf(1), max(mp.mpf(0), (x - self.lower) / (self.upper - self.lower)))

    def sf(self, x):
        return min(mp.mpf(1), max(mp.mpf(0), (self.upper - x) / (self.upper - self.lower)))

    def cf(self, t):
        t = mp.mpf(t)
        if t == 0:
            return mp.mpc(1)
        return (mp.exp(1j * t * self.upper) - mp.exp(1j * t * self.lower)) / (1j * t * (self.upper - self.lower))

    def moments(self):
        a, b = self.lower, self.upper
        return {"Mean": (a + b) / 2, "Variance": (b - a) ** 2 / 12, "Skewness": 0,
                "ExcessKurtosis": -mp.mpf(6) / 5, "Median": (a + b) / 2}


class FisherF(Continuous):
    lower = 0

    def __init__(self, d1, d2):
        self.d1, self.d2 = mp.mpf(d1), mp.mpf(d2)

    def logpdf(self, x):
        d1, d2 = self.d1, self.d2
        if x < 0:
            return -INF
        if x == 0:
            if d1 == 2:
                return mp.mpf(0)
            return INF if d1 < 2 else -INF
        return (d1 / 2 * mp.log(d1 / d2) + (d1 / 2 - 1) * mp.log(x) - (d1 + d2) / 2 * mp.log1p(d1 * x / d2)
                - mp.log(mp.beta(d1 / 2, d2 / 2)))

    def pdf(self, x):
        lp = self.logpdf(x)
        return mp.exp(lp) if lp != -INF else mp.mpf(0)

    def cdf(self, x):
        if x <= 0:
            return mp.mpf(0)
        return mp.betainc(self.d1 / 2, self.d2 / 2, 0, self.d1 * x / (self.d1 * x + self.d2), regularized=True)

    def sf(self, x):
        if x <= 0:
            return mp.mpf(1)
        return mp.betainc(self.d2 / 2, self.d1 / 2, 0, self.d2 / (self.d1 * x + self.d2), regularized=True)

    def cf(self, t):
        return BetaPrime(self.d1 / 2, self.d2 / 2).cf(self.d2 / self.d1 * mp.mpf(t))

    def moments(self):
        d1, d2 = self.d1, self.d2
        m = {}
        if d2 > 2:
            m["Mean"] = d2 / (d2 - 2)
        if d2 > 4:
            m["Variance"] = 2 * d2**2 * (d1 + d2 - 2) / (d1 * (d2 - 2) ** 2 * (d2 - 4))
        if d2 > 6:
            m["Skewness"] = (2 * d1 + d2 - 2) * mp.sqrt(8 * (d2 - 4)) / ((d2 - 6) * mp.sqrt(d1 * (d1 + d2 - 2)))
        if d1 > 2:
            m["Mode"] = (d1 - 2) / d1 * d2 / (d2 + 2)
        return m


class NoncentralChiSquared(Continuous):
    lower = 0

    def __init__(self, k, lam):
        self.k, self.lam = mp.mpf(k), mp.mpf(lam)

    def _poisson_sum(self, term):
        """Sum of Poisson(λ/2) weights times term(j). In the far tails the dominant terms lie far
        beyond the peak of the weights, so stop only once the terms themselves are decreasing
        and negligible."""
        half = self.lam / 2
        total = mp.mpf(0)
        previous = mp.mpf(0)
        j = 0
        eps = mp.mpf(10) ** (-mp.mp.dps - 5)
        while True:
            w = mp.exp(-half + j * mp.log(half) - mp.loggamma(j + 1))
            t = w * term(j)
            total += t
            if j > half + 10 and t <= previous and t <= eps * total:
                break
            previous = t
            j += 1
        return total

    def pdf(self, x):
        if x <= 0:
            return Gamma(self.k / 2, mp.mpf(1) / 2).pdf(x) * mp.exp(-self.lam / 2) if x == 0 else mp.mpf(0)
        return self._poisson_sum(lambda j: Gamma(self.k / 2 + j, mp.mpf(1) / 2).pdf(x))

    def cdf(self, x):
        if x <= 0:
            return mp.mpf(0)
        return self._poisson_sum(lambda j: mp.gammainc(self.k / 2 + j, 0, x / 2, regularized=True))

    def sf(self, x):
        if x <= 0:
            return mp.mpf(1)
        return self._poisson_sum(lambda j: mp.gammainc(self.k / 2 + j, x / 2, INF, regularized=True))

    def cf(self, t):
        d = 1 - 2j * mp.mpf(t)
        return mp.exp(1j * self.lam * t / d) / d ** (self.k / 2)

    def moments(self):
        k, lam = self.k, self.lam
        var = 2 * (k + 2 * lam)
        return {"Mean": k + lam, "Variance": var, "Skewness": 2 ** mp.mpf(1.5) * (k + 3 * lam) / (k + 2 * lam) ** mp.mpf(1.5),
                "ExcessKurtosis": 12 * (k + 4 * lam) / (k + 2 * lam) ** 2}


class Poisson(Discrete):
    def __init__(self, lam):
        self.lam = mp.mpf(lam)

    def logpmf(self, k):
        if k < 0:
            return -INF
        return -self.lam + k * mp.log(self.lam) - mp.loggamma(k + 1)

    def pmf(self, k):
        return mp.exp(self.logpmf(k)) if k >= 0 else mp.mpf(0)

    def cdf(self, k):
        return mp.gammainc(k + 1, self.lam, INF, regularized=True) if k >= 0 else mp.mpf(0)

    def sf(self, k):
        return mp.gammainc(k + 1, 0, self.lam, regularized=True) if k >= 0 else mp.mpf(1)

    def cf(self, t):
        return mp.exp(self.lam * (mp.exp(1j * t) - 1))

    def moments(self):
        lam = self.lam
        return {"Mean": lam, "Variance": lam, "Skewness": 1 / mp.sqrt(lam), "ExcessKurtosis": 1 / lam}


class Binomial(Discrete):
    def __init__(self, n, p):
        self.n, self.p = int(n), mp.mpf(p)
        self.upper = self.n

    def logpmf(self, k):
        if k < 0 or k > self.n:
            return -INF
        return (mp.loggamma(self.n + 1) - mp.loggamma(k + 1) - mp.loggamma(self.n - k + 1)
                + k * mp.log(self.p) + (self.n - k) * mp.log1p(-self.p))

    def pmf(self, k):
        lp = self.logpmf(k)
        return mp.exp(lp) if lp != -INF else mp.mpf(0)

    def cdf(self, k):
        if k < 0:
            return mp.mpf(0)
        if k >= self.n:
            return mp.mpf(1)
        return mp.betainc(self.n - k, k + 1, 0, 1 - self.p, regularized=True)

    def sf(self, k):
        if k < 0:
            return mp.mpf(1)
        if k >= self.n:
            return mp.mpf(0)
        return mp.betainc(k + 1, self.n - k, 0, self.p, regularized=True)

    def cf(self, t):
        return (1 - self.p + self.p * mp.exp(1j * t)) ** self.n

    def moments(self):
        n, p = self.n, self.p
        q = 1 - p
        return {"Mean": n * p, "Variance": n * p * q, "Skewness": (q - p) / mp.sqrt(n * p * q),
                "ExcessKurtosis": (1 - 6 * p * q) / (n * p * q)}


class NegativeBinomial(Discrete):
    """P(k) = Gamma(r + k) / (Gamma(r) k!) p^r (1 - p)^k for k = 0, 1, ..."""

    def __init__(self, r, p):
        self.r, self.p = mp.mpf(r), mp.mpf(p)

    def logpmf(self, k):
        if k < 0:
            return -INF
        return (mp.loggamma(self.r + k) - mp.loggamma(self.r) - mp.loggamma(k + 1)
                + self.r * mp.log(self.p) + k * mp.log1p(-self.p))

    def pmf(self, k):
        return mp.exp(self.logpmf(k)) if k >= 0 else mp.mpf(0)

    def cdf(self, k):
        return mp.betainc(self.r, k + 1, 0, self.p, regularized=True) if k >= 0 else mp.mpf(0)

    def sf(self, k):
        return mp.betainc(k + 1, self.r, 0, 1 - self.p, regularized=True) if k >= 0 else mp.mpf(1)

    def cf(self, t):
        return (self.p / (1 - (1 - self.p) * mp.exp(1j * t))) ** self.r

    def moments(self):
        r, p = self.r, self.p
        q = 1 - p
        return {"Mean": r * q / p, "Variance": r * q / p**2, "Skewness": (2 - p) / mp.sqrt(q * r),
                "ExcessKurtosis": 6 / r + p * p / (q * r)}


class Bernoulli(Binomial):
    def __init__(self, p):
        super().__init__(1, p)


class UniformDiscrete(Discrete):
    def __init__(self, a, b):
        self.lower, self.upper = int(a), int(b)
        self.count = self.upper - self.lower + 1

    def pmf(self, k):
        return mp.mpf(1) / self.count if self.lower <= k <= self.upper else mp.mpf(0)

    def cdf(self, k):
        return min(mp.mpf(1), max(mp.mpf(0), mp.mpf(k - self.lower + 1) / self.count))

    def sf(self, k):
        return 1 - self.cdf(k)

    def moments(self):
        n = mp.mpf(self.count)
        return {"Mean": mp.mpf(self.lower + self.upper) / 2, "Variance": (n * n - 1) / 12, "Skewness": 0,
                "ExcessKurtosis": -6 * (n * n + 1) / (5 * (n * n - 1))}


def head_cf(pdf, t, lo, end):
    """∫_lo^end exp(itx) pdf(x) dx with breakpoints every quarter period."""
    period = 2 * mp.pi / t
    points = [lo + period * k / 4 for k in range(0, int((end - lo) / period * 4) + 1)] + [end]
    points = sorted(set(x for x in points if x <= end))
    return mp.mpc(mp.quad(lambda x: mp.cos(t * x) * pdf(x), points), mp.quad(lambda x: mp.sin(t * x) * pdf(x), points))


def power_tail_cf(t, start, beta):
    """∫_start^∞ exp(itx) x^(-beta) dx = (-it)^(beta-1) Γ(1 - beta, -it·start) for t > 0."""
    z = -1j * t
    return z ** (beta - 1) * mp.gammainc(1 - beta, z * start)


class Logarithmic(Discrete):
    lower = 1

    def __init__(self, p):
        self.p = mp.mpf(p)
        self.c = -1 / mp.log1p(-self.p)

    def pmf(self, k):
        return self.c * self.p**k / k if k >= 1 else mp.mpf(0)

    def logpmf(self, k):
        return mp.log(self.c) + k * mp.log(self.p) - mp.log(k) if k >= 1 else -INF

    def sf(self, k):
        if k < 1:
            return mp.mpf(1)
        return self.c * self.p ** (k + 1) * mp.lerchphi(self.p, 1, k + 1)

    def cdf(self, k):
        return 1 - self.sf(k) if k >= 1 else mp.mpf(0)

    def cf(self, t):
        return mp.log(1 - self.p * mp.exp(1j * mp.mpf(t))) / mp.log1p(-self.p)

    def moments(self):
        p, c = self.p, self.c
        mean = c * p / (1 - p)
        return {"Mean": mean, "Variance": c * p * (1 - c * p) / (1 - p) ** 2}


class WignerSemicircle(Continuous):
    def __init__(self, r):
        self.r = mp.mpf(r)
        self.lower, self.upper = -self.r, self.r

    def pdf(self, x):
        return 2 * mp.sqrt(self.r**2 - x * x) / (mp.pi * self.r**2) if abs(x) <= self.r else mp.mpf(0)

    def cdf(self, x):
        if x <= -self.r:
            return mp.mpf(0)
        if x >= self.r:
            return mp.mpf(1)
        return mp.mpf(1) / 2 + x * mp.sqrt(self.r**2 - x * x) / (mp.pi * self.r**2) + mp.asin(x / self.r) / mp.pi

    def sf(self, x):
        return self.cdf(-x)

    def cf(self, t):
        z = self.r * mp.mpf(t)
        return mp.mpc(1) if z == 0 else mp.mpc(2 * mp.besselj(1, z) / z)

    def moments(self):
        return {"Mean": 0, "Variance": self.r**2 / 4, "Skewness": 0, "ExcessKurtosis": -1, "Median": 0, "Mode": 0}


class GeometricStable(Continuous):
    """Geometric stable law Y = S(W): given W = w ~ Exp(1), Y ~ S1(alpha, beta, gamma w^(1/alpha), mu w).
    The S1 scale family carries the extra alpha = 1 location shift by itself. Each value is a
    single integral over w of stable values, so this class is expensive."""

    light = True

    def __init__(self, alpha, beta, gamma, mu):
        self.a, self.b, self.g, self.mu = mp.mpf(alpha), mp.mpf(beta), mp.mpf(gamma), mp.mpf(mu)

    def _conditional(self, w):
        return Stable(self.a, self.b, self.g * w ** (1 / self.a), self.mu * w)

    def _mix(self, value):
        return mp.quad(lambda w: mp.exp(-w) * value(self._conditional(w)) if w > 0 else mp.mpf(0),
                       [0, mp.mpf(1) / 100, mp.mpf(1) / 10, 1, 4, 16, mp.inf])

    def pdf(self, x):
        return self._mix(lambda d: d.pdf(x))

    def cdf(self, x):
        return self._mix(lambda d: d.cdf(x))

    def sf(self, x):
        return self._mix(lambda d: d.sf(x))


def make(cls, args):
    # the exact doubles the C++ constructor receives, not their decimal spelling
    a = [mp.mpf(float(v)) for v in args]
    table = {
        "NormalRand": lambda: Normal(*a),
        "LogNormalRand": lambda: LogNormal(*a),
        "ExponentialRand": lambda: Gamma(1, *a),
        "GammaRand": lambda: Gamma(*a),
        "ChiSquaredRand": lambda: Gamma(a[0] / 2, mp.mpf(1) / 2),
        "ErlangRand": lambda: Gamma(*a),
        "BetaRand": lambda: Beta(*a),
        "ArcsineRand": lambda: Beta(1 - a[0], a[0], *a[1:]),
        "BaldingNicholsRand": lambda: Beta(a[1] * (1 - a[0]) / a[0], (1 - a[1]) * (1 - a[0]) / a[0]),
        "BetaPrimeRand": lambda: BetaPrime(*a),
        "StudentTRand": lambda: StudentT(*a),
        "CauchyRand": lambda: Cauchy(*a),
        "LevyRand": lambda: Levy(*a),
        "StableRand": lambda: Stable(*a),
        "HoltsmarkRand": lambda: Stable(mp.mpf(3) / 2, 0, *a),
        "LandauRand": lambda: Stable(1, 1, *a),
        "ParetoRand": lambda: Pareto(*a),
        "WeibullRand": lambda: Weibull(*a),
        "FrechetRand": lambda: Frechet(*a),
        "GumbelRand": lambda: Gumbel(*a),
        "GeneralizedParetoRand": lambda: GeneralizedPareto(*a),
        "GeneralizedExtremeValueRand": lambda: GeneralizedExtremeValue(*a),
        "LaplaceRand": lambda: AsymmetricLaplace(a[0], 1 / a[1], 1),
        "AsymmetricLaplaceRand": lambda: AsymmetricLaplace(a[0], 1 / a[1], a[2]),
        "LogisticRand": lambda: Logistic(*a),
        "InverseGammaRand": lambda: InverseGamma(*a),
        "InverseGaussianRand": lambda: InverseGaussian(*a),
        "UniformRand": lambda: Uniform(*a),
        "FisherFRand": lambda: FisherF(*a),
        "NoncentralChiSquaredRand": lambda: NoncentralChiSquared(*a),
        "PoissonRand": lambda: Poisson(*a),
        "BinomialRand": lambda: Binomial(*a),
        "BernoulliRand": lambda: Bernoulli(*a),
        "NegativeBinomialRand": lambda: NegativeBinomial(*a),
        "PascalRand": lambda: NegativeBinomial(*a),
        "GeometricRand": lambda: NegativeBinomial(1, *a),
        "UniformDiscreteRand": lambda: UniformDiscrete(*a),
        "GeometricStableRand": lambda: GeometricStable(*a),
        "LogarithmicRand": lambda: Logarithmic(*a),
        "WignerSemicircleRand": lambda: WignerSemicircle(*a),
    }
    factory = table.get(cls)
    return factory() if factory else None
