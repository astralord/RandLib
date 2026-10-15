#!/usr/bin/env python3
"""Generate deterministic datasets and SciPy MLE references for fitting tests."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from scipy import stats


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "tests" / "data" / "fitting"
PROBABILITIES = np.linspace(0.01, 0.99, 99)


def cases() -> list[tuple[str, str, np.ndarray, tuple[tuple[str, float], ...]]]:
    normal = stats.norm.ppf(PROBABILITIES, loc=1.25, scale=1.8)
    normal_location, normal_scale = stats.norm.fit(normal)

    log_normal = stats.lognorm.ppf(PROBABILITIES, 0.7, loc=0, scale=np.exp(0.4))
    log_normal_shape, _, log_normal_scale = stats.lognorm.fit(log_normal, floc=0)

    gamma = stats.gamma.ppf(PROBABILITIES, 2.4, loc=0, scale=1 / 1.7)
    gamma_shape, _, gamma_scale = stats.gamma.fit(gamma, floc=0)

    beta = stats.beta.ppf(PROBABILITIES, 2.5, 4.5)
    beta_alpha, beta_beta, _, _ = stats.beta.fit(beta, floc=0, fscale=1)

    weibull = stats.weibull_min.ppf(PROBABILITIES, 1.6, loc=0, scale=2.3)
    weibull_shape, _, weibull_scale = stats.weibull_min.fit(weibull, floc=0)

    cauchy = stats.cauchy.ppf(PROBABILITIES, loc=-0.75, scale=1.4)
    cauchy_location, cauchy_scale = stats.cauchy.fit(cauchy)

    gumbel = stats.gumbel_r.ppf(PROBABILITIES, loc=0.8, scale=1.2)
    gumbel_location, gumbel_scale = stats.gumbel_r.fit(gumbel)

    inverse_gamma = stats.invgamma.ppf(PROBABILITIES, 4.2, loc=0, scale=2.5)
    inverse_gamma_shape, _, inverse_gamma_scale = stats.invgamma.fit(inverse_gamma, floc=0)

    logistic = stats.logistic.ppf(PROBABILITIES, loc=0.3, scale=0.9)
    logistic_location, logistic_scale = stats.logistic.fit(logistic)

    return [
        (
            "normal",
            "NormalRand",
            normal,
            (("location", normal_location), ("variance", normal_scale * normal_scale)),
        ),
        (
            "log_normal",
            "LogNormalRand",
            log_normal,
            (("location", np.log(log_normal_scale)), ("variance", log_normal_shape * log_normal_shape)),
        ),
        (
            "gamma",
            "GammaRand",
            gamma,
            (("shape", gamma_shape), ("rate", 1 / gamma_scale)),
        ),
        ("beta", "BetaRand", beta, (("alpha", beta_alpha), ("beta", beta_beta))),
        (
            "weibull",
            "WeibullRand",
            weibull,
            (("scale", weibull_scale), ("shape", weibull_shape)),
        ),
        (
            "cauchy",
            "CauchyRand",
            cauchy,
            (("location", cauchy_location), ("scale", cauchy_scale)),
        ),
        (
            "gumbel",
            "GumbelRand",
            gumbel,
            (("location", gumbel_location), ("scale", gumbel_scale)),
        ),
        (
            "inverse_gamma",
            "InverseGammaRand",
            inverse_gamma,
            (("shape", inverse_gamma_shape), ("rate", inverse_gamma_scale)),
        ),
        (
            "logistic",
            "LogisticRand",
            logistic,
            (("location", logistic_location), ("scale", logistic_scale)),
        ),
    ]


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    generated = cases()
    with (OUTPUT / "scipy_samples.csv").open("w", newline="") as output:
        writer = csv.writer(output, lineterminator="\n")
        writer.writerow(("case", "value"))
        for case, _, sample, _ in generated:
            writer.writerows((case, format(float(value), ".17g")) for value in sample)
    with (OUTPUT / "scipy_reference.csv").open("w", newline="") as output:
        writer = csv.writer(output, lineterminator="\n")
        writer.writerow(("case", "distribution", "parameter", "value"))
        for case, distribution, _, parameters in generated:
            writer.writerows(
                (case, distribution, parameter, format(float(value), ".17g"))
                for parameter, value in parameters
            )


if __name__ == "__main__":
    main()
