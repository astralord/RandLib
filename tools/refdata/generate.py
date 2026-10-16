#!/usr/bin/env python3
"""Writes tests/data/reference/<Class>.csv with 40-digit mpmath values for Tier 1 cases.

Usage: tools/refdata/generate.py [--jobs N] [--only SUBSTRING]

Evaluation points come from tests/data/golden_v1 (the x grid recorded for the case),
plus tail points where F or S equals 1e-10, 1e-30, 1e-100 and 1e-300. Quantiles use the
p grid of the golden dump. Every x and p is the exact double stored in the CSV.
This is a developer tool; CI only reads its committed output.
"""

import argparse
import csv
import glob
import math
import multiprocessing
import os
import re
import sys

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import distributions  # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
GOLDEN = os.path.join(ROOT, "tests", "data", "golden_v1")
OUT = os.path.join(ROOT, "tests", "data", "reference")
DPS = 40
P_GRID = [1e-10, 1e-6, 0.001, 0.05, 0.3, 0.5, 0.7, 0.95, 0.999, 1 - 1e-6]
T_GRID = [0.1, 0.5, 1, 2, 5]
TAIL_LEVELS = [1e-10, 1e-30, 1e-100, 1e-300]


def fmt_arg(x):
    return "%.17g" % x


def fmt_value(v):
    if isinstance(v, mp.mpc):
        raise TypeError("complex value")
    if mp.isnan(v):
        return "nan"
    if mp.isinf(v):
        return "inf" if v > 0 else "-inf"
    return mp.nstr(v, 25, min_fixed=-100, max_fixed=100) if v != 0 else "0"


def parse_case(case_id):
    cls, args = re.match(r"^(\w+)\((.*)\)$", case_id).groups()
    values = [float(a) for a in args.split(",")] if args.strip() else []
    return cls, values


def golden_cases():
    cases = {}
    for path in sorted(glob.glob(os.path.join(GOLDEN, "*.csv"))):
        if os.path.basename(path) == "known_deviations.csv":
            continue
        with open(path, newline="") as f:
            for row in csv.DictReader(f):
                entry = cases.setdefault(row["case"], {"x": []})
                if row["function"] in ("f", "P"):
                    entry["x"].append(row["arg"])
    cases.update({
        "GeneralizedParetoRand(0.2, 2, 1)": {"x": ["1", "1.1", "2", "5", "20"]},
        "GeneralizedParetoRand(0, 1.5, -1)": {"x": ["-1", "-0.9", "0", "2", "10"]},
        "GeneralizedParetoRand(-0.3, 2, 0)": {"x": ["0", "0.1", "1", "4", "6.5"]},
        "GeneralizedExtremeValueRand(0.2, 1.5, -0.5)": {"x": ["-7", "-2", "-0.5", "2", "10"]},
        "GeneralizedExtremeValueRand(0, 2, 1)": {"x": ["-10", "-2", "1", "5", "15"]},
        "GeneralizedExtremeValueRand(-0.25, 1, 0)": {"x": ["-10", "-2", "0", "2", "3.9"]},
        "RiceRand(0, 1)": {"x": ["0", "0.1", "0.5", "1", "2", "5"]},
        "RiceRand(2, 0.5)": {"x": ["0", "0.5", "1", "2", "3", "5"]},
        "RiceRand(5, 2)": {"x": ["0", "1", "3", "5", "8", "15"]},
    })
    return cases


class Mapped:
    """Monotone map u -> x covering the support, so bracketing works on the whole real line."""

    def __init__(self, dist):
        self.lo, self.hi = mp.mpf(dist.lower), mp.mpf(dist.upper)

    def x(self, u):
        lo_finite, hi_finite = mp.isfinite(self.lo), mp.isfinite(self.hi)
        if lo_finite and hi_finite:
            return self.lo + (self.hi - self.lo) / (1 + mp.exp(-u))
        if lo_finite:
            return self.lo + mp.exp(u)
        if hi_finite:
            return self.hi - mp.exp(-u)
        return mp.sinh(u)


def continuous_quantile(dist, level, upper_tail):
    """x with F(x) = level (or S(x) = level when upper_tail), solved on a log scale."""
    level = mp.mpf(level)
    m = Mapped(dist)
    if upper_tail:
        def g(u):
            s = dist.sf(m.x(u))
            return mp.log(level) - mp.log(s) if s > 0 else mp.inf
    else:
        def g(u):
            p = dist.cdf(m.x(u))
            return mp.log(p) - mp.log(level) if p > 0 else -mp.inf
    a, b = mp.mpf(-1), mp.mpf(1)
    ga, gb = g(a), g(b)
    while ga > 0:
        b, gb = a, ga
        a *= 2
        if a < -1e6:
            return None
        ga = g(a)
    while gb < 0:
        a, ga = b, gb
        b *= 2
        if b > 1e6:
            return None
        gb = g(b)
    # Anderson-Bjorck regula falsi: superlinear and keeps the bracket, so it needs
    # far fewer evaluations than bisection when the distribution function is costly.
    tolerance = mp.mpf(10) ** (-DPS + 8)
    for _ in range(200):
        c = b - gb * (b - a) / (gb - ga) if gb != ga else (a + b) / 2
        if not (a < c < b):
            c = (a + b) / 2
        gc = g(c)
        if gc == 0:
            return m.x(c)
        if (gc > 0) == (gb > 0):
            factor = 1 - gc / gb
            ga = ga * (factor if factor > 0 else mp.mpf(1) / 2)
            b, gb = c, gc
        else:
            a, ga = b, gb
            b, gb = c, gc
        if abs(b - a) < tolerance * max(1, abs(b)):
            break
    return m.x((a + b) / 2)


def discrete_quantile(dist, level, upper_tail):
    """Smallest k with F(k) >= level, or with S(k) <= level when upper_tail."""
    level = mp.mpf(level)

    def ok(k):
        return dist.sf(k) <= level if upper_tail else dist.cdf(k) >= level

    lo = int(dist.lower) - 1 if mp.isfinite(dist.lower) else -1
    while ok(lo):
        lo = 2 * lo - 1
    hi = max(lo + 1, 1)
    while not ok(hi):
        hi *= 2
        if hi > 2**62:
            return None
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if ok(mid):
            hi = mid
        else:
            lo = mid
    return hi


def _confirm_tiny(fn, value):
    """Recompute a tiny value at higher precision; cancellation shows up as disagreement."""
    for extra in (60, 150):
        with mp.workdps(DPS + extra):
            again = fn()
        if abs(again - value) <= abs(again) * mp.mpf(10) ** (-DPS + 10):
            return value
        value = again
    raise ArithmeticError(f"value {mp.nstr(value, 5)} is not stable under higher precision")


def compute_case(item):
    case_id, x_args = item
    mp.mp.dps = DPS
    cls, params = parse_case(case_id)
    dist = distributions.make(cls, params)
    rows = []

    def add(function, arg, fn):
        try:
            value = fn()
            if isinstance(value, mp.mpf) and value != 0 and mp.isfinite(value) and abs(value) < mp.mpf("1e-25"):
                value = _confirm_tiny(fn, value)
            rows.append((case_id, function, arg, fmt_value(value)))
        except Exception as error:  # noqa: BLE001 - a missing value is reported, not fatal
            print(f"{case_id} {function}({arg}) failed: {error}", file=sys.stderr)

    light = getattr(dist, "light", False)
    if light:
        mp.mp.dps = 20
        x_args = sorted(x_args, key=float)[1::3]
    points = {float(x) for x in x_args}
    is_stable = cls in ("StableRand", "HoltsmarkRand", "LandauRand")
    if not dist.discrete and not is_stable and not light:
        for level in TAIL_LEVELS:
            for upper in (False, True):
                x = continuous_quantile(dist, level, upper)
                if x is not None and mp.isfinite(x) and abs(x) < 1e308:
                    xd = float(x)
                    if dist.lower < xd < dist.upper and 0 < (dist.sf(xd) if upper else dist.cdf(xd)) < 0.5:
                        points.add(xd)
    if dist.discrete:
        for level in (1e-30, 1e-100):
            k = discrete_quantile(dist, level, True)
            if k is not None and k < 2**31 - 1:
                points.add(float(k))

    for x in sorted(points):
        arg = fmt_arg(x) if not dist.discrete else str(int(x))
        xm = mp.mpf(x)
        if dist.discrete:
            add("P", arg, lambda: dist.pmf(int(x)))
            add("logP", arg, lambda: dist.logpmf(int(x)))
            add("F", arg, lambda: dist.cdf(int(x)))
            add("S", arg, lambda: dist.sf(int(x)))
        else:
            add("f", arg, lambda: dist.pdf(xm))
            add("logf", arg, lambda: dist.logpdf(xm))
            add("F", arg, lambda: dist.cdf(xm))
            add("S", arg, lambda: dist.sf(xm))

    quantile = discrete_quantile if dist.discrete else continuous_quantile
    for p in ([] if light else P_GRID):
        pm = mp.mpf(p)
        for name, upper in (("Quantile", False), ("Quantile1m", True)):
            if dist.discrete:
                k = quantile(dist, pm, upper)
                levels = [dist.sf(j) if upper else dist.cdf(j) for j in (k - 1, k)]
                if any(abs(level - pm) <= pm * mp.mpf("1e-12") for level in levels):
                    continue  # p equals F or S at an atom up to rounding: two answers are correct
            add(name, fmt_arg(p), lambda upper=upper: mp.mpf(quantile(dist, pm, upper)))

    for name, value in dist.moments().items():
        add(name, "", lambda: mp.mpf(value))
    if not dist.discrete and "Median" not in dist.moments() and not is_stable and not light:
        add("Median", "", lambda: continuous_quantile(dist, mp.mpf(1) / 2, False))

    if hasattr(dist, "cf"):
        for t in T_GRID:
            value = None

            def cf():
                nonlocal value
                if value is None:
                    value = dist.cf(mp.mpf(t))
                return value

            add("CF_re", fmt_arg(t), lambda: cf().real)
            add("CF_im", fmt_arg(t), lambda: cf().imag)
    return case_id, rows


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--jobs", type=int, default=os.cpu_count())
    parser.add_argument("--only", default="", help="regenerate only cases whose id contains this text")
    parser.add_argument("--cache", default=os.path.join(ROOT, "build", "refdata-cache"),
                        help="directory for per-case results, reused by later runs (delete it to recompute)")
    parser.add_argument("--cached-only", action="store_true", help="write the CSV files from cached cases only")
    args = parser.parse_args()

    cases = golden_cases()
    work = [(cid, entry["x"]) for cid, entry in cases.items() if args.only in cid and _supported(cid)]
    os.makedirs(OUT, exist_ok=True)
    results = {}
    os.makedirs(args.cache, exist_ok=True)
    todo = []
    for case_id, xs in work:
        cached = _cache_path(args.cache, case_id)
        if os.path.exists(cached):
            with open(cached, newline="") as f:
                results[case_id] = [tuple(r) for r in csv.reader(f)]
        else:
            todo.append((case_id, xs))
    if todo and args.cached_only:
        print(f"skipping {len(todo)} uncached cases", file=sys.stderr)
        todo = []
    if todo:
        with multiprocessing.Pool(min(args.jobs, len(todo))) as pool:
            for case_id, rows in pool.imap_unordered(compute_case, sorted(todo, key=lambda w: -_cost(w[0]))):
                results[case_id] = rows
                with open(_cache_path(args.cache, case_id), "w", newline="") as f:
                    csv.writer(f, lineterminator="\n").writerows(rows)
                print(f"done {case_id}: {len(rows)} values", file=sys.stderr, flush=True)

    by_class = {}
    for case_id, rows in results.items():
        by_class.setdefault(parse_case(case_id)[0], []).extend(rows)
    for cls, rows in by_class.items():
        path = os.path.join(OUT, f"{cls}.csv")
        existing = []
        if args.only and os.path.exists(path):
            regenerated = {r[0] for r in rows}
            with open(path, newline="") as f:
                existing = [tuple(r.values()) for r in csv.DictReader(f) if r["case"] not in regenerated]
        rows = sorted(existing + rows, key=lambda r: (r[0], r[1], _arg_key(r[2])))
        with open(path, "w", newline="") as f:
            writer = csv.writer(f, lineterminator="\n")
            writer.writerow(["case", "function", "arg", "value"])
            writer.writerows(rows)


def _supported(case_id):
    try:
        cls, params = parse_case(case_id)
        return distributions.make(cls, params) is not None
    except (TypeError, ValueError):
        return False


def _cache_path(cache, case_id):
    return os.path.join(cache, re.sub(r"[^A-Za-z0-9.,-]+", "_", case_id) + ".csv")


def _cost(case_id):
    return 10 if case_id.startswith(("StableRand", "HoltsmarkRand", "LandauRand", "GeometricStableRand")) else 1


def _arg_key(arg):
    try:
        return (0, float(arg))
    except ValueError:
        return (1, arg)


if __name__ == "__main__":
    main()
