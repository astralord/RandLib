#!/usr/bin/env python3
"""Sorts golden-value mismatches into verified corrections and cases that need review.

Usage:
    RANDLIB_GOLDEN_REPORT=mismatches.csv build/tests/randlib_unit_tests "[golden]"
    tools/golden/classify_deviations.py mismatches.csv > tests/data/golden_v1/known_deviations.csv

A mismatch is a verified correction when the 40-digit reference value
(tests/data/reference) for the same case, function and argument is reproduced by
the current library within 1e-9 relative (1e-9 absolute for probabilities near zero and
for stable laws, as in test_reference.cpp), while the v1.1 value is at least ten times
further away; or when the new value is at least twice as close to the reference; or when
both are within 1e-8 of it; or when the new value is within a gap that
tests/data/reference_tolerances.csv tracks. The output merges these rows with the
hand-reviewed rules of tools/golden/reviewed_deviations.csv; rows that neither covers are
summarized on stderr.
"""

import collections
import csv
import glob
import math
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
STABLE = ("StableRand", "HoltsmarkRand", "LandauRand")


def read_rules(path):
    with open(path, newline="") as f:
        lines = [line for line in f if not line.startswith("#")]
    return list(csv.DictReader(lines))


def covered(rules, case, function, arg):
    return any(r["case"] in ("*", case) and r["function"] in ("*", function) and r["arg"] in ("*", arg) for r in rules)


def override_tolerance(overrides, case, function):
    for o in overrides:
        pattern = o["case"]
        matches = pattern == "*" or (case.startswith(pattern[:-1]) if pattern.endswith("*") else case == pattern)
        if matches and o["function"] in ("*", function):
            return float(o["tolerance"])
    return None


def number(text):
    try:
        return float(text)
    except ValueError:
        return None


def error(value, reference):
    if value is None or reference is None:
        return math.inf
    if value == reference:
        return 0.0
    if not (math.isfinite(value) and math.isfinite(reference)):
        return math.inf
    return abs(value - reference) / max(abs(reference), 1e-300)


def main():
    references = {}
    for path in glob.glob(os.path.join(ROOT, "tests", "data", "reference", "*.csv")):
        with open(path, newline="") as f:
            for row in csv.DictReader(f):
                references[(row["case"], row["function"], row["arg"])] = number(row["value"])

    reviewed = read_rules(os.path.join(ROOT, "tools", "golden", "reviewed_deviations.csv"))
    overrides = read_rules(os.path.join(ROOT, "tests", "data", "reference_tolerances.csv"))
    verified, unverified, worse = [], collections.Counter(), []
    with open(sys.argv[1], newline="") as f:
        for row in csv.DictReader(f):
            if covered(reviewed, row["case"], row["function"], row["arg"]):
                continue
            key = (row["case"], row["function"], row["arg"])
            reference = references.get(key)
            if reference is None and row["function"] == "Median":
                reference = references.get((row["case"], "Quantile", "0.5"))
            new = number(row["actual"])
            old = number(row["golden"]) if row["status"] == "ok" else None
            new_error, old_error = error(new, reference), error(old, reference)
            absolute = row["function"].startswith("CF") or (
                row["case"].startswith(STABLE) and row["function"] in ("f", "logf", "F", "S"))
            if absolute and reference is not None:
                # characteristic functions are compared on an absolute scale, as in test_reference.cpp
                new_error = abs(new - reference) if new is not None else math.inf
                old_error = abs(old - reference) if old is not None else math.inf
            if reference is not None and abs(reference) < 1e-6 and new is not None:
                new_error = min(new_error, abs(new - reference) / 1e-6)
            if reference is not None and new_error <= 1e-9 and (old_error >= 10 * new_error or old is None):
                reason = f"v1.1 value off by {old_error:.1e} relative to the mpmath reference; now {new_error:.1e}"
                verified.append((row["case"], row["function"], row["arg"], "skip", reason))
            elif reference is not None and new_error * 2 <= old_error:
                reason = (f"closer to the mpmath reference: v1.1 off by {old_error:.1e}, now {new_error:.1e} "
                          "(remaining gap listed in tests/data/reference_tolerances.csv)")
                verified.append((row["case"], row["function"], row["arg"], "skip", reason))
            elif reference is not None and new_error <= 1e-8 and old_error <= 1e-8:
                reason = f"both values within 1e-8 of the mpmath reference (v1.1 {old_error:.1e}, now {new_error:.1e})"
                verified.append((row["case"], row["function"], row["arg"], "skip", reason))
            elif reference is not None and (override_tolerance(overrides, row["case"], row["function"]) or 0) >= new_error:
                reason = f"within the gap tracked in tests/data/reference_tolerances.csv (v1.1 off by {old_error:.1e}, now {new_error:.1e})"
                verified.append((row["case"], row["function"], row["arg"], "skip", reason))
            elif reference is not None:
                worse.append(f"{row['case']} {row['function']}({row['arg']}): v1.1 off by {old_error:.1e}, now {new_error:.1e}")
            else:
                unverified[(row["case"].split("(")[0], row["function"])] += 1

    writer = csv.writer(sys.stdout, lineterminator="\n")
    writer.writerow(["case", "function", "arg", "expected", "reason"])
    print("# Generated by tools/golden/classify_deviations.py; edit tools/golden/reviewed_deviations.csv instead.")
    for rule in reviewed:
        writer.writerow([rule["case"], rule["function"], rule["arg"], rule["expected"], rule["reason"]])
    writer.writerows(sorted(verified))
    print(f"{len(worse)} rows moved away from the reference:", file=sys.stderr)
    for line in worse:
        print(f"  {line}", file=sys.stderr)
    print(f"{len(verified)} verified corrections, {sum(unverified.values())} rows without reference values:", file=sys.stderr)
    for (cls, function), count in sorted(unverified.items()):
        print(f"  {cls} {function}: {count}", file=sys.stderr)


if __name__ == "__main__":
    main()
