# Numerical accuracy and evidence policy

RandLib separates distributions into two evidence tiers. Tier 1 distributions have full probability,
sampling, and fitting expectations; Tier 2 distributions are tested for coherence and sampling and
use documented best effort where no robust reference algorithm is available.

| Quantity | Tier 1 target | Tier 2 target |
|---|---:|---:|
| PDF/PMF, CDF, survival in the bulk | `1e-12` relative | `1e-8` relative |
| Tails down to probability `1e-300` | `1e-9` relative | Best effort, documented |
| Quantile round trip | `1e-10` relative | `1e-7` relative |
| Closed-form moments | Correct to floating-point rounding | Correct to floating-point rounding |
| Characteristic function | `1e-8` absolute | Best effort, documented |
| Stable PDF/CDF | `1e-9` absolute | Not applicable |
| Sampling goodness of fit | Required | Required |
| Fitting | Required | Optional |

“Bulk” means `1e-6 < F(x) < 1-1e-6`. A target is not a blanket guarantee for every parameter and
platform: committed exceptions in `tests/data/reference_tolerances.csv` name their region, tolerance,
and tracking issue.

## What counts as evidence

Probability values are generated with mpmath at 40 decimal digits. Exact identities and closed-form
special cases supplement those tables. Another double-precision library is useful for investigation,
but it is not accepted as the sole reference for changing a golden value.

Every output difference from the v1.1 golden corpus must have one of:

- a committed mpmath reference;
- an exact formula recorded in the deviation reason; or
- a focused regression test establishing the corrected invariant.

`tools/refdata/generate.py` regenerates high-precision reference CSVs. The procedure for reviewing and
classifying golden deviations is in `NEXT_STEPS.md`. Reference generators and their dependency
versions are part of the evidence and are committed with the data.

## Test layers

The suite combines exact regression tests, high-precision point references, PDF/CDF/survival and
quantile coherence, parameter validation, deterministic engine advancement, seeded parameter
recovery, and statistical sampling tests. Statistical checks use fixed seeds and predeclared
thresholds; a single attractive histogram is not validation.

Survival functions should be implemented directly when subtraction would lose tail precision.
Tolerance relaxations require a nearby explanation and an issue link. Improvements must first remove
or tighten the corresponding registry override, so the measured contract and documentation stay in
sync.
