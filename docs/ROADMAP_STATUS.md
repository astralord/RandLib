# Roadmap status

Status of each task in [MODERNIZATION_PLAN.md](MODERNIZATION_PLAN.md). This file stands in for
the GitHub milestones and issues of plan task 0.5 until they are created.

Legend: ✅ done · 🟡 partly done (see note) · ⏳ not started · 🚫 needs the maintainer

The plan was revised on 2026-10-11 to do the hard numerics last; the `2.x` branch has the current
plan and status.

## Phase 0: Preparation

| Task | Status | Notes |
|---|---|---|
| 0.1 Coordinate upstream | 🚫 | Tracking issue on `astralord/RandLib` and agreement on the branch model are for the maintainer. |
| 0.2 Tag the baseline | ✅ | `v1.1.0-baseline` at 6536dad. |
| 0.3 Golden outputs | ✅ | `tools/golden/` records 14.6k values for 160 cases from the baseline library; calls that hang in the baseline are recorded as `timeout`. |
| 0.4 Reference data | ✅ | `tools/refdata/` (mpmath, 40 digits). Stable laws use CF inversion, convergent or asymptotic series and rotated contours, never SciPy; every method is cross-checked. Geometric stable is not covered yet. |
| 0.5 Project tracking | 🚫 | Milestones and issues need repository access; this file tracks progress meanwhile. |

## Phase 1: Stabilize (v1.2)

| Task | Status | Notes |
|---|---|---|
| 1.1 Repository hygiene | ✅ | |
| 1.2 CMake build | ✅ | Presets, install/export, `find_package`, consumer example built in CI. |
| 1.3 Compiler portability (#18) | 🟡 | GCC 11–14, Clang 15–18 and AppleClang build warning-free in CI; the MSVC 14.44 duplicate-instantiation fix awaits CI. |
| 1.4 Correctness fixes B1–B10 | ✅ | B10 is deferred to v2 by design. 27 further bugs (B11–B37) found by the test suite and fixed; see CHANGELOG.md. |
| 1.5 Continuous integration | 🟡 | GCC, Clang, sanitizer, coverage and static-analysis jobs pass; the MSVC 14.44 build fix and final AppleClang statistical result await CI. |
| 1.6 Test suite (#12) | ✅ | Golden, reference, coherence, mode/median, reproducibility, engines, parameters, regressions, sampling (KS/χ², moments, empirical CF) and fitting. |
| 1.7 Release v1.2 | 🟡 | CHANGELOG written; `v1.2.0` is tagged once CI is green on all compilers. |

## Phase 2: Modern architecture (v2.0-alpha)

| Task | Status | Notes |
|---|---|---|
| 2.1 Header-only layout and namespace | ⏳ | |
| 2.2 URBG interoperability | ⏳ | |
| 2.3 Remove hidden state | 🟡 | The static mutable state found in v1.1 is gone in v1.2 (B1, B11, B12); the CI grep check comes with 2.1. |
| 2.4 API cleanup | ⏳ | |
| 2.5 Numerical kernel (#10) | 🟡 | Adaptive Gauss–Kronrod replaced adaptive Simpson in v1.2; quantiles use a bracketing root finder. Tanh-sinh, oscillatory quadrature, TOMS 748 and optimizers remain. |
| 2.6 Tooling | ⏳ | |

## Phase 3: The niche features (v2.0)

| Task | Status | Notes |
|---|---|---|
| 3.1 Fitting API | ⏳ | |
| 3.2 Heavy-tail toolkit | ⏳ | |
| 3.3 Stable distributions (#1) | 🟡 | v1.2 fixes the CF sign, the α = 1 location, the B8 floor and the tails; S0 and the near-α = 1 region remain. |
| 3.4 Noncentral t (#7) | ⏳ | |
| 3.5 Marcum Q (#4) | 🟡 | v1.2 fixes the noncentral χ² density (Bessel bug); with the new incomplete gamma and quadrature the distribution function meets 1e-12 in the bulk, but the upper tail is still 1 − Marcum P. |
| 3.6 Log-normal CF (#6) | 🟡 | With Gauss–Kronrod the CF meets the 1e-8 target for σ² ≤ 1; for σ² = 4 the error is 6e-5 (was 4e-3). |
| 3.7 Log-stable (#5) | ⏳ | |
| 3.8 Close #11 | 🚫 | Closing the issue needs the maintainer. |

## Phase 4: Documentation, packaging, release

| Task | Status | Notes |
|---|---|---|
| 4.1 Documentation | ⏳ | |
| 4.2 Migration guide | ⏳ | |
| 4.3 Packaging | ⏳ | |
| 4.4 Release v2.0.0 | 🚫 | Release candidates, announcements and the upstream vcpkg PR need the maintainer. |
