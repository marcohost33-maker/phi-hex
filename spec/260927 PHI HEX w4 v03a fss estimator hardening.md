# W4-v03a - Pre-data hardening of the correlation-ratio FSS estimator

> **Status:** AMENDMENT / VORREGISTRIERT 2026-09-27, before any PHY049
> production data were generated or inspected.
>
> This addendum supersedes only the estimator mechanics in W4-v03 Sec.5.2/Sec.5.3
> and the concrete G4 implementation. The observable, lattice ladder,
> temperature grid, seeds, sampler and claim ceiling of v03 remain unchanged.

## 1. Why an amendment is necessary

W4-v03 correctly selected the BKT scaling form

`R(T,L) = F(X),  X = L / exp(c / sqrt(T-T_BKT))`, for `T > T_BKT`,

but its phrase "cubic spline with a common smoothing parameter" did not pin
knots, smoothing selection, boundary handling or optimizer. That leaves enough
analysis freedom to overfit a visual/data collapse. Because no PHY049
production data exist yet, the estimator is hardened now rather than tuned
post hoc.

The primary literature basis remains the correlation-ratio FSS used by
Okabe/Otsuka. The implementation principle follows quantitative data-collapse
practice: turn collapse quality into a deterministic scalar objective instead
of judging a plot by eye.

## 2. Primary estimator: symmetric cross-size collapse score

For candidate `(T0,c)`, define

`X(L,T;T0,c) = L * exp(-c / sqrt(T-T0))`, with `T > T0`.

No free representation of `F` is fitted. For every unordered size pair
`La<Lb`, each observed point of one size is linearly interpolated on the other
size in `X`, and the operation is repeated in the reverse direction. A
comparison contributes

`q = (R_a - R_b,interp)^2 / (SE_a^2 + SE_b,interp^2)`.

Interpolation uncertainty is propagated from the two bracketing standard
errors with the same linear weights. A size pair is counted only if it has at
least 4 overlapping comparisons; rejected pairs contribute neither numerator
nor denominator. At least 60 accepted comparisons and at least 6 size pairs
are required.

The objective is the mean accepted mismatch `Q = sum(q)/N_comp`. Because the
interpolated comparisons are correlated and symmetric, `Q` is an engineering
collapse score, **not** a chi-square statistic and carries no p-value claim.

## 3. Fixed search domain and deterministic optimizer

Primary high-temperature fit window: `[0.585, 0.610]`.

The discriminating hypotheses fixed in v03 lie below 0.580. Therefore the
search is deliberately bounded to

- `T0 in [0.550, 0.582]`;
- `c in [0.15, 2.50]`.

A boundary hit is non-quotable and yields INCONCLUSIVE rather than extrapolating
outside the preregistered discrimination domain.

Search algorithm:

1. coarse `T0` grid spacing `0.001`;
2. 80 logarithmically spaced `c` values over the full c-domain;
3. local refinement within `T0 +/- 0.0025` at spacing `0.0002` and
   `c = 0.70..1.30` times the coarse optimum on 31 linear points, clipped to
   the fixed c-domain.

Tie handling follows deterministic iteration order; there is no random
optimizer state.

## 4. Identifiability gate

A low collapse score alone is insufficient: perfectly size-independent or
weakly structured curves can collapse for almost any `T0`.

For every coarse `T0`, profile over the same 80-point c-grid and record the
minimum `Q_profile(T0)`. Let `Q_profile,min` be the minimum on **that same
profile grid**. Acceptable profile points satisfy

`Q_profile(T0) <= Q_profile,min + 0.20`.

The span of acceptable `T0` values is the profile width. The primary fit is
quotable only when

- profile width `<= 0.012`;
- the optimum does not hit a T0 or c boundary; and
- refined collapse score `Q <= 2.50`.

The same-grid requirement is mandatory. Comparing a coarse profile against a
later refined global minimum can create an empty accepted profile and is a
known implementation error caught during G4 development.

## 5. G4 synthetic recovery - positive and adversarial cases

Before any production interpretation, PHY050 must pass all of the following:

### Positive recovery
Synthetic data are generated from one shared smooth collapse function with the
same ladder and T-grid as PHY049, independent Gaussian noise `SE=0.001`, and
known parameters:

- `(T0,c)=(0.565,0.90)`;
- `(0.573,0.90)`;
- `(0.579,0.75)`.

Each recovered T0 must be quotable and satisfy `|T_hat-T0| <= 0.003`.

### Negative controls

1. remove one required lattice size -> estimator must fail closed;
2. size-independent flat curves -> T0 must be non-identifiable/non-quotable;
3. deliberately size-dependent non-collapse curves -> T0 must be
   non-identifiable/non-quotable.

G4 validates estimator mechanics only. It is not evidence for a PHI-Hex
transition temperature.

## 6. Production uncertainty and robustness

The seed-level pair `(g(L/4), g(L/2))` is the bootstrap unit inside each
independent `(L,T)` chain family. Numerator and denominator are always resampled
together; the estimator remains a ratio of means, never a mean of ratios.
Independent `(L,T)` chains are resampled independently. The fixed original
sample delete-one-seed jackknife SE is used as the collapse weight in every
bootstrap replicate; no nested bootstrap is performed.

Primary bootstrap: 1000 replicates, RNG seed `50049`. At least 90% of
replicates must yield a quotable collapse fit, otherwise the uncertainty gate
fails closed.

Predeclared robustness variants:

- high-window trim `[0.585,0.6075]`;
- low edge raised to `[0.5875,0.610]`;
- leave each L out once.

Every variant must be quotable and the maximum `|T_variant-T_primary|` must be
`<=0.008`; otherwise G6 fails.

For the final v03 decision:

- `sigma_boot` = bootstrap standard deviation of quotable T0 replicates;
- `sigma_model` = half-span of quotable preregistered variants, floor `0.002`;
- `sigma_tot = hypot(sigma_boot, sigma_model)`;
- G5 requires `sigma_tot <= 0.010`.

## 7. Gate semantics after this amendment

- G0: complete production inputs and finite seed-level moments;
- G1/G2/G3: geometry, aligned limit and seed uniqueness from PHY049;
- G4: deterministic PHY050 synthetic recovery + adversarial nulls;
- G5: production power (`sigma_tot <= 0.010`);
- G6: production robustness (`max delta <= 0.008`).

A green G4 changes the state only from **PRE-FSS** to
**PREPRODUCTION_G4_VALIDATED**. It authorizes production measurement; it does
**not** enable physics interpretation. A physics result requires G0-G6 on real
PHY049 data and committed gate evidence in `results/`.

## 8. Literature anchors

- Y. Okabe and H. Otsuka, J. Phys. A: Math. Theor. 58 (2025) 065003,
  arXiv:2501.07388 - correlation ratio and BKT FSS variable.
- M. Tomita and Y. Okabe, Phys. Rev. B 66, 180401(R) (2002) - correlation-ratio
  finite-size scaling for KT/BKT systems.
- S. M. Bhattacharjee and F. Seno, J. Phys. A 34, 6375 (2001) - quantitative
  objective functions for data-collapse assessment.
- Y.-D. Hsieh, Y.-J. Kao and A. W. Sandvik, J. Stat. Mech. P09001 (2013),
  arXiv:1302.2900 - high-precision BKT FSS showing that subleading logarithmic
  corrections can materially shift extrapolated T_BKT values. This does not
  replace the correlation-ratio estimator, but it motivates the strict
  leave-one-L/window robustness gate and the claim ceiling.

## 9. Preproduction review hardening - 2026-09-27

Automated review before any PHY049 production data identified additional
fail-open paths. They are now part of the binding contract:

- **VAL-BIT:** if Numba is installed, a same-seed tiny Python/Numba trajectory
  must be bit-identical before the Numba backend is production-eligible.
- **G0 row identity:** every expected row must match exact `(L,t_idx,T,s,seed)`;
  duplicate seeds, wrong temperatures, wrong RNG seeds, non-finite moments,
  incomplete products, and wall-budget omissions fail closed.
- **Wall budget:** production is committed only in complete-L blocks. Once the
  24 h deadline has been reached, no further L is started and every skipped job
  is explicitly recorded as `WALL_BUDGET_STOP`. A started L is allowed to
  finish so a partial-L block cannot masquerade as complete evidence.
- **G1 geometry:** periodic translation by L in both primitive directions and
  same-sublattice preservation for L/4 and L/2 are checked directly; the
  aligned-state oracle alone is insufficient.
- **Splay direction:** for ordered L1<L2 the high-T sign is fixed before data as
  `R_L2-R_L1 < 0`; NaN/Inf, invalid errors, or the opposite sign break
  persistence.
- **Decision overlap:** H_A=[0.557,0.573] and H_B=[0.570,0.580] overlap on
  [0.570,0.573]. `SUPPORTED_LITERATURE` therefore requires the lower 95%
  bound to exceed 0.573; any 95% interval intersecting [0.570,0.573] is
  `OVERLAP`.
- **Single adjudicator:** PHY050 `assess_production()` is the only production
  path allowed to combine G0/G4/G5/G6, bootstrap uncertainty, model spread,
  and the internal discrimination label. It remains fail-closed when any
  prerequisite is non-quotable.
- **Exact production contract:** G0 requires n_seeds=12, n_therm=1000,
  n_meas=4000, max_workers=4, wall_budget_h=24.0, the matching persisted
  `campaign_contract`, exactly 5*29*12 raw rows, `complete is True`,
  `unmeasured == []`, exact row identity, finite wall_s >= 0, and finite
  correlation means constrained to [-1,1]. It also requires the exact
  persisted PHY049 preflight gate map with every VAL-BIT/G1-G4 verdict
  literally true.
  Extra rows, altered budgets, missing metadata and malformed moments fail
  closed.
- **Bootstrap contract:** `assess_production()` authorizes production
  adjudication only with exactly 1000 replicates. Reduced bootstrap counts may
  be used in isolated unit tests but must force the production decision to
  INCONCLUSIVE.
- **Runtime backend gate:** a worker/direct `_job()` path must itself consult
  VAL-BIT before selecting Numba; Numba availability alone is insufficient.
- **Production entrypoint gate:** public `produce()` must execute and require
  the complete VAL-BIT/G1-G4 preflight before submitting any measurement job,
  and persist those gate verdicts alongside the raw product.

These changes are preproduction hardening, not post-hoc tuning: no PHY049
production curve or transition estimate had been generated when they were
committed.
