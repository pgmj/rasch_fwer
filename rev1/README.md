# Revision 1 analyses (APM-26-08-243)

Analyses for the first revision of the manuscript. The review and the
response letter are in `paper_writing/Johansson2026fwer_rev1/`.

Render a file with `quarto render rev1/<file>.qmd` from the project folder.
`_metadata.yml` gives every file here a plain Typst PDF with the code shown.

## Files

| File | Status | Answers |
|---|---|---|
| `cached_analyses.qmd` | done | R1.1, R1.3, R1.4, R1.6, E3, E4, E5 |
| `rerun_plan.md`, `fresh_rerun.R`, `fresh_rerun_smoke.R` | tiers 1 and 2 run (13.2 h, no failures), tier 3 not run | E3, E4, R1.2 |
| `rerun_results.qmd` | done | E3, E4, R1.1, R1.2, R1.3 |

## cached_analyses.qmd

It runs no simulation. It reads the `.rds` caches in `simulations/`, and it
asserts that it reproduces Table 4, a Table 9 cell and the pooled Equation 5
values of the submission before reporting anything new.

Main results:

- **R1.1, strong control.** Under partial nulls at n = 1000, every rule exceeds
  5% familywise error over the fitting items (Westfall-Young 25 to 30% at the
  anchored level for phq9, with Holm and Bonferroni a few points lower). The
  uncorrected marginal p already flags 13.6 to 19.8% of the fitting phq9 items,
  so the per-item tests are miscalibrated before any correction is applied.
- **R1.3, FDR.** Benjamini-Hochberg exceeds its 5% FDR at n = 1000 (18 to 20% for
  phq9, up to 67% with strong dichotomous misfit) and is between 1.7 and 7.0% at
  n of 250 or below. The mechanism is the same as in R1.1.
- **E5, Holm.** Holm has the same reachability floor as Benjamini-Hochberg (0%
  detection at B = 100). Above the floor it is within about 2 points of
  Westfall-Young in detection. The advantage of Westfall-Young is the absence of
  a floor, not power.
- **E4, the decline with B.** Studentizing by moments of the simulated values
  makes the observed statistic t-distributed with B - 1 df. That predicts about
  6% at B = 100, which matches the fresh-data studies. Beyond that there is a
  mild conservatism (about 4.3% at B = 2000 pooled) that is present in fresh
  data too.
- **E3, master datasets.** The dichotomous null rate at n = 1000 is high at
  every B (5.6 to 9.2%), which the fresh-data polytomous studies do not show from
  B = 400. This is consistent with a master-dataset artefact, and the rerun will
  test it.
- **R1.4, disagreement.** Disagreement measured directly from seed pairs equals
  Equation 5 after a finite-seed correction when pooled (23, 10, 6, 4%). Per
  dataset, Equation 5 overstates disagreement for Westfall-Young in 44 of 480
  cases, by up to 37 points, which is the coupling through the step-down.

## rerun_results.qmd

- **E3.** The n = 1000 rise is gone in fresh data (Westfall-Young 3.8 to 4.2%
  from B = 400, against 5.6 to 9.2% from the master dataset). Convergence,
  the section 4.3 trade-off and the reachability sequences replicate.
- **E4 and R1.2.** Westfall-Young under the null: 5.6, 4.2, 3.8, 3.4, 3.5, 3.5,
  3.4% at B = 100 to 10 000. Joint studentization removes the whole decline
  (3.5% at B = 100, paired difference 2.05 points, SE 0.32), confirming the
  explanation. A flat conservatism of about 3.4% remains and does not depend
  on B, so it belongs to the parametric bootstrap, not to Monte Carlo error.
