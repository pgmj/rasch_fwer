# Fresh-data rerun: plan

For revision 1 of APM-26-08-243. Code in `fresh_rerun.R`, smoke test in
`fresh_rerun_smoke.R` (log and `.rds` beside them). Status: tiers 1 and 2
run (13.2 h, no failures), results in `rerun_results.pdf`.
Tier 3 not run.

## What it answers

- **E3.** Studies 1 and 2 drew every sample from one 10 000-person master
  dataset per condition. The editor names both studies. In the cached results
  the Westfall-Young null rate at n = 1000 is 5.6 to 9.2% at every B, against
  2.4 to 6.8% at smaller n, a rise the fresh-data polytomous studies do not show
  from B = 400 (see `cached_analyses.pdf`).
- **E4.** Table 4 stops at B = 2000. The editor asks for 10 000, and for the
  decline with B to be discussed.
- **R1.2.** The explanation for the decline (studentizing by the moments of
  the simulated values alone) is tested directly, by a second set of p-values
  studentized jointly over the observed and simulated statistics.

## Design

Item parameters are those of the submission: the same 20 dichotomous
locations from seed 20250727, the same three misfitting items, person SD 1.5,
rho = .15 for misfit. Only the respondents change, drawn fresh per replication.

| Tier | Data | k misfitting | n | B (nested) | Replications | Replaces |
|---|---|---|---|---|---|---|
| 1 | 20 dichotomous | 0 | 150, 250, 500, 1000 | 100, 200, 400, 1000, 2000, 5000, 10 000 | 500 per n | Study 1 (Fig 1, Table 7), Study 2 null (Fig 2 dichotomous panel, Table 4, Table 8 rows) |
| 2 | 20 dichotomous | 1, 3 | 150, 250, 500, 1000 | 100, 200, 400, 1000, 2000 | 250 per cell | Study 2 misfit (section 4.3, the BH and WY sequences in 5.1) |
| 3, optional | phq9, sim15 polytomous | 0 | as in Study 3 | as tier 1 | 250 per cell | Adds B = 5000 and 10 000 to the polytomous null |

**Nested B.** Each replication runs one bootstrap at the largest B, and every
smaller B uses its first B draws. The first B draws of an iid sequence are an
iid sample of size B, so the rate at each B has the same distribution as with a
separate bootstrap, and comparisons across B become paired within the
bootstrap as well as within the dataset. It also makes B = 10 000 affordable: a
replication costs 10 000 refits instead of the 18 700 that separate bootstraps
at all seven values would cost. The sentence in the Method that says every B
uses an independent bootstrap must change for the rerun studies.

**Output per replication and B.** Observed infit, marginal p, Westfall-Young
and Benjamini-Hochberg adjusted p (easyRasch2's own routine), marginal p and
Westfall-Young adjusted p with joint studentization, and the interval bounds
and flag at the five widths. Holm and Bonferroni are computed afterwards from
the marginal p. Run time in seconds is stored, so cost can be audited.

**Precision.** Tier 1 gives 2 000 replications per B pooled over n, a binomial
standard error of 0.49 points at 5%, and 0.97 per sample size. The paired
contrast between B = 100 and large B is more precise than either rate. The
cached results put the small-B excess at about 1 to 1.8 points and the
large-B shortfall at about 0.7, so both are resolved at the pooled level.
Raising tier 1 to 1 000 per n halves the variance at twice the cost (see the
projection).

**Seeds.** Deterministic per cell and disjoint from the seeds of the original
studies (data 7e8 + ..., bootstrap 8e8 + ...). Grid order places one cell of
each condition on every core in each chunk (`block_order()`), and chunks are
blocks of ten replications, so an interrupted run stays balanced.

## Checks done in the smoke test

- **Reproduction under easyRasch2 1.4.0.** One cached cell of the submission
  (rep 1, n = 150, k = 0) was rebuilt from its master dataset and seeds. Its
  marginal, Westfall-Young and Benjamini-Hochberg p-values and its interval
  flags at B = 100 and 400 are identical to the cached values. The paper's
  results were produced under 1.1.0 and 1.1.1, so this shows that nothing the
  rerun uses has changed since.
- **Nested evaluation.** `eval_nested()` reproduces `RMitemInfit()` on the same
  bootstrap.
- **Generator.** With n = 1000 and three misfitting items, mean infit was 1.39
  for the misfitting items and 0.93 for the fitting ones, as expected at
  rho = .15.

## Second-round objections considered

| A reviewer could say | How the design answers it |
|---|---|
| Only the null was rerun, the misfit results still rest on masters | Tier 2 reruns every misfit condition of Study 2. |
| The convergence study also used masters, with only 30 replications | Tier 1 stores the interval bounds at every B, which replaces Study 1 with 500 fresh replications per n. The reference is B = 10 000 of the same bootstrap rather than a separate one, which is stated. |
| The rerun used a newer package version | The smoke test reproduces a cached cell exactly under 1.4.0. |
| Nested B conditions are not independent | The rate at each B is unaffected (first B of iid draws). The original design was already paired on the dataset. Pairing makes the B contrasts more precise, which is what E4 asks about. |
| The studentization explanation is post hoc | The joint-studentization p-values are computed on the same draws. If the B = 100 excess disappears there, the explanation is confirmed. If not, it is refuted and the letter changes. Both outcomes are reportable. |
| Is the conservatism at large B real? | 2 000 replications per B up to 10 000 iterations. Tier 3 checks it in polytomous data. The cause (the composite null) is described, not identified. |
| The n = 1000 rise disappears, so was the submission wrong? | Report both, explain the master-dataset dependence, and state which conclusions move. Pooled conclusions were not contradicted in the cached analyses. |
| Strong control in fresh data | Tier 2 gives partial-null familywise error and FDR from fresh data. |
| Why stop at 10 000? | The editor asked for it. At w = .999 it still gives only five draws per tail, which the text already says. |

**Deliberately out of scope:** the conditional DGP (experimental in the
package, would test a different null), other alpha levels, more item counts,
and the iterative removal workflow. The Limitations already name the last
three.

## Run-time projection

From the smoke test, one replication per condition, run on ten cores at once.

Seconds per replication, measured with 13 replications sharing ten cores,
which is close to the load of the real run:

| Tier | Condition | n = 150 | n = 250 | n = 500 | n = 1000 |
|---|---|---|---|---|---|
| 1 | dichotomous null, B to 10 000 | 130 | 151 | 201 | 286 |
| 2 | dichotomous, k = 3, B to 2000 | 26 | 31 | 40 | 58 |
| 3 | phq9 null, B to 10 000 | 189 | 241 | | 603 |
| 3 | sim15 null, B to 10 000 | | 510 | | 1238 |

Projected wall time on ten cores (sum of core-seconds divided by ten; the grid
order keeps the cores evenly loaded, plus about 3% for checkpoints). Tier 2
assumes k = 1 costs the same as k = 3, which it should, since cost is set by
n, items and B.

| Run | Replications | Hours |
|---|---|---|
| Tier 1 | 500 per n | 11 |
| Tier 1 | 1 000 per n | 22 |
| Tier 2 | 250 per cell | 2.2 |
| Tier 3 | 250 per cell | 20 |
| Tier 3 | 120 per cell | 9.5 |
| **Tiers 1 and 2 (recommended)** | 500 / 250 | **about 13** |
| Tiers 1, 2 and 3 | 500 / 250 / 120 | about 23 |

The smoke test itself took 20.6 minutes of wall time. Stored results are
about 83 kB per tier-1 replication in memory, so the tier-1 cache will be
roughly 170 MB in memory and much less on disk.

Tier 3 is optional. Its cost is dominated by the 15-item set at n = 1000. It
answers a question nobody has asked yet (whether the large-B conservatism holds
for polytomous items). I would run tiers 1 and 2 first and decide on tier 3 from
their results.

## How to run

From `rev1/`, in a terminal:

```bash
nohup caffeinate -i Rscript -e 'source("fresh_rerun.R")' > fresh_rerun.log 2>&1 &
```

```bash
tail -f fresh_rerun.log
```

Tier 3 runs only with `FRESH_TIER3=1` set. The run is resumable: rerunning
skips cells already in the `.rds` files, and raising `R_NULL` or `R_MIS`
computes only the new replications.
