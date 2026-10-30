# Smoke test for fresh_rerun.R. Run from rev1/:
#   Rscript fresh_rerun_smoke.R > fresh_rerun_smoke.log 2>&1
# Writes fresh_rerun_smoke.rds (checks and per-replication timings).
#
# Version note. This script passed under easyRasch2 1.4.0 before the full run
# (fresh_rerun_smoke.log, "A/B passed: 1.4.0 reproduces the cached cell"). It
# does not apply from 1.4.0.9000. There RMitemInfit() studentizes by the joint
# moments of observed and simulated values, so checks A (reproduce the cached
# submission cell) and B (nested evaluation matches RMitemInfit()) fail by
# design. Meanwhile fresh_rerun.R, which this script sources, now requires
# 1.4.0.9000 for its moments = "simulated" pin. The script is kept as the
# record of the 1.4.0 check and is not meant to be rerun.

Sys.setenv(FRESH_SMOKE = "1")
source("fresh_rerun.R")
out <- list(version = as.character(packageVersion("easyRasch2")))

# ---- A. easyRasch2 1.4.0 reproduces a cached cell of the submission -----------
# Cell rep = 1, n = 150, k = 0 of simulations/infit_width_iterations_S2.rds,
# rebuilt with the code of that study (master dataset, seeds), at B = 100, 400.
s2 <- readRDS("../simulations/infit_width_iterations_S2.rds")
ci <- which(s2$grid$rep == 1 & s2$grid$n == 150 & s2$grid$k == 0)
stopifnot(length(ci) == 1)
cached <- s2$res[[ci]]

master0 <- local({
  seed <- 20250727
  set.seed(seed)
  locs <- runif(20, -2, 2)
  W <- cbind(rep(1, 20), rep(0, 20))
  Sigma <- matrix(c(1.5^2, .15 * 1.5^2, .15 * 1.5^2, 1.5^2), 2, 2)
  set.seed(seed + 0)
  th <- MASS::mvrnorm(10000, mu = c(0, 0), Sigma = Sigma)
  d <- as.data.frame(eRm::sim.xdim(persons = th, items = locs,
                                   Sigma = Sigma, weightmat = W))
  names(d) <- paste0("V", 1:20)
  d
})
stopifnot(all.equal(dich_items$locs, local({ set.seed(20250727); runif(20, -2, 2) })))

set.seed(6e7 + 1e6 * 0 + 1e3 * 150 + 1)
d <- master0[sample(10000, 150), ]
obs <- RMitemInfit(d, output = "dataframe")
msq <- obs$Infit_MSQ[match(names(d), obs$Item)]
B2 <- c(100, 200, 400, 1000, 2000)
repro <- lapply(c(1, 3), function(j) {
  co <- RMitemInfitCutoff(d, iterations = B2[j], parallel = FALSE,
                          seed = 1 * 100 + j, verbose = FALSE)
  fw <- RMitemInfit(d, cutoff = co, p_value = TRUE, correction = "fwer",
                    alpha = ALPHA, output = "dataframe")
  bh <- RMitemInfit(d, cutoff = co, p_value = TRUE, correction = "fdr_bh",
                    alpha = ALPHA, output = "dataframe")
  i1 <- match(names(d), fw$Item)
  now <- data.frame(p_marg = fw$p_infit[i1], padj_fwer = fw$padj_infit[i1],
                    padj_bh = bh$padj_infit[match(names(d), bh$Item)])
  old <- cached$pval[cached$pval$B == B2[j], c("p_marg", "padj_fwer", "padj_bh")]
  codes_now <- unlist(lapply(WIDTHS, function(w) {
    b <- bounds_at(co$results, w, names(d)); flag_code(msq, b[1, ], b[2, ])
  }))
  codes_old <- cached$hdci$code[cached$hdci$B == B2[j]]
  # B. the nested evaluation reproduces RMitemInfit() on the same bootstrap
  nest <- eval_nested(d, co, B2[j])$pval
  list(B = B2[j],
       pval_identical = isTRUE(all.equal(as.matrix(now), as.matrix(old),
                                         check.attributes = FALSE)),
       hdci_identical = identical(as.integer(codes_now), as.integer(codes_old)),
       nested_matches = isTRUE(all.equal(nest[, c("p_marg", "padj_fwer", "padj_bh")],
                                         now, check.attributes = FALSE)),
       pooled_vs_standard_maxdiff = max(abs(nest$padj_fwer_pool - nest$padj_fwer)))
})
out$repro <- repro
print(do.call(rbind, lapply(repro, as.data.frame)))
stopifnot(all(vapply(repro, function(r) r$pval_identical && r$hdci_identical &&
                       r$nested_matches, logical(1))))
message("A/B passed: 1.4.0 reproduces the cached cell, nested evaluation matches")

# ---- C. fresh generator sanity ------------------------------------------------
d3 <- sim_dich(1000, 3, seed = 1)
f3 <- RMitemInfit(d3, output = "dataframe")
out$gen <- data.frame(item = f3$Item, msq = f3$Infit_MSQ,
                      misfit = f3$Item %in% paste0("V", dich_items$mis))
print(aggregate(msq ~ misfit, out$gen, mean))
dp <- sim_poly_null("phq9", 600, seed = 1)
stopifnot(ncol(dp) == 9, all(unlist(dp) %in% 0:3))

# ---- D. timing: one replication per condition, run as in the real job ---------
jobs_t1 <- lapply(NS, function(n)
  list(tier = 1, g = data.frame(rep = 9999, k = 0L, n = n)))
jobs_t2 <- lapply(NS, function(n)
  list(tier = 2, g = data.frame(rep = 9999, k = 3L, n = n)))
jobs_t3 <- list()
for (a in names(poly_arms)) for (n in poly_arms[[a]]$ns)
  jobs_t3[[length(jobs_t3) + 1L]] <-
    list(tier = 3, g = data.frame(arm = a, rep = 9999, n = n,
                                  stringsAsFactors = FALSE))
jobs <- c(jobs_t1, jobs_t2, jobs_t3)
t0 <- Sys.time()
tim <- parallel::mclapply(jobs, function(j) {
  r <- switch(j$tier,
              one_rep_dich(j$g, B_NULL),
              one_rep_dich(j$g, B_MIS),
              one_rep_poly(j$g, B_NULL))
  list(job = j, ok = !is.null(r), sec = if (is.null(r)) NA else r$sec,
       actual_B = if (is.null(r)) NA else r$actual_B,
       size = as.numeric(object.size(r)))
}, mc.cores = NCORES, mc.preschedule = FALSE)
out$wall_min <- as.numeric(Sys.time() - t0, units = "mins")
out$timing <- do.call(rbind, lapply(tim, function(x) data.frame(
  tier = x$job$tier, arm = if (is.null(x$job$g$arm)) "dich" else x$job$g$arm,
  n = x$job$g$n, k = x$job$g$k %||% 0L, ok = x$ok, sec = x$sec,
  actual_B = x$actual_B, kb = round(x$size / 1024))))
print(out$timing)
saveRDS(out, "fresh_rerun_smoke.rds")
message("smoke test done, wall ", round(out$wall_min, 1), " min")
