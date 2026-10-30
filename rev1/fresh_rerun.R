# Fresh-data rerun for revision 1 (APM-26-08-243). Design and rationale are in
# rev1/rerun_plan.md.
#
# Replaces the master-dataset sampling of the two dichotomous studies with data
# generated independently in every replication, and extends B to 10 000 under
# the null. One bootstrap per replication at the largest B, with every smaller
# B read from its first B draws (nested), so comparisons across B are paired
# within the bootstrap as well as within the dataset.
#
# Run from rev1/ (needs ~10 cores, see the projection in rerun_plan.md):
#   nohup caffeinate -i Rscript -e 'source("fresh_rerun.R")' > fresh_rerun.log 2>&1 &
# Resumable: rerunning skips cells already in the .rds files.
#
# Set FRESH_SMOKE=1 to only define the functions (fresh_rerun_smoke.R does).

suppressPackageStartupMessages({
  library(easyRasch2)
  library(eRm)
  library(MASS)
  library(parallel)
})
source("../R/helpers.R")   # run_resumable(), bounds_at(), flag_code(), sim_partial_score

# The standard arm pins moments = "simulated" (studentize by the B simulated
# values only, the rule of easyRasch2 <= 1.4.0 and of the submission). From
# 1.4.0.9000 the package default is joint studentization, and the `moments`
# argument exists only from that version, so the CRAN release 1.4.0 cannot run
# this script. Install the development version from GitHub (pgmj/easyRasch2).
# The cached .rds files were computed under 1.4.0, where the default was the
# simulated-only rule, so they agree with the pinned arm.
stopifnot(packageVersion("easyRasch2") >= "1.4.0.9000",
          is.function(easyRasch2:::.bootstrap_pvalues),
          "moments" %in% names(formals(easyRasch2:::.bootstrap_pvalues)),
          is.function(easyRasch2:::.wy_stepdown))

WIDTHS  <- c(0.90, 0.95, 0.99, 0.995, 0.999)
ALPHA   <- 0.05
NCORES  <- 10

# Tier 1: dichotomous null, B up to 10 000 (E3, E4, and Study 1 convergence)
B_NULL  <- c(100, 200, 400, 1000, 2000, 5000, 10000)
NS      <- c(150, 250, 500, 1000)
R_NULL  <- 500
RES_T1  <- "fresh_dich_null.rds"

# Tier 2: dichotomous misfit (k = 1, 3), B up to 2000 (rest of Study 2)
B_MIS   <- c(100, 200, 400, 1000, 2000)
R_MIS   <- 250
RES_T2  <- "fresh_dich_misfit.rds"

# Tier 3 (optional): polytomous null, B up to 10 000
R_POLY  <- 250
RES_T3  <- "fresh_poly_null.rds"

# ---- item parameters ----------------------------------------------------------
# Identical to make_master() in simulations/infit_width_iterations.qmd: the same
# seed gives the same 20 locations and the same three misfitting items. Only the
# respondents change, drawn fresh per replication instead of from a master.
dich_items <- local({
  set.seed(20250727)
  locs <- runif(20, -2, 2)
  mis <- integer(0)
  for (tg in c(0, -1, -2)) mis <- c(mis, setdiff(order(abs(locs - tg)), mis)[1])
  list(locs = locs, mis = mis)
})

sim_dich <- function(n, k, seed) {
  W <- cbind(rep(1, 20), rep(0, 20))
  use <- dich_items$mis[seq_len(k)]
  if (k > 0) W[use, ] <- matrix(c(0, 1), nrow = k, ncol = 2, byrow = TRUE)
  Sigma <- matrix(c(1.5^2, .15 * 1.5^2, .15 * 1.5^2, 1.5^2), 2, 2)
  set.seed(seed)
  th <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = Sigma)
  d <- as.data.frame(eRm::sim.xdim(persons = th, items = dich_items$locs,
                                   Sigma = Sigma, weightmat = W))
  names(d) <- paste0("V", 1:20)
  d
}

# Polytomous null arms, as in simulations/infit_width_iterations_poly_weak.qmd
poly_arms <- local({
  s15 <- seq(-2, 2, length.out = 15)
  ph <- load_phq9_arm()
  list(sim15 = list(thr = lapply(s15, function(L) L + seq(-1, 1, length.out = 3)),
                    sd = 1.5, nm = paste0("V", 1:15), ns = c(250, 1000)),
       phq9  = list(thr = ph$thr, sd = ph$sd, nm = paste0("q", 1:9),
                    ns = c(150, 250, 1000)))
})

sim_poly_null <- function(arm, n, seed) {
  a <- poly_arms[[arm]]
  set.seed(seed)
  th <- rnorm(n, 0, a$sd)
  d <- as.data.frame(sim_partial_score(a$thr, th))
  names(d) <- a$nm
  d
}

# ---- p-values on the first B draws -------------------------------------------
# Standard: easyRasch2's own routine with moments = "simulated" (studentized
# by the moments of the B simulated values). Pooled: the same statistic
# studentized by moments over the observed and the B simulated values
# together, which makes the first step of the step-down a symmetric function
# of all B + 1 datasets (R1.2). Pooled equals the package default from 1.4.0.9000.
pvals_pooled <- function(obs, sim) {
  B <- nrow(sim)
  all <- rbind(obs, sim)
  m <- colMeans(all)
  s <- apply(all, 2L, stats::sd)
  s[!is.finite(s) | s < 1e-8] <- NA_real_
  t_obs <- abs((obs - m) / s)
  t_sim <- abs(sweep(sweep(sim, 2L, m, `-`), 2L, s, `/`))
  p <- vapply(seq_along(obs), function(i) {
    if (is.na(t_obs[i])) NA_real_ else (1 + sum(t_sim[, i] >= t_obs[i])) / (B + 1)
  }, numeric(1))
  list(p = p, padj = easyRasch2:::.wy_stepdown(t_obs, t_sim, B))
}

# Everything one replication reports, at each nested B
eval_nested <- function(d, co, Bs) {
  items <- names(d)
  obs_tab <- RMitemInfit(d, output = "dataframe")
  obs <- setNames(obs_tab$Infit_MSQ[match(items, obs_tab$Item)], items)
  res <- co$results
  sim_all <- tapply(res$InfitMSQ, list(res$iteration, res$Item),
                    function(x) x[1L])[, items, drop = FALSE]
  # An occasional failed refit leaves slightly fewer draws than requested. The
  # nominal B is kept as the label and the draws actually used are recorded.
  out <- lapply(Bs, function(B_nom) {
    B <- min(B_nom, nrow(sim_all))
    sim <- sim_all[seq_len(B), , drop = FALSE]
    fw <- easyRasch2:::.bootstrap_pvalues(obs, sim, correction = "fwer",
                                          moments = "simulated")
    bh <- easyRasch2:::.bootstrap_pvalues(obs, sim, correction = "fdr_bh",
                                          moments = "simulated")
    pl <- pvals_pooled(obs, sim)
    sub <- res[res$iteration <= B, ]
    bnd <- do.call(rbind, lapply(WIDTHS, function(w) {
      b <- bounds_at(sub, w, items)
      data.frame(B = B_nom, w = w, item = seq_along(items), lo = b[1, ], hi = b[2, ],
                 code = flag_code(obs, b[1, ], b[2, ]))
    }))
    list(
      pval = data.frame(B = B_nom, B_used = B, item = seq_along(items), msq = unname(obs),
                        p_marg = fw$p, padj_fwer = fw$padj, padj_bh = bh$padj,
                        p_marg_pool = pl$p, padj_fwer_pool = pl$padj),
      hdci = bnd)
  })
  list(pval = do.call(rbind, lapply(out, `[[`, "pval")),
       hdci = do.call(rbind, lapply(out, `[[`, "hdci")),
       actual_B = nrow(sim_all))
}

# ---- one replication per tier ----------------------------------------------------
# Seeds are deterministic in the design cell and disjoint from the seeds of the
# original studies (6e7 + ... and 9e7 + ...).
one_rep_dich <- function(g, Bs) {
  out <- try({
    t0 <- Sys.time()
    d <- sim_dich(g$n, g$k, seed = 7e8 + 1e6 * g$k + 1e3 * g$n + g$rep)
    co <- RMitemInfitCutoff(d, iterations = max(Bs), parallel = FALSE,
                            seed = 8e8 + 1e6 * g$k + 1e3 * g$n + g$rep,
                            verbose = FALSE)
    r <- eval_nested(d, co, Bs)
    r$sec <- as.numeric(Sys.time() - t0, units = "secs")
    r
  }, silent = TRUE)
  if (inherits(out, "try-error")) NULL else out
}

one_rep_poly <- function(g, Bs) {
  out <- try({
    t0 <- Sys.time()
    ai <- match(g$arm, names(poly_arms))
    d <- sim_poly_null(g$arm, g$n, seed = 7.5e8 + 1e7 * ai + 1e3 * g$n + g$rep)
    co <- RMitemInfitCutoff(d, iterations = max(Bs), parallel = FALSE,
                            seed = 8.5e8 + 1e7 * ai + 1e3 * g$n + g$rep,
                            verbose = FALSE)
    r <- eval_nested(d, co, Bs)
    r$sec <- as.numeric(Sys.time() - t0, units = "secs")
    r
  }, silent = TRUE)
  if (inherits(out, "try-error")) NULL else out
}

# ---- grids --------------------------------------------------------------------
# mclapply() preschedules, giving core j the cells j, j + 10, j + 20, ... of a
# chunk. Each chunk is therefore 10 replications times all conditions, ordered
# by condition first, so every core gets one cell of each condition and the
# slow n = 1000 cells are spread evenly. Chunks are replication blocks, so an
# interrupted run stays balanced over conditions.
block_order <- function(g, cond_order) {
  g[order(ceiling(g$rep / NCORES), cond_order(g), g$rep), , drop = FALSE]
}
grid_t1 <- function(R = R_NULL) {
  g <- expand.grid(rep = seq_len(R), k = 0L, n = NS)
  block_order(g, function(g) -g$n)[, c("rep", "k", "n")]
}
grid_t2 <- function(R = R_MIS) {
  g <- expand.grid(rep = seq_len(R), k = c(1L, 3L), n = NS)
  block_order(g, function(g) -g$n * 10 + g$k)[, c("rep", "k", "n")]
}
grid_t3 <- function(R = R_POLY) {
  g <- do.call(rbind, lapply(names(poly_arms), function(a)
    expand.grid(rep = seq_len(R), arm = a, n = poly_arms[[a]]$ns,
                stringsAsFactors = FALSE)))
  block_order(g, function(g) -g$n * 10 + match(g$arm, names(poly_arms)))[
    , c("arm", "rep", "n")]
}
CHUNK_T1 <- NCORES * length(NS)        # 40
CHUNK_T2 <- NCORES * 2 * length(NS)    # 80
CHUNK_T3 <- NCORES * 5                 # 50

if (!identical(Sys.getenv("FRESH_SMOKE"), "1")) {
  message("Tier 1: dichotomous null"); t0 <- Sys.time()
  run_resumable(RES_T1, grid_t1(), c("rep", "k", "n"),
                function(g) one_rep_dich(g, B_NULL), chunk_size = CHUNK_T1,
                extra = list(items = dich_items, B = B_NULL))
  message("Tier 2: dichotomous misfit")
  run_resumable(RES_T2, grid_t2(), c("rep", "k", "n"),
                function(g) one_rep_dich(g, B_MIS), chunk_size = CHUNK_T2,
                extra = list(items = dich_items, B = B_MIS))
  if (identical(Sys.getenv("FRESH_TIER3"), "1")) {
    message("Tier 3: polytomous null")
    run_resumable(RES_T3, grid_t3(), c("arm", "rep", "n"),
                  function(g) one_rep_poly(g, B_NULL), chunk_size = CHUNK_T3,
                  extra = list(B = B_NULL))
  }
  message("done in ", round(as.numeric(Sys.time() - t0, units = "hours"), 2), " h")
}
