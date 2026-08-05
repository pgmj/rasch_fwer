# Calibration: how the inter-dimensional correlation (rho) maps to conditional
# infit MSQ, separately for dichotomous and polytomous items.
#
# Used in the Method section to justify calibrating misfit on the infit scale
# rather than on rho. At identical item locations, theta distribution and rho,
# polytomous items (4 categories, thresholds spread +/- 1 logit) show roughly
# 2.8 times the infit excess of dichotomous items, because they carry more
# information per item. The ratio is specific to this category count and
# threshold spacing.
#
# Run from this directory:  Rscript calibration_rho_infit.R
# Results are cached to calibration_rho_infit.rds.

source("../R/helpers.R")
suppressMessages({
  library(easyRasch2)
  library(eRm)
  library(MASS)
  library(parallel)
})

RES    <- "calibration_rho_infit.rds"
RHOS   <- c(0.15, 0.50, 0.75, 0.85, 0.90)
N      <- 500
R_REPS <- 60
NCORES <- 10
SD     <- 1.5

# Item locations identical to the dichotomous (w, B) grid study
set.seed(20250727)
LOCS <- runif(20, -2, 2)
MIS <- integer(0)
for (tg in c(0, -1, -2)) MIS <- c(MIS, setdiff(order(abs(LOCS - tg)), MIS)[1])
THR_POLY <- lapply(LOCS, function(L) L + seq(-1, 1, length.out = 3))

sim_dich <- function(n, rho, seed) {
  set.seed(seed)
  S <- matrix(c(SD^2, rho * SD^2, rho * SD^2, SD^2), 2, 2)
  th <- MASS::mvrnorm(n, c(0, 0), S)
  W <- cbind(rep(1, 20), rep(0, 20))
  W[MIS[1], ] <- c(0, 1)
  d <- as.data.frame(eRm::sim.xdim(persons = th, items = LOCS,
                                   Sigma = S, weightmat = W))
  names(d) <- paste0("V", 1:20)
  d
}

sim_poly <- function(n, rho, seed) {
  set.seed(seed)
  S <- matrix(c(SD^2, rho * SD^2, rho * SD^2, SD^2), 2, 2)
  th <- MASS::mvrnorm(n, c(0, 0), S)
  d <- matrix(0L, n, 20)
  fit_items <- setdiff(1:20, MIS[1])
  d[, fit_items] <- sim_partial_score(THR_POLY[fit_items], th[, 1])
  d[, MIS[1]] <- sim_partial_score(THR_POLY[MIS[1]], th[, 2])[, 1]
  d <- as.data.frame(d)
  names(d) <- paste0("V", 1:20)
  d
}

# ---- part 2: anchoring the phq9 arm against real item fit -------------------
# (a) Conditional infit of the real phq9 responses, with B = 1000 cutoffs and
#     Westfall-Young corrected p-values. These observed values anchor the
#     "realistic" misfit level.
# (b) Mean infit of a single misfitting item (q5) simulated from the phq9 arm
#     parameters across rho, which identifies the rho that reproduces the
#     observed misfit magnitude.

PH <- load_phq9_arm()
PH_N <- 600
ANCHOR_RHOS <- c(0.75, 0.85, 0.90)
ANCHOR_ITEM <- 5L   # on-target, location -0.15

sim_ph <- function(n, rho, seed) {
  set.seed(seed)
  S <- PH$sd^2 * matrix(c(1, rho, rho, 1), 2, 2)
  th <- MASS::mvrnorm(n, c(0, 0), S)
  d <- sim_partial_score(PH$thr, th[, 1])
  d[, ANCHOR_ITEM] <- sim_partial_score(PH$thr[ANCHOR_ITEM], th[, 2])[, 1]
  d <- as.data.frame(d)
  names(d) <- paste0("q", 1:9)
  d
}

if (!file.exists(RES)) {
  t0 <- Sys.time()

  # format comparison (part 1)
  grid <- expand.grid(rep = 1:R_REPS, rho = RHOS,
                      format = c("dichotomous", "polytomous"),
                      stringsAsFactors = FALSE)
  vals <- mclapply(seq_len(nrow(grid)), function(i) {
    g <- grid[i, ]
    d <- if (g$format == "dichotomous") {
      sim_dich(N, g$rho, g$rep * 13 + round(g$rho * 100))
    } else {
      sim_poly(N, g$rho, g$rep * 13 + round(g$rho * 100))
    }
    x <- try(RMitemInfit(d, output = "dataframe"), silent = TRUE)
    if (inherits(x, "try-error")) NA_real_ else x$Infit_MSQ[MIS[1]]
  }, mc.cores = NCORES)
  grid$infit <- unlist(vals)

  # real phq9 item fit (part 2a)
  d_real <- as.data.frame(PH$obs)
  co <- RMitemInfitCutoff(d_real, iterations = 1000, parallel = TRUE,
                          n_cores = NCORES, seed = 1, verbose = FALSE)
  real_fit <- suppressWarnings(
    RMitemInfit(d_real, cutoff = co, p_value = TRUE, correction = "fwer",
                output = "dataframe"))

  # phq9 arm anchoring scan (part 2b)
  agrid <- expand.grid(rep = 1:R_REPS, rho = ANCHOR_RHOS,
                       stringsAsFactors = FALSE)
  avals <- mclapply(seq_len(nrow(agrid)), function(i) {
    g <- agrid[i, ]
    d <- sim_ph(PH_N, g$rho, g$rep * 17 + round(g$rho * 100))
    x <- try(RMitemInfit(d, output = "dataframe"), silent = TRUE)
    if (inherits(x, "try-error")) NA_real_ else x$Infit_MSQ[ANCHOR_ITEM]
  }, mc.cores = NCORES)
  agrid$infit <- unlist(avals)

  saveRDS(list(grid = grid, locs = LOCS, misfit_item = MIS[1], n = N,
               real_fit = real_fit, anchor = agrid, anchor_item = ANCHOR_ITEM,
               anchor_n = PH_N, phq9_sd = PH$sd, phq9_locs = PH$locs,
               elapsed_min = as.numeric(Sys.time() - t0, "mins")), RES)
}

cal <- readRDS(RES)
agg <- aggregate(infit ~ format + rho, data = cal$grid, FUN = mean,
                 na.rm = TRUE)
print(reshape(agg, idvar = "rho", timevar = "format", direction = "wide"))
cat("\nreal phq9 item fit:\n")
print(cal$real_fit[, c("Item", "Infit_MSQ", "padj_infit", "Flagged")], digits = 3)
cat("\nphq9 arm anchoring (mean infit of q5):\n")
print(aggregate(infit ~ rho, data = cal$anchor, FUN = mean, na.rm = TRUE))
