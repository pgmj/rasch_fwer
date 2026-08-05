# Calibration for infit_short_scale.qmd: which rho injects a realistic infit
# level into a five-item polytomous scale.
#
# calibration_rho_infit.R showed that rho is not portable across item formats.
# This shows it is not portable across item counts either. At 15 polytomous
# items rho = .85 injects infit around 1.32, while at 5 items the same rho
# injects only 1.20. A short scale needs stronger misfit for the same infit,
# because the misfitting item is a larger share of the total score it is
# conditioned on and therefore attenuates its own statistic more.
#
# Result: rho = 0.75 injects infit ~1.31 at five items, inside the 1.23 to 1.32
# band that the real phq9 items reach.
#
# Run from this directory:  Rscript calibration_rho_short.R
# Results are cached to calibration_rho_short.rds.

source("../R/helpers.R")
suppressMessages({
  library(easyRasch2)
  library(MASS)
  library(parallel)
})

RES    <- "calibration_rho_short.rds"
RHOS   <- c(0.75, 0.85, 0.90, 0.95)
NS     <- c(250, 600)
R_REPS <- 80
NCORES <- 10

# Same arm definition as infit_short_scale.qmd
LOCS <- seq(-2, 2, length.out = 5)
THR  <- lapply(LOCS, function(L) L + seq(-1, 1, length.out = 3))
MIS  <- which.min(abs(LOCS))   # item 3, at 0 logits
SD   <- 1.5

sim_short <- function(n, rho, seed) {
  set.seed(seed)
  S <- SD^2 * matrix(c(1, rho, rho, 1), 2, 2)
  th <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = S)
  d <- sim_partial_score(THR, th[, 1])
  d[, MIS] <- sim_partial_score(THR[MIS], th[, 2])[, 1]
  d <- as.data.frame(d)
  names(d) <- paste0("V", seq_along(LOCS))
  d
}

if (!file.exists(RES)) {
  t0 <- Sys.time()
  grid <- expand.grid(rep = 1:R_REPS, rho = RHOS, n = NS)
  vals <- mclapply(seq_len(nrow(grid)), function(i) {
    d <- sim_short(grid$n[i], grid$rho[i], 5e5 + i)
    x <- try(RMitemInfit(d, output = "dataframe"), silent = TRUE)
    if (inherits(x, "try-error")) NA_real_ else x$Infit_MSQ[MIS]
  }, mc.cores = NCORES)
  grid$infit <- unlist(vals)
  saveRDS(list(grid = grid, locs = LOCS, misfit_item = MIS, sd = SD,
               reps = R_REPS,
               elapsed_min = as.numeric(Sys.time() - t0, "mins")), RES)
}

cal <- readRDS(RES)
agg <- aggregate(infit ~ rho + n, data = cal$grid, FUN = mean, na.rm = TRUE)
agg$infit <- round(agg$infit, 3)
print(agg[order(agg$n, agg$rho), ], row.names = FALSE)
