# Calibration for infit_outfit.qmd: which rho puts the dichotomous arm at the
# same infit level as the polytomous phq9 arm, separately for an on-target and
# an off-target misfitting item.
#
# The main calibration script (calibration_rho_infit.R) anchors a single
# on-target item. This study also misfits an off-target item, and a misfitting
# item far from the person mean injects less infit at the same rho, so the
# anchoring has to be checked at both targeting levels.
#
# Result: rho = 0.60 for the dichotomous arm and rho = 0.85 for phq9 give a
# mean injected infit of 1.274 and 1.268 respectively, both inside the range
# the real phq9 items reach (1.23 to 1.32). Outfit is reported alongside, and
# shows that the same infit injection produces a much larger outfit excess in
# the dichotomous arm. That is a property of the arms, not a prediction about
# detection, since outfit's null variance is larger too.
#
# Run from this directory:  Rscript calibration_rho_targeting.R
# Results are cached to calibration_rho_targeting.rds.
#
# Requires easyRasch2 with the `statistic` argument in RMitemInfit()
# (added after 1.1.0).

source("../R/helpers.R")
suppressMessages({
  library(easyRasch2)
  library(eRm)
  library(MASS)
  library(parallel)
})

RES    <- "calibration_rho_targeting.rds"
N      <- 600
R_REPS <- 80
NCORES <- 10
D_RHOS <- c(0.50, 0.60, 0.65, 0.75)
P_RHOS <- 0.85

if (!"statistic" %in% names(formals(easyRasch2::RMitemInfit))) {
  stop("Installed easyRasch2 (", packageVersion("easyRasch2"), ") has no ",
       "`statistic` argument in RMitemInfit().", call. = FALSE)
}

# Same arm definitions as infit_outfit.qmd
set.seed(20250727)
D_LOCS <- runif(20, -2, 2)
ph <- load_phq9_arm()

ARMS <- list(
  dich20 = list(locs = D_LOCS, sd = 1.5, nm = paste0("V", 1:20),
                on = which.min(abs(D_LOCS)), off = which.max(D_LOCS),
                thr = NULL),
  phq9   = list(locs = ph$locs, sd = ph$sd, nm = paste0("q", 1:9),
                on = 5L, off = 8L, thr = ph$thr)
)

sim_arm <- function(arm, n, target, rho, seed) {
  a <- ARMS[[arm]]
  ni <- length(a$locs)
  mis <- switch(target, none = integer(0), on = a$on, off = a$off)
  set.seed(seed)
  S <- a$sd^2 * matrix(c(1, rho, rho, 1), 2, 2)
  th <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = S)
  if (is.null(a$thr)) {
    W <- cbind(rep(1, ni), rep(0, ni))
    if (length(mis)) W[mis, ] <- c(0, 1)
    d <- as.data.frame(eRm::sim.xdim(persons = th, items = a$locs,
                                     Sigma = S, weightmat = W))
  } else {
    d <- sim_partial_score(a$thr, th[, 1])
    if (length(mis)) d[, mis] <- sim_partial_score(a$thr[mis], th[, 2])[, 1]
    d <- as.data.frame(d)
  }
  names(d) <- a$nm
  list(data = d, mis = if (length(mis)) mis else 0L)
}

if (!file.exists(RES)) {
  t0 <- Sys.time()
  grid <- rbind(
    expand.grid(rep = 1:R_REPS, target = c("on", "off"), rho = D_RHOS,
                arm = "dich20", stringsAsFactors = FALSE),
    expand.grid(rep = 1:R_REPS, target = c("on", "off"), rho = P_RHOS,
                arm = "phq9", stringsAsFactors = FALSE)
  )
  vals <- mclapply(seq_len(nrow(grid)), function(i) {
    g <- grid[i, ]
    s <- sim_arm(g$arm, N, g$target, g$rho, 7e6 + i)
    x <- try(RMitemInfit(s$data, output = "dataframe"), silent = TRUE)
    y <- try(RMitemInfit(s$data, statistic = "outfit", output = "dataframe"),
             silent = TRUE)
    if (inherits(x, "try-error") || inherits(y, "try-error")) {
      return(c(NA_real_, NA_real_))
    }
    c(x$Infit_MSQ[s$mis], y$Outfit_MSQ[s$mis])
  }, mc.cores = NCORES)
  grid$infit <- vapply(vals, `[`, numeric(1), 1)
  grid$outfit <- vapply(vals, `[`, numeric(1), 2)
  saveRDS(list(grid = grid, n = N, reps = R_REPS, locs = D_LOCS,
               phq9_locs = ph$locs, phq9_sd = ph$sd,
               elapsed_min = as.numeric(Sys.time() - t0, "mins")), RES)
}

cal <- readRDS(RES)
agg <- aggregate(cbind(infit, outfit) ~ arm + rho + target, data = cal$grid,
                 FUN = mean, na.rm = TRUE)
agg[, c("infit", "outfit")] <- round(agg[, c("infit", "outfit")], 3)
print(agg[order(agg$arm, agg$rho, agg$target), ], row.names = FALSE)

cat("\nmean across targeting levels:\n")
m <- aggregate(infit ~ arm + rho, data = cal$grid, FUN = mean, na.rm = TRUE)
m$infit <- round(m$infit, 3)
print(m, row.names = FALSE)
