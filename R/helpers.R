# Shared helper functions for the simulation studies in simulations/ and for
# the analysis chunks in index.qmd.
#
# The simulation studies were originally developed in the easyRasch2 source
# tree (easyRasch2/dev/) and used pkgload::load_all() to reach package
# internals. This file instead accesses the installed easyRasch2 package
# (>= 1.1.0), which exports the phq9 dataset and contains the internal
# functions used by the data generators. That makes this folder self-contained
# apart from installed packages.

# ---- package internals used by the data generators --------------------------
sim_partial_score <- easyRasch2:::sim_partial_score
.fit_cml_thresholds <- easyRasch2:::.fit_cml_thresholds
.estimate_prior_sd <- easyRasch2:::.estimate_prior_sd
.grid_loglik <- easyRasch2:::.grid_loglik
.logp_tables <- easyRasch2:::.logp_tables

# ---- phq9 arm: CML item parameters from the packaged real responses ---------
# Returns thresholds, item mean locations, and the estimated latent SD.
# Identical to the setup block used in the original dev/ studies.
load_phq9_arm <- function() {
  e <- new.env()
  data("phq9", package = "easyRasch2", envir = e)
  obs <- as.matrix(e$phq9[, paste0("q", 1:9)])
  obs <- obs[complete.cases(obs), ]
  thr <- .fit_cml_thresholds(obs)
  locs <- vapply(thr, mean, numeric(1))
  g <- seq(-6, 6, length.out = 81)
  sd_est <- .estimate_prior_sd(.grid_loglik(obs, .logp_tables(thr, g), g), g, 0)
  list(obs = obs, thr = thr, locs = locs, sd = sd_est)
}

# ---- resumable simulation runner --------------------------------------------
# Cells already present in the cached grid (keyed on `keys`) are skipped, so
# raising a replication constant and re-rendering computes only the new cells.
# Results are checkpointed to disk after every chunk.
run_resumable <- function(res_file, grid, keys, fn,
                          ncores = NCORES, chunk_size = CHUNK,
                          extra = list()) {
  key_of <- function(g) {
    if (nrow(g) == 0L) return(character(0))
    do.call(paste, c(unname(as.list(g[keys])), sep = "|"))
  }
  acc <- if (file.exists(res_file)) readRDS(res_file) else
    c(list(res = list(), grid = grid[0, , drop = FALSE], elapsed_min = 0), extra)
  todo <- which(!(key_of(grid) %in% key_of(acc$grid)))
  if (length(todo) == 0L) {
    message("cached: ", nrow(acc$grid), " cells, nothing to run"); return(acc)
  }
  message("running ", length(todo), " new cells (", nrow(acc$grid), " cached)")
  for (s in seq(1, length(todo), by = chunk_size)) {
    idx <- todo[s:min(s + chunk_size - 1L, length(todo))]
    t0 <- Sys.time()
    out <- parallel::mclapply(idx, function(i) fn(grid[i, , drop = FALSE]),
                              mc.cores = ncores)
    # The per-cell functions swallow errors and return NULL so that an
    # occasional degenerate cell does not kill a long run. A whole chunk of
    # NULLs is not that: it means the cells cannot run at all (a missing
    # package attachment, say). Stop rather than silently caching failures,
    # which the resumable key would then treat as done.
    if (all(vapply(out, is.null, logical(1)))) {
      stop("run_resumable(): every cell in this chunk failed. ",
           "The results so far are unchanged on disk. Check that the study's ",
           "setup chunk attaches everything the cell function needs.",
           call. = FALSE)
    }
    acc$res <- c(acc$res, out)
    acc$grid <- rbind(acc$grid, grid[idx, , drop = FALSE])
    acc$elapsed_min <- acc$elapsed_min + as.numeric(Sys.time() - t0, "mins")
    acc <- utils::modifyList(acc, extra)
    saveRDS(acc, res_file)
    message("  checkpoint ", nrow(acc$grid), " / ", nrow(grid),
            "  (", round(acc$elapsed_min, 1), " min)")
  }
  acc
}

# ---- HDCI bounds and flag coding --------------------------------------------
# Bounds per item at a given width, from a cutoff object's stored draws.
bounds_at <- function(results, width, items) {
  vapply(items, function(it) {
    v <- results$InfitMSQ[results$Item == it]
    as.numeric(ggdist::hdci(v, .width = width))[1:2]
  }, numeric(2))
}

# 0 = not flagged, 1 = overfit (below interval), 2 = underfit (above)
flag_code <- function(msq, lo, hi) {
  ifelse(msq > hi, 2L, ifelse(msq < lo, 1L, 0L))
}

# ---- shared plot theme -------------------------------------------------------
# Used by every figure in the manuscript. Requires the DM Sans font to be
# installed on the system.
theme_paper <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size, base_family = "DM Sans") +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 12)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 12))
    )
}

# ---- table display -----------------------------------------------------------
# Blank a key column wherever it repeats within the group defined by the key
# columns to its left, so each label is printed once per group instead of on
# every row. Display only, and it assumes the data are already sorted by `cols`.
# Numeric key columns should be formatted to strings first, since the blanking
# converts the column to character and knitr::kable(digits=) will no longer
# apply to it.
collapse_keys <- function(df, cols) {
  df <- as.data.frame(df)
  # Keys must come from the unmodified values. Building them from columns that
  # earlier iterations have already blanked collapses distinct groups together
  # and silently drops labels.
  orig <- df[cols]
  for (j in seq_along(cols)) {
    key <- do.call(paste, c(unname(as.list(orig[cols[seq_len(j)]])), sep = "\r"))
    df[[cols[j]]] <- ifelse(duplicated(key), "", as.character(orig[[cols[j]]]))
  }
  df
}

# ---- result-binding helpers used by the paper's analysis chunks -------------
# Bind one stored part ("hdci" or "pval") of a cached study into a long
# data.frame, prefixing the grid columns of each cell.
bind_part <- function(sim, part, grid_cols) {
  do.call(rbind, lapply(seq_len(nrow(sim$grid)), function(i) {
    r <- sim$res[[i]][[part]]
    if (is.null(r)) return(NULL)
    cbind(sim$grid[i, grid_cols, drop = FALSE], r)
  }))
}
