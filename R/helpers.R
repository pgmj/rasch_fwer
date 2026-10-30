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

# ---- checkpoint writing ------------------------------------------------------
# Write to a temporary file in the same directory, then rename it into place.
# The rename is atomic, so a write that is interrupted or blocked can never
# leave a truncated cache behind, and it is far less likely to collide with a
# cloud-sync client holding the target file open. These folders live under
# OneDrive, where a plain saveRDS() of a large cache intermittently fails with
# "cannot open the connection", so the write is also retried with a backoff.
save_checkpoint <- function(obj, path, tries = 5) {
  for (i in seq_len(tries)) {
    tmp <- paste0(path, ".tmp", Sys.getpid())
    ok <- tryCatch({
      saveRDS(obj, tmp)
      file.rename(tmp, path)
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (isTRUE(ok)) return(invisible(TRUE))
    unlink(tmp)
    Sys.sleep(2 * i)
  }
  stop("run_resumable(): could not write the checkpoint to ", path,
       " after ", tries, " attempts. The previous checkpoint is intact.",
       call. = FALSE)
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
    save_checkpoint(acc, res_file)
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

# ---- figure encodings: colour online, line type and shape in print -----------
# The journal prints in black and white (editor's comment E6, revision 1), so
# every figure distinguishes conditions by line type and point shape. Colour is
# kept as a redundant cue for the online version, from colour-blind-safe
# palettes checked with the dataviz skill's validator (OKLab, Machado et al.
# 2009 CVD simulation):
# - ordered variables (interval width, correlation) use one-hue blue ramps,
#   monotone in lightness, so they read as an order and print as a grey ramp;
# - the nominal decision rules use four Okabe-Ito hues, which pass the CVD and
#   normal-vision floors over all pairs (worst CVD dE 11.0, normal 18.4).
# Give colour, shape and linetype the same title so ggplot merges the legends.
bw_shapes <- c(16, 1, 17, 2, 15, 0)
bw_lines  <- c("solid", "42", "dotted", "dotdash", "longdash", "twodash")
pal_ramp5 <- c("#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#0d366b")  # light -> dark
pal_ramp3 <- c("#5598e7", "#256abf", "#0d366b")
pal_rules <- c("Westfall-Young" = "#0072B2", "Benjamini-Hochberg" = "#D55E00",
               "Interval, w = .999" = "#009E73", "Marginal p" = "#56B4E9")
scale_bw <- function(name = ggplot2::waiver(), colours = NULL, shape = TRUE,
                     linetype = TRUE) {
  out <- list(if (is.null(colours)) {
    ggplot2::scale_colour_grey(name = name, start = 0, end = 0.55)
  } else {
    ggplot2::scale_colour_manual(name = name, values = colours)
  })
  if (shape) out <- c(out, list(ggplot2::scale_shape_manual(name = name,
                                                            values = bw_shapes)))
  if (linetype) out <- c(out, list(ggplot2::scale_linetype_manual(name = name,
                                                                  values = bw_lines)))
  out
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

# ---- fresh-data dichotomous studies (revision 1) ------------------------------
# Studies 1 and 2 were rerun for the revision with respondents generated anew in
# every replication (rev1/fresh_rerun.R) instead of sampled from one master
# dataset per condition. The null tier serves both studies: one bootstrap of
# 10 000 iterations per replication, every smaller B read from its first draws.
# This returns the two caches in the structure of the submitted
# infit_width_iterations_S1.rds and _S2.rds, so the analysis chunks need no
# change. own_cells / own_min count each run once in the Reproducibility totals.
load_fresh_dich <- function(dir = "rev1") {
  t1 <- readRDS(file.path(dir, "fresh_dich_null.rds"))
  t2 <- readRDS(file.path(dir, "fresh_dich_misfit.rds"))
  stopifnot(!any(vapply(t1$res, is.null, logical(1))),
            !any(vapply(t2$res, is.null, logical(1))))
  pcols <- c("B", "item", "p_marg", "padj_fwer", "padj_bh",
             "p_marg_pool", "padj_fwer_pool")
  slim <- function(r) list(hdci = r$hdci[, c("B", "w", "item", "code")],
                           pval = r$pval[, pcols])
  conv <- list(
    grid = t1$grid[, c("rep", "n")],
    res = lapply(t1$res, function(r) {
      h <- r$hdci
      data.frame(B = h$B, w = h$w, item = h$item, width = h$hi - h$lo)
    }),
    elapsed_min = t1$elapsed_min, own_cells = nrow(t1$grid),
    own_min = t1$elapsed_min)
  grid <- list(
    grid = rbind(t1$grid[, c("rep", "k", "n")], t2$grid[, c("rep", "k", "n")]),
    res = c(lapply(t1$res, slim), lapply(t2$res, slim)),
    misfit = t1$items$mis,
    elapsed_min = t1$elapsed_min + t2$elapsed_min,
    own_cells = nrow(t2$grid), own_min = t2$elapsed_min)
  list(conv = conv, grid = grid)
}
