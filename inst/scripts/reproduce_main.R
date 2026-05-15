###############################################################################
# reproduce_main.R
#
# Reproduces the main-text simulation tables and figures of the manuscript:
#   Table 3 : Correct OBD selection probability (%, MCSE), by scenario x design
#   Table S3: Mean utility loss, near-optimal selection rate, no-selection rate
#   Table S4: Patient allocation by exposure category, mean toxicity burden
#   Figure 1: Dose-level selection probability bars for S4 block
#   Figure 2: Dose-level selection probability bars for S7 block
#
# Designs compared (12 scenarios x 4 designs = 48 cells):
#   joint    : posterior predictive utility under joint Gaussian model
#   marginal : posterior predictive utility evaluating T, E independently
#   plugin   : utility at posterior mean (no integration)
#   boin12   : BOIN12 (Lin et al. 2020) with binary endpoints
#
# Runtime: about 30 to 60 minutes total on a recent laptop with
# parallel::mclapply. Set N_REPS = 20 and SCEN_KEYS = "S4b" for a
# fast smoke test.
###############################################################################

suppressPackageStartupMessages({
  library(ppod)
  library(parallel)
})

# ---- Configuration ----
N_REPS    <- 400L
N_CORES   <- max(1L, parallel::detectCores() - 1L)
OUT_DIR   <- "output_main"
SCEN_KEYS <- names(ppod_scenarios)
DESIGNS   <- c("joint", "marginal", "plugin", "boin12")

# Allow command-line pilot.
args <- commandArgs(trailingOnly = TRUE)
if ("--pilot" %in% args || identical(Sys.getenv("PILOT"), "1")) {
  N_REPS    <- 20L
  SCEN_KEYS <- c("S1", "S4b", "S4d")
  cat("=== PILOT MODE: 20 reps x 3 scenarios ===\n")
}

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
set.seed(20260515)

cat("ppod main-text reproducibility script\n")
cat("Scenarios:", paste(SCEN_KEYS, collapse = ", "), "\n")
cat("Designs:  ", paste(DESIGNS, collapse = ", "), "\n")
cat("Replicates per cell:", N_REPS, "\n")
cat("Parallel workers:   ", N_CORES, "\n\n")

# ---- Run all (scenario, design) cells ----
all_results <- vector("list", length(SCEN_KEYS))
names(all_results) <- SCEN_KEYS

for (scen_key in SCEN_KEYS) {
  cat("=== Scenario", scen_key, "===\n")
  scen <- ppod_scenarios[[scen_key]]
  cell <- vector("list", length(DESIGNS))
  names(cell) <- DESIGNS
  for (design in DESIGNS) {
    t0 <- Sys.time()
    cell[[design]] <- run_replicates(
      scen, design = design, n_reps = N_REPS,
      seed_base = 1000L * which(SCEN_KEYS == scen_key) +
                  100L * which(DESIGNS == design),
      parallel = .Platform$OS.type != "windows",
      n_cores  = N_CORES
    )
    dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    oc <- operating_characteristics(cell[[design]], scen)
    cat(sprintf("  %-9s  p_correct = %5.3f (SE %5.3f)  loss = %6.4f  %5.1fs\n",
                design, oc$p_correct, oc$se_correct, oc$mean_loss, dt))
  }
  all_results[[scen_key]] <- cell
}

# ---- Save raw results ----
saveRDS(all_results, file.path(OUT_DIR, "all_results.rds"))
cat("\nRaw results saved to", file.path(OUT_DIR, "all_results.rds"), "\n")

# ---- Summarize ----
oc_table <- do.call(rbind, lapply(SCEN_KEYS, function(sk) {
  scen <- ppod_scenarios[[sk]]
  do.call(rbind, lapply(DESIGNS, function(d) {
    oc <- operating_characteristics(all_results[[sk]][[d]], scen)
    data.frame(
      scenario   = sk,
      design     = d,
      p_correct  = oc$p_correct,
      se_correct = oc$se_correct,
      mean_loss  = oc$mean_loss,
      near_05    = unname(oc$p_near["eps_0.05"]),
      near_10    = unname(oc$p_near["eps_0.1"]),
      p_no_sel   = oc$p_no_select,
      alloc_over = unname(oc$alloc["over"]),
      alloc_opt  = unname(oc$alloc["optimal"]),
      alloc_sub  = unname(oc$alloc["suboptimal"]),
      tox_burden = oc$mean_tox_burden,
      stringsAsFactors = FALSE
    )
  }))
}))

write.csv(oc_table, file.path(OUT_DIR, "operating_characteristics.csv"),
           row.names = FALSE)
cat("Operating characteristics summary saved to",
    file.path(OUT_DIR, "operating_characteristics.csv"), "\n")

# ---- Table 3: correct OBD selection by design x scenario ----
tab3 <- reshape(
  oc_table[, c("scenario", "design", "p_correct", "se_correct")],
  idvar = "scenario", timevar = "design", direction = "wide"
)
write.csv(tab3, file.path(OUT_DIR, "table3_correct_selection.csv"),
           row.names = FALSE)
cat("Table 3 saved to",
    file.path(OUT_DIR, "table3_correct_selection.csv"), "\n")

# ---- Figures 1 and 2: dose-level selection bars ----
dose_selection_table <- function(cell, J) {
  sels <- vapply(cell, function(r) {
    if (is.null(r$selected) || is.na(r$selected)) NA_integer_ else r$selected
  }, integer(1L))
  tab <- tabulate(sels, nbins = J)
  c(tab, no_sel = sum(is.na(sels)))
}

make_panel_fig <- function(scen_keys, filename, mfrow = c(2, 2)) {
  pdf(filename, width = 8, height = 6)
  par(mfrow = mfrow, mar = c(4, 4, 2, 1))
  for (sk in scen_keys) {
    scen <- ppod_scenarios[[sk]]
    J <- length(scen$mu_T)
    mat <- sapply(DESIGNS, function(d) {
      dose_selection_table(all_results[[sk]][[d]], J) / N_REPS
    })
    rownames(mat) <- c(paste0("d", seq_len(J)), "no_sel")
    barplot(t(mat), beside = TRUE, ylim = c(0, 1),
            main = sk, ylab = "Selection probability",
            legend = TRUE, args.legend = list(cex = 0.6, x = "topright"))
  }
  dev.off()
}

S4_keys <- intersect(SCEN_KEYS, c("S4a", "S4b", "S4c", "S4d"))
if (length(S4_keys) > 0L) {
  make_panel_fig(S4_keys, file.path(OUT_DIR, "fig_S4.pdf"),
                  mfrow = c(2, 2))
  cat("Figure 1 saved to", file.path(OUT_DIR, "fig_S4.pdf"), "\n")
}

S7_keys <- intersect(SCEN_KEYS, c("S7a", "S7b", "S7c"))
if (length(S7_keys) > 0L) {
  make_panel_fig(S7_keys, file.path(OUT_DIR, "fig_S7.pdf"),
                  mfrow = c(1, 3))
  cat("Figure 2 saved to", file.path(OUT_DIR, "fig_S7.pdf"), "\n")
}

cat("\nDone.\n")
