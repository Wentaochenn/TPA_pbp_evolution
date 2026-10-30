#!/usr/bin/env Rscript
# p3_02_homoplasy.R
# Purpose:
#   1. Consistency Index (CI) and Retention Index (RI) per locus
#   2. Parsimony minimum number of changes
#   3. Plot CI distribution and parsimony vs frequency
#
# CI = 1 / parsimony_score for a binary character (minimum changes m = 1).
# Lower CI = stronger homoplasy.
# Parsimony changes = lower bound on the number of state changes.
#
# Interpretation (must combine three lines of evidence):
#   low CI + high SIMMAP gain + scattered on tree -> multiple independent origins (e.g. A506T)
#   low CI + SIMMAP gain ~1 + concentrated in one clade -> single origin + clonal expansion (e.g. P564L)
#   A low CI alone does NOT prove multiple origins.
#
# Usage: Rscript p3_02_homoplasy.R

suppressMessages({
    library(ape)
    library(phangorn)
    library(phytools)
    library(readxl)
})

TREE_FILE   <- "results/03_iqtree/tpa_ss14_rooted.treefile"
MATRIX_FILE <- Sys.getenv("PBP_MATRIX_FILE", unset = "results/05_extracted_variants/sample_locus_binary_matrix.tsv")
META_FILE   <- "TPA_master_metadata.xlsx"
OUTDIR      <- "results/07_homoplasy"

# CI classification thresholds (consistent with plot below)
CI_EXTREME  <- 0.1   # CI < 0.1  : extreme homoplasy
CI_STRONG   <- 0.3   # CI < 0.3  : strong
CI_MODERATE <- 0.6   # CI < 0.6  : moderate; >= 0.6 weak

dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)

cat("Reading data...\n")
tree_full <- read.tree(TREE_FILE)
mat_full  <- read.table(MATRIX_FILE, sep="\t", header=TRUE,
                        row.names=1, check.names=FALSE)
meta      <- as.data.frame(read_excel(META_FILE))
rownames(meta) <- as.character(meta[["SRA_ID (Key)"]])

common <- intersect(tree_full$tip.label, rownames(mat_full))
tips_to_drop <- setdiff(tree_full$tip.label, common)
tree <- if (length(tips_to_drop) > 0)
    drop.tip(tree_full, tips_to_drop) else tree_full
cat(sprintf("Pruned tree: %d tips\n", length(tree$tip.label)))

# Analyze loci with >= 3 mutant samples
loci_all <- colnames(mat_full)
loci_ok  <- loci_all[sapply(loci_all, function(l) {
    v <- mat_full[tree$tip.label, l]
    sum(v==1, na.rm=TRUE) >= 3
})]
cat(sprintf("Loci analyzed: %d\n\n", length(loci_ok)))

# CI / RI / parsimony
results <- list()

for (locus in loci_ok) {
    states <- as.integer(mat_full[tree$tip.label, locus])
    states[is.na(states)] <- NA
    n_mut <- sum(states==1, na.rm=TRUE)
    n_wt  <- sum(states==0, na.rm=TRUE)

    states_char <- ifelse(is.na(states), NA, ifelse(states==1,"1","0"))
    names(states_char) <- tree$tip.label

    # Drop NA tips (parsimony does not allow missing data)
    valid_tips <- names(states_char)[!is.na(states_char)]
    if (length(valid_tips) < length(tree$tip.label)) {
        tree_pars <- drop.tip(tree, setdiff(tree$tip.label, valid_tips))
    } else {
        tree_pars <- tree
    }
    states_valid <- states_char[tree_pars$tip.label]

    pdat <- phyDat(as.matrix(states_valid), type="USER", levels=c("0","1"))

    # Parsimony score (minimum changes)
    pars_score <- parsimony(tree_pars, pdat)

    # CI = m / s; for binary character m = 1
    m_min <- 1L
    CI    <- if (pars_score > 0) m_min / pars_score else NA

    # RI = (g - s) / (g - m); g = max possible changes = min(n_mut, n_wt)
    g  <- min(n_mut, n_wt)
    RI <- tryCatch({
        if (g > m_min && g > pars_score)
            (g - pars_score) / (g - m_min)
        else if (g == pars_score)
            0
        else
            NA
    }, error=function(e) NA)

    cat(sprintf("  %-40s n_mut=%4d  pars=%3d  CI=%.3f  RI=%.3f\n",
                locus, n_mut, pars_score,
                ifelse(is.na(CI), 0, CI), ifelse(is.na(RI), 0, RI)))

    results[[locus]] <- data.frame(
        locus        = locus,
        n_mutant     = n_mut,
        n_wildtype   = n_wt,
        mut_freq     = round(n_mut/(n_mut+n_wt), 4),
        parsimony_changes = pars_score,
        CI           = round(CI, 4),
        RI           = round(RI, 4),
        homoplasy_level = ifelse(is.na(CI), "unknown",
                          ifelse(CI < CI_EXTREME,  "extreme",
                          ifelse(CI < CI_STRONG,   "strong",
                          ifelse(CI < CI_MODERATE, "moderate", "weak")))),
        stringsAsFactors = FALSE
    )
}

df_out <- do.call(rbind, results)
df_out <- df_out[order(df_out$CI, na.last=TRUE), ]

out_path <- file.path(OUTDIR, "homoplasy_ci_ri.tsv")
write.table(df_out, out_path, sep="\t", row.names=FALSE, quote=FALSE)
cat(sprintf("\n[OUT] CI/RI table: %s\n", out_path))

# Summary
cat("\n======= Homoplasy summary =======\n")
cat(sprintf("%-40s %7s %6s %6s %6s %12s\n",
            "locus", "mut%", "pars", "CI", "RI", "homoplasy"))
cat(paste(rep("-",85), collapse=""), "\n")
for (i in seq_len(nrow(df_out))) {
    r <- df_out[i,]
    cat(sprintf("%-40s %7.1f%% %6d %6.3f %6.3f %12s\n",
                r$locus, r$mut_freq*100, r$parsimony_changes,
                ifelse(is.na(r$CI), 0, r$CI),
                ifelse(is.na(r$RI), 0, r$RI), r$homoplasy_level))
}

# Plots
pdf(file.path(OUTDIR, "homoplasy_summary.pdf"), width=12, height=6)
par(mfrow=c(1,2), mar=c(8,4,3,2))

# Plot 1: CI barplot (colors use the same thresholds as classification)
ci_vals <- df_out$CI
ci_vals[is.na(ci_vals)] <- 0
col_ci  <- ifelse(ci_vals < CI_EXTREME,  "#E74C3C",
           ifelse(ci_vals < CI_STRONG,   "#E67E22",
           ifelse(ci_vals < CI_MODERATE, "#F1C40F", "#27AE60")))
barplot(ci_vals,
        names.arg = sub(".*:p\\.", "", df_out$locus),
        las=2, col=col_ci, cex.names=0.6,
        ylab="Consistency Index (CI)",
        main="Per-locus CI (lower = stronger homoplasy)",
        ylim=c(0,1))
abline(h=c(CI_EXTREME, CI_STRONG, CI_MODERATE), lty=2,
       col=c("#E74C3C","#E67E22","#27AE60"))
legend("topright", bty="n", cex=0.7,
       legend=c("CI<0.1: extreme homoplasy",
                "0.1-0.3: strong",
                "0.3-0.6: moderate",
                ">0.6: weak"),
       fill=c("#E74C3C","#E67E22","#F1C40F","#27AE60"))

# Plot 2: parsimony changes vs mutant frequency
plot(df_out$mut_freq*100, df_out$parsimony_changes,
     pch=16, col="#3498DB", cex=1.0,
     xlab="Mutant frequency (%)",
     ylab="Parsimony minimum changes",
     main="Mutant frequency vs minimum changes")
text(df_out$mut_freq*100, df_out$parsimony_changes,
     labels=sub(".*:p\\.", "", df_out$locus), cex=0.5, pos=3)
abline(a=0, b=max(df_out$parsimony_changes, na.rm=TRUE)/
           max(df_out$mut_freq*100, na.rm=TRUE),
       lty=2, col="gray")

dev.off()
cat(sprintf("[OUT] plot: %s\n", file.path(OUTDIR, "homoplasy_summary.pdf")))

cat("\nDone.\n")
cat("Interpretation reminder:\n")
cat("  Parsimony changes = lower bound on number of state changes.\n")
cat("  Combine CI with SIMMAP gain and tree distribution before\n")
cat("  claiming multiple independent origins (low CI alone is not enough).\n")
