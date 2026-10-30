#!/usr/bin/env Rscript
# p3_01_asr.R
#
# Ancestral state reconstruction for TPA PBP-related variant sites.
#
# Scope and interpretation:
#   Analyze all matrix loci with >=3 carriers and >=1 noncarrier.
#   Retain unknown calls as unknown. Count labels use the less frequent state.
#   These count categories do not certify inferential reliability.
#   A complete inventory and 1-2-carrier report are produced by p3_03.
#
# Framework:
#   1. States: 0 = wildtype (wt), 1 = mutant (mut)
#   2. Model: Mk ER, equal rates, conservative and robust
#   3. ASR method: marginal reconstruction via castor::asr_mk_model,
#      using the transition matrix fitted by castor::fit_mk, not re-optimized.
#   4. Independent-origin estimation:
#      Method A: posterior threshold crossing.
#                These are transition-boundary indicators, NOT strict
#                independent origins.
#      Method B: SIMMAP stochastic character mapping.
#                This is the reportable estimate of gain/loss counts, but
#                low-count loci are marked as exploratory.
#   5. These variants are NOT labelled resistant/susceptible.
#
# Usage:
#   Rscript SC/p3_01_asr.R

suppressMessages({
    library(ape)
    library(castor)
    library(phytools)
    library(readxl)
})

# -----------------------------
# Config
# -----------------------------

TREE_FILE   <- "results/03_iqtree/tpa_ss14_rooted.treefile"
MATRIX_FILE <- Sys.getenv("PBP_MATRIX_FILE", unset = "results/05_extracted_variants/sample_locus_binary_matrix.tsv")
META_FILE   <- "TPA_master_metadata.xlsx"
OUTDIR      <- "results/06_asr"

source("script/pbp_mk_missing_helpers.R")
source("script/pbp_variant_scope_helpers.R")

# ASR threshold settings
THRESHOLD_HIGH <- 0.8   # P(mut) > 0.8 = high-confidence mutant node
THRESHOLD_LOW  <- 0.2   # P(mut) < 0.2 = high-confidence wildtype node

# Minimum mutant samples for ASR.
# Loci with fewer mutants are skipped entirely.
MIN_MUTANT <- 3

# Map every eligible ASR locus. Count categories are defined in
# pbp_variant_scope_helpers.R and use min(carriers, noncarriers).
SIMMAP_MIN_MUTANT <- MIN_MUTANT
SIMMAP_N <- 100
USE_SIMMAP <- TRUE

# Optional seed for reproducibility of stochastic maps.
SEED <- 42
set.seed(SEED)

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Helper functions
# -----------------------------

stop_if_missing <- function(path, label) {
    if (!file.exists(path)) {
        stop(sprintf("Missing %s: %s", label, path), call. = FALSE)
    }
}

safe_numeric_vector <- function(x, label) {
    y <- suppressWarnings(as.numeric(as.character(x)))

    bad <- (!is.na(x) & is.na(y)) | (!is.na(y) & !(y %in% c(0, 1)))
    if (any(bad)) {
        stop(
            sprintf(
                "Non-integer values detected in %s. Values must be 0, 1, or NA.",
                label
            ),
            call. = FALSE
        )
    }

    y
}

get_clade_info <- function(tree, node_idx, meta) {
    n_tips <- length(tree$tip.label)

    if (node_idx <= n_tips) {
        tips <- tree$tip.label[node_idx]
    } else {
        queue <- c(node_idx)
        all_desc <- c()

        while (length(queue) > 0) {
            cur <- queue[1]
            queue <- queue[-1]

            children <- tree$edge[tree$edge[, 1] == cur, 2]
            all_desc <- c(all_desc, children)
            queue <- c(queue, children[children > n_tips])
        }

        tips <- tree$tip.label[all_desc[all_desc <= n_tips]]
    }

    tips <- tips[tips %in% rownames(meta)]

    if (length(tips) == 0) {
        return(list(
            n = 0,
            countries = "NA",
            clades = "NA",
            year_min = NA,
            year_max = NA
        ))
    }

    country_col <- "Country"
    clade_col <- "Clade (Lineage)"
    year_col <- "Isolation Year"

    countries <- "NA"
    clades <- "NA"

    if (country_col %in% colnames(meta)) {
        country_tab <- sort(table(meta[tips, country_col]), decreasing = TRUE)
        if (length(country_tab) > 0) {
            countries <- paste(
                names(country_tab)[seq_len(min(3, length(country_tab)))],
                collapse = ","
            )
        }
    }

    if (clade_col %in% colnames(meta)) {
        clade_tab <- sort(table(meta[tips, clade_col]), decreasing = TRUE)
        if (length(clade_tab) > 0) {
            clades <- paste(names(clade_tab), collapse = ",")
        }
    }

    years <- numeric(0)
    if (year_col %in% colnames(meta)) {
        years <- suppressWarnings(as.numeric(meta[tips, year_col]))
        years <- years[!is.na(years)]
    }

    list(
        n = length(tips),
        countries = countries,
        clades = clades,
        year_min = if (length(years) > 0) min(years) else NA,
        year_max = if (length(years) > 0) max(years) else NA
    )
}

get_aic <- function(f) {
    if (!is.null(f$AICc)) {
        return(f$AICc)
    }

    if (!is.null(f$AIC)) {
        return(f$AIC)
    }

    NA_real_
}

# These labels describe counts; they do not establish inferential reliability.
get_simmap_reliability <- function(n_mut, n_wt) {
    pbp_tip_count_category(n_mut, n_wt)
}
get_simmap_interpretation <- function(n_mut, n_wt) {
    pbp_count_interpretation(n_mut, n_wt)
}

count_simmap_transitions <- function(simmaps, nsim) {
    gain_counts <- numeric(nsim)
    loss_counts <- numeric(nsim)

    for (sim_i in seq_len(nsim)) {
        sm <- simmaps[[sim_i]]

        for (e_i in seq_len(nrow(sm$edge))) {
            hist_i <- sm$maps[[e_i]]
            states_i <- names(hist_i)

            if (length(states_i) < 2) {
                next
            }

            for (k in seq_len(length(states_i) - 1)) {
                if (states_i[k] == "0" && states_i[k + 1] == "1") {
                    gain_counts[sim_i] <- gain_counts[sim_i] + 1
                }

                if (states_i[k] == "1" && states_i[k + 1] == "0") {
                    loss_counts[sim_i] <- loss_counts[sim_i] + 1
                }
            }
        }
    }

    list(
        gain_counts = gain_counts,
        loss_counts = loss_counts
    )
}

run_simmap_for_locus <- function(tree_pruned, tip_states_01, fit_ER, locus) {
    tip_states_01_named <- tip_states_01
    names(tip_states_01_named) <- tree_pruned$tip.label
    tip_states_01_named <- tip_states_01_named[!is.na(tip_states_01_named)]

    if (length(unique(tip_states_01_named)) < 2) {
        cat("    [SIMMAP] skipped: only one observed state after removing NA\n")
        return(NULL)
    }

    tree_simmap <- drop.tip(
        tree_pruned,
        setdiff(tree_pruned$tip.label, names(tip_states_01_named))
    )

    if (length(tree_simmap$tip.label) < 3) {
        cat("    [SIMMAP] skipped: fewer than 3 tips after removing NA\n")
        return(NULL)
    }

    # Q matrix from fitted ER rate.
    q_rate <- fit_ER$transition_matrix[1, 2]

    if (!is.finite(q_rate) || q_rate <= 0) {
        cat("    [SIMMAP] skipped: invalid fitted ER rate\n")
        return(NULL)
    }

    Q <- matrix(
        c(
            -q_rate, q_rate,
            q_rate, -q_rate
        ),
        2,
        2,
        dimnames = list(c("0", "1"), c("0", "1"))
    )

    tip_factor <- factor(
        as.character(tip_states_01_named[tree_simmap$tip.label]),
        levels = c("0", "1")
    )
    names(tip_factor) <- tree_simmap$tip.label

    simmaps <- tryCatch(
        suppressMessages(
            make.simmap(
                tree_simmap,
                tip_factor,
                Q = Q,
                nsim = SIMMAP_N,
                message = FALSE
            )
        ),
        error = function(e) {
            cat(sprintf("    [SIMMAP] failed for %s: %s\n", locus, e$message))
            NULL
        }
    )

    if (is.null(simmaps)) {
        return(NULL)
    }

    counts <- count_simmap_transitions(simmaps, SIMMAP_N)

    gain_counts <- counts$gain_counts
    loss_counts <- counts$loss_counts

    list(
        gain_mean = round(mean(gain_counts), 2),
        gain_sd = round(sd(gain_counts), 2),
        gain_q025 = round(as.numeric(quantile(gain_counts, 0.025)), 2),
        gain_q500 = round(as.numeric(quantile(gain_counts, 0.500)), 2),
        gain_q975 = round(as.numeric(quantile(gain_counts, 0.975)), 2),
        loss_mean = round(mean(loss_counts), 2),
        loss_sd = round(sd(loss_counts), 2),
        loss_q025 = round(as.numeric(quantile(loss_counts, 0.025)), 2),
        loss_q500 = round(as.numeric(quantile(loss_counts, 0.500)), 2),
        loss_q975 = round(as.numeric(quantile(loss_counts, 0.975)), 2)
    )
}

# -----------------------------
# Read data
# -----------------------------

stop_if_missing(TREE_FILE, "TREE_FILE")
stop_if_missing(MATRIX_FILE, "MATRIX_FILE")
stop_if_missing(META_FILE, "META_FILE")

cat("Reading tree...\n")
tree_full <- read.tree(TREE_FILE)
cat(sprintf(
    "  raw tree: %d tips, %d internal nodes\n",
    length(tree_full$tip.label),
    tree_full$Nnode
))

cat("Reading variant genotype matrix...\n")
mat_full <- read.table(
    MATRIX_FILE,
    sep = "\t",
    header = TRUE,
    row.names = 1,
    check.names = FALSE
)
cat(sprintf(
    "  matrix: %d samples x %d loci\n",
    nrow(mat_full),
    ncol(mat_full)
))

cat("Reading metadata...\n")
meta <- as.data.frame(read_excel(META_FILE))

if (!"SRA_ID (Key)" %in% colnames(meta)) {
    stop("META_FILE must contain column: SRA_ID (Key)", call. = FALSE)
}

rownames(meta) <- as.character(meta[["SRA_ID (Key)"]])

# -----------------------------
# Align tree and matrix
# -----------------------------

common_samples <- intersect(tree_full$tip.label, rownames(mat_full))
cat(sprintf("Samples shared by tree and matrix: %d\n", length(common_samples)))

if (length(common_samples) == 0) {
    stop("No overlapping sample IDs between tree and matrix.", call. = FALSE)
}

tips_to_drop <- setdiff(tree_full$tip.label, common_samples)

if (length(tips_to_drop) > 0) {
    cat(sprintf("Dropping %d tips not in matrix\n", length(tips_to_drop)))
    tree_pruned <- drop.tip(tree_full, tips_to_drop)
} else {
    tree_pruned <- tree_full
}

cat(sprintf("Pruned tree: %d tips\n\n", length(tree_pruned$tip.label)))

matrix_only <- setdiff(rownames(mat_full), tree_pruned$tip.label)
if (length(matrix_only) > 0) {
    cat(sprintf(
        "[Info] Matrix contains %d samples not in tree; they are ignored for ASR.\n",
        length(matrix_only)
    ))
}

# -----------------------------
# Main loop
# -----------------------------

cat("Running ASR...\n")
cat(sprintf("Total loci in matrix: %d\n", ncol(mat_full)))
cat(sprintf("MIN_MUTANT for ASR: %d\n", MIN_MUTANT))
cat(sprintf("SIMMAP_MIN_MUTANT: %d\n", SIMMAP_MIN_MUTANT))
cat("Count groups use the less frequent binary state.\n")
cat(sprintf("SIMMAP_N: %d\n\n", SIMMAP_N))

all_events <- list()
rate_summary <- data.frame()
fit_failures <- list()
node_prob_list <- list()

# Checkpoints are valid only for identical inputs, code, settings and runtime.
signature_files <- c(TREE_FILE, MATRIX_FILE, META_FILE,
    "script/p3_01_asr_simmap_binary_loci.R", "script/pbp_mk_missing_helpers.R",
    "script/pbp_variant_scope_helpers.R")
signature <- c(paste(names(tools::md5sum(signature_files)), tools::md5sum(signature_files), sep="="),
    paste("SIMMAP_N", SIMMAP_N, "MIN_MUTANT", MIN_MUTANT, "SEED", SEED),
    R.version.string, vapply(c("ape","castor","phytools","readxl"),
        function(pkg) paste(pkg,packageVersion(pkg)), character(1)))
signature_tmp <- tempfile(); writeLines(signature, signature_tmp)
signature_id <- unname(tools::md5sum(signature_tmp)); unlink(signature_tmp)
checkpoint_dir <- file.path(OUTDIR, "recovery_checkpoints", signature_id)
dir.create(checkpoint_dir, recursive=TRUE, showWarnings=FALSE)
writeLines(signature, file.path(checkpoint_dir,"input_signature.txt"))
atomic_tsv <- function(x, path) {
    temporary <- paste0(path,".tmp")
    write.table(x,temporary,sep="\t",row.names=FALSE,quote=FALSE,na="NA")
    if (!file.rename(temporary,path)) stop("Cannot atomically save ",path)
}

for (locus in colnames(mat_full)) {
    cat(sprintf("  [%s]\n", locus))

    tip_states_01 <- safe_numeric_vector(
        mat_full[tree_pruned$tip.label, locus],
        locus
    )
    names(tip_states_01) <- tree_pruned$tip.label

    if (!all(tip_states_01 %in% c(0L, 1L, NA))) {
        stop(
            sprintf(
                "Locus %s contains values other than 0, 1, or NA.",
                locus
            ),
            call. = FALSE
        )
    }

    # castor uses states 1 and 2, not 0 and 1.
    # Original biological encoding:
    #   0 = wt  -> castor state 1
    #   1 = mut -> castor state 2
    tip_states_12 <- tip_states_01 + 1L
    names(tip_states_12) <- tree_pruned$tip.label

    stopifnot(all(names(tip_states_12) == tree_pruned$tip.label))

    n_mut <- sum(tip_states_12 == 2, na.rm = TRUE)
    n_wt <- sum(tip_states_12 == 1, na.rm = TRUE)
    n_na <- sum(is.na(tip_states_12))

    if ((n_mut + n_wt) == 0) {
        cat("    [SKIP] no non-NA states\n")
        next
    }

    freq <- n_mut / (n_mut + n_wt)

    cat(sprintf(
        "    mutant: %d (%.1f%%) | wildtype: %d | NA: %d\n",
        n_mut,
        freq * 100,
        n_wt,
        n_na
    ))

    if (n_mut < MIN_MUTANT) {
        cat(sprintf("    [SKIP] mutant count < %d\n", MIN_MUTANT))
        next
    }

    if (n_wt == 0) {
        cat("    [SKIP] no wildtype samples; binary Mk model is not informative\n")
        next
    }

    checkpoint_path <- file.path(checkpoint_dir,
        paste0(gsub("[^A-Za-z0-9]", "_", locus),".rds"))
    if (file.exists(checkpoint_path) && !length(readRDS(checkpoint_path)$failures)) {
        saved <- readRDS(checkpoint_path)
        if (!identical(saved$locus,locus) || !identical(saved$signature_id,signature_id))
            stop("Checkpoint identity mismatch: ",locus)
        rate_summary <- rbind(rate_summary,saved$row)
        all_events <- c(all_events,saved$events)
        node_prob_list[[locus]] <- saved$probs
        if (length(saved$failures)) fit_failures[[locus]] <- saved$failures
        assign(".Random.seed",saved$rng_after,envir=.GlobalEnv)
        cat("    [CHECKPOINT] reused saved locus result\n")
        atomic_tsv(rate_summary,file.path(OUTDIR,"asr_evolutionary_rates.partial.tsv"))
        next
    }

    simmap_reliability <- get_simmap_reliability(n_mut, n_wt)
    low_count_flag <- min(n_mut, n_wt) < 10
    simmap_interpretation <- get_simmap_interpretation(n_mut, n_wt)

    if (low_count_flag) {
        cat(sprintf(
            "    [LOW-COUNT] n_mut=%d; SIMMAP will be reported as exploratory/background only\n",
            n_mut
        ))
    }

    # -----------------------------
    # Step 1: fit ER model
    # -----------------------------

    fit_ER <- tryCatch(
        pbp_fit_mk(tree_pruned, tip_states_01, rate_model = "ER"),
        error = function(e) {
            cat(sprintf("    [ERROR] fit_mk failed: %s\n", e$message))
            NULL
        }
    )

    if (is.null(fit_ER) || !isTRUE(fit_ER$success)) {
        msg <- if (!is.null(fit_ER$error)) fit_ER$error else "fit failed"
        fit_failures[[locus]] <- data.frame(locus = locus, stage = "fit_ER", error = msg)
        cat(sprintf("    [ERROR] fit_mk did not converge: %s\n", msg))
        next
    }

    # -----------------------------
    # Step 2: ASR using fitted transition matrix
    # -----------------------------

    asr_result <- tryCatch(
        pbp_asr_mk(
            tree = tree_pruned, states01 = tip_states_01,
            transition_matrix = fit_ER$transition_matrix, reroot = TRUE
        ),
        error = function(e) {
            cat(sprintf("    [ERROR] ASR failed: %s\n", e$message))
            NULL
        }
    )

    if (is.null(asr_result)) {
        fit_failures[[locus]] <- data.frame(locus = locus, stage = "asr_ER", error = "See ASR error above")
        next
    }

    internal_probs <- asr_result$ancestral_likelihoods
    n_tips <- length(tree_pruned$tip.label)
    n_internal <- tree_pruned$Nnode

    # -----------------------------
    # Method A: posterior threshold crossing
    # -----------------------------

    threshold_events <- list()

    n_gain_internal <- 0
    n_loss_internal <- 0
    n_gain_terminal <- 0
    n_loss_terminal <- 0

    for (edge_i in seq_len(nrow(tree_pruned$edge))) {
        parent_node <- tree_pruned$edge[edge_i, 1]
        child_node <- tree_pruned$edge[edge_i, 2]

        branch_len <- NA_real_
        if (!is.null(tree_pruned$edge.length)) {
            branch_len <- tree_pruned$edge.length[edge_i]
        }

        is_terminal <- child_node <= n_tips

        parent_internal <- parent_node - n_tips
        if (parent_internal < 1 || parent_internal > n_internal) {
            next
        }

        p_parent_mut <- internal_probs[parent_internal, 2]

        if (is_terminal) {
            s <- tip_states_12[tree_pruned$tip.label[child_node]]
            if (is.na(s)) {
                next
            }

            p_child_mut <- if (s == 2) 1.0 else 0.0
        } else {
            child_internal <- child_node - n_tips
            if (child_internal < 1 || child_internal > n_internal) {
                next
            }

            p_child_mut <- internal_probs[child_internal, 2]
        }

        if (is.na(p_parent_mut) || is.na(p_child_mut)) {
            next
        }

        event_type <- NULL

        if (p_parent_mut < THRESHOLD_LOW && p_child_mut > THRESHOLD_HIGH) {
            event_type <- "gain"

            if (is_terminal) {
                n_gain_terminal <- n_gain_terminal + 1
            } else {
                n_gain_internal <- n_gain_internal + 1
            }
        } else if (p_parent_mut > THRESHOLD_HIGH && p_child_mut < THRESHOLD_LOW) {
            event_type <- "loss"

            if (is_terminal) {
                n_loss_terminal <- n_loss_terminal + 1
            } else {
                n_loss_internal <- n_loss_internal + 1
            }
        }

        if (!is.null(event_type)) {
            clade <- get_clade_info(tree_pruned, child_node, meta)

            threshold_events[[length(threshold_events) + 1]] <- data.frame(
                locus = locus,
                event_type = event_type,
                node = child_node,
                parent_node = parent_node,
                p_parent_mut = round(p_parent_mut, 4),
                p_node_mut = round(p_child_mut, 4),
                is_terminal = is_terminal,
                n_descendants = clade$n,
                top_countries = clade$countries,
                clades = clade$clades,
                year_min = clade$year_min,
                year_max = clade$year_max,
                branch_len = round(branch_len, 8),
                stringsAsFactors = FALSE
            )
        }
    }

    cat(sprintf(
        "    [threshold] internal gain: %d | internal loss: %d | terminal gain: %d | terminal loss: %d\n",
        n_gain_internal,
        n_loss_internal,
        n_gain_terminal,
        n_loss_terminal
    ))

    all_events <- c(all_events, threshold_events)

    # -----------------------------
    # Method B: SIMMAP
    # -----------------------------

    simmap_run <- FALSE
    simmap_gain_mean <- NA_real_
    simmap_gain_sd <- NA_real_
    simmap_gain_q025 <- NA_real_
    simmap_gain_q500 <- NA_real_
    simmap_gain_q975 <- NA_real_
    simmap_loss_mean <- NA_real_
    simmap_loss_sd <- NA_real_
    simmap_loss_q025 <- NA_real_
    simmap_loss_q500 <- NA_real_
    simmap_loss_q975 <- NA_real_

    if (USE_SIMMAP && n_mut >= SIMMAP_MIN_MUTANT) {
        cat(sprintf(
            "    [SIMMAP] running %d samples; tip_count_category=%s\n",
            SIMMAP_N,
            simmap_reliability
        ))

        simmap_result <- run_simmap_for_locus(
            tree_pruned = tree_pruned,
            tip_states_01 = tip_states_01,
            fit_ER = fit_ER,
            locus = locus
        )

        if (!is.null(simmap_result)) {
            simmap_run <- TRUE

            simmap_gain_mean <- simmap_result$gain_mean
            simmap_gain_sd <- simmap_result$gain_sd
            simmap_gain_q025 <- simmap_result$gain_q025
            simmap_gain_q500 <- simmap_result$gain_q500
            simmap_gain_q975 <- simmap_result$gain_q975

            simmap_loss_mean <- simmap_result$loss_mean
            simmap_loss_sd <- simmap_result$loss_sd
            simmap_loss_q025 <- simmap_result$loss_q025
            simmap_loss_q500 <- simmap_result$loss_q500
            simmap_loss_q975 <- simmap_result$loss_q975

            cat(sprintf(
                "    [SIMMAP] gain: %.2f +/- %.2f [%.2f, %.2f] | loss: %.2f +/- %.2f [%.2f, %.2f]\n",
                simmap_gain_mean,
                simmap_gain_sd,
                simmap_gain_q025,
                simmap_gain_q975,
                simmap_loss_mean,
                simmap_loss_sd,
                simmap_loss_q025,
                simmap_loss_q975
            ))
        } else {
            cat("    [SIMMAP] failed or skipped internally\n")
            fit_failures[[locus]] <- data.frame(locus=locus,stage="simmap_ER",
                error="SIMMAP failed; see per-locus console output")
        }
    } else {
        cat(sprintf(
            "    [SIMMAP] skipped: USE_SIMMAP=%s, n_mut=%d, SIMMAP_MIN_MUTANT=%d\n",
            USE_SIMMAP,
            n_mut,
            SIMMAP_MIN_MUTANT
        ))
    }

    # -----------------------------
    # Save summary row
    # -----------------------------

    q_gain <- fit_ER$transition_matrix[1, 2]
    q_loss <- fit_ER$transition_matrix[2, 1]

    rate_summary <- rbind(
        rate_summary,
        data.frame(
            locus = locus,
            n_mutant = n_mut,
            n_wildtype = n_wt,
            n_missing = n_na,
            minor_state_count = min(n_mut, n_wt),
            tip_count_category = simmap_reliability,
            n_total_nonmissing = n_mut + n_wt,
            mut_freq = round(freq, 6),
            mut_freq_percent = round(freq * 100, 3),

            AIC_ER = round(get_aic(fit_ER), 2),
            q_wt_to_mut = round(q_gain, 8),
            q_mut_to_wt = round(q_loss, 8),

            methodA_gain_internal = n_gain_internal,
            methodA_loss_internal = n_loss_internal,
            methodA_gain_terminal = n_gain_terminal,
            methodA_loss_terminal = n_loss_terminal,

            simmap_run = simmap_run,
            simmap_n = ifelse(simmap_run, SIMMAP_N, NA_integer_),
            simmap_gain_mean = simmap_gain_mean,
            simmap_gain_sd = simmap_gain_sd,
            simmap_gain_q025 = simmap_gain_q025,
            simmap_gain_q500 = simmap_gain_q500,
            simmap_gain_q975 = simmap_gain_q975,
            simmap_loss_mean = simmap_loss_mean,
            simmap_loss_sd = simmap_loss_sd,
            simmap_loss_q025 = simmap_loss_q025,
            simmap_loss_q500 = simmap_loss_q500,
            simmap_loss_q975 = simmap_loss_q975,

            low_count_flag = low_count_flag,
            simmap_reliability = simmap_reliability,
            simmap_interpretation = simmap_interpretation,

            stringsAsFactors = FALSE
        )
    )

    node_prob_list[[locus]] <- internal_probs
    # Persist a full-precision result and RNG state before the next locus.
    saved <- list(locus=locus,signature_id=signature_id,
        row=rate_summary[nrow(rate_summary),,drop=FALSE],events=threshold_events,
        probs=internal_probs,failures=fit_failures[[locus]],rng_after=.Random.seed)
    saveRDS(saved,paste0(checkpoint_path,".tmp"))
    if (!file.rename(paste0(checkpoint_path,".tmp"),checkpoint_path))
        stop("Cannot save locus checkpoint: ",locus)
    atomic_tsv(rate_summary,file.path(OUTDIR,"asr_evolutionary_rates.partial.tsv"))
}

# -----------------------------
# Save results
# -----------------------------

cat("\nSaving results...\n")

failure_file <- file.path(OUTDIR, "asr_fit_failures.tsv")
fail_df <- if (length(fit_failures)) do.call(rbind, fit_failures) else
    data.frame(locus = character(), stage = character(), error = character())
write.table(fail_df, failure_file, sep = "\t", row.names = FALSE, quote = FALSE)
if (!nrow(rate_summary)) stop("All eligible ASR fits failed; no valid ASR result generated")
# Write successful outputs before applying the focal-locus completion gate.
if (nrow(fail_df)) cat(sprintf("[UNRESOLVED] %d loci; inspect %s\n",nrow(fail_df),failure_file))

rates_file <- file.path(OUTDIR, "asr_evolutionary_rates.tsv")
write.table(
    rate_summary,
    rates_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)
cat(sprintf(
    "[OUT] rates + event counts: evolutionary_rates.tsv (%d loci)\n",
    nrow(rate_summary)
))

if (length(all_events) > 0) {
    events_df <- do.call(rbind, all_events)

    events_file <- file.path(OUTDIR, "asr_threshold_events.tsv")
    write.table(
        events_df,
        events_file,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )

    cat(sprintf(
        "[OUT] threshold events: threshold_events.tsv (%d rows)\n",
        nrow(events_df)
    ))
    cat("  Note: these are posterior threshold crossings, not strict independent origins.\n")
    cat("  For independent origins, use SIMMAP gain results.\n")
} else {
    cat("[OUT] no threshold events detected; threshold_events.tsv not written\n")
}

# Per-locus node probabilities
for (locus in names(node_prob_list)) {
    probs <- node_prob_list[[locus]]

    df <- data.frame(
        node_internal = seq_len(nrow(probs)),
        node_global = seq(
            length(tree_pruned$tip.label) + 1,
            length(tree_pruned$tip.label) + nrow(probs)
        ),
        P_wildtype = round(probs[, 1], 6),
        P_mutant = round(probs[, 2], 6),
        stringsAsFactors = FALSE
    )

    locus_safe <- gsub("[^A-Za-z0-9]", "_", locus)

    write.table(
        df,
        file.path(OUTDIR, paste0(locus_safe, "_node_probs.tsv")),
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )
}

cat(sprintf("[OUT] node probability files: %d\n", length(node_prob_list)))

# Save reliability summary.
if (nrow(rate_summary) > 0) {
    reliability_summary <- as.data.frame(table(rate_summary$simmap_reliability))
    colnames(reliability_summary) <- c("simmap_reliability", "n_loci")

    write.table(
        reliability_summary,
        file.path(OUTDIR, "simmap_reliability_summary.tsv"),
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )

    cat("[OUT] SIMMAP reliability summary: simmap_reliability_summary.tsv\n")
}

# Successful outputs above remain available even if a focal locus failed.
focal_missing <- setdiff(PBP_FOCAL_LOCI,rate_summary$locus)
focal_mapping_failed <- rate_summary$locus[
    rate_summary$locus %in% PBP_FOCAL_LOCI & !rate_summary$simmap_run]
focal_failed <- intersect(PBP_FOCAL_LOCI,fail_df$locus)
focal_unresolved <- unique(c(focal_missing,focal_mapping_failed,focal_failed))
if (length(focal_unresolved)) stop(paste("Focal ASR/SIMMAP unresolved:",
    paste(focal_unresolved,collapse=", "),"; successful outputs saved"))
cat(sprintf("[ASR COMPLETION] %d result rows; %d explicitly unresolved loci; all six focal loci complete\n",
    nrow(rate_summary),nrow(fail_df)))

# -----------------------------
# Console summary
# -----------------------------

cat("\n======= ASR summary =======\n")
cat("Method notes:\n")
cat("  methodA_gain_internal: posterior threshold crossings, not independent origins.\n")
cat("  simmap_gain_mean: SIMMAP mean gain count; now run for all n_mutant >= MIN_MUTANT.\n")
cat(sprintf("  MIN_MUTANT for ASR/SIMMAP: %d\n", MIN_MUTANT))
cat("  Count categories use min(carriers, noncarriers); they do not establish reliability.\n")
cat("  Count categories describe sample counts; assess model and mapping sensitivity separately.\n\n")

if (nrow(rate_summary) > 0) {
    cat(sprintf(
        "%-38s %7s %6s %6s %12s %-24s\n",
        "locus",
        "mut_pct",
        "A_gain",
        "A_loss",
        "SIMMAP_gain",
        "SIMMAP_count_category"
    ))

    cat(paste(rep("-", 105), collapse = ""), "\n")

    for (i in seq_len(nrow(rate_summary))) {
        r <- rate_summary[i, ]

        simmap_str <- if (!is.na(r$simmap_gain_mean)) {
            sprintf("%.2f+/-%.2f", r$simmap_gain_mean, r$simmap_gain_sd)
        } else {
            "N/A"
        }

        cat(sprintf(
            "%-38s %7.1f%% %6d %6d %12s %-24s\n",
            r$locus,
            r$mut_freq * 100,
            r$methodA_gain_internal,
            r$methodA_loss_internal,
            simmap_str,
            r$simmap_reliability
        ))
    }

    cat("\nTip-count category counts:\n")
    print(table(rate_summary$simmap_reliability))

    cat("\nLow-count exploratory loci:\n")
    low_count_loci <- rate_summary$locus[rate_summary$low_count_flag]

    if (length(low_count_loci) > 0) {
        for (x in low_count_loci) {
            cat(sprintf("  - %s\n", x))
        }
    } else {
        cat("  none\n")
    }
}

cat("\nDone.\n")
cat(sprintf("Loci analyzed: %d\n", length(loci_done)))
cat(sprintf("Results written to: %s\n", OUTDIR))
