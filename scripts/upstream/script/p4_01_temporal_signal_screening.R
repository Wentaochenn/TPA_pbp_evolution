#!/usr/bin/env Rscript
# p4_01_temporal_signal_screening.R
#
# Integrated phylogeny-informed subsampling + date randomization test.
#
# This v2 script replaces / integrates:
#   - p4_01_subsample.R
#   - p4_02_date_randomization.R
#
# Main upgrades:
#   1. Extract and test:
#        - Full dataset
#        - SS14-only subtree
#        - Nichols-only subtree
#        - Multiple enlarged phylo-informed stratified SS14 subsamples
#
#   2. Improved stratified subsampling:
#        - Keeps all early anchors.
#        - Uses year-bin x continent stratification.
#        - Preserves lineage diversity within each bin when possible.
#        - Prefers larger root-to-tip distances, but avoids choosing only
#          near-duplicate representatives from the same local lineage.
#        - Adds focal PBP carrier rescue so key focal variants are represented.
#
#   3. Integrated strict date randomization:
#        - Root-to-tip regression is exploratory.
#        - Date randomization p-value is the main temporal-signal test.
#        - Positive slope is required.
#        - Automatically recommends the largest subsample retaining temporal signal.
#
# Inputs:
#   results/03_iqtree/tpa_ss14_rooted.treefile
#   sample_dates.tsv
#   TPA_master_metadata.xlsx
#   results/05_extracted_variants/sample_locus_binary_matrix.tsv
#
# Outputs:
#   results/11_temporal_signal/full_tree_used.nwk
#   results/11_temporal_signal/ss14_midpoint_rooted.nwk
#   results/11_temporal_signal/nichols_midpoint_rooted.nwk
#   results/11_temporal_signal/subsampled_phylo_informed_<profile>.tsv
#   results/11_temporal_signal/subsampled_tree_<profile>.nwk
#   results/11_temporal_signal/subsample_pbp_coverage.tsv
#   results/11_temporal_signal/temporal_signal_summary.tsv
#   results/11_temporal_signal/date_randomization_summary.pdf
#
# Compatibility outputs for downstream scripts:
#   results/11_temporal_signal/recommended_temporal_subset.tsv
#   results/11_temporal_signal/recommended_temporal_tree.nwk
#
# Usage:
#   Rscript SC/p4_01_temporal_signal_screening.R

suppressMessages({
    library(ape)
    library(phytools)
})

# -----------------------------
# Config
# -----------------------------

TREE_FILE   <- "results/03_iqtree/tpa_ss14_rooted.treefile"
DATES_FILE  <- "sample_dates.tsv"
META_FILE   <- "TPA_master_metadata.xlsx"
MATRIX_FILE <- Sys.getenv("PBP_MATRIX_FILE", unset = "results/05_extracted_variants/sample_locus_binary_matrix.tsv")

OUTDIR      <- "results/11_temporal_signal"

N_PERM      <- 1000
SEED        <- 42
set.seed(SEED)

# If TRUE, only samples present in the PBP binary matrix are eligible
# for stratified temporal subsampling.
REQUIRE_PBP <- TRUE

# Focal loci for coverage rescue.
# A506V is included because Figure 2 now treats it as a focal, model-sensitive
# site-506 variant.
PBP_FOCAL <- c(
    "mrcA:p.Met625Val",
    "mrcA:p.Gly708Ser",
    "mrcA:p.Ala506Thr",
    "mrcA:p.Ala506Val",
    "pbp1:p.Pro564Leu",
    "tp47:p.Ser394Arg"
)

# Temporal tiers
EARLY_CUTOFF     <- 1990
EXPANSION_END    <- 2010

# Candidate enlarged subsampling profiles.
# These names are approximate targets, not hard sample-size guarantees.
# The actual sample number depends on temporal/geographic/lineage distribution.
# After running this script, check subsample_profile_files.tsv for the
# realized n_samples per profile and rename profile labels below to match
# the actual achieved sizes if they drift from these target labels.
PROFILES <- data.frame(
    profile = c("pickup_200", "pickup_250", "pickup_400", "pickup_500", "pickup_600"),
    year_bin_width = c(2, 2, 1, 1, 1),
    max_per_bin_expansion = c(5, 6, 8, 10, 12),
    max_per_bin_recent = c(3, 4, 6, 7, 9),
    min_focal_carriers = c(12, 15, 25, 30, 35),
    stringsAsFactors = FALSE
)

# If TRUE, early samples are always retained in all candidate subsamples.
KEEP_ALL_EARLY <- TRUE

# Minimum samples required to test temporal signal.
MIN_N_TEMPORAL_TEST <- 10

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Logging helpers
# -----------------------------

info <- function(...) {
    cat(sprintf(...), "\n")
}

warn <- function(...) {
    warning(sprintf(...), call. = FALSE)
}

stop_if_missing <- function(path, label) {
    if (!file.exists(path)) {
        stop(sprintf("Missing %s: %s", label, path), call. = FALSE)
    }
}

# -----------------------------
# Basic helper functions
# -----------------------------

normalize_country <- function(x) {
    x <- as.character(x)
    x <- trimws(x)
    x <- gsub(" ", "_", x)
    x
}

get_continent <- function(country) {
    country <- normalize_country(country)

    continent_map <- c(
        USA = "NorthAmerica",
        United_States = "NorthAmerica",
        Canada = "NorthAmerica",
        Mexico = "NorthAmerica",

        UK = "Europe",
        United_Kingdom = "Europe",
        England = "Europe",
        France = "Europe",
        Germany = "Europe",
        Spain = "Europe",
        Netherlands = "Europe",
        Belgium = "Europe",
        Sweden = "Europe",
        Norway = "Europe",
        Denmark = "Europe",
        Switzerland = "Europe",
        Austria = "Europe",
        Italy = "Europe",
        Poland = "Europe",
        Czech_Republic = "Europe",
        Ireland = "Europe",
        Portugal = "Europe",
        Finland = "Europe",
        Hungary = "Europe",
        Romania = "Europe",
        Slovakia = "Europe",
        Greece = "Europe",

        Australia = "Oceania",
        New_Zealand = "Oceania",

        China = "Asia",
        Japan = "Asia",
        South_Korea = "Asia",
        Korea = "Asia",
        Taiwan = "Asia",
        Thailand = "Asia",
        Vietnam = "Asia",
        Singapore = "Asia",
        India = "Asia",
        Israel = "Asia",
        Turkey = "Asia",

        South_Africa = "Africa",
        Malawi = "Africa",
        Nigeria = "Africa",
        Kenya = "Africa",
        Uganda = "Africa",

        Brazil = "SouthAmerica",
        Argentina = "SouthAmerica",
        Colombia = "SouthAmerica",
        Peru = "SouthAmerica",
        Chile = "SouthAmerica"
    )

    out <- continent_map[country]
    out[is.na(out)] <- "Unknown"
    as.character(out)
}

assign_tier <- function(year) {
    ifelse(
        year < EARLY_CUTOFF,
        "early",
        ifelse(year <= EXPANSION_END, "expansion", "recent")
    )
}

calc_rtd <- function(tree) {
    if (is.null(tree$edge.length)) {
        stop("Tree has no branch lengths; cannot compute root-to-tip distances.", call. = FALSE)
    }

    depths <- node.depth.edgelength(tree)
    rtd_vec <- depths[seq_along(tree$tip.label)]
    names(rtd_vec) <- tree$tip.label

    rtd_vec
}

safe_read_dates <- function(path) {
    stop_if_missing(path, "DATES_FILE")

    dates_df <- read.table(
        path,
        sep = "\t",
        header = FALSE,
        col.names = c("sample", "year"),
        stringsAsFactors = FALSE
    )

    dates_df$sample <- as.character(dates_df$sample)
    dates_df$year <- suppressWarnings(as.numeric(dates_df$year))

    dates_df <- dates_df[!is.na(dates_df$sample) & !is.na(dates_df$year), ]

    if (nrow(dates_df) == 0) {
        stop("No valid dates found in sample_dates.tsv.", call. = FALSE)
    }

    setNames(dates_df$year, dates_df$sample)
}

safe_read_metadata <- function(path) {
    stop_if_missing(path, "META_FILE")

    suppressMessages({
        library(readxl)
    })

    meta_raw <- as.data.frame(read_excel(path), stringsAsFactors = FALSE)

    required_cols <- c("SRA_ID (Key)", "Clade (Lineage)", "Country", "Isolation Year")
    missing_cols <- setdiff(required_cols, colnames(meta_raw))

    if (length(missing_cols) > 0) {
        stop(
            sprintf(
                "Metadata missing required columns:\n%s",
                paste(missing_cols, collapse = "\n")
            ),
            call. = FALSE
        )
    }

    meta <- data.frame(
        sample = as.character(meta_raw[["SRA_ID (Key)"]]),
        lineage = as.character(meta_raw[["Clade (Lineage)"]]),
        country = normalize_country(meta_raw[["Country"]]),
        year_meta = suppressWarnings(as.numeric(meta_raw[["Isolation Year"]])),
        stringsAsFactors = FALSE
    )

    meta$lineage[is.na(meta$lineage) | meta$lineage == ""] <- "Unknown"
    meta$country[is.na(meta$country) | meta$country == ""] <- "Unknown"
    meta$continent <- get_continent(meta$country)

    meta
}

safe_read_matrix <- function(path) {
    if (!file.exists(path)) {
        warn("PBP matrix not found: %s", path)
        return(NULL)
    }

    mat <- read.table(
        path,
        sep = "\t",
        header = TRUE,
        row.names = 1,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    mat
}

# -----------------------------
# Tree extraction
# -----------------------------

extract_ss14_and_nichols <- function(tree, meta) {
    lineage_norm <- tolower(trimws(meta$lineage))

    nichols_tips <- meta$sample[
        lineage_norm == "nichols" &
            meta$sample %in% tree$tip.label
    ]

    ss14_tips <- setdiff(tree$tip.label, nichols_tips)

    if (length(nichols_tips) == 0) {
        warn("No Nichols tips were identified from metadata.")
    }

    if (length(ss14_tips) < 3) {
        stop("Too few SS14 tips after removing Nichols.", call. = FALSE)
    }

    tree_ss14 <- drop.tip(tree, setdiff(tree$tip.label, ss14_tips))
    tree_ss14_mp <- midpoint.root(tree_ss14)

    tree_nichols_mp <- NULL
    if (length(nichols_tips) >= 3) {
        tree_nichols <- drop.tip(tree, setdiff(tree$tip.label, nichols_tips))
        tree_nichols_mp <- midpoint.root(tree_nichols)
    } else {
        warn(
            "Nichols-only subtree has fewer than 3 tips (%d); skipping Nichols midpoint tree.",
            length(nichols_tips)
        )
    }

    list(
        ss14 = tree_ss14_mp,
        nichols = tree_nichols_mp,
        nichols_tips = nichols_tips,
        ss14_tips = ss14_tips
    )
}

# -----------------------------
# Build candidate sample table
# -----------------------------

build_sample_df <- function(tree, dates, meta, mat = NULL, require_pbp = TRUE) {
    rtd_vec <- calc_rtd(tree)

    sample_df <- data.frame(
        sample = tree$tip.label,
        rtd = as.numeric(rtd_vec[tree$tip.label]),
        year = as.numeric(dates[tree$tip.label]),
        stringsAsFactors = FALSE
    )

    idx <- match(sample_df$sample, meta$sample)

    sample_df$lineage <- meta$lineage[idx]
    sample_df$country <- meta$country[idx]
    sample_df$continent <- meta$continent[idx]

    sample_df$lineage[is.na(sample_df$lineage) | sample_df$lineage == ""] <- "Unknown"
    sample_df$country[is.na(sample_df$country) | sample_df$country == ""] <- "Unknown"
    sample_df$continent[is.na(sample_df$continent) | sample_df$continent == ""] <- "Unknown"

    before <- nrow(sample_df)
    sample_df <- sample_df[
        !is.na(sample_df$year) &
            is.finite(sample_df$year) &
            !is.na(sample_df$rtd) &
            is.finite(sample_df$rtd),
    ]

    info("  Valid dated tips: %d / %d", nrow(sample_df), before)

    if (require_pbp && !is.null(mat)) {
        before_pbp <- nrow(sample_df)
        sample_df <- sample_df[sample_df$sample %in% rownames(mat), ]
        info("  PBP genotype filter: %d -> %d", before_pbp, nrow(sample_df))
    }

    sample_df
}

# -----------------------------
# Improved stratified sampling
# -----------------------------

choose_diverse_by_lineage_rtd <- function(df, n_keep) {
    if (nrow(df) == 0 || n_keep <= 0) {
        return(df[0, , drop = FALSE])
    }

    if (is.infinite(n_keep) || n_keep >= nrow(df)) {
        return(df)
    }

    df$lineage_tmp <- df$lineage
    df$lineage_tmp[is.na(df$lineage_tmp) | df$lineage_tmp == ""] <- "Unknown"

    selected <- df[0, , drop = FALSE]

    # First pass: one high-RTD representative per lineage.
    lineage_levels <- unique(df$lineage_tmp)

    lineage_best <- lapply(lineage_levels, function(lin) {
        sub <- df[df$lineage_tmp == lin, , drop = FALSE]
        sub <- sub[order(-sub$rtd, sub$year, sub$sample), , drop = FALSE]
        sub[1, , drop = FALSE]
    })

    lineage_best <- do.call(rbind, lineage_best)
    lineage_best <- lineage_best[order(-lineage_best$rtd, lineage_best$year, lineage_best$sample), , drop = FALSE]

    n_first <- min(n_keep, nrow(lineage_best))
    selected <- lineage_best[seq_len(n_first), , drop = FALSE]

    # Fill remaining slots by high RTD.
    if (nrow(selected) < n_keep) {
        remaining <- df[!(df$sample %in% selected$sample), , drop = FALSE]
        remaining <- remaining[order(-remaining$rtd, remaining$year, remaining$sample), , drop = FALSE]

        n_more <- min(n_keep - nrow(selected), nrow(remaining))

        if (n_more > 0) {
            selected <- rbind(
                selected,
                remaining[seq_len(n_more), , drop = FALSE]
            )
        }
    }

    selected$lineage_tmp <- NULL
    rownames(selected) <- NULL

    selected
}

choose_temporally_diverse <- function(df, n_keep) {
    if (nrow(df) == 0 || n_keep <= 0) {
        return(df[0, , drop = FALSE])
    }

    if (is.infinite(n_keep) || n_keep >= nrow(df)) {
        return(df)
    }

    df <- df[order(df$year, -df$rtd, df$sample), , drop = FALSE]

    # Pick evenly across sorted years, then fill by RTD if duplicated indices occur.
    idx <- unique(round(seq(1, nrow(df), length.out = n_keep)))
    selected <- df[idx, , drop = FALSE]

    if (nrow(selected) < n_keep) {
        remaining <- df[!(df$sample %in% selected$sample), , drop = FALSE]
        remaining <- remaining[order(-remaining$rtd, remaining$year, remaining$sample), , drop = FALSE]
        n_more <- min(n_keep - nrow(selected), nrow(remaining))

        if (n_more > 0) {
            selected <- rbind(
                selected,
                remaining[seq_len(n_more), , drop = FALSE]
            )
        }
    }

    rownames(selected) <- NULL
    selected
}

add_focal_carrier_rescue <- function(selected_df, candidate_df, mat, focal_loci, min_carriers) {
    if (is.null(mat) || is.na(min_carriers) || min_carriers <= 0) {
        return(selected_df)
    }

    selected_samples <- selected_df$sample

    for (locus in focal_loci) {
        if (!(locus %in% colnames(mat))) {
            next
        }

        locus_values <- mat[, locus]
        carriers <- rownames(mat)[!is.na(locus_values) & locus_values == 1]
        carriers <- intersect(carriers, candidate_df$sample)

        current <- sum(selected_samples %in% carriers)

        if (current >= min_carriers) {
            next
        }

        pool <- candidate_df[
            candidate_df$sample %in% carriers &
                !(candidate_df$sample %in% selected_samples),
            ,
            drop = FALSE
        ]

        need <- min(min_carriers - current, nrow(pool))

        if (need <= 0) {
            next
        }

        add_df <- choose_temporally_diverse(pool, need)

        if (nrow(add_df) > 0) {
            selected_df <- rbind(selected_df, add_df)
            selected_samples <- selected_df$sample
            info(
                "    Rescue %-28s: %d -> %d carriers (+%d samples)",
                locus,
                current,
                sum(selected_samples %in% carriers),
                nrow(add_df)
            )
        }
    }

    selected_df <- selected_df[!duplicated(selected_df$sample), , drop = FALSE]
    rownames(selected_df) <- NULL

    selected_df
}

make_subsample <- function(sample_df, profile_row, mat = NULL, focal_loci = character(0)) {
    profile <- profile_row$profile
    year_bin_width <- profile_row$year_bin_width
    max_expansion <- profile_row$max_per_bin_expansion
    max_recent <- profile_row$max_per_bin_recent
    min_focal_carriers <- profile_row$min_focal_carriers

    info("")
    info("========================================")
    info("Subsampling profile: %s", profile)
    info("========================================")
    info(
        "  year_bin_width=%s | max expansion=%s | max recent=%s | min focal carriers=%s",
        year_bin_width,
        max_expansion,
        max_recent,
        min_focal_carriers
    )

    df <- sample_df
    df$year_bin <- floor(df$year / year_bin_width) * year_bin_width
    df$tier <- assign_tier(df$year)

    selected_list <- list()
    counter <- 1

    for (yb in sort(unique(df$year_bin))) {
        sub_y <- df[df$year_bin == yb, , drop = FALSE]

        if (nrow(sub_y) == 0) {
            next
        }

        tier <- sub_y$tier[1]

        if (tier == "early" && KEEP_ALL_EARLY) {
            selected_list[[counter]] <- sub_y
            counter <- counter + 1

            info(
                "  %s-%s [early]: ALL %d kept",
                yb,
                yb + year_bin_width,
                nrow(sub_y)
            )

            next
        }

        max_per <- ifelse(tier == "expansion", max_expansion, max_recent)

        parts <- list()
        part_counter <- 1

        for (cont in sort(unique(sub_y$continent))) {
            sub_c <- sub_y[sub_y$continent == cont, , drop = FALSE]
            n_keep <- min(max_per, nrow(sub_c))

            chosen_c <- choose_diverse_by_lineage_rtd(sub_c, n_keep)

            parts[[part_counter]] <- chosen_c
            part_counter <- part_counter + 1
        }

        bin_df <- do.call(rbind, parts)

        selected_list[[counter]] <- bin_df
        counter <- counter + 1

        info(
            "  %s-%s [%s]: %d / %d kept (%d continents)",
            yb,
            yb + year_bin_width,
            tier,
            nrow(bin_df),
            nrow(sub_y),
            length(unique(sub_y$continent))
        )
    }

    selected_df <- do.call(rbind, selected_list)
    selected_df <- selected_df[!duplicated(selected_df$sample), , drop = FALSE]
    rownames(selected_df) <- NULL

    info("  Initial selected: %d samples", nrow(selected_df))

    if (!is.null(mat) && length(focal_loci) > 0) {
        info("  Focal PBP carrier rescue:")
        selected_df <- add_focal_carrier_rescue(
            selected_df = selected_df,
            candidate_df = df,
            mat = mat,
            focal_loci = focal_loci,
            min_carriers = min_focal_carriers
        )
    }

    selected_df <- selected_df[order(selected_df$year, selected_df$continent, selected_df$sample), ]
    rownames(selected_df) <- NULL

    info("  Final selected: %d samples", nrow(selected_df))
    info(
        "  Year range: %.1f - %.1f",
        min(selected_df$year, na.rm = TRUE),
        max(selected_df$year, na.rm = TRUE)
    )
    info(
        "  Tier counts: early=%d, expansion=%d, recent=%d",
        sum(selected_df$tier == "early"),
        sum(selected_df$tier == "expansion"),
        sum(selected_df$tier == "recent")
    )

    selected_df
}

summarize_pbp_coverage <- function(selected_df, mat, profile) {
    if (is.null(mat)) {
        return(data.frame())
    }

    out <- list()
    k <- 1

    for (locus in PBP_FOCAL) {
        if (!(locus %in% colnames(mat))) {
            next
        }

        msub <- mat[rownames(mat) %in% selected_df$sample, , drop = FALSE]
        n <- sum(msub[, locus] == 1, na.rm = TRUE)

        out[[k]] <- data.frame(
            profile = profile,
            locus = locus,
            n_carrier = n,
            n_total = nrow(msub),
            freq_percent = ifelse(nrow(msub) > 0, n / nrow(msub) * 100, NA_real_),
            stringsAsFactors = FALSE
        )
        k <- k + 1
    }

    if (length(out) == 0) {
        return(data.frame())
    }

    do.call(rbind, out)
}

# -----------------------------
# Date randomization functions
# -----------------------------

get_rtd_year <- function(tree, dates) {
    tips_with_dates <- intersect(tree$tip.label, names(dates))
    tips_with_dates <- tips_with_dates[!is.na(dates[tips_with_dates])]

    if (length(tips_with_dates) < MIN_N_TEMPORAL_TEST) {
        return(NULL)
    }

    if (length(tree$tip.label) > length(tips_with_dates)) {
        drop <- setdiff(tree$tip.label, tips_with_dates)
        tree <- drop.tip(tree, drop)
    }

    rtd_raw <- node.depth.edgelength(tree)
    rtd_vec <- rtd_raw[seq_along(tree$tip.label)]
    names(rtd_vec) <- tree$tip.label
    yr_vec <- as.numeric(dates[tree$tip.label])

    ok <- is.finite(rtd_vec) & is.finite(yr_vec)

    if (sum(ok) < MIN_N_TEMPORAL_TEST) {
        return(NULL)
    }

    list(
        rtd = rtd_vec[ok],
        year = yr_vec[ok],
        tips = tree$tip.label[ok],
        tree = tree
    )
}

calc_regression <- function(rtd, year) {
    if (length(rtd) < MIN_N_TEMPORAL_TEST) {
        return(list(
            r2 = NA_real_,
            slope = NA_real_,
            intercept = NA_real_,
            slope_se = NA_real_,
            p_slope = NA_real_
        ))
    }

    fit <- lm(rtd ~ year)
    sm <- summary(fit)

    list(
        r2 = sm$r.squared,
        slope = coef(fit)[2],
        intercept = coef(fit)[1],
        slope_se = sm$coefficients[2, 2],
        p_slope = sm$coefficients[2, 4]
    )
}

date_randomization <- function(rtd, year, n_perm = 1000, seed = 42) {
    set.seed(seed)

    reg <- calc_regression(rtd, year)
    real_r2 <- reg$r2

    if (is.na(real_r2)) {
        return(NULL)
    }

    perm_r2 <- numeric(n_perm)

    for (i in seq_len(n_perm)) {
        yr_perm <- sample(year)
        fit_p <- lm(rtd ~ yr_perm)
        perm_r2[i] <- summary(fit_p)$r.squared
    }

    list(
        real_r2 = real_r2,
        slope = reg$slope,
        intercept = reg$intercept,
        slope_se = reg$slope_se,
        p_slope = reg$p_slope,
        perm_r2 = perm_r2,
        p_value = mean(perm_r2 >= real_r2),
        perm_mean = mean(perm_r2),
        perm_sd = sd(perm_r2),
        perm_95pct = as.numeric(quantile(perm_r2, 0.95)),
        n_samples = length(rtd),
        n_perm = n_perm,
        rtd = rtd,
        year = year
    )
}

interpret_signal <- function(res) {
    if (is.null(res)) {
        return(list(signal = "INSUFFICIENT DATA", decision = "SKIP"))
    }

    p <- res$p_value
    slp <- res$slope

    if (!is.na(p) && p < 0.05 && !is.na(slp) && slp > 0) {
        return(list(
            signal = "SIGNIFICANT (p<0.05, positive slope)",
            decision = "PROCEED_TREETIME"
        ))
    }

    if (!is.na(p) && p < 0.10 && !is.na(slp) && slp > 0) {
        return(list(
            signal = "MARGINAL (p<0.10, positive slope)",
            decision = "MARGINAL_EXPLORATORY"
        ))
    }

    if (!is.na(slp) && slp <= 0) {
        return(list(
            signal = "NO SIGNAL (negative or zero slope)",
            decision = "STOP_NO_SIGNAL"
        ))
    }

    list(
        signal = "NOT SIGNIFICANT (p>=0.10)",
        decision = "STOP_NO_SIGNAL"
    )
}

run_temporal_test <- function(dataset, label, tree, dates, n_perm = 1000, seed = 42, outdir = OUTDIR) {
    info("")
    info("----------------------------------------")
    info("Temporal test: %s", label)
    info("----------------------------------------")

    data <- get_rtd_year(tree, dates)

    if (is.null(data)) {
        info("  Insufficient dated tips")
        return(list(
            dataset = dataset,
            label = label,
            result = NULL
        ))
    }

    rtd_points <- data.frame(
        sample = data$tips,
        year = data$year,
        rtd = data$rtd,
        stringsAsFactors = FALSE
    )
    rtd_points_file <- file.path(outdir, sprintf("rtd_points_%s.tsv", dataset))
    write.table(rtd_points, rtd_points_file, sep = "\t", row.names = FALSE, quote = FALSE)
    info("  [Output] %s (%d points)", rtd_points_file, nrow(rtd_points))

    info("  Running %d date permutations on n=%d tips...", n_perm, length(data$rtd))

    res <- date_randomization(
        rtd = data$rtd,
        year = data$year,
        n_perm = n_perm,
        seed = seed
    )

    interp <- interpret_signal(res)

    info("  n samples:    %d", res$n_samples)
    info("  Real R2:      %.6f", res$real_r2)
    info("  Slope:        %.3e sub/site/year", res$slope)
    info("  Perm R2 mean: %.6f", res$perm_mean)
    info("  Perm R2 95th: %.6f", res$perm_95pct)
    info("  p-value:      %.4f", res$p_value)
    info("  Signal:       %s", interp$signal)
    info("  Decision:     %s", interp$decision)

    list(
        dataset = dataset,
        label = label,
        result = res,
        signal = interp$signal,
        decision = interp$decision
    )
}

result_to_summary_row <- function(x) {
    if (is.null(x$result)) {
        return(data.frame(
            dataset = x$dataset,
            label = x$label,
            n_samples = NA_integer_,
            real_r2 = NA_real_,
            slope = NA_real_,
            slope_se = NA_real_,
            p_slope_exploratory = NA_real_,
            perm_r2_mean = NA_real_,
            perm_r2_95 = NA_real_,
            p_value = NA_real_,
            signal = "INSUFFICIENT DATA",
            decision = "SKIP",
            stringsAsFactors = FALSE
        ))
    }

    res <- x$result

    data.frame(
        dataset = x$dataset,
        label = x$label,
        n_samples = res$n_samples,
        real_r2 = round(res$real_r2, 6),
        slope = signif(res$slope, 4),
        slope_se = signif(res$slope_se, 4),
        p_slope_exploratory = signif(res$p_slope, 4),
        perm_r2_mean = round(res$perm_mean, 6),
        perm_r2_95 = round(res$perm_95pct, 6),
        p_value = round(res$p_value, 4),
        signal = x$signal,
        decision = x$decision,
        stringsAsFactors = FALSE
    )
}

plot_date_randomization <- function(results, out_pdf) {
    valid <- results[!sapply(results, function(x) is.null(x$result))]

    if (length(valid) == 0) {
        warn("No valid temporal test results to plot.")
        return(invisible(NULL))
    }

    pdf(out_pdf, width = 8, height = 8)

    for (x in valid) {
        res <- x$result

        par(mfrow = c(2, 1), mar = c(5, 4, 3, 2))

        plot(
            res$year,
            res$rtd * 1e6,
            pch = 16,
            col = rgb(0.2, 0.4, 0.8, 0.30),
            cex = 0.6,
            xlab = "Sampling year",
            ylab = "Root-to-tip distance (sub/site x1e6)",
            main = sprintf(
                "%s\nn=%d, R2=%.4f, p=%.4f",
                x$label,
                res$n_samples,
                res$real_r2,
                res$p_value
            )
        )

        if (!is.na(res$slope) && is.finite(res$slope)) {
            abline(lm(res$rtd * 1e6 ~ res$year), col = "red", lwd = 2)
        }

        mtext(
            sprintf(
                "slope=%.2e sub/site/year | decision=%s",
                res$slope,
                x$decision
            ),
            side = 3,
            line = 0,
            cex = 0.75
        )

        hist(
            res$perm_r2,
            breaks = 50,
            col = "lightgray",
            border = "white",
            xlab = "R2 after date randomization",
            main = sprintf("Permutation distribution (%d permutations)", res$n_perm),
            xlim = range(c(res$perm_r2, res$real_r2, res$perm_95pct))
        )

        abline(v = res$real_r2, col = "red", lwd = 2, lty = 1)
        abline(v = res$perm_95pct, col = "blue", lwd = 1.5, lty = 2)

        legend(
            "topright",
            bty = "n",
            lwd = c(2, 1.5),
            col = c("red", "blue"),
            legend = c(
                sprintf("Real R2=%.4f", res$real_r2),
                sprintf("Perm 95th=%.4f", res$perm_95pct)
            ),
            cex = 0.8
        )
    }

    dev.off()

    info("[Output] %s", out_pdf)
}

choose_recommended_profile <- function(summary_df) {
    sub <- summary_df[grepl("^subsample_", summary_df$dataset), , drop = FALSE]

    if (nrow(sub) == 0) {
        return(NA_character_)
    }

    sig <- sub[
        sub$decision == "PROCEED_TREETIME" &
            !is.na(sub$slope) &
            sub$slope > 0,
        ,
        drop = FALSE
    ]

    if (nrow(sig) > 0) {
        sig <- sig[order(-sig$n_samples, sig$p_value), , drop = FALSE]
        return(sub("^subsample_", "", sig$dataset[1]))
    }

    marginal <- sub[
        sub$decision == "MARGINAL_EXPLORATORY" &
            !is.na(sub$slope) &
            sub$slope > 0,
        ,
        drop = FALSE
    ]

    if (nrow(marginal) > 0) {
        marginal <- marginal[order(-marginal$n_samples, marginal$p_value), , drop = FALSE]
        return(sub("^subsample_", "", marginal$dataset[1]))
    }

    sub <- sub[order(-sub$n_samples), , drop = FALSE]
    sub("^subsample_", "", sub$dataset[1])
}

# -----------------------------
# Main
# -----------------------------

cat("========================================\n")
cat("p4_01_temporal_signal_screening\n")
cat("========================================\n")

cat("\n[Step 0] Loading data...\n")

stop_if_missing(TREE_FILE, "TREE_FILE")
tree_full <- read.tree(TREE_FILE)
write.tree(tree_full, file.path(OUTDIR, "full_tree_used.nwk"))

dates <- safe_read_dates(DATES_FILE)
meta <- safe_read_metadata(META_FILE)

mat <- NULL
if (REQUIRE_PBP) {
    mat <- safe_read_matrix(MATRIX_FILE)
    if (is.null(mat)) {
        stop("REQUIRE_PBP=TRUE but MATRIX_FILE could not be read.", call. = FALSE)
    }
    info("PBP matrix: %d samples x %d loci", nrow(mat), ncol(mat))
}

info("Full tree: %d tips", length(tree_full$tip.label))
info("Dates: %d samples with dates", length(dates))
info("Metadata: %d samples", nrow(meta))

cat("\n[Step 1] Extracting SS14 and Nichols subtrees...\n")

trees <- extract_ss14_and_nichols(tree_full, meta)

tree_ss14 <- trees$ss14
tree_nichols <- trees$nichols

info("SS14 subtree: %d tips", length(tree_ss14$tip.label))
info("Nichols tips identified: %d", length(trees$nichols_tips))

write.tree(tree_ss14, file.path(OUTDIR, "ss14_midpoint_rooted.nwk"))
info("[Output] %s", file.path(OUTDIR, "ss14_midpoint_rooted.nwk"))

if (!is.null(tree_nichols)) {
    info("Nichols subtree: %d tips", length(tree_nichols$tip.label))
    write.tree(tree_nichols, file.path(OUTDIR, "nichols_midpoint_rooted.nwk"))
    info("[Output] %s", file.path(OUTDIR, "nichols_midpoint_rooted.nwk"))
}

cat("\n[Step 2] Building SS14 candidate table...\n")

sample_df <- build_sample_df(
    tree = tree_ss14,
    dates = dates,
    meta = meta,
    mat = mat,
    require_pbp = REQUIRE_PBP
)

if (nrow(sample_df) < MIN_N_TEMPORAL_TEST) {
    stop("Too few SS14 dated samples after filtering.", call. = FALSE)
}

write.table(
    sample_df,
    file.path(OUTDIR, "ss14_candidate_samples.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)
info("[Output] %s", file.path(OUTDIR, "ss14_candidate_samples.tsv"))

cat("\n[Step 3] Creating enlarged stratified subsamples...\n")

coverage_rows <- list()
coverage_k <- 1

profile_selected_files <- data.frame(
    profile = character(0),
    n_samples = integer(0),
    tsv_file = character(0),
    tree_file = character(0),
    stringsAsFactors = FALSE
)

for (i in seq_len(nrow(PROFILES))) {
    prof <- PROFILES[i, ]

    selected_df <- make_subsample(
        sample_df = sample_df,
        profile_row = prof,
        mat = mat,
        focal_loci = PBP_FOCAL
    )

    out_tsv <- file.path(
        OUTDIR,
        paste0("subsampled_phylo_informed_", prof$profile, ".tsv")
    )

    write.table(
        selected_df[, c(
            "sample",
            "year",
            "continent",
            "country",
            "lineage",
            "tier",
            "year_bin",
            "rtd"
        )],
        out_tsv,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )

    tips_keep <- selected_df$sample
    tips_to_drop <- setdiff(tree_ss14$tip.label, tips_keep)
    tree_sub <- drop.tip(tree_ss14, tips_to_drop)

    out_tree <- file.path(
        OUTDIR,
        paste0("subsampled_tree_", prof$profile, ".nwk")
    )

    write.tree(tree_sub, out_tree)

    info("[Output] %s", out_tsv)
    info("[Output] %s (%d tips)", out_tree, length(tree_sub$tip.label))

    cov <- summarize_pbp_coverage(
        selected_df = selected_df,
        mat = mat,
        profile = prof$profile
    )

    if (nrow(cov) > 0) {
        coverage_rows[[coverage_k]] <- cov
        coverage_k <- coverage_k + 1

        cat("\n  PBP focal coverage:\n")
        print(cov)
    }

    profile_selected_files <- rbind(
        profile_selected_files,
        data.frame(
            profile = prof$profile,
            n_samples = nrow(selected_df),
            tsv_file = out_tsv,
            tree_file = out_tree,
            stringsAsFactors = FALSE
        )
    )
}

if (length(coverage_rows) > 0) {
    coverage_df <- do.call(rbind, coverage_rows)

    write.table(
        coverage_df,
        file.path(OUTDIR, "subsample_pbp_coverage.tsv"),
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )

    info("[Output] %s", file.path(OUTDIR, "subsample_pbp_coverage.tsv"))
}

write.table(
    profile_selected_files,
    file.path(OUTDIR, "subsample_profile_files.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)
info("[Output] %s", file.path(OUTDIR, "subsample_profile_files.tsv"))

cat("\n[Step 4] Running date randomization tests...\n")

temporal_results <- list()

temporal_results[[length(temporal_results) + 1]] <- run_temporal_test(
    dataset = "full",
    label = "Full dataset",
    tree = tree_full,
    dates = dates,
    n_perm = N_PERM,
    seed = SEED
)

temporal_results[[length(temporal_results) + 1]] <- run_temporal_test(
    dataset = "ss14_midpoint",
    label = "SS14 subtree (midpoint-rooted)",
    tree = tree_ss14,
    dates = dates,
    n_perm = N_PERM,
    seed = SEED
)

if (!is.null(tree_nichols)) {
    temporal_results[[length(temporal_results) + 1]] <- run_temporal_test(
        dataset = "nichols_midpoint",
        label = "Nichols subtree (midpoint-rooted)",
        tree = tree_nichols,
        dates = dates,
        n_perm = N_PERM,
        seed = SEED
    )
}

for (i in seq_len(nrow(profile_selected_files))) {
    prof <- profile_selected_files$profile[i]

    tree_sub <- read.tree(profile_selected_files$tree_file[i])
    sub_df <- read.table(
        profile_selected_files$tsv_file[i],
        sep = "\t",
        header = TRUE,
        stringsAsFactors = FALSE
    )
    dates_sub <- setNames(sub_df$year, sub_df$sample)

    temporal_results[[length(temporal_results) + 1]] <- run_temporal_test(
        dataset = paste0("subsample_", prof),
        label = paste0("Stratified SS14 subsample ", prof),
        tree = tree_sub,
        dates = dates_sub,
        n_perm = N_PERM,
        seed = SEED
    )
}

summary_df <- do.call(
    rbind,
    lapply(temporal_results, result_to_summary_row)
)

summary_file <- file.path(OUTDIR, "temporal_signal_summary.tsv")

write.table(
    summary_df,
    summary_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)

info("[Output] %s", summary_file)

plot_date_randomization(
    results = temporal_results,
    out_pdf = file.path(OUTDIR, "date_randomization_summary.pdf")
)

cat("\n[Step 5] Choosing recommended subsample profile...\n")

info("Candidate profiles: %s", paste(PROFILES$profile, collapse = ", "))

recommended_profile <- choose_recommended_profile(summary_df)

if (is.na(recommended_profile)) {
    warn("No recommended profile could be chosen.")
} else {
    info("Recommended profile: %s (selected from candidates: %s)",
         recommended_profile, paste(PROFILES$profile, collapse = ", "))

    chosen <- profile_selected_files[
        profile_selected_files$profile == recommended_profile,
        ,
        drop = FALSE
    ]

    if (nrow(chosen) == 1) {
        generic_tsv <- file.path(OUTDIR, "recommended_temporal_subset.tsv")
        generic_tree <- file.path(OUTDIR, "recommended_temporal_tree.nwk")

        file.copy(chosen$tsv_file, generic_tsv, overwrite = TRUE)
        file.copy(chosen$tree_file, generic_tree, overwrite = TRUE)

        info("[Output compatibility] %s", generic_tsv)
        info("[Output] %s", generic_tree)

        writeLines(recommended_profile, file.path(OUTDIR, "recommended_profile.txt"))
        info("[Output] %s", file.path(OUTDIR, "recommended_profile.txt"))
    }
}

cat("\n========================================\n")
cat("FINAL TEMPORAL SIGNAL SUMMARY\n")
cat("========================================\n")
print(summary_df)

cat("
Done.
")
