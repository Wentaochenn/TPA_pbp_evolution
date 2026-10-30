#!/usr/bin/env Rscript
# Reconstruct focal PBP acquisition dates from one clock analysis.
# Called by clock_histories/run.py with a sampled posterior and explicit inputs.

suppressPackageStartupMessages({
    library(ape)
    library(phytools)
})
source(Sys.getenv("PBP_MK_HELPER"))

# =========================================================
# 0. Configuration
# =========================================================

MATRIX_FILE <- Sys.getenv("PBP_MATRIX_FILE", unset = "results/05_extracted_variants/sample_locus_binary_matrix.tsv")
META_FILE   <- "TPA_master_metadata.xlsx"

MCC_TREE_CANDIDATES <- c(
    Sys.getenv("PBP_MCC_TREE_FILE", unset = NA_character_),
    "results/13_beast2_selected_subset/run/strict/mcc_strict.tree",
    "results/13_beast2_selected_subset/run/strict/mcc_strict.tree",
    "results/13_beast2_selected_subset/run/strict/mcc.tree",
    "results/13_beast2_selected_subset/run/strict/mcc_strict.tree",
    "results/11_temporal_signal/beast2_run/strict/mcc_strict.tree",
    "results/11_temporal_signal/beast2_run/strict/mcc_strict.tree"
)

POSTERIOR_TREE_CANDIDATES <- c(
    Sys.getenv("PBP_POSTERIOR_TREES_FILE", unset = NA_character_),
    "results/13_beast2_selected_subset/run/strict/*.trees",
    "results/13_beast2_selected_subset/run/strict/*.treefile",
    "results/13_beast2_selected_subset/run/*.trees",
    "results/13_beast2_selected_subset/run/*.treefile",
    "results/11_temporal_signal/beast2_run/strict/*.trees",
    "results/11_temporal_signal/beast2_run/strict/*.treefile"
)

DATES_FILE_CANDIDATES <- c(
    Sys.getenv("PBP_DATES_FILE",unset=NA_character_),
    "results/13_beast2_selected_subset/beast2_input/beast2_dates.tsv",
    "results/13_beast2_selected_subset/input/beast2_dates.tsv",
    "results/13_beast2_selected_subset/run/beast2_dates.tsv",
    "results/11_temporal_signal/beast2_input/beast2_dates.tsv",
    "results/11_temporal_signal/beast2_input/beast2_dates.tsv"
)

OUTDIR <- Sys.getenv("PBP_OUTDIR", unset = "results/13_beast2_selected_subset/pbp_temporal")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

FOCAL_LOCI <- c(
    "mrcA:p.Ala506Thr",
    "mrcA:p.Ala506Val",
    "mrcA:p.Gly708Ser",
    "mrcA:p.Met625Val",
    "pbp1:p.Pro564Leu",
    "tp47:p.Ser394Arg"
)

FOCAL_LABELS <- c(
    "mrcA:p.Ala506Thr" = "A506T",
    "mrcA:p.Ala506Val" = "A506V",
    "mrcA:p.Gly708Ser" = "G708S",
    "mrcA:p.Met625Val" = "M625V",
    "pbp1:p.Pro564Leu" = "P564L",
    "tp47:p.Ser394Arg" = "S394R"
)

ROOT_STATE <- 0L
ROOT_PI <- c("0" = 0.999, "1" = 0.001)
PENICILLIN_YEAR <- 1943

# Primary dated origins exclude terminal singleton gains, as in v3.
COUNT_TERMINAL_GAINS_AS_ORIGINS <- FALSE

# Posterior HPD settings.
POSTERIOR_BURNIN <- as.numeric(Sys.getenv("PBP_POSTERIOR_BURNIN", unset = "0.10"))
POSTERIOR_THIN <- as.integer(Sys.getenv("PBP_POSTERIOR_THIN", unset = "10"))
POSTERIOR_SIMMAPS_PER_TREE <- as.integer(Sys.getenv("PBP_SIMMAPS_PER_TREE", unset = "1"))
MAX_POSTERIOR_TREES <- as.integer(Sys.getenv("PBP_MAX_POSTERIOR_TREES", unset = "500"))
MATCH_JACCARD_THRESHOLD <- as.numeric(Sys.getenv("PBP_MATCH_JACCARD", unset = "0.50"))
MIN_POSTERIOR_MATCH_SUPPORT <- as.numeric(Sys.getenv("PBP_MIN_MATCH_SUPPORT", unset = "0.50"))
STOP_IF_NO_POSTERIOR <- toupper(Sys.getenv("PBP_STOP_IF_NO_POSTERIOR", unset = "TRUE")) %in% c("TRUE", "T", "1", "YES", "Y")
MIN_CARRIERS_FOR_SIMMAP <- 3L
MIN_WT_FOR_SIMMAP <- 3L
SIMMAP_MODEL <- "ER"
SEED <- 20250616L
set.seed(SEED)

# =========================================================
# 1. Helper functions
# =========================================================

info <- function(...) cat(sprintf(...), "\n")

empty_df <- function() data.frame(stringsAsFactors = FALSE)

to_num <- function(x) suppressWarnings(as.numeric(x))

read_table_auto <- function(path) {
    read.table(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
               check.names = FALSE, quote = "", comment.char = "")
}

first_existing <- function(paths) {
    paths <- paths[!is.na(paths) & nzchar(paths)]
    hit <- paths[file.exists(paths)]
    if (length(hit) == 0) return(NA_character_)
    hit[1]
}

expand_candidate_files <- function(patterns, exclude_mcc = FALSE) {
    patterns <- patterns[!is.na(patterns) & nzchar(patterns)]
    out <- character()
    for (p in patterns) {
        if (grepl("[*?]", p)) {
            out <- c(out, Sys.glob(p))
        } else {
            out <- c(out, p)
        }
    }
    out <- unique(out[file.exists(out)])
    if (exclude_mcc) {
        out <- out[!grepl("mcc", basename(out), ignore.case = TRUE)]
    }
    out
}

find_col <- function(df, candidates, required = FALSE, label = "column") {
    nms <- names(df)
    low <- tolower(nms)
    cand <- tolower(candidates)
    exact <- match(cand, low)
    exact <- exact[!is.na(exact)]
    if (length(exact) > 0) return(nms[exact[1]])
    for (cc in cand) {
        hit <- which(grepl(cc, low, fixed = TRUE))
        if (length(hit) > 0) return(nms[hit[1]])
    }
    if (required) {
        stop(sprintf("Cannot find %s. Tried: %s. Available: %s",
                     label, paste(candidates, collapse = ", "), paste(nms, collapse = ", ")),
             call. = FALSE)
    }
    NA_character_
}

short_label <- function(locus) {
    y <- sub("^.*:p\\.", "", locus)
    aa3 <- c(Ala="A", Arg="R", Asn="N", Asp="D", Cys="C", Gln="Q", Glu="E",
             Gly="G", His="H", Ile="I", Leu="L", Lys="K", Met="M", Phe="F",
             Pro="P", Ser="S", Thr="T", Trp="W", Tyr="Y", Val="V")
    for (nm in names(aa3)) y <- gsub(nm, aa3[[nm]], y, fixed = TRUE)
    y
}

locus_label <- function(locus) {
    if (locus %in% names(FOCAL_LABELS)) return(unname(FOCAL_LABELS[[locus]]))
    short_label(locus)
}

hpd_interval <- function(x, prob = 0.95) {
    x <- sort(x[is.finite(x)])
    n <- length(x)
    if (n < 2) return(c(low = NA_real_, high = NA_real_, width = NA_real_))
    k <- ceiling(prob * n)
    if (k >= n) return(c(low = x[1], high = x[n], width = x[n] - x[1]))
    starts <- seq_len(n - k + 1)
    ends <- starts + k - 1
    widths <- x[ends] - x[starts]
    j <- which.min(widths)
    c(low = x[starts[j]], high = x[ends[j]], width = widths[j])
}

summarise_dates <- function(x, threshold = PENICILLIN_YEAR, prob = 0.95) {
    x <- x[is.finite(x)]
    n <- length(x)
    if (n == 0) {
        return(data.frame(n = 0L, median = NA_real_, mean = NA_real_,
                          q025 = NA_real_, q975 = NA_real_,
                          hpd_low = NA_real_, hpd_high = NA_real_, hpd_width = NA_real_,
                          p_after_1943 = NA_real_, stringsAsFactors = FALSE))
    }
    hp <- hpd_interval(x, prob = prob)
    data.frame(
        n = n,
        median = median(x),
        mean = mean(x),
        q025 = unname(quantile(x, 0.025, names = FALSE, type = 8)),
        q975 = unname(quantile(x, 0.975, names = FALSE, type = 8)),
        hpd_low = unname(hp["low"]),
        hpd_high = unname(hp["high"]),
        hpd_width = unname(hp["width"]),
        p_after_1943 = mean(x > threshold),
        stringsAsFactors = FALSE
    )
}

split_set <- function(x) {
    if (length(x) == 0 || is.na(x) || !nzchar(x)) return(character())
    unique(strsplit(x, ";", fixed = TRUE)[[1]])
}

collapse_set <- function(x) {
    x <- sort(unique(x[!is.na(x) & nzchar(x)]))
    paste(x, collapse = ";")
}

jaccard_set <- function(a, b) {
    a <- unique(a); b <- unique(b)
    if (length(a) == 0 && length(b) == 0) return(1)
    u <- union(a, b)
    if (length(u) == 0) return(0)
    length(intersect(a, b)) / length(u)
}

read_tree_any <- function(paths, allow_glob = TRUE, exclude_mcc = FALSE, label = "tree") {
    files <- if (allow_glob) expand_candidate_files(paths, exclude_mcc = exclude_mcc) else paths[file.exists(paths)]
    if (length(files) == 0) {
        stop("No ", label, " found. Checked:\n", paste(paths, collapse = "\n"), call. = FALSE)
    }
    f <- files[1]
    tr <- tryCatch(read.nexus(f), error = function(e) NULL)
    if (is.null(tr)) tr <- tryCatch(read.tree(f), error = function(e) NULL)
    if (is.null(tr)) stop("Could not read ", label, ": ", f, call. = FALSE)
    if (inherits(tr, "multiPhylo")) tr <- tr[[1]]
    if (is.null(tr$edge.length)) stop("Tree has no edge lengths: ", f, call. = FALSE)
    if(any(!is.finite(tr$edge.length))||any(tr$edge.length < -1e-10))stop("Invalid MCC: use common-ancestor heights")
    tr$edge.length[tr$edge.length<=0]<-1e-10
    tr <- reorder.phylo(tr, "cladewise")
    list(file = f, tree = tr)
}

read_posterior_trees <- function(paths) {
    if(POSTERIOR_BURNIN!=0 || POSTERIOR_THIN!=1)stop("Use externally sampled posterior file with burnin=0,thin=1")
    files <- expand_candidate_files(paths, exclude_mcc = TRUE)
    if (length(files) == 0) return(list(file = NA_character_, trees = NULL))
    f <- files[1]
    info("Reading posterior trees: %s", f)
    trs <- tryCatch(read.nexus(f), error = function(e) NULL)
    if (is.null(trs)) trs <- tryCatch(read.tree(f), error = function(e) NULL)
    if (is.null(trs)) stop("Could not read posterior trees: ", f, call. = FALSE)
    if (!inherits(trs, "multiPhylo")) trs <- structure(list(trs), class = "multiPhylo")
    list(file = f, trees = trs)
}

read_tip_dates <- function() {
    f <- first_existing(DATES_FILE_CANDIDATES)
    if (!is.na(f)) {
        df <- read_table_auto(f)
        tax_col <- find_col(df, c("taxon", "sample", "sample_id", "strain", "name"), TRUE, "date taxon column")
        date_col <- find_col(df, c("date", "year", "decimal_date", "decimal_year"), TRUE, "date column")
        dates <- setNames(to_num(df[[date_col]]), as.character(df[[tax_col]]))
        return(list(file = f, dates = dates))
    }
    if (file.exists(META_FILE) && requireNamespace("readxl", quietly = TRUE)) {
        df <- as.data.frame(readxl::read_excel(META_FILE), stringsAsFactors = FALSE)
        tax_col <- find_col(df, c("SRA_ID (Key)", "sample", "sample_id", "taxon", "strain"), TRUE, "metadata sample column")
        date_col <- find_col(df, c("Isolation Year", "Sample_Year", "year", "sampling_year"), TRUE, "metadata year column")
        dates <- setNames(to_num(df[[date_col]]), as.character(df[[tax_col]]))
        return(list(file = META_FILE, dates = dates))
    }
    stop("No dates file found and metadata fallback failed.", call. = FALSE)
}

prepare_tree_inputs <- function(tree, mat, tip_dates_all, min_date_fraction = 0.50) {
    common <- intersect(tree$tip.label, rownames(mat))
    if(length(common)!=length(tree$tip.label)||anyDuplicated(tree$tip.label))
        stop("Posterior/MCC sample set does not match the matrix")
    tip_dates <- tip_dates_all[tree$tip.label]
    date_ok <- is.finite(tip_dates)
    if (sum(date_ok) != length(tree$tip.label)) {
        stop(sprintf("Too few dated tree tips: %d/%d", sum(date_ok), length(tree$tip.label)), call. = FALSE)
    }
    tree$edge.length[!is.finite(tree$edge.length) | tree$edge.length <= 0] <- 1e-8
    tree <- reorder.phylo(tree, "cladewise")
    mat2 <- mat[tree$tip.label, , drop = FALSE]
    list(tree = tree, mat = mat2, tip_dates = tip_dates)
}

estimate_node_dates <- function(tree, tip_dates) {
    n_tip <- length(tree$tip.label)
    n_total <- n_tip + tree$Nnode
    depths <- node.depth.edgelength(tree)
    date_ok <- is.finite(tip_dates)
    root_candidates <- tip_dates[date_ok] - depths[seq_len(n_tip)][date_ok]
    root_date <- median(root_candidates, na.rm = TRUE)
    node_dates <- root_date + depths[seq_len(n_total)]
    list(root_date = root_date,
         root_date_min = min(root_candidates, na.rm = TRUE),
         root_date_max = max(root_candidates, na.rm = TRUE),
         root_date_mad = mad(root_candidates, na.rm = TRUE),
         node_dates = node_dates)
}

get_descendant_nodes <- function(tree, node) {
    edges <- tree$edge
    desc <- integer(0)
    queue <- as.integer(node)
    while (length(queue) > 0) {
        cur <- queue[1]
        queue <- queue[-1]
        kids <- edges[edges[, 1] == cur, 2]
        if (length(kids) > 0) {
            desc <- c(desc, kids)
            queue <- c(queue, kids)
        }
    }
    unique(desc)
}

get_descendant_tips <- function(tree, node) {
    n_tip <- length(tree$tip.label)
    if (node <= n_tip) return(tree$tip.label[node])
    dd <- get_descendant_nodes(tree, node)
    tip_nodes <- dd[dd <= n_tip]
    tree$tip.label[tip_nodes]
}

fitch_rooted_binary <- function(tree, tip_state_01, root_state = 0L) {
    n_tip <- length(tree$tip.label)
    n_total <- n_tip + tree$Nnode
    depths <- node.depth.edgelength(tree)
    state_sets <- vector("list", n_total)
    for (i in seq_len(n_tip)) {
        st <- tip_state_01[i]
        if (is.na(st)) state_sets[[i]] <- c(0L, 1L) else state_sets[[i]] <- as.integer(st)
    }
    internal_nodes <- (n_tip + 1L):n_total
    internal_nodes <- internal_nodes[order(depths[internal_nodes], decreasing = TRUE)]
    for (node in internal_nodes) {
        kids <- tree$edge[tree$edge[, 1] == node, 2]
        if (length(kids) == 0) {
            state_sets[[node]] <- c(0L, 1L)
            next
        }
        child_sets <- state_sets[kids]
        inter <- Reduce(intersect, child_sets)
        if (length(inter) > 0) state_sets[[node]] <- sort(as.integer(inter))
        else state_sets[[node]] <- sort(unique(as.integer(unlist(child_sets))))
    }
    node_state <- rep(NA_integer_, n_total)
    root_node <- n_tip + 1L
    if (root_state %in% state_sets[[root_node]]) node_state[root_node] <- as.integer(root_state)
    else node_state[root_node] <- as.integer(state_sets[[root_node]][1])
    edge_order <- order(depths[tree$edge[, 1]], depths[tree$edge[, 2]])
    for (edge_i in edge_order) {
        parent <- tree$edge[edge_i, 1]
        child <- tree$edge[edge_i, 2]
        if (is.na(node_state[parent])) next
        child_set <- state_sets[[child]]
        if (node_state[parent] %in% child_set) node_state[child] <- node_state[parent]
        else node_state[child] <- as.integer(child_set[1])
    }
    list(node_state = node_state, state_sets = state_sets)
}

classify_temporal_mode <- function(n_gain_primary, freq) {
    if (is.na(n_gain_primary)) return("unresolved")
    if (n_gain_primary <= 1) return("single_origin_clonal_expansion")
    if (n_gain_primary >= 5) return("multiple_independent_origins")
    if (n_gain_primary >= 2 && n_gain_primary <= 4) return("oligoclonal_or_intermediate")
    if (freq < 0.02) return("rare_or_underpowered")
    "intermediate_or_unresolved"
}

infer_mcc_parsimony_events <- function(tree, mat, locus, tip_dates, root_state = ROOT_STATE) {
    n_tip <- length(tree$tip.label)
    n_total <- n_tip + tree$Nnode
    nd <- estimate_node_dates(tree, tip_dates)
    node_dates <- nd$node_dates
    tip_state <- to_num(mat[tree$tip.label, locus])
    tip_state[!(tip_state %in% c(0, 1))] <- NA_real_
    n_known <- sum(!is.na(tip_state))
    n_carriers <- sum(tip_state == 1, na.rm = TRUE)
    n_wt <- sum(tip_state == 0, na.rm = TRUE)
    freq <- ifelse(n_known > 0, n_carriers / n_known, NA_real_)

    rec <- fitch_rooted_binary(tree, as.integer(tip_state), root_state = root_state)
    node_state <- rec$node_state
    rows <- list()
    gain_idx <- 0L
    loss_idx <- 0L

    for (edge_i in seq_len(nrow(tree$edge))) {
        parent <- tree$edge[edge_i, 1]
        child <- tree$edge[edge_i, 2]
        ps <- node_state[parent]
        cs <- node_state[child]
        if (is.na(ps) || is.na(cs) || ps == cs) next
        event_type <- ifelse(ps == 0 && cs == 1, "gain", "loss")
        if (event_type == "gain") gain_idx <- gain_idx + 1L else loss_idx <- loss_idx + 1L
        parent_date <- node_dates[parent]
        child_date <- node_dates[child]
        event_date <- mean(c(parent_date, child_date), na.rm = TRUE)
        desc_tips <- get_descendant_tips(tree, child)
        desc_states <- tip_state[match(desc_tips, tree$tip.label)]
        carrier_tips <- desc_tips[which(desc_states == 1)]
        rows[[length(rows) + 1L]] <- data.frame(
            locus = locus,
            label = locus_label(locus),
            event = event_type,
            independent_origin_index = ifelse(event_type == "gain", gain_idx, NA_integer_),
            parent_node = parent,
            node = child,
            child_node = child,
            parent_state = ps,
            child_state = cs,
            parent_date = parent_date,
            child_date = child_date,
            event_date = event_date,
            origin_lower = min(parent_date, child_date),
            origin_upper = max(parent_date, child_date),
            branch_length_years = abs(child_date - parent_date),
            n_descendants = length(desc_tips),
            n_carrier_descendants = length(carrier_tips),
            n_known_descendants = sum(!is.na(desc_states)),
            is_terminal = child <= n_tip,
            origin_type = ifelse(child <= n_tip, "terminal_singleton", "internal_clade"),
            carrier_tip_set = collapse_set(carrier_tips),
            support = "parsimony_root0",
            method = "root_noncarrier_binary_parsimony_on_CA_MCC",
            stringsAsFactors = FALSE
        )
    }
    events <- if (length(rows) > 0) do.call(rbind, rows) else empty_df()
    if (nrow(events) > 0) {
        gains <- events[events$event == "gain", , drop = FALSE]
        if (nrow(gains) > 0) {
            gains <- gains[order(gains$event_date, gains$node), , drop = FALSE]
            gains$independent_origin_index <- seq_len(nrow(gains))
            events[events$event == "gain", ] <- gains
        }
    }
    summary <- data.frame(
        locus = locus,
        label = locus_label(locus),
        n_known = n_known,
        n_carriers = n_carriers,
        n_wt = n_wt,
        freq_pct = round(freq * 100, 4),
        stringsAsFactors = FALSE
    )
    list(events = events, summary = summary, node_dates = node_dates, root_calibration = nd)
}

fit_q_matrix <- function(tree, tip_state, model = SIMMAP_MODEL) {
    known <- !is.na(tip_state)
    tr <- drop.tip(tree,tree$tip.label[!known])
    states <- setNames(tip_state[known],tree$tip.label[known])[tr$tip.label]
    if(length(unique(states))<2L)return(NULL)
    tr$edge.length[tr$edge.length<=0] <- 1e-10
    fit <- pbp_fit_mk(tr,states,rate_model="ER",root_prior="flat")
    if(!pbp_mk_fit_ok(fit))stop("Temporal ER fitting failed: ",fit$error)
    Q <- fit$transition_matrix;dimnames(Q)<-list(c("0","1"),c("0","1"));Q
}

extract_simmap_gain_events <- function(sm, locus, tip_state_named, node_dates,
                                       posterior_tree_index, simmap_index,
                                       internal_only = FALSE) {
    n_tip <- length(sm$tip.label)
    rows <- list()
    for (edge_i in seq_len(nrow(sm$edge))) {
        parent <- sm$edge[edge_i, 1]
        child <- sm$edge[edge_i, 2]
        is_terminal <- child <= n_tip
        if (internal_only && is_terminal) next
        maps <- sm$maps[[edge_i]]
        if (is.null(maps) || length(maps) < 2) next
        states <- names(maps)
        dur <- as.numeric(maps)
        if (any(!is.finite(dur)) || sum(dur) <= 0) next
        cum <- cumsum(dur)
        parent_date <- node_dates[parent]
        child_date <- node_dates[child]
        for (k in seq_len(length(states) - 1L)) {
            if (states[k] == "0" && states[k + 1L] == "1") {
                prop <- cum[k] / sum(dur)
                event_date <- parent_date + prop * (child_date - parent_date)
                desc_tips <- get_descendant_tips(sm, child)
                desc_states <- tip_state_named[desc_tips]
                carrier_tips <- desc_tips[which(desc_states == 1)]
                rows[[length(rows) + 1L]] <- data.frame(
                    locus = locus,
                    label = locus_label(locus),
                    posterior_tree_index = posterior_tree_index,
                    simmap_index = simmap_index,
                    edge_index = edge_i,
                    parent_node = parent,
                    child_node = child,
                    parent_date = parent_date,
                    child_date = child_date,
                    event_date = event_date,
                    branch_lower = min(parent_date, child_date),
                    branch_upper = max(parent_date, child_date),
                    transition_prop_from_parent = prop,
                    n_descendants = length(desc_tips),
                    n_carrier_descendants = length(carrier_tips),
                    is_terminal = is_terminal,
                    origin_type = ifelse(is_terminal, "terminal_singleton", "internal_clade"),
                    carrier_tip_set = collapse_set(carrier_tips),
                    method = "posterior_tree_SIMMAP_transition",
                    stringsAsFactors = FALSE
                )
            }
        }
    }
    if (length(rows) > 0) do.call(rbind, rows) else empty_df()
}

run_posterior_simmap_for_tree <- function(tree, mat, tip_dates_all, locus, Q, tree_index,
                                          nsim = POSTERIOR_SIMMAPS_PER_TREE) {
    prep <- prepare_tree_inputs(tree,mat,tip_dates_all,min_date_fraction=1)
    tree2<-prep$tree;mat2<-prep$mat
    x<-to_num(mat2[tree2$tip.label,locus]);x[!(x%in%c(0,1))]<-NA_real_
    if(sum(x==1,na.rm=TRUE)<MIN_CARRIERS_FOR_SIMMAP || sum(x==0,na.rm=TRUE)<MIN_WT_FOR_SIMMAP)
        return(list(events=empty_df(),success=0L,status="not_identifiable_tip_states"))
    if(is.null(Q))stop("Missing Q for eligible temporal locus: ",locus)
    known<-!is.na(x);named<-setNames(x[known],tree2$tip.label[known])
    tr<-drop.tip(tree2,setdiff(tree2$tip.label,names(named)));named<-named[tr$tip.label]
    if(any(tr$edge.length < -1e-10))stop("Negative posterior branch")
    tr$edge.length[tr$edge.length<=0]<-1e-10
    expected_jumps<-sum(tr$edge.length)*max(-diag(Q))
    if(!is.finite(expected_jumps)||expected_jumps>1e5)
        stop("ER rate implies excessive mapping jumps; inspect Q rather than silently cap the rate")
    dates<-prep$tip_dates[tr$tip.label]
    nd<-estimate_node_dates(tr,dates)
    if(nd$root_date_max-nd$root_date_min>1e-5)stop("Posterior tree/date calibration mismatch")
    fac<-factor(as.character(as.integer(named)),levels=c("0","1"));names(fac)<-tr$tip.label
    # Deterministic per-tree/per-locus seeds; independent of iteration failures.
    set.seed(SEED+as.integer(tree_index)*100L+match(locus,FOCAL_LOCI))
    sm<-suppressMessages(make.simmap(tr,fac,Q=Q,nsim=nsim,pi=ROOT_PI,message=FALSE))
    if(!inherits(sm,"multiSimmap"))sm<-structure(list(sm),class="multiSimmap")
    if(length(sm)!=nsim)stop("Incomplete stochastic maps")
    events<-lapply(seq_along(sm),function(i) extract_simmap_gain_events(sm[[i]],locus,
       setNames(as.integer(named),names(named)),nd$node_dates,tree_index,i,
       internal_only=!COUNT_TERMINAL_GAINS_AS_ORIGINS))
    events<-Filter(function(e)nrow(e)>0,events)
    list(events=if(length(events))do.call(rbind,events)else empty_df(),
         success=length(sm),status="complete")
}

match_posterior_events_to_mcc <- function(mcc_events, posterior_events, threshold = MATCH_JACCARD_THRESHOLD) {
    if (nrow(mcc_events) == 0 || nrow(posterior_events) == 0) return(empty_df())
    out <- list()
    # Split by locus and posterior tree/map so matching is one-to-one within each posterior realization.
    group_key <- paste(posterior_events$locus, posterior_events$posterior_tree_index,
                       posterior_events$simmap_index, sep = "||")
    groups <- split(posterior_events, group_key)
    for (gname in names(groups)) {
        pe <- groups[[gname]]
        locus <- unique(pe$locus)
        if (length(locus) != 1) next
        me <- mcc_events[mcc_events$locus == locus, , drop = FALSE]
        if (nrow(me) == 0 || nrow(pe) == 0) next
        pairs <- list()
        idx <- 0L
        for (i in seq_len(nrow(me))) {
            mset <- split_set(me$carrier_tip_set[i])
            for (j in seq_len(nrow(pe))) {
                pset <- split_set(pe$carrier_tip_set[j])
                jac <- jaccard_set(mset, pset)
                if (is.finite(jac) && jac >= threshold) {
                    idx <- idx + 1L
                    pairs[[idx]] <- data.frame(mcc_i = i, post_j = j, jaccard = jac, stringsAsFactors = FALSE)
                }
            }
        }
        if (length(pairs) == 0) next
        pair_df <- do.call(rbind, pairs)
        pair_df <- pair_df[order(-pair_df$jaccard), , drop = FALSE]
        used_m <- integer(0); used_p <- integer(0)
        for (r in seq_len(nrow(pair_df))) {
            mi <- pair_df$mcc_i[r]
            pj <- pair_df$post_j[r]
            if (mi %in% used_m || pj %in% used_p) next
            used_m <- c(used_m, mi); used_p <- c(used_p, pj)
            mrow <- me[mi, , drop = FALSE]
            prow <- pe[pj, , drop = FALSE]
            out[[length(out) + 1L]] <- data.frame(
                event_uid = mrow$event_uid,
                locus = mrow$locus,
                label = mrow$label,
                mcc_independent_origin_index = mrow$independent_origin_index,
                mcc_event_date = mrow$event_date,
                mcc_origin_lower = mrow$origin_lower,
                mcc_origin_upper = mrow$origin_upper,
                posterior_tree_index = prow$posterior_tree_index,
                simmap_index = prow$simmap_index,
                posterior_event_date = prow$event_date,
                posterior_branch_lower = prow$branch_lower,
                posterior_branch_upper = prow$branch_upper,
                posterior_transition_prop_from_parent = prow$transition_prop_from_parent,
                posterior_n_carrier_descendants = prow$n_carrier_descendants,
                posterior_origin_type = prow$origin_type,
                jaccard = pair_df$jaccard[r],
                mcc_carrier_tip_set = mrow$carrier_tip_set,
                posterior_carrier_tip_set = prow$carrier_tip_set,
                stringsAsFactors = FALSE
            )
        }
    }
    if (length(out) > 0) do.call(rbind, out) else empty_df()
}

safe_merge <- function(x, y, by, all.x = TRUE) {
    if (nrow(y) == 0) return(x)
    merge(x, y, by = by, all.x = all.x, sort = FALSE)
}

# =========================================================
# 2. Load MCC tree, dates, and matrix
# =========================================================

info("========================================")
info("PBP temporal origin analysis: posterior SIMMAP HPD")
info("========================================")
info("Seed: %d", SEED)
info("Output directory: %s", OUTDIR)

if (!file.exists(MATRIX_FILE)) stop("Missing MATRIX_FILE: ", MATRIX_FILE, call. = FALSE)
mat <- read.table(MATRIX_FILE, sep = "\t", header = TRUE, row.names = 1,
                  check.names = FALSE, stringsAsFactors = FALSE)
for (j in seq_len(ncol(mat))) mat[[j]] <- to_num(mat[[j]])

date_obj <- read_tip_dates()
tip_dates_all <- date_obj$dates
info("Date source: %s", date_obj$file)

mcc_obj <- read_tree_any(MCC_TREE_CANDIDATES, allow_glob = TRUE, exclude_mcc = FALSE, label = "MCC tree")
info("MCC tree: %s", mcc_obj$file)
prep_mcc <- prepare_tree_inputs(mcc_obj$tree, mat, tip_dates_all)
mcc_tree <- prep_mcc$tree
mcc_mat <- prep_mcc$mat
mcc_tip_dates <- prep_mcc$tip_dates
mcc_cal <- estimate_node_dates(mcc_tree, mcc_tip_dates)
info("MCC tree after filtering: %d tips, %d internal nodes", length(mcc_tree$tip.label), mcc_tree$Nnode)
info("MCC root-date median %.2f; root-date MAD %.3f", mcc_cal$root_date, mcc_cal$root_date_mad)
if(mcc_cal$root_date_max-mcc_cal$root_date_min>1e-5)stop("MCC dates not chronologically calibrated")

calib <- data.frame(
    mcc_tree_file = mcc_obj$file,
    date_source = date_obj$file,
    n_mcc_tips = length(mcc_tree$tip.label),
    n_mcc_tips_with_dates = sum(is.finite(mcc_tip_dates)),
    mcc_root_date_median = mcc_cal$root_date,
    mcc_root_date_min = mcc_cal$root_date_min,
    mcc_root_date_max = mcc_cal$root_date_max,
    mcc_root_date_mad = mcc_cal$root_date_mad,
    latest_tip_date = max(mcc_tip_dates, na.rm = TRUE),
    earliest_tip_date = min(mcc_tip_dates, na.rm = TRUE),
    stringsAsFactors = FALSE
)
write.table(calib, file.path(OUTDIR, "pbp_temporal_tree_calibration.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# =========================================================
# 3. MCC parsimony event dating and Q matrices
# =========================================================

info("")
info("========================================")
info("MCC parsimony gains and locus-level summaries")
info("========================================")

mcc_event_rows <- list()
summary_rows <- list()
Q_BY_LOCUS <- list()

for (locus in FOCAL_LOCI) {
    if (!locus %in% colnames(mcc_mat)) {
        info("[SKIP] %s not found in matrix", locus)
        next
    }
    res <- infer_mcc_parsimony_events(mcc_tree, mcc_mat, locus, mcc_tip_dates, ROOT_STATE)
    ev <- res$events
    sumrow <- res$summary
    if (nrow(ev) > 0) mcc_event_rows[[locus]] <- ev

    gains <- if (nrow(ev) > 0) ev[ev$event == "gain", , drop = FALSE] else empty_df()
    n_gain <- nrow(gains)
    n_internal <- if (n_gain > 0) sum(!gains$is_terminal) else 0L
    n_terminal <- if (n_gain > 0) sum(gains$is_terminal) else 0L
    n_primary <- if (COUNT_TERMINAL_GAINS_AS_ORIGINS) n_gain else n_internal
    freq <- sumrow$freq_pct / 100

    sumrow$n_gain_events <- n_gain
    sumrow$n_internal_gain_events <- n_internal
    sumrow$n_terminal_gain_events <- n_terminal
    sumrow$n_gain_primary <- n_primary
    sumrow$earliest_mcc_gain_date <- if (n_gain > 0) min(gains$event_date, na.rm = TRUE) else NA_real_
    sumrow$latest_mcc_gain_date <- if (n_gain > 0) max(gains$event_date, na.rm = TRUE) else NA_real_
    sumrow$temporal_mode <- if(sumrow$n_wt==0L) "all_known_carriers_origin_unresolved" else if(n_primary==0L) "no_reconstructed_gain" else classify_temporal_mode(n_primary, freq)
    sumrow$method <- "CA_MCC_root_noncarrier_parsimony_plus_posterior_SIMMAP_HPD"
    summary_rows[[locus]] <- sumrow

    tip_state <- to_num(mcc_mat[mcc_tree$tip.label, locus])
    tip_state[!(tip_state %in% c(0, 1))] <- NA_real_
    Q <- fit_q_matrix(mcc_tree, tip_state, SIMMAP_MODEL)
    Q_BY_LOCUS[[locus]] <- Q

    info("%s: carriers %d/%d (%.2f%%); MCC gains %d total (%d internal, %d terminal); Q=%s",
         locus_label(locus), sumrow$n_carriers, sumrow$n_known, sumrow$freq_pct,
         n_gain, n_internal, n_terminal, ifelse(is.null(Q), "NOT_IDENTIFIABLE", "OK"))
}

mcc_events_all <- if (length(mcc_event_rows) > 0) do.call(rbind, mcc_event_rows) else empty_df()
summary_df <- if (length(summary_rows) > 0) do.call(rbind, summary_rows) else empty_df()

# Primary events for HPD matching: gain events, excluding terminal singleton gains unless configured otherwise.
mcc_gains <- if (nrow(mcc_events_all) > 0) mcc_events_all[mcc_events_all$event == "gain", , drop = FALSE] else empty_df()
if (nrow(mcc_gains) > 0 && !COUNT_TERMINAL_GAINS_AS_ORIGINS) {
    mcc_primary_gains <- mcc_gains[!mcc_gains$is_terminal, , drop = FALSE]
} else {
    mcc_primary_gains <- mcc_gains
}

if (nrow(mcc_primary_gains) > 0) {
    mcc_primary_gains <- mcc_primary_gains[order(mcc_primary_gains$locus, mcc_primary_gains$event_date, mcc_primary_gains$node), , drop = FALSE]
    mcc_primary_gains$event_uid <- sprintf("%s_origin_%02d",
                                           gsub("[^A-Za-z0-9]+", "_", mcc_primary_gains$label),
                                           ave(mcc_primary_gains$event_date, mcc_primary_gains$locus, FUN = seq_along))
}

write.table(mcc_events_all, file.path(OUTDIR, "pbp_origin_events_all.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary_df, file.path(OUTDIR, "pbp_temporal_summary_mcc_only.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mcc_primary_gains, file.path(OUTDIR, "pbp_dated_independent_origins_mcc_only.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# =========================================================
# 4. Posterior-tree SIMMAP gain dating
# =========================================================

info("")
info("========================================")
info("Posterior-tree SIMMAP origin-date sampling")
info("========================================")

post_obj <- read_posterior_trees(POSTERIOR_TREE_CANDIDATES)
if (is.null(post_obj$trees)) {
    msg <- paste0("No posterior .trees file found. True HPD cannot be estimated. ",
                  "Set PBP_POSTERIOR_TREES_FILE=/path/to/posterior.trees or place a .trees file in the BEAST2 run directory.")
    if (STOP_IF_NO_POSTERIOR) stop(msg, call. = FALSE) else warning(msg)
}

posterior_gain_events <- empty_df()
posterior_assignments <- empty_df()
locus_map_counts <- data.frame(locus = FOCAL_LOCI, n_posterior_tree_maps_attempted = 0L,
                               n_posterior_tree_maps_with_simmap = 0L,
                               stringsAsFactors = FALSE)

if (!is.null(post_obj$trees)) {
    posterior_trees <- post_obj$trees
    n_total_post <- length(posterior_trees)
    burnin_n <- if (POSTERIOR_BURNIN < 1) floor(n_total_post * POSTERIOR_BURNIN) else as.integer(POSTERIOR_BURNIN)
    keep_idx <- seq.int(burnin_n + 1L, n_total_post, by = max(1L, POSTERIOR_THIN))
    if (length(keep_idx) > MAX_POSTERIOR_TREES) {
        keep_idx <- keep_idx[seq_len(MAX_POSTERIOR_TREES)]
    }
    posterior_trees <- posterior_trees[keep_idx]
    info("Posterior file: %s", post_obj$file)
    info("Total posterior trees: %d; burn-in removed: %d; thin: %d; retained: %d",
         n_total_post, burnin_n, POSTERIOR_THIN, length(posterior_trees))
    info("SIMMAPs per retained tree per locus: %d", POSTERIOR_SIMMAPS_PER_TREE)
    info("Event matching Jaccard threshold: %.2f", MATCH_JACCARD_THRESHOLD)

    all_post_events <- list()
    map_audit <- list()
    for (ti in seq_along(posterior_trees)) {
        if (ti %% 25 == 0 || ti == 1) info("  Processing posterior tree %d/%d", ti, length(posterior_trees))
        tr <- posterior_trees[[ti]]
        if (is.null(tr$edge.length)) next
        tr$edge.length[!is.finite(tr$edge.length) | tr$edge.length <= 0] <- 1e-8
        tr <- reorder.phylo(tr, "cladewise")
        for (locus in FOCAL_LOCI) {
            Q <- Q_BY_LOCUS[[locus]]
            locus_map_counts$n_posterior_tree_maps_attempted[locus_map_counts$locus == locus] <-
                locus_map_counts$n_posterior_tree_maps_attempted[locus_map_counts$locus == locus] + POSTERIOR_SIMMAPS_PER_TREE
            result <- run_posterior_simmap_for_tree(tr, mat, tip_dates_all, locus, Q,
                                                tree_index = keep_idx[ti],
                                                nsim = POSTERIOR_SIMMAPS_PER_TREE)
            ev <- result$events
            map_audit[[length(map_audit)+1L]] <- data.frame(locus=locus,
                posterior_tree_index=keep_idx[ti],status=result$status,
                successful_maps=result$success,internal_gain_events=nrow(ev))
            locus_map_counts$n_posterior_tree_maps_with_simmap[locus_map_counts$locus == locus] <-
                locus_map_counts$n_posterior_tree_maps_with_simmap[locus_map_counts$locus == locus] + result$success
            if(nrow(ev)>0)all_post_events[[length(all_post_events)+1L]]<-ev
        }
    }
    write.table(do.call(rbind,map_audit),file.path(OUTDIR,"pbp_mapping_execution_audit.tsv"),
                sep="\t",quote=FALSE,row.names=FALSE)
    posterior_gain_events <- if (length(all_post_events) > 0) do.call(rbind, all_post_events) else empty_df()

    write.table(posterior_gain_events, file.path(OUTDIR, "pbp_posterior_simmap_gain_events.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(locus_map_counts, file.path(OUTDIR, "pbp_posterior_simmap_map_counts.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    info("Posterior SIMMAP gain transitions extracted: %d", nrow(posterior_gain_events))

    if (nrow(mcc_primary_gains) > 0 && nrow(posterior_gain_events) > 0) {
        posterior_assignments <- match_posterior_events_to_mcc(mcc_primary_gains, posterior_gain_events,
                                                               threshold = MATCH_JACCARD_THRESHOLD)
    }
    write.table(posterior_assignments, file.path(OUTDIR, "pbp_posterior_simmap_origin_assignments.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    info("Posterior SIMMAP transitions matched to MCC origin events: %d", nrow(posterior_assignments))
}

# =========================================================
# 5. Summarise posterior HPD intervals
# =========================================================

info("")
info("========================================")
info("Summarising origin-date posterior HPD intervals")
info("========================================")

origin_hpd_rows <- list()
if (nrow(mcc_primary_gains) > 0) {
    for (i in seq_len(nrow(mcc_primary_gains))) {
        ev <- mcc_primary_gains[i, , drop = FALSE]
        ass <- if (nrow(posterior_assignments) > 0) {
            posterior_assignments[posterior_assignments$event_uid == ev$event_uid, , drop = FALSE]
        } else empty_df()
        s <- summarise_dates(if (nrow(ass) > 0) ass$posterior_event_date else numeric(0),
                             threshold = PENICILLIN_YEAR)
        n_possible <- locus_map_counts$n_posterior_tree_maps_attempted[locus_map_counts$locus == ev$locus]
        if (length(n_possible) == 0) n_possible <- NA_integer_
        n_possible_success <- locus_map_counts$n_posterior_tree_maps_with_simmap[locus_map_counts$locus == ev$locus]
        if (length(n_possible_success) == 0) n_possible_success <- NA_integer_
        origin_hpd_rows[[length(origin_hpd_rows) + 1L]] <- data.frame(
            event_uid = ev$event_uid,
            locus = ev$locus,
            label = ev$label,
            mcc_independent_origin_index = ev$independent_origin_index,
            mcc_event_date = ev$event_date,
            mcc_branch_lower = ev$origin_lower,
            mcc_branch_upper = ev$origin_upper,
            mcc_branch_width = ev$origin_upper - ev$origin_lower,
            mcc_n_carrier_descendants = ev$n_carrier_descendants,
            mcc_origin_type = ev$origin_type,
            posterior_n_assigned = s$n,
            posterior_n_possible_tree_maps = n_possible,
            posterior_n_tree_maps_with_any_simmap_gain = n_possible_success,
            posterior_event_matching_support = ifelse(is.finite(n_possible) && n_possible > 0, s$n / n_possible, NA_real_),
            posterior_simmap_median = s$median,
            posterior_simmap_mean = s$mean,
            posterior_simmap_q025 = s$q025,
            posterior_simmap_q975 = s$q975,
            posterior_simmap_hpd_low = s$hpd_low,
            posterior_simmap_hpd_high = s$hpd_high,
            posterior_simmap_hpd_width = s$hpd_width,
            posterior_simmap_p_after_1943 = s$p_after_1943,
            mean_jaccard = if (nrow(ass) > 0) mean(ass$jaccard, na.rm = TRUE) else NA_real_,
            min_jaccard = if (nrow(ass) > 0) min(ass$jaccard, na.rm = TRUE) else NA_real_,
            match_jaccard_threshold = MATCH_JACCARD_THRESHOLD,
            min_posterior_match_support = MIN_POSTERIOR_MATCH_SUPPORT,
            is_robust_origin = is.finite(n_possible) && n_possible > 0 && (s$n / n_possible) >= MIN_POSTERIOR_MATCH_SUPPORT,
            hpd_method = "BEAST_posterior_trees_plus_SIMMAP_transition_dates_matched_by_descendant_carrier_sets",
            stringsAsFactors = FALSE
        )
    }
}
origin_hpd <- if (length(origin_hpd_rows) > 0) do.call(rbind, origin_hpd_rows) else empty_df()
write.table(origin_hpd, file.path(OUTDIR, "pbp_origin_hpd.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Per-locus HPD for all internal posterior SIMMAP gains, useful when individual event matching is unstable.
locus_hpd_rows <- list()
if (nrow(posterior_gain_events) > 0) {
    for (locus in unique(posterior_gain_events$locus)) {
        pe <- posterior_gain_events[posterior_gain_events$locus == locus, , drop = FALSE]
        s <- summarise_dates(pe$event_date, threshold = PENICILLIN_YEAR)
        n_possible <- locus_map_counts$n_posterior_tree_maps_attempted[locus_map_counts$locus == locus]
        if (length(n_possible) == 0) n_possible <- NA_integer_
        locus_hpd_rows[[length(locus_hpd_rows) + 1L]] <- data.frame(
            locus = locus,
            label = locus_label(locus),
            posterior_n_gain_transitions = s$n,
            posterior_n_possible_tree_maps = n_possible,
            posterior_mean_gains_per_tree_map = ifelse(is.finite(n_possible) && n_possible > 0, s$n / n_possible, NA_real_),
            posterior_all_gain_median = s$median,
            posterior_all_gain_q025 = s$q025,
            posterior_all_gain_q975 = s$q975,
            posterior_all_gain_hpd_low = s$hpd_low,
            posterior_all_gain_hpd_high = s$hpd_high,
            posterior_all_gain_hpd_width = s$hpd_width,
            posterior_all_gain_p_after_1943 = s$p_after_1943,
            hpd_method = "all_internal_posterior_SIMMAP_0to1_gain_dates_by_locus",
            stringsAsFactors = FALSE
        )
    }
}
locus_hpd <- if (length(locus_hpd_rows) > 0) do.call(rbind, locus_hpd_rows) else empty_df()
write.table(locus_hpd, file.path(OUTDIR, "pbp_posterior_simmap_locus_hpd.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Join HPD to primary MCC origins.
primary_final <- if (nrow(mcc_primary_gains) > 0) {
    hpd_cols <- intersect(c("event_uid", "posterior_n_assigned", "posterior_n_possible_tree_maps",
                            "posterior_event_matching_support", "posterior_simmap_median",
                            "posterior_simmap_q025", "posterior_simmap_q975",
                            "posterior_simmap_hpd_low", "posterior_simmap_hpd_high",
                            "posterior_simmap_hpd_width", "posterior_simmap_p_after_1943",
                            "mean_jaccard", "min_jaccard", "hpd_method"), names(origin_hpd))
    safe_merge(mcc_primary_gains, origin_hpd[, hpd_cols, drop = FALSE], by = "event_uid", all.x = TRUE)
} else empty_df()
write.table(primary_final, file.path(OUTDIR, "pbp_dated_independent_origins.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
# Backward-compatible aliases.
write.table(primary_final, file.path(OUTDIR, "pbp_origin_events_dated.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(primary_final, file.path(OUTDIR, "pbp_gain_events_dated.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Locus summary joined with posterior per-locus gain HPD.
summary_final <- summary_df
if (nrow(locus_hpd) > 0) {
    summary_final <- safe_merge(summary_final, locus_hpd, by = c("locus", "label"), all.x = TRUE)
}
write.table(summary_final, file.path(OUTDIR, "pbp_temporal_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Input check.
input_check <- data.frame(
    item = c("mcc_tree", "posterior_trees", "matrix_file", "date_source", "output_dir",
             "posterior_burnin", "posterior_thin", "posterior_simmaps_per_tree",
             "max_posterior_trees", "match_jaccard_threshold", "root_pi", "simmap_model"),
    value = c(mcc_obj$file, ifelse(is.null(post_obj$trees), NA_character_, post_obj$file), MATRIX_FILE, date_obj$file, OUTDIR,
              POSTERIOR_BURNIN, POSTERIOR_THIN, POSTERIOR_SIMMAPS_PER_TREE,
              MAX_POSTERIOR_TREES, MATCH_JACCARD_THRESHOLD,
              paste(names(ROOT_PI), ROOT_PI, sep = "=", collapse = ";"), SIMMAP_MODEL),
    stringsAsFactors = FALSE
)
write.table(input_check, file.path(OUTDIR, "pbp_temporal_input_check.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# =========================================================
# 6. Plot
# =========================================================

plot_file <- file.path(OUTDIR, "pbp_dated_independent_origins.pdf")
if (requireNamespace("ggplot2", quietly = TRUE) && nrow(primary_final) > 0) {
    suppressPackageStartupMessages(library(ggplot2))
    plot_df <- primary_final
    plot_df$label <- factor(plot_df$label, levels = rev(unname(FOCAL_LABELS[FOCAL_LOCI])))
    plot_df$plot_x <- ifelse(is.finite(plot_df$posterior_simmap_median),
                             plot_df$posterior_simmap_median, plot_df$event_date)
    plot_df$plot_low <- ifelse(is.finite(plot_df$posterior_simmap_hpd_low),
                               plot_df$posterior_simmap_hpd_low, plot_df$origin_lower)
    plot_df$plot_high <- ifelse(is.finite(plot_df$posterior_simmap_hpd_high),
                                plot_df$posterior_simmap_hpd_high, plot_df$origin_upper)
    plot_df$uncertainty_type <- ifelse(is.finite(plot_df$posterior_simmap_hpd_low),
                                       "posterior_SIMMAP_95HPD", "MCC_branch_interval_only")
    p <- ggplot(plot_df, aes(x = plot_x, y = label)) +
        geom_vline(xintercept = PENICILLIN_YEAR, linetype = "dotted", color = "grey45") +
        annotate("text", x = PENICILLIN_YEAR, y = Inf, label = "Penicillin\nintroduced",
                 hjust = 1.1, vjust = 1.1, size = 2.8, color = "grey35") +
        geom_errorbar(aes(xmin = plot_low, xmax = plot_high, linetype = uncertainty_type),
                      orientation = "y", width = 0.12, linewidth = 0.42, alpha = 0.65, na.rm = TRUE) +
        geom_point(aes(size = n_carrier_descendants, shape = origin_type), alpha = 0.88, na.rm = TRUE) +
        scale_shape_manual(values = c("internal_clade" = 16, "terminal_singleton" = 1), name = "Origin type") +
        scale_linetype_manual(values = c("posterior_SIMMAP_95HPD" = "solid",
                                         "MCC_branch_interval_only" = "dashed"), name = "Uncertainty") +
        scale_size_area(max_size = 6, name = "Carrier descendants") +
        labs(title = "Dated independent PBP origins with posterior SIMMAP HPD intervals",
             subtitle = "Points are posterior SIMMAP medians where matched; otherwise MCC branch midpoints. Terminal singleton gains excluded by default.",
             x = "Calendar year", y = "Focal PBP variant") +
        theme_bw(base_size = 10) +
        theme(legend.position = "bottom", axis.text.y = element_text(face = "bold"))
    ggsave(plot_file, p, width = 11.5, height = 6.3)
    info("[Output] %s", plot_file)
}

# =========================================================
# 7. Console summary
# =========================================================

info("")
info("========================================")
info("Output summary")
info("========================================")
info("Primary MCC origins + HPD: %s", file.path(OUTDIR, "pbp_dated_independent_origins.tsv"))
info("Origin HPD table:       %s", file.path(OUTDIR, "pbp_origin_hpd.tsv"))
info("Posterior assignments:  %s", file.path(OUTDIR, "pbp_posterior_simmap_origin_assignments.tsv"))
info("Per-locus HPD:          %s", file.path(OUTDIR, "pbp_posterior_simmap_locus_hpd.tsv"))

if (nrow(origin_hpd) > 0) {
    show_cols <- intersect(c("label", "mcc_independent_origin_index", "mcc_event_date",
                             "posterior_n_assigned", "posterior_simmap_median",
                             "posterior_simmap_hpd_low", "posterior_simmap_hpd_high",
                             "posterior_simmap_p_after_1943", "posterior_event_matching_support"),
                           names(origin_hpd))
    print(origin_hpd[, show_cols, drop = FALSE], row.names = FALSE)
} else {
    info("No matched posterior HPD intervals were produced. Check posterior .trees path, SIMMAP settings, and Jaccard threshold.")
}

info("")
info("Recommended manuscript wording:")
info("  'Points denote MCC or posterior-SIMMAP median origin dates; horizontal bars denote 95%% HPD intervals from BEAST posterior trees and stochastic maps when event matching was possible, otherwise MCC branch spans are shown as conservative point-estimate intervals.'")
info("Origins with posterior_event_matching_support below %.2f are flagged is_robust_origin=FALSE in pbp_origin_hpd.tsv.", MIN_POSTERIOR_MATCH_SUPPORT)
