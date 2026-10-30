# Prespecified biological scope and descriptive count labels.
PBP_SCOPE_GENES <- c("mrcA", "pbp1", "pbp2", "tp47")
PBP_FOCAL_LOCI <- c("mrcA:p.Ala506Thr", "mrcA:p.Ala506Val",
                    "mrcA:p.Met625Val", "mrcA:p.Gly708Ser",
                    "pbp1:p.Pro564Leu", "tp47:p.Ser394Arg")
PBP_TP_CODES <- c(mrcA="TP0705", pbp1="TP0500", pbp2="TP0760", tp47="TP0574")

pbp_tip_count_category <- function(n1, n0) {
    minor <- min(n1, n0)
    if (minor == 0) "single_observed_state" else if (minor < 3)
        "minor_state_count_1to2" else if (minor < 10)
        "minor_state_count_3to9" else "minor_state_count_ge10"
}

pbp_count_interpretation <- function(n1, n0) {
    sprintf("Tip-count category only: carriers=%d; noncarriers=%d; minor state=%d. Assess model and mapping sensitivity separately.",
            n1, n0, min(n1, n0))
}

pbp_locus_seed <- function(locus, offset=0) {
    as.integer((42 + sum(as.double(utf8ToInt(locus)) * seq_along(utf8ToInt(locus))) + offset) %% 2147483646)
}

pbp_read_inputs <- function(tree_file, matrix_file) {
    tree <- ape::read.tree(tree_file)
    mat <- read.table(matrix_file, sep="\t", header=TRUE, row.names=1,
                      check.names=FALSE, quote="", comment.char="")
    if (anyDuplicated(tree$tip.label) || anyDuplicated(colnames(mat)))
        stop("Duplicate tree tips or matrix locus names")
    shared <- intersect(tree$tip.label, rownames(mat))
    if (length(shared) < 3) stop("Fewer than three shared samples")
    if (length(shared) < length(tree$tip.label))
        tree <- ape::drop.tip(tree, setdiff(tree$tip.label, shared))
    mat <- mat[tree$tip.label, , drop=FALSE]
    for (locus in names(mat)) {
        x <- mat[[locus]]
        y <- suppressWarnings(as.numeric(as.character(x)))
        if (any(!is.na(x) & is.na(y)) || any(!is.na(y) & !(y %in% c(0,1))))
            stop(sprintf("Invalid ternary matrix value at %s", locus))
        mat[[locus]] <- y
    }
    if (is.null(tree$edge.length) || any(!is.finite(tree$edge.length)) ||
        any(tree$edge.length < 0) || sum(tree$edge.length) <= 0)
        stop("Tree requires finite nonnegative branch lengths with positive total length")
    list(tree=tree, matrix=mat)
}

pbp_make_inventory <- function(tree, mat) {
    rows <- lapply(names(mat), function(locus) {
        x <- mat[[locus]]; names(x) <- tree$tip.label
        n1 <- sum(x==1,na.rm=TRUE); n0 <- sum(x==0,na.rm=TRUE); missing <- sum(is.na(x))
        gene <- sub(":.*$", "", locus); in_scope <- gene %in% PBP_SCOPE_GENES
        status <- if (!in_scope) "outside_PBP_scope" else if (n1+n0==0)
            "no_known_states" else if (n1==0) "no_observed_carriers" else if (n0==0)
            "all_known_samples_are_carriers" else if (n1<3)
            "low_carrier_count_1to2_descriptive" else "eligible_ER_ARD"
        score <- if (n1+n0==0) NA_real_ else 0
        mono <- NA
        if (n1>0 && n0>0) {
            tr <- if (missing) ape::drop.tip(tree, names(x)[is.na(x)]) else tree
            z <- x[tr$tip.label]
            dat <- phangorn::phyDat(matrix(as.character(z),ncol=1,
                       dimnames=list(names(z),NULL)),type="USER",levels=c("0","1"))
            score <- as.numeric(phangorn::parsimony(tr,dat))
            mono <- if (n1==1) TRUE else ape::is.monophyletic(tr,names(z)[z==1])
        }
        code <- unname(PBP_TP_CODES[gene]); if (!length(code)) code <- NA_character_
        data.frame(locus=locus,gene=gene,TP_code=code,is_focal=locus %in% PBP_FOCAL_LOCI,
                   in_PBP_scope=in_scope,n_mutant=n1,n_noncarrier=n0,n_wildtype=n0,
                   n_missing=missing,n_total=n1+n0+missing,n_known=n1+n0,
                   carrier_frequency=if (n1+n0>0) n1/(n1+n0) else NA_real_,
                   minor_state_count=min(n1,n0),tip_count_category=pbp_tip_count_category(n1,n0),
                   scope_status=status,eligible_ER_ARD=status=="eligible_ER_ARD",
                   parsimony_steps=score,carriers_monophyletic_among_known=mono,
                   carrier_tip_ids=paste(names(x)[!is.na(x) & x==1],collapse=";"),
                   stringsAsFactors=FALSE)
    })
    do.call(rbind,rows)
}

pbp_map_summary <- function(tree, states, Q, nsim, seed) {
    set.seed(seed)
    dimnames(Q) <- list(c("0","1"),c("0","1"))
    fac <- factor(as.character(states[tree$tip.label]),levels=c("0","1"))
    names(fac) <- tree$tip.label
    maps <- suppressMessages(phytools::make.simmap(tree,fac,Q=Q,nsim=nsim,
                              pi=c(0.5,0.5),message=FALSE))
    if (nsim==1) maps <- list(maps)
    counts <- matrix(0,nrow=nsim,ncol=2,dimnames=list(NULL,c("gain","loss")))
    n_tips <- length(tree$tip.label)
    node_mut <- matrix(NA_real_,nrow=nsim,ncol=tree$Nnode)
    for (i in seq_len(nsim)) {
        sm <- maps[[i]]; parents <- unique(sm$edge[,1])
        for (parent in parents) {
            edges <- which(sm$edge[,1]==parent)
            start <- vapply(edges,function(e) names(sm$maps[[e]])[1],character(1))
            if (length(unique(start))!=1) stop("Inconsistent mapped states at a shared parent")
            node_mut[i,parent-n_tips] <- as.integer(start[1]=="1")
        }
        for (history in sm$maps) {
            z <- names(history)
            if (length(z)>1) {
                counts[i,"gain"] <- counts[i,"gain"] + sum(head(z,-1)=="0" & tail(z,-1)=="1")
                counts[i,"loss"] <- counts[i,"loss"] + sum(head(z,-1)=="1" & tail(z,-1)=="0")
            }
        }
    }
    if (anyNA(node_mut)) stop("Missing mapped internal-node states")
    summarize <- function(x) c(mean=mean(x),sd=sd(x),q025=unname(quantile(x,.025)),
                               median=median(x),q975=unname(quantile(x,.975)))
    list(gain=summarize(counts[,"gain"]),loss=summarize(counts[,"loss"]),
         counts=counts,node_probability=colMeans(node_mut))
}
