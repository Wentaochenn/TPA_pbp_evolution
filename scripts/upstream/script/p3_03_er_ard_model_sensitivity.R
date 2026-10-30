#!/usr/bin/env Rscript
# ER/ARD sensitivity for every eligible variant in the prespecified PBP genes.
suppressPackageStartupMessages({library(ape);library(castor);library(phytools);library(phangorn)})
source("script/pbp_mk_missing_helpers.R")
source("script/pbp_variant_scope_helpers.R")
TREE_FILE <- "results/03_iqtree/tpa_ss14_rooted.treefile"
MATRIX_FILE <- Sys.getenv("PBP_MATRIX_FILE", "results/05_extracted_variants/callability/sample_locus_ternary_matrix.tsv")
OUTDIR <- "results/08_model_sensitivity"
NSIM <- suppressWarnings(as.integer(Sys.getenv("PBP_SIMMAP_N", "100")))
if (length(NSIM)!=1 || is.na(NSIM) || NSIM<2) stop("PBP_SIMMAP_N must be an integer >=2")
dir.create(OUTDIR,recursive=TRUE,showWarnings=FALSE)
cat("R:",R.version.string,"castor:",as.character(packageVersion("castor")),
    "phytools:",as.character(packageVersion("phytools")),"\n")
inputs <- pbp_read_inputs(TREE_FILE,MATRIX_FILE)
tree <- inputs$tree; mat <- inputs$matrix
cat("Shared tree tips:",length(tree$tip.label),"matrix loci:",ncol(mat),"\n")
inventory <- pbp_make_inventory(tree,mat)
write_tsv <- function(df,name) write.table(df,file.path(OUTDIR,name),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
write_tsv(inventory,"pbp_variant_inventory.tsv")
write_tsv(inventory[inventory$in_PBP_scope & inventory$n_mutant %in% 1:2,,drop=FALSE],
          "pbp_low_carrier_count_1to2.tsv")
targets <- inventory$locus[inventory$eligible_ER_ARD]
cat("Automatically eligible PBP loci:",length(targets),"; maps/model/locus:",NSIM,"\n")
if (!length(targets)) stop("No eligible PBP loci; inventory was written")

signature_files <- c(TREE_FILE,MATRIX_FILE,"script/pbp_mk_missing_helpers.R",
                     "script/pbp_variant_scope_helpers.R","script/p3_03_er_ard_model_sensitivity.R",
                     "script/pbp_binary_fit.py")
signature <- c(paste(names(tools::md5sum(signature_files)),tools::md5sum(signature_files),sep="="),
               paste0("maps_per_model=",NSIM),"root_prior=flat; mapping_pi=0.5,0.5",
               "scope=mrcA,pbp1,pbp2,tp47; carrier_min=3; noncarrier_min=1",
               R.version.string,paste0("castor=",packageVersion("castor")),
               paste0("phytools=",packageVersion("phytools")),paste0("ape=",packageVersion("ape")))
signature_tmp <- tempfile();writeLines(signature,signature_tmp)
signature_id <- unname(tools::md5sum(signature_tmp));unlink(signature_tmp)
writeLines(signature,file.path(OUTDIR,"expanded_input_signature.txt"))
cache_dir <- file.path(OUTDIR,"expanded_checkpoints",signature_id)
dir.create(cache_dir,recursive=TRUE,showWarnings=FALSE)
tree_dir <- file.path(OUTDIR,"node_comparison_trees")
dir.create(tree_dir,recursive=TRUE,showWarnings=FALSE)

get_root_probability <- function(tr,x,fit) {
    if (!pbp_mk_fit_ok(fit)) return(NA_real_)
    ans <- pbp_asr_mk(tr,x,fit$transition_matrix,reroot=FALSE,root_prior="flat")
    root <- setdiff(unique(tr$edge[,1]),tr$edge[,2])
    if (length(root)!=1) stop("Expected one root")
    ans$ancestral_likelihoods[root-length(tr$tip.label),2]
}
state_class <- function(p) ifelse(is.na(p),NA_character_,
                    ifelse(p<=.2,"noncarrier",ifelse(p>=.8,"carrier","uncertain")))
fit_error <- function(f) if (pbp_mk_fit_ok(f)) "" else
                    if (!is.null(f$error)) f$error else "invalid or non-converged fit"
boundary_flag <- function(f) {
    if (!pbp_mk_fit_ok(f)) return(NA)
    q <- f$transition_matrix[cbind(c(1,2),c(2,1))]
    if (max(q)==0) TRUE else min(q)/max(q)<1e-6
}

# Numerically fit the same ER/ARD likelihood if either primary optimizer fails.
fit_binary_pair <- function(tr,x,locus) {
    folder <- file.path(OUTDIR,"numerical_fits",gsub("[^A-Za-z0-9]","_",locus))
    dir.create(folder,recursive=TRUE,showWarnings=FALSE)
    tf <- file.path(folder,"known_tips.tree")
    sf <- file.path(folder,"known_states.tsv")
    result <- file.path(folder,"rates.tsv")
    ape::write.tree(tr,tf)
    write.table(data.frame(sample=tr$tip.label,state=as.integer(x[tr$tip.label])),
                sf,sep="\t",quote=FALSE,row.names=FALSE)
    python <- Sys.getenv("PBP_ANALYSIS_PY","python3")
    code <- system2(python,c(shQuote("script/pbp_binary_fit.py"),
                   "--tree",shQuote(tf),"--states",shQuote(sf),"--out",shQuote(result)),
                   stdout=file.path(folder,"fit.log"),stderr=file.path(folder,"fit.log"))
    if (code!=0L || !file.exists(result)) stop("Numerical ER/ARD fitting failed: ",folder)
    rates <- read.delim(result,stringsAsFactors=FALSE)
    if (nrow(rates)!=2L || !setequal(rates$model,c("ER","ARD"))) stop("Invalid rate table")
    ans <- setNames(lapply(c("ER","ARD"),function(model) {
        r <- rates[rates$model==model,,drop=FALSE]
        Q <- matrix(c(-r$q01,r$q01,r$q10,-r$q10),2,2,byrow=TRUE)
        f <- list(success=TRUE,converged=TRUE,transition_matrix=Q,
                  loglikelihood=r$loglikelihood,AIC=r$AIC,
                  optimizer_used="scaled_binary_likelihood_multiple_starts")
        if (!pbp_mk_fit_ok(f) || abs(r$independent_loglikelihood-r$loglikelihood)>1e-8)
            stop("Invalid fitted likelihood")
        f
    }),c("ER","ARD"))
    ans
}

analyze_locus <- function(locus,tr,x,inv) {
    cat("[LOCUS]",locus,"carriers=",inv$n_mutant,"noncarriers=",inv$n_noncarrier,"\n")
    fits <- list();maps <- list();errors <- character();root <- c(ER=NA_real_,ARD=NA_real_)
    for (model in c("ER","ARD")) {
        set.seed(pbp_locus_seed(locus,if(model=="ER") 0 else 1000))
        fits[[model]] <- pbp_fit_mk(tr,x,model,root_prior="flat")
    }
    if (!all(vapply(fits,pbp_mk_fit_ok,logical(1)))) {
        recovered <- tryCatch(fit_binary_pair(tr,x,locus),error=function(e) {
            errors <<- c(errors,paste("numerical likelihood fit:",conditionMessage(e)))
            NULL
        })
        if (!is.null(recovered)) fits <- recovered
    }
    for (model in c("ER","ARD")) {
        fit <- fits[[model]]
        if (!pbp_mk_fit_ok(fit)) {errors <- c(errors,paste(model,fit_error(fit)));next}
        root[model] <- tryCatch(get_root_probability(tr,x,fit),error=function(e) {
            errors <<- c(errors,paste(model,"root:",conditionMessage(e)));NA_real_
        })
        maps[[model]] <- tryCatch(
            pbp_map_summary(tr,x,fit$transition_matrix,NSIM,
                           pbp_locus_seed(locus,if(model=="ER") 2000 else 3000)),
            error=function(e) {errors <<- c(errors,paste(model,"mapping:",conditionMessage(e)));NULL})
        cat("  ",model,"fit_ok=",pbp_mk_fit_ok(fit),"mapping_ok=",!is.null(maps[[model]]),"\n")
    }
    ok <- vapply(fits,pbp_mk_fit_ok,logical(1))
    aic <- vapply(fits,function(f) if(pbp_mk_fit_ok(f)) f$AIC else NA_real_,numeric(1))
    q <- function(model,i,j) if(ok[[model]]) fits[[model]]$transition_matrix[i,j] else NA_real_
    delta <- aic[["ARD"]]-aic[["ER"]]
    root_class <- state_class(root)
    root_agree <- if (anyNA(root_class) || any(root_class=="uncertain")) NA else root_class[1]==root_class[2]
    status <- if (!all(ok)) "fit_failed" else if(anyNA(root)) "root_reconstruction_failed" else
                  if (is.null(maps$ER) || is.null(maps$ARD)) "mapping_failed" else "complete"
    row <- data.frame(locus=locus,gene=inv$gene,TP_code=inv$TP_code,is_focal=inv$is_focal,
        n_mutant=inv$n_mutant,n_wildtype=inv$n_noncarrier,n_missing=inv$n_missing,
        mut_freq=inv$carrier_frequency,minor_state_count=inv$minor_state_count,
        tip_count_category=inv$tip_count_category,parsimony_changes=inv$parsimony_steps,
        tree_scope="known_state_tips_for_this_variant",root_prior="flat",
        ER_q=q("ER",1,2),ER_AIC=aic[["ER"]],ER_loglikelihood=if(ok[["ER"]]) fits$ER$loglikelihood else NA_real_,
        ER_root_Pmut=unname(root[["ER"]]),ARD_q_gain=q("ARD",1,2),ARD_q_loss=q("ARD",2,1),
        ARD_ratio=if(is.finite(q("ARD",2,1)) && q("ARD",2,1)>0) q("ARD",1,2)/q("ARD",2,1) else NA_real_,
        ARD_AIC=aic[["ARD"]],ARD_loglikelihood=if(ok[["ARD"]]) fits$ARD$loglikelihood else NA_real_,
        ARD_root_Pmut=unname(root[["ARD"]]),delta_AIC_ARD_ER=delta,
        preferred_model=if(is.finite(delta)) {if(delta< -2) "ARD" else "ER"} else "unresolved",
        AIC_evidence=if(!is.finite(delta)) "unresolved" else if(abs(delta)<=2) "weak_model_separation" else
                         if(delta<0) "ARD_lower_AIC" else "ER_lower_AIC",
        ER_fit_ok=unname(ok[["ER"]]),ARD_fit_ok=unname(ok[["ARD"]]),
        ER_fit_error=fit_error(fits$ER),ARD_fit_error=fit_error(fits$ARD),
        ARD_near_rate_boundary=boundary_flag(fits$ARD),
        ER_root_class=root_class[1],ARD_root_class=root_class[2],root_state_agree=unname(root_agree),
        mapping_n=NSIM,node_comparison_nnodes=tr$Nnode,analysis_status=status,analysis_error=paste(errors,collapse="; "),stringsAsFactors=FALSE)
    for (model in c("ER","ARD")) for (direction in c("gain","loss")) {
        for (stat in c("mean","sd","q025","median","q975"))
            row[[paste(model,"simmap",direction,stat,sep="_")]] <-
                if(is.null(maps[[model]])) NA_real_ else unname(maps[[model]][[direction]][stat])
    }
    nodes <- data.frame();counts <- data.frame()
    for (model in c("ER","ARD")) if (!is.null(maps[[model]])) {
        cc <- data.frame(locus=locus,model=model,map_index=seq_len(NSIM),
                         gain=maps[[model]]$counts[,"gain"],loss=maps[[model]]$counts[,"loss"])
        counts <- rbind(counts,cc)
    }
    if (!is.null(maps$ER) && !is.null(maps$ARD)) {
        nt <- length(tr$tip.label);root_id <- setdiff(unique(tr$edge[,1]),tr$edge[,2])
        pER <- maps$ER$node_probability;pARD <- maps$ARD$node_probability
        nodes <- data.frame(locus=locus,node=nt+seq_len(tr$Nnode),is_root=(nt+seq_len(tr$Nnode))==root_id,
                  n_known_tips=nt,ER_Pcarrier=pER,ARD_Pcarrier=pARD,
                  absolute_probability_difference=abs(pER-pARD),
                  ER_state_class=state_class(pER),ARD_state_class=state_class(pARD),
                  maps_per_model=NSIM)
        row$max_node_probability_difference <- max(nodes$absolute_probability_difference)
        row$n_nodes_with_changed_state_class <- sum(nodes$ER_state_class!=nodes$ARD_state_class)
    } else {
        row$max_node_probability_difference <- NA_real_
        row$n_nodes_with_changed_state_class <- NA_integer_
    }
    list(row=row,nodes=nodes,counts=counts)
}

results <- list()
for (locus in targets) {
    inv <- inventory[inventory$locus==locus,,drop=FALSE]
    x <- mat[[locus]];names(x) <- tree$tip.label
    tr <- if(anyNA(x)) ape::drop.tip(tree,names(x)[is.na(x)]) else tree
    x <- x[tr$tip.label]
    safe <- gsub("[^A-Za-z0-9]","_",locus)
    ape::write.tree(tr,file.path(tree_dir,paste0(safe,".tree")))
    cache <- file.path(cache_dir,paste0(safe,".rds"))
    if (file.exists(cache)) {
        cat("[CACHE]",locus,"\n");answer <- readRDS(cache)
    } else {
        answer <- analyze_locus(locus,tr,x,inv)
        if (answer$row$analysis_status=="complete") {
            tmp <- paste0(cache,".tmp");saveRDS(answer,tmp);if(!file.rename(tmp,cache)) stop("Checkpoint rename failed")
        }
    }
    results[[locus]] <- answer
}
combined <- do.call(rbind,lapply(results,`[[`,"row"));rownames(combined) <- NULL
write_tsv(combined,"er_ard_model_sensitivity.tsv")
bind_nonempty <- function(key) {
    rows <- lapply(results,`[[`,key);rows <- rows[vapply(rows,nrow,integer(1))>0]
    if(length(rows)) do.call(rbind,rows) else data.frame(locus=character())
}
write_tsv(bind_nonempty("nodes"),"er_ard_node_state_comparison.tsv")
write_tsv(bind_nonempty("counts"),"er_ard_simmap_counts.tsv")
audit <- inventory
m <- match(audit$locus,combined$locus)
audit$analysis_status <- ifelse(is.na(m),audit$scope_status,combined$analysis_status[m])
audit$analysis_error <- ifelse(is.na(m),"",combined$analysis_error[m])
write_tsv(audit,"pbp_variant_analysis_audit.tsv")
write_tsv(combined[combined$analysis_status!="complete",,drop=FALSE],"er_ard_unresolved_loci.tsv")
writeLines(c("Binary 0/1 denote absence/presence of this particular variant; NA is unobserved.",
 "All eligible loci in mrcA,pbp1,pbp2,tp47 are attempted; carriers>=3 and noncarriers>=1.",
 "Low-count loci remain in inventory; count categories do not establish reliability.",
 "Each model is fitted and mapped on the same tree of known-state tips for that variant.",
 "Both models use a flat root prior. The original ASR numeric results are retained separately.",
 "Mappings condition on fitted Q; empirical map quantiles do not integrate Q-estimation uncertainty.",
 "Node probabilities are map proportions on each saved known-tip tree; they have Monte Carlo uncertainty.",
 "root_state_agree only compares root states when both probabilities are <=0.2 or >=0.8.",
 "No binary label of overall evolutionary-conclusion robustness is assigned.",
 "AIC differences within +/-2 are recorded as weak separation; ER is the simpler reporting default.",
 "A relative ARD rate ratio below 1e-6 is a numerical boundary flag, not proof of biological irreversibility."),
 file.path(OUTDIR,"expanded_methods.txt"))
if (any(combined$analysis_status!="complete"))
    stop("ER/ARD analysis is incomplete; inspect er_ard_unresolved_loci.tsv")

good <- combined[combined$analysis_status=="complete",,drop=FALSE]
pdf(file.path(OUTDIR,"er_ard_gain_loss_comparison.pdf"),width=8,height=7)
for (direction in c("gain","loss")) {
    x <- good[[paste0("ER_simmap_",direction,"_mean")]]
    y <- good[[paste0("ARD_simmap_",direction,"_mean")]]
    lim <- range(c(x,y),finite=TRUE);if(diff(lim)==0)lim<-lim+c(-.5,.5)
    plot(x,y,xlim=lim,ylim=lim,pch=ifelse(good$is_focal,19,1),
         xlab=paste("ER: mean",direction,"count"),ylab=paste("ARD: mean",direction,"count"),
         main=paste("Model sensitivity:",direction,"counts"))
    abline(0,1,col="grey60",lty=2)
    text(x[good$is_focal],y[good$is_focal],labels=sub(":p[.]",":",good$locus[good$is_focal]),pos=3,cex=.6)
}
dev.off()
cat("[DONE] eligible=",nrow(combined)," complete=",sum(combined$analysis_status=="complete"),
    " unresolved=",sum(combined$analysis_status!="complete"),"\n")
