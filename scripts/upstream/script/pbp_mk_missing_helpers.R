# Missing-data likelihoods for binary Mk models; unknown tips remain unknown.
pbp_tip_likelihoods <- function(states01, tip_labels) {
    if (length(states01) != length(tip_labels)) stop("Tip-state length mismatch")
    if (!is.null(names(states01)) && !identical(names(states01), tip_labels)) {
        stop("Tip-state order does not match tree$tip.label")
    }
    x <- suppressWarnings(as.numeric(as.character(states01)))
    if (any(!is.na(states01) & is.na(x)) ||
        any(!is.na(x) & (!is.finite(x) | !(x %in% c(0, 1))))) {
        stop("Binary states must be exactly 0, 1, or NA")
    }
    L <- matrix(1, nrow = length(x), ncol = 2,
                dimnames = list(tip_labels, NULL))
    known <- which(!is.na(x))
    if (length(known)) {
        L[known, ] <- 0
        L[cbind(known, as.integer(x[known]) + 1L)] <- 1
    }
    L
}

pbp_mk_fit_ok <- function(fit) {
    if (is.null(fit) || !isTRUE(fit$success)) return(FALSE)
    if (!is.null(fit$converged) && !isTRUE(fit$converged)) return(FALSE)
    Q <- fit$transition_matrix
    is.matrix(Q) && identical(dim(Q), c(2L, 2L)) &&
        all(is.finite(Q)) && all(Q[row(Q) != col(Q)] >= 0) &&
        max(abs(rowSums(Q))) < 1e-7 && is.finite(fit$loglikelihood)
}

pbp_fit_mk <- function(tree, states01, rate_model = "ER", root_prior = "auto") {
    # Validate order with names, then pass unnamed likelihood rows in that order.
    # castor 1.8.5 checks matrix row names with vector-valued && (missing any()).
    # Removing labels after our own validation avoids that defect while keeping
    # castor's remaining input checks enabled.
    L <- unname(pbp_tip_likelihoods(states01, tree$tip.label))
    errors <- character()
    attempts <- list()
    for (algorithm in c("nlminb", "optim")) {
        fit <- tryCatch(
            castor::fit_mk(tree, Nstates = 2, tip_states = NULL, tip_priors = L,
                           rate_model = rate_model, root_prior = root_prior, Ntrials = 3,
                           optim_algorithm = algorithm,
                           optim_max_iterations = 3000),
            error = function(e) list(success = FALSE, error = conditionMessage(e))
        )
        attempts[[algorithm]] <- fit
        if (pbp_mk_fit_ok(fit)) { fit$optimizer_used <- algorithm; return(fit) }
        msg <- if (!is.null(fit$error)) fit$error else "Non-converged or non-finite Mk fit"
        errors <- c(errors, paste(algorithm, msg, sep = ": "))
    }
    if (identical(rate_model,"ER")) {
        recovered <- pbp_retry_ER(tree, L, root_prior, attempts)
        if (pbp_mk_fit_ok(recovered)) return(recovered)
        errors <- c(errors,recovered$error)
    }
    fit$success <- FALSE
    fit$error <- paste(errors, collapse = "; ")
    fit
}

pbp_asr_mk <- function(tree, states01, transition_matrix, reroot = TRUE, root_prior = "auto") {
    result <- castor::asr_mk_model(
        tree = tree, tip_states = NULL, Nstates = 2,
        tip_priors = unname(pbp_tip_likelihoods(states01, tree$tip.label)),
        transition_matrix = transition_matrix, reroot = reroot, root_prior = root_prior,
        include_ancestral_likelihoods = TRUE
    )
    if (!isTRUE(result$success) ||
        !identical(dim(result$ancestral_likelihoods), c(as.integer(tree$Nnode), 2L)) ||
        any(!is.finite(result$ancestral_likelihoods))) {
        msg <- if (!is.null(result$error)) result$error else "Invalid ancestral likelihoods"
        stop(msg)
    }
    result
}

# Retry a one-rate ER fit only after both normal optimizers failed validation.
# Keep model, root prior and input likelihoods unchanged. Never force convergence.
pbp_retry_ER <- function(tree, L, root_prior, initial_attempts) {
    structurally_valid <- function(f) {
        if (is.null(f) || !isTRUE(f$success)) return(FALSE)
        Q <- f$transition_matrix
        is.matrix(Q) && identical(dim(Q),c(2L,2L)) && all(is.finite(Q)) &&
            all(Q[row(Q)!=col(Q)]>=0) && max(abs(rowSums(Q)))<1e-7 &&
            is.finite(f$loglikelihood)
    }
    valid <- Filter(structurally_valid, initial_attempts)
    if (!length(valid)) return(list(success=FALSE,error="ER retry has no finite initial fit for a scale/reference"))
    best_initial <- valid[[which.max(vapply(valid,function(f) f$loglikelihood,numeric(1)))]]
    base_q <- best_initial$transition_matrix[1,2]
    if (!is.finite(base_q) || base_q<=0)
        return(list(success=FALSE,error="ER retry requires a positive rate scale"))
    all_fits <- initial_attempts
    labels <- names(initial_attempts)
    for (factor in c(0.01,0.1,0.5,2,10,100)) {
        q <- base_q*factor
        Q <- matrix(c(-q,q,q,-q),2,2)
        for (algorithm in c("nlminb","optim")) {
            label <- paste0(algorithm,"_single_start_x",factor)
            candidate <- tryCatch(castor::fit_mk(tree,Nstates=2,tip_states=NULL,
                tip_priors=L,rate_model="ER",root_prior=root_prior,Ntrials=1,
                guess_transition_matrix=Q,optim_algorithm=algorithm,optim_max_iterations=3000),
                error=function(e) list(success=FALSE,error=conditionMessage(e)))
            all_fits[[label]] <- candidate
            labels <- c(labels,label)
        }
    }
    audit <- do.call(rbind,lapply(names(all_fits),function(label) {
        f <- all_fits[[label]]
        data.frame(attempt=label,success=isTRUE(f$success),converged=isTRUE(f$converged),
            loglikelihood=if(is.null(f$loglikelihood)) NA_real_ else f$loglikelihood,
            q=if(is.null(f$transition_matrix)) NA_real_ else f$transition_matrix[1,2],
            accepted=pbp_mk_fit_ok(f),stringsAsFactors=FALSE)
    }))
    cat("[ER RETRY] root_prior=",paste(root_prior,collapse=","),"\n",sep="")
    print(audit,row.names=FALSE)
    passing <- Filter(pbp_mk_fit_ok,all_fits)
    finite_fits <- Filter(structurally_valid,all_fits)
    if (!length(passing)) return(list(success=FALSE,error="No ER retry reported valid convergence",retry_audit=audit))
    best <- passing[[which.max(vapply(passing,function(f) f$loglikelihood,numeric(1)))]]
    highest <- max(vapply(finite_fits,function(f) f$loglikelihood,numeric(1)))
    gap <- highest-best$loglikelihood
    if (!is.finite(gap) || gap>1e-6)
        return(list(success=FALSE,error=paste("ER converged retry falls below best finite likelihood by",gap),retry_audit=audit))
    fixed <- tryCatch(castor::asr_mk_model(tree,Nstates=2,tip_states=NULL,tip_priors=L,
        transition_matrix=best$transition_matrix,root_prior=root_prior,reroot=FALSE,
        include_ancestral_likelihoods=FALSE),error=function(e) NULL)
    if (is.null(fixed) || !isTRUE(fixed$success) || !is.finite(fixed$loglikelihood) ||
        abs(fixed$loglikelihood-best$loglikelihood)>1e-5)
        return(list(success=FALSE,error="Fixed-Q likelihood did not validate the converged ER retry",retry_audit=audit))
    best$optimizer_used <- "validated_ER_single_start_retries"
    best$retry_audit <- audit
    best$retry_likelihood_gap <- gap
    cat(sprintf("[ER RETRY PASS] logL=%.12f; q=%.12g; best-finite gap=%.3g\n",
        best$loglikelihood,best$transition_matrix[1,2],gap))
    best
}
