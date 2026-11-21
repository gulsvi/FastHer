###################################################################
###        Regularization of a reference LD matrix              ###
###################################################################
 
# Validates and completes the control list of the regularization
.reg_control <- function(reg_control) {
    ctl <- list(min_fraction = 0.95, min_rank = 5L, wald = 1.96)
    if (!is.null(reg_control)) {
        if (!is.list(reg_control) || !all(names(reg_control) %in% names(ctl))) {
            stop("'reg_control' must be a list with elements among: ",
                 paste(names(ctl), collapse = ", "), ".", call. = FALSE)
        }
        ctl[names(reg_control)] <- reg_control
    }
    ok <- is.numeric(ctl$min_fraction) && length(ctl$min_fraction) == 1L &&
        ctl$min_fraction >= 0 && ctl$min_fraction < 1 &&
        is.numeric(ctl$min_rank) && length(ctl$min_rank) == 1L && ctl$min_rank >= 1 &&
        is.numeric(ctl$wald) && length(ctl$wald) == 1L && ctl$wald > 0
    if (!isTRUE(ok)) {
        stop("Invalid 'reg_control': need min_fraction in [0, 1), min_rank >= 1, wald > 0.",
             call. = FALSE)
    }
    ctl$min_rank <- as.integer(ctl$min_rank)
    ctl
}
 
# Regularization of the reference LD by truncation of its spectrum.
# Vz2, Uval: squared projections and eigenvalues of the first components
# (descending order); cum: cumulative fraction of the LD variance retained by
# the first k components (k = 1, 2, ...). Starting from reg_tau, components are
# dropped one at a time from the tail while the Wald statistic for d = 1
# exceeds ctl$wald; dropping stops when the statistic is not larger than
# ctl$wald, when the retained fraction is not larger than ctl$min_fraction, or
# when the rank would become not larger than ctl$min_rank.
.reml_regularized <- function(Vz2, Uval, cum, n, m, reml, control_list, reg_tau, ctl) {
    rank <- min(length(Uval), max(sum(cum < reg_tau), ctl$min_rank + 1L))
    path <- NULL
    repeat {
        fit <- .reml_fit(Vz2[seq_len(rank)], Uval[seq_len(rank)], n, m, rank,
                         reml = reml, control_list = control_list)
        d <- fit$estimate[["d"]]
        wald <- abs(d - 1) / fit$se[["d"]]
        frac <- cum[rank]
        path <- rbind(path, data.frame(rank = rank, fraction = frac,
                                       h2 = fit$estimate[["h2"]], d = d,
                                       se_d = fit$se[["d"]], wald = wald))
        converged <- is.finite(wald) && wald <= ctl$wald
        if (converged || rank - 1L <= ctl$min_rank || frac <= ctl$min_fraction) break
        rank <- rank - 1L
    }
    list(fit = fit, rank = rank, fraction = frac, converged = converged, path = path)
}
