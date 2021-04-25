#' Plot Interval Probabilities
#'
#' Plot Interval Probabilities
#'
#' @export
#'
#'
predsu_plot_pintervals <- function(pintervals, lambda_surv = NULL) {
    x  <- pintervals[, "T1"]
    pl <- cumsum(pintervals[, "PL"])
    pu <- cumsum(pintervals[, "PU"])

    plot(x, pl, type = "l")
    lines(x, pu, col = "blue", type = "l")

    if (!is.null(lambda_surv)) {
        ptruth <- pexp(x, lambda_surv)
        lines(x, ptruth, col = "red", type = "l")
    }
}
