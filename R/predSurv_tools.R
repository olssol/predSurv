#' Test R function
#'
#' @export
#'
#'
predsu_test <- function() {
    print("Hello World!");
}

#' Simulate Survical Data
#'
#' Simulate a survival dataset
#'
#' @param n Sample size
#' @param t_ana Time from study starts to the (interim) anlysis
#' @param t_enrollment Time from study starts to the enrollment finishes
#'
#' @return A dataframe with enrollment time, time of survival, time of
#'     censoring, time of event (survival or censoring), and the indicator of
#'     censoring
#'
#' @export
#'
#'
predsu_simu <- function(n = 300,
                        t_ana         = 1,
                        t_enrollment  = 0.25,
                        lambda_surv   = 1,
                        lambda_censor = 0.5) {

    time_enroll <- runif(n, 0, min(t_ana, t_enrollment))
    y           <- rexp(n, lambda_surv)
    x           <- rexp(n, lambda_censor)
    censor      <- as.integer(x < y)

    t_event     <- censor * x + (1 - censor) * y
    t_tmp       <- t_event + time_enroll
    inx         <- which(t_tmp > t_ana)

    if (length(inx) > 0) {
        t_event[inx] <- t_ana - time_enroll[inx]
        censor[inx]  <- 2
    }

    data.frame(Time_Enroll = time_enroll,
               T_Survival  = y,
               T_Censor    = x,
               T_Event     = t_event,
               Ind_Censor  = factor(censor,
                                    levels = 0:2,
                                    label  = c("Dead",
                                               "Lost to follow up",
                                               "Censored at analysis")))
}

#' Partial censoring information
#'
#'
#' @export
#'
#'
predsu_pci <- function(vec_tcensor, vec_tsurv, option = c("upper", "lower")) {
    f_upper <- function(tc) {
        inx <- which(vec_tsurv < tc)
        vec_tsurv[max(inx)]
    }

    f_lower <- function(tc) {
        inx <- which(vec_tsurv > tc)
        vec_tsurv[min(inx)]
    }

    option <- match.arg(option)
    f      <- swicth(option,
                     upper = f_upper,
                     lower = f_lower)

    vec_tsurv <- c(0, sort(vec_tsurv), Inf)

    sapply(vec_tcensor, f)
}

#' Get probability for each interval
#'
#'
#' @export
#'
#'
predsu_pinterval <- function(vec_tsurv, vec_tcensor, t_max = 10) {

    n_surv <- length(vec_tsurv)
    n_tot  <- n_surv + length(vec_tcensor)
    tsurv  <- c(0, unname(sort(vec_tsurv)), t_max)

    rst <- NULL
    for (i in 0:n_surv) {
        t0_i <- tsurv[i + 1]
        t1_i <- tsurv[i + 2]

        ## C(i)
        c_i  <- length(which(vec_tcensor <= t1_i))

        if (0 == i) {
            c_i_1 <- 0
        } else {
            c_i_1 <- rst[i, "C"]
        }

        ## lambda(i)
        lu_i <- 1 / (n_tot - (i - 1) - c_i)
        ll_i <- 1 / (n_tot - (i - 1) - c_i_1)

        ## P(i)
        if (0 == i) {
            pu_i <- lu_i
            pl_i <- ll_i
        } else {
            pu_i <- prod(1 - rst[, "LU"]) * lu_i
            pl_i <- prod(1 - rst[, "LL"]) * ll_i
        }

        rst <- rbind(rst,
                     c(i  = i,
                       T0 = t0_i,
                       T1 = t1_i,
                       C  = c_i,
                       LU = lu_i,
                       LL = ll_i,
                       PU = pu_i,
                       PL = pl_i))
    }

    ## fix P(i) based on lambda_L
    rst[n_surv + 1, "PL"] <- 1 - sum(rst[1:n_surv, "PL"])

    ## return
    rst
}


#' Get probability of survival
#'
#' Get probability T > t
#'
#'
#' @export
#'
#'
predsu_psurv <- function(time, prob_ints, option = c("upper", "lower")) {
    option <- match.arg(option)

    ## survival probability
    prob <- switch(option,
                   upper = prob_ints[, "PU"],
                   lower = prob_ints[, "PL"])

    nint <- nrow(prob_ints)

    ## interval
    inx  <- max(which(prob_ints[, "T0"] < time))
    pinx <- prob_ints[inx, ]

    ## within interval probability > time
    p_int <- pinx["T1"] - time
    p_int <- p_int / (pinx["T1"] - pinx["T0"])
    p_int <- p_int * prob[inx]

    ## other interval probabilities
    if (inx == nint) {
        p_other <- 0
    } else {
        p_other <- sum(prob[(inx + 1):nint])
    }

    p_int + p_other
}

#' Predict survival time
#'
#' Predict survival time conditioning on censoring time based on interval
#' probabilities
#'
#'
#' @export
#'
#'
predsu_pred <- function(prob_ints, condition_t = 0,
                        option = c("upper", "lower")) {

    ## survival probability
    option <- match.arg(option)
    prob   <- switch(option,
                     upper = prob_ints[, "PU"],
                     lower = prob_ints[, "PL"])

    prob_ints <- cbind(prob_ints, P = prob)

    ## interval
    inx  <- min(which(prob_ints[, "T1"] > condition_t))
    if (inx > 1) {
        prob_ints <- prob_ints[-(1:(inx - 1)), , drop = FALSE]
    }

    ## within interval probability > time
    p_int <- prob_ints[1, "T1"] - condition_t
    p_int <- p_int / (prob_ints[1, "T1"] - prob_ints[1, "T0"])
    p_int <- p_int * prob_ints[1, "P"]

    prob_ints[1, "T0"] <- condition_t
    prob_ints[1, "P"]  <- p_int

    ## standardize
    prob <- prob_ints[, "P"] / sum(prob_ints[, "P"])
    prob <- cumsum(prob)

    ## random
    rand_p <- runif(1)
    inx    <- min(which(prob > rand_p))
    t0     <- prob_ints[inx, "T0"]
    t1     <- prob_ints[inx, "T1"]

    if (1 == inx) {
        p_prev <- 0
    } else {
        p_prev <- prob[inx - 1]
    }

    ## predicted time
    rst <- (rand_p - p_prev) / (prob[inx] - p_prev)
    rst <- rst * (t1 - t0)
    rst <- rst + t0

    attr(rst, "prob") <- 1 - rand_p
    rst
}



#' Predict All
#'
#' Predict survival time conditioning on censoring time based on interval
#' probabilities
#'
#' @param dta Survival dataset
#' @param v_enroll Column name for enrollment time (starting from 0 when the
#'     study started)
#' @param v_event Column name for event or censoring time
#' @param v_censor Column name for censoring indicator
#' @param censor_event Level of censor indicator corresponding to event
#' @param censor_lfu Level of censor indicator corresponding to censored due to
#'     lost to follow up
#' @param censor_ana Level of censor indicator corresponding to censored due
#'     to analysis time occur first
#' @param t_ana Time of the current (interim) analysis (starting from 0 when the study
#'     stared)
#' @param t_enrollment Duration of enrollment (starting from 0 when the study
#'     stared). Bigger than t_ana if need to enroll more patients
#' @param t_next_ana Time for the next analysis (starting from 0 when the
#'     study stared)
#' @param n_to_enroll Number of patients to be enrolled from the time of the
#'     interim analysis to the time enrollment will finish
#'
#' @return A data frame with the following new columns. Pred_T_Event: predicted
#'     event or censoring time; Pred_Ind_Censor: predicted censoring status;
#'     Pred_T_Surv: predicted survival time
#'
#' @export
#'
#'
predsu_pred_all <- function(dta,
                            censor_event    = "Dead",
                            censor_lfu      = "Lost to follow up",
                            censor_ana      = "Censored at analysis",
                            v_enroll        = "Time_Enroll",
                            v_event         = "T_Event",
                            v_censor        = "Ind_Censor",
                            t_ana           = 0.75,
                            t_enrollment    = 1,
                            t_next_ana      = 2,
                            n_to_enroll     = 15,
                            ...) {

    ## add column
    dta$Pred_Ind_Censor <- dta[[v_censor]]
    dta$Pred_T_Event    <- dta[[v_event]]
    dta$Pred_T_Survival <- NA

    ## to be predicted
    inx_censored_ana <- which(dta[[v_censor]] == censor_ana)
    inx_all          <- c(inx_censored_ana, rep(Inf, n_to_enroll))

    if (0 == length(inx_all)) {
        return(dta)
    }

    ## probabilities of intervals
    inx_dead         <- which(dta$Pred_Ind_Censor == censor_event)
    inx_censored     <- which(dta$Pred_Ind_Censor != censor_event)
    vec_tsurv        <- dta[inx_dead,     "T_Event"]
    vec_tcensor      <- dta[inx_censored, "T_Event"]
    prob_ints        <- predsu_pinterval(vec_tsurv, vec_tcensor, ...)

    ## predict all
    for (j in inx_all) {

        ## new patient
        if (Inf == j) {
            t_enroll <- runif(1, t_ana, t_enrollment)
            t_cond   <- 0
        } else {
            t_enroll <- dta[j, v_enroll]
            t_cond   <- dta[j, v_event]
        }

        ## predict time and censoring
        pred_surv <- predsu_pred(prob_ints, condition_t = t_cond)
        if (pred_surv + t_enroll > t_next_ana) {
            new_censor <- censor_ana
            new_event  <- t_next_ana - t_enroll
        } else if (pred_surv > max(prob_ints[, "T0"])) {
            new_censor <- censor_lfu
            new_event  <- max(prob_ints[, "T0"])
        } else {
            new_censor <- censor_event
            new_event  <- pred_surv
        }

        ## append patient
        if (Inf == j) {
            tmp                <- nrow(dta) + 1
            dta[tmp, ]         <- NA
            dta[tmp, v_enroll] <- t_enroll
        } else {
            tmp <- j
        }

        dta[tmp, "Pred_Ind_Censor"] <- new_censor
        dta[tmp, "Pred_T_Event"]    <- new_event
        dta[tmp, "Pred_T_Survival"] <- pred_surv
    }

    dta
}
