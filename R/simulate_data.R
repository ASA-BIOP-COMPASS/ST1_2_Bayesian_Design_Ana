#' Simulate patient-level binary trial data for analyze_logistic()
#' @param n Named recruited sample sizes, control first.
#' @param response_rates Named marginal-ish baseline stratum response rates
#'   at reference covariates (all dummy variables zero).
#' @param strata Named list of level-probability vectors, e.g.
#'   list(region=c(North=.4,South=.35,West=.25), severity=c(Mild=.6,Severe=.4)).
#' @param stratum_log_or Named numeric list of log-odds effects, one per
#'   nonreference level of each factor. Reference is first listed level.
#' @param seed Random seed.
#' @return Data frame with id, arm, y, and categorical stratification columns.
#' @export
simulate_logistic_data <- function(n, response_rates, strata=list(),
                                   stratum_log_or=list(), seed=123) {
  if (is.null(names(n)) || !identical(names(n),names(response_rates)) ||
      anyDuplicated(names(n)) || length(n)<2) stop('Use matching named arms.',call.=FALSE)
  .check_counts(rep(0,length(n)),n)
  .prob(response_rates,'response_rates', strict=TRUE)
  if (!is.list(strata) || (length(strata) && (is.null(names(strata)) ||
      any(!nzchar(names(strata))) || anyDuplicated(names(strata)))))
    stop('strata must be a named list.',call.=FALSE)
  set.seed(seed)
  arm <- rep(names(n),n)
  dat <- data.frame(id=seq_along(arm),arm=arm,stringsAsFactors=FALSE)
  eta <- stats::qlogis(as.numeric(response_rates[match(arm,names(response_rates))]))
  for (nm in names(strata)) {
    pr <- strata[[nm]]
    if (is.null(names(pr)) || length(pr)<2 || anyDuplicated(names(pr)) ||
        any(!is.finite(pr)) || any(pr<0) || abs(sum(pr)-1)>1e-8)
      stop('Every stratum requires >=2 named levels with probabilities summing to one.',call.=FALSE)
    lev <- names(pr)
    z <- sample(lev,length(arm),replace=TRUE,prob=pr)
    dat[[nm]] <- factor(z,levels=lev)
    effects <- stratum_log_or[[nm]]
    if (is.null(effects)) effects <- rep(0,length(lev)-1L)
    if (!is.numeric(effects) || length(effects)!=length(lev)-1L || any(!is.finite(effects)))
      stop(paste('stratum_log_or for',nm,'must contain one value per nonreference level.'),call.=FALSE)
    eta <- eta + c(0,effects)[match(z,lev)]
  }
  dat$y <- stats::rbinom(length(arm),1,stats::plogis(eta))
  dat
}
