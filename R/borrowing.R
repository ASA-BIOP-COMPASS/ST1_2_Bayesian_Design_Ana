#' Create a conjugate historical control borrowing prior
#'
#' @param historical_rate Historical or anticipated control response rate in (0,1).
#' @param ess Borrowed information (pseudo-patient total), nonnegative.
#' @param base_prior Baseline Beta prior shapes. Defaults to Jeffreys' prior.
#' @return Named list with prior shapes and borrowing ESS.
#' @export
borrowing_prior <- function(historical_rate, ess = 0, base_prior = c(0.5, 0.5)) {
  .prob(historical_rate, 'historical_rate', strict = TRUE)
  .nonnegative(ess, 'ess')
  if (length(base_prior) != 2L || any(!is.finite(base_prior)) || any(base_prior <= 0))
    stop('base_prior must be two positive shapes.', call. = FALSE)
  list(alpha = base_prior[1] + ess * historical_rate,
       beta = base_prior[2] + ess * (1 - historical_rate),
       borrowed_ess = ess, beta_total = sum(base_prior) + ess,
       historical_rate = historical_rate)
}

#' Compute a prespecified cap on borrowing from historical controls
#'
#' @param n_control Planned concurrent control count.
#' @param historical_ess Available borrowing ESS from the historical prior.
#' @param max_fraction Maximum fraction of combined control information
#'   attributable to borrowing, between 0 and 1 (exclusive of 1).
#' @param max_absolute Maximum borrowed pseudo-subjects.
#' @param max_ratio Maximum borrowed pseudo-subjects per concurrent control.
#' @return One-row data frame. These are policy caps, not a statistically
#'   validated dynamic-borrowing rule. For n_control=0, fraction/ratio caps are 0.
#' @export
borrowing_limits <- function(n_control, historical_ess, max_fraction = 0.5,
                             max_absolute = Inf, max_ratio = Inf) {
  .nonnegative(n_control, 'n_control')
  .nonnegative(historical_ess, 'historical_ess')
  if (length(max_fraction) != 1 || !is.finite(max_fraction) || max_fraction < 0 || max_fraction >= 1)
    stop('max_fraction must be in [0,1).', call. = FALSE)
  if (length(max_absolute) != 1 || is.na(max_absolute) || max_absolute < 0)
    stop('max_absolute must be nonnegative.', call. = FALSE)
  if (length(max_ratio) != 1 || is.na(max_ratio) || max_ratio < 0)
    stop('max_ratio must be nonnegative.', call. = FALSE)
  cap <- min(historical_ess, max_absolute, n_control * max_ratio,
             n_control * max_fraction / (1 - max_fraction))
  data.frame(available_ess = historical_ess, allowed_ess = cap,
             restricted_ess = historical_ess - cap,
             borrowed_fraction = if (n_control + cap > 0) cap / (n_control + cap) else NA_real_)
}

#' Estimate historical control MAP prior ESS via RBesT
#'
#' @param data Data frame with columns `study`, `n`, `r` (responders).
#' @param seed Random seed.
#' @param ... Passed to RBesT::gMAP.
#' @return List of MAP model, mixture approximation, and three ESS estimates.
#' @export
historical_map_ess <- function(data, seed = 123, ...) {
  if (!requireNamespace('RBesT', quietly = TRUE)) stop('Install RBesT first.', call. = FALSE)
  if (!all(c('study', 'n', 'r') %in% names(data))) stop('Need columns study, n, r.', call. = FALSE)
  if (any(data$n <= 0) || any(data$r < 0 | data$r > data$n)) stop('Invalid historical counts.', call. = FALSE)
  set.seed(seed)
  fit <- RBesT::gMAP(cbind(r, n - r) ~ 1 | study, data = data,
                     tau.dist = 'HalfNormal', tau.prior = 1, beta.prior = 2,
                     family = binomial, ...)
  mix <- RBesT::automixfit(fit)
  vals <- vapply(c('elir', 'moment', 'morita'),
                 function(x) as.numeric(RBesT::ess(mix, method = x)), numeric(1))
  list(map = fit, mixture = mix, ess = vals)
}
