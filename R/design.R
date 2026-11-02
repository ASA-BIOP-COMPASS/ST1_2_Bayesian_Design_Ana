# Validate integer randomization weights, reducing common factors.
# A block of sum(weights / gcd(weights)) participants must be recruited
# to preserve an exact allocation ratio.
.ratio_units <- function(ratios) {
  if (!is.numeric(ratios) || length(ratios) < 2L ||
      any(!is.finite(ratios)) || any(ratios <= 0) ||
      any(ratios != round(ratios)) || any(ratios > .Machine$integer.max))
    stop('ratios must contain at least two positive integer weights.', call. = FALSE)
  if (is.null(names(ratios)) || anyNA(names(ratios)) ||
      any(!nzchar(names(ratios))) || anyDuplicated(names(ratios)))
    stop('ratios must have distinct, non-empty arm names.', call. = FALSE)
  gcd_two <- function(a, b) {
    while (b != 0) { tmp <- b; b <- a %% b; a <- tmp }
    a
  }
  divisor <- Reduce(gcd_two, as.numeric(ratios))
  units <- ratios / divisor
  if (sum(units) > .Machine$integer.max)
    stop('The exact allocation block is too large.', call. = FALSE)
  setNames(as.integer(units), names(ratios))
}

#' Smallest total sample size that preserves exact allocation weights
#' @param ratios Named positive integer allocation weights, placebo first.
#' @return A single integer: the sum of reduced allocation weights.
#' @export
allocation_unit <- function(ratios) sum(.ratio_units(ratios))

#' Generate feasible integer total sample sizes for an allocation ratio
#' @param ratios Named positive integer allocation weights, placebo first.
#' @param min_total Smallest randomized total to consider (inclusive).
#' @param max_total Largest randomized total to consider (inclusive).
#' @param step_blocks Gap between feasible allocation blocks, default 1.
#' @return Integer candidate totals, all exactly divisible by the allocation unit.
#' @examples
#' candidate_totals(c(Placebo=1, Low=1, Mid=1, High=1), 160, 180)
#' candidate_totals(c(Placebo=1, Low=2, Mid=2, High=1), 120, 150)
#' @export
candidate_totals <- function(ratios, min_total, max_total, step_blocks = 1L) {
  block <- allocation_unit(ratios)
  .nonnegative(min_total, 'min_total')
  .nonnegative(max_total, 'max_total')
  if (max_total < min_total) stop('max_total must be >= min_total.', call. = FALSE)
  if (length(step_blocks) != 1L || !is.numeric(step_blocks) ||
      !is.finite(step_blocks) || step_blocks < 1 || step_blocks != floor(step_blocks))
    stop('step_blocks must be a positive integer.', call. = FALSE)
  first <- max(1, ceiling(min_total / block))
  last <- floor(max_total / block)
  if (first > last) return(integer())
  vals <- seq.int(first, last, by = step_blocks) * block
  as.integer(vals)
}

#' Allocate an EXACT integer ratio across placebo and all doses
#' @param total_n Total randomized integer sample size, divisible by allocation_unit(ratios).
#' @param ratios Named positive integer allocation weights (placebo first).
#' @return Named integer sizes, with exact specified allocation ratio.
#' @examples
#' allocate_arms(168, c(Placebo=1, Low=1, Mid=1, High=1))
#' @export
allocate_arms <- function(total_n, ratios = c(Control = 1, Dose1 = 1)) {
  .nonnegative(total_n, 'total_n')
  if (total_n != floor(total_n)) stop('total_n must be an integer.', call. = FALSE)
  units <- .ratio_units(ratios)
  block <- sum(units)
  if (total_n < block || total_n %% block != 0)
    stop(sprintf('total_n must be a positive multiple of %s for the exact allocation ratio. Use candidate_totals().',
                 block), call. = FALSE)
  setNames(as.integer((total_n / block) * units), names(ratios))
}

#' Analyze binary treatment arms against a shared control
#'
#' @param responses Named response counts with control first.
#' @param n Named sample sizes, same order as `responses`.
#' @param control_prior Beta prior shapes for control.
#' @param treatment_prior Common Beta prior shapes for all doses.
#' @param n_post Number of joint posterior Monte Carlo samples.
#' @param margin Clinical superiority margin in risk difference units.
#' @param alpha One-sided posterior tail threshold for decisions.
#' @param multiplicity Adjustment applied to posterior tail areas as a
#'   decision convention, not a frequentist p-value adjustment guarantee.
#' @param seed Random seed, or NULL to preserve current stream.
#' @return Data frame with posterior RD summaries, superiority probability,
#'   adjusted tail score and dose-level decisions.
#' @export
analyze_binary <- function(responses, n, control_prior = c(0.5, 0.5),
                           treatment_prior = c(0.5, 0.5), n_post = 4000,
                           margin = 0, alpha = 0.025,
                           multiplicity = 'hochberg', seed = NULL) {
  .check_counts(responses, n)
  if (length(n) < 2) stop('At least one control and one treatment required.', call. = FALSE)
  if (is.null(names(n)) || !identical(names(n), names(responses))) stop('Use matching named arms.', call. = FALSE)
  if (any(c(control_prior, treatment_prior) <= 0) || length(control_prior) != 2 || length(treatment_prior) != 2)
    stop('Beta prior shapes must be two positive numbers.', call. = FALSE)
  .nonnegative(n_post, 'n_post')
  if (n_post < 100 || n_post != floor(n_post)) stop('n_post must be an integer >=100.', call. = FALSE)
  .prob(alpha, 'alpha', strict = TRUE)
  multiplicity <- .method(multiplicity)
  if (!is.null(seed)) set.seed(seed)
  control <- stats::rbeta(n_post, control_prior[1] + responses[1],
                         control_prior[2] + n[1] - responses[1])
  k <- length(n) - 1L
  out <- vector('list', k)
  for (j in seq_len(k)) {
    idx <- j + 1L
    trt <- stats::rbeta(n_post, treatment_prior[1] + responses[idx],
                        treatment_prior[2] + n[idx] - responses[idx])
    rd <- trt - control
    out[[j]] <- data.frame(arm = names(n)[idx], n = n[idx], responses = responses[idx],
                           rd_mean = mean(rd), rd_lower = unname(stats::quantile(rd, .025)),
                           rd_upper = unname(stats::quantile(rd, .975)),
                           prob_superiority = mean(rd > margin))
  }
  res <- do.call(rbind, out)
  res$tail_score <- 1 - res$prob_superiority
  adj <- .adjust(res$tail_score, multiplicity, alpha)
  res$adjusted_tail_score <- adj$adjusted
  res$success <- adj$success
  rownames(res) <- NULL
  attr(res, 'control_posterior_mean') <- mean(control)
  attr(res, 'multiplicity') <- multiplicity
  res
}

#' Simulate multi-arm design operating characteristics
#'
#' @param n Named sample sizes, first control then dose arms.
#' @param response_rates Named true probabilities aligned with n.
#' @param borrow_ess Historical-control pseudo-sample size.
#' @param historical_rate Prior center for control borrowing, defaults to true control rate.
#' @param base_control_prior Beta baseline prior shapes for placebo, before
#'   any historical borrowing. Default c(0.5, 0.5) (Jeffreys' prior).
#' @param treatment_prior Common Beta prior shapes for active doses. Default
#'   c(0.5, 0.5) (Jeffreys' prior).
#' @param n_sim Number of trial replicates.
#' @param n_post Posterior draws per trial; increase for final calibration.
#' @param multiplicity Adjustment: none, bonferroni, holm, hochberg or BH.
#' @param alpha One-sided decision threshold for posterior tail score.
#' @param margin Risk-difference threshold of clinical relevance.
#' @param seed Seed for reproducibility.
#' @return List: per-arm power, any-dose power, all-dose power, and scenario.
#' @export
simulate_power <- function(n, response_rates, borrow_ess = 0, historical_rate = NULL,
                           n_sim = 500, n_post = 2000, multiplicity = 'hochberg',
                           alpha = 0.025, margin = 0, seed = 123,
                           base_control_prior = c(0.5, 0.5),
                           treatment_prior = c(0.5, 0.5)) {
  if (length(n) < 2L || is.null(names(n)) || anyDuplicated(names(n)) ||
      !identical(names(n), names(response_rates))) stop('Use matching named control and treatment arms.', call. = FALSE)
  .check_counts(rep(0, length(n)), n)
  .prob(response_rates, 'response_rates')
  .nonnegative(borrow_ess, 'borrow_ess')
  .nonnegative(n_sim, 'n_sim')
  if (n_sim < 1 || n_sim != floor(n_sim)) stop('n_sim must be a positive integer.', call. = FALSE)
  if (is.null(historical_rate)) historical_rate <- response_rates[1]
  if (borrow_ess > 0) .prob(historical_rate, 'historical_rate', strict = TRUE)
  multiplicity <- .method(multiplicity)
  if (!is.numeric(base_control_prior) || length(base_control_prior) != 2L ||
      any(!is.finite(base_control_prior)) || any(base_control_prior <= 0))
    stop('base_control_prior must contain two positive finite Beta shapes.', call. = FALSE)
  if (!is.numeric(treatment_prior) || length(treatment_prior) != 2L ||
      any(!is.finite(treatment_prior)) || any(treatment_prior <= 0))
    stop('treatment_prior must contain two positive finite Beta shapes.', call. = FALSE)
  cp <- borrowing_prior(if (borrow_ess > 0) historical_rate else 0.5,
                        borrow_ess, base_prior = base_control_prior)
  set.seed(seed)
  results <- matrix(FALSE, nrow = n_sim, ncol = length(n) - 1L,
                    dimnames = list(NULL, names(n)[-1]))
  for (i in seq_len(n_sim)) {
    y <- stats::rbinom(length(n), n, response_rates)
    names(y) <- names(n)
    a <- analyze_binary(y, n, control_prior = c(cp$alpha, cp$beta),
                        treatment_prior = treatment_prior,
                        n_post = n_post, margin = margin, alpha = alpha,
                        multiplicity = multiplicity)
    results[i, ] <- a$success
  }
  p <- colMeans(results)
  data <- data.frame(arm = names(p), power = as.numeric(p),
                     mc_se = sqrt(as.numeric(p) * (1 - as.numeric(p)) / n_sim))
  list(per_arm = data, any_power = mean(rowSums(results) > 0),
       all_power = mean(rowSums(results) == ncol(results)),
       scenario = list(n = n, response_rates = response_rates, borrow_ess = borrow_ess,
                       historical_rate = historical_rate, multiplicity = multiplicity,
                       base_control_prior = base_control_prior,
                       treatment_prior = treatment_prior,
                       actual_control_prior = c(alpha = cp$alpha, beta = cp$beta),
                       alpha = alpha, margin = margin, n_sim = n_sim, n_post = n_post))
}

#' Search candidate total sample sizes for target operating characteristics
#'
#' @param total_n Candidate overall sample sizes.
#' @param ratios Named allocation ratios.
#' @param response_rates Named true response rates in same order.
#' @param target_power Minimum probability of declaring success.
#' @param target Which power: any, all, or a specified treatment name.
#' @param ... Passed to simulate_power.
#' @return List containing scenario table and smallest passing *candidate*.
#' @export
find_sample_size <- function(total_n, ratios, response_rates, target_power = 0.8,
                             target = 'any', ...) {
  .prob(target_power, 'target_power')
  if (!length(total_n) || any(total_n < 0) || any(total_n != floor(total_n)))
    stop('total_n must contain nonnegative integers.', call. = FALSE)
  if (!identical(names(ratios), names(response_rates))) stop('Match ratio and response names.', call. = FALSE)
  if (!target %in% c('any', 'all', names(ratios)[-1])) stop('Unknown power target.', call. = FALSE)
  unit <- allocation_unit(ratios)
  candidates <- sort(unique(total_n[total_n %% unit == 0 & total_n >= unit]))
  if (!length(candidates)) stop('No feasible total_n; use candidate_totals().', call. = FALSE)
  ans <- lapply(candidates, function(size) {
    n <- allocate_arms(size, ratios)
    out <- simulate_power(n = n, response_rates = response_rates, ...)
    score <- if (target == 'any') out$any_power else if (target == 'all') out$all_power else
      out$per_arm$power[match(target, out$per_arm$arm)]
    data.frame(total_n = size, power = score, passes = score >= target_power,
               as.list(n), check.names = FALSE)
  })
  tab <- do.call(rbind, ans)
  list(results = tab, minimum_candidate = if (any(tab$passes))
    tab[which(tab$passes)[1], , drop = FALSE] else NULL,
    target = target, target_power = target_power)
}
