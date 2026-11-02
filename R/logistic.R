# Shortest contiguous interval containing the requested fraction of MCMC draws.
# This is a simulation-based approximation to a unimodal HPD interval.
.hpd_draw_interval <- function(draws, mass = 0.95) {
  if (!is.numeric(draws) || length(draws) < 2L || any(!is.finite(draws)))
    stop('HPD interval requires at least two finite posterior draws.', call. = FALSE)
  z <- sort(as.numeric(draws))
  n <- length(z)
  gap <- ceiling(mass * n) - 1L
  gap <- max(1L, min(n - 1L, gap))
  left <- seq_len(n - gap)
  right <- left + gap
  idx <- which.min(z[right] - z[left])
  unname(c(z[left[idx]], z[right[idx]]))
}

#' Bayesian logistic analysis with shared covariate effects and control borrowing
#'
#' @param data Patient-level data frame with binary `outcome`, arm and covariates.
#' @param outcome Name of 0/1 outcome column.
#' @param arm Name of treatment arm column.
#' @param control Character control arm label.
#' @param covariates Character vector of covariate names. Numeric and categorical (factor/character) columns; factors may have multiple levels.
#' @param center `overall`, `arm`, or `none` numeric covariate centering.
#' @param control_prior Beta prior shapes for the reference control response
#'   probability. Defaults to c(0.5, 0.5) (Jeffreys' Bernoulli prior),
#'   without automatic historical borrowing. Supply explicit positive shapes
#'   to use an informative prior.
#' @param intercept_sd Prior SD on each active-dose log-odds intercept.
#'   Default 100 gives a proper, very diffuse Normal(0, 100^2) prior.
#' @param covariate_sd Prior SD on shared covariate coefficients.
#'   Default 100 gives a proper, very diffuse Normal(0, 100^2) prior.
#' @param n_iter Total Metropolis iterations per chain.
#' @param burnin Burn-in iterations per chain.
#' @param thin Retain every nth iteration after burn-in.
#' @param n_chains Independent chains initialized separately.
#' @param proposal_sd Scalar or parameter-length random-walk proposal SD.
#' @param margin RD threshold for superiority.
#' @param alpha One-sided posterior tail threshold.
#' @param multiplicity Posterior-tail multiplicity decision convention.
#' @param seed Seed for reproducibility.
#' @param interval_method Credible-interval reporting method: 'both' (default),
#'   'HPD' (approximate shortest posterior-draw interval), or 'equal_tail'.
#'   When 'both', rd_lower/rd_upper retain equal-tailed bounds for backward
#'   compatibility; separate HPD and equal-tailed columns are also returned.
#' @param credible_level Posterior probability contained in each credible
#'   interval; default 0.95.
#' @param sas_compatible_summary If TRUE, also provide a one-row
#'   `sas_summary` table with SAS-style pP_ref_mean, pD1_ref_mean,
#'   RD_D1_mean, RD_D1_LCL, RD_D1_UCL and Pr_Sup_D1 columns for all doses.
#'   The LCL and UCL columns always use HPD intervals, matching the
#'   historical SAS output convention; Pr_Sup means Pr(RD > 0).
#' @return A bbd_logistic object with summaries, retained posterior and MCMC diagnostics.
#' @details Direct independent Metropolis random-walk implementation; designed as a
#'   reference implementation for moderate-sized datasets, not a validated replacement
#'   for SAS PROC MCMC or production HMC. Posterior risk differences are calculated
#'   at covariate values of zero on the *centered* design matrix scale. When
#'   center='arm', the reference populations differ by arm and the RD is not a
#'   standardized marginal contrast. HPD intervals use the shortest contiguous
#'   interval among retained posterior draws; when posteriors are multimodal,
#'   this is not necessarily a disjoint highest-density region. SAS PROC MCMC
#'   and this custom sampler may produce different Monte Carlo approximations;
#'   compare posterior means, superiority probabilities and convergence before
#'   attributing differences to interval definitions. MCMC diagnostics must be
#'   reviewed.
#' @export
analyze_logistic <- function(data, outcome = 'y', arm = 'arm', control = 'Placebo',
                             covariates = character(), center = 'overall',
                             control_prior = c(0.5, 0.5),
                             intercept_sd = 100, covariate_sd = 100,
                             n_iter = 12000, burnin = 2000, thin = 5,
                             n_chains = 2, proposal_sd = 0.15,
                             margin = 0, alpha = 0.025, multiplicity = 'hochberg',
                             seed = 123, interval_method = 'both',
                             credible_level = 0.95,
                             sas_compatible_summary = FALSE) {
  needed <- c(outcome, arm, covariates)
  if (!all(needed %in% names(data)) || anyNA(data[, needed, drop = FALSE]))
    stop('Missing requested columns or values.', call. = FALSE)
  y <- data[[outcome]]
  if (!is.numeric(y) || any(!y %in% c(0, 1))) stop('Outcome must be numeric 0/1.', call. = FALSE)
  center <- match.arg(center, c('overall', 'arm', 'none'))
  multiplicity <- .method(multiplicity)
  interval_method <- match.arg(interval_method, c('both', 'HPD', 'equal_tail'))
  if (!is.numeric(credible_level) || length(credible_level) != 1L ||
      !is.finite(credible_level) || credible_level <= 0 || credible_level >= 1)
    stop('credible_level must be a single number strictly between 0 and 1.', call. = FALSE)
  if (!is.logical(sas_compatible_summary) || length(sas_compatible_summary) != 1L ||
      is.na(sas_compatible_summary))
    stop('sas_compatible_summary must be TRUE or FALSE.', call. = FALSE)
  if (length(control_prior) != 2 || any(!is.finite(control_prior)) || any(control_prior <= 0))
    stop('control_prior must be two positive Beta shapes.', call. = FALSE)
  if (!is.numeric(intercept_sd) || length(intercept_sd) != 1 || intercept_sd <= 0 || !is.finite(intercept_sd) ||
      !is.numeric(covariate_sd) || length(covariate_sd) != 1 || covariate_sd <= 0 || !is.finite(covariate_sd))
    stop('Prior SD must be positive and finite.', call. = FALSE)
  if (n_iter <= burnin || burnin < 0 || thin < 1 || n_chains < 1 ||
      any(c(n_iter, burnin, thin, n_chains) != floor(c(n_iter, burnin, thin, n_chains))))
    stop('Invalid MCMC iteration settings.', call. = FALSE)
  .prob(alpha, 'alpha', strict = TRUE)
  arms <- c(control, setdiff(unique(as.character(data[[arm]])), control))
  if (length(arms) < 2 || !control %in% data[[arm]]) stop('Control and >=1 dose required.', call. = FALSE)
  arm_index <- match(as.character(data[[arm]]), arms)
  xx <- matrix(numeric(0), nrow = nrow(data), ncol = 0)
  if (length(covariates)) {
    covdf <- data[, covariates, drop=FALSE]
    for (nm in covariates) {
      if (is.character(covdf[[nm]]) || is.logical(covdf[[nm]]))
        covdf[[nm]] <- factor(covdf[[nm]])
    }
    # model.matrix expands every multi-level factor to k-1 indicators.
    xx <- stats::model.matrix(~ ., data=covdf)[,-1,drop=FALSE]
    if (center == 'overall' && ncol(xx)) xx <- sweep(xx,2,colMeans(xx),'-')
    if (center == 'arm' && ncol(xx)) for (i in unique(arm_index)) {
      z <- arm_index == i
      xx[z,] <- sweep(xx[z,,drop=FALSE],2,colMeans(xx[z,,drop=FALSE]),'-')
    }
  }
  design <- cbind(stats::model.matrix(~ factor(arm_index, levels = seq_along(arms)) - 1), xx)
  k <- length(arms)
  d <- ncol(design)
  if (length(proposal_sd) == 1) proposal_sd <- rep(proposal_sd, d)
  if (length(proposal_sd) != d || any(!is.finite(proposal_sd)) || any(proposal_sd <= 0))
    stop('proposal_sd must be positive scalar or vector matching parameter count.', call. = FALSE)
  # Stable likelihood uses dbinom to avoid log(0) on numerical extremes.
  logpost <- function(theta) {
    p0 <- stats::plogis(theta[1])
    prior <- stats::dbeta(p0, control_prior[1], control_prior[2], log = TRUE) +
      log(p0) + log1p(-p0) # Jacobian of logit transform
    if (k > 1) prior <- prior + sum(stats::dnorm(theta[2:k], 0, intercept_sd, log = TRUE))
    if (d > k) prior <- prior + sum(stats::dnorm(theta[(k + 1):d], 0, covariate_sd, log = TRUE))
    prior + sum(stats::dbinom(y, size = 1, prob = stats::plogis(drop(design %*% theta)), log = TRUE))
  }
  set.seed(seed)
  saved_indices <- seq.int(burnin + 1L, n_iter)
  saved_indices <- saved_indices[(seq_along(saved_indices) - 1L) %% thin == 0]
  traces <- vector('list', n_chains)
  acceptance <- numeric(n_chains)
  for (chain in seq_len(n_chains)) {
    theta <- c(stats::qlogis((control_prior[1] + sum(y[arm_index == 1])) /
                                   (sum(control_prior) + sum(arm_index == 1))),
               rep(0, d - 1L)) + stats::rnorm(d, sd = 0.05)
    lp <- logpost(theta)
    draws <- matrix(NA_real_, nrow = length(saved_indices), ncol = d)
    accepted <- 0L
    s <- 0L
    for (i in seq_len(n_iter)) {
      # Component-wise MH tends to mix better than a single high-dimensional proposal.
      for (j in seq_len(d)) {
        test <- theta
        test[j] <- test[j] + stats::rnorm(1, sd = proposal_sd[j])
        lp_new <- logpost(test)
        if (is.finite(lp_new) && log(stats::runif(1)) < lp_new - lp) {
          theta <- test
          lp <- lp_new
          accepted <- accepted + 1L
        }
      }
      if (i > burnin && (i - burnin - 1L) %% thin == 0L) {
        s <- s + 1L
        draws[s, ] <- theta
      }
    }
    traces[[chain]] <- draws
    acceptance[chain] <- accepted / (n_iter * d)
  }
  # Use expanded model-matrix column names: a factor with K levels creates K-1 coefficients.
  parameter_names <- c(paste0('logit_', make.names(arms)),
                       paste0('beta_', colnames(xx)))
  traces <- lapply(traces, function(z) {colnames(z) <- parameter_names; z})
  pool <- do.call(rbind, traces)
  reference_p <- stats::plogis(pool[, seq_len(k), drop = FALSE])
  colnames(reference_p) <- arms
  diffs <- sweep(reference_p[, -1, drop = FALSE], 1, reference_p[, 1], '-')
  res <- lapply(seq_len(k - 1L), function(j) {
    val <- diffs[, j]
    tails <- (1 - credible_level) / 2
    et <- unname(stats::quantile(val, c(tails, 1 - tails)))
    hp <- .hpd_draw_interval(val, credible_level)
    # Keep the longstanding rd_lower/rd_upper interface for existing users.
    primary <- if (interval_method == 'HPD') hp else et
    out <- data.frame(arm = arms[j + 1L], rd_mean = mean(val),
                      rd_lower = primary[1], rd_upper = primary[2],
                      prob_superiority = mean(val > margin))
    if (interval_method %in% c('both', 'HPD')) {
      out$hpd_lower <- hp[1]
      out$hpd_upper <- hp[2]
    }
    if (interval_method %in% c('both', 'equal_tail')) {
      out$equal_tail_lower <- et[1]
      out$equal_tail_upper <- et[2]
    }
    out
  })
  res <- do.call(rbind, res)
  res$tail_score <- 1 - res$prob_superiority
  a <- .adjust(res$tail_score, multiplicity, alpha)
  res$adjusted_tail_score <- a$adjusted
  res$success <- a$success
  rownames(res) <- NULL
  diagnostic <- mcmc_diagnostics(traces)
  sas_table <- NULL
  if (sas_compatible_summary) {
    # SAS-style wide output: dose indices follow the observed arm order.
    sas_values <- c(pP_ref_mean = unname(mean(reference_p[, 1])))
    for (j in seq_len(k - 1L)) {
      sas_values[paste0('pD', j, '_ref_mean')] <- mean(reference_p[, j + 1L])
    }
    for (j in seq_len(k - 1L)) {
      v <- diffs[, j]
      hp <- .hpd_draw_interval(v, credible_level)
      sas_values[paste0('RD_D', j, '_mean')] <- mean(v)
      sas_values[paste0('RD_D', j, '_LCL')] <- hp[1]
      sas_values[paste0('RD_D', j, '_UCL')] <- hp[2]
    }
    for (j in seq_len(k - 1L)) {
      sas_values[paste0('Pr_Sup_D', j)] <- mean(diffs[, j] > 0)
    }
    sas_table <- as.data.frame(as.list(sas_values), check.names = FALSE)
  }
  out <- list(summary = res, reference_means = colMeans(reference_p),
              posterior = pool, reference_prob_draws = reference_p,
              chains = traces, mcmc = diagnostic, acceptance = acceptance,
              sas_summary = sas_table, interval_method = interval_method,
              credible_level = credible_level,
              arms = arms, covariates = covariates, center = center,
              control_prior = control_prior, intercept_sd = intercept_sd,
              covariate_sd = covariate_sd, multiplicity = multiplicity)
  class(out) <- 'bbd_logistic'
  out
}

#' Basic diagnostic summaries for multiple independent MCMC chains
#'
#' @param chains List of parameter matrices (iteration x parameter).
#' @return Data frame with approximate split-Rhat and chain-averaged ESS.
#' @details Basic Geyer initial-positive-pairs autocorrelation ESS and split Rhat.
#'   Not a rank-normalized diagnostic; for final work use posterior or bayesplot.
#' @export
mcmc_diagnostics <- function(chains) {
  m <- length(chains)
  p <- ncol(chains[[1]])
  n <- min(vapply(chains, nrow, integer(1)))
  ans <- vector('list', p)
  for (j in seq_len(p)) {
    x <- sapply(chains, function(z) z[seq_len(n), j])
    if (is.null(dim(x))) x <- matrix(x, ncol = 1)
    # split chain for Rhat
    h <- floor(n / 2)
    halves <- if (h >= 2) do.call(cbind, lapply(seq_len(m), function(i)
      cbind(x[seq_len(h), i], x[(n - h + 1L):n, i]))) else matrix(numeric(0), 0, 0)
    rhat <- NA_real_
    if (nrow(halves) >= 2 && ncol(halves) >= 2) {
      W <- mean(apply(halves, 2, stats::var))
      B <- nrow(halves) * stats::var(colMeans(halves))
      if (W > 0) rhat <- sqrt(((nrow(halves) - 1) / nrow(halves) * W + B / nrow(halves)) / W)
    }
    # Conservative chain-wise autocorrelation truncation.
    per_chain <- vapply(seq_len(m), function(i) {
      if (n < 5 || stats::sd(x[, i]) == 0) return(as.numeric(n))
      ac <- as.numeric(stats::acf(x[, i], plot = FALSE, lag.max = min(500L, n - 1L))$acf)[-1]
      pairs <- if (length(ac) >= 2) ac[seq(1, length(ac) - 1L, 2L)] +
        ac[seq(2, length(ac), 2L)] else numeric()
      z <- which(pairs <= 0)
      positive <- if (length(z)) pairs[seq_len(z[1] - 1L)] else pairs
      max(1, n / (1 + 2 * sum(positive)))
    }, numeric(1))
    ans[[j]] <- data.frame(parameter = colnames(chains[[1]])[j],
                           rhat = rhat, ess_mcmc = sum(per_chain))
  }
  do.call(rbind, ans)
}

#' @export
print.bbd_logistic <- function(x, ...) {
  cat('Bayesian logistic analysis (reference Metropolis implementation)\n')
  print(x$summary, row.names = FALSE)
  cat('Acceptance fractions:', paste(round(x$acceptance, 3), collapse = ', '), '\n')
  invisible(x)
}
