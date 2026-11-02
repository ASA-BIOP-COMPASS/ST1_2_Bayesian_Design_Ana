.prob <- function(x, name, strict = FALSE) {
  if (!is.numeric(x) || !length(x) || any(!is.finite(x)) ||
      any(if (strict) x <= 0 | x >= 1 else x < 0 | x > 1))
    stop(name, ' must contain probabilities ', if (strict) 'strictly within (0,1).' else 'in [0,1].', call. = FALSE)
}
.nonnegative <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 0 || !is.finite(x))
    stop(name, ' must be one finite nonnegative number.', call. = FALSE)
}
.check_counts <- function(y, n) {
  if (length(y) != length(n) || any(!is.finite(y)) || any(!is.finite(n)) ||
      any(y < 0 | y > n) || any(y != floor(y)) || any(n < 0 | n != floor(n)))
    stop('Responses and sample sizes must be nonnegative integers, with y <= n.', call. = FALSE)
}
.method <- function(x) match.arg(x, c('none', 'bonferroni', 'holm', 'hochberg', 'BH'))
.adjust <- function(tail_prob, method, alpha) {
  adjusted <- if (method == 'none') tail_prob else stats::p.adjust(tail_prob, method = method)
  list(adjusted = adjusted, success = adjusted <= alpha)
}
