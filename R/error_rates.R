#' Estimate locus-specific genotyping error rates
#'
#' Estimates a separate genotyping error rate for each locus from known
#' parent-offspring pairs, following the approach of Amiri Roudbar et al.
#' (2025), section 2.6: mismatches between offspring and parents are counted
#' across the population, and pairs carrying more than `max_pair_mismatch`
#' mismatching loci are set aside on the grounds that they are more likely to be
#' wrongly assigned relationships than repeatedly mistyped ones.
#'
#' @section Estimators:
#' The published description gives the numerator (a count of mismatching loci)
#' but not the denominator that converts it into a rate, so two estimators are
#' offered and the choice is explicit.
#'
#' `"naive"` returns the raw proportion of scored pairs that mismatch at the
#' locus. This is biased low, because a genotyping error only becomes visible
#' when the replacement genotype happens to be Mendelian-incompatible with the
#' other member of the pair, and at a locus with few common alleles most errors
#' stay hidden.
#'
#' `"corrected"`, the default, inverts that detection probability. Under the
#' random-genotype-replacement model exactly one member of a pair is mistyped
#' with probability \eqn{2e(1-e)}, and such an error is detected with
#' probability
#' \deqn{d = 1 - \sum_{g_p} \sum_{g_o} P(g_p) P(g_o) 1\{T(g_o \mid g_p) > 0\},}
#' computed exactly from `freqs`. The observed mismatch rate \eqn{r} then
#' satisfies \eqn{r = 2 e (1-e) d}, which is solved for the smaller root
#' \eqn{e = [1 - \sqrt{1 - 2r/d}]/2}. This inversion is the package's own
#' derivation from the published error model, not a published equation.
#'
#' @param g A `"genotypes"` object.
#' @param pairs A `data.frame` or `matrix` with two columns giving known
#'   parent and offspring identifiers, one row per pair. Trios should be supplied
#'   as two rows.
#' @param freqs An `"allele_freqs"` object. Required for
#'   `method = "corrected"`.
#' @param method `"corrected"` or `"naive"`.
#' @param max_pair_mismatch Pairs mismatching at more than this many loci are
#'   excluded as probable non-relatives. `NULL`, the default, keeps every pair
#'   supplied, on the assumption that the pedigree given is correct. Amiri
#'   Roudbar et al. (2025) describe retaining pairs with a single typing error,
#'   which corresponds to `max_pair_mismatch = 1`; note that this discards
#'   exactly the pairs carrying the most errors and so biases the estimate
#'   downward, increasingly so as the true rate rises. Use it when the pedigree
#'   itself is uncertain and the bias is the lesser problem.
#' @param min_error,max_error Bounds applied to every estimate. The lower bound
#'   keeps the likelihood well behaved at loci where no mismatch was seen, since
#'   an error rate of exactly zero makes a single mismatch an absolute exclusion.
#'
#' @return A `data.frame` of class `"error_rates"` with one row per locus:
#'   `locus`, `n_pairs` (retained pairs scored at the locus), `n_mismatch`,
#'   `mismatch_rate`, `detectable` (the probability \eqn{d}, `NA` for the naive
#'   method) and `error`. The `error` column can be passed directly as the
#'   `error` argument of [parentage_lod()].
#'
#' @references
#' Amiri Roudbar, M., Mousavi, S.F., Akbarzadeh, M., Brounts, S.H. & Momen, M.
#' (2025) Pairwise paternity assignment with forward-backward simulations:
#' refining CERVUS using trio-based likelihood and locus-specific error rates.
#' \emph{Ecology and Evolution}, \strong{15}(10), e72230.
#' \doi{10.1002/ece3.72230}
#'
#' Kalinowski, S.T., Taper, M.L. & Marshall, T.C. (2007) Revising how the
#' computer program CERVUS accommodates genotyping error increases success in
#' paternity assignment. \emph{Molecular Ecology}, \strong{16}, 1099--1106.
#' \doi{10.1111/j.1365-294X.2007.03089.x}
#'
#' @examples
#' sim <- simulate_population(n = 120, n_loci = 8, n_alleles = 6,
#'                            n_offspring = 120, error = 0.03, seed = 11)
#' p <- allele_freqs(sim$genotypes)
#' pairs <- rbind(
#'   data.frame(parent = sim$pedigree$sire, offspring = sim$pedigree$offspring),
#'   data.frame(parent = sim$pedigree$dam,  offspring = sim$pedigree$offspring)
#' )
#' er <- estimate_error_rates(sim$genotypes, pairs, p)
#' er
#' @export
estimate_error_rates <- function(g, pairs, freqs = NULL,
                                 method = c("corrected", "naive"),
                                 max_pair_mismatch = NULL,
                                 min_error = 1e-4, max_error = 0.5) {
  stopifnot(inherits(g, "genotypes"))
  method <- match.arg(method)
  if (method == "corrected" && is.null(freqs)) {
    stop("`freqs` is required for method = 'corrected'.", call. = FALSE)
  }
  if (is.matrix(pairs)) pairs <- as.data.frame(pairs, stringsAsFactors = FALSE)
  if (!is.data.frame(pairs) || ncol(pairs) < 2L) {
    stop("`pairs` must be a two-column table of parent and offspring identifiers.", call. = FALSE)
  }
  par_id <- .as_ids(g, as.character(pairs[[1L]]), "pairs[, 1]")
  off_id <- .as_ids(g, as.character(pairs[[2L]]), "pairs[, 2]")
  keep <- !is.na(par_id) & !is.na(off_id)
  par_id <- par_id[keep]; off_id <- off_id[keep]
  if (length(par_id) == 0L) stop("No usable parent-offspring pairs.", call. = FALSE)

  pi_ <- match(par_id, g$ids); oi <- match(off_id, g$ids)
  nl <- n_loci(g)
  fr <- if (is.null(freqs)) allele_freqs(g) else freqs
  .check_freqs(fr, g)

  mm <- matrix(NA, length(pi_), nl)
  for (l in seq_len(nl)) {
    mm[, l] <- .mismatch_locus(g$a1[oi, l], g$a2[oi, l], g$a1[pi_, l], g$a2[pi_, l],
                               fr[[l]], type = "both_unknown")
  }
  per_pair <- rowSums(mm, na.rm = TRUE)
  use <- if (is.null(max_pair_mismatch)) rep(TRUE, length(per_pair)) else per_pair <= max_pair_mismatch
  if (!any(use)) {
    stop("No pair had at most ", max_pair_mismatch,
         " mismatching loci; check the pedigree or raise `max_pair_mismatch`.", call. = FALSE)
  }
  mm <- mm[use, , drop = FALSE]

  n_pairs <- colSums(!is.na(mm))
  n_mis <- colSums(mm, na.rm = TRUE)
  rate <- ifelse(n_pairs > 0L, n_mis / n_pairs, NA_real_)

  if (method == "naive") {
    det <- rep(NA_real_, nl)
    est <- rate
  } else {
    det <- vapply(seq_len(nl), function(l) .detect_prob(fr[[l]]), numeric(1))
    est <- vapply(seq_len(nl), function(l) {
      r <- rate[l]; d <- det[l]
      if (is.na(r) || is.na(d) || d <= 0) return(NA_real_)
      z <- 1 - 2 * r / d
      if (z <= 0) max_error else (1 - sqrt(z)) / 2
    }, numeric(1))
  }
  est[is.na(est)] <- min_error
  est <- pmin(pmax(est, min_error), max_error)

  out <- data.frame(locus = g$loci, n_pairs = n_pairs, n_mismatch = n_mis,
                    mismatch_rate = rate, detectable = det, error = est,
                    stringsAsFactors = FALSE)
  attr(out, "method") <- method
  attr(out, "n_pairs_used") <- sum(use)
  attr(out, "n_pairs_dropped") <- sum(!use)
  class(out) <- c("error_rates", "data.frame")
  out
}

# Probability that a random-genotype replacement in one member of a
# parent-offspring pair produces a visible Mendelian incompatibility.
.detect_prob <- function(p) {
  k <- length(p)
  gi <- unlist(lapply(seq_len(k), function(i) rep(i, k - i + 1L)))
  gj <- unlist(lapply(seq_len(k), function(i) seq.int(i, k)))
  Pg <- geno_freq(gi, gj, p)
  compat <- 0
  for (m in seq_along(gi)) {
    t <- trans_prob(gi, gj, gi[m], gj[m], p)
    compat <- compat + Pg[m] * sum(Pg[t > 0])
  }
  1 - compat
}

#' @param x An `"error_rates"` object.
#' @param ... Passed to `print.data.frame`.
#' @rdname estimate_error_rates
#' @export
print.error_rates <- function(x, ...) {
  cat("<error_rates>", nrow(x), "loci, method =", attr(x, "method"), "\n")
  cat("  pairs used: ", attr(x, "n_pairs_used"),
      " (dropped ", attr(x, "n_pairs_dropped"), " over the mismatch limit)\n", sep = "")
  cat("  error rate: median ", format(stats::median(x$error), digits = 3),
      ", range ", paste(format(range(x$error), digits = 3), collapse = " - "), "\n", sep = "")
  print.data.frame(as.data.frame(x), digits = 3, ...)
  invisible(x)
}
