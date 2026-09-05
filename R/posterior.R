#' Posterior probability of parentage
#'
#' Converts LOD scores into the posterior probability that each candidate is the
#' true parent, given how many candidates were offered and how likely it is that
#' the true parent is among them.
#'
#' @section Why this exists:
#' The two threshold criteria in this package control different things, and
#' neither reports the quantity most studies want. The trio-specific rule of
#' [lod_critical()] bounds the chance of rejecting a true parent, so its recall
#' is stable but its precision falls as the candidate pool grows. The
#' population-level rule of [delta_critical()] bounds the proportion of
#' assignments that are wrong, so its precision is stable but its recall falls
#' instead. Simulation across pool sizes and sampling fractions shows the two
#' failing in exactly opposite directions.
#'
#' A posterior avoids the trade because it is a statement about one assignment
#' rather than about a decision rule. It is also trio-specific and pool-aware at
#' once: the evidence enters through that offspring's own LOD scores, while the
#' size of the pool enters through the prior.
#'
#' @section The model:
#' For one offspring with `n` candidates, either the true parent is among them,
#' with probability `prop_sampled`, or it is not. Given that it is, each
#' candidate is equally likely a priori, so
#' \deqn{\Pr(H_j) = \pi / n, \qquad \Pr(H_0) = 1 - \pi.}
#' The genotypes of the other candidates are uninformative under every
#' hypothesis and cancel, leaving the likelihood ratio of candidate `j` against
#' the unsampled-parent hypothesis equal to \eqn{\exp(\mathrm{LOD}_j)}. Hence
#' \deqn{\Pr(H_j \mid g) = \frac{(\pi/n)\,e^{L_j}}
#'   {(1-\pi) + (\pi/n)\sum_k e^{L_k}}.}
#' The sum in the denominator is what makes the posterior pool-aware: adding
#' candidates dilutes the prior, so the same LOD score buys less confidence in a
#' larger pool. `post_unsampled` reports the complementary probability that no
#' candidate offered is the parent.
#'
#' Computation is on the log scale, so LOD scores of any magnitude are safe.
#'
#' @section Assumptions:
#' Candidates are assumed unrelated to each other and to the offspring except
#' under the hypothesis being tested, and exchangeable a priori. Close relatives
#' of the true parent in the pool break both, and will inflate the posterior of
#' whichever relative scores highest. `prop_sampled` is a genuine prior and the
#' result depends on it; treat it as an assumption to be varied, not a constant
#' to be guessed once.
#'
#' @param x A `"parentage_lod_matrix"` from `parentage_lod(output = "matrix")`,
#'   or a `"parentage_lod"` long table.
#' @param prop_sampled Prior probability that the true parent is among the
#'   candidates offered: one value, or one per offspring.
#' @param n_candidates Number of candidates the prior should be spread over.
#'   Defaults to the number actually scored for each offspring. Set it larger
#'   than the number scored if candidates were screened out before this point,
#'   so that the prior still reflects the pool that was really searched.
#' @param full If `TRUE`, also return the posterior for every candidate rather
#'   than only the best one.
#'
#' @return A `data.frame` of class `"parentage_posterior"`, one row per
#'   offspring: `offspring`, `candidate` (the highest-posterior candidate),
#'   `lod`, `posterior`, `post_unsampled`, `n_candidates`. With `full = TRUE`,
#'   the complete offspring-by-candidate posterior matrix is attached as the
#'   attribute `"posterior_matrix"`.
#'
#' @references
#' Marshall, T.C., Slate, J., Kruuk, L.E.B. & Pemberton, J.M. (1998)
#' Statistical confidence for likelihood-based paternity inference in natural
#' populations. \emph{Molecular Ecology}, \strong{7}, 639--655.
#' \doi{10.1046/j.1365-294x.1998.00374.x}
#'
#' Kalinowski, S.T., Taper, M.L. & Marshall, T.C. (2007) Revising how the
#' computer program CERVUS accommodates genotyping error increases success in
#' paternity assignment. \emph{Molecular Ecology}, \strong{16}, 1099--1106.
#' \doi{10.1111/j.1365-294X.2007.03089.x}
#'
#' Amiri Roudbar, M., Mousavi, S.F., Akbarzadeh, M., Brounts, S.H. & Momen, M.
#' (2025) Pairwise paternity assignment with forward-backward simulations:
#' refining CERVUS using trio-based likelihood and locus-specific error rates.
#' \emph{Ecology and Evolution}, \strong{15}(10), e72230.
#' \doi{10.1002/ece3.72230}
#'
#' @examples
#' sim <- simulate_population(n = 200, n_loci = 12, n_alleles = 8,
#'                            n_offspring = 30, error = 0.01, seed = 8)
#' p <- allele_freqs(sim$genotypes)
#' m <- parentage_lod(sim$genotypes, sim$pedigree$offspring, sim$sires, p,
#'                    error = 0.01, output = "matrix")
#' post <- parentage_posterior(m, prop_sampled = 0.9)
#' head(post)
#' parentage_fdr(post, target = 0.05)
#' @export
parentage_posterior <- function(x, prop_sampled, n_candidates = NULL,
                                full = FALSE) {
  L <- .as_lod_matrix(x)
  no <- nrow(L); nc <- ncol(L)
  pi_ <- as.numeric(prop_sampled)
  if (length(pi_) == 1L) pi_ <- rep(pi_, no)
  if (length(pi_) != no) stop("`prop_sampled` must have length 1 or one value per offspring.", call. = FALSE)
  if (any(is.na(pi_) | pi_ < 0 | pi_ > 1)) stop("`prop_sampled` must lie in [0, 1].", call. = FALSE)

  PM <- if (full) matrix(NA_real_, no, nc, dimnames = dimnames(L)) else NULL
  best <- integer(no); bp <- numeric(no); bl <- numeric(no)
  pu <- numeric(no); nn <- integer(no)

  for (i in seq_len(no)) {
    v <- L[i, ]
    ok <- which(!is.na(v))
    n_i <- if (is.null(n_candidates)) length(ok) else as.integer(n_candidates)
    nn[i] <- n_i
    if (length(ok) == 0L || n_i < 1L) {
      best[i] <- NA_integer_; bp[i] <- NA_real_; bl[i] <- NA_real_; pu[i] <- NA_real_
      next
    }
    lp <- log(pi_[i]) - log(n_i) + v[ok]
    l0 <- log1p(-pi_[i])
    ldenom <- .logsumexp(c(l0, .logsumexp(lp)))
    post <- exp(lp - ldenom)
    pu[i] <- exp(l0 - ldenom)
    j <- which.max(post)
    best[i] <- ok[j]; bp[i] <- post[j]; bl[i] <- v[ok[j]]
    if (full) PM[i, ok] <- post
  }

  out <- data.frame(offspring = rownames(L),
                    candidate = colnames(L)[best],
                    lod = bl, posterior = bp, post_unsampled = pu,
                    n_candidates = nn, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  if (full) attr(out, "posterior_matrix") <- PM
  attr(out, "prop_sampled") <- prop_sampled
  class(out) <- c("parentage_posterior", "data.frame")
  out
}

.logsumexp <- function(v) {
  v <- v[!is.na(v)]
  if (length(v) == 0L) return(-Inf)
  m <- max(v)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(v - m)))
}

.as_lod_matrix <- function(x) {
  if (inherits(x, "parentage_lod_matrix")) return(x$lod)
  if (inherits(x, "parentage_lod")) {
    off <- unique(x$offspring)
    cnd <- unique(if (is.null(x$candidate2)) x$candidate
                  else paste(x$candidate, x$candidate2, sep = " x "))
    L <- matrix(NA_real_, length(off), length(cnd), dimnames = list(off, cnd))
    ci <- if (is.null(x$candidate2)) x$candidate
          else paste(x$candidate, x$candidate2, sep = " x ")
    L[cbind(match(x$offspring, off), match(ci, cnd))] <- x$lod
    return(L)
  }
  if (is.matrix(x)) return(x)
  stop("`x` must be a parentage_lod result or a LOD matrix.", call. = FALSE)
}

#' Choose a posterior cutoff that controls the false discovery rate
#'
#' Given posterior probabilities from [parentage_posterior()], the expected
#' number of wrong assignments among a set of accepted ones is the sum of
#' `1 - posterior` over that set. Accepting offspring in decreasing order of
#' posterior and stopping when that expected proportion reaches `target` gives
#' the largest set whose expected false discovery rate is at most `target`.
#'
#' This is the practical payoff of working with posteriors rather than
#' thresholds: the error rate being controlled is stated directly, in the units
#' a results section reports, and it does not drift when the candidate pool
#' changes size.
#'
#' @param x A `"parentage_posterior"` object, or a numeric vector of posteriors.
#' @param target Target expected false discovery rate.
#'
#' @return A list with the posterior `cutoff`, the number `n_assigned` at that
#'   cutoff, the `expected_false` count, and the realised `fdr`.
#'
#' @examples
#' sim <- simulate_population(n = 200, n_loci = 12, n_alleles = 8,
#'                            n_offspring = 30, error = 0.01, seed = 8)
#' p <- allele_freqs(sim$genotypes)
#' m <- parentage_lod(sim$genotypes, sim$pedigree$offspring, sim$sires, p,
#'                    error = 0.01, output = "matrix")
#' parentage_fdr(parentage_posterior(m, prop_sampled = 0.9), target = 0.05)
#' @export
parentage_fdr <- function(x, target = 0.05) {
  p <- if (inherits(x, "parentage_posterior")) x$posterior else as.numeric(x)
  p <- p[!is.na(p)]
  if (length(p) == 0L) {
    return(list(cutoff = NA_real_, n_assigned = 0L, expected_false = 0, fdr = NA_real_))
  }
  o <- sort(p, decreasing = TRUE)
  fdr <- cumsum(1 - o) / seq_along(o)
  k <- which(fdr <= target)
  if (length(k) == 0L) {
    return(list(cutoff = NA_real_, n_assigned = 0L, expected_false = 0, fdr = NA_real_))
  }
  k <- max(k)
  list(cutoff = o[k], n_assigned = k, expected_false = sum(1 - o[seq_len(k)]),
       fdr = fdr[k])
}

#' @param x A `"parentage_posterior"` object.
#' @param ... Passed to `print.data.frame`.
#' @rdname parentage_posterior
#' @export
print.parentage_posterior <- function(x, ...) {
  cat("<parentage_posterior>", nrow(x), "offspring, prop_sampled =",
      format(attr(x, "prop_sampled")), "\n")
  cat("  posterior >= 0.95: ", sum(x$posterior >= 0.95, na.rm = TRUE),
      ";  >= 0.99: ", sum(x$posterior >= 0.99, na.rm = TRUE), "\n", sep = "")
  cat("  expected wrong among all top candidates: ",
      format(sum(1 - x$posterior, na.rm = TRUE), digits = 4), "\n", sep = "")
  print.data.frame(utils::head(as.data.frame(x), 6L), digits = 4, row.names = FALSE, ...)
  if (nrow(x) > 6L) cat("... ", nrow(x) - 6L, " more rows\n", sep = "")
  invisible(x)
}
