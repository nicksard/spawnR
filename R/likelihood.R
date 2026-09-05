#' Per-locus likelihood ratio under the genotyping-error model
#'
#' Computes the single-locus LOD contribution
#' \eqn{\ln[L(H_1) / L(H_2)]} for the three parentage configurations, under the
#' random-genotype-replacement model of genotyping error with a locus-specific
#' error rate.
#'
#' @section Error model:
#' Following Marshall et al. (1998) and Kalinowski et al. (2007), an observed
#' genotype equals the true genotype with probability \eqn{1 - e}, and with
#' probability \eqn{e} is replaced by a genotype drawn at random under
#' Hardy-Weinberg proportions:
#' \deqn{\Pr(g \mid G) = (1 - e)\,1\{G = g\} + e\,P(g).}
#' Amiri Roudbar et al. (2025) retain this model but allow \eqn{e} to differ
#' among loci, which is what `error` permits here.
#'
#' @section Expansion:
#' Marginalising the unobserved true genotypes gives one term per subset of
#' individuals that was mistyped. With three individuals involved (a known or
#' alleged second parent, the candidate, and the offspring):
#' \deqn{L(H_1) \propto (1-e)^3 T(o \mid m, a)
#'   + e(1-e)^2 [ T(o \mid m) + T(o \mid a) + P(o) ]
#'   + 3 e^2 (1-e) P(o) + e^3 P(o)}
#' \deqn{L(H_2) \propto (1-e)^3 T(o \mid m)
#'   + e(1-e)^2 [ T(o \mid m) + 2 P(o) ] + 3 e^2 (1-e) P(o) + e^3 P(o)}
#' for `type = "one_known"`; for `type = "pair"` the numerator is the same with
#' the alleged mother in place of the known mother and the denominator collapses
#' to \eqn{P(o)}; and for `type = "both_unknown"`, with only two individuals,
#' \deqn{L(H_1) \propto (1-e)^2 T(o \mid a) + 2 e(1-e) P(o) + e^2 P(o),
#'   \qquad L(H_2) \propto P(o).}
#' The common factors \eqn{P(g_m) P(g_a)} cancel in the ratio and are omitted.
#' Every weight set sums to one over the offspring genotype, so both hypotheses
#' are proper distributions; this is asserted in the package's regression tests.
#'
#' @section Relationship to the published coefficients:
#' The multinomial multipliers above (`3` on the \eqn{e^2(1-e)} term, `2` on the
#' \eqn{e(1-e)} term, and the second \eqn{P(o)} in the `"one_known"`
#' denominator) are re-derived from the error model rather than transcribed. The
#' corresponding coefficients printed by Amiri Roudbar et al. (2025) do not sum
#' to one and cannot be a probability decomposition; a corrigendum to Kalinowski
#' et al. (2007) likewise records typesetting errors in that paper's appendix
#' without any change to the underlying method. See the derivation note shipped
#' with the package, `system.file("notes", "derivation.md", package =
#' "spawnR")`, for the full argument.
#'
#' @param o1,o2 Integer allele codes of the offspring genotype at the locus.
#' @param c1,c2 Integer allele codes of the candidate parent's genotype. For
#'   `type = "pair"` this is the candidate father.
#' @param k1,k2 Integer allele codes of the known parent (`type = "one_known"`)
#'   or of the candidate mother (`type = "pair"`). Ignored when
#'   `type = "both_unknown"`.
#' @param p Numeric vector of allele frequencies at the locus.
#' @param error Genotyping error rate at the locus, in `[0, 1]`.
#' @param type One of `"one_known"`, `"both_unknown"` or `"pair"`.
#'
#' @return A numeric vector of per-locus LOD contributions. `NA` where a
#'   required genotype is missing; `-Inf` at a locus that excludes the candidate
#'   outright (only reachable when `error` is exactly zero).
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
#' p <- c(0.4, 0.3, 0.2, 0.1)
#' # Candidate 1/2, offspring 1/3, no second parent, 1% error
#' lod_locus(1L, 3L, 1L, 2L, p = p, error = 0.01, type = "both_unknown")
#' # The same comparison with the mother 3/3 known
#' lod_locus(1L, 3L, 1L, 2L, k1 = 3L, k2 = 3L, p = p,
#'           error = 0.01, type = "one_known")
#' @export
lod_locus <- function(o1, o2, c1, c2, p, error, k1 = NULL, k2 = NULL,
                      type = c("one_known", "both_unknown", "pair")) {
  lk <- .locus_lik(o1, o2, c1, c2, p, error, k1, k2, match.arg(type))
  .safe_log_ratio(lk$num, lk$den)
}

# Unnormalised L(H1) and L(H2) at one locus, with the common factors
# P(g_k) P(g_c) divided out. Both are proper probability distributions over the
# offspring genotype; `test-likelihood.R` asserts that they sum to one.
.locus_lik <- function(o1, o2, c1, c2, p, error, k1 = NULL, k2 = NULL,
                       type = "one_known") {
  e <- error[1L]
  if (is.na(e) || e < 0 || e > 1) stop("`error` must lie in [0, 1].", call. = FALSE)

  Po  <- geno_freq(o1, o2, p)
  T_c <- trans_prob(o1, o2, c1, c2, p)

  if (type == "both_unknown") {
    u0 <- (1 - e)^2; u1 <- e * (1 - e); u2 <- e^2
    num <- u0 * T_c + (2 * u1 + u2) * Po
    den <- Po
  } else {
    if (is.null(k1) || is.null(k2)) {
      stop("`k1` and `k2` are required for type = \'", type, "\'.", call. = FALSE)
    }
    w0 <- (1 - e)^3; w1 <- e * (1 - e)^2; w2 <- e^2 * (1 - e); w3 <- e^3
    T_k  <- trans_prob(o1, o2, k1, k2, p)
    T_kc <- trans_prob_pair(o1, o2, k1, k2, c1, c2, p)
    num <- w0 * T_kc + w1 * (T_k + T_c + Po) + (3 * w2 + w3) * Po
    den <- if (type == "pair") Po else w0 * T_k + w1 * (T_k + 2 * Po) + (3 * w2 + w3) * Po
  }
  list(num = num, den = den)
}

.safe_log_ratio <- function(num, den) {
  n <- max(length(num), length(den))
  num <- rep_len(num, n)
  den <- rep_len(den, n)
  out <- rep(NA_real_, n)
  ok <- !is.na(num) & !is.na(den)
  both0 <- ok & num == 0 & den == 0
  n0    <- ok & num == 0 & den > 0
  d0    <- ok & num > 0  & den == 0
  good  <- ok & num > 0  & den > 0
  out[both0] <- 0
  out[n0]    <- -Inf
  out[d0]    <- Inf
  out[good]  <- log(num[good]) - log(den[good])
  out
}

# Mendelian incompatibility of a candidate with an offspring, ignoring error.
# Used for mismatch counting and for exclusion-based screening.
.mismatch_locus <- function(o1, o2, c1, c2, p, k1 = NULL, k2 = NULL,
                            type = "one_known") {
  z <- if (type == "both_unknown") {
    trans_prob(o1, o2, c1, c2, p)
  } else {
    trans_prob_pair(o1, o2, k1, k2, c1, c2, p)
  }
  out <- z == 0
  out[is.na(z)] <- NA
  out
}
