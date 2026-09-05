#' Mendelian transition and genotype probabilities
#'
#' Elementary probabilities underlying every likelihood ratio in the package,
#' vectorised over individuals at a single locus. Alleles are supplied as
#' integer codes indexing `p`; `NA` propagates.
#'
#' `geno_freq()` is the Hardy-Weinberg genotype frequency
#' \eqn{P(o) = p_x^2} for a homozygote and \eqn{2 p_x p_y} for a heterozygote.
#'
#' `trans_prob()` is \eqn{T(o \mid f)}, the probability of the offspring
#' genotype given one parent's genotype when the other parent is an unsampled
#' individual drawn from the population. The known parent transmits each of its
#' alleles with probability one half and the unsampled parent contributes an
#' allele at its population frequency, giving the standard table
#' \eqn{p_i}, \eqn{p_j}, \eqn{(p_i + p_j)/2}, \eqn{p_k/2}, and zero when no
#' allele is shared.
#'
#' `trans_prob_pair()` is \eqn{T(o \mid m, f)}, the probability of the offspring
#' genotype given both parents, equal to one quarter of the number of the four
#' allele pairings that reconstruct the offspring genotype.
#'
#' Two identities that these functions satisfy exactly, and that the package's
#' regression tests verify numerically, are
#' \eqn{\sum_F P(F) T(o \mid m, F) = T(o \mid m)} and
#' \eqn{\sum_M \sum_F P(M) P(F) T(o \mid M, F) = P(o)}.
#'
#' @param o1,o2 Integer allele codes of the offspring genotype.
#' @param f1,f2 Integer allele codes of the (candidate) parent genotype.
#' @param m1,m2 Integer allele codes of the second parent genotype.
#' @param p Numeric vector of allele frequencies at the locus, indexed by allele
#'   code.
#'
#' @return A numeric vector of probabilities, `NA` where any input allele is
#'   `NA`.
#'
#' @references
#' Marshall, T.C., Slate, J., Kruuk, L.E.B. & Pemberton, J.M. (1998)
#' Statistical confidence for likelihood-based paternity inference in natural
#' populations. \emph{Molecular Ecology}, \strong{7}, 639--655.
#' \doi{10.1046/j.1365-294x.1998.00374.x}
#'
#' @examples
#' p <- c(0.5, 0.3, 0.2)
#' # heterozygous parent 1/2, heterozygous offspring 1/2
#' trans_prob(1L, 2L, 1L, 2L, p)          # (p1 + p2) / 2 = 0.4
#' trans_prob_pair(1L, 2L, 1L, 1L, 2L, 2L, p)  # certain: 1
#' geno_freq(1L, 2L, p)                   # 2 * 0.5 * 0.3 = 0.3
#' @export
geno_freq <- function(o1, o2, p) {
  n <- max(length(o1), length(o2))
  o1 <- rep_len(o1, n); o2 <- rep_len(o2, n)
  out <- rep(NA_real_, n)
  ok <- !is.na(o1) & !is.na(o2)
  hom <- ok & o1 == o2
  het <- ok & o1 != o2
  out[hom] <- p[o1[hom]]^2
  out[het] <- 2 * p[o1[het]] * p[o2[het]]
  out
}

#' @rdname geno_freq
#' @export
trans_prob <- function(o1, o2, f1, f2, p) {
  n <- max(length(o1), length(o2), length(f1), length(f2))
  o1 <- rep_len(o1, n); o2 <- rep_len(o2, n)
  f1 <- rep_len(f1, n); f2 <- rep_len(f2, n)
  ok <- !is.na(o1) & !is.na(o2) & !is.na(f1) & !is.na(f2)
  hom <- ok & o1 == o2
  het <- ok & o1 != o2
  h <- function(a) {
    out <- numeric(n)
    i <- which(hom & a == o1); out[i] <- p[o1[i]]
    i <- which(het & a == o1); out[i] <- p[o2[i]]
    i <- which(het & a == o2); out[i] <- p[o1[i]]
    out
  }
  out <- 0.5 * (h(f1) + h(f2))
  out[!ok] <- NA_real_
  out
}

#' @rdname geno_freq
#' @export
trans_prob_pair <- function(o1, o2, m1, m2, f1, f2, p) {
  n <- max(length(o1), length(o2), length(m1), length(m2), length(f1), length(f2))
  o1 <- rep_len(o1, n); o2 <- rep_len(o2, n)
  m1 <- rep_len(m1, n); m2 <- rep_len(m2, n)
  f1 <- rep_len(f1, n); f2 <- rep_len(f2, n)
  ok <- !is.na(o1) & !is.na(o2) & !is.na(m1) & !is.na(m2) & !is.na(f1) & !is.na(f2)
  hit <- function(a, b) {
    v <- (a == o1 & b == o2) | (a == o2 & b == o1)
    v[is.na(v)] <- FALSE
    as.numeric(v)
  }
  out <- 0.25 * (hit(m1, f1) + hit(m1, f2) + hit(m2, f1) + hit(m2, f2))
  out[!ok] <- NA_real_
  out
}
