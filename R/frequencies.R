#' Estimate allele frequencies
#'
#' Allele frequencies are estimated by simple allele counting over all typed
#' individuals at each locus. Because a likelihood-ratio calculation divides by
#' Hardy-Weinberg genotype frequencies, an allele that is observed in a candidate
#' or offspring but has estimated frequency zero would make the likelihood
#' undefined; `min_freq` guards against that by flooring every frequency before
#' renormalising.
#'
#' Marshall et al. (1998) discuss the sensitivity of LOD scores to the reference
#' allele frequencies and recommend estimating them from the study population
#' itself, which is what `reference` allows when the individuals used for
#' frequency estimation are not the same as those being tested.
#'
#' @param g A `"genotypes"` object.
#' @param reference Optional subset of `g` used to estimate frequencies: a
#'   character vector of identifiers, or integer/logical positions. Defaults to
#'   all individuals.
#' @param min_freq Lower bound applied to every estimated frequency before
#'   renormalisation. The default is `min(1 / (2 * n), 0.01)` with `n` the
#'   number of typed individuals at the locus: the frequency of a single
#'   unobserved allele copy, capped at 1%. The cap matters for small reference
#'   sets. Uncapped, six typed parents would floor every allele at 1/12, about
#'   0.083, so an allele never seen in the reference would count as common as
#'   one seen once, which distorts LOD scores and the Delta simulation, since
#'   [sim_delta()] draws genotypes from these frequencies. A message is shown
#'   when fewer than 25 individuals are typed at a locus, because frequencies
#'   from so few allele copies are noisy whatever the floor; supply `min_freq`
#'   explicitly to choose the floor and silence it. Set to `0` for unmodified
#'   counting.
#'
#' @return An object of class `"allele_freqs"`: a list of numeric vectors, one
#'   per locus, each summing to one and named by allele label.
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
#' @examples
#' g <- simulate_population(n = 60, n_loci = 5, n_alleles = 6, seed = 1)$genotypes
#' p <- allele_freqs(g)
#' round(p[[1]], 3)
#' @export
allele_freqs <- function(g, reference = NULL, min_freq = NULL) {
  stopifnot(inherits(g, "genotypes"))
  gg <- if (is.null(reference)) g else subset_ind(g, reference)
  out <- vector("list", n_loci(g))
  names(out) <- g$loci
  inv2n <- rep(NA_real_, n_loci(g))
  for (l in seq_len(n_loci(g))) {
    lab <- g$alleles[[l]]
    k <- length(lab)
    if (k == 0L) stop("Locus '", g$loci[l], "' has no typed genotypes.", call. = FALSE)
    cnt <- tabulate(c(gg$a1[, l], gg$a2[, l]), nbins = k)
    ntyped <- sum(!is.na(gg$a1[, l]))
    mf <- if (is.null(min_freq)) {
      if (ntyped > 0L) min(1 / (2 * ntyped), 0.01) else 1 / (2 * k)
    } else min_freq
    inv2n[l] <- if (ntyped > 0L) 1 / (2 * ntyped) else Inf
    p <- if (sum(cnt) > 0) cnt / sum(cnt) else rep(1 / k, k)
    p <- pmax(p, mf)
    p <- p / sum(p)
    names(p) <- lab
    out[[l]] <- p
  }
  if (is.null(min_freq) && any(inv2n > 0.02)) {
    big <- inv2n > 0.02
    nmin <- if (is.finite(max(inv2n))) round(1 / (2 * max(inv2n))) else 0L
    message("allele_freqs(): ", sum(big), " of ", length(inv2n), " loci have ",
            "fewer than 25 typed individuals in the reference set (as few as ",
            nmin, "). Frequencies from so few allele copies are noisy, which ",
            "affects LOD scores and the Delta simulation; unobserved alleles are ",
            "floored at 0.01. Consider a larger reference set or external ",
            "frequencies, or set `min_freq` explicitly to silence this.")
  }
  structure(out, class = c("allele_freqs", "list"))
}

#' @param x An `"allele_freqs"` object.
#' @param ... Ignored.
#' @rdname allele_freqs
#' @export
print.allele_freqs <- function(x, ...) {
  k <- vapply(x, length, integer(1))
  cat("<allele_freqs>", length(x), "loci,", sum(k), "alleles total\n")
  cat("  alleles per locus: ", paste0(range(k), collapse = "-"), "\n", sep = "")
  cat("  loci: ", paste(utils::head(names(x), 6L), collapse = ", "),
      if (length(x) > 6L) ", ..." else "", "\n", sep = "")
  invisible(x)
}

.check_freqs <- function(p, g) {
  if (!inherits(p, "allele_freqs") && !is.list(p)) {
    stop("`freqs` must be an `allele_freqs` object (see ?allele_freqs).", call. = FALSE)
  }
  if (length(p) != n_loci(g)) {
    stop("`freqs` covers ", length(p), " loci but the genotypes have ", n_loci(g), ".", call. = FALSE)
  }
  for (l in seq_along(p)) {
    if (length(p[[l]]) < length(g$alleles[[l]])) {
      stop("`freqs` for locus '", g$loci[l], "' has fewer alleles than the genotype data. ",
           "Estimate frequencies from the same `genotypes` object.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

.expand_error <- function(error, nl, what = "error") {
  if (is.null(error)) stop("`", what, "` must be supplied.", call. = FALSE)
  e <- as.numeric(error)
  if (length(e) == 1L) e <- rep(e, nl)
  if (length(e) != nl) {
    stop("`", what, "` must be a single value or one value per locus (", nl, ").", call. = FALSE)
  }
  if (anyNA(e) || any(e < 0) || any(e > 1)) {
    stop("`", what, "` values must lie in [0, 1].", call. = FALSE)
  }
  e
}
