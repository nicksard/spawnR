#' Marker informativeness and non-exclusion probabilities
#'
#' The per-locus summary that CERVUS reports alongside its allele frequency
#' analysis: number of alleles, observed and expected heterozygosity,
#' polymorphic information content, and five non-exclusion probabilities, per
#' locus and combined over loci. It says how much power a panel has before any
#' parentage is attempted, and what a count of mismatching loci from
#' [parentage_lod()] or [assign_trios()] can be expected to achieve.
#'
#' @section Non-exclusion probabilities:
#' Each is the probability that an individual who is *not* the parent (or not
#' the same individual) is nevertheless *not excluded* by Mendelian
#' inheritance at the locus, for codominant markers without null alleles,
#' mutation or genotyping error, under Hardy-Weinberg proportions. Writing
#' \eqn{a_i = \sum_j p_j^i} for the power sums of the allele frequencies, the
#' exclusion probabilities are those of Jamieson & Taylor (1997):
#' \describe{
#'   \item{`ne_1p`}{First parent: an unrelated candidate tested against the
#'     offspring alone, the other parent unknown (their eq. 2a),
#'     \eqn{1 - 4a_2 + 2a_2^2 + 4a_3 - 3a_4}.}
#'   \item{`ne_2p`}{Second parent: an unrelated candidate tested when the other
#'     parent's genotype is known (their eq. 1a, the classical paternity
#'     exclusion), \eqn{1 - 2a_2 + a_3 + 2a_4 - 3a_5 - 2a_2^2 + 3a_2 a_3}.}
#'   \item{`ne_pp`}{Parent pair: an unrelated mother-father pair tested jointly
#'     (their eq. 3a), \eqn{1 + 4a_4 - 4a_5 - 3a_6 - 8a_2^2 + 8a_2 a_3 +
#'     2a_3^2}.}
#' }
#' and the non-exclusion probability is one minus each. The two identity
#' probabilities are the chance that two individuals share a genotype:
#' \describe{
#'   \item{`ne_i`}{Two unrelated individuals, \eqn{2a_2^2 - a_4}.}
#'   \item{`ne_si`}{Two full sibs, \eqn{1/4 + a_2/2 + a_2^2/2 - a_4/4}.}
#' }
#' Every one of these is checked in the test suite against direct enumeration
#' over all genotypes. Combined values multiply across loci, which assumes the
#' loci are independent: linked SNPs overstate the combined power.
#'
#' The naming follows CERVUS: NE-1P, NE-2P, NE-PP, NE-I and NE-SI. Note that
#' Jamieson & Taylor number their equations the other way round: their eq. 1a
#' (paternity, mother known) is NE-2P here, and their eq. 2a (one parent
#' unavailable) is NE-1P.
#'
#' @section Other columns:
#' `h_obs` is the proportion of typed individuals that are heterozygous, and
#' `h_exp` the unbiased expected heterozygosity of Nei (1987),
#' \eqn{2n/(2n-1)\,(1 - a_2)} for \eqn{n} typed individuals. `pic` is the
#' polymorphic information content of Botstein et al. (1980),
#' \eqn{1 - a_2 - \sum_{i<j} 2 p_i^2 p_j^2}. `f_is` is \eqn{1 - h_{obs} /
#' h_{exp}}: strongly positive values suggest null alleles or allelic dropout,
#' and strongly negative values, excess heterozygosity, which in species with
#' a recent genome duplication (salmonids among them) often marks a paralogous
#' locus that will cause spurious mismatches.
#'
#' @param g A `"genotypes"` object.
#' @param reference Optional subset of individuals from which frequencies and
#'   heterozygosities are computed: identifiers, or integer/logical positions.
#'   Typically the candidate parents. Defaults to everyone in `g`.
#' @param freqs Optional `"allele_freqs"` object to use instead of counting
#'   alleles in `reference`. By default frequencies are counted directly with
#'   no floor, as CERVUS does; the floor applied by [allele_freqs()] protects
#'   the likelihood from zero frequencies but would bias these summaries.
#'
#' @return A `data.frame` of class `"marker_summary"` with one row per locus,
#'   columns `locus`, `n_alleles`, `n_typed`, `h_obs`, `h_exp`, `pic`, `f_is`,
#'   `ne_1p`, `ne_2p`, `ne_pp`, `ne_i` and `ne_si`, and attribute `combined`:
#'   the mean number of alleles, mean heterozygosities and PIC, and the
#'   combined non-exclusion probabilities, the product over loci.
#'
#' @references
#' Jamieson, A. & Taylor, St.C.S. (1997) Comparisons of three probability
#' formulae for parentage exclusion. \emph{Animal Genetics}, \strong{28},
#' 397--400.
#'
#' Marshall, T.C., Slate, J., Kruuk, L.E.B. & Pemberton, J.M. (1998)
#' Statistical confidence for likelihood-based paternity inference in natural
#' populations. \emph{Molecular Ecology}, \strong{7}, 639--655.
#' \doi{10.1046/j.1365-294x.1998.00374.x}
#'
#' Botstein, D., White, R.L., Skolnick, M. & Davis, R.W. (1980) Construction
#' of a genetic linkage map in man using restriction fragment length
#' polymorphisms. \emph{American Journal of Human Genetics}, \strong{32},
#' 314--331.
#'
#' Nei, M. (1987) \emph{Molecular Evolutionary Genetics}. Columbia University
#' Press, New York.
#'
#' @examples
#' sim <- simulate_population(n = 60, n_loci = 8, n_alleles = 6, seed = 1)
#' ms <- marker_summary(sim$genotypes, reference = c(sim$dams, sim$sires))
#' ms
#' attr(ms, "combined")
#' @export
marker_summary <- function(g, reference = NULL, freqs = NULL) {
  stopifnot(inherits(g, "genotypes"))
  gg <- if (is.null(reference)) g else subset_ind(g, reference)
  nl <- n_loci(g)
  if (!is.null(freqs)) .check_freqs(freqs, g)

  rows <- lapply(seq_len(nl), function(l) {
    k <- length(g$alleles[[l]])
    typed <- !is.na(gg$a1[, l])
    n <- sum(typed)
    p <- if (!is.null(freqs)) {
      as.numeric(freqs[[l]])
    } else {
      cnt <- tabulate(c(gg$a1[typed, l], gg$a2[typed, l]), nbins = k)
      if (n > 0L) cnt / sum(cnt) else rep(NA_real_, k)
    }
    p <- p[!is.na(p) & p > 0]
    if (!length(p)) {
      return(data.frame(locus = g$loci[l], n_alleles = 0L, n_typed = n,
                        h_obs = NA_real_, h_exp = NA_real_, pic = NA_real_,
                        f_is = NA_real_, ne_1p = NA_real_, ne_2p = NA_real_,
                        ne_pp = NA_real_, ne_i = NA_real_, ne_si = NA_real_))
    }
    a2 <- sum(p^2); a3 <- sum(p^3); a4 <- sum(p^4); a5 <- sum(p^5); a6 <- sum(p^6)
    h_obs <- if (n > 0L) mean(gg$a1[typed, l] != gg$a2[typed, l]) else NA_real_
    h_exp <- if (n > 1L) 2 * n / (2 * n - 1) * (1 - a2) else 1 - a2
    pic <- 1 - a2 - (a2^2 - a4)        # sum_{i<j} 2 p_i^2 p_j^2 = a2^2 - a4
    data.frame(
      locus = g$loci[l], n_alleles = length(p), n_typed = n,
      h_obs = h_obs, h_exp = h_exp, pic = pic,
      f_is = if (is.na(h_obs) || h_exp <= 0) NA_real_ else 1 - h_obs / h_exp,
      ne_1p = 1 - (1 - 4 * a2 + 2 * a2^2 + 4 * a3 - 3 * a4),
      ne_2p = 1 - (1 - 2 * a2 + a3 + 2 * a4 - 3 * a5 - 2 * a2^2 + 3 * a2 * a3),
      ne_pp = 1 - (1 + 4 * a4 - 4 * a5 - 3 * a6 - 8 * a2^2 + 8 * a2 * a3 + 2 * a3^2),
      ne_i  = 2 * a2^2 - a4,
      ne_si = 0.25 + 0.5 * a2 + 0.5 * a2^2 - 0.25 * a4)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL

  ne <- c("ne_1p", "ne_2p", "ne_pp", "ne_i", "ne_si")
  comb <- c(mean_alleles = mean(out$n_alleles),
            mean_h_obs = mean(out$h_obs, na.rm = TRUE),
            mean_h_exp = mean(out$h_exp, na.rm = TRUE),
            mean_pic = mean(out$pic, na.rm = TRUE),
            vapply(ne, function(v) prod(out[[v]], na.rm = TRUE), numeric(1)))
  attr(out, "combined") <- comb
  class(out) <- c("marker_summary", "data.frame")
  out
}

#' @param x A `"marker_summary"` object.
#' @param ... Passed to `print.data.frame`.
#' @rdname marker_summary
#' @export
print.marker_summary <- function(x, ...) {
  cat("<marker_summary>", nrow(x), "loci\n")
  print.data.frame(utils::head(as.data.frame(x), 20L), digits = 4,
                   row.names = FALSE, ...)
  if (nrow(x) > 20L) cat("... ", nrow(x) - 20L, " more loci\n", sep = "")
  cmb <- attr(x, "combined")
  if (!is.null(cmb)) {
    cat("\nCombined over loci (non-exclusion = product over loci):\n")
    cat(sprintf("  mean alleles %.2f, mean H_obs %.3f, mean H_exp %.3f, mean PIC %.3f\n",
                cmb[["mean_alleles"]], cmb[["mean_h_obs"]], cmb[["mean_h_exp"]],
                cmb[["mean_pic"]]))
    cat(sprintf("  NE-1P %.3g   NE-2P %.3g   NE-PP %.3g   NE-I %.3g   NE-SI %.3g\n",
                cmb[["ne_1p"]], cmb[["ne_2p"]], cmb[["ne_pp"]], cmb[["ne_i"]],
                cmb[["ne_si"]]))
  }
  invisible(x)
}
