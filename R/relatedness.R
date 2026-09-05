#' Average IBD coefficients for a mixture of relationship classes
#'
#' Builds the `ibd` argument used by [parentage_lod()] and [lod_locus()] from a
#' description of who is in the candidate pool. The coefficients describe how a
#' candidate that is *not* the parent is related to the offspring.
#'
#' @section Why an average is enough:
#' The likelihood under the non-parent hypothesis is linear in `(k0, k1, k2)`,
#' so a mixture over relationship classes has exactly the likelihood of the
#' mixture-averaged coefficient vector. You therefore never need to enumerate
#' the classes at run time; you only need their average.
#'
#' @section What the classes mean:
#' Coefficients are for the **candidate-to-offspring** relationship, which is
#' one step further than the relationship people usually have in mind. A
#' candidate who is a full sib of the true father is the offspring's uncle, so
#' the relevant sharing is avuncular, not full-sib.
#'
#' \describe{
#'   \item{`unrelated`}{`c(1, 0, 0)`}
#'   \item{`cousin`}{`c(0.75, 0.25, 0)` -- also a half sib of the true parent}
#'   \item{`avuncular`}{`c(0.5, 0.5, 0)` -- a full sib of the true parent, a
#'     grandparent of the offspring, or a half sib of the offspring}
#'   \item{`full_sib`}{`c(0.25, 0.5, 0.25)` -- another offspring of the same
#'     two parents}
#'   \item{`parent`}{`c(0, 1, 0)` -- given for completeness; using it makes
#'     every LOD score zero, since it is the numerator hypothesis}
#' }
#'
#' @param ... Named weights for relationship classes, for example
#'   `unrelated = 0.9, avuncular = 0.1`. Weights are normalised to sum to one.
#'
#' @return A numeric vector of length three, `c(k0, k1, k2)`.
#'
#' @seealso [parentage_lod()], [pool_relatedness()]
#'
#' @examples
#' ibd_mixture(unrelated = 1)
#' ibd_mixture(unrelated = 0.9, avuncular = 0.1)
#' ibd_mixture(unrelated = 0.8, cousin = 0.15, avuncular = 0.05)
#' @export
ibd_mixture <- function(...) {
  w <- unlist(list(...))
  tab <- list(unrelated = c(1, 0, 0), cousin = c(0.75, 0.25, 0),
              avuncular = c(0.5, 0.5, 0), full_sib = c(0.25, 0.5, 0.25),
              parent = c(0, 1, 0))
  if (length(w) == 0L) return(tab$unrelated)
  nm <- names(w)
  if (is.null(nm) || any(!nzchar(nm))) stop("All weights must be named.", call. = FALSE)
  bad <- setdiff(nm, names(tab))
  if (length(bad)) {
    stop("Unknown relationship class(es): ", paste(bad, collapse = ", "),
         ". Available: ", paste(names(tab), collapse = ", "), ".", call. = FALSE)
  }
  if (any(w < 0)) stop("Weights must be non-negative.", call. = FALSE)
  if (sum(w) <= 0) stop("Weights must not all be zero.", call. = FALSE)
  w <- w / sum(w)
  out <- colSums(do.call(rbind, lapply(nm, function(k) tab[[k]])) * w)
  names(out) <- c("k0", "k1", "k2")
  out
}

#' Screen a candidate pool for close relatives
#'
#' Estimates how many candidates carry more allele sharing with each offspring
#' than an unrelated individual plausibly would, which is the quantity that
#' decides whether the `relatedness` argument of [parentage_lod()] is needed.
#'
#' For each offspring a reference distribution of LOD scores for genuinely
#' unrelated candidates is generated analytically from the allele frequencies,
#' conditioned on that offspring's own genotypes and missing-data pattern. The
#' number of candidates exceeding the `1 - alpha` quantile of that distribution
#' is compared with the `alpha * n` expected by chance, and the excess is
#' reported as an estimated count of relatives.
#'
#' Read the result as a **lower bound**, not an estimate. Only relatives extreme
#' enough to clear the threshold are counted, and many are not: in a pool seeded
#' with ten full sibs of the true sire plus the sire himself, 5.5% of the pool by
#' construction, this returns about 3%. It also cannot distinguish a full sib of
#' the parent from a grandparent, and the true parent counts as an excess, so
#' subtract one where you expect the parent to be present.
#'
#' Use it to decide *whether* a relatedness-aware alternative is warranted, then
#' set the classes in [parentage_lod()] from what you know about how the pool was
#' collected, not from this number.
#'
#' @param g A `"genotypes"` object.
#' @param offspring Identifiers of offspring to screen.
#' @param candidates Identifiers of the candidate pool.
#' @param freqs An `"allele_freqs"` object.
#' @param error Genotyping error rate, one value or one per locus.
#' @param alpha Upper-tail probability defining "more sharing than expected".
#' @param nsim Replicates for the unrelated reference distribution.
#'
#' @return A `data.frame` with one row per offspring: `offspring`, the
#'   `threshold` LOD, `n_above`, the `expected` count under no relatedness, and
#'   `excess_frac`, the estimated fraction of the pool that is related.
#'
#' @seealso [ibd_mixture()], [parentage_lod()]
#'
#' @examples
#' sim <- simulate_population(n = 200, n_loci = 12, n_alleles = 8,
#'                            n_offspring = 10, error = 0.01, seed = 12)
#' p <- allele_freqs(sim$genotypes)
#' pool_relatedness(sim$genotypes, sim$pedigree$offspring[1:5], sim$sires,
#'                  p, error = 0.01, nsim = 500)
#' @export
pool_relatedness <- function(g, offspring, candidates, freqs, error,
                             alpha = 0.01, nsim = 2000L) {
  stopifnot(inherits(g, "genotypes"))
  off <- .as_ids(g, offspring, "offspring")
  L <- parentage_lod(g, off, candidates, freqs, error, output = "matrix")$lod
  nc <- ncol(L)
  out <- lapply(seq_along(off), function(i) {
    nd <- .sim_unrelated(.row_of(g, off[i], "offspring"), freqs,
                         .expand_error(error, n_loci(g)), nsim)
    thr <- stats::quantile(nd, 1 - alpha, names = FALSE)
    v <- L[i, ]; v <- v[!is.na(v)]
    n_above <- sum(v > thr)
    data.frame(offspring = off[i], threshold = thr, n_above = n_above,
               expected = alpha * length(v),
               excess_frac = max(0, (n_above - alpha * length(v)) / length(v)),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

# LOD distribution for candidates drawn at random from the population, i.e.
# genuinely unrelated to this offspring.
.sim_unrelated <- function(og, freqs, e, nsim) {
  nl <- length(freqs); parts <- vector("list", nl); m <- 0L
  for (l in seq_len(nl)) {
    o1 <- og$a1[l]; o2 <- og$a2[l]
    if (is.na(o1)) next
    p <- freqs[[l]]; tb <- .locus_tab(length(p))
    Phwe <- geno_freq(tb$gi, tb$gj, p)
    vals <- lod_locus(o1, o2, tb$gi, tb$gj, p = p, error = e[l],
                      type = "both_unknown")
    vals[!is.finite(vals)] <- 0
    m <- m + 1L
    parts[[m]] <- list(vals = vals, probs = Phwe)
  }
  .draw_sum(parts[seq_len(m)], nsim)
}
