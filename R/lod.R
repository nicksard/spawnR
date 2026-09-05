#' Multilocus LOD scores for candidate parents
#'
#' Computes the total LOD score, summed over independent loci, for every
#' offspring by candidate combination requested, in any of the three parentage
#' configurations of Kalinowski et al. (2007): a candidate parent evaluated with
#' the other parent's genotype known, a candidate parent evaluated alone, and a
#' candidate parent pair evaluated jointly.
#'
#' Loci at which any genotype required by the configuration is missing contribute
#' nothing and are excluded from `n_compared`. A positive LOD means the candidate
#' is more likely to be a true parent of the offspring than a random individual
#' drawn from the population; the score alone carries no significance, which is
#' supplied by [lod_critical()] or [delta_critical()].
#'
#' `mismatches` counts loci that are Mendelian-incompatible ignoring genotyping
#' error, and is reported for screening and quality control rather than used in
#' the likelihood.
#'
#' @param g A `"genotypes"` object holding all individuals referred to below.
#' @param offspring Identifiers (or positions) of the offspring to test.
#' @param candidates Identifiers of the candidate parents. For `type = "pair"`
#'   this may instead be a two-column `data.frame` or `matrix` whose columns give
#'   candidate father and candidate mother identifiers respectively, in which
#'   case `mothers` is ignored.
#' @param freqs An `"allele_freqs"` object, from [allele_freqs()].
#' @param error Genotyping error rate: a single value applied to every locus, or
#'   one value per locus (see [estimate_error_rates()]).
#' @param known For `type = "one_known"`, the identifier of the known parent: a
#'   single value recycled over all offspring, or one value per offspring, with
#'   `NA` allowed for offspring whose second parent is unknown.
#' @param type One of `"both_unknown"` (candidate parent alone),
#'   `"one_known"` (candidate parent with the other parent's genotype supplied
#'   through `known`), or `"pair"` (candidate parent pair jointly).
#' @param mothers For `type = "pair"`, identifiers of candidate mothers; every
#'   combination of `candidates` and `mothers` is evaluated.
#' @param exclude_self Drop combinations in which a candidate is the offspring
#'   itself or the known parent.
#' @param max_mismatch If not `NULL`, rows with more than this many mismatching
#'   loci are dropped from the result. Useful to keep output small when
#'   screening many candidates.
#'
#' @return A `data.frame` of class `"parentage_lod"` with one row per evaluated
#'   combination and columns `offspring`, `candidate` (and `candidate2` for
#'   `type = "pair"`), `n_compared`, `mismatches` and `lod`, sorted by offspring
#'   and then decreasing LOD.
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
#' sim <- simulate_population(n = 40, n_loci = 8, n_alleles = 6, seed = 42)
#' g <- sim$genotypes
#' p <- allele_freqs(g)
#' res <- parentage_lod(g, offspring = sim$pedigree$offspring[1:3],
#'                      candidates = sim$sires, freqs = p, error = 0.01,
#'                      type = "both_unknown")
#' head(res)
#' @export
parentage_lod <- function(g, offspring, candidates, freqs, error,
                          known = NULL,
                          type = c("both_unknown", "one_known", "pair"),
                          mothers = NULL, exclude_self = TRUE,
                          max_mismatch = NULL) {
  stopifnot(inherits(g, "genotypes"))
  type <- match.arg(type)
  .check_freqs(freqs, g)
  nl <- n_loci(g)
  e <- .expand_error(error, nl)

  off <- .as_ids(g, offspring, "offspring")
  no <- length(off)

  pair_df <- NULL
  if (type == "pair") {
    if (is.data.frame(candidates) || (is.matrix(candidates) && ncol(candidates) == 2L)) {
      pair_df <- data.frame(candidate = as.character(candidates[[1L]]),
                            candidate2 = as.character(candidates[[2L]]),
                            stringsAsFactors = FALSE)
    } else {
      if (is.null(mothers)) stop("`mothers` is required for type = 'pair' unless `candidates` is a two-column table.", call. = FALSE)
      cf <- .as_ids(g, candidates, "candidates")
      cm <- .as_ids(g, mothers, "mothers")
      pair_df <- expand.grid(candidate = cf, candidate2 = cm,
                             KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    }
    .as_ids(g, pair_df$candidate, "candidates")
    .as_ids(g, pair_df$candidate2, "mothers")
    cand <- pair_df$candidate
    cand2 <- pair_df$candidate2
  } else {
    cand <- .as_ids(g, candidates, "candidates")
    cand2 <- NULL
  }
  nc <- length(cand)
  if (nc == 0L) stop("No candidates supplied.", call. = FALSE)

  kn <- if (type == "one_known") {
    if (is.null(known)) stop("`known` is required for type = 'one_known'.", call. = FALSE)
    k <- as.character(known)
    if (length(k) == 1L) k <- rep(k, no)
    if (length(k) != no) stop("`known` must have length 1 or length(offspring).", call. = FALSE)
    .as_ids(g, k[!is.na(k)], "known", allow_na = TRUE)
    k
  } else rep(NA_character_, no)

  ci  <- match(cand, g$ids)
  ci2 <- if (is.null(cand2)) NULL else match(cand2, g$ids)
  ca1 <- g$a1[ci, , drop = FALSE]; ca2 <- g$a2[ci, , drop = FALSE]
  cb1 <- if (is.null(ci2)) NULL else g$a1[ci2, , drop = FALSE]
  cb2 <- if (is.null(ci2)) NULL else g$a2[ci2, , drop = FALSE]

  chunks <- vector("list", no)
  for (i in seq_len(no)) {
    oi <- match(off[i], g$ids)
    o1 <- g$a1[oi, ]; o2 <- g$a2[oi, ]

    if (type == "one_known") {
      if (is.na(kn[i])) {
        kk1 <- rep(NA_integer_, nl); kk2 <- kk1
      } else {
        kj <- match(kn[i], g$ids)
        kk1 <- g$a1[kj, ]; kk2 <- g$a2[kj, ]
      }
    }

    lodm <- matrix(NA_real_, nc, nl)
    mism <- matrix(NA, nc, nl)
    for (l in seq_len(nl)) {
      p <- freqs[[l]]
      if (type == "both_unknown") {
        lodm[, l] <- lod_locus(o1[l], o2[l], ca1[, l], ca2[, l], p = p,
                               error = e[l], type = "both_unknown")
        mism[, l] <- .mismatch_locus(o1[l], o2[l], ca1[, l], ca2[, l], p,
                                     type = "both_unknown")
      } else if (type == "one_known") {
        lodm[, l] <- lod_locus(o1[l], o2[l], ca1[, l], ca2[, l], p = p,
                               error = e[l], k1 = kk1[l], k2 = kk2[l],
                               type = "one_known")
        mism[, l] <- .mismatch_locus(o1[l], o2[l], ca1[, l], ca2[, l], p,
                                     k1 = kk1[l], k2 = kk2[l], type = "one_known")
      } else {
        lodm[, l] <- lod_locus(o1[l], o2[l], ca1[, l], ca2[, l], p = p,
                               error = e[l], k1 = cb1[, l], k2 = cb2[, l],
                               type = "pair")
        mism[, l] <- .mismatch_locus(o1[l], o2[l], ca1[, l], ca2[, l], p,
                                     k1 = cb1[, l], k2 = cb2[, l], type = "pair")
      }
    }
    ok <- !is.na(lodm)
    d <- data.frame(offspring = off[i], candidate = cand,
                    stringsAsFactors = FALSE)
    if (!is.null(cand2)) d$candidate2 <- cand2
    if (type == "one_known") d$known <- kn[i]
    d$n_compared <- rowSums(ok)
    d$mismatches <- rowSums(mism, na.rm = TRUE)
    lodm[!ok] <- 0
    d$lod <- rowSums(lodm)
    d$lod[d$n_compared == 0L] <- NA_real_

    drop <- rep(FALSE, nc)
    if (exclude_self) {
      drop <- drop | cand == off[i]
      if (!is.null(cand2)) drop <- drop | cand2 == off[i] | cand2 == cand
      if (type == "one_known" && !is.na(kn[i])) drop <- drop | cand == kn[i]
    }
    if (!is.null(max_mismatch)) drop <- drop | d$mismatches > max_mismatch
    chunks[[i]] <- d[!drop, , drop = FALSE]
  }

  out <- do.call(rbind, chunks)
  rownames(out) <- NULL
  out <- out[order(match(out$offspring, off), -.na_last(out$lod)), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "type") <- type
  attr(out, "error") <- e
  class(out) <- c("parentage_lod", "data.frame")
  out
}

.na_last <- function(x) {
  x[is.na(x)] <- -Inf
  x
}

.as_ids <- function(g, ids, what, allow_na = FALSE) {
  if (is.numeric(ids) || is.logical(ids)) return(g$ids[ids])
  ids <- as.character(ids)
  k <- match(ids, g$ids)
  bad <- is.na(k) & !(allow_na & is.na(ids))
  if (any(bad)) {
    stop("`", what, "` contains identifiers not present in the genotype data: ",
         paste(utils::head(unique(ids[bad]), 5L), collapse = ", "), call. = FALSE)
  }
  ids
}

#' @param x A `"parentage_lod"` object.
#' @param ... Passed to `print.data.frame`.
#' @rdname parentage_lod
#' @export
print.parentage_lod <- function(x, ...) {
  cat("<parentage_lod>", nrow(x), "comparisons,",
      length(unique(x$offspring)), "offspring, type =", attr(x, "type"), "\n")
  print.data.frame(utils::head(as.data.frame(x), 10L), ...)
  if (nrow(x) > 10L) cat("... ", nrow(x) - 10L, " more rows\n", sep = "")
  invisible(x)
}

#' Delta: the gap between the best and second-best candidate
#'
#' The discriminant statistic of Marshall et al. (1998): the difference in LOD
#' score between the most likely candidate parent and the second most likely.
#' When a single candidate remains after any exclusion-based screening, Delta is
#' defined as that candidate's LOD score.
#'
#' Only candidates with a positive LOD are treated as contenders, matching the
#' convention that a candidate less likely than a random individual is not a
#' meaningful runner-up.
#'
#' @param x A `"parentage_lod"` object from [parentage_lod()].
#'
#' @return A `data.frame` with one row per offspring: `offspring`, the best
#'   candidate (`candidate`, and `candidate2` where the configuration is a pair),
#'   `lod`, `lod2` (the runner-up's score, `NA` if there is none), `delta`,
#'   `n_candidates` (contenders with positive LOD) and `mismatches`.
#'
#' @references
#' Marshall, T.C., Slate, J., Kruuk, L.E.B. & Pemberton, J.M. (1998)
#' Statistical confidence for likelihood-based paternity inference in natural
#' populations. \emph{Molecular Ecology}, \strong{7}, 639--655.
#' \doi{10.1046/j.1365-294x.1998.00374.x}
#'
#' @examples
#' sim <- simulate_population(n = 40, n_loci = 8, n_alleles = 6, seed = 42)
#' p <- allele_freqs(sim$genotypes)
#' res <- parentage_lod(sim$genotypes, sim$pedigree$offspring[1:3],
#'                      sim$sires, p, error = 0.01, type = "both_unknown")
#' delta_stat(res)
#' @export
delta_stat <- function(x) {
  stopifnot(inherits(x, "parentage_lod"))
  has_pair <- !is.null(x$candidate2)
  sp <- split(seq_len(nrow(x)), factor(x$offspring, levels = unique(x$offspring)))
  rows <- lapply(sp, function(idx) {
    d <- x[idx, , drop = FALSE]
    d <- d[!is.na(d$lod), , drop = FALSE]
    if (nrow(d) == 0L) return(NULL)
    d <- d[order(-d$lod), , drop = FALSE]
    pos <- sum(d$lod > 0)
    top <- d[1L, ]
    lod2 <- if (pos >= 2L) d$lod[2L] else NA_real_
    delta <- if (pos >= 2L) top$lod - d$lod[2L] else if (pos == 1L) top$lod else NA_real_
    out <- data.frame(offspring = top$offspring, candidate = top$candidate,
                      stringsAsFactors = FALSE)
    if (has_pair) out$candidate2 <- top$candidate2
    out$lod <- top$lod
    out$lod2 <- lod2
    out$delta <- delta
    out$n_candidates <- pos
    out$mismatches <- top$mismatches
    out
  })
  out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  rownames(out) <- NULL
  out
}
