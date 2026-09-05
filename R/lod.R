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
#' @param ibd Identity-by-descent coefficients `c(k0, k1, k2)` describing how a
#'   *non-parent* candidate is related to the offspring, changing what the LOD
#'   score is a ratio against. The default `NULL` means unrelated, `c(1, 0, 0)`,
#'   which is the classical assumption. Use [ibd_mixture()] to build this from
#'   relationship classes. Not yet supported for `type = "pair"`.
#' @param relatedness Named weights describing the relationship classes present
#'   among *non-parent* candidates, for example
#'   `c(unrelated = 0.95, avuncular = 0.05)`; see [ibd_mixture()] for the class
#'   names. The LOD then becomes a ratio against the best-fitting alternative
#'   rather than against an unrelated stranger, which is what keeps a close
#'   relative of an unsampled parent from being mistaken for the parent. Mutually
#'   exclusive with `ibd`. Not yet supported for `type = "pair"`.
#' @param output `"long"` returns one row per comparison. `"matrix"` returns
#'   the offspring-by-candidate matrices directly, which is much cheaper for
#'   large candidate sets: a thousand offspring against several thousand
#'   candidates is millions of rows in long form, and assembling that table
#'   costs more than computing the scores.
#'
#' @return A `data.frame` of class `"parentage_lod"` with one row per evaluated
#'   combination and columns `offspring`, `candidate` (and `candidate2` for
#'   `type = "pair"`), `n_compared`, `mismatches` and `lod`, sorted by offspring
#'   and then decreasing LOD. With `output = "matrix"`, an object of class
#'   `"parentage_lod_matrix"`: a list of `lod`, `mismatches` and `n_compared`
#'   matrices with offspring in rows and candidates in columns.
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
                          max_mismatch = NULL, ibd = NULL, relatedness = NULL,
                          output = c("long", "matrix")) {
  stopifnot(inherits(g, "genotypes"))
  type <- match.arg(type)
  output <- match.arg(output)
  if (!is.null(relatedness)) {
    if (!is.null(ibd)) stop("Give either `ibd` or `relatedness`, not both.", call. = FALSE)
    if (type == "pair") stop("`relatedness` is not yet supported for type = 'pair'.", call. = FALSE)
    rw <- relatedness / sum(relatedness)
    if (any(rw < 0) || is.null(names(rw))) {
      stop("`relatedness` must be a named vector of non-negative class weights; ",
           "see ?ibd_mixture for the class names.", call. = FALSE)
    }
  }
  if (!is.null(ibd)) {
    ibd <- as.numeric(ibd)
    if (length(ibd) != 3L || any(ibd < 0) || abs(sum(ibd) - 1) > 1e-8) {
      stop("`ibd` must be three non-negative coefficients summing to one; ",
           "see ?ibd_mixture.", call. = FALSE)
    }
    if (type == "pair") stop("`ibd` is not yet supported for type = 'pair'.", call. = FALSE)
  }
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

  gi <- .geno_index(g)
  oi <- match(off, g$ids)
  ci <- match(cand, g$ids)

  MIS <- NULL; NCM <- NULL
  if (type == "pair") {
    c2i <- match(cand2, g$ids)
    LOD <- MIS <- NCM <- matrix(0, no, nc)
    for (m in unique(cand2)) {
      j <- which(cand2 == m)
      M <- .lod_matrix(gi, oi, ci[j], rep(match(m, g$ids), no), freqs, e, "pair")
      LOD[, j] <- M$lod; MIS[, j] <- M$mism; NCM[, j] <- M$ncmp
    }
  } else {
    ki <- if (type == "one_known") match(kn, g$ids) else NULL
    if (is.null(relatedness)) {
      M <- .lod_matrix(gi, oi, ci, ki, freqs, e, type, ibd = ibd)
      LOD <- M$lod; MIS <- M$mism; NCM <- M$ncmp
    } else {
      # Relatedness is a property of the whole genome, not an independent draw
      # at every locus, so the classes must be mixed over multilocus
      # likelihoods, not inside the per-locus denominator. Writing L_c for the
      # LOD against class c alone,
      #     LOD = -log( sum_c w_c exp(-L_c) ),
      # which is dominated by whichever alternative explains the pair best.
      # Mixing per locus instead shifts every candidate almost equally and
      # leaves true parents and their relatives just as hard to tell apart.
      acc <- NULL
      for (cl in names(rw)) {
        Mc <- .lod_matrix(gi, oi, ci, ki, freqs, e, type,
                          ibd = ibd_mixture(structure(1, names = cl)))
        term <- log(rw[[cl]]) - Mc$lod
        acc <- if (is.null(acc)) term else {
          mx <- pmax(acc, term)
          mx + log(exp(acc - mx) + exp(term - mx))
        }
        if (is.null(MIS)) { MIS <- Mc$mism; NCM <- Mc$ncmp }
      }
      LOD <- -acc
    }
  }

  if (output == "matrix") {
    dn <- list(off, if (is.null(cand2)) cand else paste(cand, cand2, sep = " x "))
    lod_m <- LOD; lod_m[NCM == 0L] <- NA_real_
    return(structure(list(lod = `dimnames<-`(lod_m, dn),
                          mismatches = `dimnames<-`(MIS, dn),
                          n_compared = `dimnames<-`(NCM, dn),
                          type = type, error = e),
                     class = "parentage_lod_matrix"))
  }

  # Long form, assembled in one pass rather than one data.frame per offspring.
  d <- data.frame(offspring = rep(off, each = nc),
                  candidate = rep(cand, times = no),
                  stringsAsFactors = FALSE)
  if (!is.null(cand2)) d$candidate2 <- rep(cand2, times = no)
  if (type == "one_known") d$known <- rep(kn, each = nc)
  d$n_compared <- as.integer(t(NCM))
  d$mismatches <- as.integer(t(MIS))
  lodv <- as.vector(t(LOD))
  lodv[d$n_compared == 0L] <- NA_real_
  d$lod <- lodv

  drop <- rep(FALSE, nrow(d))
  if (exclude_self) {
    drop <- drop | d$candidate == d$offspring
    if (!is.null(cand2)) drop <- drop | d$candidate2 == d$offspring | d$candidate2 == d$candidate
    if (type == "one_known") drop <- drop | (!is.na(d$known) & d$candidate == d$known)
  }
  if (!is.null(max_mismatch)) drop <- drop | d$mismatches > max_mismatch
  chunks <- list(d[!drop, , drop = FALSE])

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

#' @param x A `"parentage_lod_matrix"` object.
#' @param ... Ignored.
#' @rdname parentage_lod
#' @export
print.parentage_lod_matrix <- function(x, ...) {
  cat("<parentage_lod_matrix>", nrow(x$lod), "offspring x", ncol(x$lod),
      "candidates, type =", x$type, "\n")
  cat("  LOD range: ", paste(format(range(x$lod, na.rm = TRUE), digits = 4),
                             collapse = " to "), "\n", sep = "")
  invisible(x)
}
